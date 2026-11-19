# bnatt500 当前基线的稀疏基分解失败

生产参数下，最后 64 次换基中的第 25 次把原本非奇异的基变成了精确奇异的基。这次发生在 LP 迭代计数为 3865 时，步长是零，在基位置 1906 用列 8 换出列 14。生产求解沿用换基乘子更新后的方向，算出的主元为 `1.2300387375833728e-9`，刚超过原 `1e-9` 门槛。同一换基前的基从原矩阵重分解后，主元算为 `-0.0`。独立 SuperLU 用两种列序分别算出约 `-8.05e-11` 和 `1.21e-12`。生产求解直到第 64 次换基后的周期重构才发现基坏了，在因子位置 4950 看到候选绝对值为零，累计报告 3904 次迭代。前一份证据中的 2688 步故障有另一条换基轨迹。

证书在 [verify.py](verify.py)。19 个非零整数行权重作为左向量，作用到换基后的 5058 个基列时，按保存系数的精确二进制有理数计算全部为零。作用到换出的列 14 时，点积为 `7205759403792785/18014398509481984`。因此对换基前的基方程 `B α = A[:,8]` 左乘这个向量，基位置 1906 上的真实方向分量必为零。在线方向该分量的误差把一次本不满足主元门槛的零步换基放了过去。生产重分解重放显示前 24 次换基后仍能完成分解，换入列 8 后至最终失败基均不能。独立 SuperLU 换基后以 COLAMD 列序报奇异，用 NATURAL 列序时最小绝对 U 主元约 `6.01e-16`。

默认构建与实际开启采集的构建各运行一次，全部非耗时字段精确相同，均为 `LimitReached`、1 个节点、3904 次 LP 迭代、无可行解和界。单独开启已有的实验特性 `risk_pivot_refresh`，仍用原生产参数，30 秒到限，未捕获 LU 分解失败，做了 105584 次 LP 迭代、保留 50 条根割、返回界 0。该次运行未获得可行解或完成证明，也未逐次核对是否经过第 3865 次换基对应的基。这里只记录候选机制的定点现象。

## 独立复核

首个失败的原始捕获只保存一次。压缩包包含最后 64 次换基及其涉及的全部矩阵列，SHA-256 钉住原始未压缩内容。用标准库复核证书，加 `--check-prior` 可进一步用模素数消元验证失效前的 25 个基在保存系数下都非奇异。先决条件只是矩阵系数为精确二进制分数，其分母均可在该素数域中求逆。

```sh
evidence=bench/large_scale_rounds/results/bnatt-current-lu
(cd "$evidence" && shasum -a 256 -c SHA256SUMS)
python3 "$evidence/verify.py"
python3 "$evidence/verify.py" --check-prior
gzip -dc "$evidence/capture.json.gz" > /tmp/mic-bnatt-current-offline.json
python3 bench/inspect_lu_failure.py /tmp/mic-bnatt-current-offline.json --scipy
```

最后一项使用独立环境中的 SciPy 1.16.3 与 NumPy 2.5.5。旧的生产算术历史重放以及阈值都随捕获保存。证书只涉及一个基的秩与该方向分量，不提供 LP 或整数模型可行性的判定。

## 现场重采

源码为 `67dd6ee43bca13eaa063a8b549648461636a3e58`，MILP 树为 `23d628b71db35e3e5b846dcc66457842e5926402`。构建使用 aarch64-apple-darwin 的 rustc 1.97.0 (`2d8144b78`)。冻结转换件 `bnatt500.milp` 的 SHA-256 为 `80f3109a1af2fe652334df4e3f0dbb9d7c597fc790697ede2bb51fc46e3de77c`，外部基准原模型导出的摘要为 `7d48f7d203013102ca06fca727c334b200c32d5cbbbfe87ccaa35a265bbd7359`。

在该提交的独立工作树应用 [探针补丁](capture-probe.patch)。补丁只接通原 `round_probe` 的诊断，在原生产求解结束后做历史重放，保持所有参数与判断。将下列 `model` 替换为恢复的冻结转换件。

```sh
model=/tmp/mic-rounds-miplib2017/bnatt500.milp
cargo build --manifest-path milp/Cargo.toml --release --example round_probe
milp/target/release/examples/round_probe solve "$model" --early-stop off --probe-reuse off --binary-propagation auto --threads 1 --time-limit-secs 30 --node-limit 2000000
cargo build --manifest-path milp/Cargo.toml --release --features lu_failure_diagnostics --example round_probe
MIC_LU_CAPTURE=/tmp/mic-bnatt-current-capture.json milp/target/release/examples/round_probe solve "$model" --early-stop off --probe-reuse off --binary-propagation auto --threads 1 --time-limit-secs 30 --node-limit 2000000
python3 /path/to/evidence/verify.py /tmp/mic-bnatt-current-capture.json --check-prior
```

对照的三个 release 可执行文件摘要分别为默认 `9fca64e33c640560f85a741e5e162a53db4910f5feb991cbf30fe05a75321844`、采集 `bb55fbd05cb4be23f0eb88b9472d0634a470faa305d36f4fc2d69384376cd112`、采集加 `risk_pivot_refresh` `42215295fa242b2afb63d9ac90076338bc0e9b87dd2ef2b2d7fd774d422108b2`。原始报告中的耗时仅保留为现场记录。
