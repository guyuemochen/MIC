# PyBoolNet 原规则定点对照

17 个网络中，当前 MIC 的专用路径命中 1 个。旧 MIC、当前 MIC、Gurobi 10.0.3 都证明了 17 个模型的最优值。直接逐项执行原规则，穷举 87,032 个状态，共找到 61 个固定点。每个模型都至少有一个固定点，本组没有无解样例。

**样本来源。** 固定 PyBoolNet 的提交 `3394310d07e3e25810a828d8da4f6d7241115113`。30 个 BNet 文件里，有 29 个与所在目录同名的主模型。取 2–16 个变量的全部 19 个主模型，再剔除上游明确标注为多值编码的 `multivalued`，以及与 `dinwoodie_stomatal` 的规则完全相同的 `saadatpour_guardcell`。所有未入选文件和理由均记录在结果的 `skipped` 字段。[上游固定快照](https://github.com/hklarner/pyboolnet/tree/3394310d07e3e25810a828d8da4f6d7241115113/pyboolnet/repository)

**模型与核验。** 原文件给出各变量的布尔更新规则。令每个变量等于自己的下一步取值，就得到固定点问题。上游没有优化目标，本轮统一最小化第一条规则的输出位。对每个合取或析取门，以两条整数线性不等式表示其真值关系。原规则内的嵌套门使用辅助变量，复制和取反直接写成等式，恒等规则无需约束。穷举器独立执行源规则。每次求解都检查目标与界，检查返回解及回调解的全部取值、MILP 行和源规则。首轮的 51 个最终向量、52 个回调向量及原文件和转换模型的 SHA-256 均已归档。向量按变量声明顺序记录，开头为全部原始状态位。[逐例数据](results.json)

**独立命中。** `randomnet_n7k3` 是上游用 BoolNet 生成的随机网络。27 个 MILP 变量、54 行被压缩为 7 个状态位和 30 个子句。原规则共有 10 个固定点，最小目标为 0。旧 MIC 用了 261 次 LP 迭代，新 MIC 用了 0 次。每方预热一次，再轮换顺序串行运行 9 次。求解耗时中位数为旧 MIC 1.277 ms、当前 MIC 0.045 ms、Gurobi 0.428 ms。时间从各求解入口开始，含解回调，不计进程启动、读文件或 Gurobi 模型装载。其余 16 个样例均走其他路径，其中 15 个含专用识别器不接收的两列约束。剩余的 `randomnet_n15k3` 有输出位直接进入另一个析取门，如 `Gene1` 的输入 `Gene8`，超出了当前两层结构的范围。[命中模型原文](https://github.com/hklarner/pyboolnet/blob/3394310d07e3e25810a828d8da4f6d7241115113/pyboolnet/repository/randomnet_n7k3/randomnet_n7k3.bnet) · [15 位对照模型](https://github.com/hklarner/pyboolnet/blob/3394310d07e3e25810a828d8da4f6d7241115113/pyboolnet/repository/randomnet_n15k3/randomnet_n15k3.bnet)

**复现。** 在 MIC 仓库根目录执行以下命令。脚本直接读取固定版本的 Git 对象，原网络与转换后的 MILP 都保存在临时目录。旧 MIC 取 `31fbb74786f20a8d00a16ea65d9db77f7b2addbc`，当前 MIC 取 `6f8be8e7f9dcac311c514e8a81f2cfda3da3651a`。两版均使用 release 构建，单线程，每例限时 60 秒、节点上限 200 万。若本机没有 Gurobi，可省略 `--gurobi` 参数重放两版 MIC 对照。

```sh
work=$(mktemp -d)
audit=bench/large_scale_rounds/results/bnatt-external
git clone https://github.com/hklarner/pyboolnet.git "$work/pyboolnet"
for rev in 31fbb74 6f8be8e; do
  mkdir "$work/$rev"
  git archive "$rev" milp | tar -x -C "$work/$rev"
  cargo build --manifest-path "$work/$rev/milp/Cargo.toml" --locked --release --example round_probe
done
# 将观察脚本置于指定版本，固定识别器与当前 MIC 的源码。
mkdir -p "$work/6f8be8e/$audit"
cp "$audit/build_probe.py" "$work/6f8be8e/$audit/"
python3 "$work/6f8be8e/$audit/build_probe.py" "$work/detector"
python3 "$audit/independent/check.py" --repo "$work/pyboolnet" --models "$work/models" \
  --baseline "$work/31fbb74/milp/target/release/examples/round_probe" \
  --current "$work/6f8be8e/milp/target/release/examples/round_probe" \
  --detector "$work/detector" --gurobi /Library/gurobi1003/macos_universal2/lib/libgurobi100.dylib \
  --output "$work/results.json"
```
