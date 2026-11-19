# bnatt500 稀疏第一阶段冷恢复定点

默认关闭的 `sparse_phase1_recovery` 捕获了同一个坏基，并让初始根 LP 正常结束。完整求解仍在根节点用完 30 秒，返回 `LimitReached`、界 0，无整数解。随后首批根割试算仍出现另一处 LU 失败。

两组共用冻结源码 `ab838a9`，其求解代码是 `b234f51`，继承了 `14d9a22` 已正式开启的根割复查。只给候选组开启恢复特性。单线程，30 秒与 2000000 节点，自动选 LP 引擎和二元传播，停滞早停与探测复用关闭。两个探针都把原模型重新导出，与外部基准归档逐字节匹配。完整来源及命令参数见 [manifest.json](manifest.json)。

- **生产默认**在初始根 LP 第 3904 次迭代处因 LU 主元不足退出。候选绝对值为 0，因子位置 4950。最终 1 节点，未产生对偶界、割或可行解。
- **恢复候选**从相同首错重启，初始根 LP 和同参数根重启都在 7211 次迭代后正常完成，界为 0。首批 10 条割另遇热启动中的 LU 第 5052 列主元不足，原回滚路径恢复了根界，最终保留 0 条。根启发式恢复 LP 明确返回 `time_limit`，仍有第一阶段欠账。最终 1 节点、90675 次累计 LP 迭代、界 0。冻结原 MPS 最小化非负的 x1（转换件 x0），可直接由变量下界验证 0。
- 两组的生产与只读诊断输出分别对账，各自 23 个非耗时字段完全一致。两个当前首错捕获还与 [67dd6ee 的旧定点捕获](../bnatt-current-lu/capture.json.gz) 逐字节相同，共用其矩阵和换基日志。旧报告中全程开启风险刷新曾保留 50 条根割，其实验范围和生产默认均不同。
- 健康小例 `assign-4x4` 显式选稀疏引擎，默认与候选的 23 个非耗时字段和完整向量相同，36 次 LP 迭代，最优值与界均为 10。四个行最小成本 1、2、3、4 分别位于不同列，报告向量同时达到这四个下界。

恢复仅用于直接冷启动第一阶段中已执行换基所触发的有限 LU 主元不足。坏基所在尝试的迭代与工作步、绝对期限、当前界、原行、割和反循环记账跨重建保留。第二遍第一阶段复用风险主元刷新，每次 LP 调用只可冷重建一次，人工清理和热路径仍沿原数值处理。已运行但无可用因子的 resolve 在读旧因子前直接冷建，同样受此约束。风险刷新内重分解即便失败，也保留数值原因并遵守重分解后期限检查。候选保持默认关闭，本次为定点功能证据。

## 复核和重采

```sh
archive=bench/large_scale_rounds/results/bnatt-phase1-recovery
(cd "$archive" && shasum -a 256 -c SHA256SUMS)
python3 "$archive/verify.py"
```

在冻结源码工作树用 `cargo build --manifest-path milp/Cargo.toml --release --example round_probe` 构建各组，按清单增加特性。诊断组额外应用 [现场采集补丁](../bnatt-current-lu/capture-probe.patch)，并通过 `MIC_LU_CAPTURE` 写出原始现场。诊断 JSON 的外层耗时包含求解后的离线重放。生产探针不加诊断补丁，求解报告耗时保留用于还原现场。

重采时用清单参数执行 `round_probe solve`。将诊断组输出命名为 `bnatt-default-capture.json`、`bnatt-candidate-capture.json`。分别用两组生产二进制执行 `round_probe export`，输出 `default-original.json`、`candidate-original.json`。四个文件放在同一目录，可用 `verify.py --live 目录` 复核其字节。健康小例沿用生产参数，额外指定 `--lp-engine revised`。

候选单开专项 8 项、叠加工作限额及诊断等特性的专项 9 项、默认完整 MILP 539 项、组合完整 MILP 596 项、候选 release 验收 16 项、planner 127 项通过。默认和组合完整测试各保留原 6 项忽略，planner 保留原 2 项忽略。

```sh
combined=sparse_phase1_recovery,root_heuristic_work_limit,lu_failure_diagnostics,phase1_progress_diagnostics,risk_pivot_refresh,proof_profile
cargo test --manifest-path milp/Cargo.toml --features sparse_phase1_recovery --lib phase1_cold_recovery_tests
cargo test --manifest-path milp/Cargo.toml --features "$combined" --lib phase1_cold_recovery_tests
cargo test --manifest-path milp/Cargo.toml
cargo test --manifest-path milp/Cargo.toml --features "$combined"
cargo test --manifest-path milp/Cargo.toml --release --features sparse_phase1_recovery --test acceptance -- --include-ignored
cargo test --manifest-path planner/Cargo.toml
```
