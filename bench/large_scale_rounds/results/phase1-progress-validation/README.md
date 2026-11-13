# 第一阶段只读诊断验证

在 `53ab2c4` 的独立工作树验证默认关闭的 `phase1_progress_diagnostics`。诊断在线程内显式启用，只认领首个实际运行的第一阶段。引擎构造来源登记最多八份，采集最多四份完整快照和 4096 条摘要事件。浮点数据以原始 IEEE 754 位串输出。`estimated_bytes` 是诊断数据的估算预算，计入预留的事件与快照容量、来源和捕获数据。256 MiB 不代表分配器实际占用或进程常驻内存的硬上限。

核对了稠密与稀疏路径。观察点位于第一阶段入口、人工变量正部和首次到达可行容差且负人工未超差时、首次失败以及人工清理前。定价记录来自当轮已有结果，按状态世代标记适用范围。观察读取求解器状态，额外扫描只计算诊断摘要，不调用定价、重构或求解，也不新增求解提前退出。显式阶段界、第二阶段残基人工固定为零、主元执行前安全复核和全基本变量界守卫沿用基线。工作树没有运行公开模型性能基准或 `bnatt500`。

初次编译的八项专项均通过。静态复核发现诊断原有估算漏计摘要事件向量容量，修正后也显式计入四份快照的固定容量，并增加预算断言。稀疏重算检测到非有限值时，同步记录首次故障原因和局部重算值，避免恢复后终态覆盖该证据。另补两种混合行序断言，验证原有汇总在任意行出现非有限基本值时始终保留 NaN。验证过程没有失败的测试。

| 日志 | 执行内容 | 结果 |
| --- | --- | --- |
| `01-specialized-initial.log` | `cargo test --manifest-path milp/Cargo.toml --features phase1_progress_diagnostics --lib phase1_progress_diagnostics::tests` | 初版 8 通过 |
| `02-specialized-budget.log` | 同上 | 诊断预算修正后 8 通过 |
| `03-feature-debug.log` | `cargo test --manifest-path milp/Cargo.toml --features phase1_progress_diagnostics` | 库 478 通过、4 忽略，CLI 29 通过，验收 14 通过、2 忽略，冷热对账 1 通过 |
| `04-feature-release-acceptance.log` | `cargo test --manifest-path milp/Cargo.toml --features phase1_progress_diagnostics --release --test acceptance -- --include-ignored` | 验收 16 通过，含常规 debug 忽略的 2 项 |
| `07-default-example-build.log` | `cargo build --manifest-path milp/Cargo.toml --example phase1_trajectory_probe` | 默认特性构建成功 |
| `09-example-validation.log` | 对账两种构建的轨迹文本并解析两个 JSON | 轨迹逐字节相同，两个引擎均捕获入口、首个近零状态、清理前状态 |
| `10-specialized-nonfinite-order.log` | 同首行命令 | 补充混合行序断言后 8 通过，生产代码与全量验证时一致 |

例程仅使用自造的 `two-equalities.milp`。在带特性的 debug 构建后运行 `phase1_trajectory_probe 输入.milp 1000`，随后在默认构建运行同一命令，比较 `05-trajectory-enabled.txt` 与 `08-trajectory-default.txt`。`phase1_progress_probe 输入.milp 输出.json 30 dense` 和 `phase1_progress_probe 输入.milp 输出.json 30 revised` 各运行一次，分别保存为 `06-capture-dense.json` 与 `06-capture-revised.json`。
