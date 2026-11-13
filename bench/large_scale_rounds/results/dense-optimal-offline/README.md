# p0201 已完成的离线复核归档

本目录仅整理已完成计算，没有重跑。主报告见 `mic-dense-optimal-audit-report.md`；完整结果见 `mic-dense-optimal-exact-audit.json`；计算脚本见 `mic_dense_optimal_exact_audit.py`。

输入 artifact：`/tmp/mic-dense-failure-profile/bench/large_scale_rounds/results/dense-optimal-profile/p0201-root.stdout`

SHA256：`19350aa7eaeae3177390ae103a0966dbe39adc9375149e82eeb6a5544b27d7d5`。诊断源码 d200934，证据提交 f90e33e。只验证 `normalized_current_lp` 保存的 binary64 输入；不假定外部 MIP 坐标映射。

- 完整样本只有 Optimal 序号 88–91；没有 104、183、184 的完整快照。
- 精确核验确认保存松弛值不一致；最大异常对应原不等式仍有明显余量，不能称原行不可行。
- 微小负松弛检验数导致负无穷箱下界，这是合法无穷端点结果，不是非有限算术错误；它不足以证明返回界错误。
- 已做有限的逐分量有理化候选尝试，未得到严格可行点。目标为 7227 的候选仍违反若干行，不能用作反证。
- 后续精确可行点修复/反证、对偶符号投影、183/184 定点补采均未执行；留待 root 后续 CPU 窗口。没有矩阵求解、优化器调用或生产代码修改。

脚本保留计算时的绝对输入/输出路径；如搬迁归档，运行前需调整路径。R11 计时期间不要重跑。
