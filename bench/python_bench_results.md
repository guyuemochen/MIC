# micpython × 开源求解器（HiGHS / SCIP）同机对比

环境：conda `mic-bench`（py3.11，pip 的 scipy/ortools 官方轮子），threads=1，时限见正文；机器 15W 笔记本（墙钟噪声 ±20-60%，见 docs/parallel-profiling.md §1.3）。solve 单位秒；`watchdog_kill` = 求解器时限+120s 宽限内未返回被强杀。CBC 列已按决策移除（2026-10-15，理由见脚本头注释：挂死不受时限约束 + 三例错误 Unbounded 判定）。

| 算例 | micpython | highs | scip | 目标一致 |
|---|---|---|---|---|
| knapsack-200 | 0.002/optimal | 0.03/optimal | 0.003/optimal | ✓ |
| knapsack-500 | 0.02/optimal | 0.054/optimal | 0.009/optimal | ✓ |
| knapsack-2000 | 0.011/optimal | 0.121/optimal | 0.158/optimal | ✓ |
| knapsack-5000 | 0.026/optimal | 0.395/optimal | 0.042/optimal | △ |
| mknap-3x100 | 0.082/optimal | 0.655/optimal | 0.219/optimal | ✓ |
| mknap-5x300 | 13.655/optimal | 2.853/optimal | 2.091/optimal | ✓ |
| setcover-100x150 | 0.282/optimal | 0.206/optimal | 0.248/optimal | ✓ |
| setcover-300x600 | 70.348/optimal | 1.714/optimal | 1.698/optimal | ✓ |
| milp-mixed-250i40c | 0.012/optimal | 0.46/optimal | 0.295/optimal | ✓ |
| marketsplit-4x30 | 0.01/optimal | 7.01/optimal | 5.915/optimal | ✓ |
| marketsplit-4x40 | 9.078/optimal | 60.007/limit | 59.817/status6 | ✓ |
| assign-40x40 | 0.077/optimal | 0.023/optimal | 0.022/optimal | ✓ |
| lp-wide-200x300 | 0.009/optimal | 0.026/optimal | 0.106/optimal | ✓ |
| lp-sparse-4000x800 | 9.323/optimal | 0.784/optimal | 1.872/optimal | ✓ |

## MISMATCH（目标不一致，需人工核查）

- knapsack-5000: [('micpython', 88526.0, 'optimal'), ('highs', 88526.0, 'optimal'), ('scip', 88525.0, 'optimal')]
- knapsack-5000 [DP]: [('scip', 88525.0, 88526.0)]
