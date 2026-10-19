# micpython × 开源求解器（HiGHS / SCIP）同机对比

环境：conda `mic-bench`（py3.11，pip 的 scipy/ortools 官方轮子），threads=1，时限见正文；机器 15W 笔记本（墙钟噪声 ±20-60%，见 docs/parallel-profiling.md §1.3）。solve 单位秒；`watchdog_kill` = 求解器时限+120s 宽限内未返回被强杀。CBC 列已按决策移除（2026-10-15，理由见脚本头注释：挂死不受时限约束 + 三例错误 Unbounded 判定）。

| 算例 | micpython | highs | scip | 目标一致 |
|---|---|---|---|---|
| knapsack-200 | 0.002/optimal | 0.034/optimal | 0.012/optimal | ✓ |
| knapsack-500 | 0.011/optimal | 0.054/optimal | 0.011/optimal | ✓ |
| knapsack-2000 | 0.011/optimal | 0.128/optimal | 0.127/optimal | ✓ |
| knapsack-5000 | 0.02/optimal | 0.361/optimal | 0.044/optimal | △ |
| mknap-3x100 | 0.077/optimal | 0.624/optimal | 0.205/optimal | ✓ |
| mknap-5x300 | 7.093/optimal | 2.954/optimal | 2.113/optimal | ✓ |
| setcover-100x150 | 0.228/optimal | 0.212/optimal | 0.248/optimal | ✓ |
| setcover-300x600 | 21.95/optimal | 1.707/optimal | 1.715/optimal | ✓ |
| milp-mixed-250i40c | 0.006/optimal | 0.462/optimal | 0.284/optimal | ✓ |
| marketsplit-4x30 | 0.01/optimal | 6.135/optimal | 5.281/optimal | ✓ |
| marketsplit-4x40 | 8.488/optimal | 60.01/limit | 59.528/status6 | ✓ |
| assign-40x40 | 0.08/optimal | 0.024/optimal | 0.021/optimal | ✓ |
| lp-wide-200x300 | 0.008/optimal | 0.027/optimal | 0.124/optimal | ✓ |
| lp-sparse-4000x800 | 5.291/optimal | 0.887/optimal | 1.912/optimal | ✓ |

## MISMATCH（目标不一致，需人工核查）

- knapsack-5000: [('micpython', 88526.0, 'optimal'), ('highs', 88526.0, 'optimal'), ('scip', 88525.0, 'optimal')]
- knapsack-5000 [DP]: [('scip', 88525.0, 88526.0)]
