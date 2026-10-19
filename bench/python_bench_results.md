# micpython × 开源求解器（HiGHS / SCIP）同机对比

环境：conda `mic-bench`（py3.11，pip 的 scipy/ortools 官方轮子），threads=1，时限见正文；机器 15W 笔记本（墙钟噪声 ±20-60%，见 docs/parallel-profiling.md §1.3）。solve 单位秒；`watchdog_kill` = 求解器时限+120s 宽限内未返回被强杀。CBC 列已按决策移除（2026-10-15，理由见脚本头注释：挂死不受时限约束 + 三例错误 Unbounded 判定）。

| 算例 | micpython | highs | scip | 目标一致 |
|---|---|---|---|---|
| knapsack-200 | 0.007/optimal | 0.109/optimal | 0.019/optimal | ✓ |
| knapsack-500 | 0.016/optimal | 0.153/optimal | 0.024/optimal | ✓ |
| knapsack-2000 | 0.024/optimal | 0.286/optimal | 0.293/optimal | ✓ |
| knapsack-5000 | 0.034/optimal | 0.808/optimal | 0.125/optimal | △ |
| mknap-3x100 | 0.134/optimal | 1.21/optimal | 0.493/optimal | ✓ |
| mknap-5x300 | 9.784/optimal | 6.148/optimal | 4.776/optimal | ✓ |
| setcover-100x150 | 0.45/optimal | 0.426/optimal | 0.601/optimal | ✓ |
| setcover-300x600 | 37.687/optimal | 3.441/optimal | 3.739/optimal | ✓ |
| milp-mixed-250i40c | 0.013/optimal | 0.897/optimal | 0.728/optimal | ✓ |
| marketsplit-4x30 | 0.02/optimal | 13.938/optimal | 12.004/optimal | ✓ |
| marketsplit-4x40 | 18.119/optimal | 60.016/limit | 59.948/status6 | ✓ |
| assign-40x40 | 0.215/optimal | 0.052/optimal | 0.045/optimal | ✓ |
| lp-wide-200x300 | 0.021/optimal | 0.057/optimal | 0.283/optimal | ✓ |
| lp-sparse-4000x800 | 10.161/optimal | 1.965/optimal | 4.204/optimal | ✓ |

## MISMATCH（目标不一致，需人工核查）

- knapsack-5000: [('micpython', 88526.0, 'optimal'), ('highs', 88526.0, 'optimal'), ('scip', 88525.0, 'optimal')]
- knapsack-5000 [DP]: [('scip', 88525.0, 88526.0)]
