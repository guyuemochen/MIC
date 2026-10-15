# micpython × 开源求解器（CBC / HiGHS / SCIP）同机对比

环境：conda `mic-bench`（py3.11，pip 的 pulp/scipy/ortools 官方轮子），threads=1，时限见正文；机器 15W 笔记本（墙钟噪声 ±20-60%，见 docs/parallel-profiling.md §1.3）。solve 单位秒；`watchdog_kill` = 求解器时限+120s 宽限内未返回被强杀（实测 CBC 在 mknap/marketsplit 类会进入不受时限约束的阶段）。

| 算例 | micpython | CBC | HiGHS | SCIP | 目标一致 |
|---|---|---|---|---|---|
| knapsack-200 | 0.003/optimal | 0.167/optimal | 0.034/optimal | 0.013/optimal | ✓ |
| knapsack-500 | 0.019/optimal | 0.095/optimal | 0.058/optimal | 0.009/optimal | ✓ |
| knapsack-2000 | 0.012/optimal | 0.114/optimal | 0.127/optimal | 0.137/optimal | ✓ |
| knapsack-5000 | 0.023/optimal | 0.172/optimal | 0.398/optimal | 0.062/optimal | △ |
| mknap-3x100 | 0.086/optimal | 180.0/watchdog_kill | 0.662/optimal | 0.233/optimal | ✓ |
| mknap-5x300 | 10.357/optimal | 240.0/watchdog_kill | 3.123/optimal | 2.457/optimal | ✓ |
| setcover-100x150 | 0.313/optimal | 0.239/optimal | 0.229/optimal | 0.286/optimal | ✓ |
| setcover-300x600 | 57.71/optimal | 240.0/watchdog_kill | 1.796/optimal | 1.782/optimal | ✓ |
| milp-mixed-250i40c | 0.013/optimal | 0.157/Unbounded | 0.5/optimal | 0.311/optimal | △ |
| marketsplit-4x30 | 0.01/optimal | 180.0/watchdog_kill | 7.195/optimal | 6.232/optimal | ✓ |
| marketsplit-4x40 | 9.082/optimal | 180.0/watchdog_kill | 60.011/limit | 59.405/status6 | ✓ |
| assign-40x40 | 0.111/optimal | 0.179/optimal | 0.027/optimal | 0.024/optimal | ✓ |
| lp-wide-200x300 | 0.01/optimal | 0.152/Unbounded | 0.031/optimal | 0.126/optimal | △ |
| lp-sparse-4000x800 | 9.156/optimal | 0.832/Unbounded | 0.895/optimal | 2.095/optimal | △ |

## MISMATCH（目标不一致，需人工核查）

- knapsack-5000: [('micpython', 88526.0, 'optimal'), ('cbc', 88526.0, 'optimal'), ('highs', 88526.0, 'optimal'), ('scip', 88525.0, 'optimal')]
- knapsack-5000 [DP]: [('scip', 88525.0, 88526.0)]
- milp-mixed-250i40c: [('micpython', 1275.791751746084, 'optimal'), ('cbc', 445078867824000.0, 'Unbounded'), ('highs', 1275.7917517463113, 'optimal'), ('scip', 1275.791751746084, 'optimal')]
- lp-wide-200x300: [('micpython', 112.66861478047609, 'optimal'), ('cbc', 23475.870994099987, 'Unbounded'), ('highs', 112.66861478047613, 'optimal'), ('scip', 112.66861478047622, 'optimal')]
- lp-sparse-4000x800: [('micpython', 5893.705756618241, 'optimal'), ('cbc', 360949672.7894396, 'Unbounded'), ('highs', 5893.705755757289, 'optimal'), ('scip', 5893.705756618257, 'optimal')]
