# micpython × 开源求解器（CBC / HiGHS / SCIP）同机对比

环境：conda `mic-bench`（py3.11，pip 的 pulp/scipy/ortools 官方轮子），threads=1，时限见正文；机器 15W 笔记本（墙钟噪声 ±20-60%，见 docs/parallel-profiling.md §1.3）。solve 单位秒；`watchdog_kill` = 求解器时限+120s 宽限内未返回被强杀（实测 CBC 在 mknap/marketsplit 类会进入不受时限约束的阶段）。

| 算例 | micpython | CBC | HiGHS | SCIP | 目标一致 |
|---|---|---|---|---|---|
| knapsack-200 | 0.002/optimal | 0.161/optimal | 0.034/optimal | 0.005/optimal | ✓ |
| knapsack-500 | 0.02/optimal | 0.1/optimal | 0.065/optimal | 0.012/optimal | ✓ |
| knapsack-2000 | 0.015/optimal | 0.117/optimal | 0.144/optimal | 0.137/optimal | ✓ |
| knapsack-5000 | 0.026/optimal | 0.176/optimal | 0.421/optimal | 0.043/optimal | △ |
| mknap-3x100 | 0.081/optimal | 180.0/watchdog_kill | 0.791/optimal | 0.216/optimal | ✓ |
| mknap-5x300 | 10.849/optimal | 240.0/watchdog_kill | 2.947/optimal | 2.225/optimal | ✓ |
| setcover-100x150 | 0.365/optimal | 0.232/optimal | 0.208/optimal | 0.314/optimal | ✓ |
| setcover-300x600 | 59.102/optimal | 240.0/watchdog_kill | 1.838/optimal | 1.957/optimal | ✓ |
| milp-mixed-250i40c | 0.013/optimal | 0.161/Unbounded | 0.523/optimal | 0.287/optimal | △ |
| marketsplit-4x30 | 0.014/optimal | 180.0/watchdog_kill | 6.544/optimal | 19.995/optimal | ✓ |
| marketsplit-4x40 | 8.555/optimal | 180.0/watchdog_kill | 60.008/limit | 59.476/status6 | ✓ |
| assign-40x40 | 0.102/optimal | 0.212/optimal | 0.035/optimal | 0.033/optimal | ✓ |
| lp-wide-200x300 | 0.011/optimal | 0.187/Unbounded | 0.035/optimal | 0.145/optimal | △ |
| lp-sparse-4000x800 | 34.757/optimal | 1.002/Unbounded | 1.036/optimal | 2.872/optimal | △ |

## MISMATCH（目标不一致，需人工核查）

- knapsack-5000: [('micpython', 88526.0, 'optimal'), ('cbc', 88526.0, 'optimal'), ('highs', 88526.0, 'optimal'), ('scip', 88525.0, 'optimal')]
- knapsack-5000 [DP]: [('scip', 88525.0, 88526.0)]
- milp-mixed-250i40c: [('micpython', 1275.791751746084, 'optimal'), ('cbc', 445078867824000.0, 'Unbounded'), ('highs', 1275.7917517463113, 'optimal'), ('scip', 1275.791751746084, 'optimal')]
- lp-wide-200x300: [('micpython', 112.66861478047609, 'optimal'), ('cbc', 23475.870994099987, 'Unbounded'), ('highs', 112.66861478047613, 'optimal'), ('scip', 112.66861478047622, 'optimal')]
- lp-sparse-4000x800: [('micpython', 5893.705756617977, 'optimal'), ('cbc', 360949672.7894396, 'Unbounded'), ('highs', 5893.705755757289, 'optimal'), ('scip', 5893.705756618257, 'optimal')]
