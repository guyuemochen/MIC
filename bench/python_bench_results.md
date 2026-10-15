# micpython × 开源求解器（CBC / HiGHS / SCIP）同机对比

环境：conda `mic-bench`（py3.11，pip 的 pulp/scipy/ortools 官方轮子），threads=1，时限见正文；机器 15W 笔记本（墙钟噪声 ±20-60%，见 docs/parallel-profiling.md §1.3）。solve 单位秒；`watchdog_kill` = 求解器时限+120s 宽限内未返回被强杀（实测 CBC 在 mknap/marketsplit 类会进入不受时限约束的阶段）。

| 算例 | micpython | CBC | HiGHS | SCIP | 目标一致 |
|---|---|---|---|---|---|
| knapsack-200 | 0.003/optimal | 0.275/optimal | 0.057/optimal | 0.012/optimal | ✓ |
| knapsack-500 | 0.021/optimal | 0.196/optimal | 0.07/optimal | 0.014/optimal | ✓ |
| knapsack-2000 | 0.011/optimal | 0.109/optimal | 0.129/optimal | 0.161/optimal | ✓ |
| knapsack-5000 | 0.026/optimal | 0.178/optimal | 0.404/optimal | 0.066/optimal | △ |
| mknap-3x100 | 0.093/optimal | 180.0/watchdog_kill | 0.743/optimal | 0.262/optimal | ✓ |
| mknap-5x300 | 13.069/optimal | 240.0/watchdog_kill | 10.538/optimal | 5.559/optimal | ✓ |
| setcover-100x150 | 0.709/optimal | 0.344/optimal | 0.602/optimal | 0.667/optimal | ✓ |
| setcover-300x600 | 82.574/optimal | 240.0/watchdog_kill | 1.865/optimal | 1.942/optimal | ✓ |
| milp-mixed-250i40c | 0.014/optimal | 0.16/Unbounded | 0.52/optimal | 0.322/optimal | △ |
| marketsplit-4x30 | 0.01/optimal | 180.0/watchdog_kill | 7.269/optimal | 6.236/optimal | ✓ |
| marketsplit-4x40 | 9.873/optimal | 180.0/watchdog_kill | 60.01/limit | 59.047/status6 | ✓ |
| assign-40x40 | 0.087/optimal | 0.193/optimal | 0.032/optimal | 0.025/optimal | ✓ |
| lp-wide-200x300 | 0.009/optimal | 0.161/Unbounded | 0.029/optimal | 0.126/optimal | △ |
| lp-sparse-4000x800 | 10.046/optimal | 1.06/Unbounded | 1.183/optimal | 2.84/optimal | △ |

## MISMATCH（目标不一致，需人工核查）

- knapsack-5000: [('micpython', 88526.0, 'optimal'), ('cbc', 88526.0, 'optimal'), ('highs', 88526.0, 'optimal'), ('scip', 88525.0, 'optimal')]
- knapsack-5000 [DP]: [('scip', 88525.0, 88526.0)]
- milp-mixed-250i40c: [('micpython', 1275.791751746084, 'optimal'), ('cbc', 445078867824000.0, 'Unbounded'), ('highs', 1275.7917517463113, 'optimal'), ('scip', 1275.791751746084, 'optimal')]
- lp-wide-200x300: [('micpython', 112.66861478047609, 'optimal'), ('cbc', 23475.870994099987, 'Unbounded'), ('highs', 112.66861478047613, 'optimal'), ('scip', 112.66861478047622, 'optimal')]
- lp-sparse-4000x800: [('micpython', 5893.705756618241, 'optimal'), ('cbc', 360949672.7894396, 'Unbounded'), ('highs', 5893.705755757289, 'optimal'), ('scip', 5893.705756618257, 'optimal')]
