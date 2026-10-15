# MIC — Mixed-Integer & Constraint Solvers

MIC 是一组完全独立的 Rust 求解器 crate：**M**ixed-**I**nteger 精确求解（`milp/`）与 **C**onstraint 启发式求解（`planner/`）。两个 crate 由并行工作流同时开发并各自通过全部验证门禁，不共享任何代码、互不依赖，可单独使用、单独测试。

| 目录 | 求解器 | 核心算法 | 验证状态 |
|---|---|---|---|
| [`milp/`](milp/README.md) | 混合整数线性规划（MILP）**精确**求解器 | 单纯形法解 LP 松弛**双引擎**（稠密表格 / 稀疏修正单纯形，`--lp-engine auto` 类自适应默认：表格 ≥10^6 格走稀疏引擎；另有大表稀疏簿记、退化停滞的 Bland 释放）+ 分支定界（类自适应分支）+ Gomory/聚合 MIR 割（循环预算回滚、树内局部割）+ 数值缩放 + 预处理/潜水/子集和 MITM 启发式；多线程（工作栈按类选栈/本地亲和/全局 DFS，2026-10 W1 起多线程默认按类选栈）、Gurobi 式进度日志、自适应启发式与策略级自动降级、t1 预算制确定性模式（`--determinism budget`） | 6 个已知最优解算例全部精确命中；CLI 端到端一致；大规模基准（bench gate/full 计时矩阵，18 算例 t1/t4）九轮逐轮 A/B 验证（八轮保留 + 一轮空集回退）并经两轮独立复核（无正确性阻断项）；稀疏核心与批量接口轮 S1–S5 全保留（lp-sparse-4000x800 默认引擎约 3.5×、目标恒同），micpython 与 CBC/HiGHS/SCIP 外部交叉验证 14 算例全对账；5000 变量背包与独立 DP 真值一致 |
| [`planner/`](planner/README.md) | OptaPlanner 式**启发式**约束求解器 | HardSoftScore 硬软评分 + First/Best Fit 构造启发式 + 禁忌/模拟退火/迟到接受局部搜索 | 6 个算例的可行性（hard=0）、搜索改进、固定种子可复现全部满足 |

## 目录结构

```
MIC/
├── milp/               MILP 精确求解器（独立 crate）
│   ├── src/            建模/解析、单纯形双引擎（稠密表格 + 稀疏修正）、预处理、割平面、数值缩放、启发式、分支定界、编排
│   ├── instances/      8 个 .milp 算例（LP / 背包 / 指派 / 一般 MILP / 不可行 / 无界）
│   ├── tests/          验收测试（已知最优解断言）
│   └── README.md       算例格式、CLI、编程 API
├── planner/            OptaPlanner 式启发式求解器（独立 crate）
│   ├── src/            域模型、评分、RNG、解析、构造启发式、移动、局部搜索、编排
│   ├── instances/      6 个算例（排班 / 云分配 / 课程表 / 指派）
│   ├── tests/          验收测试 + 基准测试
│   └── README.md       算例格式、CLI、编程 API
├── milp-ffi/           milp 的 C ABI 绑定 crate（cdylib，std-only）
├── micpython/          gurobipy 风格 Python 包（ctypes，零第三方依赖）
│   ├── examples/       knapsack.py 端到端示例
│   ├── tests/          Python unittest 套件（已知最优/状态/参数/句柄）
│   └── README.md       安装前置与完整 API 参考
├── bench/              大规模基准算例与生成器（gen_instances.py，种子确定；计时矩阵
│   └── instances/      gate=6 / full=18 算例：背包 200–5000（附 DP 核对最优 .opt）、多背包、
│                       指派、混合 MILP、宽 LP、稀疏 LP 4000×800、setcover 400×800、marketsplit，共 22 个 .milp）
├── docs/
│   ├── 交付报告.md      验证结果、评审发现（含证据与状态）、未覆盖范围
│   └── large-scale-report.md   milp 大规模优化轮（R1–R8 保留 + R9 回退 + P1–P3 + W1–W7 + S1–S5 稀疏核心与批量接口）逐轮记录、实测数据与遗留问题
└── README.md           本文件
```

## 快速开始

两个 crate 均为纯标准库实现（无外部依赖），从工作区根目录运行：

### milp（精确求解，保证最优）

```bash
cargo run --manifest-path milp/Cargo.toml --release -- milp/instances/knapsack-10.milp
```

常用参数：`--verbose`（输出求解过程）、`--time-limit-secs N`（限时）。
求得最优时最后一行输出 `objective = <数值>`；不可行/无界有独立状态与退出码。

大规模算例示例（`bench/instances/`，2026-10 大规模优化轮 R1–R8 与 S1–S5 验证）——
稀疏纯 LP 800 行 × 4000 列（密度 8%），默认引擎 auto（≥10^6 格自动走稀疏修正
单纯形）单线程约 9–10s 证最优（tableau 稠密引擎约 32–36s，约 3.5×）：

```bash
cargo run --manifest-path milp/Cargo.toml --release -- bench/instances/lp-sparse-4000x800.milp
# status = optimal / objective = 5893.705757（1 节点，7774 次单纯形迭代，两引擎迭代序列一致）
```

集合覆盖 100×150（R8 割循环退化修复后 gate t1 约 0.24s，树 45→13 节点）：

```bash
cargo run --manifest-path milp/Cargo.toml --release -- bench/instances/setcover-100x150.milp
# status = optimal / objective = 276.000000
```

（更快的核验路径：`bench/instances/knapsack-5000.milp` 单线程约 2.2s 证最优
88526，与生成器独立 DP 真值一致。）

### planner（启发式求解，可行 + 快速）

```bash
cargo run --manifest-path planner/Cargo.toml --release -- planner/instances/roster-small.txt \
  --steps 2000 --seed 7 --algo tabu --construct bf
```

常用参数：`--steps N` 或 `--time-secs F`（二选一）、`--seed N`（固定种子结果可复现）、
`--algo tabu|sa|la`（禁忌/模拟退火/迟到接受）、`--construct ff|bf`（First/Best Fit）、
`--samples K`（每步邻域采样数）、`--verbose`。
结束时输出 `hard = <整数>`（0 表示可行）与 `objective = <数值>`（软惩罚，越小越好）。

## 运行测试

```bash
# 单元 + 集成测试
cargo test --manifest-path milp/Cargo.toml
cargo test --manifest-path planner/Cargo.toml

# 验收测试（release）：milp 断言已知最优；planner 断言可行性/改进/可复现/必达算例
cargo test --manifest-path milp/Cargo.toml --release --test acceptance
cargo test --manifest-path planner/Cargo.toml --release --test acceptance
```

## API 文档

```bash
cargo doc --manifest-path milp/Cargo.toml --no-deps --open
cargo doc --manifest-path planner/Cargo.toml --no-deps --open
```

编程 API 的详细说明见各自的 `src/lib.rs` 文档注释与 README。

## 从 Python 调用（gurobipy 风格）

```bash
cargo build --release --manifest-path milp-ffi/Cargo.toml   # 先产出 DLL
python micpython/examples/knapsack.py                        # 端到端示例
python -m unittest discover -s micpython/tests -t .          # Python 测试套件
```

```python
from micpython import Model, GRB
m = Model("knapsack")
x = m.addVars(range(10), vtype=GRB.BINARY)
m.setObjective(quicksum(v[i] * x[i] for i in range(10)), GRB.MAXIMIZE)
m.addConstr(quicksum(w[i] * x[i] for i in range(10)) <= capacity)
m.setParam("Threads", 4)
m.optimize()
if m.Status == GRB.OPTIMAL:
    print(m.ObjVal)
```

完整 API（TimeLimit/MIPGap/HeuristicMode 参数、状态映射、句柄管理）见 [micpython/README.md](micpython/README.md)。

大规模建模性能（S1）：`quicksum` 为共享项字典的就地累加（O(n)；旧 O(n²)
算子链在 5000 项时约 1.3s），`addVars` 单次 FFI 跨界建整批变量（任一元素
非法整批拒绝、零突变）——5000 变量模型构建从约 2.35s 降到约 0.05s。

与开源求解器交叉验证（S5，conda 环境 `mic-bench`）：

```bash
conda run --no-capture-output -n mic-bench python -u bench/python_bench.py
# 单算例核验（仅打印不落盘）：python -u bench/python_bench.py --pair <算例名> micpython
```

14 算例 micpython 全部 optimal，与 HiGHS/SCIP 全部一致（纯 LP 位级吻合）。

## 验证状态与已知问题

开发过程中的全部验证结果与独立代码评审结论见 [docs/交付报告.md](docs/交付报告.md)。要点：

- 所有验收门禁（编译、全部测试、已知最优/可行性/可复现断言、CLI 端到端）均通过；
- **milp 大规模基准**（2026-10，分支 `feat/gurobi-large-scale`）：九轮优化
  （八轮保留 + 一轮空集回退）——稀疏混合表格（R1）、聚合 MIR 割（R2）、
  分支类自适应默认（R3）、数值缩放（R4）、子集和 MITM 启发式（R5）、全局
  DFS 工作栈（R6）、复核遗留修复（R7）、setcover 割循环退化修复（R8：
  setcover-100x150 树 45→13 节点、迭代 −72%）、并行证据轮空集回退（R9）
  ——全部以同窗交错 A/B 中位数为口径实测，`threads=1` 行为字段
  （status/nodes/迭代/目标）对账或逐位一致、目标值与已知最优/DP 真值全部
  精确；bench 计时矩阵扩为 gate（6 算例、约 3 分钟）/full（18 算例 ×
  t1/t4）两集（gate 单次约 2–3 分钟）；两轮独立复核均通过（无正确性阻断
  项，遗留发现见报告清单）。逐轮数据（含各轮内部被回退的尝试）与遗留问题
  见 [docs/large-scale-report.md](docs/large-scale-report.md)；
- **milp 稀疏核心与批量接口轮**（2026-10-14/15，S1–S5 五轮全保留）：S1
  批量 addVars FFI 符号（22→23，ABI 恒 1）+ micpython quicksum O(n) 修复
  （5000 变量建模 2.35→0.05s）；S2 稀疏修正单纯形核心（CSC + 稀疏 LU +
  ETA + 有界原单纯形，lp-sparse-4000x800 约 3.5×、迭代序列与稠密引擎逐次
  一致）；S3 接入主求解路径（LpEngine trait、热启动/加割、`--lp-engine`
  旗标，14 bench 实例双引擎目标恒同对账）；S4 类自适应默认 auto（≥10^6 格
  走 revised，阈下留 tableau，acceptance 默认输出逐字节不变）；S5 数值稳健
  性常量对齐 + python_bench 全量外部复跑（micpython 14/14 optimal、与
  HiGHS/SCIP 全部一致、MISMATCH 零新增）。终局独立复核通过（无正确性阻断
  项，5 条发现均为文档/死代码类，见报告 §8.7）；
- 评审发现的 **2 个高危问题已修复并通过回归验证**（milp 单纯形热启动误报 LimitReached；planner 解析器 CAP/RES 指令交错的静默建模错误）；
- 另有 **6 个中低危发现未修复**（数值语义、容差边界、文档陈旧等），每条都附证据与位置，处理优先级由您决定。
