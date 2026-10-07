# MIC — Mixed-Integer & Constraint Solvers

MIC 是一组完全独立的 Rust 求解器 crate：**M**ixed-**I**nteger 精确求解（`milp/`）与 **C**onstraint 启发式求解（`planner/`）。两个 crate 由并行工作流同时开发并各自通过全部验证门禁，不共享任何代码、互不依赖，可单独使用、单独测试。

| 目录 | 求解器 | 核心算法 | 验证状态 |
|---|---|---|---|
| [`milp/`](milp/README.md) | 混合整数线性规划（MILP）**精确**求解器 | 单纯形法解 LP 松弛 + 分支定界 + Gomory 混合整数割 + 预处理/潜水启发式 | 5 个已知最优解算例全部精确命中；CLI 端到端一致 |
| [`planner/`](planner/README.md) | OptaPlanner 式**启发式**约束求解器 | HardSoftScore 硬软评分 + First/Best Fit 构造启发式 + 禁忌/模拟退火/迟到接受局部搜索 | 6 个算例的可行性（hard=0）、搜索改进、固定种子可复现全部满足 |

## 目录结构

```
MIC/
├── milp/               MILP 精确求解器（独立 crate）
│   ├── src/            建模/解析、单纯形、预处理、割平面、启发式、分支定界、编排
│   ├── instances/      7 个 .milp 算例（LP / 背包 / 指派 / 一般 MILP / 不可行 / 无界）
│   ├── tests/          验收测试（已知最优解断言）
│   └── README.md       算例格式、CLI、编程 API
├── planner/            OptaPlanner 式启发式求解器（独立 crate）
│   ├── src/            域模型、评分、RNG、解析、构造启发式、移动、局部搜索、编排
│   ├── instances/      6 个算例（排班 / 云分配 / 课程表 / 指派）
│   ├── tests/          验收测试 + 基准测试
│   └── README.md       算例格式、CLI、编程 API
├── docs/
│   └── 交付报告.md      验证结果、8 项评审发现（含证据与状态）、未覆盖范围
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

## 验证状态与已知问题

开发过程中的全部验证结果与独立代码评审结论见 [docs/交付报告.md](docs/交付报告.md)。要点：

- 所有验收门禁（编译、全部测试、已知最优/可行性/可复现断言、CLI 端到端）均通过；
- 评审发现的 **2 个高危问题已修复并通过回归验证**（milp 单纯形热启动误报 LimitReached；planner 解析器 CAP/RES 指令交错的静默建模错误）；
- 另有 **6 个中低危发现未修复**（数值语义、容差边界、文档陈旧等），每条都附证据与位置，处理优先级由您决定。
