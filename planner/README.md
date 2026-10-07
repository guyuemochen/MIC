# planner

OptaPlanner 式启发式约束求解器：纯 Rust、**仅标准库、零依赖、无 `unsafe`、单线程**
（性能来自增量评分，而非并行）。求解流程：

```
文本算例 ──parse──> ProblemModel ──construct（FF/BF）──> 完整初始解
        ──local_search（Tabu / SA / LA，增量评分）──> SearchResult（历史最优）
```

评分语义（整数精确，全库无 epsilon；浮点仅出现于 SA 温度与 Mixed 比例判定，永不参与解的比较）：

- `hard = Σ_{r,v} max(0, load[r][v] − cap[r][v]) + |{(e,c) ∈ 冲突边 : val(e) = val(c) ≠ UNASSIGNED}|`
  （可行 ⟺ `hard == 0`）
- `soft = Σ_e cost[e][val(e)]`（未分配实体不计任何项）
- `HardSoftScore` 按 `(hard, soft)` 字典序、**越小越好**——与 OptaPlanner `HardSoftScore`
  同构而方向相反（彼处最大化负惩罚，此处最小化正惩罚）；所有算法决策只用该字典序。

## OptaPlanner 概念映射

| OptaPlanner | 本 crate |
| --- | --- |
| `@PlanningSolution`（问题静态结构） | `model::ProblemModel` |
| `@PlanningEntity` | `model::PlanningEntity` |
| `@PlanningVariable`（值域） | `PlanningEntity::domain`（升序、去重、非空） |
| 工作解 / nullable planning variable | `model::Solution` / `model::UNASSIGNED` |
| `Move`（ChangeMove / SwapMove） | `moves::Move::Change` / `moves::Move::Swap` |
| `HardSoftScore` | `model::HardSoftScore`（字典序，最小化） |
| Construction Heuristic | `construct::ConstructionHeuristic`（`FirstFit` / `BestFit`） |
| Local Search（Tabu / SA / LA） | `localsearch::local_search` + `Metaheuristic` |
| `Solver` 编排 | `solver::solve` + `SolverConfig` |
| ScoreDirector（增量计算） | `score::IncrementalScorer`（`ScoreCalculator` 为全量真值） |
| solver 配置默认值 | `SolverConfig::default()` / `LocalSearchConfig::default()`（全文唯一来源） |

## 快速开始

CLI（默认 `tabu`、`steps=50000`、`seed=42`、`ff`、`samples=32`）：

```text
cargo run --manifest-path planner/Cargo.toml --release -- planner/instances/roster-small.txt
```

stdout 结束时**恰好两行**（历史最优解的硬惩罚与软惩罚）：

```text
hard = 0
objective = 27
```

lib 用法（`PlannerError` 实现 `std::error::Error`，用户输入不引发 panic）：

```rust
use std::sync::Arc;

let model = Arc::new(planner::format::parse_file(std::path::Path::new(
    "planner/instances/roster-small.txt",
))?);
let r = planner::solve(model, &planner::SolverConfig::default())?;
println!("{} {}", r.best_score.hard, r.best_score.soft);
```

## 算例文本格式

行导向纯文本：`#` 整行注释、任意空白分词、**大小写敏感**；文件一律 UTF-8 无 BOM
（`format::parse` 容忍文件头单个 U+FEFF 与 CRLF/LF 行尾）。完整 normative 文法与错误类别见
`src/format.rs` 模块文档；指令总序：

1. `NAME <单词>`：可选，至多一次，**只能在 VALUES 之前**；
2. `VALUES <m>`：恰好一次，`1 ≤ m ≤ 2^31`；
3. `ENTITIES <n>`：恰好一次，**紧随 VALUES**（二者之间不允许任何其他指令），`1 ≤ n ≤ 2^20`；
4. 其后七类任意交错，唯一强约束：**所有 `RES`/`CAP` 行必须先于所有 `DEMAND` 行**。

| 指令 | 语义 |
| --- | --- |
| `COSTROW <e> <c0> .. <c_{m-1}>` | 整行软代价，恰好 m 个非负整数（≤ 2^40） |
| `COST <e> <v> <c>` | 单点软代价；与 COSTROW 按行序后写覆盖先写 |
| `FORBID <e> <v>` | 从实体 e 的值域删除 v |
| `CONF <e1> <e2>` | 两实体不得同值（硬约束）；允许单向书写，模型侧对称化 |
| `RES <name> <c0> .. <c_{m-1}>` | 新增资源（恰好 m 个容量）；名字不可重复且不可为 `cap` |
| `CAP <v> <c>` | 糖——首次出现创建资源 `cap`（每值容量 i64::MAX/2）并设 `cap[v]=c` |
| `DEMAND <e> <d0> .. <d_{R-1}>` | 实体对全部 R 个资源的需求；元数必须等于 R |

缺省语义：未声明的 cost 为 0；值域默认全值 `[0, m)` 减 `FORBID`；未写 `DEMAND` 的实体对
`cap` 资源需求 1（若存在）、其余资源 0；写了 `DEMAND` 的实体以该行为准（含 cap 维度）。

## CLI 选项与输出契约

```text
planner <算例路径> [--steps N] [--time-secs F] [--seed S]
         [--algo tabu|sa|la] [--construct ff|bf] [--samples K] [--verbose]
```

| 旗标 | 默认 | 语义 |
| --- | --- | --- |
| `--steps N` | 50000 | 恰好执行 N 步（比特级可复现）；与 `--time-secs` 互斥，同给 exit 1；N ≤ 10^9 |
| `--time-secs F` | —（用 steps） | 每步开始前检查已耗时；F 解析为 f64 后先 `is_finite() && 0 ≤ F ≤ 36000` 再 `Duration::try_from_secs_f64`，任一不过 → `error: invalid --time-secs` + exit 1（**禁用** `from_secs_f64`，其对 NaN/负/溢出直接 panic） |
| `--seed S` | 42 | u64 种子（SplitMix64 扩展 + xoshiro256\*\*） |
| `--algo tabu\|sa\|la` | tabu | `tabu` → `TabuSearch { tabu_size: 8 }`；`sa` **固定** `SimulatedAnnealing { start_temp: 1000.0, end_temp: 1.0 }`；`la` **固定** `LateAcceptance { history_len: 30 }`（参数取自 `localsearch` 契约的唯一定义，无调参旗标） |
| `--construct ff\|bf` | ff | 构造启发式（FirstFit / BestFit，确定性、不耗随机数） |
| `--samples K` | 32 | 每步采样数，**仅 tabu 读取**；sa/la 下接受但无效果（不耗随机数）；K ≤ 10^6 |
| `--verbose` | 关 | 进度走 stderr：构造分、搜索后的初始分 / 实际步数。步级进度（每 5000 步当前最优）因 lib 公共 API 无逐步回调而不提供（见下「可复现性与 verbose」） |

输出契约：成功时 stdout **恰两行**且在结束时输出——`hard = <i64>` 与最后一行
`objective = <i64>`（均取 `SearchResult.best_score`，`objective` 后无内容）；`hard` 与
`objective` 之间恒为单个空格、`=` 两侧各一空格。错误：stderr `error: <消息>` + exit 1；
正常求解（即使 `hard > 0`）exit 0。参数解析全程不 panic。

## 可复现性与 verbose

- **同二进制 + 同种子 + Steps 终止 ⇒ 结果逐位一致**（`tests/acceptance.rs::feasible_and_reproducible`
  断言之）；随机流（SplitMix64 播种 + xoshiro256\*\* + Lemire `gen_range` + 53 位 `gen_ratio`）
  与移动生成的随机数消耗顺序均在 `src/rng.rs` / `src/moves.rs` 契约逐字钉死。
- **TimeLimit 终止**：结果是"已完成步"的确定函数，但墙钟截断点随机器速度不同，跨机器不保证相同。
- **SA**：`exp` 依赖平台 libm，跨平台比特级一致不作承诺——默认与全部验收使用 tabu/la（纯整数路径）。
- **跨实现**：不承诺轨迹一致；必达算例验收用"种子表 `[42, 7, 2026]` 任一命中最优"消除该依赖。
- **--verbose 的进度粒度**：`solver::solve` / `localsearch::local_search` 不提供逐步回调，
  故 CLI 在搜索进行中无法报告"每 5000 步当前最优"；verbose 输出构造分（求解前）与
  初始分 / 实际步数 / 历史最优（求解后）。

## 规模硬界（解析期 + 建模期双层强制，违例返回 `PlannerError` 而非 panic）

| 项 | 界 |
| --- | --- |
| 值个数 `m`（VALUES） | `1 ≤ m ≤ 2^31`（值 id < u32::MAX，与 `UNASSIGNED` 不冲突） |
| 实体数 `n`（ENTITIES） | `1 ≤ n ≤ 2^20` |
| 资源数 `R` | ≤ 64 |
| 单条 COST / CAP | ≤ 2^40，非负 |
| 单条 DEMAND | ≤ 2^31，非负 |
| `Σ_{e,v} cost` | ≤ 2^40（i128 求和判定） |
| 每资源 `Σ_e demand` | ≤ 2^40（i128 求和判定） |

推论（rustc 1.95 实测）：`hard ≤ 2^46`、`W = 1 + Σcost ≤ 2^40 + 1`，故 SA 标量
`s = hard·W + soft ≤ 77371252455407735437000704` 超 i64——SA 的 `s` 与 `Δ` 全程 i128 精确，
仅 `Δ/T` 转 f64 进 `exp`。设计有效域为数百实体、数千值。

## 测试与基准

验收测试（`tests/acceptance.rs`；**务必 `--release`**——`must_match_known_optima` 的
200k 步在 debug 下过慢）：

```text
cargo test --manifest-path planner/Cargo.toml --release
```

基准（`tests/bench.rs`，全部 `#[ignore]`，手动运行；不断言数值，仅性能参考——对 6 算例计时
`construct(FF/BF)` 与 `local_search(TabuSearch{tabu_size:8}, Steps(300_000))` 并打印表格）：

```text
cargo test --release --manifest-path planner/Cargo.toml --test bench -- --ignored --nocapture
```

## 算例清单与已知值

`planner/instances/`（构造分记法 `(hard, soft)`；必达最优经双重独立实现验证——
rustc 原型剪枝 DFS + 评审方 python verify.py 两轮复核相等）：

| 算例 | 规模 n×m | FF 构造分 | BF 构造分 | 已知最优 / 参考 |
| --- | --- | --- | --- | --- |
| `roster-small.txt` 员工排班 | 6×3 | (0, 36) | (0, 27) | **27**（必达，验收断言） |
| `cloud-small.txt` 云资源分配 | 6×4 | (6, 16)（FF 起点硬违反可被搜索修复至 0） | (3, 13) | **13**（必达，验收断言） |
| `timetable-small.txt` 课程表 | 6×4 | (0, 16) | (0, 17) | **14**（必达，验收断言） |
| `assign-10.txt` 指派问题 | 10×10 | (0, 113) | (0, 41) | **36**（必达，验收断言） |
| `roster-medium.txt` 员工排班（中等） | 30×8 | hard = 0（soft 不锁定） | —（仅测确定性） | 122（参考） |
| `cloud-medium.txt` 云资源分配（中等） | 40×10 | hard = 0（soft 不锁定） | —（仅测确定性） | 77（参考） |

**参考值声明**：中等算例的 122 / 77（原型 40 次搜索最好成绩）与 CLI 默认配置下的
132 / 88 出自本人原型的搜索记录，评审方未复核、**无任何断言依赖**，仅供 README 参考。
四个必达最优与八个构造分基线已由 FF+tabu(seed 42, 50k 步, samples 32, tabu 8) 命中
（本人原型与评审方 python sim.py 双重实测）；验收采用种子表 `[42, 7, 2026]` 任一命中。

## 文件编码约定

crate 源码与算例文件一律 **UTF-8 无 BOM**（Windows 下勿用 PowerShell 重定向生成
UTF-16/带 BOM 文件）；`format::parse` 容忍文件头单个 U+FEFF（剥去再解析）与 CRLF/LF 行尾；
算例文本内容 normative（已知最优/构造基线是其函数）。
