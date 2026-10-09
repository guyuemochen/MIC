//! # planner
//!
//! 纯标准库、零依赖、无 unsafe 的 OptaPlanner 式启发式约束求解 crate
//! （threads=1 时单线程求解路径保持不变；[`parallel`] 提供多起点并行编排）。
//!
//! ## 模块总览
//! - [`rng`]：确定性伪随机数（SplitMix64 播种 + xoshiro256\*\*，随机流逐位钉死；
//!   多起点种子派生 [`rng::derive_seed`]）
//! - [`model`]：域模型（ProblemModel / PlanningEntity / 值域 / Solution / HardSoftScore）
//! - [`score`]：全量评分器（真值）+ 增量评分器（preview/apply 共用同一增量内核）
//! - [`format`]：行导向文本算例解析（指令总序与数值界在解析期强制）
//! - [`construct`]：构造启发式（First Fit / Best Fit，返回 Result 不 panic）
//! - [`moves`]：Change/Swap 移动与随机邻域（随机数消耗顺序逐字固定、跨类回退为单一无条件规则）
//! - [`localsearch`]：禁忌搜索 / 模拟退火 / 迟到接受 + 步数/时间终止
//! - [`solver`]：构造 → 局部搜索编排
//! - [`parallel`]：多起点并行编排（std scoped threads、零共享可变状态、聚合确定性）
//!
//! ## 求解流程与数据流
//! 解析并校验 → 构造完整初始解（未分配变量记 [`model::UNASSIGNED`]）→ 元启发式搜索
//! （每步增量预演，字典序 `(hard, soft)` 比较，跟踪历史最优快照）→ [`localsearch::SearchResult`]。
//! 核心是"单一变更路径"：解只能经 [`moves::Move::apply`] 同步修改 Solution 与
//! [`score::IncrementalScorer`] 影子状态，增量与全量重算一致。
//!
//! ## CLI 输出契约（见 src/main.rs）
//! stdout 结束时恰两行：`hard = <整数>`（硬惩罚，0 = 可行）与最后一行 `objective = <整数>`
//! （软惩罚，越小越好）。
//!
//! ## 实现状态
//! 全模块均已完整实现（非骨架）：解析与建模、FF/BF 构造、禁忌/模拟退火/迟到接受局部
//! 搜索、单线程编排与多起点并行（[`parallel`]）全部可用；随机流与种子派生契约在
//! [`rng`] / [`crate::moves`] 逐位钉死，验收测试见 `tests/acceptance.rs`
//! （单线程六项 + 并行五项）。

#![forbid(unsafe_code)]

pub mod rng;
pub mod model;
pub mod score;
pub mod format;
pub mod construct;
pub mod moves;
pub mod localsearch;
pub mod solver;
pub mod parallel;

pub use model::{
    EntityId, ValueId, UNASSIGNED, HardSoftScore, ProblemModel, PlanningEntity, Resource,
    Solution, PlannerError,
};
pub use score::{ScoreCalculator, IncrementalScorer};
pub use format::{parse, parse_file};
pub use construct::{ConstructionHeuristic, construct};
pub use moves::{Move, MoveKind, random_move};
pub use localsearch::{Termination, Metaheuristic, LocalSearchConfig, SearchResult, local_search};
pub use solver::{SolverConfig, solve};
pub use parallel::{solve_parallel, solve_parallel_detailed, ParallelReport};
