//! milp —— 纯 Rust 标准库实现的混合整数线性规划（MILP）精确求解器。
//!
//! 模块划分（工作包归属见设计文档）：
//! - [`error`]   ：统一错误类型 `MilpError` 与 `Result` 别名（WP1）；
//! - [`model`]   ：问题建模（变量/约束/目标、容差集中定义）（WP1）；
//! - [`instance`]：.milp 人类可读算例格式的解析与序列化（WP1）；
//! - [`simplex`] ：有界变量两阶段表格单纯形 LP 求解器（WP2）；
//! - [`presolve`]：保守预处理与还原（WP3）；
//! - [`cuts`]    ：根节点 Gomory 混合整数割 + 覆盖割（WP3）；
//! - [`heuristics`]：取整/潜水启发式（WP3）；
//! - [`branch_bound`]：深度优先分支定界与受控搜索入口（WP4）；
//! - [`solver`]  ：统一求解编排与终检（WP4）；
//! - [`search_ctrl`]：搜索控制层（进度行/gap/自适应启发式/自动降级，纯逻辑）；
//! - [`parallel_bb`]：并行分支定界引擎（threads ≥ 2 时由 branch_bound 委托）；
//! - CLI 入口在 `src/main.rs`，验收测试在 `tests/acceptance.rs`，
//!   热启动差分在 `tests/warm_cold_parity.rs`（WP5）。
//!
//! 约定：纯标准库、无 unsafe、单线程确定性（并列一律按最小下标破键）；
//! 用户输入路径禁止 unwrap/expect/panic；诚实性铁律——只有满足最优性终止
//! 条件才报 Optimal，触发上限一律报 LimitReached。
//!
//! 当前进度：全部模块已交付实现（WP1–WP5），含并行引擎与搜索控制层；
//! 行为由 crate 单元测试与验收/差分测试钉死。

// 注：模块声明顺序为 WP1 契约逐字内容（非字母序），勿被 rustfmt
// （reorder_modules 默认开启）重排。
pub mod error;
pub mod model;
pub mod instance;
pub mod simplex;
pub mod presolve;
pub mod cuts;
pub mod heuristics;
pub mod branch_bound;
pub mod solver;
// 追加声明（并行分支定界引擎与搜索控制层）：仅追加，勿重排上方既有声明。
pub mod search_ctrl;
pub mod parallel_bb;
