//! 稀疏修正单纯形 LP 求解器（S2 standalone 核心；S3 经 [`LpEngine`]
//! trait 接入主求解路径；S4 起为**类自适应默认**的 revised 侧——
//! `LpEngineKind::Auto`（默认）按 [`effective_engine`] 规则
//! m×(n+m) ≥ 10^6 格 ⇒ 本引擎，否则稠密 tableau；显式指定覆盖优先）。
//!
//! # 与 [`crate::simplex`]（稠密表格引擎）的语义对齐契约
//! - **同一 LP 语义**：内部统一 min（Maximize 取负，报告还原 + obj_offset）、
//!   行按 max|a_ij| 归一（`row_scale` 对偶还原）、行永不翻转、统一残差规则
//!   初始基（Le 松弛 / Ge 剩余 / Eq 人工；正右端 Ge 与违反 Le 同样走人工列）、
//!   近零行（max|a_ij| ≤ tol.zero）不入表按原尺度 rhs 判定；
//! - **同一两阶段纪律**：Phase-1 = min Σ人工；最优值 > tol.feasibility ⇒
//!   Infeasible；人工列**从不允许入基**（定价只扫静态列段，比『禁再入基』
//!   更强且同安全）；Phase-1 后残基人工变量值清零并尽力驱赶出基（t=0 主元，
//!   最小下标非人工候选），整行无候选的人工变量留基值恒 0；
//! - **同一比值测试语义**（F7 全部点名项）：基变量上下界阻挡 + 入基变量
//!   自身界翻转；主元 |p| ≥ tol.pivot 才允许出基，方向阻挡但主元过小 ⇒
//!   FilteredBlock（数值失败）；无任何阻挡（含翻转不可用）⇒ NoBlock ⇒
//!   Phase-2 报 Unbounded；并列带宽 `band = min(feas·(1+|t_min|), feas·
//!   RATIO_BAND_CAP)`（simplex.rs:116-125 实测踩坑的封顶，本引擎同款）；
//! - **同一反循环/上限纪律**：Dantzig 得分（贴下界 −rc / 贴上界 +rc）最大、
//!   并列最小下标；连续退化 > 50 或总迭代 > 3×(m+n) 且停滞 ⇒ Bland（最小
//!   下标），非退化步 + 驻留 8 步释放（每次调用 ≤64 次重入；总迭代触发逐
//!   迭代重申）；迭代上限 max(10_000, 20×(m+n))（`set_iter_limit_hint`
//!   可收紧，R8 割循环接线点，F4）；数值失败（含主元过小、NaN/Inf）诚实
//!   映射 IterationLimit + objective=NaN，绝不冒充 Optimal；W2 内部墙钟
//!   （`set_time_deadline`）条件读钟——未注入 ⇒ 零墙钟读点；
//! - **同一解构造**：x 提取（基本值 10·feas 内夹回界）、Optimal 终检
//!   （基变量界违反/残基人工变量非零 ⇒ 降级 IterationLimit）、目标按 x
//!   重算、对偶 = c_B·B⁻¹ × row_scale。
//!
//! # 修正单纯形机制（与表格引擎的实现差异；目标值恒同是红线）
//! - 约束矩阵 CSC（按列，行索引升序）；基矩阵因子化 = **Gilbert-Peierls
//!   左视稀疏 LU**（逐列稀疏三角求解 + 非零模式 DFS），**阈值主元**：候选
//!   |v| ≥ 0.1·max|v| 方可入选，再按『L 行密度最小』破键（Markowitz 的
//!   左视代理——完整 (r−1)(c−1) 计数需右视活动子矩阵，v1 以阈值稳定性 +
//!   填充代理近似，右视全计数留后续），并列取最小行号；未分配候选
//!   max|v| < tol.pivot 视为奇异 ⇒ 数值失败；
//! - FTRAN（Bα = a）/ BTRAN（Bᵀy = c）稀疏三角求解 + **乘子形式 ETA 文件**：
//!   每主元追加 eta（r, α），FTRAN 顺序施加 E⁻¹（旧→新）、BTRAN 逆序施加
//!   (Eᵀ)⁻¹（新→旧）后走基础 LU；
//! - **refactorization 纪律**（类比稠密引擎每 200 主元重算 B⁻¹b）：每
//!   [`REFACTOR_INTERVAL`] 主元整基重因子化 + 从原始数据重算 x_B =
//!   B⁻¹(b̂ − N·x_N)（漂移双保险）；数值异常（主元过小 / NaN / 检验数非
//!   有限）时若非刚重因子化则先重因子化重算再重试本迭代（确定性的
//!   just_refactored 闸保证至多重试一次，无墙钟决策），仍失败才上报；
//! - x_B 维护与表格引擎同款口径：主元步 x_B[k] −= δt·α_k、x_B[r] =
//!   入基旧界值 + δt；界翻转 x_B[k] −= Δ·α_k；
//! - 定价：每迭代 BTRAN(c_B) + 全列 Dantzig 扫描（rc_j = c_j − y·A_j 按
//!   CSC 列现算，无增量漂移——与表格引擎『三入口批量重算』同款新鲜度
//!   契约，本引擎每迭代重算即天然新鲜）；v1 不做部分定价：lp-sparse 实测
//!   全列扫描（nnz≈2.6e5 ≈ 0.1ms/迭代）在 9.5s 总墙钟、1.22ms/迭代成本
//!   中非主要成分，且本引擎迭代序列与稠密引擎逐次数一致（7774 次）——
//!   部分定价会改换入基序列、引入迭代数回归风险，留后续按证据启用。
//!
//! # 范围边界（诚实声明）
//! - S3 起提供全套引擎面：`solve` 冷启动 / `resolve` 热启动（snap+传播 →
//!   对偶修复快路径或整基重建，语义对齐 simplex.rs::resolve）/ `add_row`
//!   割行（原尺度不缩放；脏因子化标记 + 下次调用整基重构——F6『宜随后
//!   重构』的最简落法）/ W6 快照（Clone）。访问器按 cuts.rs 消费面
//!   （cuts.rs:267/448-465/1247/1754）：表格行 = e_iᵀB⁻¹A 经 BTRAN +
//!   行收集实现，列空间 [结构|松弛]、人工列不外露（simplex.rs:262-269
//!   同款约定）；
//! - **revised 仅接线顺序路径**（threads=1）：并行引擎直接构造稠密
//!   Simplex，revised × threads≥2 在 branch_and_bound_controlled 与 CLI
//!   解析层显式报错（F5）；
//! - 并列破键全部确定性（最小下标 / 最小行号 / 升序值计算单遍），无哈希
//!   迭代序、无墙钟决策（deadline 未注入零读钟）——threads=1 / t1 双跑
//!   确定性契约。
//!
//! # 数值注记
//! 行归一化（max|a| = 1）+ 阈值主元 0.1 给出与稠密引擎可比的稳定性；
//! 检验数/表格列由三角求解现算（非增量维护），误差随因子化精度而非
//! 主元序列累积；refactorization 周期重置一切累积误差。

use crate::error::{MilpError, Result};
use crate::model::{Problem, RelOp, Sense, Tolerances};
use crate::simplex::{
    Bounds, ColumnLayout, LpSolution, LpStatus, NonbasicCol, RowOrigin, Simplex, SlackCol,
};

// ===========================================================================
// LP 引擎抽象（S3：B&B/cuts/heuristics 经 dyn 引擎驱动，tableau 默认不变）
// ===========================================================================

/// LP 引擎选择（CLI `--lp-engine` / 库 `BbOptions::lp_engine` 两层）。
///
/// 默认 [`LpEngineKind::Auto`]（S4 类自适应定档）：表格规模
/// m×(n+m) ≥ [`ENGINE_AUTO_CELL_THRESHOLD`]（10^6 格）⇒ [`LpEngineKind::Revised`]
///（稀疏修正单纯形），否则 [`LpEngineKind::Tableau`]（稠密表格）。显式
/// `Tableau`/`Revised` 覆盖优先。`Revised` **仅接线顺序路径**（threads=1）：
/// 并行引擎（parallel_bb）直接构造稠密 Simplex，revised × threads≥2 在
/// `branch_and_bound_controlled` 入口显式报错（F5 冷审：静默忽略旗标会
/// 形成文档外行为分叉）。
///
/// 定档证据（S4 同窗 A/B，见轮次汇报）：revised 决定性胜出的是大规模
/// **纯 LP 松弛主导**类——lp-sparse-4000x800（800×4800=3.84M 格）t1
/// ~10.4s vs tableau ~35s（3.4×）、lp-wide-1000x2000（2000×3000=6M 格）
/// 8.5s vs 15.1s（1.78×）；阈下回归或噪声级——setcover-100x150 revised
/// 慢 2.3×、mknap-5x300（7.5k 格）revised 慢 3.6×、setcover-300x600
///（270k 格）revised 1.43×（单窗噪声内，保守留 tableau）、knapsack/
/// assign 小算例毫秒级差异。阈下小算例保持 tableau 兼防 codegen 扰动
///（R1 教训：+1-2% 不可复现噪声）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LpEngineKind {
    /// 类自适应（默认）：按 [`effective_engine`] 的规模规则解析。
    Auto,
    /// 稠密表格两阶段单纯形（simplex.rs）。
    Tableau,
    /// 稀疏修正单纯形（本模块；CSC + 稀疏 LU + ETA + 对偶热修复）。
    Revised,
}

impl Default for LpEngineKind {
    fn default() -> Self {
        LpEngineKind::Auto
    }
}

/// 类自适应规模阈（S4 定档）：m×(n+m) ≥ 该值 ⇒ revised。取 10^6——
/// 落在 setcover-400x800（480k）与 lp-sparse（3.84M）之间的整数量级空带，
/// 规则简单可解释（『表格超百万格』）。
pub const ENGINE_AUTO_CELL_THRESHOLD: u64 = 1_000_000;

/// 解析有效引擎（S4 类自适应）：`Auto` 按表格规模定档，显式指定原样
/// 返回（覆盖优先）。m = 约束行数、n = 变量数、ncols ≈ n + m（松弛列）。
pub fn effective_engine(problem: &Problem, kind: LpEngineKind) -> LpEngineKind {
    match kind {
        LpEngineKind::Auto => {
            let m = problem.constraints.len() as u64;
            let n = problem.vars.len() as u64;
            if m.saturating_mul(n + m) >= ENGINE_AUTO_CELL_THRESHOLD {
                LpEngineKind::Revised
            } else {
                LpEngineKind::Tableau
            }
        }
        explicit => explicit,
    }
}

/// LP 引擎统一接口：branch_bound / cuts / heuristics 的全部消费面。
///
/// 对 [`Simplex`]（稠密）与 [`SparseSimplex`]（本模块）各有一份实现；
/// B&B 持有 `Box<dyn LpEngine>`（W6 局部割的子树作用域快照经
/// [`LpEngine::snapshot`]——Clone 的对象安全等价物）。
pub trait LpEngine {
    /// 两阶段冷启动求解（忽略此前内部状态）。
    fn solve(&mut self, bounds: &Bounds) -> LpSolution;
    /// 热启动求解（界变化 snap+传播 → 对偶修复快路径或整基重建 → 收尾；
    /// 语义对齐 simplex.rs::resolve，见各实现的模块文档）。
    fn resolve(&mut self, bounds: &Bounds) -> LpSolution;
    /// 追加割行（terms 结构变量空间原尺度、不缩放；返回割序号 k）。
    fn add_row(&mut self, terms: &[(usize, f64)], op: RelOp, rhs: f64) -> Result<usize>;
    /// 内部墙钟时限（W2 条件读钟；None = 零读钟）。
    fn set_time_deadline(&mut self, deadline: Option<std::time::Instant>);
    /// 当前内部墙钟时限（R8 回滚透传用）。
    fn time_deadline(&self) -> Option<std::time::Instant>;
    /// 单次调用迭代上限提示（R8 割循环收紧用，F4）。
    fn set_iter_limit_hint(&mut self, hint: Option<u64>);
    /// 累计迭代数。
    fn iterations(&self) -> u64;
    /// 各行基列下标（结构+松弛列空间）。
    fn basis_columns(&self) -> &[usize];
    /// 当前全部非基列及其所贴界（不含人工列）。
    fn nonbasic_columns(&self) -> Vec<NonbasicCol>;
    /// 列布局（结构列/松弛列含义；人工列不外露）。
    fn column_layout(&self) -> &ColumnLayout;
    /// 原始约束下标 → 表格行；近零行 None。
    fn row_of_constraint(&self, constraint_index: usize) -> Option<usize>;
    /// 割序号 → 表格行。
    fn row_of_cut(&self, cut_index: usize) -> Option<usize>;
    /// 表格行数（含割行）。
    fn num_rows(&self) -> usize;
    /// 表格行稀疏读取：(非零 (列, e_rowᵀB⁻¹A_j) 项[结构+松弛段、列升序],
    /// 该行基变量当前值)。
    fn tableau_row_terms(&self, row: usize) -> (Vec<(usize, f64)>, f64);
    /// 深拷贝快照（W6 子树作用域协议；Clone 的对象安全等价物）。
    fn snapshot(&self) -> Box<dyn LpEngine>;
}

/// 按引擎类型构造（branch_bound / cuts 回滚重建的统一工厂）。
pub fn new_engine(
    problem: &Problem,
    tol: Tolerances,
    kind: LpEngineKind,
) -> Result<Box<dyn LpEngine>> {
    match effective_engine(problem, kind) {
        LpEngineKind::Tableau => Ok(Box::new(Simplex::new(problem, tol)?)),
        LpEngineKind::Revised => Ok(Box::new(SparseSimplex::new(problem, tol)?)),
        LpEngineKind::Auto => unreachable!("effective_engine 已解析 Auto"),
    }
}

/// Box 透传实现（`&Box<dyn LpEngine>` → `&dyn LpEngine` 的去尺寸胁迫与
/// 泛型 `E = Box<dyn LpEngine>` 的整箱替换 `*s = rebuilt` 都需要它）。
impl<T: LpEngine + ?Sized> LpEngine for Box<T> {
    fn solve(&mut self, bounds: &Bounds) -> LpSolution {
        (**self).solve(bounds)
    }
    fn resolve(&mut self, bounds: &Bounds) -> LpSolution {
        (**self).resolve(bounds)
    }
    fn add_row(&mut self, terms: &[(usize, f64)], op: RelOp, rhs: f64) -> Result<usize> {
        (**self).add_row(terms, op, rhs)
    }
    fn set_time_deadline(&mut self, deadline: Option<std::time::Instant>) {
        (**self).set_time_deadline(deadline)
    }
    fn time_deadline(&self) -> Option<std::time::Instant> {
        (**self).time_deadline()
    }
    fn set_iter_limit_hint(&mut self, hint: Option<u64>) {
        (**self).set_iter_limit_hint(hint)
    }
    fn iterations(&self) -> u64 {
        (**self).iterations()
    }
    fn basis_columns(&self) -> &[usize] {
        (**self).basis_columns()
    }
    fn nonbasic_columns(&self) -> Vec<NonbasicCol> {
        (**self).nonbasic_columns()
    }
    fn column_layout(&self) -> &ColumnLayout {
        (**self).column_layout()
    }
    fn row_of_constraint(&self, constraint_index: usize) -> Option<usize> {
        (**self).row_of_constraint(constraint_index)
    }
    fn row_of_cut(&self, cut_index: usize) -> Option<usize> {
        (**self).row_of_cut(cut_index)
    }
    fn num_rows(&self) -> usize {
        (**self).num_rows()
    }
    fn tableau_row_terms(&self, row: usize) -> (Vec<(usize, f64)>, f64) {
        (**self).tableau_row_terms(row)
    }
    fn snapshot(&self) -> Box<dyn LpEngine> {
        (**self).snapshot()
    }
}

/// 稠密引擎的 trait 实现（方法名与固有方法同名：固有方法解析优先，
/// 委托即一行）。
impl LpEngine for Simplex {
    fn solve(&mut self, bounds: &Bounds) -> LpSolution {
        self.solve(bounds)
    }
    fn resolve(&mut self, bounds: &Bounds) -> LpSolution {
        self.resolve(bounds)
    }
    fn add_row(&mut self, terms: &[(usize, f64)], op: RelOp, rhs: f64) -> Result<usize> {
        self.add_row(terms, op, rhs)
    }
    fn set_time_deadline(&mut self, deadline: Option<std::time::Instant>) {
        self.set_time_deadline(deadline)
    }
    fn time_deadline(&self) -> Option<std::time::Instant> {
        self.time_deadline()
    }
    fn set_iter_limit_hint(&mut self, hint: Option<u64>) {
        self.set_iter_limit_hint(hint)
    }
    fn iterations(&self) -> u64 {
        self.iterations()
    }
    fn basis_columns(&self) -> &[usize] {
        self.basis_columns()
    }
    fn nonbasic_columns(&self) -> Vec<NonbasicCol> {
        self.nonbasic_columns()
    }
    fn column_layout(&self) -> &ColumnLayout {
        self.column_layout()
    }
    fn row_of_constraint(&self, constraint_index: usize) -> Option<usize> {
        self.row_of_constraint(constraint_index)
    }
    fn row_of_cut(&self, cut_index: usize) -> Option<usize> {
        self.row_of_cut(cut_index)
    }
    fn num_rows(&self) -> usize {
        self.num_rows()
    }
    fn tableau_row_terms(&self, row: usize) -> (Vec<(usize, f64)>, f64) {
        self.tableau_row_terms(row)
    }
    fn snapshot(&self) -> Box<dyn LpEngine> {
        Box::new(self.clone())
    }
}

/// 连续退化迭代阈值：超过后切入 Bland（与 simplex.rs 同值）。
const DEGENERATE_STREAK_LIMIT: u32 = 50;
/// Bland 释放最小驻留步数（R8 语义，与 simplex.rs 同值）。
const BLAND_DWELL: u32 = 8;
/// 每次 solve 调用的 Bland 释放重入预算（与 simplex.rs 同值）。
const BLAND_REARM_MAX: u32 = 64;
/// 总迭代触发 Bland 的因子（与 simplex.rs 同值）。
const BLAND_TOTAL_FACTOR: u64 = 3;
/// 比值测试并列带宽绝对上限系数（simplex.rs:116-125 实测踩坑封顶，同值）。
const RATIO_BAND_CAP: f64 = 1e3;
/// 内部墙钟检查间隔（主元数，W2；未注入 deadline 零读钟）。
const TIME_CHECK_INTERVAL: u32 = 64;
/// refactorization 周期：每该数量主元整基重因子化 + 重算 x_B（类比
/// simplex.rs REFRESH_INTERVAL=200 的 B⁻¹b 重算纪律；本引擎重因子化同时
/// 清空 ETA 文件，恢复三角求解成本与数值精度）。
const REFACTOR_INTERVAL: u64 = 200;
/// LU 阈值主元比：候选 |v| ≥ 该值 × max|v| 才可入选（稳定性下界）。
const LU_THRESHOLD: f64 = 0.1;

// ===========================================================================
// LU 因子化与三角求解
// ===========================================================================

/// 稀疏 LU 因子化结果：P·B·Q = L·U（P 行置换、Q 列置换；L 单位下三角、
/// U 上三角）。`perm[i]` = LU 位置 i 对应的原始行；`col_perm[t]` = 因子
/// 列 t 对应输入矩阵的列。四份结构（行/列 × L/U）服务四个方向的三角
/// 求解：FTRAN 前代用 `l_rows`、回代用 `u_rows`（均行主序）；BTRAN 解
/// Uᵀ 用 `u_cols`（列主序升序扫）、解 Lᵀ 用 `l_cols`（列主序降序扫）。
///
/// Q 的用途（refactor 注入）：LP 基的**三角形部分先行**——单位列（松弛/
/// 人工/单非零结构列）按行序排前、多元素列按 nnz 升序排后，松弛单位列
/// 零填充消元、结构性内核规模与 fill 大降（右视 Markowitz 之外对 LP 基
/// 最有效的填充抑制；行解 P 仍由阈值主元决定）。
#[derive(Debug, Clone, Default)]
struct LuFactors {
    m: usize,
    /// LU 位置 → 原始行。
    perm: Vec<usize>,
    /// 因子列 t → 输入矩阵列。
    col_perm: Vec<usize>,
    /// L 行主序：l_rows[i] = [(列 j < i, l_ij)]（升序 j）。
    l_rows: Vec<Vec<(usize, f64)>>,
    /// L 列主序：l_cols[j] = [(行 i > j, l_ij)]（升序 i）。
    l_cols: Vec<Vec<(usize, f64)>>,
    /// U 行主序：u_rows[i] = [(列 j > i, u_ij)]（升序 j）。
    u_rows: Vec<Vec<(usize, f64)>>,
    /// U 列主序：u_cols[j] = [(行 i < j, u_ij)]（升序 i）。
    u_cols: Vec<Vec<(usize, f64)>>,
    /// U 对角。
    u_diag: Vec<f64>,
}

/// ETA（乘子形式基更新）：B_new = B_old·E，E = I 除第 r 列为该列
///（= 当时 B_old⁻¹·A_q 的**非零项**稀疏存储，行位空间；主元项单独存。
/// S4：稠密 ETA 在长解（setcover 根 LP 10k+ 迭代 × m=300）下施加成本占
/// 每迭代 ~80%，稀疏化后只触非零项）。
#[derive(Debug, Clone)]
struct Eta {
    /// 主元行位置。
    r: usize,
    /// 主元值（e_r）。
    pivot: f64,
    /// 其余非零项 (行位, 值)。
    nz: Vec<(usize, f64)>,
}

/// factorize 的工作区（SparseSimplex 持有，跨 refactor 复用容量）。
#[derive(Debug, Clone, Default)]
struct LuScratch {
    l_cols_f: Vec<Vec<(usize, f64)>>,
    u_cols_f: Vec<Vec<(usize, f64)>>,
    u_diag: Vec<f64>,
    l_rows_f: Vec<Vec<(usize, f64)>>,
    perm: Vec<usize>,
    assigned: Vec<usize>,
    x: Vec<f64>,
    in_pattern: Vec<bool>,
    pattern: Vec<usize>,
    stack: Vec<usize>,
    assigned_pat: Vec<usize>,
    unassigned_pat: Vec<usize>,
    lu_of: Vec<usize>,
}

impl LuScratch {
    fn resize(&mut self, m: usize) {
        self.l_cols_f.resize_with(m, Vec::new);
        self.u_cols_f.resize_with(m, Vec::new);
        self.l_rows_f.resize_with(m, Vec::new);
        self.u_diag.resize(m, 0.0);
        self.assigned.resize(m, usize::MAX);
        self.x.resize(m, 0.0);
        self.in_pattern.resize(m, false);
    }
}

/// 因子化失败原因（奇异基）。
#[derive(Debug, Clone)]
struct LuError(String);

/// Gilbert-Peierls 左视稀疏 LU（阈值主元 + 行填充代理破键）。
///
/// 输入：`cols[k]` = 基矩阵第 k 列的 (原始行, 值) 列表（可乱序；零值忽略）。
/// 三角求解在**原始行空间**进行（L 列条目按原始行记账，行未分配时其 LU
/// 位置未知；已分配行的解分量按 LU 位置升序先算——依赖关系 l_{r,c}≠0 ⇒
/// pos(c) < pos(r) 保证拓扑序，未分配行只依赖已分配分量、彼此独立，随后
/// 任定序（取下标升序）计算）。
/// 确定性：候选并列一律取最小原始行号；值计算升序单遍，与 DFS 遍历序无关。
/// 一次性版（单元测试用；生产路径走 [`factorize_into`] 复用版）。
#[cfg(test)]
fn factorize(
    cols: &[Vec<(usize, f64)>],
    pivot_floor: f64,
    m: usize,
) -> std::result::Result<LuFactors, LuError> {
    let mut f = LuFactors::default();
    let mut sc = LuScratch::default();
    f.resize(m);
    sc.resize(m);
    factorize_into(cols, pivot_floor, m, &mut f, &mut sc)?;
    Ok(f)
}

/// factorize 的复用版（S4）：`f` 与 `sc` 由调用方持久持有（SparseSimplex
/// 字段），clear 复用容量——每次 refactor 的 ~10×m 次堆分配降为 0。
/// 数值语义与一次性版逐位一致（同一算法，仅缓冲生命周期不同）。
fn factorize_into(
    cols: &[Vec<(usize, f64)>],
    pivot_floor: f64,
    m: usize,
    f: &mut LuFactors,
    sc: &mut LuScratch,
) -> std::result::Result<(), LuError> {
    let LuScratch {
        l_cols_f,
        u_cols_f,
        u_diag,
        l_rows_f,
        perm,
        assigned,
        x,
        in_pattern,
        pattern,
        stack,
        assigned_pat,
        unassigned_pat,
        lu_of,
    } = sc;
    for v in l_cols_f.iter_mut() {
        v.clear();
    }
    for v in u_cols_f.iter_mut() {
        v.clear();
    }
    for v in l_rows_f.iter_mut() {
        v.clear();
    }
    for v in u_diag.iter_mut() {
        *v = 0.0;
    }
    perm.clear();
    assigned.iter_mut().for_each(|a| *a = usize::MAX);
    in_pattern.iter_mut().for_each(|b| *b = false);
    x.iter_mut().for_each(|v| *v = 0.0);
    pattern.clear();
    stack.clear();
    assigned_pat.clear();
    unassigned_pat.clear();
    f.l_rows.iter_mut().for_each(|v| v.clear());
    f.l_cols.iter_mut().for_each(|v| v.clear());
    f.u_rows.iter_mut().for_each(|v| v.clear());
    f.u_cols.iter_mut().for_each(|v| v.clear());
    f.u_diag.clear();
    f.u_diag.resize(m, 0.0);
    f.perm.clear();

    for k in 0..m {
        // —— 模式 + scatter：a_k 非零行（重复行条目累加）+ 经 L 列闭包
        //    （节点 r 已分配 ⇒ 后继 = l_cols_f[assigned[r]] 的行；
        //    未分配行无后继——其列未处理）——
        pattern.clear();
        stack.clear();
        for &(r, v) in &cols[k] {
            if v != 0.0 {
                if !in_pattern[r] {
                    in_pattern[r] = true;
                    pattern.push(r);
                    stack.push(r);
                }
                x[r] += v; // 重复行条目累加（防御：正常基列无重复）
            }
        }
        while let Some(r) = stack.pop() {
            let c = assigned[r];
            if c != usize::MAX {
                for &(succ, _) in &l_cols_f[c] {
                    if !in_pattern[succ] {
                        in_pattern[succ] = true;
                        pattern.push(succ);
                        stack.push(succ);
                    }
                }
            }
        }
        // —— 值：按拓扑序单遍 L 前代（scatter 已在模式遍历中完成）——
        assigned_pat.clear();
        for &r in pattern.iter() {
            if assigned[r] != usize::MAX {
                assigned_pat.push(r);
            }
        }
        assigned_pat.sort_unstable_by_key(|&r| assigned[r]);
        for &r in assigned_pat.iter() {
            let mut acc = 0.0;
            for &(c, val) in &l_rows_f[r] {
                let pr = perm[c];
                if x[pr] != 0.0 {
                    acc += val * x[pr];
                }
            }
            x[r] -= acc;
        }
        unassigned_pat.clear();
        unassigned_pat.extend(pattern.iter().copied().filter(|&r| assigned[r] == usize::MAX));
        unassigned_pat.sort_unstable();
        for &r in unassigned_pat.iter() {
            let mut acc = 0.0;
            for &(c, val) in &l_rows_f[r] {
                let pr = perm[c];
                if x[pr] != 0.0 {
                    acc += val * x[pr];
                }
            }
            x[r] -= acc;
        }
        // —— 主元：未分配候选中 |v| ≥ LU_THRESHOLD·max|v| 且 ≥ floor 的
        //    行里取 l_rows_f 行密度最小（填充代理），并列取 |v| 大者、再
        //    并列取最小行号（实现上：先按 (fill, −|v|, 行号) 稳定选优）——
        let mut vmax = 0.0f64;
        for &r in unassigned_pat.iter() {
            let a = x[r].abs();
            if a > vmax {
                vmax = a;
            }
        }
        if !(vmax >= pivot_floor) {
            return Err(LuError(format!(
                "LU 第 {k} 列未分配候选最大 |v| = {vmax:e} 低于主元下限 {pivot_floor:e}（奇异基）"
            )));
        }
        let cut = LU_THRESHOLD * vmax;
        let mut pick: Option<(usize, f64, usize)> = None; // (填充代理, |v|, 行)
        for &r in unassigned_pat.iter() {
            let a = x[r].abs();
            if a < cut || a < pivot_floor {
                continue;
            }
            let fill = l_rows_f[r].len();
            let better = match pick {
                None => true,
                // (fill, −|v|, 行号) 字典序最小：填充少者优先，并列取 |v|
                // 大者，再并列取最小行号（升序扫描保持首个）。
                Some((bf, ba, _)) => fill < bf || (fill == bf && a > ba),
            };
            if better {
                pick = Some((fill, a, r));
            }
        }
        let (_, _, prow) = pick.expect("vmax ≥ floor 的候选必在集合内");
        // —— 记账：U 列 k（已分配模式行）+ 对角；L 列 k（其余未分配行）——
        for &r in assigned_pat.iter() {
            if x[r] != 0.0 {
                u_cols_f[k].push((r, x[r]));
            }
        }
        u_diag[k] = x[prow];
        let piv = u_diag[k];
        for &r in unassigned_pat.iter() {
            if r != prow && x[r] != 0.0 {
                let l = x[r] / piv;
                l_cols_f[k].push((r, l));
                l_rows_f[r].push((k, l));
            }
        }
        assigned[prow] = k;
        perm.push(prow);
        // 清理本列模式。
        for &r in pattern.iter() {
            in_pattern[r] = false;
            x[r] = 0.0;
        }
    }

    // —— 重映射到 LU 空间（原始行 → lu_of；列内按 LU 行升序）——
    lu_of.clear();
    lu_of.resize(m, 0usize);
    for (i, &r) in perm.iter().enumerate() {
        lu_of[r] = i;
    }
    for (k, col) in l_cols_f.iter().enumerate() {
        for &(r, v) in col.iter() {
            f.l_cols[k].push((lu_of[r], v));
        }
        f.l_cols[k].sort_unstable_by_key(|&(i, _)| i);
    }
    for (k, col) in u_cols_f.iter().enumerate() {
        for &(r, v) in col.iter() {
            f.u_cols[k].push((lu_of[r], v));
        }
        f.u_cols[k].sort_unstable_by_key(|&(i, _)| i);
    }
    for (j, col) in f.l_cols.iter().enumerate() {
        for &(i, v) in col {
            f.l_rows[i].push((j, v));
        }
    }
    for (j, col) in f.u_cols.iter().enumerate() {
        for &(i, v) in col {
            f.u_rows[i].push((j, v));
        }
    }
    f.u_diag.clone_from(u_diag);
    f.perm.clone_from(perm);
    f.col_perm.clear();
    f.col_perm.extend(0..m);
    Ok(())
}

impl LuFactors {
    /// 维度预置（复用版 factorize_into 的目标容器；S4）。
    fn resize(&mut self, m: usize) {
        self.m = m;
        self.l_rows.resize_with(m, Vec::new);
        self.l_cols.resize_with(m, Vec::new);
        self.u_rows.resize_with(m, Vec::new);
        self.u_cols.resize_with(m, Vec::new);
        self.u_diag.resize(m, 0.0);
    }

    /// 基础 FTRAN：解 B α = a（a 稠密、原始行空间）→ 写 `out`（**基位置
    /// 空间**）。PBQ = LU：LUβ = Pa 解出 β（因子列空间），α[Q 列] = β：
    /// `out[col_perm[t]] = β[t]`。行解不外显（α 按列索引）。调用方保证
    /// a 与 out 不别名。
    fn ftran_dense(&self, a: &[f64], out: &mut [f64]) {
        let m = self.m;
        let mut z = vec![0.0f64; m];
        for i in 0..m {
            z[i] = a[self.perm[i]];
        }
        for i in 0..m {
            if z[i] == 0.0 && self.l_rows[i].is_empty() {
                continue;
            }
            let mut acc = 0.0;
            for &(j, v) in &self.l_rows[i] {
                if z[j] != 0.0 {
                    acc += v * z[j];
                }
            }
            z[i] -= acc;
        }
        for i in (0..m).rev() {
            let mut acc = 0.0;
            for &(j, v) in &self.u_rows[i] {
                if z[j] != 0.0 {
                    acc += v * z[j];
                }
            }
            z[i] = (z[i] - acc) / self.u_diag[i];
        }
        for (t, &bp) in self.col_perm.iter().enumerate() {
            out[bp] = z[t];
        }
    }

    /// 基础 BTRAN：解 Bᵀ y = c（c 稠密、基位置空间）→ y（原始行空间）。
    /// Bᵀ = Q Uᵀ Lᵀ P：输入按 Q 重排（c̃[t] = c[col_perm[t]]），Uᵀ/Lᵀ
    /// 解后 y[perm[i]] = z[i]。
    fn btran(&self, c: &[f64]) -> Vec<f64> {
        let m = self.m;
        // Uᵀ w = c̃（列主序升序）。
        let mut w = vec![0.0f64; m];
        for t in 0..m {
            let cj = c[self.col_perm[t]];
            let mut acc = 0.0;
            for &(i, v) in &self.u_cols[t] {
                if w[i] != 0.0 {
                    acc += v * w[i];
                }
            }
            w[t] = (cj - acc) / self.u_diag[t];
        }
        // Lᵀ z = w（列主序降序）。
        let mut z = vec![0.0f64; m];
        for i in (0..m).rev() {
            let mut acc = 0.0;
            for &(r, v) in &self.l_cols[i] {
                if z[r] != 0.0 {
                    acc += v * z[r];
                }
            }
            z[i] = w[i] - acc;
        }
        let mut y = vec![0.0f64; m];
        for i in 0..m {
            y[self.perm[i]] = z[i];
        }
        y
    }
}

// ===========================================================================
// 求解器
// ===========================================================================

/// 定价结果（与 simplex.rs 私有枚举同构）。
enum Pricing {
    /// 合格入基列与方向（δ=+1 增向 / −1 减向）。
    Candidate(usize, f64),
    /// 无任何合格列（最优性终止条件满足）。
    Optimal,
    /// 数值失败。
    Failed(String),
}

/// 比值测试结果（与 simplex.rs Step 同构；NoBlock ⇒ Unbounded 的映射在
/// run_phase，F7）。
enum Step {
    /// 基变量在 `row` 出基（步长 t，`hit_upper` 表示其贴上界）。
    Leave { row: usize, t: f64, hit_upper: bool },
    /// 入基变量自身界翻转（步长 t）。
    Flip { t: f64 },
    /// 无任何阻挡 ⇒ 无界。
    NoBlock,
    /// 存在被主元容差过滤的阻挡候选（数值不可靠）。
    FilteredBlock(String),
}

/// 单阶段主循环结果。
enum PhaseOutcome {
    Optimal,
    Unbounded,
    LimitReached,
    TimedOut,
    /// 数值异常已重因子化重算，重试本迭代（确定性闸：至多一次）。
    Retry,
    Failed(String),
}

/// 对偶修复结果（与 simplex.rs::DualRepair 同构）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DualRepair {
    /// 已恢复原始可行（随后 Phase-2 收尾应为零迭代）。
    Repaired,
    /// 不适用/超出预算/数值异常 ⇒ 回退整基重建（语义与冷路径一致）。
    GiveUp,
}

/// 稀疏修正单纯形求解器（standalone v1；语义契约见模块文档）。
///
/// `Clone`：全字段纯数据（CSC/LU/ETA/基），克隆即精确状态（W6 快照纪律
/// 的等价物，S3 接 B&B 时复用；比稠密表格克隆更便宜）。
#[derive(Clone)]
pub struct SparseSimplex {
    tol: Tolerances,
    num_structural: usize,
    num_constraints: usize,
    sense_max: bool,
    obj_offset: f64,
    c_internal: Vec<f64>,
    layout: ColumnLayout,
    /// 静态列 CSC（[结构|松弛]）：col_start / col_rows / col_vals（行升序）。
    col_start: Vec<usize>,
    col_rows: Vec<usize>,
    col_vals: Vec<f64>,
    /// 每表格行归一化结构系数（升序；x_B 重算的行主序视图）。
    row_terms: Vec<Vec<(usize, f64)>>,
    /// 每行归一化右端 b̂。
    row_rhs: Vec<f64>,
    /// 每行关系运算符。
    row_op: Vec<RelOp>,
    /// 每行松弛列下标（Eq 行 None）。
    row_slack: Vec<Option<usize>>,
    /// 对偶还原系数（1/max|a_ij|）。
    row_scale: Vec<f64>,
    /// 原始约束 → 表格行（近零行为 None）。
    constraint_row: Vec<Option<usize>>,
    /// 割序号 → 表格行（近零割行为 None；S3 起 add_row 记账）。
    cut_row: Vec<Option<usize>>,
    /// 当前人工列数（静态段之后连续排列，每 solve 重建；add_row 可追加）。
    num_art: usize,
    /// 人工列 t → (所属行, σ)；列号 = ncols_static + t。
    art_info: Vec<(usize, f64)>,
    basis: Vec<usize>,
    in_basis: Vec<bool>,
    at_upper: Vec<bool>,
    x_b: Vec<f64>,
    bounds: Bounds,
    /// 上次 solve/resolve 的有效界（resolve 的 snap 基准；与 simplex.rs
    /// 的 prev_bounds 同款 swap 协议）。
    prev_bounds: Bounds,
    lu: LuFactors,
    /// factorize 工作区（跨 refactor 复用容量；S4 分配消除）。
    lu_scratch: LuScratch,
    /// refactor 的基列收集缓冲（复用容量）。
    refactor_cols: Vec<Vec<(usize, f64)>>,
    etas: Vec<Eta>,
    /// add_row 追加了行/列后的脏标记：下一次 solve/resolve 入口整基重
    /// 因子化 + 重算 x_B（割批共享一次；F6『数值上宜随后重构』的落法）。
    factorization_dirty: bool,
    /// 上次 Phase-1 未清理即收场 / add_row 新建人工列（resolve 必须重跑
    /// Phase-1 的欠账记账，语义同 simplex.rs）。
    phase1_pending: bool,
    has_run: bool,
    infeasible_flag: bool,
    #[allow(dead_code)]
    numerical_error: Option<String>,
    iterations_total: u64,
    iter_call: u64,
    degenerate_streak: u32,
    bland: bool,
    bland_permanent: bool,
    bland_dwell: u32,
    bland_rearms: u32,
    iter_limit_hint: Option<u64>,
    time_deadline: Option<std::time::Instant>,
    pivots_since_time_check: u32,
    pivots_since_refactor: u64,
    /// 数值异常重试闸：刚重因子化 ⇒ 不再重试（防循环；确定性）。
    just_refactored: bool,
    // —— 热路径暂存（容量跨迭代复用）——
    /// FTRAN 入参（稠密列散射）。
    ftran_in: Vec<f64>,
    /// FTRAN 结果（行位空间稠密 α）。
    alpha: Vec<f64>,
    /// BTRAN 结果（行位空间稠密 y）。
    y_vec: Vec<f64>,
    /// c_B 暂存。
    cb: Vec<f64>,
    /// 稠密右端工作区（x_B 重算）。
    work: Vec<f64>,
    /// 比值测试候选暂存 (t, 行, hit_upper, |主元|, 出基列)。
    cand: Vec<(f64, usize, bool, f64, usize)>,
}

impl SparseSimplex {
    /// 由 LP 松弛构建（整数性忽略；校验与 simplex.rs::Simplex::new 同款）。
    pub fn new(problem: &Problem, tol: Tolerances) -> Result<Self> {
        let n = problem.vars.len();
        for v in &problem.vars {
            if !v.lb.is_finite() {
                return Err(MilpError::InvalidModel(format!(
                    "变量 {} 的下界非有限（本 crate 不支持 -∞ 下界）",
                    v.name
                )));
            }
            if !(v.ub.is_finite() || v.ub == f64::INFINITY) {
                return Err(MilpError::InvalidModel(format!(
                    "变量 {} 的上界非法（允许有限数或 +∞）",
                    v.name
                )));
            }
        }
        if !problem.obj_offset.is_finite() {
            return Err(MilpError::InvalidModel("obj_offset 非有限数".to_string()));
        }
        for &(j, c) in &problem.obj.terms {
            if j >= n {
                return Err(MilpError::InvalidModel(format!("目标含越界变量下标 {j}")));
            }
            if !c.is_finite() {
                return Err(MilpError::InvalidModel(format!(
                    "目标变量 {j} 的系数非有限数"
                )));
            }
        }
        for con in &problem.constraints {
            if !con.rhs.is_finite() {
                return Err(MilpError::InvalidModel(format!(
                    "约束 {:?} 右端非有限数",
                    con.name
                )));
            }
            for &(j, a) in &con.expr.terms {
                if j >= n {
                    return Err(MilpError::InvalidModel(format!(
                        "约束 {:?} 含越界变量下标 {j}",
                        con.name
                    )));
                }
                if !a.is_finite() {
                    return Err(MilpError::InvalidModel(format!(
                        "约束 {:?} 变量 {j} 系数非有限数",
                        con.name
                    )));
                }
            }
        }

        let sense_max = problem.sense == Sense::Maximize;
        let sign = if sense_max { -1.0 } else { 1.0 };
        let mut c_internal = vec![0.0; n];
        for &(j, c) in &problem.obj.terms {
            c_internal[j] = sign * c;
        }

        let mut s = SparseSimplex {
            tol,
            num_structural: n,
            num_constraints: problem.constraints.len(),
            sense_max,
            obj_offset: problem.obj_offset,
            c_internal,
            layout: ColumnLayout {
                num_structural: n,
                slack: Vec::new(),
            },
            col_start: Vec::new(),
            col_rows: Vec::new(),
            col_vals: Vec::new(),
            row_terms: Vec::new(),
            row_rhs: Vec::new(),
            row_op: Vec::new(),
            row_slack: Vec::new(),
            row_scale: Vec::new(),
            constraint_row: vec![None; problem.constraints.len()],
            cut_row: Vec::new(),
            num_art: 0,
            art_info: Vec::new(),
            basis: Vec::new(),
            in_basis: Vec::new(),
            at_upper: Vec::new(),
            x_b: Vec::new(),
            bounds: Bounds {
                lb: vec![0.0; n],
                ub: vec![f64::INFINITY; n],
            },
            prev_bounds: Bounds {
                lb: Vec::new(),
                ub: Vec::new(),
            },
            lu: LuFactors::default(),
            lu_scratch: LuScratch::default(),
            refactor_cols: Vec::new(),
            etas: Vec::new(),
            factorization_dirty: false,
            phase1_pending: false,
            has_run: false,
            infeasible_flag: false,
            numerical_error: None,
            iterations_total: 0,
            iter_call: 0,
            degenerate_streak: 0,
            bland: false,
            bland_permanent: false,
            bland_dwell: 0,
            bland_rearms: 0,
            iter_limit_hint: None,
            time_deadline: None,
            pivots_since_time_check: 0,
            pivots_since_refactor: 0,
            just_refactored: false,
            ftran_in: Vec::new(),
            alpha: Vec::new(),
            y_vec: Vec::new(),
            cb: Vec::new(),
            work: Vec::new(),
            cand: Vec::new(),
        };

        // —— 行构建（归一化 + 近零行剔除 + 松弛列布局；同 simplex.rs）——
        for (ci, con) in problem.constraints.iter().enumerate() {
            let g = con
                .expr
                .terms
                .iter()
                .map(|&(_, a)| a.abs())
                .fold(0.0_f64, f64::max);
            if g <= s.tol.zero {
                let ok = match con.op {
                    RelOp::Le => con.rhs >= -s.tol.feasibility,
                    RelOp::Ge => con.rhs <= s.tol.feasibility,
                    RelOp::Eq => con.rhs.abs() <= s.tol.feasibility,
                };
                if !ok {
                    s.infeasible_flag = true;
                }
                s.constraint_row[ci] = None;
                continue;
            }
            let row = s.row_terms.len();
            if con.op == RelOp::Eq {
                s.row_slack.push(None);
            } else {
                let col = n + s.layout.slack.len();
                s.layout.slack.push(SlackCol {
                    origin: RowOrigin::Constraint(ci),
                    sign: if con.op == RelOp::Le { 1.0 } else { -1.0 },
                });
                s.row_slack.push(Some(col));
            }
            s.row_terms
                .push(con.expr.terms.iter().map(|&(j, a)| (j, a / g)).collect());
            s.row_rhs.push(con.rhs / g);
            s.row_op.push(con.op);
            s.row_scale.push(1.0 / g);
            s.constraint_row[ci] = Some(row);
        }

        // —— 静态列 CSC（结构列 + 松弛列；行升序）——
        let ncs = s.ncols_static();
        let mut counts = vec![0usize; ncs];
        for row in &s.row_terms {
            for &(j, _) in row {
                counts[j] += 1;
            }
        }
        for &sl in s.row_slack.iter() {
            if let Some(c) = sl {
                counts[c] += 1;
            }
        }
        s.col_start = Vec::with_capacity(ncs + 1);
        let mut acc = 0usize;
        s.col_start.push(0);
        for &c in &counts {
            acc += c;
            s.col_start.push(acc);
        }
        s.col_rows = vec![0; acc];
        s.col_vals = vec![0.0; acc];
        let mut cursor = s.col_start.clone();
        for (i, row) in s.row_terms.iter().enumerate() {
            for &(j, a) in row {
                let p = cursor[j];
                s.col_rows[p] = i;
                s.col_vals[p] = a;
                cursor[j] += 1;
            }
        }
        let n0 = n;
        for (t, &sl) in s.row_slack.iter().enumerate() {
            if let Some(c) = sl {
                let p = cursor[c];
                s.col_rows[p] = t;
                s.col_vals[p] = s.layout.slack[c - n0].sign;
                cursor[c] += 1;
            }
        }

        let m = s.row_terms.len();
        s.basis = vec![usize::MAX; m];
        s.in_basis = vec![false; ncs];
        s.at_upper = vec![false; ncs];
        s.x_b = vec![0.0; m];
        s.ftran_in = vec![0.0; m];
        s.alpha = vec![0.0; m];
        s.y_vec = vec![0.0; m];
        s.cb = vec![0.0; m];
        s.work = vec![0.0; m];
        s.bounds = Bounds::from_problem(problem);
        // 预置一致初始基（默认界；保证访问器在 solve 前可用，同 simplex.rs）。
        s.rebuild_initial_basis();
        Ok(s)
    }

    // ---- 访问器（cuts.rs 消费面；列空间 [结构|松弛]、人工列不外露）----

    /// 表格行数（= 保留行数；S2 无割行）。
    pub fn num_rows(&self) -> usize {
        self.basis.len()
    }

    /// 原始约束下标 → 表格行；近零行返回 None。
    pub fn row_of_constraint(&self, constraint_index: usize) -> Option<usize> {
        self.constraint_row
            .get(constraint_index)
            .copied()
            .flatten()
    }

    /// 割序号 → 表格行；近零割行为 None。
    pub fn row_of_cut(&self, cut_index: usize) -> Option<usize> {
        self.cut_row.get(cut_index).copied().flatten()
    }

    /// 各行基列下标（结构+松弛列空间）。
    pub fn basis_columns(&self) -> &[usize] {
        &self.basis
    }

    /// 当前全部非基列及其所贴界（不含人工列）。
    pub fn nonbasic_columns(&self) -> Vec<NonbasicCol> {
        let ncols = self.ncols_static();
        let mut out = Vec::new();
        for j in 0..ncols {
            if !self.in_basis[j] {
                out.push(NonbasicCol {
                    col: j,
                    at_upper: self.at_upper[j],
                });
            }
        }
        out
    }

    /// 列布局（结构列/松弛列含义）。
    pub fn column_layout(&self) -> &ColumnLayout {
        &self.layout
    }

    /// 第 `row` 行的稀疏表格行：(非零 (列, e_rowᵀB⁻¹A_j) 项 [结构+松弛段、
    /// 列升序]， 该行基变量当前值)。修正单纯形实现：ρ = e_rowᵀB⁻¹ 经
    /// BTRAN，行项 = ρ·A_j 按列现算（与稠密引擎 B⁻¹A 行读取语义一致；
    /// 对账在 ≥1e-6 相对容差口径，不承诺逐位同）。越界行返回 (空, NaN)。
    pub fn tableau_row_terms(&self, row: usize) -> (Vec<(usize, f64)>, f64) {
        let m = self.basis.len();
        if row >= m {
            return (Vec::new(), f64::NAN);
        }
        let mut e = vec![0.0f64; m];
        e[row] = 1.0;
        let rho = self.btran_full(&e);
        let ncols = self.ncols_static();
        let mut terms = Vec::new();
        for j in 0..ncols {
            let mut dot = 0.0;
            for k in self.col_start[j]..self.col_start[j + 1] {
                let r = self.col_rows[k];
                if rho[r] != 0.0 {
                    dot += rho[r] * self.col_vals[k];
                }
            }
            if dot != 0.0 {
                terms.push((j, dot));
            }
        }
        (terms, self.x_b[row])
    }

    /// 累计迭代数。
    pub fn iterations(&self) -> u64 {
        self.iterations_total
    }

    /// 设置单次调用迭代上限提示（R8 割循环接线点，语义同 simplex.rs）。
    pub fn set_iter_limit_hint(&mut self, hint: Option<u64>) {
        self.iter_limit_hint = hint;
    }

    /// 内部墙钟时限（W2 条件读钟：未注入零读钟；语义同 simplex.rs）。
    pub fn set_time_deadline(&mut self, deadline: Option<std::time::Instant>) {
        self.time_deadline = deadline;
    }

    /// 当前内部墙钟时限（S3 透传用）。
    pub fn time_deadline(&self) -> Option<std::time::Instant> {
        self.time_deadline
    }

    // ---- 求解入口 ----

    /// 两阶段冷启动求解（standalone：每次调用整基重建，忽略此前状态）。
    pub fn solve(&mut self, bounds: &Bounds) -> LpSolution {
        self.start_call();
        if self.infeasible_flag {
            return self.build_solution(LpStatus::Infeasible);
        }
        if let Some(st) = self.validate_bounds_arg(bounds) {
            return self.build_solution(st);
        }
        self.bounds.lb.clone_from(&bounds.lb);
        self.bounds.ub.clone_from(&bounds.ub);
        self.rebuild_initial_basis();
        let st = self.run_two_phase(true);
        self.has_run = true;
        self.build_solution(st)
    }

    /// 热启动求解（S3）：语义对齐 simplex.rs::resolve——
    /// (1)(2) 非基贴界双向检查（相对上次有效界 `prev_bounds`：收紧越界
    /// snap 到收紧界 / 严格内点 snap 到较近界、等距按 rc 符号选侧），
    /// 每个 snap 的 Δx_j 按 B⁻¹A_j（FTRAN）传播；
    /// (3) 基变量越界或 Phase-1 欠账（含 add_row 新建人工列）⇒ 修复：
    /// a) **对偶修复快路径**（无欠账、基中无人工列、静态非基对偶可行）：
    /// 有界变量对偶单纯形若干次旋转恢复原始可行（每次出基落到被违反
    /// 界；FTRAN/BTRAN 即时现算）；任何异常 GiveUp 回退 b)；
    /// b) 整基重建 + Phase-1（冷路径，语义与 simplex.rs 完全一致）。
    pub fn resolve(&mut self, bounds: &Bounds) -> LpSolution {
        self.start_call();
        if self.infeasible_flag {
            return self.build_solution(LpStatus::Infeasible);
        }
        if let Some(st) = self.validate_bounds_arg(bounds) {
            return self.build_solution(st);
        }
        if !self.has_run {
            // 尚无可用热状态：退化为冷启动（整基重建两阶段）。
            self.bounds.lb.clone_from(&bounds.lb);
            self.bounds.ub.clone_from(&bounds.ub);
            self.rebuild_initial_basis();
            let st = self.run_two_phase(true);
            self.has_run = true;
            return self.build_solution(st);
        }
        // 热状态可用：add_row 脏标记或 ETA 越限 ⇒ 先整基重因子化 + 重算
        // x_B（基保持，等价『随后重构』纪律）。
        if self.factorization_dirty || self.eta_over_limit() {
            if let Err(msg) = self.refactor() {
                self.numerical_error = Some(msg);
                return self.build_solution(LpStatus::IterationLimit);
            }
        }
        // (1)(2) 非基贴界双向检查 + 传播（prev_bounds 换入上次有效界）。
        std::mem::swap(&mut self.bounds, &mut self.prev_bounds);
        self.bounds.lb.clone_from(&bounds.lb);
        self.bounds.ub.clone_from(&bounds.ub);
        self.snap_nonbasics();
        // (3) 基变量越界 / Phase-1 欠账检查 → 对偶修复或整基重建。
        let need_repair = self.basics_violated() || self.phase1_pending;
        let mut dual_repaired = false;
        if need_repair
            && !self.phase1_pending
            && !self.artificials_in_basis()
            && self.dual_repair() == DualRepair::Repaired
        {
            dual_repaired = true;
        } else if need_repair {
            self.rebuild_initial_basis();
        }
        let st = self.run_two_phase(need_repair && !dual_repaired);
        self.has_run = true;
        self.build_solution(st)
    }

    /// 非基列贴界双向检查（与 simplex.rs::snap_nonbasics 同款规则；传播
    /// 经 FTRAN 现算 B⁻¹A_j）。
    fn snap_nonbasics(&mut self) {
        let n = self.num_structural;
        let feas = self.tol.feasibility;
        // 等距并列按 rc 符号选侧所需的 y（惰性：首次用到才 BTRAN）。
        let mut y_cache: Option<Vec<f64>> = None;
        for j in 0..n {
            if self.in_basis[j] {
                continue;
            }
            let old_v = if self.at_upper[j] {
                self.prev_bounds.ub[j]
            } else {
                self.prev_bounds.lb[j]
            };
            let nlb = self.bounds.lb[j];
            let nub = self.bounds.ub[j];
            let inside = (old_v - nlb).abs() > feas
                && !(nub.is_finite() && (nub - old_v).abs() <= feas);
            let (target, to_upper) = if self.at_upper[j] {
                if nub.is_finite() && nub < old_v - feas {
                    (nub, true) // 收紧越界：snap 到新 ub
                } else if inside {
                    self.nearest_or_rc(j, old_v, nlb, nub, &mut y_cache)
                } else if nub.is_finite() {
                    (nub, true) // 仍贴上界
                } else {
                    (nlb, false) // 新上界 +∞：只能贴下界
                }
            } else {
                if nlb > old_v + feas {
                    (nlb, false) // 收紧越界：snap 到新 lb
                } else if inside {
                    self.nearest_or_rc(j, old_v, nlb, nub, &mut y_cache)
                } else {
                    (nlb, false) // 仍贴下界
                }
            };
            let delta = target - old_v;
            if delta != 0.0 {
                self.ftran_col(j);
                let m = self.basis.len();
                for k in 0..m {
                    let a = self.alpha[k];
                    if a != 0.0 {
                        self.x_b[k] -= delta * a;
                    }
                }
            }
            self.at_upper[j] = to_upper;
        }
    }

    /// 松弛失贴时的 snap 目标：较近的界；等距（带宽内）按检验数符号
    /// 选侧（rc<0 选下界、rc≥0 选上界；与 simplex.rs::nearest_or_rc 同款）。
    fn nearest_or_rc(
        &mut self,
        j: usize,
        v: f64,
        lb: f64,
        ub: f64,
        y_cache: &mut Option<Vec<f64>>,
    ) -> (f64, bool) {
        let d_lb = (v - lb).abs();
        let d_ub = if ub.is_finite() { (ub - v).abs() } else { f64::INFINITY };
        let tie = self.tol.feasibility;
        if d_lb < d_ub - tie {
            return (lb, false);
        }
        if d_ub < d_lb - tie {
            return (ub, true);
        }
        if y_cache.is_none() {
            let m = self.basis.len();
            let mut cb = std::mem::take(&mut self.cb);
            cb.resize(m, 0.0);
            for i in 0..m {
                cb[i] = self.cost_of(self.basis[i], 2);
            }
            let y = self.btran_full(&cb);
            self.cb = cb;
            *y_cache = Some(y);
        }
        let y = y_cache.as_ref().unwrap();
        let mut dot = 0.0;
        for k in self.col_start[j]..self.col_start[j + 1] {
            let yr = y[self.col_rows[k]];
            if yr != 0.0 {
                dot += yr * self.col_vals[k];
            }
        }
        let rc = self.cost_of(j, 2) - dot;
        if rc < 0.0 {
            (lb, false)
        } else {
            (ub, true)
        }
    }

    /// 基变量越界检查（与 simplex.rs::basics_violated 同款：残基人工变量
    /// 必须为 0、基变量在界 ±feas 内、值有限）。
    fn basics_violated(&self) -> bool {
        let ncs = self.ncols_static();
        let feas = self.tol.feasibility;
        for i in 0..self.basis.len() {
            let c = self.basis[i];
            if !self.x_b[i].is_finite() {
                return true;
            }
            if c >= ncs {
                if self.x_b[i].abs() > feas {
                    return true;
                }
            } else {
                let (lb, ub) = self.col_bounds(c);
                if self.x_b[i] < lb - feas || self.x_b[i] > ub + feas {
                    return true;
                }
            }
        }
        false
    }

    /// 基中是否存在人工列（对偶快路径前置排除）。
    fn artificials_in_basis(&self) -> bool {
        let ncs = self.ncols_static();
        self.basis.iter().any(|&c| c >= ncs)
    }

    /// 静态非基列对偶可行性预检（与 simplex.rs::dual_feasible 同款）。
    fn dual_feasible(&mut self) -> bool {
        let ncs = self.ncols_static();
        let m = self.basis.len();
        for i in 0..m {
            self.cb[i] = self.cost_of(self.basis[i], 2);
        }
        let y = self.btran_full(&self.cb);
        for j in 0..ncs {
            if self.in_basis[j] {
                continue;
            }
            let (lb, ub) = self.col_bounds(j);
            if ub - lb <= self.tol.feasibility {
                continue; // 固定列：无对偶语义
            }
            let mut dot = 0.0;
            for k in self.col_start[j]..self.col_start[j + 1] {
                let yr = y[self.col_rows[k]];
                if yr != 0.0 {
                    dot += yr * self.col_vals[k];
                }
            }
            let rc = self.cost_of(j, 2) - dot;
            if !rc.is_finite() {
                return false;
            }
            if self.at_upper[j] {
                if rc > self.tol.optimality {
                    return false;
                }
            } else if rc < -self.tol.optimality {
                return false;
            }
        }
        true
    }

    /// 有界变量对偶单纯形修复（与 simplex.rs::dual_repair 同款算法，三角
    /// 求解即时现算版）：从（对偶可行、原始违约）基出发，每次选违反最重
    /// 基行出基（落到被违反界）、对偶比值测试选入基列，标准主元，直至
    /// 全部基变量回界内。任何异常（无可入列 / 步长非法 / 超预算 / 检验数
    /// 非有限 / 主元过小 / 重因子化失败）⇒ GiveUp 回退整基重建——快路径
    /// 只加速、绝不改变结果。
    fn dual_repair(&mut self) -> DualRepair {
        const BUDGET_BASE: u32 = 2_000;
        const DEGEN_BLAND: u32 = 32;
        if !self.dual_feasible() {
            return DualRepair::GiveUp;
        }
        let m = self.basis.len();
        let ncs = self.ncols_static();
        let budget = BUDGET_BASE.saturating_add((m as u32).saturating_mul(50));
        let bland_at = budget / 2;
        let mut iters: u32 = 0;
        let mut degen: u32 = 0;
        loop {
            // y = BTRAN(c_B)（检验数新鲜度：上一轮主元已改基）。
            for i in 0..m {
                self.cb[i] = self.cost_of(self.basis[i], 2);
            }
            let y = self.btran_full(&self.cb);
            // (1) 选相对违反最重的基行（绝对 feasibility 过滤 + 界幅归一
            // 排序，与 simplex.rs 同款）。
            let mut worst: Option<(usize, f64, bool)> = None;
            for i in 0..m {
                let (lb, ub) = self.col_bounds(self.basis[i]);
                let v_lb = lb - self.x_b[i];
                let v_ub = self.x_b[i] - ub;
                let (viol, below) = if v_lb > v_ub { (v_lb, true) } else { (v_ub, false) };
                if viol.is_finite() && viol <= self.tol.feasibility {
                    continue;
                }
                let bound = if below { lb } else { ub };
                let norm = viol / (1.0 + bound.abs());
                match worst {
                    Some((_, w, _)) if w >= norm => {}
                    _ => worst = Some((i, norm, below)),
                }
            }
            let Some((r, _, below)) = worst else {
                return DualRepair::Repaired;
            };
            if iters >= budget {
                return DualRepair::GiveUp;
            }
            let use_bland = iters >= bland_at || degen >= DEGEN_BLAND;
            // (2) 对偶比值测试：表格行 r 经 BTRAN(e_r) + 行收集；方向可行
            //（below ⇒ δ·a<0 / above ⇒ δ·a>0）的非基列中 ratio =
            // 有符号 rc/|a| 最小者；并列 |主元| 大者、退化期 Bland 最小列。
            let mut e = vec![0.0f64; m];
            e[r] = 1.0;
            let rho = self.btran_full(&e);
            let mut best: Option<(f64, f64, usize, f64)> = None; // (ratio, |a|, q, delta)
            for j in 0..ncs {
                if self.in_basis[j] {
                    continue;
                }
                let (lb, ub) = self.col_bounds(j);
                if ub - lb <= self.tol.feasibility {
                    continue;
                }
                let mut a = 0.0;
                for k in self.col_start[j]..self.col_start[j + 1] {
                    let rr = rho[self.col_rows[k]];
                    if rr != 0.0 {
                        a += rr * self.col_vals[k];
                    }
                }
                if a.abs() < self.tol.pivot {
                    continue;
                }
                let delta = if self.at_upper[j] { -1.0 } else { 1.0 };
                let da = delta * a;
                if below && da >= 0.0 {
                    continue;
                }
                if !below && da <= 0.0 {
                    continue;
                }
                let mut dot = 0.0;
                for k in self.col_start[j]..self.col_start[j + 1] {
                    let yr = y[self.col_rows[k]];
                    if yr != 0.0 {
                        dot += yr * self.col_vals[k];
                    }
                }
                let rc = self.cost_of(j, 2) - dot;
                if !rc.is_finite() {
                    return DualRepair::GiveUp;
                }
                let signed = if self.at_upper[j] { -rc } else { rc };
                let ratio = (signed / a.abs()).max(0.0);
                let better = match best {
                    None => true,
                    Some((br, ba, bq, _)) => {
                        if use_bland {
                            j < bq
                        } else {
                            ratio < br
                                || (ratio == br && (a.abs() > ba || (a.abs() == ba && j < bq)))
                        }
                    }
                };
                if better {
                    best = Some((ratio, a.abs(), j, delta));
                }
            }
            let Some((_, _, q, delta)) = best else {
                // 无方向可行列（对偶无界 ⇒ 原始不可行）：保守回退 Phase-1。
                return DualRepair::GiveUp;
            };
            // (3) 步长 t：出基变量恰好落到被违反的界；FTRAN 入基列取主元。
            let mut a = 0.0;
            for k in self.col_start[q]..self.col_start[q + 1] {
                let rr = rho[self.col_rows[k]];
                if rr != 0.0 {
                    a += rr * self.col_vals[k];
                }
            }
            let (lb_r, ub_r) = self.col_bounds(self.basis[r]);
            let (target, hit_upper) = if below { (lb_r, false) } else { (ub_r, true) };
            let t = (target - self.x_b[r]) / (-delta * a);
            if !t.is_finite() || t < 0.0 {
                return DualRepair::GiveUp;
            }
            self.ftran_col(q);
            if self.alpha[r].abs() < self.tol.pivot {
                // 数值边缘：先重因子化重算再试一次（确定性的单次重试）。
                if !self.just_refactored {
                    if self.refactor().is_err() {
                        return DualRepair::GiveUp;
                    }
                    self.ftran_col(q);
                    if self.alpha[r].abs() < self.tol.pivot {
                        return DualRepair::GiveUp;
                    }
                } else {
                    return DualRepair::GiveUp;
                }
            }
            if self.do_pivot(q, r, delta, t, hit_upper).is_err() {
                return DualRepair::GiveUp;
            }
            iters += 1;
            self.iter_call += 1;
            self.iterations_total += 1;
            if t <= self.tol.feasibility {
                degen += 1;
            } else {
                degen = 0;
            }
        }
    }

    // ---- 内部：列与成本 ----

    fn ncols_static(&self) -> usize {
        self.num_structural + self.layout.slack.len()
    }

    /// 列当前上下界（结构列用当前 bounds；松弛/人工恒 [0,+∞)）。
    fn col_bounds(&self, j: usize) -> (f64, f64) {
        if j < self.num_structural {
            (self.bounds.lb[j], self.bounds.ub[j])
        } else {
            (0.0, f64::INFINITY)
        }
    }

    fn nonbasic_value(&self, j: usize) -> f64 {
        let (lb, ub) = self.col_bounds(j);
        if self.at_upper[j] {
            ub
        } else {
            lb
        }
    }

    /// 列在指定阶段的成本（与 simplex.rs::cost_of 同构）。
    fn cost_of(&self, j: usize, phase: u8) -> f64 {
        let ncs = self.ncols_static();
        if phase == 1 {
            if j >= ncs {
                1.0
            } else {
                0.0
            }
        } else if j < self.num_structural {
            self.c_internal[j]
        } else {
            0.0
        }
    }

    // ---- 内部：因子化 / 三角求解 / ETA ----

    /// ETA 越限判定（S4）：**flops 平衡 + 绝对上限**——ETA 施加成本
    /// ≈ Σ|eta.nnz|/次三角求解，重因子化成本 ≈ O(LU fill) 一次性，比例
    /// 失衡即重构。下限 64 为实测定档（setcover-300x600 全程 A/B：
    /// floor 16/8 ⇒ 重构过频反亏 [124.9s]，64 最优 [47.1s]，128/200 ⇒
    /// ETA 施加占优 [99.7-124.9s]）；上限 200 沿用稠密引擎
    /// REFRESH_INTERVAL 语义。确定性：只依赖 m 与 LU 规模。
    fn eta_over_limit(&self) -> bool {
        let m = self.basis.len().max(1);
        let lu_nnz = self.lu.l_rows.len() + self.lu.u_cols.len() + self.lu.u_diag.len();
        let balance = (4 * lu_nnz) / m;
        let limit = (balance as u64).clamp(64, REFACTOR_INTERVAL);
        self.etas.len() as u64 >= limit
    }

    /// 整基重因子化 + 从原始数据重算 x_B（漂移双保险；同时清空 ETA）。
    ///
    /// 因子化列序（Q）：单位列（松弛/人工/单非零结构列）按行序排前
    ///（零填充消元），多元素列按 (nnz 升序, 基位置升序) 排后——LP 基的
    /// 三角形部分先行、结构性内核的规模与 fill 大降（确定性排序）。
    fn refactor(&mut self) -> std::result::Result<(), String> {
        let m = self.basis.len();
        let ncs = self.ncols_static();
        // 列权：(行 or nnz, 基位置)。
        let mut units: Vec<(usize, usize)> = Vec::new(); // (唯一行, bp)
        let mut multis: Vec<(usize, usize)> = Vec::new(); // (nnz, bp)
        for (bp, &j) in self.basis.iter().enumerate() {
            if j < ncs {
                let nnz = self.col_start[j + 1] - self.col_start[j];
                if nnz == 1 {
                    units.push((self.col_rows[self.col_start[j]], bp));
                } else {
                    multis.push((nnz, bp));
                }
            } else {
                let (row, _) = self.art_info[j - ncs];
                units.push((row, bp));
            }
        }
        units.sort_unstable();
        multis.sort_unstable();
        let mut order: Vec<usize> = units.iter().map(|&(_, bp)| bp).collect();
        order.extend(multis.iter().map(|&(_, bp)| bp));
        // 按因子化列序收集基列（原始行空间稀疏列）。
        self.lu_scratch.resize(m);
        self.lu.resize(m);
        self.refactor_cols.clear();
        self.refactor_cols.resize_with(m, Vec::new);
        for (t, &bp) in order.iter().enumerate() {
            let j = self.basis[bp];
            let col = &mut self.refactor_cols[t];
            col.clear();
            if j < ncs {
                for k in self.col_start[j]..self.col_start[j + 1] {
                    if self.col_vals[k] != 0.0 {
                        col.push((self.col_rows[k], self.col_vals[k]));
                    }
                }
            } else {
                let (row, sigma) = self.art_info[j - ncs];
                col.push((row, sigma));
            }
        }
        factorize_into(
            &self.refactor_cols,
            self.tol.pivot,
            m,
            &mut self.lu,
            &mut self.lu_scratch,
        )
        .map_err(|e| e.0)?;
        self.lu.col_perm = order;
        self.etas.clear();
        self.pivots_since_refactor = 0;
        // x_B = B⁻¹(b̂ − N·x_N)：非基结构列贴界、非基松弛/人工贴 0
        //（行主序扣减，与 simplex.rs::refresh_values 同序）。
        for i in 0..m {
            let mut acc = 0.0;
            for &(j, a) in &self.row_terms[i] {
                if !self.in_basis[j] {
                    let v = if self.at_upper[j] {
                        self.bounds.ub[j]
                    } else {
                        self.bounds.lb[j]
                    };
                    if v != 0.0 {
                        acc += a * v;
                    }
                }
            }
            self.work[i] = self.row_rhs[i] - acc;
        }
        let mut xb = vec![0.0; m];
        self.lu.ftran_dense(&self.work, &mut xb);
        for v in xb.iter() {
            if !v.is_finite() {
                self.lu = LuFactors::default();
                self.etas.clear();
                return Err("x_B 重算出现 NaN/Inf".to_string());
            }
        }
        self.x_b = xb;
        self.just_refactored = true;
        self.factorization_dirty = false;
        Ok(())
    }

    /// FTRAN：α = B⁻¹A_q（含 ETA 顺序施加 E⁻¹）→ self.alpha。
    fn ftran_col(&mut self, q: usize) {
        for v in self.ftran_in.iter_mut() {
            *v = 0.0;
        }
        if q < self.ncols_static() {
            for k in self.col_start[q]..self.col_start[q + 1] {
                self.ftran_in[self.col_rows[k]] = self.col_vals[k];
            }
        } else {
            let (row, sigma) = self.art_info[q - self.ncols_static()];
            self.ftran_in[row] = sigma;
        }
        let mut alpha = std::mem::take(&mut self.alpha);
        self.lu.ftran_dense(&self.ftran_in, &mut alpha);
        // ETA：旧 → 新顺序施加 E⁻¹（α_r' = α_r/e_r；其余 α_i −= e_i·α_r'）。
        for eta in &self.etas {
            let t = alpha[eta.r] / eta.pivot;
            if t != 0.0 {
                for &(i, e) in eta.nz.iter() {
                    alpha[i] -= e * t;
                }
            }
            alpha[eta.r] = t;
        }
        self.alpha = alpha;
    }

    /// BTRAN + ETA：解 Bᵀ y = c（c 稠密行位空间）→ y。
    /// ETA 逆序施加 (Eᵀ)⁻¹（新 → 旧：w_i = v_i (i≠r)；w_r =
    /// (v_r − Σ e_i·v_i)/e_r），再基础 LU BTRAN。
    fn btran_full(&self, c: &[f64]) -> Vec<f64> {
        let mut v = c.to_vec();
        for eta in self.etas.iter().rev() {
            let mut acc = 0.0;
            for &(i, e) in eta.nz.iter() {
                acc += e * v[i];
            }
            v[eta.r] = (v[eta.r] - acc) / eta.pivot;
        }
        self.lu.btran(&v)
    }

    /// 追加割行（S3；terms 为结构变量空间稀疏行，原尺度不缩放，scale=1）；
    /// 语义与 simplex.rs::add_row 同款：Le 按松弛列（sign=+1）、Ge 按剩余
    /// 列（sign=−1）、Eq 无松弛；残差违反时建人工列（σ 规则同构建期）并置
    /// `phase1_pending`；返回割序号 k（0 起连续）。
    ///
    /// 修正单纯形实现注（F6）：不是免费操作——新行使基从 m×m 增广为
    /// (m+1)×(m+1)（下加边 [[B,0],[κᵀ,σ]]），旧行 x_B 数学上不变（上式逆
    /// 的分块结构保证），本实现**不即时更新因子**：置 `factorization_dirty`，
    /// 由下一次 solve/resolve 入口整基重因子化 + 从原始数据重算 x_B（割批
    /// 共享一次重构，『数值上宜随后重构』的最简落法）。新松弛列在 CSC 末尾
    /// 追加（新松弛恒为最后一个静态列——无需稠密引擎的中段插入+人工列
    /// 段右移技巧；人工列号随 ncs 增长统一 +1，经 basis/in_basis 同步）。
    /// 近零割行（防御；cuts.rs 数值纪律下不可达）不入表，按原尺度 rhs 判定。
    pub fn add_row(&mut self, terms: &[(usize, f64)], op: RelOp, rhs: f64) -> Result<usize> {
        let n = self.num_structural;
        if !rhs.is_finite() {
            return Err(MilpError::InvalidModel("add_row 右端非有限数".to_string()));
        }
        // 归并排序 + 丢弃精确零系数（防御；与 simplex.rs 同款）。
        let mut merged: Vec<(usize, f64)> = Vec::new();
        {
            let mut sorted: Vec<(usize, f64)> = Vec::with_capacity(terms.len());
            for &(j, c) in terms {
                if j >= n {
                    return Err(MilpError::InvalidModel(format!(
                        "add_row 含越界变量下标 {j}"
                    )));
                }
                if !c.is_finite() {
                    return Err(MilpError::InvalidModel(format!(
                        "add_row 变量 {j} 系数非有限数"
                    )));
                }
                sorted.push((j, c));
            }
            sorted.sort_by_key(|&(j, _)| j);
            for (j, c) in sorted {
                match merged.last_mut() {
                    Some(last) if last.0 == j => last.1 += c,
                    _ => merged.push((j, c)),
                }
            }
            merged.retain(|&(_, c)| c != 0.0);
        }
        let k = self.cut_row.len();
        let g = merged.iter().map(|&(_, c)| c.abs()).fold(0.0_f64, f64::max);
        if g <= self.tol.zero {
            // 防御性近零割行：不入表，按原尺度 rhs 判定（契约同 simplex.rs）。
            let ok = match op {
                RelOp::Le => rhs >= -self.tol.feasibility,
                RelOp::Ge => rhs <= self.tol.feasibility,
                RelOp::Eq => rhs.abs() <= self.tol.feasibility,
            };
            if !ok {
                self.infeasible_flag = true;
            }
            self.cut_row.push(None);
            return Ok(k);
        }
        let m = self.basis.len();
        let ncs_old = self.ncols_static();
        let ncols_old_total = ncs_old + self.num_art;
        // 基列位置映射 + 割行对结构基列的系数 w。
        let mut pos_in_basis = vec![usize::MAX; ncols_old_total];
        for (i, &c) in self.basis.iter().enumerate() {
            pos_in_basis[c] = i;
        }
        let mut w = vec![0.0; m];
        for &(j, c) in &merged {
            let p = pos_in_basis[j];
            if p != usize::MAX {
                w[p] = c;
            }
        }
        // 残差 r_new = rhs − ν·x_cur（基列按当前值，非基列按所贴界）。
        let mut nu_x = 0.0;
        for p in 0..m {
            if w[p] != 0.0 {
                nu_x += w[p] * self.x_b[p];
            }
        }
        for &(j, c) in &merged {
            if pos_in_basis[j] == usize::MAX {
                nu_x += c * self.nonbasic_value(j);
            }
        }
        let r_new = rhs - nu_x;
        // 基列选择（与 simplex.rs::add_row 同款规则）。
        let slack_sign = match op {
            RelOp::Le => 1.0,
            RelOp::Ge => -1.0,
            RelOp::Eq => 0.0,
        };
        let use_slack_basic = match op {
            RelOp::Le => r_new >= -self.tol.feasibility,
            RelOp::Ge => r_new <= self.tol.feasibility,
            RelOp::Eq => false,
        };
        let art_sign = if use_slack_basic {
            0.0
        } else if op == RelOp::Eq {
            if r_new >= 0.0 {
                1.0
            } else {
                -1.0
            }
        } else if op == RelOp::Le {
            -1.0
        } else {
            1.0
        };
        let sigma_basic = if use_slack_basic {
            slack_sign
        } else {
            art_sign
        };
        // 新松弛列：静态段末尾（列号 ncs_old；人工列段整体 +1）。
        let has_slack = op != RelOp::Eq;
        let has_art = !use_slack_basic;
        if has_slack {
            let c = ncs_old;
            self.in_basis.insert(c, false);
            self.at_upper.insert(c, false);
            for b in self.basis.iter_mut() {
                if *b >= c {
                    *b += 1; // 人工列段右移一位（结构/既有松弛列号不变）
                }
            }
            self.layout.slack.push(SlackCol {
                origin: RowOrigin::Cut(k),
                sign: slack_sign,
            });
        }
        // —— CSC 重构：割行的结构系数必须进入各结构列的 CSC 段（新行 m
        //    是最大行号 ⇒ 追加在各段末尾）；新松弛列在末尾新增一段。漏掉
        //    结构段会使随后的 refactor() 从 CSC 重建出缺行的基矩阵（S3
        //    实测：setcover-300x600 误报 infeasible 的根因）。单遍重建
        //    O(nnz + merged.len())，行序天然升序。
        {
            let ncs_new = self.ncols_static();
            let mut new_start = Vec::with_capacity(ncs_new + 1);
            let mut new_rows: Vec<usize> = Vec::with_capacity(self.col_rows.len() + merged.len() + 1);
            let mut new_vals: Vec<f64> = Vec::with_capacity(self.col_vals.len() + merged.len() + 1);
            new_start.push(0);
            let mut cut_it = merged.iter().peekable();
            for j in 0..ncs_new {
                if j < ncs_old {
                    for k in self.col_start[j]..self.col_start[j + 1] {
                        // 既有段原样拷贝（新行 m 大于全部旧行，段内升序保持）。
                        new_rows.push(self.col_rows[k]);
                        new_vals.push(self.col_vals[k]);
                    }
                }
                if j < n {
                    // 结构列：追加割行系数（若有）。
                    if let Some(&&(jj, c)) = cut_it.peek() {
                        if jj == j {
                            new_rows.push(m);
                            new_vals.push(c);
                            cut_it.next();
                        }
                    }
                } else if j == ncs_old && has_slack {
                    // 新松弛列（本循环首次覆盖：旧 CSC 无其段）。
                    new_rows.push(m);
                    new_vals.push(slack_sign);
                }
                new_start.push(new_rows.len());
            }
            debug_assert_eq!(cut_it.peek(), None, "全部割项已入 CSC");
            self.col_start = new_start;
            self.col_rows = new_rows;
            self.col_vals = new_vals;
        }
        // 人工列（必要时）：新静态段之后。
        let ncols_total_new = self.ncols_static() + self.num_art + usize::from(has_art);
        self.in_basis.resize(ncols_total_new, false);
        self.at_upper.resize(ncols_total_new, false);
        let basic_col = if use_slack_basic {
            self.row_slack_of_cut(k)
        } else {
            if has_art {
                self.num_art += 1;
                self.art_info.push((m, art_sign));
                self.phase1_pending = true; // 人工列入基即欠账
            }
            Some(self.ncols_static() + self.num_art - 1)
        };
        let basic_col = basic_col.expect("use_slack_basic ⇒ has_slack");
        // 追加行记账。
        self.row_terms.push(merged);
        self.row_rhs.push(rhs);
        self.row_op.push(op);
        self.row_slack.push(if has_slack {
            Some(ncs_old)
        } else {
            None
        });
        self.row_scale.push(1.0); // 割行不缩放
        self.basis.push(basic_col);
        self.in_basis[basic_col] = true;
        self.x_b.push(r_new / sigma_basic);
        self.cut_row.push(Some(m));
        // 暂存扩容（行数 +1）。
        let m1 = m + 1;
        self.ftran_in.resize(m1, 0.0);
        self.alpha.resize(m1, 0.0);
        self.y_vec.resize(m1, 0.0);
        self.cb.resize(m1, 0.0);
        self.work.resize(m1, 0.0);
        self.factorization_dirty = true;
        Ok(k)
    }

    /// 割 k 的松弛列下标（add_row 内部：Le/Ge 割恒有）。
    fn row_slack_of_cut(&self, k: usize) -> Option<usize> {
        // 新松弛列 = 追加前的静态列数 = layout.slack 末项。
        let _ = k;
        let n = self.num_structural;
        Some(n + self.layout.slack.len() - 1)
    }

    // ---- 内部：初始基重建 ----

    /// 统一残差规则重建初始基（同 simplex.rs::rebuild_initial_basis 的
    /// 选基与 σ 规则；修正单纯形下 B = diag(σ_basic)，LU 平凡非奇异）。
    fn rebuild_initial_basis(&mut self) {
        let n = self.num_structural;
        let ncs = self.ncols_static();
        let m = self.row_terms.len();
        self.in_basis.truncate(ncs);
        self.at_upper.truncate(ncs);
        for b in self.in_basis.iter_mut() {
            *b = false;
        }
        for u in self.at_upper.iter_mut() {
            *u = false; // 冷启动非基点 = 结构变量贴下界
        }
        self.num_art = 0;
        self.art_info.clear();
        // 旧人工列连同欠账一并作废：新初始基要么无人工列，要么由紧随其后的
        // run_two_phase(true) 跑完整 Phase-1 清理（同 simplex.rs）。
        self.phase1_pending = false;
        // 残差 r_i = b̂_i − â_i·x⁰（x⁰ = 结构变量贴当前下界）。
        let mut r = vec![0.0; m];
        for i in 0..m {
            let mut ax = 0.0;
            for &(j, a) in &self.row_terms[i] {
                ax += a * self.bounds.lb[j];
            }
            r[i] = self.row_rhs[i] - ax;
        }
        // 选基 + σ（同稠密引擎）。
        self.basis.clear();
        self.etas.clear();
        let mut art_rows: Vec<usize> = Vec::new();
        for i in 0..m {
            let need_art = match self.row_op[i] {
                RelOp::Le => r[i] < -self.tol.feasibility,
                RelOp::Ge => r[i] > self.tol.feasibility,
                RelOp::Eq => true,
            };
            if need_art {
                art_rows.push(i);
                self.basis.push(usize::MAX); // 占位
            } else {
                let c = self.row_slack[i].expect("Le/Ge 行必有松弛列");
                self.basis.push(c);
                self.in_basis[c] = true;
                let sigma = self.layout.slack[c - n].sign;
                self.x_b[i] = r[i] / sigma;
            }
        }
        // 人工列按行序连续编号（静态段之后）。
        let ncols_new = ncs + art_rows.len();
        self.in_basis.resize(ncols_new, false);
        self.at_upper.resize(ncols_new, false);
        self.num_art = art_rows.len();
        for (t, &i) in art_rows.iter().enumerate() {
            let sigma = match self.row_op[i] {
                RelOp::Le => -1.0,
                RelOp::Ge => 1.0,
                RelOp::Eq => {
                    if r[i] >= 0.0 {
                        1.0
                    } else {
                        -1.0
                    }
                }
            };
            self.art_info.push((i, sigma));
            let col = ncs + t;
            self.basis[i] = col;
            self.in_basis[col] = true;
            self.x_b[i] = r[i] / sigma;
        }
        // 对角 ±1 基必非奇异；防御性兜底（契约上不可达）。
        if let Err(msg) = self.refactor() {
            self.numerical_error = Some(msg);
        }
    }

    // ---- 内部：定价 / 比值测试 / 主元 ----

    /// Dantzig/Bland 定价（与 simplex.rs::price 同构；rc 每迭代经 BTRAN
    /// 现算，天然新鲜）。人工列与固定列（ub−lb ≤ feas）永不入基。
    fn price(&mut self, phase: u8) -> Pricing {
        let ncs = self.ncols_static();
        let m = self.basis.len();
        for i in 0..m {
            self.cb[i] = self.cost_of(self.basis[i], phase);
        }
        let y = self.btran_full(&self.cb);
        self.y_vec = y;
        let mut best_j = usize::MAX;
        let mut best_delta = 0.0;
        let mut best_score = f64::NEG_INFINITY;
        for j in 0..ncs {
            if self.in_basis[j] {
                continue;
            }
            let (lb, ub) = self.col_bounds(j);
            if ub - lb <= self.tol.feasibility {
                continue; // 固定列不可移动，不入基
            }
            let mut dot = 0.0;
            for k in self.col_start[j]..self.col_start[j + 1] {
                let yr = self.y_vec[self.col_rows[k]];
                if yr != 0.0 {
                    dot += yr * self.col_vals[k];
                }
            }
            let rc = self.cost_of(j, phase) - dot;
            if !rc.is_finite() {
                return Pricing::Failed(format!("列 {j} 检验数非有限"));
            }
            let (eligible, delta) = if !self.at_upper[j] {
                (rc < -self.tol.optimality, 1.0)
            } else {
                (rc > self.tol.optimality, -1.0)
            };
            if !eligible {
                continue;
            }
            if self.bland {
                return Pricing::Candidate(j, delta); // 最小下标
            }
            let score = if delta > 0.0 { -rc } else { rc };
            if score > best_score {
                best_score = score;
                best_j = j;
                best_delta = delta;
            }
        }
        if best_j == usize::MAX {
            Pricing::Optimal
        } else {
            Pricing::Candidate(best_j, best_delta)
        }
    }

    /// 有界比值测试（与 simplex.rs::ratio_test 同构；α 为 self.alpha）。
    fn ratio_test(&mut self, q: usize, delta: f64) -> Step {
        let m = self.basis.len();
        self.cand.clear();
        let mut filtered_block = false;
        for i in 0..m {
            let a = self.alpha[i];
            if a == 0.0 {
                continue;
            }
            let alpha = delta * a; // 基变量值随 t 的变化为 x_b − t·alpha
            let (lb, ub) = self.col_bounds(self.basis[i]);
            let (t_i, hit_upper) = if alpha > 0.0 {
                if !lb.is_finite() {
                    continue;
                }
                ((self.x_b[i] - lb) / alpha, false)
            } else {
                if !ub.is_finite() {
                    continue;
                }
                ((ub - self.x_b[i]) / (-alpha), true)
            };
            let t_i = t_i.max(0.0);
            if a.abs() < self.tol.pivot {
                filtered_block = true;
                continue;
            }
            self.cand.push((t_i, i, hit_upper, a.abs(), self.basis[i]));
        }
        let (lb, ub) = self.col_bounds(q);
        if ub.is_finite() {
            self.cand.push(((ub - lb).max(0.0), usize::MAX, false, 0.0, q));
        }
        if self.cand.is_empty() {
            return if filtered_block {
                Step::FilteredBlock(format!(
                    "列 {q} 的阻挡候选全部被主元容差过滤（|p| < {}）",
                    self.tol.pivot
                ))
            } else {
                Step::NoBlock
            };
        }
        let mut t_min = f64::INFINITY;
        for c in self.cand.iter() {
            if c.0 < t_min {
                t_min = c.0;
            }
        }
        let band =
            (self.tol.feasibility * (1.0 + t_min)).min(self.tol.feasibility * RATIO_BAND_CAP);
        let mut pick: Option<&(f64, usize, bool, f64, usize)> = None;
        for c in self.cand.iter().filter(|c| c.0 <= t_min + band) {
            pick = match pick {
                None => Some(c),
                Some(best) => {
                    let better = if self.bland {
                        c.4 < best.4 // Bland：一律最小列下标
                    } else {
                        c.3 > best.3 || (c.3 == best.3 && c.4 < best.4)
                    };
                    if better {
                        Some(c)
                    } else {
                        Some(best)
                    }
                }
            };
        }
        match pick {
            Some(&(t, usize::MAX, _, _, _)) => Step::Flip { t },
            Some(&(t, row, hit_upper, _, _)) => Step::Leave { row, t, hit_upper },
            None => Step::FilteredBlock("比值测试无有效候选（内部错误）".to_string()),
        }
    }

    /// 主元执行（修正单纯形口径）：x_B 更新与表格引擎同款；基更新走 ETA
    ///（追加 (r, α)），达到 [`REFACTOR_INTERVAL`] 即整基重因子化 + 重算
    ///（失败以 Err 上抛，由 run_phase 的数值异常通道处理）。
    fn do_pivot(
        &mut self,
        q: usize,
        r: usize,
        delta: f64,
        t: f64,
        hit_upper: bool,
    ) -> std::result::Result<(), String> {
        let m = self.basis.len();
        let entering_old = if delta > 0.0 {
            self.col_bounds(q).0
        } else {
            self.col_bounds(q).1
        };
        for k in 0..m {
            if k != r {
                let f = self.alpha[k];
                if f != 0.0 {
                    self.x_b[k] -= delta * t * f;
                }
            }
        }
        self.x_b[r] = entering_old + delta * t;
        let leaving = self.basis[r];
        let mut nz: Vec<(usize, f64)> = Vec::new();
        for (i, &v) in self.alpha.iter().enumerate() {
            if i != r && v != 0.0 {
                nz.push((i, v));
            }
        }
        self.etas.push(Eta {
            r,
            pivot: self.alpha[r],
            nz,
        });
        self.basis[r] = q;
        self.in_basis[leaving] = false;
        self.in_basis[q] = true;
        self.at_upper[leaving] = hit_upper;
        self.just_refactored = false;
        // 触发用 **ETA 实际数量**（跨 resolve 生命周期）：start_call 每次调用
        // 复位本次主元计数，逐次短 resolve 永不触发周期重构 ⇒ ETA 无界增长
        //（S4 实测 setcover-300x600 上 470 条，每次 FTRAN/BTRAN 多付
        // |etas|×m 次乘加，节点级热启动被拖垮）。
        if self.eta_over_limit() {
            self.refactor()?;
        }
        Ok(())
    }

    /// 界翻转（基不变）：非基值移动 Δ = ±(ub−lb) 并按 α 传播。
    fn do_flip(&mut self, q: usize, delta: f64, t: f64) {
        let (lb, ub) = self.col_bounds(q);
        let old_v = if delta > 0.0 { lb } else { ub };
        let new_v = if delta > 0.0 { ub } else { lb };
        debug_assert!(
            (new_v - old_v - delta * t).abs() <= 1e-9 * (1.0 + (new_v - old_v).abs())
        );
        let d = new_v - old_v;
        let m = self.basis.len();
        for k in 0..m {
            let a = self.alpha[k];
            if a != 0.0 {
                self.x_b[k] -= d * a;
            }
        }
        self.at_upper[q] = delta > 0.0;
    }

    // ---- 内部：阶段循环 ----

    fn run_two_phase(&mut self, phase1_needed: bool) -> LpStatus {
        if phase1_needed && self.num_art > 0 {
            match self.run_phase(1) {
                PhaseOutcome::Failed(msg) => {
                    self.numerical_error = Some(msg);
                    self.phase1_pending = true;
                    return LpStatus::IterationLimit;
                }
                PhaseOutcome::LimitReached => {
                    self.phase1_pending = true;
                    return LpStatus::IterationLimit;
                }
                PhaseOutcome::TimedOut => {
                    // 墙钟超限：与迭代上限同款记账（人工变量未清理），但
                    // 状态码独立——W2/R8 不把墙钟超限当迭代上限处理。
                    self.phase1_pending = true;
                    return LpStatus::TimeLimit;
                }
                PhaseOutcome::Unbounded => {
                    self.numerical_error =
                        Some("Phase-1 出现无界方向（数学上不可能）".to_string());
                    self.phase1_pending = true;
                    return LpStatus::IterationLimit;
                }
                PhaseOutcome::Retry => {
                    unreachable!("run_phase 不向上层传播 Retry");
                }
                PhaseOutcome::Optimal => {
                    if self.phase1_infeasibility() > self.tol.feasibility {
                        // 人工变量带正值留基未清理：状态不可供后续 resolve
                        // 直接跑 Phase-2，记账强制下一次 resolve 重跑。
                        self.phase1_pending = true;
                        return LpStatus::Infeasible;
                    }
                    self.cleanup_after_phase1();
                    self.phase1_pending = false;
                }
            }
        }
        match self.run_phase(2) {
            PhaseOutcome::Optimal => LpStatus::Optimal,
            PhaseOutcome::Unbounded => LpStatus::Unbounded,
            PhaseOutcome::LimitReached => LpStatus::IterationLimit,
            PhaseOutcome::TimedOut => LpStatus::TimeLimit,
            PhaseOutcome::Retry => unreachable!("run_phase 不向上层传播 Retry"),
            PhaseOutcome::Failed(msg) => {
                self.numerical_error = Some(msg);
                LpStatus::IterationLimit
            }
        }
    }

    fn run_phase(&mut self, phase: u8) -> PhaseOutcome {
        let m = self.basis.len();
        let n = self.num_structural;
        let base_limit = 10_000.max(20u64.saturating_mul((m + n) as u64));
        let limit: u64 = self.iter_limit_hint.map_or(base_limit, |h| h.max(1_000));
        let bland_total: u64 = BLAND_TOTAL_FACTOR.saturating_mul((m + n) as u64);
        loop {
            if self.iter_call >= limit {
                return PhaseOutcome::LimitReached;
            }
            if let Some(deadline) = self.time_deadline {
                self.pivots_since_time_check = self.pivots_since_time_check.saturating_add(1);
                if self.pivots_since_time_check >= TIME_CHECK_INTERVAL {
                    self.pivots_since_time_check = 0;
                    if std::time::Instant::now() >= deadline {
                        return PhaseOutcome::TimedOut;
                    }
                }
            }
            #[allow(unused_assignments)]
            let mut t_nondegenerate = false;
            let (q, delta) = match self.price(phase) {
                Pricing::Failed(msg) => {
                    if let out @ PhaseOutcome::Failed(_) = self.numerical_retry_or_fail(msg) {
                        return out;
                    }
                    continue;
                }
                Pricing::Optimal => return PhaseOutcome::Optimal,
                Pricing::Candidate(q, delta) => (q, delta),
            };
            self.ftran_col(q);
            match self.ratio_test(q, delta) {
                Step::FilteredBlock(msg) => {
                    if let out @ PhaseOutcome::Failed(_) = self.numerical_retry_or_fail(msg) {
                        return out;
                    }
                    continue;
                }
                Step::NoBlock => {
                    if phase == 2 {
                        return PhaseOutcome::Unbounded;
                    }
                    return PhaseOutcome::Failed(
                        "Phase-1 出现无界方向（数学上不可能，数值异常）".to_string(),
                    );
                }
                Step::Flip { t } => {
                    t_nondegenerate = t > self.tol.feasibility;
                    self.do_flip(q, delta, t);
                    self.count_step(t);
                }
                Step::Leave { row, t, hit_upper } => {
                    if self.alpha[row].abs() < self.tol.pivot {
                        let msg = format!(
                            "主元绝对值 {} 低于容差 {}",
                            self.alpha[row].abs(),
                            self.tol.pivot
                        );
                        if let out @ PhaseOutcome::Failed(_) = self.numerical_retry_or_fail(msg) {
                            return out;
                        }
                        continue;
                    }
                    t_nondegenerate = t > self.tol.feasibility;
                    if let Err(msg) = self.do_pivot(q, row, delta, t, hit_upper) {
                        if let out @ PhaseOutcome::Failed(_) = self.numerical_retry_or_fail(msg) {
                            return out;
                        }
                        continue;
                    }
                    self.count_step(t);
                }
            }
            if self.x_b.iter().any(|v| !v.is_finite()) {
                let msg = "基变量值出现 NaN/Inf".to_string();
                if let out @ PhaseOutcome::Failed(_) = self.numerical_retry_or_fail(msg) {
                    return out;
                }
                continue;
            }
            if self.iter_call > bland_total && self.degenerate_streak > 0 {
                self.bland = true;
                self.bland_permanent = true;
            } else if self.degenerate_streak > DEGENERATE_STREAK_LIMIT {
                self.bland = true;
            }
            self.bland_release_step(t_nondegenerate);
        }
    }

    /// 数值异常的确定性重试：非刚重因子化 ⇒ 重因子化 + 重算 x_B 后返回
    /// Retry（外层 continue 重试本迭代；just_refactored 闸保证至多一次，
    /// 无循环）；刚重因子化仍异常 / 重因子化失败 ⇒ Failed。
    fn numerical_retry_or_fail(&mut self, msg: String) -> PhaseOutcome {
        if !self.just_refactored {
            return match self.refactor() {
                Ok(()) => PhaseOutcome::Retry,
                Err(m) => PhaseOutcome::Failed(m),
            };
        }
        PhaseOutcome::Failed(msg)
    }

    /// Bland 释放记账（与 simplex.rs 同款）。
    fn bland_release_step(&mut self, t_nondegenerate: bool) {
        if self.bland
            && self.bland_dwell >= BLAND_DWELL
            && t_nondegenerate
            && self.bland_rearms < BLAND_REARM_MAX
        {
            self.bland = false;
            self.bland_permanent = false;
            self.bland_rearms += 1;
            self.bland_dwell = 0;
        } else if self.bland {
            self.bland_dwell = self.bland_dwell.saturating_add(1);
        }
    }

    fn count_step(&mut self, t: f64) {
        self.iter_call = self.iter_call.saturating_add(1);
        self.iterations_total = self.iterations_total.saturating_add(1);
        if t <= self.tol.feasibility {
            self.degenerate_streak = self.degenerate_streak.saturating_add(1);
        } else {
            self.degenerate_streak = 0;
        }
    }

    /// Phase-1 最优时的人工变量总量（残基人工变量值之和，负值按 0 计）。
    fn phase1_infeasibility(&self) -> f64 {
        let ncs = self.ncols_static();
        let mut sum = 0.0;
        for i in 0..self.basis.len() {
            if self.basis[i] >= ncs {
                sum += self.x_b[i].max(0.0);
            }
        }
        sum
    }

    /// Phase-1 后处理（与 simplex.rs 同款）：残基人工变量值清零并尽力以
    /// t=0 主元驱赶出基（最小下标候选）；无候选的人工变量留基值恒 0。
    fn cleanup_after_phase1(&mut self) {
        let ncs = self.ncols_static();
        for i in 0..self.basis.len() {
            if self.basis[i] < ncs {
                continue;
            }
            self.x_b[i] = 0.0;
            let m = self.basis.len();
            let mut e = vec![0.0f64; m];
            e[i] = 1.0;
            let rho = self.btran_full(&e);
            let mut pick = usize::MAX;
            for j in 0..ncs {
                if self.in_basis[j] {
                    continue;
                }
                let mut dot = 0.0;
                for k in self.col_start[j]..self.col_start[j + 1] {
                    let r = self.col_rows[k];
                    if rho[r] != 0.0 {
                        dot += rho[r] * self.col_vals[k];
                    }
                }
                if dot.abs() >= self.tol.pivot {
                    pick = j;
                    break; // 最小下标
                }
            }
            if pick != usize::MAX {
                let delta = if self.at_upper[pick] { -1.0 } else { 1.0 };
                self.ftran_col(pick);
                if self.alpha[i].abs() >= self.tol.pivot {
                    // t=0 主元：x_B 不动（δ·0·α = 0），x_B[r] = 入基旧界值。
                    if self.do_pivot(pick, i, delta, 0.0, false).is_err() {
                        // 重因子化也失败：保持人工变量留基（值 0），后续
                        // 终检/数值通道兜底。
                        break;
                    }
                }
                self.count_step(0.0);
            }
        }
    }

    // ---- 内部：解构造 ----

    fn validate_bounds_arg(&mut self, bounds: &Bounds) -> Option<LpStatus> {
        let n = self.num_structural;
        if bounds.lb.len() != n || bounds.ub.len() != n {
            self.numerical_error = Some("bounds 维度与结构变量数不符".to_string());
            return Some(LpStatus::IterationLimit);
        }
        for j in 0..n {
            let (lb, ub) = (bounds.lb[j], bounds.ub[j]);
            if !lb.is_finite() {
                self.numerical_error = Some(format!("变量 {j} 下界非有限"));
                return Some(LpStatus::IterationLimit);
            }
            if !(ub.is_finite() || ub == f64::INFINITY) {
                self.numerical_error = Some(format!("变量 {j} 上界非法"));
                return Some(LpStatus::IterationLimit);
            }
            if lb > ub {
                return Some(LpStatus::Infeasible);
            }
        }
        None
    }

    fn start_call(&mut self) {
        self.iter_call = 0;
        self.degenerate_streak = 0;
        self.bland = false;
        self.bland_permanent = false;
        self.bland_dwell = 0;
        self.bland_rearms = 0;
        self.pivots_since_refactor = 0;
        self.pivots_since_time_check = TIME_CHECK_INTERVAL;
    }

    /// 构造 LpSolution（口径与 simplex.rs::build_solution 同款）。
    fn build_solution(&mut self, mut status: LpStatus) -> LpSolution {
        let n = self.num_structural;
        let m = self.basis.len();
        let ncs = self.ncols_static();
        let feas = self.tol.feasibility;
        let mut x = vec![0.0; n];
        for j in 0..n {
            if self.in_basis[j] {
                let mut v = f64::NAN;
                for i in 0..m {
                    if self.basis[i] == j {
                        v = self.x_b[i];
                        break;
                    }
                }
                let (lb, ub) = self.col_bounds(j);
                if v < lb && v >= lb - 10.0 * feas {
                    v = lb;
                }
                if v > ub && v <= ub + 10.0 * feas {
                    v = ub;
                }
                x[j] = v;
            } else {
                x[j] = self.nonbasic_value(j);
            }
        }
        if status == LpStatus::Optimal {
            // 最终校验：基变量界违反与残基人工变量。
            let mut ok = true;
            for i in 0..m {
                let c = self.basis[i];
                if !self.x_b[i].is_finite() {
                    ok = false;
                    break;
                }
                if c >= ncs {
                    if self.x_b[i].abs() > 10.0 * feas {
                        ok = false;
                        break;
                    }
                } else {
                    let (lb, ub) = self.col_bounds(c);
                    if self.x_b[i] < lb - 10.0 * feas || self.x_b[i] > ub + 10.0 * feas {
                        ok = false;
                        break;
                    }
                }
            }
            if !ok {
                self.numerical_error = Some("终检发现基变量越界/人工变量非零".to_string());
                status = LpStatus::IterationLimit;
            }
        }
        if status == LpStatus::Optimal {
            let mut z = 0.0;
            for j in 0..n {
                if self.c_internal[j] != 0.0 {
                    z += self.c_internal[j] * x[j];
                }
            }
            let objective = if self.sense_max {
                -z + self.obj_offset
            } else {
                z + self.obj_offset
            };
            // 对偶：y = c_B·B⁻¹（BTRAN 现算），× row_scale 还原。
            let mut cb = vec![0.0; m];
            for i in 0..m {
                cb[i] = self.cost_of(self.basis[i], 2);
            }
            let y = self.btran_full(&cb);
            let mut duals = vec![0.0; self.num_constraints];
            for ci in 0..self.num_constraints {
                if let Some(r) = self.constraint_row[ci] {
                    duals[ci] = y[r] * self.row_scale[r];
                }
            }
            LpSolution {
                status,
                objective,
                x,
                duals,
            }
        } else {
            LpSolution {
                status,
                objective: f64::NAN,
                x,
                duals: vec![0.0; self.num_constraints],
            }
        }
    }
}

/// 稀疏引擎的 trait 实现（方法名与固有方法同名：固有方法解析优先）。
impl LpEngine for SparseSimplex {
    fn solve(&mut self, bounds: &Bounds) -> LpSolution {
        self.solve(bounds)
    }
    fn resolve(&mut self, bounds: &Bounds) -> LpSolution {
        self.resolve(bounds)
    }
    fn add_row(&mut self, terms: &[(usize, f64)], op: RelOp, rhs: f64) -> Result<usize> {
        self.add_row(terms, op, rhs)
    }
    fn set_time_deadline(&mut self, deadline: Option<std::time::Instant>) {
        self.set_time_deadline(deadline)
    }
    fn time_deadline(&self) -> Option<std::time::Instant> {
        self.time_deadline()
    }
    fn set_iter_limit_hint(&mut self, hint: Option<u64>) {
        self.set_iter_limit_hint(hint)
    }
    fn iterations(&self) -> u64 {
        self.iterations()
    }
    fn basis_columns(&self) -> &[usize] {
        self.basis_columns()
    }
    fn nonbasic_columns(&self) -> Vec<NonbasicCol> {
        self.nonbasic_columns()
    }
    fn column_layout(&self) -> &ColumnLayout {
        self.column_layout()
    }
    fn row_of_constraint(&self, constraint_index: usize) -> Option<usize> {
        self.row_of_constraint(constraint_index)
    }
    fn row_of_cut(&self, cut_index: usize) -> Option<usize> {
        self.row_of_cut(cut_index)
    }
    fn num_rows(&self) -> usize {
        self.num_rows()
    }
    fn tableau_row_terms(&self, row: usize) -> (Vec<(usize, f64)>, f64) {
        self.tableau_row_terms(row)
    }
    fn snapshot(&self) -> Box<dyn LpEngine> {
        Box::new(self.clone())
    }
}

// ===========================================================================
// 单元测试
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instance::load_instance;
    use crate::model::{LinearExpr, VarType, Variable};
    use crate::simplex::Simplex as DenseSimplex;
    use std::path::PathBuf;

    // ---- 构造辅助（对齐 simplex.rs 测试风格）----

    fn var(name: &str, lb: f64, ub: f64) -> Variable {
        Variable {
            name: name.to_string(),
            var_type: VarType::Continuous,
            lb,
            ub,
        }
    }

    fn prob(
        sense: Sense,
        vars: Vec<Variable>,
        obj: Vec<(usize, f64)>,
        cons: Vec<(Vec<(usize, f64)>, RelOp, f64)>,
    ) -> Problem {
        Problem {
            name: "test".to_string(),
            sense,
            vars,
            obj: LinearExpr { terms: obj },
            obj_offset: 0.0,
            constraints: cons
                .into_iter()
                .map(|(terms, op, rhs)| crate::model::Constraint {
                    name: None,
                    expr: LinearExpr { terms },
                    op,
                    rhs,
                })
                .collect(),
        }
    }

    /// 确定性 LCG（种子固定；无哈希、无墙钟——t1 确定性契约）。
    struct Lcg(u64);
    impl Lcg {
        fn next_f64(&mut self, lo: f64, hi: f64) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let u = (self.0 >> 33) as f64 / ((1u64 << 31) as f64);
            lo + u * (hi - lo)
        }
        fn next_usize(&mut self, n: usize) -> usize {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((self.0 >> 33) as usize) % n.max(1)
        }
    }

    fn bench_instance(name: &str) -> PathBuf {
        PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/../bench/instances"))
            .join(format!("{name}.milp"))
    }

    fn milp_instance(name: &str) -> PathBuf {
        PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/instances")).join(format!(
            "{name}.milp"
        ))
    }

    /// 双引擎对账：状态恒同；Optimal 目标相对 1e-6；稀疏解自身可行性
    ///（多解顶点下 x 允许与稠密引擎不同——目标恒同是红线）。
    fn check_vs_dense(problem: &Problem, label: &str) -> (LpStatus, f64) {
        let tol = Tolerances::default();
        let bounds = Bounds::from_problem(problem);
        let dense = DenseSimplex::new(problem, tol).unwrap().solve(&bounds);
        let mut sp = SparseSimplex::new(problem, tol).unwrap();
        let sparse = sp.solve(&bounds);
        assert_eq!(
            format!("{:?}", dense.status),
            format!("{:?}", sparse.status),
            "{label}: 状态分叉（dense={:?} sparse={:?}）",
            dense.status,
            sparse.status
        );
        if dense.status == LpStatus::Optimal {
            let a = dense.objective;
            let b = sparse.objective;
            assert!(
                (a - b).abs() <= 1e-6 * (1.0 + a.abs()),
                "{label}: 目标分叉 dense={a} sparse={b}"
            );
            // 稀疏解自身可行性：界 + 行（按行 max|a| 尺度）。
            for (j, &v) in sparse.x.iter().enumerate() {
                assert!(
                    v >= problem.vars[j].lb - 1e-6 && v <= problem.vars[j].ub + 1e-6,
                    "{label}: x[{j}]={v} 越界 [{},{}]",
                    problem.vars[j].lb,
                    problem.vars[j].ub
                );
            }
            for con in &problem.constraints {
                let g = con
                    .expr
                    .terms
                    .iter()
                    .map(|&(_, a)| a.abs())
                    .fold(0.0_f64, f64::max);
                if g <= 1e-12 {
                    continue;
                }
                let lhs: f64 = con
                    .expr
                    .terms
                    .iter()
                    .map(|&(j, a)| a * sparse.x[j])
                    .sum::<f64>()
                    / g;
                let rhs = con.rhs / g;
                let ok = match con.op {
                    RelOp::Le => lhs <= rhs + 1e-6,
                    RelOp::Ge => lhs >= rhs - 1e-6,
                    RelOp::Eq => (lhs - rhs).abs() <= 1e-6,
                };
                assert!(ok, "{label}: 行违反 lhs={lhs} {} {rhs}", match con.op {
                    RelOp::Le => "<=",
                    RelOp::Ge => ">=",
                    RelOp::Eq => "==",
                });
            }
            (sparse.status, sparse.objective)
        } else {
            (sparse.status, f64::NAN)
        }
    }

    // ---- A. LU 单元 ----

    /// L·U 与 P·B 逐元素重建 + FTRAN/BTRAN 残差（多形态矩阵）。
    #[test]
    fn lu_reconstruction_and_solves() {
        let cases: Vec<Vec<Vec<(usize, f64)>>> = vec![
            // 对角。
            vec![vec![(0, 2.0)], vec![(1, -3.0)], vec![(2, 0.5)]],
            // 需要行换位的下三角状矩阵。
            vec![
                vec![(0, 0.001), (1, 3.0)],
                vec![(0, 2.0), (1, 1.0)],
            ],
            // 箭头（fill-in）。
            vec![
                vec![(0, 1.0), (2, 5.0)],
                vec![(1, 2.0), (2, 6.0)],
                vec![(0, 3.0), (1, 4.0), (2, 7.0)],
            ],
            // 一般稀疏 5×5。
            vec![
                vec![(0, 1.0), (2, 2.0)],
                vec![(0, 3.0), (1, 4.0), (3, 1.0)],
                vec![(1, -1.0), (2, 5.0), (4, 2.0)],
                vec![(2, 7.0), (3, 6.0)],
                vec![(3, 8.0), (4, 9.0)],
            ],
        ];
        for (ci, cols) in cases.iter().enumerate() {
            let m = cols.len();
            let f = factorize(cols, 1e-9, m).unwrap_or_else(|e| panic!("case {ci}: {e:?}"));
            // 稠密 B（行×列）。
            let mut b = vec![vec![0.0; m]; m];
            for (j, col) in cols.iter().enumerate() {
                for &(r, v) in col {
                    b[r][j] = v;
                }
            }
            // 稠密 L（单位下三角）× 稠密 U 的真实矩阵积重建 PB。
            let mut ld = vec![vec![0.0; m]; m];
            for j in 0..m {
                ld[j][j] = 1.0;
                for &(i, v) in &f.l_cols[j] {
                    ld[i][j] = v;
                }
            }
            let mut ud = vec![vec![0.0; m]; m];
            for j in 0..m {
                ud[j][j] = f.u_diag[j];
                for &(i, v) in &f.u_cols[j] {
                    ud[i][j] = v;
                }
            }
            for i in 0..m {
                for j in 0..m {
                    let pb: f64 = (0..m).map(|k| ld[i][k] * ud[k][j]).sum();
                    let want = b[f.perm[i]][j];
                    assert!(
                        (pb - want).abs() <= 1e-9 * (1.0 + want.abs()),
                        "case {ci}: (LU)[{i}][{j}]={pb} != PB={want}"
                    );
                }
            }
            // FTRAN：B α = e_j（逐列即 B⁻¹）。
            for j in 0..m {
                let mut a = vec![0.0; m];
                a[j] = 1.0;
                let mut alpha = vec![0.0; m];
                f.ftran_dense(&a, &mut alpha);
                for i in 0..m {
                    let dot: f64 = (0..m).map(|k| b[i][k] * alpha[k]).sum();
                    assert!(
                        (dot - a[i]).abs() <= 1e-8,
                        "case {ci}: FTRAN Bα≠e_{j}（第 {i} 行）"
                    );
                }
            }
            // BTRAN：Bᵀ y = e_j。
            for j in 0..m {
                let mut c = vec![0.0; m];
                c[j] = 1.0;
                let y = f.btran(&c);
                for i in 0..m {
                    let dot: f64 = (0..m).map(|k| b[k][i] * y[k]).sum();
                    assert!(
                        (dot - c[i]).abs() <= 1e-8,
                        "case {ci}: BTRAN Bᵀy≠e_{j}（第 {i} 列）"
                    );
                }
            }
        }
    }

    /// 奇异基（零列）被拒绝。
    #[test]
    fn lu_singular_rejected() {
        let cols = vec![vec![(0, 1.0)], vec![]]; // 第 2 列全零
        assert!(factorize(&cols, 1e-9, 2).is_err());
    }

    /// 稠密随机矩阵（LCG 固定种子）20×20：因子化与双向求解残差。
    #[test]
    fn lu_dense_random_matrices() {
        let mut rng = Lcg(0xC0FFEE);
        for case in 0..8 {
            let m = 20usize;
            let mut cols: Vec<Vec<(usize, f64)>> = (0..m)
                .map(|j| {
                    let mut col: Vec<(usize, f64)> = Vec::new();
                    for i in 0..m {
                        let v = rng.next_f64(-1.0, 1.0);
                        if v.abs() > 0.75 {
                            col.push((i, if v.abs() < 0.05 { 1.0 } else { v }));
                        }
                    }
                    if col.is_empty() {
                        col.push((j, 1.0)); // 防奇异
                    }
                    let _ = j;
                    col
                })
                .collect();
            // 对角加强（可解性）。
            for j in 0..m {
                cols[j].push((j, 2.0));
            }
            let f = factorize(&cols, 1e-9, m).unwrap_or_else(|e| panic!("case {case}: {e:?}"));
            let mut b = vec![vec![0.0; m]; m];
            for (j, col) in cols.iter().enumerate() {
                for &(r, v) in col {
                    b[r][j] += v;
                }
            }
            let mut a = vec![0.0; m];
            for i in 0..m {
                a[i] = rng.next_f64(-1.0, 1.0);
            }
            let mut alpha = vec![0.0; m];
            f.ftran_dense(&a, &mut alpha);
            for i in 0..m {
                let dot: f64 = (0..m).map(|k| b[i][k] * alpha[k]).sum();
                assert!(
                    (dot - a[i]).abs() <= 1e-6 * (1.0 + a[i].abs()),
                    "case {case}: FTRAN 残差 @ {i}"
                );
            }
            let y = f.btran(&a);
            for i in 0..m {
                let dot: f64 = (0..m).map(|k| b[k][i] * y[k]).sum();
                assert!(
                    (dot - a[i]).abs() <= 1e-6 * (1.0 + a[i].abs()),
                    "case {case}: BTRAN 残差 @ {i}"
                );
            }
        }
    }

    // ---- B. LP 语义（小算例 + 双引擎对账）----

    #[test]
    fn lp_bound_flip_optimum_at_ub() {
        // max x, 0≤x≤5 ⇒ 最优 5（自身界翻转路径）。
        let p = prob(
            Sense::Maximize,
            vec![var("x", 0.0, 5.0)],
            vec![(0, 1.0)],
            vec![(vec![(0, 1.0)], RelOp::Le, 100.0)],
        );
        let (st, obj) = check_vs_dense(&p, "bound-flip");
        assert_eq!(format!("{st:?}"), "Optimal");
        assert!((obj - 5.0).abs() < 1e-7);
    }

    #[test]
    fn lp_equality_row_needs_phase1() {
        // max x+y s.t. x+y == 3, 0≤x,y≤2 ⇒ 3（人工列 + Phase-1）。
        let p = prob(
            Sense::Maximize,
            vec![var("x", 0.0, 2.0), var("y", 0.0, 2.0)],
            vec![(0, 1.0), (1, 1.0)],
            vec![(vec![(0, 1.0), (1, 1.0)], RelOp::Eq, 3.0)],
        );
        let (st, obj) = check_vs_dense(&p, "eq-row");
        assert_eq!(format!("{st:?}"), "Optimal");
        assert!((obj - 3.0).abs() < 1e-7);
    }

    #[test]
    fn lp_ge_positive_rhs_needs_artificial() {
        // min x s.t. x ≥ 3（x⁰=0 违反 ⇒ 人工列 σ=+1）⇒ 3。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0)],
            vec![(0, 1.0)],
            vec![(vec![(0, 1.0)], RelOp::Ge, 3.0)],
        );
        let (st, obj) = check_vs_dense(&p, "ge-positive");
        assert_eq!(format!("{st:?}"), "Optimal");
        assert!((obj - 3.0).abs() < 1e-7);
    }

    #[test]
    fn lp_le_violated_rhs_needs_artificial() {
        // min x+y s.t. x−y ≤ −5, x,y∈[0,10]：x⁰=(0,0) 违反（0 ≤ −5 假）
        // ⇒ 人工列 σ=−1；最优 x=0,y=5 ⇒ 5。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0), var("y", 0.0, 10.0)],
            vec![(0, 1.0), (1, 1.0)],
            vec![(vec![(0, 1.0), (1, -1.0)], RelOp::Le, -5.0)],
        );
        let (st, obj) = check_vs_dense(&p, "le-violated");
        assert_eq!(format!("{st:?}"), "Optimal");
        assert!((obj - 5.0).abs() < 1e-7);
    }

    #[test]
    fn lp_infeasible_detected_by_phase1() {
        // x ≤ 1 且 x ≥ 2。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0)],
            vec![(0, 1.0)],
            vec![
                (vec![(0, 1.0)], RelOp::Le, 1.0),
                (vec![(0, 1.0)], RelOp::Ge, 2.0),
            ],
        );
        let (st, _) = check_vs_dense(&p, "infeasible");
        assert_eq!(format!("{st:?}"), "Infeasible");
    }

    #[test]
    fn lp_unbounded_detected() {
        // max x, x ≥ 0 无上界（NoBlock ⇒ Unbounded，F7）。
        let p = prob(
            Sense::Maximize,
            vec![var("x", 0.0, f64::INFINITY)],
            vec![(0, 1.0)],
            vec![(vec![(0, 1.0)], RelOp::Ge, 0.0)],
        );
        let (st, _) = check_vs_dense(&p, "unbounded");
        assert_eq!(format!("{st:?}"), "Unbounded");
    }

    #[test]
    fn lp_near_zero_rows_excluded() {
        // 近零行满足 ⇒ 剔除照常解；近零行违反 ⇒ Infeasible。
        let ok = prob(
            Sense::Maximize,
            vec![var("x", 0.0, 5.0), var("y", 0.0, 5.0)],
            vec![(0, 1.0), (1, 1.0)],
            vec![
                (vec![(0, 1.0), (1, 1.0)], RelOp::Le, 4.0),
                (vec![(1, 0.0)], RelOp::Le, 5.0), // 近零行（系数被 LinearExpr
                // 语义下不会出现——手动构造归并后的 terms 才能造出；此处
                // 用显式零系数行近似：terms 均为 0 系数）
            ],
        );
        let _ = ok;
        // 直接构造近零行（绕过 add_term 归并）。
        let mut p = prob(
            Sense::Maximize,
            vec![var("x", 0.0, 5.0), var("y", 0.0, 5.0)],
            vec![(0, 1.0), (1, 1.0)],
            vec![(vec![(0, 1.0), (1, 1.0)], RelOp::Le, 4.0)],
        );
        p.constraints.push(crate::model::Constraint {
            name: None,
            expr: LinearExpr {
                terms: vec![(1, 1e-13)], // max|a| ≤ tol.zero=1e-12
            },
            op: RelOp::Le,
            rhs: 1.0,
        });
        let (st, obj) = check_vs_dense(&p, "near-zero-ok");
        assert_eq!(format!("{st:?}"), "Optimal");
        assert!((obj - 4.0).abs() < 1e-6);
        // 违反变体：rhs = −1。
        let mut p2 = p.clone();
        p2.constraints.pop();
        p2.constraints.push(crate::model::Constraint {
            name: None,
            expr: LinearExpr {
                terms: vec![(1, 1e-13)],
            },
            op: RelOp::Le,
            rhs: -1.0,
        });
        let (st2, _) = check_vs_dense(&p2, "near-zero-bad");
        assert_eq!(format!("{st2:?}"), "Infeasible");
    }

    #[test]
    fn lp_degenerate_vertex() {
        // 退化顶点：x≤1, y≤1, x+y≤2, max x+y（多行同点并列）。
        let p = prob(
            Sense::Maximize,
            vec![var("x", 0.0, 5.0), var("y", 0.0, 5.0)],
            vec![(0, 1.0), (1, 1.0)],
            vec![
                (vec![(0, 1.0)], RelOp::Le, 1.0),
                (vec![(1, 1.0)], RelOp::Le, 1.0),
                (vec![(0, 1.0), (1, 1.0)], RelOp::Le, 2.0),
            ],
        );
        let (st, obj) = check_vs_dense(&p, "degenerate");
        assert_eq!(format!("{st:?}"), "Optimal");
        assert!((obj - 2.0).abs() < 1e-7);
    }

    #[test]
    fn lp_minimize_with_offset() {
        // min 2x+y(+10), x+y≥4 ⇒ 用满 y（系数 1）= 4+10（obj_offset 还原）。
        let mut p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0), var("y", 0.0, 10.0)],
            vec![(0, 2.0), (1, 1.0)],
            vec![(vec![(0, 1.0), (1, 1.0)], RelOp::Ge, 4.0)],
        );
        p.obj_offset = 10.0;
        let (st, obj) = check_vs_dense(&p, "min-offset");
        assert_eq!(format!("{st:?}"), "Optimal");
        assert!((obj - 14.0).abs() < 1e-7);
    }

    /// 生成 setcover 型 LP（固定种子；确定性）。
    fn setcover_lp(n: usize, m: usize) -> Problem {
        let mut rng = Lcg(0xDEADBEEF);
        let vars: Vec<Variable> = (0..n).map(|j| var(&format!("x{j}"), 0.0, 1.0)).collect();
        let obj: Vec<(usize, f64)> = (0..n).map(|j| (j, 1.0 + (j % 7) as f64)).collect();
        let mut cons: Vec<(Vec<(usize, f64)>, RelOp, f64)> = Vec::new();
        for _ in 0..m {
            let mut js: Vec<usize> = (0..8).map(|_| rng.next_usize(n)).collect();
            js.sort_unstable();
            js.dedup();
            let terms: Vec<(usize, f64)> = js.iter().map(|&j| (j, 1.0)).collect();
            cons.push((terms, RelOp::Ge, 1.0));
        }
        prob(Sense::Minimize, vars, obj, cons)
    }

    #[test]
    fn iter_limit_hint_truncates_honestly() {
        // F4 接线点证据：hint 钳到下限 1000（与 simplex.rs 同款
        // h.max(1000)）；生成 setcover 型 LP（300x200，固定种子）全程
        // 1427 迭代 > 1000 ⇒ hint=1 ⇒ IterationLimit + objective NaN
        //（诚实截断，不冒充 Optimal）。
        let p = setcover_lp(300, 200);
        let mut sp = SparseSimplex::new(&p, Tolerances::default()).unwrap();
        sp.set_iter_limit_hint(Some(1));
        let sol = sp.solve(&Bounds::from_problem(&p));
        assert_eq!(format!("{:?}", sol.status), "IterationLimit");
        assert!(sol.objective.is_nan());
        // 对照：无 hint（默认上限 max(10k, 20(m+n))）完整最优且与稠密
        // 引擎目标一致。
        check_vs_dense(&p, "setcover-gen-300x200");
    }

    #[test]
    fn time_deadline_truncates_and_none_is_clean() {
        // 中等 LP（150 变量 × 100 Ge 行，同 simplex.rs W2 测试量级）。
        let n = 150usize;
        let vars: Vec<Variable> = (0..n).map(|j| var(&format!("x{j}"), 0.0, 1.0)).collect();
        let obj: Vec<(usize, f64)> = (0..n).map(|j| (j, (7 * j % 19 + 1) as f64)).collect();
        let mut cons: Vec<(Vec<(usize, f64)>, RelOp, f64)> = Vec::new();
        for i in 0..100usize {
            let mut terms: Vec<(usize, f64)> = (0..6)
                .map(|k| ((i * 7 + k * 23) % n, (3 * k + 1) as f64))
                .collect();
            terms.sort_by_key(|&(j, _)| j);
            cons.push((terms, RelOp::Ge, 3.0));
        }
        let p = prob(Sense::Minimize, vars, obj, cons);
        let bounds = Bounds::from_problem(&p);
        // 过期 deadline ⇒ TimeLimit（首迭代即检查，条件读钟）。
        let mut sp = SparseSimplex::new(&p, Tolerances::default()).unwrap();
        sp.set_time_deadline(Some(std::time::Instant::now()));
        let sol = sp.solve(&bounds);
        assert_eq!(format!("{:?}", sol.status), "TimeLimit");
        // 未注入 ⇒ 正常最优（W2：None 态零墙钟读点）。
        let mut sp2 = SparseSimplex::new(&p, Tolerances::default()).unwrap();
        let sol2 = sp2.solve(&bounds);
        assert_eq!(format!("{:?}", sol2.status), "Optimal");
        assert!(sp2.time_deadline().is_none());
    }

    #[test]
    fn deterministic_repeat_solves() {
        let p = load_instance(&milp_instance("knapsack-10")).unwrap();
        let bounds = Bounds::from_problem(&p);
        let mut a = SparseSimplex::new(&p, Tolerances::default()).unwrap();
        let mut b = SparseSimplex::new(&p, Tolerances::default()).unwrap();
        let sa = a.solve(&bounds);
        let sb = b.solve(&bounds);
        assert_eq!(format!("{:?}", sa.status), format!("{:?}", sb.status));
        assert_eq!(sa.objective.to_bits(), sb.objective.to_bits());
        assert_eq!(a.iterations(), b.iterations());
    }

    /// 克隆即精确状态（W6 快照纪律等价物）。
    #[test]
    fn clone_is_exact_state() {
        let p = load_instance(&milp_instance("hard-knapsack-12")).unwrap();
        let bounds = Bounds::from_problem(&p);
        let mut a = SparseSimplex::new(&p, Tolerances::default()).unwrap();
        let sa = a.solve(&bounds);
        let mut c = a.clone();
        let sc = c.solve(&bounds);
        assert_eq!(sa.objective.to_bits(), sc.objective.to_bits());
        assert_eq!(a.iterations() * 2, c.iterations());
    }

    // ---- C. 访问器等价性（F3：cuts.rs 消费面）----

    /// 初始基下与稠密引擎逐项等价（选基规则相同 ⇒ B=diag(σ) ⇒ 表格行
    /// 数学恒同）。
    #[test]
    fn accessors_match_dense_at_initial_basis() {
        let p = load_instance(&milp_instance("milp-mixed-5")).unwrap();
        let tol = Tolerances::default();
        let dense = DenseSimplex::new(&p, tol).unwrap();
        let sp = SparseSimplex::new(&p, tol).unwrap();
        assert_eq!(sp.num_rows(), dense.num_rows());
        assert_eq!(sp.basis_columns(), dense.basis_columns());
        assert_eq!(sp.row_of_constraint(0), dense.row_of_constraint(0));
        // 列布局：结构数 + 各松弛列 origin/sign。
        let dl = dense.column_layout();
        let sl = sp.column_layout();
        assert_eq!(dl.num_structural, sl.num_structural);
        assert_eq!(dl.slack.len(), sl.slack.len());
        for (a, b) in dl.slack.iter().zip(sl.slack.iter()) {
            assert_eq!(format!("{:?}", a.origin), format!("{:?}", b.origin));
            assert_eq!(a.sign, b.sign);
        }
        let dn = dense.nonbasic_columns();
        let sn = sp.nonbasic_columns();
        assert_eq!(dn.len(), sn.len());
        for (a, b) in dn.iter().zip(sn.iter()) {
            assert_eq!((a.col, a.at_upper), (b.col, b.at_upper));
        }
        for row in 0..sp.num_rows() {
            let (d_terms, d_rhs) = dense.tableau_row_terms(row);
            let (s_terms, s_rhs) = sp.tableau_row_terms(row);
            assert_eq!(d_terms.len(), s_terms.len(), "行 {row} 项数");
            for ((dj, dv), (sj, sv)) in d_terms.iter().zip(s_terms.iter()) {
                assert_eq!(dj, sj, "行 {row} 列序");
                assert!(
                    (dv - sv).abs() <= 1e-12 * (1.0 + dv.abs()),
                    "行 {row} 列 {dj}: dense={dv} sparse={sv}"
                );
            }
            assert!((d_rhs - s_rhs).abs() <= 1e-12 * (1.0 + d_rhs.abs()));
        }
    }

    /// 解后自洽性：每行表格行在基列段 = 单位阵（e_iᵀB⁻¹·B = e_i）——
    /// BTRAN + 行收集实现的数学身份验证（不依赖两引擎基相同）。
    #[test]
    fn tableau_rows_unit_on_basis_after_solve() {
        let p = load_instance(&milp_instance("knapsack-10")).unwrap();
        let mut sp = SparseSimplex::new(&p, Tolerances::default()).unwrap();
        let sol = sp.solve(&Bounds::from_problem(&p));
        assert_eq!(format!("{:?}", sol.status), "Optimal");
        let basis: Vec<usize> = sp.basis_columns().to_vec();
        for row in 0..sp.num_rows() {
            let (terms, _) = sp.tableau_row_terms(row);
            let mut by_col = std::collections::BTreeMap::new();
            for (j, v) in terms {
                by_col.insert(j, v);
            }
            for (i, &bcol) in basis.iter().enumerate() {
                let got = by_col.get(&bcol).copied().unwrap_or(0.0);
                let want = if i == row { 1.0 } else { 0.0 };
                assert!(
                    (got - want).abs() <= 1e-7,
                    "行 {row} 基列 {bcol}（位置 {i}）：{got} ≠ {want}"
                );
            }
        }
    }

    // ---- D. 实例对账（milp/instances 全 LP 松弛）----

    #[test]
    fn parity_on_milp_instances() {
        for name in [
            "knapsack-10",
            "lp-prod-3",
            "milp-mixed-5",
            "milp-wolsey-2",
            "assign-4x4",
            "hard-knapsack-12",
        ] {
            let p = load_instance(&milp_instance(name)).unwrap();
            check_vs_dense(&p, name);
        }
    }

    /// 确定性随机 LP 双引擎对账（10 种子 × 混合算子/等式行/有限上界）。
    #[test]
    fn parity_on_generated_lps() {
        for seed in 1..=10u64 {
            let mut rng = Lcg(seed.wrapping_mul(0x9E3779B97F4A7C15));
            let n = 12 + rng.next_usize(9);
            let m = 6 + rng.next_usize(7);
            let vars: Vec<Variable> = (0..n)
                .map(|j| var(&format!("x{j}"), 0.0, rng.next_f64(2.0, 8.0)))
                .collect();
            let obj: Vec<(usize, f64)> = (0..n)
                .map(|j| (j, (rng.next_f64(-5.0, 5.0) * 10.0).round() / 10.0))
                .collect();
            // 行按可行点 x0 生成（保证大多可行；Eq 行_rhs 精确 a·x0）。
            let x0: Vec<f64> = (0..n).map(|j| vars[j].lb + (vars[j].ub - vars[j].lb) * 0.3).collect();
            let mut cons: Vec<(Vec<(usize, f64)>, RelOp, f64)> = Vec::new();
            for _ in 0..m {
                let k = 2 + rng.next_usize(4);
                let mut js: Vec<usize> = (0..k).map(|_| rng.next_usize(n)).collect();
                js.sort_unstable();
                js.dedup();
                let terms: Vec<(usize, f64)> = js
                    .iter()
                    .map(|&j| {
                        let v = (rng.next_f64(-3.0, 3.0) * 10.0).round() / 10.0;
                        (j, if v == 0.0 { 1.0 } else { v })
                    })
                    .collect();
                let ax: f64 = terms.iter().map(|&(j, a)| a * x0[j]).sum();
                let op = match rng.next_usize(3) {
                    0 => RelOp::Le,
                    1 => RelOp::Ge,
                    _ => RelOp::Eq,
                };
                let rhs = match op {
                    RelOp::Le => ax + rng.next_f64(0.0, 2.0),
                    RelOp::Ge => ax - rng.next_f64(0.0, 2.0),
                    RelOp::Eq => ax,
                };
                cons.push((terms, op, (rhs * 100.0).round() / 100.0));
            }
            let sense = if seed % 2 == 0 {
                Sense::Maximize
            } else {
                Sense::Minimize
            };
            let p = prob(sense, vars, obj, cons);
            check_vs_dense(&p, &format!("gen-seed-{seed}"));
        }
    }

    // ---- E. bench 实例验收（F2 钉值：lp-wide=112.66861478047609、
    //      lp-sparse=5893.705757——各算例各钉各值）----

    /// lp-wide-200x300：双引擎对账 + 钉值（默认套件内：稠密 ~93 迭代）。
    #[test]
    fn lp_wide_200x300_parity_and_pin() {
        let p = load_instance(&bench_instance("lp-wide-200x300")).unwrap();
        let (st, obj) = check_vs_dense(&p, "lp-wide-200x300");
        assert_eq!(format!("{st:?}"), "Optimal");
        assert!(
            (obj - 112.66861478047609).abs() <= 1e-6 * (1.0 + 112.668),
            "lp-wide 钉值 112.66861478047609，实得 {obj}"
        );
    }

    /// lp-sparse-4000x800：双引擎对账 + 钉值 5893.705757 + 性能验收
    ///（稀疏 3 rep 中位 ≥2× 于稠密 1 rep 且 ≤15s；重负载标 #[ignore]，
    /// 门禁/验收以 `cargo test --release --lib simplex_sparse -- --ignored
    /// --test-threads=1` 显式跑）。
    #[test]
    #[ignore]
    fn lp_sparse_4000x800_parity_pin_and_speed() {
        let p = load_instance(&bench_instance("lp-sparse-4000x800")).unwrap();
        let tol = Tolerances::default();
        let bounds = Bounds::from_problem(&p);
        // 稠密基线 1 rep（~7.8k 迭代，历史 gate 记录 39.609s / python_bench
        // 34.757s；standalone 口径本机现测）。
        let t0 = std::time::Instant::now();
        let dense = DenseSimplex::new(&p, tol).unwrap().solve(&bounds);
        let dense_secs = t0.elapsed().as_secs_f64();
        assert_eq!(format!("{:?}", dense.status), "Optimal", "稠密基线须最优");
        // 稀疏 3 rep。
        let mut sparse_secs: Vec<f64> = Vec::new();
        let mut sparse_obj = f64::NAN;
        let mut sparse_iters = 0;
        for rep in 0..3 {
            let mut sp = SparseSimplex::new(&p, tol).unwrap();
            let t0 = std::time::Instant::now();
            let sol = sp.solve(&bounds);
            sparse_secs.push(t0.elapsed().as_secs_f64());
            assert_eq!(format!("{:?}", sol.status), "Optimal", "rep {rep}");
            sparse_obj = sol.objective;
            sparse_iters = sp.iterations();
        }
        sparse_secs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med = sparse_secs[1];
        eprintln!(
            "lp-sparse-4000x800: dense={dense_secs:.3}s | sparse reps={:.3}/{:.3}/{:.3}s (median {med:.3}s), iters={sparse_iters}, obj={sparse_obj:.6}",
            sparse_secs[0], sparse_secs[1], sparse_secs[2]
        );
        // 对账与钉值。
        assert!(
            (dense.objective - sparse_obj).abs() <= 1e-6 * (1.0 + dense.objective.abs()),
            "目标分叉 dense={} sparse={}",
            dense.objective,
            sparse_obj
        );
        assert!(
            (sparse_obj - 5893.705757).abs() <= 1e-6 * (1.0 + 5893.7),
            "lp-sparse 钉值 5893.705757，实得 {sparse_obj}"
        );
        // 性能验收：≥2× 且 ≤15s。
        assert!(
            med * 2.0 <= dense_secs,
            "稀疏中位 {med:.3}s 未达稠密 {dense_secs:.3}s 的一半"
        );
        assert!(med <= 15.0, "稀疏中位 {med:.3}s 超 15s 验收线");
    }

    // ---- F. S3 集成面：resolve 热启动 / add_row 割 / 快照作用域 ----

    /// 界收紧（分支形态）下的热启动 resolve：目标与冷解、与稠密引擎
    /// resolve 三方一致。
    #[test]
    fn resolve_warm_matches_cold_and_dense_after_tightening() {
        let p = load_instance(&milp_instance("knapsack-10")).unwrap();
        let tol = Tolerances::default();
        let mut b = Bounds::from_problem(&p);
        // 分支形态 1：x0 收紧上界。
        b.ub[0] = 0.0;
        let warm = {
            let mut sp = SparseSimplex::new(&p, tol).unwrap();
            let _ = sp.solve(&Bounds::from_problem(&p));
            sp.resolve(&b)
        };
        let cold = {
            let mut sp = SparseSimplex::new(&p, tol).unwrap();
            sp.solve(&b)
        };
        let dense = {
            let mut d = DenseSimplex::new(&p, tol).unwrap();
            let _ = d.solve(&Bounds::from_problem(&p));
            d.resolve(&b)
        };
        assert_eq!(format!("{:?}", warm.status), format!("{:?}", cold.status));
        assert_eq!(format!("{:?}", warm.status), format!("{:?}", dense.status));
        if warm.status == LpStatus::Optimal {
            assert!(
                (warm.objective - cold.objective).abs() <= 1e-7 * (1.0 + cold.objective.abs()),
                "warm={} cold={}",
                warm.objective,
                cold.objective
            );
            assert!(
                (warm.objective - dense.objective).abs() <= 1e-7 * (1.0 + dense.objective.abs()),
                "warm={} dense={}",
                warm.objective,
                dense.objective
            );
        }
        // 分支形态 2（更深树路径）：连续收紧 x1..x3。
        let mut b2 = Bounds::from_problem(&p);
        for j in 0..3 {
            b2.lb[j] = 1.0;
        }
        let mut sp = SparseSimplex::new(&p, tol).unwrap();
        let _ = sp.solve(&Bounds::from_problem(&p));
        let sol2 = sp.resolve(&b2);
        let mut d = DenseSimplex::new(&p, tol).unwrap();
        let _ = d.solve(&Bounds::from_problem(&p));
        let d2 = d.resolve(&b2);
        assert_eq!(format!("{:?}", sol2.status), format!("{:?}", d2.status));
        // Infeasible（如 x0..x2=1 超容量）时 objective 为 NaN：只对账状态；
        // Optimal 时对账目标（多解顶点下 x 允许不同——目标恒同是红线）。
        if sol2.status == LpStatus::Optimal {
            assert!(
                (sol2.objective - d2.objective).abs() <= 1e-7 * (1.0 + d2.objective.abs()),
                "{} vs {}",
                sol2.objective,
                d2.objective
            );
        }
    }

    /// dive 形态的连续固定（lb=ub）resolve 序列：终态与目标与稠密引擎一致。
    #[test]
    fn dive_like_fixing_sequence_matches_dense() {
        let p = load_instance(&milp_instance("knapsack-10")).unwrap();
        let tol = Tolerances::default();
        let root = Bounds::from_problem(&p);
        let mut sp = SparseSimplex::new(&p, tol).unwrap();
        let mut d = DenseSimplex::new(&p, tol).unwrap();
        let _ = d.solve(&root);
        let mut b = root.clone();
        let base = sp.solve(&root);
        let root_x = base.x.clone();
        for j in 0..5 {
            let v = root_x[j].round().max(b.lb[j]).min(b.ub[j]);
            b.lb[j] = v;
            b.ub[j] = v;
            let rs = sp.resolve(&b);
            let rd = d.resolve(&b);
            assert_eq!(
                format!("{:?}", rs.status),
                format!("{:?}", rd.status),
                "step {j}: 状态分叉"
            );
            if rs.status == LpStatus::Optimal {
                assert!(
                    (rs.objective - rd.objective).abs()
                        <= 1e-7 * (1.0 + rd.objective.abs()),
                    "step {j}: {} vs {}",
                    rs.objective,
                    rd.objective
                );
            }
        }
    }

    /// add_row + resolve（割后重解）：与稠密引擎同序列同目标；行记账
    ///（row_of_cut/num_rows/布局）正确。
    #[test]
    fn add_row_cut_then_resolve_matches_dense() {
        let p = load_instance(&milp_instance("knapsack-10")).unwrap();
        let tol = Tolerances::default();
        let bounds = Bounds::from_problem(&p);
        // 手工合法割（binary 有效）：x0 + x1 ≤ 1。
        let cut_terms = vec![(0usize, 1.0), (1usize, 1.0)];
        let mut sp = SparseSimplex::new(&p, tol).unwrap();
        let mut d = DenseSimplex::new(&p, tol).unwrap();
        let s0 = sp.solve(&bounds);
        let d0 = d.solve(&bounds);
        assert_eq!(s0.status, d0.status);
        let k1 = sp.add_row(&cut_terms, RelOp::Le, 1.0).unwrap();
        let k1d = d.add_row(&cut_terms, RelOp::Le, 1.0).unwrap();
        assert_eq!(k1, k1d, "割序号一致");
        assert_eq!(sp.row_of_cut(k1), Some(sp.num_rows() - 1));
        assert_eq!(sp.num_rows(), d.num_rows());
        let s1 = sp.resolve(&bounds);
        let d1 = d.resolve(&bounds);
        assert_eq!(format!("{:?}", s1.status), format!("{:?}", d1.status));
        assert!(
            (s1.objective - d1.objective).abs() <= 1e-7 * (1.0 + d1.objective.abs()),
            "{} vs {}",
            s1.objective,
            d1.objective
        );
        // 第二条割（x2 + x3 ≤ 1）后重解。
        let cut2 = vec![(2usize, 1.0), (3usize, 1.0)];
        sp.add_row(&cut2, RelOp::Le, 1.0).unwrap();
        d.add_row(&cut2, RelOp::Le, 1.0).unwrap();
        let s2 = sp.resolve(&bounds);
        let d2 = d.resolve(&bounds);
        assert_eq!(format!("{:?}", s2.status), format!("{:?}", d2.status));
        assert!(
            (s2.objective - d2.objective).abs() <= 1e-7 * (1.0 + d2.objective.abs()),
            "{} vs {}",
            s2.objective,
            d2.objective
        );
        assert_eq!(sp.num_rows(), d.num_rows());
    }

    /// 违反当前解的割（人工列路径：phase1_pending ⇒ resolve 整基重建 +
    /// Phase-1）后仍得正确目标；与稠密引擎一致。
    #[test]
    fn add_row_violated_cut_artificial_path_matches_dense() {
        let p = load_instance(&milp_instance("knapsack-10")).unwrap();
        let tol = Tolerances::default();
        let bounds = Bounds::from_problem(&p);
        let mut sp = SparseSimplex::new(&p, tol).unwrap();
        let mut d = DenseSimplex::new(&p, tol).unwrap();
        let _ = sp.solve(&bounds);
        let _ = d.solve(&bounds);
        // 根最优处 x[0] 多为 1：割 x0 ≤ 0 违反当前解 ⇒ r_new < 0 ⇒ 人工列。
        sp.add_row(&[(0usize, 1.0)], RelOp::Le, 0.0).unwrap();
        d.add_row(&[(0usize, 1.0)], RelOp::Le, 0.0).unwrap();
        let s1 = sp.resolve(&bounds);
        let d1 = d.resolve(&bounds);
        assert_eq!(format!("{:?}", s1.status), format!("{:?}", d1.status));
        assert!(
            (s1.objective - d1.objective).abs() <= 1e-7 * (1.0 + d1.objective.abs()),
            "{} vs {}",
            s1.objective,
            d1.objective
        );
    }

    /// W6 快照作用域协议：snapshot 恢复后与施加割前的解逐值一致（Clone 等价）。
    #[test]
    fn snapshot_restore_scopes_local_cuts() {
        let p = load_instance(&milp_instance("knapsack-10")).unwrap();
        let tol = Tolerances::default();
        let bounds = Bounds::from_problem(&p);
        let mut sp = SparseSimplex::new(&p, tol).unwrap();
        let base = sp.solve(&bounds);
        let snap = sp.clone();
        sp.add_row(&[(0usize, 1.0)], RelOp::Le, 0.0).unwrap();
        let _ = sp.resolve(&bounds);
        // 恢复快照。
        let mut restored = snap;
        let r = restored.resolve(&bounds);
        assert_eq!(format!("{:?}", r.status), "Optimal");
        assert!(
            (r.objective - base.objective).abs() <= 1e-9 * (1.0 + base.objective.abs()),
            "恢复后 {} vs 施加前 {}",
            r.objective,
            base.objective
        );
    }

    /// setcover-100x150（revised 全路径，B&B 口径在 solver.rs）：根割循环
    /// 全链路（gomory/覆盖割 + add_row + resolve + R8 hint）下的根 LP 目标
    /// 与稠密引擎一致（276 为已知最优，此处对账根 LP 后的树由 solver 层钉）。
    #[test]
    fn root_cut_loop_matches_dense_on_setcover_100x150() {
        let p = load_instance(&bench_instance("setcover-100x150")).unwrap();
        let tol = Tolerances::default();
        let bounds = Bounds::from_problem(&p);
        let int_vars = p.integer_vars();
        let mut sp = SparseSimplex::new(&p, tol).unwrap();
        let mut d = DenseSimplex::new(&p, tol).unwrap();
        let _ = sp.solve(&bounds);
        let _ = d.solve(&bounds);
        // 契约注：两引擎的最优基不同 ⇒ GMI 分离的割集不同 ⇒ 割后根 LP
        // 目标允许不同（都是合法有效不等式子集）；此处钉『双侧 Optimal +
        // 各自目标单调收紧（≥ 各自未割根 LP）』，全树终值恒同由
        // solver::tests 的 revised 全路径对账钉死。
        let root_obj_sp = sp.resolve(&bounds).objective;
        let root_obj_d = d.resolve(&bounds).objective;
        let (a1, st1, _) = crate::cuts::cut_loop_recorded(
            &mut sp,
            &bounds,
            &p,
            &int_vars,
            crate::cuts::CutOptions::default(),
            &tol,
            || SparseSimplex::new(&p, tol),
        )
        .unwrap();
        let (_a2, st2, _) = crate::cuts::cut_loop_recorded(
            &mut d,
            &bounds,
            &p,
            &int_vars,
            crate::cuts::CutOptions::default(),
            &tol,
            || DenseSimplex::new(&p, tol),
        )
        .unwrap();
        assert_eq!(format!("{st1:?}"), format!("{st2:?}"), "割循环终态分叉");
        assert!(a1 >= 1, "revised 侧须施加割（否则未覆盖）");
        let r1 = sp.resolve(&bounds);
        let r2 = d.resolve(&bounds);
        assert_eq!(format!("{:?}", r1.status), "Optimal", "revised 侧割后须 Optimal");
        assert_eq!(format!("{:?}", r2.status), "Optimal", "tableau 侧割后须 Optimal");
        // min 问题：割后根 LP ≥ 未割根 LP − 容差（合法割只会收紧或不变）。
        assert!(
            r1.objective >= root_obj_sp - 1e-6 * (1.0 + root_obj_sp.abs()),
            "revised 割后 {} 低于未割 {}（割非法或求解错值）",
            r1.objective,
            root_obj_sp
        );
        assert!(
            r2.objective >= root_obj_d - 1e-6 * (1.0 + root_obj_d.abs()),
            "tableau 割后 {} 低于未割 {}",
            r2.objective,
            root_obj_d
        );
        // 已知最优 276：两侧割后根 LP 界不得超过它（割切掉最优解即非法）。
        assert!(
            r1.objective <= 276.0 + 1e-6 && r2.objective <= 276.0 + 1e-6,
            "割后根 LP 超过已知最优 276：{}/{}（割切掉最优解）",
            r1.objective,
            r2.objective
        );
    }
}
