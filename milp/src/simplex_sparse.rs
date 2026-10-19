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
//! - **同一反循环/上限纪律**（定价规则除外——本引擎阈上 Devex 与稠密
//!   引擎 Dantzig 的**有意分叉**（阈下回 Dantzig，见下方机制节）；反循环
//!   触发条件与 Bland 路径两引擎完全同款）：连续退化 > 50 或总迭代 >
//!   3×(m+n) 且停滞 ⇒ Bland
//!   （最小下标直返，与定价权重无关），非退化步 + 驻留 8 步释放（每次
//!   调用 ≤64 次重入；总迭代触发逐迭代重申）；迭代上限
//!   max(10_000, 20×(m+n))（`set_iter_limit_hint`
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
//! - **refactorization 纪律**（类比稠密引擎每 200 主元重算 B⁻¹b）：ETA
//!   条数达到 [`ETA_REFLOOR_FLOOR`]（S4 实测定档）即整基重因子化 + 从原始
//!   数据重算 x_B = B⁻¹(b̂ − N·x_N)（漂移双保险；[`REFACTOR_INTERVAL`] 为
//!   历史上限与纪律同源标记）；数值异常（主元过小 / NaN / 检验数非
//!   有限）时若非刚重因子化则先重因子化重算再重试本迭代（确定性的
//!   just_refactored 闸保证至多重试一次，无墙钟决策），仍失败才上报；
//! - x_B 维护与表格引擎同款口径：主元步 x_B[k] −= δt·α_k、x_B[r] =
//!   入基旧界值 + δt；界翻转 x_B[k] −= Δ·α_k；
//! - 定价（**规模分档**：构造期按原始 m×(n+m) 与 [`ENGINE_AUTO_CELL_THRESHOLD`]
//!   比较定死，割行追加不翻转——阈上 = revised 的 auto 默认域 ⇒ **Devex**
//!   （Forrest–Goldfarb 参考框架近似最陡边）；阈下 ⇒ Dantzig，与本引擎
//!   Devex 化之前逐字节同入基序列（A/B 兼容锚）。分档依据（S6 三算例
//!   A/B + 确定性节点预算探针实测）：Devex 在大规模纯 LP 主导类大赢——
//!   lp-sparse-4000x800 端到端 11.27→6.59s、lp-wide-1000x2000 8.84→7.56s
//!   （3 rep min）；但阈下中规模 B&B 类 setcover-300x600（270k 格）上
//!   Devex 顶点轨迹使搜索形态劣化：无 LP 失败/Bland/上限触发的前提下
//!   2000 节点 incumbent 158，未达 Dantzig 1832 节点证优 144（~42s），
//!   端到端 240s 超时（best 156）——分档后阈下恢复 Dantzig 轨迹）。
//!   Devex 机制：入基打分 score_j = rc_j²/w_j（w_j 为参考权重，≥1；
//!   eligible 判定沿用贴下界 −rc / 贴上界 +rc 的阈值与方向逻辑），严格
//!   大于取最大 ⇒ 并列最小下标；每主元后以出基行 αr = e_rᵀB⁻¹A
//!   （BTRAN + 行收集）更新权重：非基 j w_j ← max(w_j, (αr_j/αr_q)²·w_q)、
//!   出基列 w ← max(1, w_q/αr_q²)，越过 [`DEVEX_RESET_THRESHOLD`] 即
//!   重置框架（当前非基集、全 w=1）；solve/resolve 入口与 Phase-1→
//!   Phase-2 切换均重置。rc_j = c_j − y·A_j 每迭代经 BTRAN(c_B) 按 CSC
//!   列现算，无增量漂移（与表格引擎『三入口批量重算』同款新鲜度契约）。
//!   **口径变更（用户确认）**：Devex 域内本引擎迭代序列不再与稠密引擎
//!   逐次一致（Dantzig 期曾逐次一致，7774 次），两引擎**目标恒同**是
//!   红线（多解顶点下 x 允许不同）；v1 不做部分定价（全列扫描在
//!   lp-sparse 实测成本中非主要成分，部分定价另行扰动入基序列，留后续
//!   按证据启用）。
//!
//! # 范围边界（诚实声明）
//! - S3 起提供全套引擎面：`solve` 冷启动 / `resolve` 热启动（snap+传播 →
//!   **三层修复**（S11 重实现，设计源 S10 节留档）：对偶主循环（有界
//!   变量对偶单纯形 + BFRT 长步界翻转 + 对偶 Bland，见
//!   `dual_main_resolve`）→ 对偶修复快路径 `dual_repair`（原样保留）
//!   → 整基重建 + Phase-1） / `add_row`
//!   割行（原尺度不缩放；脏因子化标记 + 下次调用整基重构——F6『宜随后
//!   重构』的最简落法）/ W6 快照（Clone）。**resolve 与稠密引擎的语义
//!   对齐点 = 目标值与状态恒同**（多解顶点下热启动迭代序列与 x 允许
//!   分叉——与 Devex 口径变更同源的既定红线）；三层化的前两层只加速、
//!   任何触发逐层回退，正确性下界 = 老路径。第一层（主循环）的启用受
//!   三闸约束——构造期规模闸 [`DUAL_MAIN_CELL_CEILING`]（大格域无红利
//!   净亏）+ 失败熔断双通道（连续/比率），闸住或熔断时跳过第一层
//!   **仍走 `dual_repair`**（等价老路径行为，不走远贵的整基重建）。
//!   访问器按 cuts.rs 消费面
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
    /// 只读对偶访问面：全部**结构列**的 Phase-2 检验数 rc_j（内部 min
    /// 口径，长度 = num_structural）。根阶段 reduced-cost 变量固定
    /// （branch_bound.rs 方向 3）消费；与定价（price）同源同公式**现算**，
    /// 不是缓存值。
    ///
    /// 契约（两引擎同源语义）：
    /// - **只读**：不改变任何求解状态（不触碰基/表格/因子化/迭代计数）。
    /// - **时序**：仅在上一次 solve/resolve 返回 Optimal 后调用（与
    ///   [`Self::tableau_row_terms`] 同款时序契约——稀疏引擎 add_row 后
    ///   `factorization_dirty` 间隙内 LU 维度与基不符，BTRAN 结果无意义；
    ///   根固定的调用点在根阶段末尾的记账性 resolve 之后，满足契约）。
    /// - **同源公式**：稠密 = 表格现算
    ///   `rc_j = cost_of(j,2) − Σ_i cb_i·(B⁻¹A)_ij`（即定价所用的
    ///   rc_of）；稀疏 = 同公式的 BTRAN 形式（cb → y = B⁻ᵀ·cb →
    ///   `rc_j = cost_of(j,2) − yᵀA_j`，CSC 列点积，与 price 的逐列
    ///   现算同式）。最优基上非基贴下界列 rc ≥ −tol.optimality、贴上界
    ///   列 rc ≤ +tol.optimality（内部 min 口径）。
    /// - **行缩放不变**：rc = c_j − c_B B⁻¹ A_j 对正行缩放不变
    ///   （B̂⁻¹·R = B⁻¹，R 为行缩放阵），引擎内部缩放口径下仍为原始系
    ///   检验数，可与原始空间阈值直接比较。
    /// - 值不跨 solve/resolve 复用（新鲜度规范与 compute_rc_batch 的
    ///   入口重算纪律同源）。
    fn structural_reduced_costs(&self) -> Vec<f64>;
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
    fn structural_reduced_costs(&self) -> Vec<f64> {
        (**self).structural_reduced_costs()
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
    fn structural_reduced_costs(&self) -> Vec<f64> {
        self.structural_reduced_costs()
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
/// 历史上限与纪律同源标记（与稠密引擎 REFRESH_INTERVAL 同源）：eta_over_limit
/// 曾以其为 clamp 上界，S5 复核后固化常数下限 [`ETA_REFLOOR_FLOOR`]，本常量
/// 仅由常量对账测试消费钉值。
#[allow(dead_code)]
const REFACTOR_INTERVAL: u64 = 200;
/// ETA 条数下限（S4 实测定档，见 `eta_over_limit` 文档）。
const ETA_REFLOOR_FLOOR: u64 = 64;
/// LU 阈值主元比：候选 |v| ≥ 该值 × max|v| 才可入选（稳定性下界）。
/// S5 定档证据（代表集防过拟合，目标值全部恒同）：lp-sparse
/// 0.05/0.1/0.25 = 9.74/9.56/10.16s（噪声带内持平）；setcover-300x600
///（revised 路径）= 59.8/41.2/72.8s——0.1 显著最优（0.05 主元质量差、
/// 0.25 填充重），维持 0.1 不改。
const LU_THRESHOLD: f64 = 0.1;
/// Devex 参考权重重置阈值：任一 w_j 越过（或失真为非有限）该值 ⇒ 参考
/// 框架重置为当前非基集（全部 w_j = 1）。取值依据（Forrest–Goldfarb 型
/// 实现通行 10^6 量级，本改动规格建议同量级）：权重只在框架内沿主元链
/// 单调放大（单次主元 (αr_j/αr_q)² 即可乘上 10^6），过大的 w_j 把对应列
/// 打分 rc²/w_j 压到趋零（等效冻结该列入基）且对真实最陡边范数的近似
/// 失真；而过小的阈值又迫使频繁重置、每次丢掉全部累积权重信息退化回
/// Dantzig。1e6 居中兼顾（重置后从全 1 起步需多轮放大才触阈）。
const DEVEX_RESET_THRESHOLD: f64 = 1e6;
/// 对偶主循环（S10 三层 resolve 第一层）单次调用预算 = m + 该基值
///（**每行 1 步 + 64 基值**，近线性于规模）。定档依据（A/B 调参，
/// 2026-10-16）：旧值 2000+50/m 在阈上大格域「慢磨」失控（2.96× 硬
/// 回归）——该域已由 [`DUAL_MAIN_CELL_CEILING`] 规模闸整体禁用；阈下
/// 启用域的两实测点：setcover-100x150 修复最长 ~100+ 步（budget=64 会
/// 掐掉长修复丢红利，节点 19→23）、setcover-300x600 平均 47 步——
/// m+64 同时覆盖两者；且慢磨浪费上限被压到 m+64 步（300 行基 ≈
/// 0.1s/触发），失败熔断（连续 8 次）最坏前期浪费 < 1s。
const DUAL_MAIN_BUDGET_BASE: u32 = 64;
/// 对偶主循环预算每行增量：1（m + 64；定档依据见上）。
const DUAL_MAIN_BUDGET_PER_ROW: u32 = 1;
/// 失败熔断阈值·连续通道（A/B 调参，2026-10-16）：主循环**连续**
/// Fallthrough 达该次数 ⇒ 本实例后续 resolve 跳过主循环、直接走
/// dual_repair（等价老路径行为，正确性下界不变）。
const DUAL_MAIN_FAIL_FUSE: u32 = 8;
/// 失败熔断阈值·比率通道（A/B 调参，2026-10-16）：入口调用数 ≥
/// [`DUAL_MAIN_FUSE_MIN_CALLS`] 且累计 Fallthrough 率 ≥ 40% ⇒ 永久
/// 熔断（不清零——浪费的成本不因后续成功退回）。定档依据（真 dive
/// 诊断）：setcover-300x600 失败率 22%（18 次中 4 次、成功修复平均
/// 47 步 + B&B 节点 591→229 的轨迹红利 ⇒ 净赚，不熔断）；
/// setcover-400x800 失败率 44%（budget=64 下 74 步修复被掐 + 成功
/// 修复 49 步也比 dual_repair 贵 + B&B 节点无减少 ⇒ 每次触发净亏
/// ~40ms，数千次触发 = 2.96× 回归的主体）。40% 线在两者之间；误熔断
/// （健康算例失败率漂移过线）只是回到老路径行为，安全侧。计数跨
/// resolve 生命周期、随快照（Clone）继承/回退，纯计数确定性。
const DUAL_MAIN_FUSE_MIN_CALLS: u64 = 32;
/// 比率通道的失败率阈值（百分数；fall×100 ≥ calls×40）。
const DUAL_MAIN_FAIL_RATE_PCT: u64 = 40;
/// 主循环规模闸（A/B 调参，2026-10-16）：构造期表格规模 m×(n+m) ≥ 该
/// 值 ⇒ 主循环构造期即禁用（resolve 直接走 dual_repair，等价老路径）。
/// 定档依据（实测）：主循环的收益机制是**修复基质量对 B&B 树形态的
/// 改善**（setcover-300x600 = 270k 格：节点 591→229，端到端 0.66-0.78×）
/// ——该红利随规模衰减（割主导树形态、基质量影响变小），且主循环修复
/// 步数（~50）恒为 dual_repair（~20，|主元|最大破键收敛快）的 2.5×，
/// 无红利即净亏（setcover-400x800 = 480k 格：每次触发净亏 ~40ms，节点
/// 反增，端到端 1.4-3×）；两者 engine 层在线指标同构（失败率 22% vs
/// 24%、成功步数 47 vs 49）——红利差异只在 B&B 层，engine 层不可见，
/// 故按 S4 类自适应（ENGINE_AUTO_CELL_THRESHOLD）先例做构造期规模闸，
/// 阈值居 270k（赢）与 480k（亏）之间取 350k。阈下域内失败熔断双通道
/// 仍是安全网。
const DUAL_MAIN_CELL_CEILING: u64 = 350_000;
/// 对偶退化步（对偶比值测试步长 θ ≤ tol.optimality，即对偶目标零改进）
/// 连续达该步数 ⇒ 切入对偶 Bland（S10 钉值 32：eligible 最小下标入基、
/// 禁长步）。量级对齐 primal 主路径 DEGENERATE_STREAK_LIMIT=50 的「数十
/// 步内介入」纪律，取 2 的幂便于人工核对。
const DUAL_BLAND_DEGEN_LIMIT: u32 = 32;
/// 复盘补丁 A（S10 复核中危意见）：对偶 Bland 介入后**仍**连续退化达该
/// 步数 ⇒ 直接 GiveUp 回退老路径，不许只靠迭代 cap（budget）兜底。
/// 取值依据：与 [`DUAL_BLAND_DEGEN_LIMIT`] 同源（2×）——Bland 已是最保守
/// 的确定性入基规则（eligible 最小下标、禁长步），其后同量级证据窗口
///（2×32 步）内仍无一次非退化步，说明当前基上对偶修复在原地打转，
/// 继续迭代只是烧预算；独立阈值让「退化打转」在预算耗尽前就触发回退，
/// 回退链（dual_repair → 整基重建）的正确性下界兜底。
const DUAL_BLAND_GIVEUP_DEGEN: u32 = 64;

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

/// 对偶主循环结果（S10 三层 resolve 第一层）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DualResolveOutcome {
    /// 已恢复原始可行（随后 Phase-2 收尾应为零迭代）。
    Repaired,
    /// 主循环放弃：回退第二层 dual_repair → 第三层整基重建（普通回退
    /// 通道；含「候选空 = 该基上对偶无界 ⇒ 交 Phase-1 判原始不可行」的
    /// 场景——主循环不下结论，只回退）。因子化完好，dual_repair 入口的
    /// dual_feasible 复查可安全执行。
    Fallthrough,
    /// 因子化已损坏（refactor 失败）：factorize_into 奇异失败会清空
    /// lu.perm（f.m 保持旧维度 ⇒ btran 的 y[self.perm[i]] 空索引 panic）、
    /// x_B 重算 NaN 会置 lu = default（m=0 ⇒ 空 y 再被索引 panic）。此态
    /// 上 dual_repair 入口的 dual_feasible 复查**自身会 panic 而非返回
    /// false**——该出口不得进 dual_repair，直达整基重建（rebuild 内部
    /// 对 ±1 对角基 refactor 必成，损坏因子化被完全覆写；老路径对同类
    /// 失败同款语义：数值失败 ⇒ rebuild + Phase-1 诚实映射，绝不 panic）。
    FactorizationBroken,
    /// 墙钟超限：独立 TimedOut 通道——跳过 dual_repair（其循环不读钟，
    /// 超限后跑它等于盲迭代）直达整基重建，run_phase 的 W2 条件读钟
    /// 接管，Phase-1/2 将诚实返回 TimeLimit。
    TimedOut,
}

/// 对偶主循环的退化/Bland 记账（S10 对偶 Bland + 复盘补丁 A）。独立
/// 结构体（而非内联计数）是为了确定性单测能直接驱动计数语义——退化
/// 打转的触发时机是正确性相关行为，值得脱离 LP 实例钉死。
#[derive(Debug, Clone, Default)]
struct DualDegenGuard {
    /// 连续对偶退化步（θ ≤ tol.optimality）计数；任一非退化步清零。
    degen: u32,
    /// 对偶 Bland 是否已介入。一旦介入保持到本次主循环结束（快路径无
    /// 释放机制——与 primal 主路径的 Bland 释放纪律有意分叉：resolve 快
    /// 路径生命周期短，持续退化的正确处置是回退而非继续抗战，由补丁 A
    /// 的独立 GiveUp 兜底）。
    bland: bool,
    /// Bland 介入后的连续退化步计数（补丁 A 的独立触发器）。
    degen_after_bland: u32,
}

impl DualDegenGuard {
    /// 记录一步（theta = 本步对偶步长，tol_opt = tol.optimality）；
    /// 返回 true = 补丁 A 触发（Bland 后仍连续退化达
    /// [`DUAL_BLAND_GIVEUP_DEGEN`] ⇒ 直接 GiveUp 回退老路径）。
    fn on_step(&mut self, theta: f64, tol_opt: f64) -> bool {
        if theta <= tol_opt {
            self.degen = self.degen.saturating_add(1);
            if self.bland {
                self.degen_after_bland = self.degen_after_bland.saturating_add(1);
                if self.degen_after_bland >= DUAL_BLAND_GIVEUP_DEGEN {
                    return true; // 补丁 A：退化打转独立 GiveUp，不等迭代 cap。
                }
            } else if self.degen >= DUAL_BLAND_DEGEN_LIMIT {
                self.bland = true; // 对偶 Bland 介入：eligible 最小下标、禁长步。
            }
        } else {
            self.degen = 0;
            self.degen_after_bland = 0;
        }
        false
    }
}

/// 对偶比值测试 + BFRT 长步选择（S10 主循环核心；纯函数，确定性单测
/// 直接驱动选列与吸收语义）。
///
/// 输入：方向可行候选 `cands`（元素 `(θ_j, j, αr_j, δ_j, range_j)`，**已
/// 按 (θ, 列号) 全序升序**——调用方排序）；`rem` = 出基行到目标界的带
/// 符号剩余违反（方向过滤保证 rem 与每个候选的 −δ_j·αr_j 同号 ⇒ 商
/// 非负）。θ_j = 该列 rc 越零所需的对偶步长（有符号 rc 取对偶可行侧、
/// 除以 |αr_j|）；range_j = 该列界幅 ub−lb（翻转自由度）。
///
/// 输出：`Some((k, t))` = 以候选 k 入基、primal 步长 t（出基变量恰落到
/// 被违反界；此前 0..k 的候选**全部界翻转吸收**——调用方逐列 FTRAN +
/// 翻转传播）；`None` = 候选耗尽（全部翻转后仍到不了目标界 ⇒ 该基上
/// 对偶无界 ⇒ 原始不可行候选，交回退链 Phase-1 判定）或步长异常
///（t 非有限/为负的数值防御）。两种 None 都发生在任何翻转执行之前 ⇒
/// 调用方状态干净。
///
/// 机制（BFRT / Fourer–Maros 型长步比值测试）：沿 θ 升序走候选——入选
/// 列以步长 θ_k 更新对偶变量后，θ_j < θ_k 的候选其 rc 已越零变号，翻到
/// 对面界即「步后 rc=0 落对偶可行侧」（对偶目标继续上升）；候选 j 若在
/// 自己界幅内走不到目标界（t_j > range_j），翻它吸收一部分剩余违反、
/// 取下一候选，直到某候选以 ≤ 其界幅的步长入基清零 rem。无界上界列
/// range = +∞ ⇒ 恒可入基 ⇒ 循环必终于某次入基（或候选耗尽）。
fn dual_ratio_bfrt(cands: &[(f64, usize, f64, f64, f64)], rem: f64) -> Option<(usize, f64)> {
    let mut rem = rem;
    for (k, &(_, _, a, delta, range)) in cands.iter().enumerate() {
        // x_B[r] 每单位 j 移动的变化量 = −δ_j·αr_j（与 rem 同号）。
        let move_per_unit = -delta * a;
        let t = rem / move_per_unit;
        if !t.is_finite() || t < 0.0 {
            return None; // 数值防御（方向过滤后不应出现）
        }
        if t <= range {
            return Some((k, t));
        }
        rem -= move_per_unit * range; // 界翻转吸收后剩余违反递减
    }
    None
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
    /// Devex 参考权重（静态列空间 0..ncols_static；人工列永不入基不设
    /// 权重）。参考框架 = 最近一次重置时的非基集，w_j ≥ 1 单调不降（框架
    /// 内）；定价打分 rc_j²/w_j。重置点与更新规则见 `devex_reset` /
    /// `devex_update`。
    devex_w: Vec<f64>,
    /// Devex 规模分档（构造期按原始 m×(n+m) 定死，割行追加不翻转）：
    /// ≥ [`ENGINE_AUTO_CELL_THRESHOLD`]（revised 的 auto 默认域）才启用
    /// Devex 定价与权重维护；阈下保持 Dantzig——与分档前的 Dantzig 版
    /// 逐字节同入基序列（A/B 兼容锚）。分档依据见模块头定价节。
    devex_active: bool,
    iter_limit_hint: Option<u64>,
    time_deadline: Option<std::time::Instant>,
    pivots_since_time_check: u32,
    /// 对偶主循环（S10 三层 resolve 第一层）累计主元数（跨调用累计，
    /// 口径同 iterations_total）：诊断/测试用——断言三层路径被真实覆盖。
    dual_main_pivots: u64,
    /// 对偶主循环累计 BFRT 界翻转数（跨调用累计）：诊断/测试用。
    dual_main_flips: u64,
    /// 对偶主循环累计入口调用数（含入口红线拒入）：调参诊断用——
    /// 平均每次调用的主元数 = pivots / calls。
    dual_main_calls: u64,
    /// 对偶主循环累计 Fallthrough 次数（含入口红线拒入与各回退触发）：
    /// 调参诊断用——失败占比与熔断阈值的依据。
    dual_main_fallthroughs: u64,
    /// 失败熔断记账：主循环**连续** Fallthrough 次数（Repaired 清零）；
    /// 达 [`DUAL_MAIN_FAIL_FUSE`] ⇒ 后续 resolve 跳过主循环直接走
    /// dual_repair（等价老路径行为）。跨 resolve 生命周期、随快照继承。
    dual_main_consecutive_fails: u32,
    /// 主循环规模闸（构造期定死，同 devex_active 模式）：m×(n+m) ≥
    /// [`DUAL_MAIN_CELL_CEILING`] ⇒ false（大格域主循环无红利净亏，
    /// 见常量 doc）。
    dual_main_enabled: bool,
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
            devex_w: Vec::new(),
            // 规模分档与 effective_engine 同阈同维（m = 约束行数、n = 变量
            // 数、原始问题口径）：revised 的类自适应默认域内才启用 Devex。
            devex_active: (problem.constraints.len() as u64)
                .saturating_mul(n as u64 + problem.constraints.len() as u64)
                >= ENGINE_AUTO_CELL_THRESHOLD,
            iter_limit_hint: None,
            time_deadline: None,
            pivots_since_time_check: 0,
            dual_main_pivots: 0,
            dual_main_flips: 0,
            dual_main_calls: 0,
            dual_main_fallthroughs: 0,
            dual_main_consecutive_fails: 0,
            // 规模闸与 effective_engine 同阈同维（m = 约束行数、n+m 列
            // 数、原始问题口径）：大格域主循环无红利净亏（见
            // DUAL_MAIN_CELL_CEILING doc），构造期即禁用。
            dual_main_enabled: (problem.constraints.len() as u64)
                .saturating_mul(n as u64 + problem.constraints.len() as u64)
                < DUAL_MAIN_CELL_CEILING,
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
        // Devex 参考权重：构造期即对齐静态列空间（全 1 = 全新框架；
        // 与 in_basis/at_upper 同步维护，add_row 插入由此恒安全）。
        s.devex_w = vec![1.0; ncs];
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

    /// 全部结构列的 Phase-2 检验数（只读访问面；契约与公式见 trait 方法
    /// [`LpEngine::structural_reduced_costs`] 文档：cb → y = Bᵀ⁻¹cb（BTRAN）
    /// → rc_j = cost_of(j,2) − yᵀA_j 按列现算，与 price 同源）。
    pub fn structural_reduced_costs(&self) -> Vec<f64> {
        // cb = 基列 Phase-2 成本（人工列 phase 2 成本为 0，退化基中基本
        // 人工列不贡献对偶）。
        let cb: Vec<f64> = self
            .basis
            .iter()
            .map(|&j| self.cost_of(j, 2))
            .collect();
        let y = self.btran_full(&cb);
        let n = self.num_structural;
        let mut out = Vec::with_capacity(n);
        for j in 0..n {
            let mut dot = 0.0;
            for k in self.col_start[j]..self.col_start[j + 1] {
                let yr = y[self.col_rows[k]];
                if yr != 0.0 {
                    dot += yr * self.col_vals[k];
                }
            }
            out.push(self.cost_of(j, 2) - dot);
        }
        out
    }

    /// 第 `row` 行的稀疏表格行：(非零 (列, e_rowᵀB⁻¹A_j) 项 [结构+松弛段、
    /// 列升序]， 该行基变量当前值)。修正单纯形实现：ρ = e_rowᵀB⁻¹ 经
    /// BTRAN，行项 = ρ·A_j 按列现算（与稠密引擎 B⁻¹A 行读取语义一致；
    /// 对账在 ≥1e-6 相对容差口径，不承诺逐位同）。越界行返回 (空, NaN)。
    ///
    /// 时序契约（与稠密引擎的差异点，独立复核指认）：`add_row` 追加割行后
    /// 置 factorization_dirty 推迟重构，该间隙内 lu 仍为旧维度而
    /// basis.len() 已 +1——此时调用本方法，行收集处 col_rows 可能含新行号
    /// 而越界 panic。现有全部消费方（cuts 的割分离、branch_bound 的局部割）
    /// 均在 resolve 之后调用，不触发；新调用方必须先 solve/resolve 再取表行。
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
        // Devex 冷启动重置：参考框架 = 初始非基集（rebuild 后、两阶段前）。
        self.devex_reset();
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
        // Devex 热启动入口重置：参考框架 = 当时非基集（reset-to-全 1 与
        // 具体时刻的非基集选择无关——框架只经权重生效；置于任何贴界/修复
        // 主元之前，其后 dual_repair / Phase-1 的主元照常更新权重）。
        self.devex_reset();
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
        // (3) 基变量越界 / Phase-1 欠账检查 → 三层修复（S10：对偶主循环
        // → dual_repair → 整基重建；正确性下界 = 老路径——前两层只加速、
        // 绝不改变结果，任何触发都逐层回退）。
        let need_repair = self.basics_violated() || self.phase1_pending;
        let mut dual_repaired = false;
        // 主循环启用闸 = 构造期规模闸（大格域无红利净亏，见
        // DUAL_MAIN_CELL_CEILING doc）AND 在线失败熔断（双通道，阈下
        // 域内的安全网）。熔断后跳过第一层**仍走 dual_repair**（等价老
        // 路径行为，不走整基重建——重建 + Phase-1 远贵于 dual_repair，
        // 调参实测该误接一次就把 setcover-300x600 拖垮 4.5×）；未熔断
        // 时的 Repaired 会清零连续计数——健康算例的零星失败不触发连续
        // 通道（比率通道不清零）。
        let main_fused = !self.dual_main_enabled
            || self.dual_main_consecutive_fails >= DUAL_MAIN_FAIL_FUSE
            || (self.dual_main_calls >= DUAL_MAIN_FUSE_MIN_CALLS
                && self.dual_main_fallthroughs.saturating_mul(100)
                    >= self.dual_main_calls.saturating_mul(DUAL_MAIN_FAIL_RATE_PCT));
        if need_repair && !self.phase1_pending && !self.artificials_in_basis() {
            // 熔断跳过时**不**计入 Fallthrough 计数（主循环未被调用，
            // 计数冻结在熔断阈值上；重开尝试的通道是上方 Repaired 清零，
            // 只对未熔断状态有效）。
            let outcome = if main_fused {
                DualResolveOutcome::Fallthrough
            } else {
                self.dual_main_resolve()
            };
            match outcome {
                DualResolveOutcome::Repaired => {
                    dual_repaired = true;
                    self.dual_main_consecutive_fails = 0;
                }
                DualResolveOutcome::TimedOut => {
                    // 墙钟独立 TimedOut 通道：跳过 dual_repair（其循环不
                    // 读钟，超限后跑它等于盲迭代）直达整基重建；run_phase
                    // 的 W2 条件读钟接管，Phase-1/2 诚实返回 TimeLimit。
                    self.rebuild_initial_basis();
                }
                DualResolveOutcome::FactorizationBroken => {
                    // 因子化损坏通道（refactor 失败）：损坏态（perm 空 /
                    // m=0）上 dual_repair 入口的 dual_feasible→btran 会索引
                    // 越界 panic 而非返回 false——跳过 dual_repair 直达整
                    // 基重建（rebuild 内部对 ±1 对角基 refactor 必成，因子
                    // 化与中间态一并覆写）。老路径对同类失败（dual_repair
                    // 内 refactor Err ⇒ GiveUp ⇒ rebuild）同款语义：数值
                    // 失败 ⇒ rebuild + Phase-1 诚实映射，绝不 panic。
                    self.rebuild_initial_basis();
                }
                DualResolveOutcome::Fallthrough => {
                    // 第二层：dual_repair（原样保留，经长期验证的稳健回
                    // 退；其入口 dual_feasible 复查同时是第一层 BFRT 中途
                    // GiveUp 后状态污染的防线，见 dual_main_loop 注释）。
                    if !main_fused {
                        self.dual_main_fallthroughs += 1;
                        self.dual_main_consecutive_fails =
                            self.dual_main_consecutive_fails.saturating_add(1);
                    }
                    if self.dual_repair() == DualRepair::Repaired {
                        dual_repaired = true;
                    } else {
                        self.rebuild_initial_basis();
                    }
                }
            }
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

    // ---- 内部：对偶主循环（S10 三层 resolve 第一层） ----

    /// 对偶主循环（S10 设计重实现，仅 revised 引擎热启动路径）：从（对偶
    /// 可行、原始违约）基出发的有界变量对偶单纯形。规则：
    /// - 离基行 = |违反| 最大，并列最小行号（绝对违反量，非 dual_repair
    ///   的界幅归一——主循环求速度，先修最痛的行）；
    /// - 对偶比值测试最小 θ 并列最小列号（(θ, 列号) 全序，确定性破键）；
    /// - BFRT 长步：升序逐列界翻转吸收剩余违反（`dual_ratio_bfrt`）；
    /// - 对偶退化步（θ ≤ tol.optimality）连续 [`DUAL_BLAND_DEGEN_LIMIT`]
    ///   步 ⇒ 对偶 Bland（eligible 最小下标、禁长步）；
    /// - 复盘补丁 A：Bland 后仍连续退化达 [`DUAL_BLAND_GIVEUP_DEGEN`] ⇒
    ///   直接 GiveUp 回退老路径（不等迭代 cap）。
    ///
    /// 回退触发表（逐条见函数内注释）：对偶不可行入口扫描不过 / x_B 或
    /// rc 非有限 / 候选空或耗尽（对偶无界 ⇒ 交 Phase-1 判原始不可行）/
    /// 主元过小（refactor 单次确定性重试后仍小）/ 步长异常 / 迭代超限
    ///（m+64，快路径量级）/ 补丁 A 的 Bland 后仍退化——以上 ⇒
    /// [`DualResolveOutcome::Fallthrough`]（因子化完好，dual_repair 可
    /// 安全接手）；do_pivot 或主元闸的 refactor 失败（因子化损坏，
    /// dual_repair 入口复查会 panic）⇒ [`DualResolveOutcome::
    /// FactorizationBroken`] 直达整基重建；墙钟超限单独走
    /// [`DualResolveOutcome::TimedOut`] 通道。连续 [`DUAL_MAIN_FAIL_FUSE`]
    /// 次 Fallthrough ⇒ resolve 接线熔断（后续跳过主循环，见 resolve）。
    /// **正确性下界 = 老路径**。
    ///
    /// 与第二层 dual_repair 的算法分叉是有意的（同为对偶单纯形、规则
    /// 不同）：主循环求速度（BFRT 长步、θ 最小列），dual_repair 是经
    /// 长期验证的稳健回退（归一化行选择、|主元| 最大破键）。两层的
    /// rc/BTRAN 现算与主元执行纪律（|主元| 闸、refactor 单次重试）同款。
    fn dual_main_resolve(&mut self) -> DualResolveOutcome {
        // 预算 = m + 64（定档依据见常量 doc）。
        let m = self.basis.len();
        let budget =
            DUAL_MAIN_BUDGET_BASE.saturating_add((m as u32).saturating_mul(DUAL_MAIN_BUDGET_PER_ROW));
        self.dual_main_loop(budget)
    }

    /// [`Self::dual_main_resolve`] 的预算参数化内核：预算是唯一的可注入
    /// 参数，确定性单测以极小预算直接驱动「预算耗尽 ⇒ Fallthrough ⇒
    /// 回退链仍给出正确解」的诚实回退路径。
    fn dual_main_loop(&mut self, budget: u32) -> DualResolveOutcome {
        self.dual_main_calls += 1;
        // 入口对偶可行性 rc 扫描红线：对偶不可行基不跑主循环（对偶单纯形
        // 的前提），回退第二层 dual_repair——其入口同款 dual_feasible 复查
        // 同时是本循环 BFRT 中途 GiveUp 后状态污染的防线（见下方两处
        // GiveUp 的【复盘补丁 B】注释）。
        if !self.dual_feasible() {
            self.dual_main_fallthroughs += 1;
            return DualResolveOutcome::Fallthrough;
        }
        let m = self.basis.len();
        let ncs = self.ncols_static();
        let mut guard = DualDegenGuard::default();
        let mut iters: u32 = 0;
        // 候选缓冲 (θ_j, j, αr_j, δ_j, range_j)，跨迭代复用容量。
        let mut cands: Vec<(f64, usize, f64, f64, f64)> = Vec::new();
        loop {
            if iters >= budget {
                return DualResolveOutcome::Fallthrough; // 迭代超限
            }
            // 墙钟：W2 条件读钟纪律（deadline 未注入 ⇒ 零墙钟读点）。
            if let Some(deadline) = self.time_deadline {
                self.pivots_since_time_check = self.pivots_since_time_check.saturating_add(1);
                if self.pivots_since_time_check >= TIME_CHECK_INTERVAL {
                    self.pivots_since_time_check = 0;
                    if std::time::Instant::now() >= deadline {
                        return DualResolveOutcome::TimedOut; // 独立通道，防超限后盲迭代
                    }
                }
            }
            // y = BTRAN(c_B)：rc 每迭代现算（与 run_phase 同款新鲜度契约，
            // 无增量漂移）。
            for i in 0..m {
                self.cb[i] = self.cost_of(self.basis[i], 2);
            }
            let y = self.btran_full(&self.cb);
            if y.iter().any(|v| !v.is_finite()) {
                return DualResolveOutcome::Fallthrough; // 对偶变量非有限
            }
            // (1) 离基行 = |违反| 最大，并列最小行号（升序扫描 + 严格大于
            // 才替换 ⇒ 并列保留先见的最小行号）。
            let mut worst: Option<(usize, f64, bool)> = None; // (行, 违反, below)
            for i in 0..m {
                let (lb, ub) = self.col_bounds(self.basis[i]);
                let v_lb = lb - self.x_b[i];
                let v_ub = self.x_b[i] - ub;
                let (viol, below) = if v_lb > v_ub { (v_lb, true) } else { (v_ub, false) };
                if !viol.is_finite() {
                    return DualResolveOutcome::Fallthrough; // x_B 非有限
                }
                if viol <= self.tol.feasibility {
                    continue;
                }
                match worst {
                    // 旧值 ≥ 当前（含并列）⇒ 保留先见者 ⇒ 并列取最小行号
                    //（与 dual_repair 的 w >= norm 同款破键方向）。
                    Some((_, w, _)) if w >= viol => {}
                    _ => worst = Some((i, viol, below)),
                }
            }
            let Some((r, _, below)) = worst else {
                return DualResolveOutcome::Repaired; // 全部基变量回界 ⇒ 修复完成
            };
            let (lb_r, ub_r) = self.col_bounds(self.basis[r]);
            let target = if below { lb_r } else { ub_r };
            // (2) BTRAN(e_r) 取表格行 ρ，扫全部静态非基列现算 αr_j 与 rc_j，
            // 收集方向可行候选（below ⇒ δ·αr < 0 / above ⇒ δ·αr > 0——
            // 只有向界移动能缩小本行违反的列才 eligible）。
            let mut e = vec![0.0f64; m];
            e[r] = 1.0;
            let rho = self.btran_full(&e);
            if rho.iter().any(|v| !v.is_finite()) {
                return DualResolveOutcome::Fallthrough; // ρ 非有限
            }
            cands.clear();
            for j in 0..ncs {
                if self.in_basis[j] {
                    continue;
                }
                let (lb, ub) = self.col_bounds(j);
                if ub - lb <= self.tol.feasibility {
                    continue; // 固定列：无对偶语义、无翻转自由度
                }
                let mut a = 0.0;
                for k in self.col_start[j]..self.col_start[j + 1] {
                    let rr = rho[self.col_rows[k]];
                    if rr != 0.0 {
                        a += rr * self.col_vals[k];
                    }
                }
                if !a.is_finite() || a.abs() < self.tol.pivot {
                    continue; // 主元过小的列不参与（含非有限数值防御）
                }
                let delta = if self.at_upper[j] { -1.0 } else { 1.0 };
                let da = delta * a;
                if below && da >= 0.0 {
                    continue; // 方向不可行
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
                    return DualResolveOutcome::Fallthrough; // rc 非有限
                }
                let signed = if self.at_upper[j] { -rc } else { rc };
                let theta = (signed / a.abs()).max(0.0);
                if !theta.is_finite() {
                    return DualResolveOutcome::Fallthrough; // θ 非有限
                }
                cands.push((theta, j, a, delta, ub - lb));
            }
            if cands.is_empty() {
                // 候选空：该基上对偶无界 ⇒ 原始不可行候选——主循环不在此
                // 下结论，回退链交 Phase-1 判定（老路径同款语义；此时未
                // 执行任何翻转/主元，状态干净）。
                return DualResolveOutcome::Fallthrough;
            }
            // (3) 选列与长步。Bland 期：eligible 最小下标、禁长步（无翻转
            // 吸收）；否则 (θ, 列号) 全序升序 + BFRT 逐列吸收。
            let (pick, t, theta_final, flip_end): (usize, f64, f64, usize);
            if guard.bland {
                let mut best = 0usize;
                for idx in 1..cands.len() {
                    if cands[idx].1 < cands[best].1 {
                        best = idx; // 并列最小下标（j 互异 ⇒ 唯一）
                    }
                }
                let (_, _, a, delta, _) = cands[best];
                let t_bland = (target - self.x_b[r]) / (-delta * a);
                if !t_bland.is_finite() || t_bland < 0.0 {
                    return DualResolveOutcome::Fallthrough; // 步长异常（Bland 期未翻转，状态干净）
                }
                pick = best;
                t = t_bland;
                theta_final = cands[best].0;
                flip_end = 0; // 禁长步
            } else {
                cands.sort_unstable_by(|x, w| {
                    // θ 已过滤为有限 ⇒ partial_cmp 恒 Some；(θ, j) 二元组
                    // 因 j 互异而唯一 ⇒ 排序结果与算法无关地确定。
                    x.0.partial_cmp(&w.0)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then(x.1.cmp(&w.1))
                });
                match dual_ratio_bfrt(&cands, target - self.x_b[r]) {
                    Some((k, tk)) => {
                        pick = k;
                        t = tk;
                        theta_final = cands[k].0;
                        flip_end = k; // 0..k 的候选全部界翻转吸收
                    }
                    None => {
                        // 候选耗尽/步长异常——纯函数阶段未执行任何翻转，
                        // 状态干净，直接回退（对偶无界 ⇒ 交 Phase-1）。
                        return DualResolveOutcome::Fallthrough;
                    }
                }
            }
            let q = cands[pick].1;
            let delta_q = cands[pick].3;
            // (4) 执行翻转（先于主元）：被翻列走到对面界、x_B 全行传播；
            // 入选列以 θ_final 更新对偶变量后这些列 rc 越零变号、「步后
            // rc=0 落对偶可行侧」（本引擎 y/rc 每迭代 BTRAN 现算，无需
            // 显式对偶更新）。
            for idx in 0..flip_end {
                let (_, j, _, delta, range) = cands[idx];
                self.ftran_col(j);
                self.do_flip(j, delta, range);
                self.dual_main_flips += 1;
            }
            // (5) 入基主元：FTRAN 列、|主元| 闸（refactor 单次确定性重试，
            // just_refactored 闸保证至多一次，与 dual_repair/run_phase 同款）。
            self.ftran_col(q);
            if self.alpha[r].abs() < self.tol.pivot {
                // 【复盘补丁 B：状态污染声明之一】走到这三个 GiveUp 出口
                // 时，(4) 的 BFRT 翻转（若有）**已执行**——留下「已翻转未
                // 主元」中间态：被翻列 x_j 已到新界侧，但其 rc 还在旧界侧
                // 符号（y 未更新 ⇒ 该列对偶不可行），x_b 已含翻转传播。
                // 出口分两类，安全依赖不同（勿混）：
                // · refactor Err 出口：因子化已损坏（perm 清空 / m=0）——
                //   dual_repair 入口的 dual_feasible 复查在损坏因子化上
                //   **必然 panic 而非返回 false**，补丁 B 的复查防线在此
                //   不成立 ⇒ 走 FactorizationBroken 直达整基重建（±1 对角
                //   基 refactor 必成，污染与损坏因子化一并覆写），不依赖
                //   任何前置复查。
                // · 两个「主元仍过小」出口：因子化完好（refactor 成功或
                //   未动），仅翻转污染 ⇒ 普通回退，安全性**完全依赖回退
                //   链**：第二层 dual_repair 入口的 dual_feasible() 复查会
                //   发现被翻列的对偶违反 ⇒ GiveUp ⇒ 第三层整基重建从零
                //   重来。⚠️ 禁止在 dual_repair 之前插入任何不经
                //   dual_feasible 复查的修复/续跑路径——那会在对偶不可行
                //   基上跑对偶算法，静默产生错误结果。
                if !self.just_refactored {
                    if self.refactor().is_err() {
                        return DualResolveOutcome::FactorizationBroken; // refactor 出错（因子化损坏态）
                    }
                    self.ftran_col(q);
                    if self.alpha[r].abs() < self.tol.pivot {
                        return DualResolveOutcome::Fallthrough; // 重试后主元仍过小（因子化完好）
                    }
                } else {
                    return DualResolveOutcome::Fallthrough; // 刚重因子化仍过小（因子化完好）
                }
            }
            if self.do_pivot(q, r, delta_q, t, !below).is_err() {
                // 【复盘补丁 B：状态污染声明之二】do_pivot 的 Err 唯一来源
                // 是其内部 ETA 越限触发的 refactor 失败——因子化已损坏
                // （perm 清空 / m=0），且发生在 (4) 翻转已执行之后。与上方
                // refactor Err 出口同类：dual_repair 入口复查在损坏态上必
                // panic ⇒ 走 FactorizationBroken 直达整基重建覆写一切，
                // 不依赖 dual_feasible 复查；同样不得插入不经复查的修复
                // 路径。
                return DualResolveOutcome::FactorizationBroken; // do_pivot 内 refactor 失败（因子化损坏态）
            }
            iters += 1;
            self.iter_call += 1;
            self.iterations_total += 1;
            self.dual_main_pivots += 1;
            if guard.on_step(theta_final, self.tol.optimality) {
                // 复盘补丁 A：Bland 后仍连续退化达独立阈值 ⇒ 直接 GiveUp
                // 回退老路径，不等迭代 cap 兜底。此刻上一步主元已完整
                // 落地（翻转与主元成对执行），基自洽、无中间态污染——
                // dual_repair 可安全接手。
                return DualResolveOutcome::Fallthrough;
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

    /// ETA 越限判定（S4 实测定档）：ETA 条数达到 [`ETA_REFLOOR_FLOOR`] 即
    /// 整基重因子化。曾设想按 LU fill 做 flops 平衡（ETA 施加成本 ∝
    /// Σ|eta.nnz|/次三角求解 vs 重构成本 O(LU fill)），但左视 LU 的
    /// l_rows/u_cols/u_diag 外层长度恒为 m、并非真实非零计数，稠密场景恒落
    /// clamp 下限——独立复核指认后按实测语义固化为常数（setcover-300x600
    /// 全程 A/B：floor 8/16 重构过频反亏 124.9s、64 最优 47.1s、128/200
    /// ETA 施加占优 99.7-124.9s）。[`REFACTOR_INTERVAL`] 保留为文档化的
    /// 历史上限与漂移控制纪律同源标记（常量对账测试钉值）。确定性：只依赖
    /// ETA 计数。
    fn eta_over_limit(&self) -> bool {
        self.etas.len() as u64 >= ETA_REFLOOR_FLOOR
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
            // Devex 权重随静态列空间同步插入（新列从全 1 起步；与
            // in_basis/at_upper 同位插入保持列号对齐）。
            self.devex_w.insert(c, 1.0);
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

    /// Devex/Dantzig/Bland 定价（**与 simplex.rs::price 的分叉点**：本
    /// 引擎构造期规模分档（m×(n+m) ≥ [`ENGINE_AUTO_CELL_THRESHOLD`]）阈上
    /// Devex 打分、稠密引擎恒 Dantzig；阈下 Dantzig 与本引擎旧版逐字节
    /// 同轨迹。用户确认口径——两引擎目标恒同是红线，迭代序列不承诺逐次
    /// 一致；eligible 判定与 Bland 路径完全同款）。rc 每迭代经 BTRAN
    /// 现算，天然新鲜；Devex 域内合格列打分 score_j = rc_j²/w_j（w_j 为
    /// 参考权重，见 `devex_update`），严格大于取最大 ⇒ 并列最小下标破键
    /// （确定性不变）。人工列与固定列（ub−lb ≤ feas）永不入基。
    fn price(&mut self, phase: u8) -> Pricing {
        let ncs = self.ncols_static();
        if self.devex_w.len() != ncs {
            // 防御：权重与静态列空间不同步（各入口重置后不应发生）——按
            // 当前非基集口径重置。
            self.devex_reset();
        }
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
                return Pricing::Candidate(j, delta); // 最小下标（与权重无关）
            }
            let score = if self.devex_active {
                // Devex 打分：rc² 方向无关（贴下界 rc<0 / 贴上界 rc>0 平方后
                // 同向可比）；w_j ≥ 1 恒有限无除零。
                rc * rc / self.devex_w[j]
            } else {
                // 阈下 Dantzig：与旧版（Dantzig 定价期）逐字节同公式——
                // 入基序列、迭代数、顶点全部不变（A/B 兼容锚）。
                if delta > 0.0 { -rc } else { rc }
            };
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

    /// Devex 参考权重重置：参考框架 = 当前非基集，全部 w_j = 1（静态列
    /// 空间；人工列永不入基不设权重）。调用点：solve/resolve 入口、
    /// Phase-1→Phase-2 切换、任一权重越过 [`DEVEX_RESET_THRESHOLD`]。
    fn devex_reset(&mut self) {
        let ncs = self.ncols_static();
        self.devex_w.clear();
        self.devex_w.resize(ncs, 1.0);
    }

    /// Devex 参考权重更新（Forrest–Goldfarb 标准形式）。**必须在基/ETA
    /// 变更之前调用**：更新所用的出基表格行 αr = e_rᵀB⁻¹A 属于主元前的
    /// 基（此时 `btran_full` 读到的 LU+ETA 恰为主元前因子化）。行收集
    /// 复用本引擎既有机制（BTRAN(e_r) + CSC 逐列点积，与 dual_repair /
    /// tableau_row_terms 同款，不另起稠密实现）；主元元 αr_q 直接取
    /// `self.alpha[r]`（= e_rᵀB⁻¹A_q，与收集值数学恒同，省一列收集）。
    ///
    /// 更新规则：非基 j（αr_j ≠ 0）w_j ← max(w_j, (αr_j/αr_q)²·w_q)；
    /// 出基变量转入非基集 w ← max(1, w_q/αr_q²)（人工列出基无权重条目，
    /// 跳过）。任一权重越过（或失真为非有限）[`DEVEX_RESET_THRESHOLD`]
    /// ⇒ 参考框架重置。Bland 模式打分不用权重，但更新照常维护——释放
    /// Bland 后权重仍与基演进一致。
    fn devex_update(&mut self, q: usize, r: usize, leaving: usize) {
        if !self.devex_active {
            return; // 阈下 Dantzig：不维护权重（打分也不消费）
        }
        let ncs = self.ncols_static();
        if self.devex_w.len() != ncs {
            self.devex_reset(); // 防御（各入口重置后不应发生）
        }
        let a_rq = self.alpha[r];
        if a_rq == 0.0 {
            // 主元为零不可除：上游三调用点（run_phase/dual_repair/
            // cleanup_after_phase1）均已按 |αr_q| ≥ tol.pivot 拦截，此处
            // 仅为确定性防御——跳过更新，打分退回近似（不谎报、不上抛）。
            return;
        }
        let w_q = self.devex_w[q]; // q 恒为静态列（定价只扫静态段）
        let m = self.basis.len();
        // ρ = e_rᵀB⁻¹（主元前基）。
        let mut e = vec![0.0f64; m];
        e[r] = 1.0;
        let rho = self.btran_full(&e);
        let mut breach = false;
        for j in 0..ncs {
            if j == q || self.in_basis[j] {
                // 入基列自身更新为恒等式（max(w_q, w_q)）；基列在行空间
                // 为单位向量（αr_j ∈ {0,1}）且不入定价，不设权重。
                continue;
            }
            let mut a = 0.0;
            for k in self.col_start[j]..self.col_start[j + 1] {
                let rr = rho[self.col_rows[k]];
                if rr != 0.0 {
                    a += rr * self.col_vals[k];
                }
            }
            if a == 0.0 {
                continue;
            }
            let ratio = a / a_rq;
            let cand = ratio * ratio * w_q;
            // cand 为 NaN（上游数值异常）时比较为假、原值保留——权重
            // 向量永不引入 NaN；inf/越阈则触发重置兜底。
            if cand > self.devex_w[j] {
                self.devex_w[j] = cand;
                if !cand.is_finite() || cand > DEVEX_RESET_THRESHOLD {
                    breach = true;
                }
            }
        }
        if leaving < ncs {
            // 出基变量转入非基集：按规范**赋值** max(1, w_q/αr_q²)（其在
            // 基期间名义范数已变，旧权重作废）。
            let cand = (w_q / (a_rq * a_rq)).max(1.0);
            self.devex_w[leaving] = cand;
            if !cand.is_finite() || cand > DEVEX_RESET_THRESHOLD {
                breach = true;
            }
        }
        if breach {
            self.devex_reset();
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
    ///（追加 (r, α)），ETA 条数达到 [`ETA_REFLOOR_FLOOR`] 即整基重因子化 +
    /// 重算（失败以 Err 上抛，由 run_phase 的数值异常通道处理）。
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
        // Devex 权重更新：须在 ETA 追加与基交换之前（αr 为主元前基的
        // 表格行，见 `devex_update` 文档）。
        self.devex_update(q, r, leaving);
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
                    // Devex：Phase-1 基演进（含清理期 t=0 主元）期间的权重
                    // 对 Phase-2 已失参考意义——切换点重置为当前非基集。
                    self.devex_reset();
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
    fn structural_reduced_costs(&self) -> Vec<f64> {
        self.structural_reduced_costs()
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

    /// S6：Devex 规模分档钉值——构造期 m×(n+m)（原始问题口径）≥
    /// [`ENGINE_AUTO_CELL_THRESHOLD`] 才启用 Devex；阈下 Dantzig 与旧版
    /// 逐字节同轨迹（A/B 兼容锚，分档依据见模块头定价节）。阈上算例取
    /// m=998、n=5（998×1003 ≈ 1.001M 格刚好过阈）——变量极少 ⇒ 两引擎
    /// 求解都只需个位数主元，默认套件可承受；同时以目标恒同对账覆盖
    /// Devex 打分 + 权重更新路径的端到端正确性。
    #[test]
    fn devex_gate_by_instance_size() {
        // 小算例：阈下 ⇒ Dantzig。
        let p = load_instance(&milp_instance("milp-mixed-5")).unwrap();
        let sp = SparseSimplex::new(&p, Tolerances::default()).unwrap();
        assert!(!sp.devex_active, "阈下须保持 Dantzig");
        // 998 行 × 5 变量（行系数轮转单变量 ⇒ 最紧行 x_j ≤ 1）：阈上 ⇒
        // Devex，且最优目标 = 5.0 与稠密引擎恒同。
        let n = 5usize;
        let m = 998usize;
        let vars: Vec<Variable> = (0..n).map(|j| var(&format!("x{j}"), 0.0, 10.0)).collect();
        let obj: Vec<(usize, f64)> = (0..n).map(|j| (j, 1.0)).collect();
        let cons: Vec<(Vec<(usize, f64)>, RelOp, f64)> = (0..m)
            .map(|i| {
                let rhs = 1.0 + 0.001 * (i / n) as f64; // 同变量的更松行冗余
                (vec![(i % n, 1.0)], RelOp::Le, rhs)
            })
            .collect();
        let p2 = prob(Sense::Maximize, vars, obj, cons);
        let sp2 = SparseSimplex::new(&p2, Tolerances::default()).unwrap();
        assert!(sp2.devex_active, "阈上须启用 Devex");
        let (st, obj2) = check_vs_dense(&p2, "devex-gate-large");
        assert_eq!(format!("{st:?}"), "Optimal");
        assert!((obj2 - 5.0).abs() <= 1e-6, "Devex 域目标 {obj2} ≠ 5.0");
    }

    #[test]
    fn iter_limit_hint_truncates_honestly() {
        // F4 接线点证据：hint 钳到下限 1000（与 simplex.rs 同款
        // h.max(1000)）；生成 setcover 型 LP（300x200，固定种子）全程
        // 1427 迭代 > 1000 ⇒ hint=1 ⇒ IterationLimit + objective NaN
        //（诚实截断，不冒充 Optimal）。锚口径注记（S6 规模分档）：本例
        // m×(n+m)=200×500=100k 阈下 ⇒ Dantzig 轨迹（1427 为 Dantzig 期
        // 实测并经分档复原逐字节重现）；Devex 全开的中途版本曾降至 951
        //（< 1000 截断不触发）并临时改锚 300x400，分档后复原原锚。
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

    /// S5：与稠密引擎（simplex.rs）的退化/重构纪律常量逐项对齐钉死。
    /// dense 侧常量为私有（不可程序读取），此处按其源码钉值
    ///（simplex.rs:97-136）：任何一侧改动使口径分叉时，本测试提醒人工
    /// 复核对齐（允许有意分叉——但须改钉值并记因）。
    #[test]
    fn bland_and_refactor_constants_parity_with_tableau() {
        // simplex.rs:97 / 本文件同值。
        assert_eq!(DEGENERATE_STREAK_LIMIT, 50);
        // simplex.rs:108/111/115。
        assert_eq!(BLAND_DWELL, 8);
        assert_eq!(BLAND_REARM_MAX, 64);
        assert_eq!(BLAND_TOTAL_FACTOR, 3);
        // simplex.rs:125（RATIO_BAND_CAP 实测踩坑封顶）。
        assert_eq!(RATIO_BAND_CAP, 1e3);
        // simplex.rs:136（W2 条件读钟间隔）。
        assert_eq!(TIME_CHECK_INTERVAL, 64);
        // simplex.rs:113 REFRESH_INTERVAL=200：本引擎 REFACTOR_INTERVAL 为
        // 历史上限与纪律同源标记（eta_over_limit 曾以其为 clamp 上界，复核
        // 后固化常数下限 ETA_REFLOOR_FLOOR），语义对应『每 200 主元用原始
        // 数据重算』的漂移控制纪律。
        assert_eq!(REFACTOR_INTERVAL, 200);
        assert_eq!(ETA_REFLOOR_FLOOR, 64);
        // S5 定档证据见 LU_THRESHOLD 文档注释；防过拟合代表集实测维持 0.1。
        assert_eq!(LU_THRESHOLD, 0.1);
    }

    /// S5：并列破键敏感的退化 LP（全对称系数 ⇒ 任一定价规则下打分全并列
    ///（Dantzig 与 Devex 同为全等 rc/全等权重）+ 退化顶点）：resolve 序列
    ///（逐步固定）与稠密引擎目标恒同——Bland 触发/释放路径的值级行为
    /// 对账。重锚注（S10 三层化）：口径 = 目标值与状态恒同（热启动迭代
    /// 序列允许分叉）。
    #[test]
    fn degenerate_tie_heavy_parity_with_dense() {
        // 8 变量全对称 max Σx s.t. 两两容量行（全等系数 ⇒ 全部比值并列）。
        let n = 8usize;
        let vars: Vec<Variable> = (0..n).map(|j| var(&format!("x{j}"), 0.0, 1.0)).collect();
        let obj: Vec<(usize, f64)> = (0..n).map(|j| (j, 1.0)).collect();
        let mut cons: Vec<(Vec<(usize, f64)>, RelOp, f64)> = Vec::new();
        for i in 0..4 {
            // 每行 4 个全等系数、rhs=2：任意 2 个变量取 1 即触界（退化）。
            let terms: Vec<(usize, f64)> = (0..4).map(|k| ((i * 2 + k) % n, 1.0)).collect();
            cons.push((terms, RelOp::Le, 2.0));
        }
        let p = prob(Sense::Maximize, vars, obj, cons);
        let tol = Tolerances::default();
        let root = Bounds::from_problem(&p);
        let mut sp = SparseSimplex::new(&p, tol).unwrap();
        let mut d = DenseSimplex::new(&p, tol).unwrap();
        let s0 = sp.solve(&root);
        let d0 = d.solve(&root);
        assert_eq!(s0.status, d0.status);
        assert!((s0.objective - d0.objective).abs() <= 1e-9);
        // 逐步固定序列（每次固定产生并列破键敏感的退化子问题）。
        let mut b = root.clone();
        for j in 0..n {
            b.lb[j] = 1.0;
            let rs = sp.resolve(&b);
            let rd = d.resolve(&b);
            assert_eq!(
                format!("{:?}", rs.status),
                format!("{:?}", rd.status),
                "step {j}"
            );
            if rs.status == LpStatus::Optimal {
                assert!(
                    (rs.objective - rd.objective).abs() <= 1e-9 * (1.0 + rd.objective.abs()),
                    "step {j}: {} vs {}",
                    rs.objective,
                    rd.objective
                );
            }
        }
    }

    // ---- G. S10 重实现：对偶主循环 resolve（三层化）----

    /// G1（S10 测试 1 对齐）：对偶主循环路径下的热启动 resolve——目标与
    /// 冷解、与稠密引擎 resolve 三方恒同，且 `dual_main_pivots > 0` 证明
    /// 第一层路径被真实覆盖（不是静默回退到 dual_repair 的假绿）。
    #[test]
    fn dual_main_resolve_matches_cold_and_dense_objective() {
        let p = load_instance(&milp_instance("knapsack-10")).unwrap();
        let tol = Tolerances::default();
        let root = Bounds::from_problem(&p);
        // 形态 1：x0 收紧上界（分支形态）。
        let mut b = root.clone();
        b.ub[0] = 0.0;
        let mut sp = SparseSimplex::new(&p, tol).unwrap();
        let _ = sp.solve(&root);
        let warm = sp.resolve(&b);
        let cold = {
            let mut sp2 = SparseSimplex::new(&p, tol).unwrap();
            sp2.solve(&b)
        };
        let mut d = DenseSimplex::new(&p, tol).unwrap();
        let _ = d.solve(&root);
        let dense = d.resolve(&b);
        assert_eq!(format!("{:?}", warm.status), format!("{:?}", cold.status));
        assert_eq!(format!("{:?}", warm.status), format!("{:?}", dense.status));
        if warm.status == LpStatus::Optimal {
            assert!((warm.objective - cold.objective).abs() <= 1e-7 * (1.0 + cold.objective.abs()));
            assert!((warm.objective - dense.objective).abs() <= 1e-7 * (1.0 + dense.objective.abs()));
        }
        // 形态 2：dive 序列连续固定（每次固定都是一次热启动修复机会）。
        let mut sp2 = SparseSimplex::new(&p, tol).unwrap();
        let mut d2 = DenseSimplex::new(&p, tol).unwrap();
        let _ = d2.solve(&root);
        let base = sp2.solve(&root);
        let root_x = base.x.clone();
        let mut b2 = root.clone();
        for j in 0..5 {
            let v = root_x[j].round().max(b2.lb[j]).min(b2.ub[j]);
            b2.lb[j] = v;
            b2.ub[j] = v;
            let rs = sp2.resolve(&b2);
            let rd = d2.resolve(&b2);
            assert_eq!(
                format!("{:?}", rs.status),
                format!("{:?}", rd.status),
                "step {j}: 状态分叉"
            );
            if rs.status == LpStatus::Optimal {
                assert!(
                    (rs.objective - rd.objective).abs() <= 1e-7 * (1.0 + rd.objective.abs()),
                    "step {j}: {} vs {}",
                    rs.objective,
                    rd.objective
                );
            }
        }
        // 路径覆盖：两个实例的 resolve 序列合计必须真的跑过对偶主循环主元。
        assert!(
            sp.dual_main_pivots + sp2.dual_main_pivots > 0,
            "对偶主循环未被覆盖（dual_main_pivots 全 0）——测试退化为旧路径对账"
        );
    }

    /// G2（S10 测试 2 对齐）：对偶比值测试选列——最小 θ 胜出；θ 并列时
    /// (θ, 列号) 全序取小列号。输入按调用方约定（升序）手工构造。
    #[test]
    fn dual_ratio_bfrt_selects_min_theta_then_min_col() {
        // rem = 4：三个候选都能以小步长直接入基（range 充裕）。
        // (θ, j, αr, δ, range)：move_per_unit = −δ·αr = +1（与 rem 同号）。
        let cands = [
            (0.1, 9usize, 1.0, -1.0, 10.0),
            (0.2, 2usize, 1.0, -1.0, 10.0),
            (0.3, 7usize, 1.0, -1.0, 10.0),
        ];
        // 最小 θ（候选 0，尽管列号 9 最大）胜出；t = rem/1 = 4。
        assert_eq!(dual_ratio_bfrt(&cands, 4.0), Some((0, 4.0)));
        // θ 并列（升序输入下小列号在前）：仍取 idx 0 = 列 3。
        let tie = [
            (0.5, 3usize, 1.0, -1.0, 10.0),
            (0.5, 7usize, 1.0, -1.0, 10.0),
        ];
        assert_eq!(dual_ratio_bfrt(&tie, 2.0), Some((0, 2.0)));
        // 步长恰等于界幅（t == range）⇒ 入基而非翻转（闭区间语义）。
        assert_eq!(dual_ratio_bfrt(&tie, 10.0), Some((0, 10.0)));
    }

    /// G3（S10 测试 3 对齐）：BFRT 长步——升序逐列界翻转吸收剩余违反，
    /// 直至某候选以 ≤ 其界幅的步长入基；全部吸收仍不够 ⇒ 候选耗尽 None；
    /// 无界上界列（range=∞）恒可入基（长步必终于入基）。
    #[test]
    fn dual_ratio_bfrt_bound_flip_absorbs_remaining() {
        // rem = 10：候选 0 界幅 3（翻吸收 rem→7）、候选 1 界幅 4（→3）、
        // 候选 2 无上界（range=∞ ⇒ 恒可入基）⇒ 入基 t = 3，翻转集 = [0,2)。
        let cands = [
            (0.1, 2usize, 1.0, -1.0, 3.0),
            (0.2, 5usize, 1.0, -1.0, 4.0),
            (0.3, 8usize, 1.0, -1.0, f64::INFINITY),
        ];
        assert_eq!(dual_ratio_bfrt(&cands, 10.0), Some((2, 3.0)));
        // 候选 2 也只有界幅 1：翻转三个候选共吸收 8，rem 仍余 2 ⇒ 耗尽
        // None（对偶无界 ⇒ 调用方回退交 Phase-1）。
        let finite = [
            (0.1, 2usize, 1.0, -1.0, 3.0),
            (0.2, 5usize, 1.0, -1.0, 4.0),
            (0.3, 8usize, 1.0, -1.0, 1.0),
        ];
        assert_eq!(dual_ratio_bfrt(&finite, 10.0), None);
        // 空候选表：候选空场景（主循环同款回退语义）。
        assert_eq!(dual_ratio_bfrt(&[], 1.0), None);
        // 负步长（rem 与 move_per_unit 反号的数值防御）⇒ None 且未翻转。
        let one = [(0.1, 2usize, 1.0, -1.0, 10.0)];
        assert_eq!(dual_ratio_bfrt(&one, -1.0), None);
    }

    /// G4（S10 测试 4 对齐）：候选空（该基上对偶无界 ⇒ 原始不可行）的
    /// GiveUp 回退——主循环不下结论、回退链交 Phase-1 判 Infeasible，
    /// 与稠密引擎 resolve 状态恒同。
    #[test]
    fn dual_main_candidate_empty_yields_infeasible_via_phase1() {
        // min x s.t. x+y ≥ 2、x−y ≤ 0、x,y ∈ [0,1]：根最优 x=y=1；
        // 收紧 y.ub=0.5 后 x+y ≤ 1.5 < 2 ⇒ 不可行。
        let vars = vec![var("x", 0.0, 1.0), var("y", 0.0, 1.0)];
        let cons = vec![
            (vec![(0usize, 1.0), (1usize, 1.0)], RelOp::Ge, 2.0),
            (vec![(0usize, 1.0), (1usize, -1.0)], RelOp::Le, 0.0),
        ];
        let p = prob(Sense::Minimize, vars, vec![(0usize, 1.0)], cons);
        let tol = Tolerances::default();
        let root = Bounds::from_problem(&p);
        let mut b = root.clone();
        b.ub[1] = 0.5;
        let mut sp = SparseSimplex::new(&p, tol).unwrap();
        let s0 = sp.solve(&root);
        assert_eq!(format!("{:?}", s0.status), "Optimal");
        let warm = sp.resolve(&b);
        let mut d = DenseSimplex::new(&p, tol).unwrap();
        let _ = d.solve(&root);
        let dense = d.resolve(&b);
        assert_eq!(
            format!("{:?}", warm.status),
            format!("{:?}", dense.status),
            "与稠密引擎 resolve 状态分叉"
        );
        assert_eq!(format!("{:?}", warm.status), "Infeasible");
    }

    /// G5（S10 测试 5 对齐）：预算耗尽诚实映射——预算参数化内核以 0 预算
    /// 直接驱动（单测规模不可稳定构造真实耗尽，0 预算 = 耗尽的极限形态：
    /// 一次主元未做即回退），回退链 dual_repair 接住后仍给出与稠密引擎
    /// 一致的正确目标；任何一层都不谎报状态。
    #[test]
    fn dual_main_budget_exhaustion_maps_honestly() {
        let p = load_instance(&milp_instance("knapsack-10")).unwrap();
        let tol = Tolerances::default();
        let root = Bounds::from_problem(&p);
        let mut b = root.clone();
        b.ub[0] = 0.0;
        // 稠密引擎参照。
        let mut d = DenseSimplex::new(&p, tol).unwrap();
        let _ = d.solve(&root);
        let dense = d.resolve(&b);
        // 手动复现 resolve 的 (1)(2) 前奏（换界 + snap），以便注入预算。
        let mut sp = SparseSimplex::new(&p, tol).unwrap();
        let _ = sp.solve(&root);
        sp.start_call();
        std::mem::swap(&mut sp.bounds, &mut sp.prev_bounds);
        sp.bounds.lb.clone_from(&b.lb);
        sp.bounds.ub.clone_from(&b.ub);
        sp.snap_nonbasics();
        assert!(
            sp.basics_violated() && !sp.phase1_pending && !sp.artificials_in_basis(),
            "前奏后应处于「对偶修复适用」状态"
        );
        // 预算 0 ⇒ 主循环立即耗尽回退（未做任何主元/翻转，状态干净）。
        assert_eq!(sp.dual_main_loop(0), DualResolveOutcome::Fallthrough);
        // 第二层 dual_repair 接住。
        assert_eq!(sp.dual_repair(), DualRepair::Repaired);
        // Phase-2 收尾 + 解构造：目标与稠密引擎 resolve 恒同。
        let st = sp.run_two_phase(false);
        let sol = sp.build_solution(st);
        assert_eq!(format!("{:?}", sol.status), format!("{:?}", dense.status));
        assert_eq!(format!("{:?}", sol.status), "Optimal");
        assert!(
            (sol.objective - dense.objective).abs() <= 1e-7 * (1.0 + dense.objective.abs()),
            "{} vs {}",
            sol.objective,
            dense.objective
        );
    }

    /// G6（S10 测试 6 对齐）：对偶退化步（θ ≤ tol.optimality 的零代价
    /// 收紧）连续 32 步 ⇒ 切入对偶 Bland；非退化步清零退化计数（Bland
    /// 本身保持——快路径无释放机制，持续退化由补丁 A 兜底）。
    #[test]
    fn dual_degen_guard_bland_threshold_and_reset() {
        let tol_opt = Tolerances::default().optimality;
        let mut g = DualDegenGuard::default();
        for _ in 0..31 {
            assert!(!g.on_step(0.0, tol_opt));
        }
        assert!(!g.bland, "31 步退化不应触发 Bland");
        assert!(!g.on_step(0.0, tol_opt));
        assert!(g.bland, "第 32 步连续退化应触发对偶 Bland");
        // 非退化步：退化计数清零，Bland 保持介入。
        assert!(!g.on_step(1.0, tol_opt));
        assert_eq!(g.degen, 0);
        assert!(g.bland);
        // θ 恰等于 tol.optimality 也是退化步（闭区间语义）。
        let mut g2 = DualDegenGuard::default();
        for _ in 0..32 {
            assert!(!g2.on_step(tol_opt, tol_opt));
        }
        assert!(g2.bland, "θ == tol.optimality 应计为退化步");
    }

    /// G7（S10 测试 7 对齐）：对偶主循环阈值钉值——S10 钉值 32（Bland
    /// 介入）与补丁 A 的 GiveUp 阈值 64（2× 同源）；预算为 A/B 调参后
    /// 的固定快路径量级 128（setcover-400x800 慢磨 2.96× 回归的修参，
    /// 见常量 doc），失败熔断 8 次。任何一侧改动使口径分叉时，本测试
    /// 提醒人工复核对齐（允许有意分叉——但须改钉值并记因）。
    #[test]
    fn dual_main_constants_pinned() {
        assert_eq!(DUAL_BLAND_DEGEN_LIMIT, 32); // S10 钉值。
        assert_eq!(DUAL_BLAND_GIVEUP_DEGEN, 64); // 复盘补丁 A：2× 同源。
        assert_eq!(DUAL_MAIN_BUDGET_BASE, 64); // A/B 调参：m + 64。
        assert_eq!(DUAL_MAIN_BUDGET_PER_ROW, 1); // 每行 1 步。
        assert_eq!(DUAL_MAIN_FAIL_FUSE, 8); // A/B 调参：连续失败熔断。
        assert_eq!(DUAL_MAIN_FUSE_MIN_CALLS, 32); // 比率熔断最小样本。
        assert_eq!(DUAL_MAIN_FAIL_RATE_PCT, 40); // 比率熔断失败率线（%）。
        assert_eq!(DUAL_MAIN_CELL_CEILING, 350_000); // 规模闸（270k 赢/480k 亏之间）。
    }

    /// G8（复盘补丁 A 用例）：对偶 Bland 介入后**仍连续**退化达 64 步 ⇒
    /// on_step 返回 true（直接 GiveUp 回退老路径，不等迭代 cap）；期间
    /// 任何一步非退化都会重置「Bland 后计数」——只有连续打转才触发。
    #[test]
    fn dual_degen_guard_gives_up_after_bland_still_degenerate() {
        let tol_opt = Tolerances::default().optimality;
        let mut g = DualDegenGuard::default();
        // 32 步退化 ⇒ Bland。
        for _ in 0..32 {
            assert!(!g.on_step(0.0, tol_opt));
        }
        assert!(g.bland);
        // Bland 后连续退化：63 步不触发，第 64 步触发 GiveUp。
        for k in 0..63 {
            assert!(!g.on_step(0.0, tol_opt), "Bland 后第 {} 步退化不应触发", k + 1);
        }
        assert!(
            g.on_step(0.0, tol_opt),
            "Bland 后连续第 64 步退化应触发补丁 A 的 GiveUp"
        );
        // 「连续」语义：Bland 后退化 63 步 → 一步非退化（重置）→ 再退化
        // 63 步不触发。
        let mut g2 = DualDegenGuard::default();
        for _ in 0..32 {
            assert!(!g2.on_step(0.0, tol_opt));
        }
        for _ in 0..63 {
            assert!(!g2.on_step(0.0, tol_opt));
        }
        assert!(!g2.on_step(1.0, tol_opt)); // 非退化步重置 after_bland。
        for _ in 0..63 {
            assert!(!g2.on_step(0.0, tol_opt));
        }
        // 重置后重新累计：第 64 步退化触发。
        assert!(g2.on_step(0.0, tol_opt));
    }

    /// G9（终审高危项回归锚）：FactorizationBroken 通道的安全支柱——
    /// refactor 失败留下的两类损坏因子化（perm 清空 / m=0）在
    /// rebuild_initial_basis 后必须被完全覆写恢复（±1 对角基 refactor 必
    /// 成），随后任何 btran 消费点不 panic。主循环两个 refactor-Err 出口
    /// 的真实触发需数值奇异基（单测无法稳定构造，如终审所注「未实机
    /// 复现」），此处直接钉住接线所依赖的恢复语义。
    #[test]
    fn rebuild_recovers_from_broken_factorization() {
        let p = load_instance(&milp_instance("knapsack-10")).unwrap();
        let tol = Tolerances::default();
        let bounds = Bounds::from_problem(&p);
        let mut sp = SparseSimplex::new(&p, tol).unwrap();
        let _ = sp.solve(&bounds);
        let m = sp.basis.len();
        // 损坏态一：x_B 重算 NaN 路径（refactor Err）——lu = default（m=0）。
        sp.lu = LuFactors::default();
        sp.rebuild_initial_basis();
        let y = sp.btran_full(&vec![1.0; m]);
        assert_eq!(y.len(), m, "m=0 损坏态 rebuild 后 btran 须恢复正常维度");
        // 损坏态二：factorize_into 奇异 Err 路径——perm 清空、f.m 保留
        // 旧维度（btran 的 y[self.perm[i]] 空索引 panic 形态）。
        sp.lu.perm.clear();
        sp.lu.u_diag.iter_mut().for_each(|v| *v = 0.0);
        sp.rebuild_initial_basis();
        let y2 = sp.btran_full(&vec![1.0; m]);
        assert_eq!(y2.len(), m);
        assert!(y2.iter().all(|v| v.is_finite()), "覆写恢复后对偶变量须有限");
        // rebuild 后引擎仍可正常求解（损坏不粘滞）。
        let mut b = bounds.clone();
        b.ub[0] = 0.0;
        let sol = sp.resolve(&b);
        assert_eq!(format!("{:?}", sol.status), "Optimal");
    }

    /// G10（终审低危项补缺）：BFRT 界翻转的**真实路径**集成覆盖——手工
    /// 构造 LP 使主循环一次迭代内「最小 θ 候选界幅不足以入基（翻转吸收）
    /// + 下一候选入基清零剩余违反」，断言 dual_main_flips/pivots 均 ≥ 1、
    /// 修复后 Optimal 且目标与手工推导值及稠密引擎 resolve 恒同。
    ///
    /// 构造（行归一化 ÷2）：max x0+x1−0.1·x2 s.t. x0−x1−2·x2 ≤ 1，
    /// x0∈[0,3]、x1∈[0,10]、x2∈[0,0.5]。根最优 = 界角点 (3,10,0)、
    /// 基 {s1}（s1=4），全部非基贴界且对偶可行（x2 正成本保贴下界）。
    /// 收紧 x1.ub→0：snap 传播 Δ=−10 ⇒ s1 = −1（below，rem=1）。
    /// 候选（θ,列号）升序：x2（θ=0.1，range=0.5，t=1>0.5 ⇒ **翻转**，
    /// rem→0.5）→ x0（θ=2，range=3，t=1≤3 ⇒ 入基）。修复后最优
    /// x0=2、x1=0、x2=0.5，目标 = 2−0.05 = 1.95。
    #[test]
    fn dual_main_bfrt_flip_real_path() {
        let vars = vec![
            var("x0", 0.0, 3.0),
            var("x1", 0.0, 10.0),
            var("x2", 0.0, 0.5),
        ];
        let cons = vec![(
            vec![(0usize, 1.0), (1usize, -1.0), (2usize, -2.0)],
            RelOp::Le,
            1.0,
        )];
        let p = prob(Sense::Maximize, vars, vec![(0, 1.0), (1, 1.0), (2, -0.1)], cons);
        let tol = Tolerances::default();
        let root = Bounds::from_problem(&p);
        let mut sp = SparseSimplex::new(&p, tol).unwrap();
        let s0 = sp.solve(&root);
        assert_eq!(format!("{:?}", s0.status), "Optimal");
        assert!(
            (s0.objective - 13.0).abs() <= 1e-9,
            "根最优须为界角点 13（构造前提），实得 {}",
            s0.objective
        );
        // 收紧 x1 上界 ⇒ 主循环 BFRT 翻转 + 入基。
        let mut b = root.clone();
        b.ub[1] = 0.0;
        let warm = sp.resolve(&b);
        let mut d = DenseSimplex::new(&p, tol).unwrap();
        let _ = d.solve(&root);
        let dense = d.resolve(&b);
        assert_eq!(
            format!("{:?}", warm.status),
            format!("{:?}", dense.status),
            "与稠密引擎 resolve 状态分叉"
        );
        assert_eq!(format!("{:?}", warm.status), "Optimal");
        assert!(
            (warm.objective - dense.objective).abs() <= 1e-7 * (1.0 + dense.objective.abs()),
            "目标恒同失败：{} vs {}",
            warm.objective,
            dense.objective
        );
        assert!(
            (warm.objective - 1.95).abs() <= 1e-7,
            "手工推导最优 1.95，实得 {}",
            warm.objective
        );
        assert!(sp.dual_main_pivots >= 1, "对偶主循环主元未被覆盖");
        assert!(sp.dual_main_flips >= 1, "BFRT 界翻转未被真实路径覆盖");
    }

    /// G11（A/B 调参诊断，重负载 #[ignore]；release 显式跑）：root 割循环 +
    /// 真 dive（每步固定当前 LP 解的最大分数基变量——B&B 分支的真实形态，
    /// 保证触发修复）序列上，对偶主循环的每次调用平均主元数（pivots/
    /// calls）与 Fallthrough 失败占比——budget/熔断阈值的定档依据。
    #[test]
    #[ignore]
    fn dual_main_diagnostics_setcover() {
        let t0 = std::time::Instant::now();
        let p = load_instance(&bench_instance("setcover-300x600")).unwrap();
        let tol = Tolerances::default();
        let bounds = Bounds::from_problem(&p);
        let int_vars = p.integer_vars();
        let mut sp = SparseSimplex::new(&p, tol).unwrap();
        let _ = sp.solve(&bounds);
        let (cuts, _st, _) = crate::cuts::cut_loop_recorded(
            &mut sp,
            &bounds,
            &p,
            &int_vars,
            crate::cuts::CutOptions::default(),
            &tol,
            || SparseSimplex::new(&p, tol),
        )
        .unwrap();
        eprintln!(
            "root+cuts({cuts}): calls={} pivots={} flips={} fall={} ({:.1}s)",
            sp.dual_main_calls,
            sp.dual_main_pivots,
            sp.dual_main_flips,
            sp.dual_main_fallthroughs,
            t0.elapsed().as_secs_f64()
        );
        // 真 dive：每步把当前解的最大分数变量固定到较近整数 ⇒ 该变量
        // 多为基变量 ⇒ 界剪出基行违反 ⇒ 触发三层修复（B&B 节点形态）。
        let mut b = bounds.clone();
        let n = p.vars.len();
        for step in 0..40 {
            let x = sp.resolve(&b).x;
            let mut bj = usize::MAX;
            let mut bfrac = 0.0f64;
            for j in 0..n {
                let f = (x[j] - x[j].round()).abs();
                if f > bfrac {
                    bfrac = f;
                    bj = j;
                }
            }
            if bj == usize::MAX {
                break; // 已整数解
            }
            let v = x[bj].round().max(b.lb[bj]).min(b.ub[bj]);
            b.lb[bj] = v;
            b.ub[bj] = v;
            if step % 8 == 7 {
                eprintln!(
                    "  step {step}: calls={} pivots={} flips={} fall={}",
                    sp.dual_main_calls,
                    sp.dual_main_pivots,
                    sp.dual_main_flips,
                    sp.dual_main_fallthroughs
                );
            }
        }
        let calls = sp.dual_main_calls.max(1);
        eprintln!(
            "TOTAL: calls={} pivots={} (avg {:.0}/call) flips={} fall={} ({:.1}s)",
            sp.dual_main_calls,
            sp.dual_main_pivots,
            sp.dual_main_pivots as f64 / calls as f64,
            sp.dual_main_flips,
            sp.dual_main_fallthroughs,
            t0.elapsed().as_secs_f64()
        );
    }

    // ---- F. S3 集成面：resolve 热启动 / add_row 割 / 快照作用域 ----

    /// 界收紧（分支形态）下的热启动 resolve：目标与冷解、与稠密引擎
    /// resolve 三方一致。
    ///
    /// 重锚注（S10 三层化：对偶主循环 → dual_repair → 整基重建）：本测试
    /// 口径 = **目标值与状态恒同**——热启动迭代序列与解顶点允许与冷解/
    /// 稠密引擎分叉（多解 LP 下 x 不作对账项）。
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
    /// 重锚注（S10 三层化）：口径 = 目标值与状态恒同（顶点/迭代序列允许
    /// 分叉，同 `resolve_warm_matches_cold_and_dense_after_tightening`）。
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
