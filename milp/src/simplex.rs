//! 有界变量两阶段表格单纯形 LP 求解器（WP2 实现）。
//!
//! # 算法结构
//! - 显式 dense 表格 + 有界变量机制（非基变量位于 lb 或 ub，记 `at_upper`），
//!   分支定界只改 [`Bounds`] 不改约束矩阵，热启动代价最低。
//! - 稀疏混合表格（轴1，**规模阈值后置**）：表格单元数 m×ncols 超过
//!   [`RC_SPARSE_CELL_THRESHOLD`]（2^18）才启用每行非零列索引 `row_nz`
//!   簿记（与 dense 表格并存）——主元消元归并游标同步、检验数按行索引
//!   散射、add_row 消基稀疏迭代、割分离走
//!   [`Simplex::tableau_row_terms`] 稀疏读取。阈内**完整走旧路径**
//!   （do_pivot_legacy：O(ncols) 计数扫描 + 纯索引消元；稠密行主序 axpy；
//!   κ 稠密消基）——小算例热路径零新增指令、成本与行为逐字节回基线。
//!   行密度阈值沿用 nnz*4 ≥ ncols：超阈行标记稠密（空索引，恒走 zip
//!   向量化循环、零簿记成本），稀疏行走索引表；fill-in 归并插入、恰零
//!   墓碑惰性删除（每 64 主元压缩）、每行密度回退。启用态语义为值级/
//!   决策级 == 等价（跳零是精确无操作）。回归门：debug 构建周期性全量
//!   重建比对。binv 保持稠密——诚实边界：这是混合表格不是修正单纯形，
//!   目标域是 1e3–1e4 行入场券，不是百万变量。
//! - 内部统一 min（Maximize 目标系数取负，报告时还原并加 `obj_offset`）；
//!   构建期每行除以 max|a_ij| 归一（割行例外，不缩放）。**行永不翻转**，
//!   因此 [`SlackCol::sign`] 只有两个合法来源（Le=+1、Ge=−1、割行恒 +1），
//!   对偶与解的还原无符号歧义。
//! - 初始基（统一残差规则）：设初始非基点 x⁰ = 各结构变量取 lb，
//!   r_i = b̂_i − â_i·x⁰：
//!   * Le 行（â·x + s = b̂）：r_i ≥ −feas ⇒ 松弛列入基（值 r_i）；
//!     否则加人工列（σ=−1，值 −r_i>0）入基，松弛列非基贴 0。
//!   * Ge 行（â·x − s = b̂）：r_i ≤ +feas ⇒ 剩余列入基（值 −r_i）；
//!     否则加人工列（σ=+1，值 r_i>0）入基，剩余列非基贴 0。
//!   * Eq 行：恒加人工列（σ = r_i≥0 ? +1 : −1，值 |r_i|）入基。
//!   人工变量绝不限于等式行：正右端 ≥ 行（如 x≥3 而 x 停在 lb=0）的剩余量
//!   恒为负，取反行救不了它；旧的“b_i<0 整行取反”规则已废除。
//! - Phase-1 目标 min Σ人工变量；最优值 > tol.feasibility ⇒ Infeasible；
//!   人工列离基后永久禁止再进基（本实现从一开始就不允许人工列入基，
//!   比契约更强且同样安全）。残基人工变量值必须为 0；Phase-2 用原目标。
//! - 近零行（max|a_ij| ≤ tol.zero）不入表，按原尺度 rhs 判定：违反则置
//!   `infeasible_flag`（后续 solve/resolve 直接 Infeasible），满足则永久剔除
//!   （duals 恒 0，`row_of_constraint` 返回 None）。
//!
//! # 定价与终止
//! - 双侧合格性：贴下界且 rc < −tol.optimality（增向入基）或 贴上界且
//!   rc > +tol.optimality（减向入基）。贴上界列的减向入基是二进制/有限上界
//!   变量的常态操作，缺它会在非最优点误报 Optimal。
//! - Dantzig 阶段得分 = 贴下界列 (−rc)、贴上界列 (+rc)，取最大，并列取最小
//!   列下标。反循环：自上次目标改善起连续退化迭代 > 50 次或本次调用总迭代
//!   超过 3×(m+n) 后进入 Bland 模式（本次调用内永久生效）：入基取合格列最小
//!   下标；出基侧全部并列候选（被阻挡基变量按其列下标、入基变量自身界翻转
//!   按入基列下标）一律取最小下标；|主元|最大规则仅 Dantzig 阶段使用。
//! - 比值测试同时考虑基变量上下界与入基变量自身界翻转；主元 |p| ≥ tol.pivot
//!   才允许；入基列无任何阻挡 ⇒ Unbounded。检验数由当前表格在三入口
//!   （price / dual_feasible / dual_repair 每轮迭代顶部）批量重算
//!   （compute_rc_batch：值级等价于逐列现算，不增量维护，无定价漂移）；
//!   每 200 次主元用原始数据重算 B⁻¹b；NaN/Inf 立即按数值失败处理。迭代
//!   上限 max(10_000, 20×(m+n))（割循环内经 set_iter_limit_hint 收紧为
//!   max(3_000, 12×(m+n))，R8）⇒ IterationLimit（诚实上报，不冒充
//!   Optimal）。Bland 释放（R8）：停滞触发的 Bland 在非退化步（目标真
//!   改善 ⇒ 已脱离触发退化区）+ 驻留满 8 步后回 Dantzig 定价（每次调用
//!   最多 64 次重入，预算耗尽后 Bland 永久驻留——无循环定理的终止路径
//!   恒在场）——Bland 的最小下标入基在宽表上单步进展极小，纯 Bland 爬行
//!   实测打满迭代上限（setcover 类根 LP）。主循环每迭代零分配：c_b/检验
//!   数/比值候选/消元系数/主元行非零下标/界交换全部走容量复用暂存字段。
//!
//! # 表格约定（实现细节）
//! - 列空间 = [结构列 | 松弛列 | 人工列]，人工列不外露（`column_layout`
//!   只含前两段）。行空间 = 未剔除的约束行 + 割行。
//! - `matrix` 维护 B⁻¹·A；`binv` 显式维护 B⁻¹（m×m，随主元同步更新），
//!   用于周期性 B⁻¹b 重算与对偶计算。`rhs[i]` 直接存第 i 行基变量的当前值
//!   （非基变量位于界上时的取值），非基变量变界后按 x_B += −Δx_j·(B⁻¹A_j)
//!   传播。
//! - 归一化记 `row_scale[i] = 1/max|a_ij|`（割行恒 1.0）：内部对偶 ŷ = ∂z/∂b̂
//!   还原为 ∂z/∂rhs 时乘以该系数（â = a/g, b̂ = b/g ⇒ ∂z/∂b = ŷ/g），
//!   满足契约“内部对偶 × scale 还原”且数学正确。结构变量的 x 无需还原
//!   （行缩放不改变变量空间）。
//!
//! # 热启动与数值失败的报告通道
//! - `resolve`：非基列贴界双向检查（收紧越界 snap / 回溯恢复后严格内点
//!   snap 到较近界、等距按 rc 符号选侧），Δ 按 B⁻¹A_j 传播；随后检查基
//!   变量（含残基人工变量）是否越界，或存在 Phase-1 欠账（上次 Phase-1
//!   未清理即收场——不可行/上限/数值失败、人工变量带正值留基，或 add_row
//!   新建人工列；其基值即便被 snap 传播推到容差内，行系数仍活跃，Phase-2
//!   会解到该行被零成本人工变量吸收的错误 LP），任一成立则重建初始基
//!   （含按需重建人工列）重跑 Phase-1 修复，否则直接 Phase-2。修复采用
//!   “重建初始基”而非复合目标 Phase-1：两者终态等价（正确的 LP 状态），
//!   且重建复用冷启动全部已测代码路径，正确性最有保障。
//! - `solve`/`resolve` 签名不可失败（返回 [`LpStatus`] 而非 `Result`），因此
//!   契约中的 `Err(Numerical)`（主元过滤后无候选、NaN/Inf 等）在库内部
//!   映射为 `LpStatus::IterationLimit`（objective=NaN，绝不冒充 Optimal），
//!   失败原因记入私有字段 `numerical_error` 供调试。

use crate::error::{MilpError, Result};
use crate::model::{Problem, RelOp, Sense, Tolerances};

/// 连续退化迭代阈值：超过后切入 Bland 模式（反循环）。
const DEGENERATE_STREAK_LIMIT: u32 = 50;
/// Bland 释放的最小驻留步数（R8）：Bland 驻留至少该步数且出现非退化步
///（目标真改善，证明已脱离当前退化区）后回到 Dantzig 定价。Bland 的
/// 最小下标入基在宽表（setcover 类 ncols ≫ m）上单步进展极小——纯退化
/// 停滞触发 Bland 后永久驻留会在 3×(m+n) 后爬行至迭代上限（setcover
/// 300x600/400x800 根 LP 实测：Phase-2 以 bland=true 打满 18000/24000
/// 上限、零 incumbent 传播）。释放规则只作用于**停滞触发**的 Bland；
/// 总迭代触发（iter_call > 3×(m+n)）的 Bland 永久生效（原有保底语义
/// 不变），且每次调用最多释放 [`BLAND_REARM_MAX`] 次——重入预算耗尽后
/// Bland 永久驻留，保证 Bland 无循环定理的终止路径恒在场。
const BLAND_DWELL: u32 = 8;
/// 每次 solve/resolve 调用的 Bland 释放重入预算（R8）：耗尽后停滞触发的
/// Bland 也永久驻留（终止保底）。
const BLAND_REARM_MAX: u32 = 64;
/// 主元计数刷新间隔：每累计本次主元次数重算一次 B⁻¹b 抑制漂移。
const REFRESH_INTERVAL: u64 = 200;
/// Bland 触发的总迭代因子：本次调用总迭代 > BLAND_TOTAL_FACTOR×(m+n) 触发。
const BLAND_TOTAL_FACTOR: u64 = 3;
/// 比值测试并列带宽的绝对上限系数：band = min(feas·(1+|t_min|), feas·CAP)。
/// 相对项在 |t_min| ≈ 1e6 时膨胀到 ≈1e-3，会把真正更小的界翻转/出基候选
/// （差距 1e-4 量级）吞进并列集合、按 |主元| 破键输给假并列的 Leave——
/// 步长越过真最小比值，基变量越界至 band·|α|（2026-10 第三轮在合成算例
/// 冷解上实测观察到，见 test_ratio_test_band_large_magnitude）。取小后：
/// |t_min| ≤ CAP−1 时与旧公式逐值一致（零扰动），超出后带宽封顶
/// feas·CAP（默认 1e-6，仍高于 |t|=1e6 处 ulp 级浮点噪声 ~2e-10 三个
/// 数量级以上，真并列不被拆散）；收窄并列集合只会让选择更贴近最小比值，
/// 真 t_min 候选恒在集合内，不引入误剪。
const RATIO_BAND_CAP: f64 = 1e3;

// —— 稀疏混合表格（轴1：行级非零索引簿记，规模阈值后置）——

/// 行非零索引表项的墓碑位（最高位）：置位 = 该列曾非零、现已**恰为 0.0**，
/// 惰性待压缩（迭代跳过、不计存活数）；置位/复活只动高 1 位，低 31 位仍为
/// 列下标、保持升序。
const NZ_TOMB: u32 = 1 << 31;
/// 行非零索引表项的列下标掩码（低 31 位）。列数 ≥ 2^31 时禁用稀疏簿记
///（单行 dense 表即需 16 GiB，实际不可达；防御性回退全稠密）。
const NZ_MASK: u32 = !NZ_TOMB;
/// 行索引簿记维护间隔（主元数，仅启用态）：压缩墓碑 + debug 全量重建比对。
const NZ_COMPACT_INTERVAL: u64 = 64;
/// 轴1 整体启用阈值（表格单元数 m×ncols）：**≤ 阈值时 do_pivot / add_row /
/// compute_rc_batch 完整走旧路径**（O(ncols) 计数扫描、纯索引消元、稠密
/// 行主序 axpy、κ 稠密消基——热路径零新增指令，小算例成本与行为逐字节回
/// 基线）；超阈才启用 row_nz 簿记 / 归并游标消元 / 检验数散射 / 稀疏读表。
/// 取 2^18 = 262144：bench 13 算例中表格单元最多的 lp-wide-200x300 为
/// 300×500 = 150000（次大 setcover-100x150 含根割 ≈1.07e5、mknap-5x300
/// 含根割 ≤9e4），全部留驻旧路径；lp-sparse-4000x800（≈3.8e6）、
/// setcover-400x800（≈4.8e5）与拉伸算例 lp-wide-1000x2000（≈6e6）、
/// setcover-300x600（≈7.2e5）超阈启用。启用态的语义契约是值级/决策级
/// == 等价（跳零是精确无操作，至多翻转个别 ±0.0 符号位），不再承诺逐字节。
const RC_SPARSE_CELL_THRESHOLD: u64 = 262_144;

/// LP 求解状态。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LpStatus {
    /// 已证明最优。
    Optimal,
    /// 不可行。
    Infeasible,
    /// 无界。
    Unbounded,
    /// 触发迭代上限（诚实上报，不冒充 Optimal）。
    IterationLimit,
}

/// 结构变量当前上下界（分支定界只改这个，不改 Problem）。
#[derive(Debug, Clone)]
pub struct Bounds {
    /// 各变量下界。
    pub lb: Vec<f64>,
    /// 各变量上界（允许 f64::INFINITY）。
    pub ub: Vec<f64>,
}

impl Bounds {
    /// 从问题默认界构建。
    pub fn from_problem(problem: &Problem) -> Self {
        Bounds {
            lb: problem.vars.iter().map(|v| v.lb).collect(),
            ub: problem.vars.iter().map(|v| v.ub).collect(),
        }
    }

    /// 只收紧一侧（下界）。
    ///
    /// 返回 `false` 表示 `var` 越界或新区间为空（调用方按剪枝处理），不 panic。
    /// 传入比现值更松的界不会放宽（“只收紧”语义）。
    pub fn tighten_lb(&mut self, var: usize, lb: f64) -> bool {
        if var >= self.lb.len() || var >= self.ub.len() {
            return false;
        }
        // f64::max 忽略 NaN 输入，天然防御非有限实参。
        let new_lb = lb.max(self.lb[var]);
        if new_lb > self.ub[var] {
            return false;
        }
        self.lb[var] = new_lb;
        true
    }

    /// 只收紧一侧（上界）。
    ///
    /// 返回 `false` 表示 `var` 越界或新区间为空（调用方按剪枝处理），不 panic。
    /// 传入比现值更松的界不会放宽（“只收紧”语义）。
    pub fn tighten_ub(&mut self, var: usize, ub: f64) -> bool {
        if var >= self.lb.len() || var >= self.ub.len() {
            return false;
        }
        let new_ub = ub.min(self.ub[var]);
        if new_ub < self.lb[var] {
            return false;
        }
        self.ub[var] = new_ub;
        true
    }
}

/// 行来源：`Constraint(i)` = `problem.constraints()[i]`；
/// `Cut(k)` = 第 k 条 [`Simplex::add_row`] 追加的割（k 即 add_row 返回值）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowOrigin {
    /// 原始约束行（下标为 problem.constraints() 下标）。
    Constraint(usize),
    /// add_row 追加的割行（k 从 0 连续递增）。
    Cut(usize),
}

/// 松弛列信息。
///
/// `sign` 仅有两来源（行永不翻转）：+1 = Le 行松弛（â·x+s=b̂），
/// −1 = Ge 行剩余（â·x−s=b̂）；割行恒 +1。
#[derive(Debug, Clone, Copy)]
pub struct SlackCol {
    /// 所属行来源。
    pub origin: RowOrigin,
    /// +1 = Le 行松弛；−1 = Ge 行剩余。
    pub sign: f64,
}

/// 标准形列布局：列空间 = [结构列 | 松弛列]，人工列不外露。
#[derive(Debug, Clone)]
pub struct ColumnLayout {
    /// 结构列数（下标 0..num_structural 对应 problem 变量）。
    pub num_structural: usize,
    /// 各松弛/剩余列（下标 num_structural.. 起顺序对应）。
    pub slack: Vec<SlackCol>,
}

/// 非基列及其所贴界。
#[derive(Debug, Clone, Copy)]
pub struct NonbasicCol {
    /// 列下标（结构+松弛列空间）。
    pub col: usize,
    /// 是否贴上界（false = 贴下界）。
    pub at_upper: bool,
}

/// LP 求解结果。
#[derive(Debug, Clone)]
pub struct LpSolution {
    /// 求解状态。
    pub status: LpStatus,
    /// 目标值：原始 sense + obj_offset 口径；非 Optimal 时为 NaN。
    pub objective: f64,
    /// 结构变量取值（长度 = num_vars）。
    pub x: Vec<f64>,
    /// 对偶值：长度 = num_constraints，按原始约束顺序；
    /// 语义 = 最优处 ∂(内部min目标)/∂rhs_i，由内部对偶 × scale_i 还原（无翻转
    /// 符号）；被剔除的近零行恒 0；仅供诊断与测试，割不依赖它。
    pub duals: Vec<f64>,
}

/// 有界变量两阶段表格单纯形求解器。
///
/// 私有状态：dense 表格（`matrix` = B⁻¹A）、显式 B⁻¹（`binv`）、基列、
/// 非基贴界状态、列布局、行缩放系数、人工列记账、迭代计数、Bland 模式
/// 标记、`infeasible_flag` 等。
pub struct Simplex {
    /// 容差集中配置。
    tol: Tolerances,
    /// 结构变量数 n。
    num_structural: usize,
    /// 原始约束数（duals 输出长度）。
    num_constraints: usize,
    /// 目标是否为 Maximize（内部统一 min，报告时取负还原）。
    sense_max: bool,
    /// 目标常数偏移（报告口径）。
    obj_offset: f64,
    /// 内部 min 目标系数（长度 n；Maximize 已取负）。
    c_internal: Vec<f64>,
    /// 列布局（结构列 + 松弛列；人工列不外露）。
    layout: ColumnLayout,
    /// 当前人工列数（人工列总在 [结构|松弛] 段之后连续排列）。
    num_art: usize,
    /// dense 表格：matrix[i][j] = (B⁻¹A)_{ij}，行 = 基位置，列 = 全列空间。
    matrix: Vec<Vec<f64>>,
    /// 行级非零列索引簿记（轴1 稀疏混合表格），与 `matrix` 并存，**仅在
    /// [`Self::nz_enabled`]（表格单元数超 [`RC_SPARSE_CELL_THRESHOLD`]）时
    /// 维护**：`row_nz[i]` = 第 i 行非零列下标的升序 u32 表。**空表 = 稠密
    /// 标记**：该行恒走稠密 zip 循环、不承担簿记义务（割/GMI 行在大 LP 上
    /// 近乎全稠密——空标记使其零簿记成本、零索引内存，密度回退是保底）。
    /// 稀疏标记行满足不变量：存活项集合 == 该行精确非零集（值 != 0.0，墓碑
    /// 项不算）且存活数低于密度阈（nnz*4 < ncols，与 do_pivot 判据同源）。
    /// 墓碑项（最高位置位）= 已恰为零的列：惰性删除，每
    /// [`NZ_COMPACT_INTERVAL`] 主元压缩一次（控制删除 memmove 成本）。
    /// 禁用态恒为全空（整基重建时跨阈翻转会全量重建/清空）。
    row_nz: Vec<Vec<u32>>,
    /// 各行索引内的墓碑数（存活数 = len − dead，O(1) 密度判定用）。
    nz_dead: Vec<u32>,
    /// 距上次行索引压缩/debug 比对的主元数（仅启用态计数）。
    pivots_since_nz: u64,
    /// 轴1 启用位：m×ncols > [`RC_SPARSE_CELL_THRESHOLD`] 时为真（在
    /// rebuild_initial_basis / add_row 收口处重判——这两处是行数/列数唯一
    /// 变化点）。假 = 全旧路径（do_pivot_legacy 等），真 = 簿记路径。
    nz_enabled: bool,
    /// 显式基逆：binv[i][k] = (B⁻¹)_{ik}。
    binv: Vec<Vec<f64>>,
    /// rhs[i] = 第 i 行基变量当前值（非基变量贴界时的取值口径）。
    rhs: Vec<f64>,
    /// 各行基列下标（全列空间）。
    basis: Vec<usize>,
    /// 列是否在基中（全列空间长度）。
    in_basis: Vec<bool>,
    /// 非基列是否贴上界（全列空间长度；基列时无意义）。
    at_upper: Vec<bool>,
    /// 每表格行的归一化结构系数（稀疏，按列升序）。
    row_terms: Vec<Vec<(usize, f64)>>,
    /// 每表格行归一化右端 b̂。
    row_rhs: Vec<f64>,
    /// 每表格行关系运算符。
    row_op: Vec<RelOp>,
    /// 每表格行自己的松弛列下标（Eq 行为 None）。
    row_slack: Vec<Option<usize>>,
    /// 每表格行对偶还原系数（约束行 = 1/max|a_ij|，割行 = 1.0）。
    row_scale: Vec<f64>,
    /// 原始约束下标 → 表格行（被剔除的近零行为 None）。
    constraint_row: Vec<Option<usize>>,
    /// 割序号 → 表格行（割行永不剔除；防御性近零割行为 None）。
    cut_row: Vec<Option<usize>>,
    /// 最近一次 solve/resolve 使用的有效界。
    bounds: Bounds,
    /// 是否已完成过一次 solve/resolve（resolve 首次调用走冷启动路径）。
    has_run: bool,
    /// 当前表格存在未经完整 Phase-1 清理的人工变量（上次 Phase-1 以不可行
    /// /中断收场且未调 cleanup，或 add_row 新建人工列）。此时人工变量可能
    /// 带正值留基且行系数活跃：其基值即便被热启动 snap 传播推到容差内，
    /// 状态也不可用（Phase-2 会在等式/不等式行被零成本人工变量吸收的错误
    /// LP 上求解），下一次 resolve 必须重跑 Phase-1 修复。
    phase1_pending: bool,
    /// 构建期近零行违反 ⇒ 整个实例不可行。
    infeasible_flag: bool,
    /// 数值失败原因（签名无错误通道，仅供诊断）。
    #[allow(dead_code)] // 保留字段：调试时可观察数值失败原因
    numerical_error: Option<String>,
    /// 累计迭代数（跨调用）。
    iterations_total: u64,
    /// 本次调用迭代数（上限与 Bland 触发口径）。
    iter_call: u64,
    /// 自上次有效改善起的连续退化迭代数。
    degenerate_streak: u32,
    /// Bland 模式。R8 起分两种：总迭代触发的永久 Bland（`bland_permanent`，
    /// 本次调用内不释放）与停滞触发的可释放 Bland（非退化步 + 驻留满
    /// [`BLAND_DWELL`] 后回 Dantzig，重入预算 [`BLAND_REARM_MAX`]）。
    bland: bool,
    /// 总迭代触发（iter_call > 3×(m+n)）的永久 Bland 置位。
    bland_permanent: bool,
    /// Bland 已驻留步数（释放的驻留门）。
    bland_dwell: u32,
    /// 单次调用迭代上限提示（R8，割循环收紧用；None = 默认公式）。
    /// 割循环每轮 resolve 是暖启动（对偶修复快路径，收敛轮实测 ~1.5k 迭代），
    /// 新割行引发的 Bland 爬行若放任到默认上限 max(10k, 20(m+n)) 会白烧
    /// （setcover-100x150 第 5 轮实测 10000 迭代打满后整树回退）。收紧后
    /// 失败轮更早中止，由 cuts.rs 的回滚语义保留已提交割。
    iter_limit_hint: Option<u64>,
    /// 本次调用已释放 Bland 的重入次数。
    bland_rearms: u32,
    /// 距上次 B⁻¹b 重算的主元数。
    pivots_since_refresh: u64,
    // —— 热路径暂存（容量跨迭代/跨调用复用，消除每迭代分配）——
    /// 批量检验数暂存（全列空间长度；三入口各自入口重算，见
    /// [`Simplex::compute_rc_batch`] 的新鲜度契约）。
    rc_scratch: Vec<f64>,
    /// 基成本 c_b 暂存（行数长度；由 compute_rc_batch 填充）。
    cb_scratch: Vec<f64>,
    /// 比值测试候选暂存（clear+push 容量复用，见 ratio_test）。
    cand_scratch: Vec<(f64, usize, bool, f64, usize)>,
    /// do_pivot 消元系数暂存（行数长度，覆盖写，见 do_pivot）。
    pivot_f: Vec<f64>,
    /// do_pivot 主元行非零下标暂存（行缩放/消元只遍历非零列，见 do_pivot）。
    pivot_nz: Vec<usize>,
    /// resolve 界交换暂存（swap 后保存“上次有效界”，clone_from 复用容量）。
    prev_bounds: Bounds,
}

/// 定价结果。
enum Pricing {
    /// 合格入基列与方向（δ=+1 增向 / −1 减向）。
    Candidate(usize, f64),
    /// 无任何合格列（最优性终止条件满足）。
    Optimal,
    /// 数值失败（检验数非有限等）。
    Failed(String),
}

/// 比值测试结果。
enum Step {
    /// 基变量在 `row` 出基（步长 t，`hit_upper` 表示其贴上界）。
    Leave { row: usize, t: f64, hit_upper: bool },
    /// 入基变量自身从所贴界翻转到对端界（步长 t）。
    Flip { t: f64 },
    /// 无任何阻挡（含自身界翻转不可用）⇒ 无界。
    NoBlock,
    /// 存在被主元容差过滤掉的阻挡候选（数值不可靠）。
    FilteredBlock(String),
}

/// 对偶修复结果。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DualRepair {
    /// 已恢复原始可行（基保持最优口径，随后 Phase-2 收尾应为零迭代）。
    Repaired,
    /// 不适用 / 超出预算 / 数值异常 ⇒ 调用方回退整基重建（语义与既有
    /// 冷路径完全一致，快路径只加速、绝不改变结果）。
    GiveUp,
}

impl Simplex {
    /// 由 LP 松弛构建（整数性忽略）。
    ///
    /// 结构变量 lb=−∞ 返回 `Err(InvalidModel)`（.milp 格式保证下界有限；
    /// ub=+∞ 合法）；系数/右端非有限数同样报错。
    pub fn new(problem: &Problem, tol: Tolerances) -> Result<Self> {
        let n = problem.vars.len();
        // —— 校验（全部通过后再构建，避免半初始化状态）——
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

        // —— 构建 ——
        let sense_max = problem.sense == Sense::Maximize;
        let sign = if sense_max { -1.0 } else { 1.0 };
        let mut c_internal = vec![0.0; n];
        for &(j, c) in &problem.obj.terms {
            c_internal[j] = sign * c;
        }

        let mut s = Simplex {
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
            num_art: 0,
            matrix: Vec::new(),
            row_nz: Vec::new(),
            nz_dead: Vec::new(),
            pivots_since_nz: 0,
            nz_enabled: false,
            binv: Vec::new(),
            rhs: Vec::new(),
            basis: Vec::new(),
            in_basis: Vec::new(),
            at_upper: Vec::new(),
            row_terms: Vec::new(),
            row_rhs: Vec::new(),
            row_op: Vec::new(),
            row_slack: Vec::new(),
            row_scale: Vec::new(),
            constraint_row: vec![None; problem.constraints.len()],
            cut_row: Vec::new(),
            bounds: Bounds {
                lb: vec![0.0; n],
                ub: vec![f64::INFINITY; n],
            },
            has_run: false,
            phase1_pending: false,
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
            pivots_since_refresh: 0,
            rc_scratch: Vec::new(),
            cb_scratch: Vec::new(),
            cand_scratch: Vec::new(),
            pivot_f: Vec::new(),
            pivot_nz: Vec::new(),
            prev_bounds: Bounds {
                lb: Vec::new(),
                ub: Vec::new(),
            },
        };

        for (ci, con) in problem.constraints.iter().enumerate() {
            // 行归一化系数 g = max|a_ij|（仅结构系数；Eq/Le/Ge 行均如此）。
            let g = con
                .expr
                .terms
                .iter()
                .map(|&(_, a)| a.abs())
                .fold(0.0_f64, f64::max);
            if g <= s.tol.zero {
                // 近零行：不入表，按原尺度 rhs 判定。
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

        // 分配表格行结构（行数以保留下来的行为准），再预置一致的初始表格
        // （默认界），保证 add_row 可在 solve 之前调用。
        let m = s.row_terms.len();
        s.matrix = vec![Vec::new(); m];
        // 行索引先置空（禁用态稠密标记），紧随其后的 rebuild_initial_basis
        // 按表格规模判定启用并重建（见 nz_update_enabled/nz_rebuild_all）。
        s.row_nz = vec![Vec::new(); m];
        s.nz_dead = vec![0; m];
        s.rhs = vec![0.0; m];
        s.basis = vec![usize::MAX; m];
        s.binv = vec![vec![0.0; m]; m];
        let ncols0 = n + s.layout.slack.len();
        s.in_basis = vec![false; ncols0];
        s.at_upper = vec![false; ncols0];
        s.bounds = Bounds::from_problem(problem);
        s.rebuild_initial_basis();
        Ok(s)
    }

    /// 两阶段冷启动求解（当前 bounds；完全忽略此前内部状态）。
    pub fn solve(&mut self, bounds: &Bounds) -> LpSolution {
        self.start_call();
        if self.infeasible_flag {
            return self.build_solution(LpStatus::Infeasible);
        }
        if let Some(st) = self.validate_bounds_arg(bounds) {
            return self.build_solution(st);
        }
        // 界写入：clone_from 容量复用（替代整体 clone 的分配）。
        self.bounds.lb.clone_from(&bounds.lb);
        self.bounds.ub.clone_from(&bounds.ub);
        self.rebuild_initial_basis();
        let st = self.run_two_phase(true);
        self.has_run = true;
        self.build_solution(st)
    }

    /// 热启动求解：非基双向 snap+传播 → 必要时 Phase-1 修复 → Phase-2。
    pub fn resolve(&mut self, bounds: &Bounds) -> LpSolution {
        self.start_call();
        if self.infeasible_flag {
            return self.build_solution(LpStatus::Infeasible);
        }
        if let Some(st) = self.validate_bounds_arg(bounds) {
            return self.build_solution(st);
        }
        if !self.has_run {
            // 尚无可用热状态：退化为冷启动。
            self.bounds.lb.clone_from(&bounds.lb);
            self.bounds.ub.clone_from(&bounds.ub);
            self.rebuild_initial_basis();
            let st = self.run_two_phase(true);
            self.has_run = true;
            return self.build_solution(st);
        }
        // (1)(2) 非基贴界双向检查 + 传播。
        // 界暂存 swap + clone_from：`prev_bounds` 换入“上次有效界”（容量跨
        // 调用复用，消除原 mem::replace 的双克隆分配），新界 clone_from 写
        // 入 self.bounds；snap_nonbasics 内部从 self.prev_bounds 读取 prev
        // （语义与原局部变量一致）。
        std::mem::swap(&mut self.bounds, &mut self.prev_bounds);
        self.bounds.lb.clone_from(&bounds.lb);
        self.bounds.ub.clone_from(&bounds.ub);
        self.snap_nonbasics();
        // (3) 基变量越界检查 + Phase-1 欠账检查：任一成立则修复。修复有
        //     两条路：
        //     a) **对偶修复快路径**（大规模主路径）：仅界变化引起（无
        //        Phase-1 欠账、基中无人工列）且静态非基对偶可行 ⇒ 有界
        //        变量对偶单纯形几次旋转恢复原始可行（基与成本未动 ⇒ 对偶
        //        可行性保持；出基变量落到其被违反的界）。任何异常一律
        //        GiveUp 回退 b)，快路径只加速、绝不改变语义。
        //     b) 整基重建 + Phase-1（原路径，语义不变）：欠账 = 上次
        //        Phase-1 未清理即收场（不可行/上限/数值失败）或 add_row
        //        新建了人工列：此时残基人工变量行系数活跃，其基值即便被
        //        (1)(2) 的 snap 传播恰好推到 ≤ tol.feasibility，Phase-2 仍
        //        会解到该行被零成本人工变量吸收的错误 LP 上（终检降级
        //        IterationLimit 或谎报 Optimal）。
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

    /// 追加割行（terms 为结构变量空间稀疏行，不做缩放，scale=1）；
    /// 以 Le 语义加松弛列 s∈[0,+∞)、sign=+1、`RowOrigin::Cut(k)`；
    /// 当前解可能违反新行，由随后 resolve 修复；返回割序号 k（0 起连续递增）。
    ///
    /// 实现说明：契约按 Le 割描述（cuts.rs 恒产 Le 割）；为防御性完整起见，
    /// op=Ge 按剩余列（sign=−1）、op=Eq 按无松弛+人工列处理，语义与构建期
    /// 行规则一致。割行永不剔除；若传入近零行（正常割纪律不会发生）则不入
    /// 表并按原尺度 rhs 判定是否置 `infeasible_flag`。
    pub fn add_row(&mut self, terms: &[(usize, f64)], op: RelOp, rhs: f64) -> Result<usize> {
        let n = self.num_structural;
        if !rhs.is_finite() {
            return Err(MilpError::InvalidModel("add_row 右端非有限数".to_string()));
        }
        // 归并排序 + 丢弃精确零系数（防御；cuts.rs 保证输入已规范）。
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
            // 防御性近零割行：不入表（契约：割行数值纪律由 cuts.rs 保证）。
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
        let ncols_old = n + self.layout.slack.len() + self.num_art;
        // 基列位置映射。
        let mut pos_in_basis = vec![usize::MAX; ncols_old];
        for i in 0..m {
            pos_in_basis[self.basis[i]] = i;
        }
        // w：新行对当前基列的系数（仅结构基列非零）。
        let mut w = vec![0.0; m];
        for &(j, c) in &merged {
            let p = pos_in_basis[j];
            if p != usize::MAX {
                w[p] = c;
            }
        }
        // κ_j = ν_j − Σ_i w_i·M[i][j]（对新行的全旧列消基；轴1 启用态对
        // 稀疏标记行按行索引迭代——跳过 M[i][j]==0 的格是精确无操作减法，
        // 逐列减法序列与稠密循环同序；κ 零元为 +0.0，减 ±0 仍得 +0，
        // 逐字节一致。禁用态走旧稠密循环）。
        let mut kappa = vec![0.0; ncols_old];
        for &(j, c) in &merged {
            kappa[j] = c;
        }
        for i in 0..m {
            if w[i] != 0.0 {
                let wi = w[i];
                let mrow = &self.matrix[i];
                if self.nz_enabled && !self.row_nz[i].is_empty() {
                    for &e in self.row_nz[i].iter() {
                        if e & NZ_TOMB != 0 {
                            continue;
                        }
                        let j = (e & NZ_MASK) as usize;
                        kappa[j] -= wi * mrow[j];
                    }
                } else {
                    for j in 0..ncols_old {
                        kappa[j] -= wi * mrow[j];
                    }
                }
            }
        }
        for i in 0..m {
            kappa[self.basis[i]] = 0.0; // 数值卫生：基列位置恒 0
        }
        // 残差 r_new = rhs − ν·x_cur（基列按当前值，非基列按所贴界）。
        let mut nu_x = 0.0;
        for i in 0..m {
            if w[i] != 0.0 {
                nu_x += w[i] * self.rhs[i];
            }
        }
        for &(j, c) in &merged {
            if pos_in_basis[j] == usize::MAX {
                nu_x += c * self.nonbasic_value(j);
            }
        }
        let r_new = rhs - nu_x;

        // 基变量选择（与构建期残差规则一致；行不缩放）。
        let slack_sign = match op {
            RelOp::Le => 1.0,
            RelOp::Ge => -1.0,
            RelOp::Eq => 0.0, // Eq 无松弛列
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

        // 新列：松弛列（Eq 无）与人工列（非基时才建）。割松弛列**插入**在
        // 静态段（结构+松弛）末尾、人工列段之前——若按旧口径直接取
        // `n + slack.len()`（= 人工列段首下标），会与存量人工列撞号：新割
        // 行基列与旧人工列共址、该列在表格中带旧行残留系数而非单位向量，
        // 表格不再是 B⁻¹A（此前靠后续整基重建偶然兜住）。插入方式：旧行在
        // 该位补 0（新松弛列对旧行的 B⁻¹A 系数恒 0——B⁻¹ 未动、新列在旧
        // 行上系数为 0），人工列段整体右移一位，基中人工列下标同步 +1；
        // layout 记账保持『slack 第 t 项 ↔ 列 n+t』的顺序对应。
        let ncs_old = n + self.layout.slack.len();
        let has_slack = op != RelOp::Eq;
        let has_art = !use_slack_basic;
        let slack_col: Option<usize> = if has_slack {
            let c = ncs_old;
            for row in self.matrix.iter_mut() {
                row.insert(c, 0.0);
            }
            // 轴1（启用态）：既有行索引同步右移一位（列 ≥ c 的存活/墓碑项
            // +1；新松弛列对旧行恒 0，不入索引。墓碑在高 1 位，+1 不进位）。
            // 禁用态索引恒空，跳过（旧路径零新增指令）。
            if self.nz_enabled {
                for idx in self.row_nz.iter_mut() {
                    for e in idx.iter_mut() {
                        if (*e & NZ_MASK) as usize >= c {
                            *e += 1;
                        }
                    }
                }
            }
            self.in_basis.insert(c, false);
            self.at_upper.insert(c, false);
            for b in self.basis.iter_mut() {
                if *b >= c {
                    *b += 1; // 人工列段整体右移一位（结构/松弛列下标不变）
                }
            }
            self.layout.slack.push(SlackCol {
                origin: RowOrigin::Cut(k),
                sign: slack_sign,
            });
            Some(c)
        } else {
            None
        };
        let art_col: Option<usize> = if has_art {
            Some(n + self.layout.slack.len() + self.num_art)
        } else {
            None
        };
        let basic_col = if use_slack_basic {
            slack_col.unwrap_or(usize::MAX)
        } else {
            art_col.unwrap_or(usize::MAX)
        };
        debug_assert!(basic_col != usize::MAX);

        let ncols_new = ncols_old + usize::from(has_slack) + usize::from(has_art);
        // 旧行补零（松弛列插入已扩一位；此处再补新人工列槽）。
        for row in self.matrix.iter_mut() {
            row.resize(ncols_new, 0.0);
        }
        self.in_basis.resize(ncols_new, false);
        self.at_upper.resize(ncols_new, false);
        // 新表格行 = κ/σ（旧列；人工列段右移一位）| 松弛列系数 slack_sign/σ
        // | 人工列系数 1。
        let shift = usize::from(has_slack);
        let mut newrow = vec![0.0; ncols_new];
        for j in 0..ncols_old {
            let dest = if j < ncs_old { j } else { j + shift };
            newrow[dest] = kappa[j] / sigma_basic;
        }
        if let Some(c) = slack_col {
            newrow[c] = slack_sign / sigma_basic;
        }
        if let Some(c) = art_col {
            newrow[c] = art_sign / sigma_basic; // = 1.0
        }

        // binv 扩一行一列：新行 = [−σ·(wᵀB⁻¹), σ]。
        let mut newbinv = vec![0.0; m + 1];
        for kk in 0..m {
            let mut acc = 0.0;
            for i in 0..m {
                if w[i] != 0.0 {
                    acc += w[i] * self.binv[i][kk];
                }
            }
            newbinv[kk] = -sigma_basic * acc;
        }
        newbinv[m] = sigma_basic;
        for row in self.binv.iter_mut() {
            row.push(0.0);
        }
        self.binv.push(newbinv);

        // 追加表格行与记账。
        let row = self.matrix.len();
        // 轴1（启用态）：新行索引（GMI 割行在大 LP 上近乎全稠密 → 多直接
        // 落稠密标记）。禁用态推空表（稠密标记语义）。
        let mut newidx: Vec<u32> = Vec::new();
        if self.nz_enabled && ncols_new < NZ_MASK as usize {
            for (j, &v) in newrow.iter().enumerate() {
                if v != 0.0 {
                    newidx.push(j as u32);
                }
            }
            if newidx.len() * 4 >= ncols_new {
                newidx.clear();
            }
        }
        self.matrix.push(newrow);
        self.row_nz.push(newidx);
        self.nz_dead.push(0);
        self.rhs.push(r_new / sigma_basic);
        self.basis.push(basic_col);
        self.in_basis[basic_col] = true;
        self.num_art += usize::from(has_art);
        if has_art {
            // 新人工列入基即欠账：其行系数活跃（等式/违反行残差吸收器），
            // 基值即使 ≤ tol.feasibility 也不允许下一次 resolve 跳过 Phase-1
            // 直接跑 Phase-2。Le 割（契约路径）基值恒 > tol.feasibility，
            // basics_violated 本可捕获；此处对 Eq/Ge 防御路径一视同仁。
            self.phase1_pending = true;
        }
        self.row_terms.push(merged);
        self.row_rhs.push(rhs);
        self.row_op.push(op);
        self.row_slack.push(slack_col);
        self.row_scale.push(1.0); // 割行不缩放（契约例外）
        self.cut_row.push(Some(row));
        // 轴1 规模门收口：行/列刚增长，重判启用位；跨阈翻转（禁用→启用）时
        // 一次性全量重建索引（含新行；此前推入的是空表）。add_row 只增维度，
        // 启用态不会翻回禁用。
        if !self.nz_enabled {
            let m_now = self.basis.len();
            let ncols_now = self.ncols_static() + self.num_art;
            if (m_now as u64) * (ncols_now as u64) > RC_SPARSE_CELL_THRESHOLD {
                self.nz_enabled = true;
                self.nz_rebuild_all();
            }
        }
        Ok(k)
    }

    /// 第 `row` 行（表格行下标）的 (B⁻¹Ã 系数[结构+松弛列], 该行基变量当前值)；
    /// 仅供基变量为结构列的行（割平面用）。
    pub fn tableau_row(&self, row: usize) -> (Vec<f64>, f64) {
        if row >= self.matrix.len() {
            return (Vec::new(), f64::NAN);
        }
        let ncols = self.num_structural + self.layout.slack.len();
        (self.matrix[row][..ncols].to_vec(), self.rhs[row])
    }

    /// 第 `row` 行的稀疏读取接口（轴1）：(非零 `(列下标, 值)` 项列表
    /// [结构+松弛列段、按列升序]，该行基变量当前值)。
    ///
    /// 与 [`Simplex::tableau_row`] 的语义关系：稠密接口中系数恰为 0.0 的列
    /// 在此省略——消费方（GMI 分离）对 α==0 的列产出的 π 恒为 0、无任何
    /// 系数/常数贡献，行为等价；NaN/Inf 非零必入表，异常行检测不受影响。
    /// 启用态稀疏标记行直接走行索引（跳过墓碑），禁用态/稠密标记行按行值
    /// 扫描收集。越界行返回 `(空, NaN)`（与稠密接口同口径）。
    pub fn tableau_row_terms(&self, row: usize) -> (Vec<(usize, f64)>, f64) {
        if row >= self.matrix.len() {
            return (Vec::new(), f64::NAN);
        }
        let ncols = self.num_structural + self.layout.slack.len();
        let mut terms: Vec<(usize, f64)> = Vec::new();
        if self.nz_enabled && !self.row_nz[row].is_empty() {
            for &e in self.row_nz[row].iter() {
                if e & NZ_TOMB != 0 {
                    continue;
                }
                let j = (e & NZ_MASK) as usize;
                if j < ncols {
                    terms.push((j, self.matrix[row][j]));
                }
            }
        } else {
            for (j, &v) in self.matrix[row][..ncols].iter().enumerate() {
                if v != 0.0 {
                    terms.push((j, v));
                }
            }
        }
        (terms, self.rhs[row])
    }

    /// 原始约束下标 → 表格行下标；被剔除的近零行返回 `None`。
    pub fn row_of_constraint(&self, constraint_index: usize) -> Option<usize> {
        self.constraint_row.get(constraint_index).copied().flatten()
    }

    /// 割序号 → 表格行下标（割行永不被剔除）。
    pub fn row_of_cut(&self, cut_index: usize) -> Option<usize> {
        self.cut_row.get(cut_index).copied().flatten()
    }

    /// 各行基列下标（结构+松弛列空间）。
    pub fn basis_columns(&self) -> &[usize] {
        &self.basis
    }

    /// 当前全部非基列及其所贴界（不含人工列）。
    pub fn nonbasic_columns(&self) -> Vec<NonbasicCol> {
        let ncols = self.num_structural + self.layout.slack.len();
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

    /// 表格行数（含割行，不含被剔除的近零行）。
    pub fn num_rows(&self) -> usize {
        self.basis.len()
    }

    /// 累计单纯形迭代次数。
    pub fn iterations(&self) -> u64 {
        self.iterations_total
    }

    // ===================== 内部：列与成本 =====================

    /// 静态列数（结构 + 松弛，不含人工列）。
    fn ncols_static(&self) -> usize {
        self.num_structural + self.layout.slack.len()
    }

    /// 列的当前上下界（结构列用当前 bounds；松弛/人工列恒 [0,+∞)）。
    fn col_bounds(&self, j: usize) -> (f64, f64) {
        if j < self.num_structural {
            (self.bounds.lb[j], self.bounds.ub[j])
        } else {
            (0.0, f64::INFINITY)
        }
    }

    /// 非基列当前取值（按所贴界；调用方保证 j 非基）。
    fn nonbasic_value(&self, j: usize) -> f64 {
        let (lb, ub) = self.col_bounds(j);
        if self.at_upper[j] {
            ub
        } else {
            lb
        }
    }

    /// 列在指定阶段的目标成本（phase 1 = min Σ人工；phase 2 = 内部 min 原目标）。
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

    /// 第 j 列的 Phase-2 检验数（用当前表格现算，无增量漂移）。
    fn rc_of(&self, j: usize) -> f64 {
        let m = self.basis.len();
        let mut rc = self.cost_of(j, 2);
        for i in 0..m {
            let cb = self.cost_of(self.basis[i], 2);
            if cb != 0.0 {
                rc -= cb * self.matrix[i][j];
            }
        }
        rc
    }

    /// 批量计算全列空间在指定阶段的检验数（写入 `rc_scratch`，同时填充
    /// 基成本暂存 `cb_scratch`）。
    ///
    /// # 值级等价
    /// 先按列填 `cost_of(j, phase)`，再按基行 i 升序、仅 c_b[i]≠0 行做
    /// 行主序 axpy——每列的减法序列与 price 原内联循环 / [`rc_of`] 的逐列
    /// 计算完全同序，f64 结果逐值相等（跳过 c_b==0 行是精确无操作：减 0
    /// 不改变任何数值，含 ±0 符号位）。
    ///
    /// 轴1（稀疏散射，规模门后置）：禁用态（m×ncols ≤
    /// [`RC_SPARSE_CELL_THRESHOLD`]）走下方纯稠密旧路径不动——小算例行为
    /// 逐字节不变；启用态稀疏标记行改按行非零索引散射 `rc[j] -= cb·row[i][j]`
    /// （跳过该行 matrix[i][j]==0 的格）：每列的减法序列仍与稠密路径同序
    /// （行 i 升序），跳过项是精确无操作（减 ±0 至多把 −0.0 翻成 +0.0，
    /// f64 == 不可分辨），结果保持 == 相等；稠密标记行（空索引）回退 zip
    /// 循环（与旧稠密路径逐字节同序）。
    ///
    /// # 新鲜度契约（三入口各自在入口无条件重算）
    /// 禁止只由 price 刷新后他处复用——do_pivot、phase 切换、add_row、
    /// 整基重建都会使缓存陈旧。三入口：
    /// (a) [`price`] 入口；(b) [`dual_feasible`] 入口（resolve 入口距上次
    /// price 可能隔 phase 切换/add_row/整基重建）；(c) [`dual_repair`]
    /// 外层循环每轮迭代顶部（上一轮 do_pivot 已改表，下一轮列扫描前必须
    /// 重算；扫描“只对过滤列读取”是省工而非新鲜度保证，二者不得混淆）。
    /// v1 不做表格版本号脏标记，三入口无条件重算。
    fn compute_rc_batch(&mut self, phase: u8) {
        let m = self.basis.len();
        let ncols = self.ncols_static() + self.num_art;
        self.cb_scratch.clear();
        self.cb_scratch.resize(m, 0.0);
        for i in 0..m {
            self.cb_scratch[i] = self.cost_of(self.basis[i], phase);
        }
        self.rc_scratch.clear();
        self.rc_scratch.resize(ncols, 0.0);
        for j in 0..ncols {
            self.rc_scratch[j] = self.cost_of(j, phase);
        }
        if self.nz_enabled {
            // 启用态走外置散射路径（轴1）：稀疏标记行按行索引散射，稠密标记
            // 行回退 zip。独立函数 + inline(never)：两路共存于同一函数体会
            // 扰动旧路径的代码生成（寄存器分配/内联预算，实测 knapsack-2000
            // +11% 墙钟），外置后阈内代码生成与原实现一致。
            self.compute_rc_batch_scatter();
            return;
        }
        for i in 0..m {
            let cb = self.cb_scratch[i];
            if cb != 0.0 {
                let row = &self.matrix[i];
                let rc = &mut self.rc_scratch;
                for j in 0..ncols {
                    rc[j] -= cb * row[j];
                }
            }
        }
    }

    /// compute_rc_batch 的启用态散射路径（轴1）：稀疏标记行按 row_nz 散射
    /// `rc[j] -= cb·row[i][j]`（跳过墓碑与零格——精确无操作，== 等价），
    /// 稠密标记行（空索引）回退 zip（与旧稠密路径逐字节同序）。
    #[inline(never)]
    fn compute_rc_batch_scatter(&mut self) {
        let m = self.basis.len();
        let ncols = self.ncols_static() + self.num_art;
        for i in 0..m {
            let cb = self.cb_scratch[i];
            if cb == 0.0 {
                continue;
            }
            let row = &self.matrix[i];
            let rc = &mut self.rc_scratch;
            if self.row_nz[i].is_empty() {
                for j in 0..ncols {
                    rc[j] -= cb * row[j];
                }
            } else {
                for &e in self.row_nz[i].iter() {
                    if e & NZ_TOMB != 0 {
                        continue;
                    }
                    let j = (e & NZ_MASK) as usize;
                    rc[j] -= cb * row[j];
                }
            }
        }
    }

    // ===================== 内部：初始基重建 =====================

    /// 按统一残差规则重建初始基：丢弃旧人工列，非基点取 x⁰ = 结构变量贴
    /// 当前下界；每行按 Le/Ge/Eq 规则选松弛或带符号人工列入基。
    /// 初始基列为 ±单位列，B⁻¹ = diag(±1)，表格行 = σ_i×归一化行（无需
    /// Gauss-Jordan，O(m·ncols)）。
    fn rebuild_initial_basis(&mut self) {
        let n = self.num_structural;
        let ncs = self.ncols_static();
        let m = self.matrix.len();
        // 截断旧人工列并清零静态区残留（重入 rebuild 时旧表格值必须归零，
        // 否则跨行松弛列位置会保留旧值污染新初始表格）。
        for row in self.matrix.iter_mut() {
            row.truncate(ncs);
            for v in row.iter_mut() {
                *v = 0.0;
            }
        }
        self.in_basis.truncate(ncs);
        self.at_upper.truncate(ncs);
        self.num_art = 0;
        // 旧人工列连同欠账一并作废：新初始基要么无人工列，要么由紧随其后的
        // run_two_phase(true) 跑完整 Phase-1 清理。
        self.phase1_pending = false;
        self.in_basis.iter_mut().for_each(|b| *b = false);
        // 冷启动非基点 x⁰ = 结构变量贴下界：at_upper 一律复位为 false。
        self.at_upper.iter_mut().for_each(|u| *u = false);

        // 残差 r_i = b̂_i − â_i·x⁰。
        let mut r = vec![0.0; m];
        for i in 0..m {
            let mut ax = 0.0;
            for &(j, a) in &self.row_terms[i] {
                ax += a * self.bounds.lb[j];
            }
            r[i] = self.row_rhs[i] - ax;
        }
        // 选基 + 符号。
        let mut sigma = vec![1.0; m];
        let mut need_art: Vec<Option<f64>> = vec![None; m];
        for i in 0..m {
            match self.row_op[i] {
                RelOp::Le => {
                    if r[i] >= -self.tol.feasibility {
                        sigma[i] = 1.0;
                    } else {
                        sigma[i] = -1.0;
                        need_art[i] = Some(-1.0);
                    }
                }
                RelOp::Ge => {
                    if r[i] <= self.tol.feasibility {
                        sigma[i] = -1.0;
                    } else {
                        sigma[i] = 1.0;
                        need_art[i] = Some(1.0);
                    }
                }
                RelOp::Eq => {
                    let s = if r[i] >= 0.0 { 1.0 } else { -1.0 };
                    sigma[i] = s;
                    need_art[i] = Some(s);
                }
            }
        }
        // 追加人工列（按行序连续编号）。
        let art_rows: Vec<usize> = (0..m).filter(|&i| need_art[i].is_some()).collect();
        let ncols_new = ncs + art_rows.len();
        for row in self.matrix.iter_mut() {
            row.resize(ncols_new, 0.0);
        }
        self.in_basis.resize(ncols_new, false);
        self.at_upper.resize(ncols_new, false);
        self.num_art = art_rows.len();
        let mut art_col_of_row = vec![usize::MAX; m];
        for (t, &i) in art_rows.iter().enumerate() {
            art_col_of_row[i] = ncs + t;
        }
        // 填表。
        for i in 0..m {
            let s = sigma[i];
            let row = &mut self.matrix[i];
            for &(j, a) in &self.row_terms[i] {
                row[j] = s * a;
            }
            if let Some(c) = self.row_slack[i] {
                row[c] = s * self.layout.slack[c - n].sign;
            }
            let basic = match need_art[i] {
                Some(_) => art_col_of_row[i],
                None => self.row_slack[i].unwrap_or(art_col_of_row[i]),
            };
            row[basic] = 1.0; // σ_i·σ_basic = 1（基列规范型）
            self.rhs[i] = s * r[i];
            self.basis[i] = basic;
            self.in_basis[basic] = true;
        }
        // B⁻¹ = diag(σ)。
        for i in 0..m {
            for k in 0..m {
                self.binv[i][k] = 0.0;
            }
            self.binv[i][i] = sigma[i];
        }
        // 轴1（④）规模门 + 索引重建：超阈（m×ncols > RC_SPARSE_CELL_THRESHOLD）
        // 启用时从填表结果一次性重建全部行非零索引（含人工列段）；阈内禁用
        // （行索引清空 = 稠密标记语义，旧路径零簿记——大表缩回阈内时同步
        // 释放索引内存）。整基重建是行数/列数变化点之一，启用位在此重判。
        self.nz_update_enabled();
        if self.nz_enabled {
            self.nz_rebuild_all();
        } else {
            for idx in self.row_nz.iter_mut() {
                idx.clear();
            }
            for d in self.nz_dead.iter_mut() {
                *d = 0;
            }
        }
    }

    // ===================== 内部：稀疏索引簿记（轴1，规模门后置） =====================

    /// 重判轴1 启用位：m×ncols > [`RC_SPARSE_CELL_THRESHOLD`] 为真。
    /// 只在维度变化点调用（rebuild_initial_basis / add_row 收口）。
    fn nz_update_enabled(&mut self) {
        let m = self.basis.len();
        let ncols = self.ncols_static() + self.num_art;
        self.nz_enabled = (m as u64) * (ncols as u64) > RC_SPARSE_CELL_THRESHOLD;
    }

    /// 行 i 索引的存活（非墓碑）项数 = 该行精确非零数（稀疏标记行不变量）。
    fn nz_alive(&self, i: usize) -> usize {
        self.row_nz[i].len() - self.nz_dead[i] as usize
    }

    /// 就地重建第 i 行索引：按 matrix[i] 现值收集非零列（写指针恒 ≤ 扫描
    /// 位置，无 memmove），存活数超密度阈（nnz*4 ≥ ncols）转稠密标记。
    /// 调用方保证 row_nz[i]/nz_dead[i] 已分配。
    fn nz_rebuild_row(&mut self, i: usize) {
        if self.matrix[i].len() >= NZ_MASK as usize {
            // 列数超出 31 位索引上限（实际不可达）：防御性禁用稀疏簿记。
            self.row_nz[i].clear();
            self.nz_dead[i] = 0;
            return;
        }
        let ncols = self.matrix[i].len();
        let sparse;
        {
            let idx = &mut self.row_nz[i];
            let row = &self.matrix[i];
            let mut w = 0;
            for (j, &v) in row.iter().enumerate() {
                if v != 0.0 {
                    if w < idx.len() {
                        idx[w] = j as u32;
                    } else {
                        idx.push(j as u32);
                    }
                    w += 1;
                }
            }
            idx.truncate(w);
            sparse = idx.len() * 4 < ncols;
        }
        self.nz_dead[i] = 0;
        if !sparse {
            self.row_nz[i].clear(); // 密度回退：稠密标记
        }
    }

    /// 全表索引重建：所有行按现值重扫。O(m×ncols)，与 rebuild_initial_basis
    /// 的填表同阶；仅在启用态（构建收口 / add_row 跨阈翻转）调用。
    fn nz_rebuild_all(&mut self) {
        for i in 0..self.matrix.len() {
            self.nz_rebuild_row(i);
        }
    }

    /// 压缩全部行索引的墓碑（每 NZ_COMPACT_INTERVAL 主元一次，控制删除
    /// memmove 成本的惰性回收）。
    fn nz_compact(&mut self) {
        for i in 0..self.row_nz.len() {
            if self.nz_dead[i] > 0 {
                self.row_nz[i].retain(|&e| e & NZ_TOMB == 0);
                self.nz_dead[i] = 0;
            }
        }
    }

    /// 行 k 索引中列 j 若为存活项则置墓碑（do_pivot 的 matrix[k][q]=0.0
    /// 卫生写同步删 q 用；二分定位，正常路径下 q 已被消元的恰零簿记处理，
    /// 此处纯防御）。
    fn nz_tomb_alive(&mut self, k: usize, j: usize) {
        let idx = &mut self.row_nz[k];
        let mut lo = 0usize;
        let mut hi = idx.len();
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if ((idx[mid] & NZ_MASK) as usize) < j {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        if lo < idx.len() && ((idx[lo] & NZ_MASK) as usize) == j && idx[lo] & NZ_TOMB == 0 {
            idx[lo] |= NZ_TOMB;
            self.nz_dead[k] += 1;
        }
    }

    /// debug 回归门：全量重建索引表与簿记结果比对相等。稀疏标记行断言
    /// ① 存活项集合 == 该行精确非零集（索引漂移会静默污染表格——值错但
    /// 无 NaN，这是兜底断言）；② 存活数低于密度阈（密度回退义务）。
    /// 稠密标记行（空索引）不承担簿记义务，跳过；禁用态全空，天然空转。
    #[cfg(debug_assertions)]
    fn nz_debug_verify(&self) {
        for i in 0..self.matrix.len() {
            let idx = &self.row_nz[i];
            if idx.is_empty() {
                continue;
            }
            let expect: Vec<u32> = self.matrix[i]
                .iter()
                .enumerate()
                .filter(|&(_, &v)| v != 0.0)
                .map(|(j, _)| j as u32)
                .collect();
            let alive: Vec<u32> = idx.iter().copied().filter(|&e| e & NZ_TOMB == 0).collect();
            assert_eq!(
                alive, expect,
                "行 {i} 非零索引簿记与全量重建不符（索引漂移会静默污染表格）"
            );
            assert!(
                alive.len() * 4 < self.matrix[i].len(),
                "行 {i} 非零数 {} 未低于密度阈（ncols = {}），应回退稠密标记",
                alive.len(),
                self.matrix[i].len()
            );
        }
    }

    // ===================== 内部：定价 =====================

    /// 定价：双侧合格性扫描。Bland 模式取最小下标合格列；Dantzig 模式取
    /// 得分（贴下界 −rc / 贴上界 +rc）最大者、并列取最小列下标。
    /// 人工列与固定列（ub−lb ≤ feasibility）永不入基。
    fn price(&mut self, phase: u8) -> Pricing {
        let ncs = self.ncols_static();
        // 新鲜度规范 (a)：入口无条件批量重算检验数（同时消除原局部 c_b
        // vec![0.0; m] 分配，c_b 由 compute_rc_batch 写入 cb_scratch）。
        self.compute_rc_batch(phase);
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
            let rc = self.rc_scratch[j];
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

    // ===================== 内部：比值测试 =====================

    /// 有界变量比值测试：考虑全部基变量上下界阻挡与入基变量自身界翻转。
    /// 主元 |p| ≥ tol.pivot 才允许出基；并列破键：Dantzig 阶段 |主元|最大、
    /// 再最小列下标；Bland 模式一律最小列下标（基变量按其列下标、界翻转按
    /// 入基列下标）。
    fn ratio_test(&mut self, q: usize, delta: f64) -> Step {
        let m = self.basis.len();
        // (t, 行, hit_upper, |主元|, 破键列下标)；行 = usize::MAX 表示界翻转。
        // 暂存字段 clear+push：容量跨迭代复用（&self→&mut self，调用点均在
        // &mut 上下文）。
        self.cand_scratch.clear();
        let mut filtered_block = false;
        for i in 0..m {
            let a = self.matrix[i][q];
            if a == 0.0 {
                continue;
            }
            let alpha = delta * a; // 基变量值随 t 的变化为 rhs − t·alpha
            let (lb, ub) = self.col_bounds(self.basis[i]);
            let (t_i, hit_upper) = if alpha > 0.0 {
                if !lb.is_finite() {
                    continue;
                }
                ((self.rhs[i] - lb) / alpha, false)
            } else {
                if !ub.is_finite() {
                    continue;
                }
                ((ub - self.rhs[i]) / (-alpha), true)
            };
            let t_i = t_i.max(0.0);
            if a.abs() < self.tol.pivot {
                filtered_block = true; // 方向阻挡但主元过小，过滤
                continue;
            }
            self.cand_scratch.push((t_i, i, hit_upper, a.abs(), self.basis[i]));
        }
        // 入基变量自身界翻转（从所贴界移动到对端界）。
        let (lb, ub) = self.col_bounds(q);
        if ub.is_finite() {
            self.cand_scratch
                .push(((ub - lb).max(0.0), usize::MAX, false, 0.0, q));
        }
        let cands = &self.cand_scratch;
        if cands.is_empty() {
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
        for c in cands.iter() {
            if c.0 < t_min {
                t_min = c.0;
            }
        }
        // 并列带宽：相对项与绝对上限取小（RATIO_BAND_CAP 文档：大数量级下
        // 相对项膨胀会吞掉真正更小的候选；|t_min| ≤ 999 时与旧公式逐值一致）。
        let band =
            (self.tol.feasibility * (1.0 + t_min)).min(self.tol.feasibility * RATIO_BAND_CAP);
        // 选择：先最小比值（带宽内视为并列），再按模式破键。
        let mut pick: Option<&(f64, usize, bool, f64, usize)> = None;
        for c in cands.iter().filter(|c| c.0 <= t_min + band) {
            pick = match pick {
                None => Some(c),
                Some(best) => {
                    let better = if self.bland {
                        c.4 < best.4 // Bland：一律最小列下标
                    } else {
                        c.3 > best.3 || (c.3 == best.3 && c.4 < best.4) // Dantzig：|主元|最大、再最小列下标
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

    // ===================== 内部：主元与翻转 =====================

    /// 执行主元：入基列 q（方向 delta、步长 t），出基行 r。
    /// rhs 采用“基变量值”口径：非出基行按 −δt·f_k 更新，出基行 = 入基变量
    /// 新值；系数与 binv 做标准表格行变换。
    ///
    /// 热路径零分配：消元系数 f 覆盖写暂存 `pivot_f`；行缩放与消元按主元
    /// 行密度自适应（O2）：稠密（非零数 ≥ ncols/4）走连续 `zip` 消元（自动
    /// 向量化，吞吐高数倍），稀疏走非零下标 `pivot_nz`（省算、间接寻址）。
    /// 跳过/乘零列是精确无操作，仅可能翻转个别 ±0.0 的符号位（f64 == 比较
    /// 行为不变，值级/决策级等价）。matrix/binv 消元内循环经 split_at_mut
    /// 同时取得可变行 k 与只读行 r（别名安全写法，可向量化）；
    /// matrix[k][q]=0 卫生写在消元循环外保留。
    /// 执行主元：入基列 q（方向 delta、步长 t），出基行 r。
    /// rhs 采用“基变量值”口径：非出基行按 −δt·f_k 更新，出基行 = 入基变量
    /// 新值；系数与 binv 做标准表格行变换。
    ///
    /// 轴1 规模门（[`Self::nz_enabled`]）：阈内（m×ncols ≤
    /// [`RC_SPARSE_CELL_THRESHOLD`]）整体走 [`Simplex::do_pivot_legacy`]
    /// ——**逐字保留的旧实现**（O(ncols) 计数扫描 + zip/纯索引消元），热路径
    /// 零新增指令，小算例成本与行为逐字节回基线；超阈走本体的簿记路径：
    /// 密度判据沿用 nnz*4 ≥ ncols 但非零数改读行索引（O(1)），空索引 = 稠密
    /// 标记（恒走 zip、零簿记成本）。消元同时同步维护行 k 的索引：旧 0 变
    /// 非零（fill-in）在归并游标处插入（同位墓碑则复活，免 memmove）；减法
    /// 结果恰为 0.0 置墓碑（惰性删除，每 [`NZ_COMPACT_INTERVAL`] 主元压缩
    /// 一次）；存活数超密度阈即转稠密标记（每行密度回退，fill-in 吞掉收益
    /// 时的保底）。启用态语义为值级/决策级 == 等价（跳零是精确无操作，至多
    /// 翻转个别 ±0.0 符号位）。binv 保持稠密 zip 消元——诚实边界：这是混合
    /// 表格不是修正单纯形，目标域是 1e3–1e4 行入场券，不是百万变量。
    fn do_pivot(&mut self, q: usize, r: usize, delta: f64, t: f64, hit_upper: bool) {
        if self.nz_enabled {
            self.do_pivot_nz(q, r, delta, t, hit_upper);
            return;
        }
        self.do_pivot_legacy_inner(q, r, delta, t, hit_upper);
    }

    /// 启用态主元（轴1 簿记路径）：密度判据读行索引、归并游标消元簿记。
    fn do_pivot_nz(&mut self, q: usize, r: usize, delta: f64, t: f64, hit_upper: bool) {
        let m = self.basis.len();
        let ncols = self.matrix[r].len();
        let p = self.matrix[r][q];
        // 密度自适应：非零数读行索引（稀疏标记行的存活数 = 精确非零数），
        // 空索引 = 稠密标记。割行/GMI 行在大 LP 上近乎全稠密，索引式内循环
        // （间接寻址）阻碍 SIMD——实测 mknap-5x300 t1 do_pivot 占墙钟 ~70%。
        let dense_elim = self.row_nz[r].is_empty() || self.nz_alive(r) * 4 >= ncols;
        if !dense_elim {
            // 主元行非零下标直接读索引（跳过墓碑），替代旧全列扫描。
            self.pivot_nz.clear();
            for &e in self.row_nz[r].iter() {
                if e & NZ_TOMB == 0 {
                    self.pivot_nz.push((e & NZ_MASK) as usize);
                }
            }
        }
        // 行 r 缩放（稠密：连续除；稀疏：只遍历非零列，除法下溢到恰 0.0 的
        // 同步转墓碑——维持「存活项 == 精确非零集」不变量）。
        if dense_elim {
            for v in self.matrix[r].iter_mut() {
                *v /= p;
            }
        } else {
            let row = &mut self.matrix[r];
            for e in self.row_nz[r].iter_mut() {
                if *e & NZ_TOMB != 0 {
                    continue;
                }
                let j = (*e & NZ_MASK) as usize;
                row[j] /= p;
                if row[j] == 0.0 {
                    *e |= NZ_TOMB;
                    self.nz_dead[r] += 1;
                }
            }
        }
        self.matrix[r][q] = 1.0;
        for v in self.binv[r].iter_mut() {
            *v /= p;
        }
        // 出基列与入基旧值。
        let leaving = self.basis[r];
        let entering_old = if delta > 0.0 {
            self.col_bounds(q).0
        } else {
            self.col_bounds(q).1
        };
        // 消元系数 f（用未更新的 matrix[k][q]）：覆盖写暂存，无分配。
        self.pivot_f.clear();
        self.pivot_f.resize(m, 0.0);
        for k in 0..m {
            if k != r {
                self.pivot_f[k] = self.matrix[k][q];
            }
        }
        // 消元（k ≠ r；f_k 用未更新的 matrix[k][q]）。
        for k in 0..m {
            if k == r {
                continue;
            }
            let fk = self.pivot_f[k];
            if fk != 0.0 {
                // matrix 行 k（可变）与行 r（只读）：k≠r，经 split_at_mut
                // 取得不相交切片。
                let (row_k, row_r): (&mut Vec<f64>, &Vec<f64>) = if r < k {
                    let (left, right) = self.matrix.split_at_mut(k);
                    (&mut right[0], &left[r])
                } else {
                    let (left, right) = self.matrix.split_at_mut(k + 1);
                    (&mut left[k], &right[r - k - 1])
                };
                if dense_elim {
                    for (x, y) in row_k.iter_mut().zip(row_r.iter()) {
                        *x -= fk * *y;
                    }
                    // —— 行 k 索引簿记（稀疏标记行才承担；matrix 与 row_nz
                    //    为不相交字段，可与上面的行切片同时持有可变借用）。
                    //    主元行稠密：行 k 全列被触，就地重建索引（写指针恒 ≤
                    //    扫描位置，无 memmove）。
                    if !self.row_nz[k].is_empty() {
                        let idx = &mut self.row_nz[k];
                        let mut w = 0;
                        for (j, &v) in row_k.iter().enumerate() {
                            if v != 0.0 {
                                if w < idx.len() {
                                    idx[w] = j as u32;
                                } else {
                                    idx.push(j as u32);
                                }
                                w += 1;
                            }
                        }
                        idx.truncate(w);
                        self.nz_dead[k] = 0;
                    }
                } else if self.row_nz[k].is_empty() {
                    // 主元行稀疏 + 行 k 稠密标记：无簿记义务，纯索引消元。
                    for &j in self.pivot_nz.iter() {
                        row_k[j] -= fk * row_r[j];
                    }
                } else {
                    // 主元行稀疏 + 行 k 稀疏标记：归并游标同步消元——算术与
                    // 索引簿记单遍完成（双方升序）：存活项减法（结果恰零转
                    // 墓碑）、墓碑位复活（旧值精确 0.0，非零结果免插入
                    // memmove）、fill-in 游标处插入（旧 0 变非零）。
                    let idx = &mut self.row_nz[k];
                    let mut pos = 0;
                    for &j in self.pivot_nz.iter() {
                        let y = row_r[j];
                        while pos < idx.len() && ((idx[pos] & NZ_MASK) as usize) < j {
                            pos += 1;
                        }
                        if pos < idx.len() && ((idx[pos] & NZ_MASK) as usize) == j {
                            if idx[pos] & NZ_TOMB == 0 {
                                // 存活：减法 + 结果恰零转墓碑。
                                let x = row_k[j] - fk * y;
                                row_k[j] = x;
                                if x == 0.0 {
                                    idx[pos] |= NZ_TOMB;
                                    self.nz_dead[k] += 1;
                                }
                            } else {
                                // 墓碑位（旧值精确 0.0）：非零结果复活。
                                let nv = -fk * y;
                                if nv != 0.0 {
                                    idx[pos] &= !NZ_TOMB;
                                    self.nz_dead[k] -= 1;
                                    row_k[j] = nv;
                                }
                            }
                        } else {
                            // fill-in：旧 0 变非零，游标处插入（0 − fk·y ≡
                            // −fk·y，== 相等，至多 ±0 符号位差异）。
                            let nv = -fk * y;
                            if nv != 0.0 {
                                idx.insert(pos, j as u32);
                                row_k[j] = nv;
                            }
                        }
                    }
                }
                // 密度回退：行 k 稀疏标记且存活数超阈（nnz*4 ≥ ncols）转稠密
                // 标记（fill-in 吞掉收益时的保底）。
                if !self.row_nz[k].is_empty() && self.nz_alive(k) * 4 >= ncols {
                    self.row_nz[k].clear();
                    self.nz_dead[k] = 0;
                }
                let (bk, br): (&mut Vec<f64>, &Vec<f64>) = if r < k {
                    let (left, right) = self.binv.split_at_mut(k);
                    (&mut right[0], &left[r])
                } else {
                    let (left, right) = self.binv.split_at_mut(k + 1);
                    (&mut left[k], &right[r - k - 1])
                };
                for (x, y) in bk.iter_mut().zip(br.iter()) {
                    *x -= fk * *y;
                }
                self.rhs[k] -= delta * t * fk;
            }
            self.matrix[k][q] = 0.0; // 数值卫生（保留在非零列循环之外）
            // 卫生写同步删 q（轴1）：正常路径下消元已把 q 精确归零
            // （fk − fk·1 = +0.0）并墓碑；此处防御性同步（稀疏标记行）。
            if !self.row_nz[k].is_empty() {
                self.nz_tomb_alive(k, q);
            }
        }
        self.rhs[r] = entering_old + delta * t;
        self.basis[r] = q;
        self.in_basis[leaving] = false;
        self.in_basis[q] = true;
        self.at_upper[leaving] = hit_upper;
        // —— 簿记维护：每 NZ_COMPACT_INTERVAL 主元压缩墓碑；debug 构建
        //    同时做全量重建比对（索引漂移会静默污染表格——值错但无 NaN，
        //    周期断言是回归门）。
        self.pivots_since_nz = self.pivots_since_nz.saturating_add(1);
        if self.pivots_since_nz >= NZ_COMPACT_INTERVAL {
            self.pivots_since_nz = 0;
            self.nz_compact();
            #[cfg(debug_assertions)]
            self.nz_debug_verify();
        }
    }

    /// do_pivot 的阈内旧路径（轴1 规模门前置前的逐字原实现）：O(ncols) 计数
    /// 扫描 + zip/纯索引消元，无任何索引簿记。保留为独立函数以保证小算例
    /// 热路径零新增指令、行为逐字节不变（bench 13 算例回基线的锚点）。
    #[inline]
    fn do_pivot_legacy_inner(&mut self, q: usize, r: usize, delta: f64, t: f64, hit_upper: bool) {
        let m = self.basis.len();
        let ncols = self.matrix[r].len();
        let p = self.matrix[r][q];
        // 密度自适应（O2）：先计数非零（O(ncols) 无分支计数），稀疏时才
        // 收集下标。割行/GMI 行在大 LP 上近乎全稠密，索引式内循环（间接
        // 寻址）阻碍 SIMD——实测 mknap-5x300 t1 do_pivot 占墙钟 ~70%。
        let nnz = self.matrix[r].iter().filter(|&&v| v != 0.0).count();
        let dense_elim = nnz * 4 >= ncols;
        if !dense_elim {
            self.pivot_nz.clear();
            for j in 0..ncols {
                if self.matrix[r][j] != 0.0 {
                    self.pivot_nz.push(j);
                }
            }
        }
        // 行 r 缩放（稠密：连续除；稀疏：只遍历非零列）。
        if dense_elim {
            for v in self.matrix[r].iter_mut() {
                *v /= p;
            }
        } else {
            let row = &mut self.matrix[r];
            for &j in self.pivot_nz.iter() {
                row[j] /= p;
            }
        }
        self.matrix[r][q] = 1.0;
        for v in self.binv[r].iter_mut() {
            *v /= p;
        }
        // 出基列与入基旧值。
        let leaving = self.basis[r];
        let entering_old = if delta > 0.0 {
            self.col_bounds(q).0
        } else {
            self.col_bounds(q).1
        };
        // 消元系数 f（用未更新的 matrix[k][q]）：覆盖写暂存，无分配。
        self.pivot_f.clear();
        self.pivot_f.resize(m, 0.0);
        for k in 0..m {
            if k != r {
                self.pivot_f[k] = self.matrix[k][q];
            }
        }
        // 消元（k ≠ r；f_k 用未更新的 matrix[k][q]）。
        for k in 0..m {
            if k == r {
                continue;
            }
            let fk = self.pivot_f[k];
            if fk != 0.0 {
                // matrix 行 k（可变）与行 r（只读）：k≠r，经 split_at_mut
                // 取得不相交切片。
                let (row_k, row_r): (&mut Vec<f64>, &Vec<f64>) = if r < k {
                    let (left, right) = self.matrix.split_at_mut(k);
                    (&mut right[0], &left[r])
                } else {
                    let (left, right) = self.matrix.split_at_mut(k + 1);
                    (&mut left[k], &right[r - k - 1])
                };
                if dense_elim {
                    for (x, y) in row_k.iter_mut().zip(row_r.iter()) {
                        *x -= fk * *y;
                    }
                } else {
                    for &j in self.pivot_nz.iter() {
                        row_k[j] -= fk * row_r[j];
                    }
                }
                let (bk, br): (&mut Vec<f64>, &Vec<f64>) = if r < k {
                    let (left, right) = self.binv.split_at_mut(k);
                    (&mut right[0], &left[r])
                } else {
                    let (left, right) = self.binv.split_at_mut(k + 1);
                    (&mut left[k], &right[r - k - 1])
                };
                for (x, y) in bk.iter_mut().zip(br.iter()) {
                    *x -= fk * *y;
                }
                self.rhs[k] -= delta * t * fk;
            }
            self.matrix[k][q] = 0.0; // 数值卫生（保留在非零列循环之外）
        }
        self.rhs[r] = entering_old + delta * t;
        self.basis[r] = q;
        self.in_basis[leaving] = false;
        self.in_basis[q] = true;
        self.at_upper[leaving] = hit_upper;
    }

    /// 执行入基变量自身界翻转（基不变）：非基值移动 Δ = ±(ub−lb) 并按
    /// B⁻¹A_q 传播全部基变量值。
    fn do_flip(&mut self, q: usize, delta: f64, t: f64) {
        let (lb, ub) = self.col_bounds(q);
        let old_v = if delta > 0.0 { lb } else { ub };
        let new_v = if delta > 0.0 { ub } else { lb };
        debug_assert!((new_v - old_v - delta * t).abs() <= 1e-9 * (1.0 + (new_v - old_v).abs()));
        let d = new_v - old_v;
        let m = self.basis.len();
        for k in 0..m {
            let a = self.matrix[k][q];
            if a != 0.0 {
                self.rhs[k] -= d * a;
            }
        }
        self.at_upper[q] = delta > 0.0;
    }

    // ===================== 内部：阶段循环 =====================

    /// 单阶段主循环。返回 Optimal / Unbounded（仅 Phase-2 合法）/ LimitReached
    /// / Failed（数值异常）。
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
            // 初始 false 仅在 Flip/Leave 两臂都被覆盖（其余臂直接 return），
            // 尾部释放检查读取的是本步结果——初始值语法上必需但语义上不达。
            #[allow(unused_assignments)]
            let mut t_nondegenerate = false;
            let (q, delta) = match self.price(phase) {
                Pricing::Failed(msg) => return PhaseOutcome::Failed(msg),
                Pricing::Optimal => return PhaseOutcome::Optimal,
                Pricing::Candidate(q, delta) => (q, delta),
            };
            match self.ratio_test(q, delta) {
                Step::FilteredBlock(msg) => return PhaseOutcome::Failed(msg),
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
                    let p = self.matrix[row][q];
                    if p.abs() < self.tol.pivot {
                        return PhaseOutcome::Failed(format!(
                            "主元绝对值 {p} 低于容差 {}",
                            self.tol.pivot
                        ));
                    }
                    t_nondegenerate = t > self.tol.feasibility;
                    self.do_pivot(q, row, delta, t, hit_upper);
                    self.count_step(t);
                    self.pivots_since_refresh = self.pivots_since_refresh.saturating_add(1);
                    if self.pivots_since_refresh >= REFRESH_INTERVAL {
                        self.pivots_since_refresh = 0;
                        if let Err(msg) = self.refresh_values() {
                            return PhaseOutcome::Failed(msg);
                        }
                    }
                }
            }
            if self.rhs.iter().any(|v| !v.is_finite()) {
                return PhaseOutcome::Failed("基变量值出现 NaN/Inf".to_string());
            }
            if self.iter_call > bland_total {
                // 总迭代触发：永久 Bland（保底终止路径，不释放）。
                self.bland = true;
                self.bland_permanent = true;
            } else if self.degenerate_streak > DEGENERATE_STREAK_LIMIT {
                // 停滞触发：可释放 Bland（见 BLAND_DWELL 文档——setcover 类
                // 宽表上纯 Bland 爬行实测打满迭代上限）。
                self.bland = true;
            }
            self.bland_release_step(t_nondegenerate);
        }
    }

    /// Bland 释放记账（R8，run_phase 每步尾部调用）：`t_nondegenerate` =
    /// 本步非退化（目标真改善）。驻留满 [`BLAND_DWELL`] 且非退化 ⇒ 回
    /// Dantzig 定价（重入预算内）；永久 Bland（总迭代触发）与预算耗尽后
    /// 不释放——Bland 无循环定理的终止路径恒在场。
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

    /// 记账一次迭代：步长 ≤ feasibility 视为退化（无有效改善）。
    fn count_step(&mut self, t: f64) {
        self.iter_call = self.iter_call.saturating_add(1);
        self.iterations_total = self.iterations_total.saturating_add(1);
        if t <= self.tol.feasibility {
            self.degenerate_streak = self.degenerate_streak.saturating_add(1);
        } else {
            self.degenerate_streak = 0;
        }
    }

    /// 用原始数据重算 B⁻¹b（x_B = B⁻¹(b̂ − Σ_非基 Â_j·x_j)），抑制漂移。
    fn refresh_values(&mut self) -> std::result::Result<(), String> {
        let m = self.basis.len();
        let n = self.num_structural;
        // 非基结构列取值（非基松弛/人工列恒贴 0，无贡献）。
        let mut val = vec![0.0; n];
        for j in 0..n {
            if !self.in_basis[j] {
                val[j] = self.nonbasic_value(j);
            }
        }
        let mut d = self.row_rhs.clone();
        for k in 0..m {
            let mut acc = 0.0;
            for &(j, a) in &self.row_terms[k] {
                if val[j] != 0.0 {
                    acc += a * val[j];
                }
            }
            d[k] -= acc;
        }
        for i in 0..m {
            let mut acc = 0.0;
            for k in 0..m {
                let b = self.binv[i][k];
                if b != 0.0 {
                    acc += b * d[k];
                }
            }
            if !acc.is_finite() {
                return Err("B⁻¹b 重算出现 NaN/Inf".to_string());
            }
            self.rhs[i] = acc;
        }
        Ok(())
    }

    /// 两阶段编排：必要时 Phase-1 → 清理人工变量 → Phase-2。
    ///
    /// Phase-1 欠账记账：凡未走到 `cleanup_after_phase1` 的 Phase-1 出口
    /// （不可行/迭代上限/数值失败/无界）都置 `phase1_pending`——人工变量
    /// 未清理留基，该状态不可供后续 resolve 直接跑 Phase-2；只有清理完成
    /// 才复位（残基人工变量值恒 0 且行系数近零，Phase-2 不会再移动它）。
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
                PhaseOutcome::Unbounded => {
                    self.numerical_error = Some("Phase-1 出现无界方向（数学上不可能）".to_string());
                    self.phase1_pending = true;
                    return LpStatus::IterationLimit;
                }
                PhaseOutcome::Optimal => {
                    if self.phase1_infeasibility() > self.tol.feasibility {
                        // 人工变量带正值留基且未清理：状态对 Phase-2 无效，
                        // 但 snap 传播可能把基值推到容差内，必须记账强制
                        // 下一次 resolve 重跑 Phase-1。
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
            PhaseOutcome::Failed(msg) => {
                self.numerical_error = Some(msg);
                LpStatus::IterationLimit
            }
        }
    }

    /// Phase-1 最优时的人工变量总量（残基人工变量值之和，负值按 0 计）。
    fn phase1_infeasibility(&self) -> f64 {
        let ncs = self.ncols_static();
        let mut sum = 0.0;
        for i in 0..self.basis.len() {
            if self.basis[i] >= ncs {
                sum += self.rhs[i].max(0.0);
            }
        }
        sum
    }

    /// Phase-1 后处理：残基人工变量值清零并尽力驱赶出基（在非人工非基列上
    /// 以 t=0 主元换出）；整行无非人工候选的视为冗余行，人工变量留基且值
    /// 恒 0（其行系数对全部可入基列近零，后续相位不会再移动它）。
    fn cleanup_after_phase1(&mut self) {
        let ncs = self.ncols_static();
        for i in 0..self.basis.len() {
            if self.basis[i] < ncs {
                continue;
            }
            // 值在容差内（infeasibility 已判 0），强制精确置 0。
            self.rhs[i] = 0.0;
            let mut pick = usize::MAX;
            for j in 0..ncs {
                if self.in_basis[j] {
                    continue;
                }
                if self.matrix[i][j].abs() >= self.tol.pivot {
                    pick = j;
                    break; // 最小下标
                }
            }
            if pick != usize::MAX {
                let delta = if self.at_upper[pick] { -1.0 } else { 1.0 };
                self.do_pivot(pick, i, delta, 0.0, false);
                self.count_step(0.0);
            }
        }
    }

    // ===================== 内部：热启动 =====================

    /// 非基列贴界双向检查（相对上次有效界，即 `self.prev_bounds`——由
    /// resolve 在调用前 swap 换入）：
    /// * 收紧越界：贴上界且新 ub < 当前值 → snap 到新 ub；贴下界且新 lb >
    ///   当前值 → snap 到新 lb。
    /// * 松弛失贴：当前值严格位于新 (lb,ub) 内部 → snap 到较近界；等距并列
    ///   按 Phase-2 检验数符号选侧（rc<0 选下界、rc≥0 选上界）。
    /// 每个 snap 的 Δx_j 按 B⁻¹A_j 传播更新全部基变量值。
    fn snap_nonbasics(&mut self) {
        let n = self.num_structural;
        let feas = self.tol.feasibility;
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
            // 是否严格位于新界内部（+∞ 上界视为距离无穷远）。
            let inside =
                (old_v - nlb).abs() > feas && !(nub.is_finite() && (nub - old_v).abs() <= feas);
            let (target, to_upper) = if self.at_upper[j] {
                if nub.is_finite() && nub < old_v - feas {
                    (nub, true) // 收紧越界：snap 到新 ub
                } else if inside {
                    self.nearest_or_rc(j, old_v, nlb, nub)
                } else if nub.is_finite() {
                    (nub, true) // 仍贴上界
                } else {
                    (nlb, false) // 新上界 +∞：只能贴下界
                }
            } else {
                if nlb > old_v + feas {
                    (nlb, false) // 收紧越界：snap 到新 lb
                } else if inside {
                    self.nearest_or_rc(j, old_v, nlb, nub)
                } else {
                    (nlb, false) // 仍贴下界
                }
            };
            let delta = target - old_v;
            if delta != 0.0 {
                let m = self.basis.len();
                for k in 0..m {
                    let a = self.matrix[k][j];
                    if a != 0.0 {
                        self.rhs[k] -= delta * a;
                    }
                }
            }
            self.at_upper[j] = to_upper;
        }
    }

    /// 松弛失贴时的 snap 目标：较近的界；等距（带宽内）按检验数符号选侧。
    fn nearest_or_rc(&self, j: usize, v: f64, lb: f64, ub: f64) -> (f64, bool) {
        let d_lb = (v - lb).abs();
        let d_ub = if ub.is_finite() {
            (ub - v).abs()
        } else {
            f64::INFINITY
        };
        let tie = self.tol.feasibility;
        if d_lb < d_ub - tie {
            (lb, false)
        } else if d_ub < d_lb - tie {
            (ub, true)
        } else if self.rc_of(j) < 0.0 {
            (lb, false)
        } else {
            (ub, true)
        }
    }

    /// 基变量越界检查：任何基变量超出其界（含残基人工变量值非 0）或出现
    /// 非有限值 ⇒ 需要 Phase-1 修复。
    fn basics_violated(&self) -> bool {
        let ncs = self.ncols_static();
        let feas = self.tol.feasibility;
        for i in 0..self.basis.len() {
            let c = self.basis[i];
            if !self.rhs[i].is_finite() {
                return true;
            }
            if c >= ncs {
                if self.rhs[i].abs() > feas {
                    return true; // 残基人工变量必须为 0
                }
            } else {
                let (lb, ub) = self.col_bounds(c);
                if self.rhs[i] < lb - feas || self.rhs[i] > ub + feas {
                    return true;
                }
            }
        }
        false
    }

    // ===================== 内部：解构造 =====================

    /// 校验 bounds 实参：lb>ub ⇒ Infeasible；长度不符/非有限界 ⇒ 数值失败
    ///（映射 IterationLimit，原因记入 numerical_error）。
    // ===================== 内部：界变化后的对偶修复快路径 =====================

    /// 基中是否存在人工列（对偶快路径前置排除：人工列在基 ⇒ 对偶可行性
    /// 语义不成立，走重建）。
    fn artificials_in_basis(&self) -> bool {
        let ncs = self.ncols_static();
        self.basis.iter().any(|&c| c >= ncs)
    }

    /// 静态非基列对偶可行性预检（min 内部口径：贴下界 rc ≥ −tol、贴上界
    /// rc ≤ +tol；固定列不参与）。界收紧后 snap 只挪非基取值、基与成本
    /// 未动 ⇒ 理论上恒可行；本预检是防御（snap 跨界翻转等罕见形态会使
    /// rc 符号与所贴界失配），失败即放弃快路径。
    fn dual_feasible(&mut self) -> bool {
        let ncs = self.ncols_static();
        // 新鲜度规范 (b)：入口一次批量重算（resolve 入口距上次 price 可能
        // 隔 phase 切换/add_row/整基重建，不得复用 price 的缓存；原逐列
        // rc_of 改为入口一次批量，&self→&mut self，调用点均在 &mut 上下文）。
        self.compute_rc_batch(2);
        for j in 0..ncs {
            if self.in_basis[j] {
                continue;
            }
            let (lb, ub) = self.col_bounds(j);
            if ub - lb <= self.tol.feasibility {
                continue; // 固定列：无对偶语义
            }
            let rc = self.rc_scratch[j];
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

    /// 有界变量对偶单纯形修复：从当前（对偶可行、原始违约）基出发，每次
    /// 选违反最重的基行出基（落到其被违反的界）、在对偶比值测试选出的
    /// 入基列上做标准表格主元，直至全部基变量回到界内。
    ///
    /// 与 [`do_pivot`] 的口径衔接：出基变量落界 ⇒ `hit_upper` = 落上界？
    /// 步长 t 由「出基变量值恰好到达被违反界」解出；入基列从其所贴界沿
    /// delta 方向移动 t。迭代计数与 Phase 路径同口径累加（统计诚实）。
    /// 任何异常（无可入列、步长非法、超预算、检验数非有限）⇒ [`DualRepair::GiveUp`]
    /// 由调用方回退整基重建——包括「对偶无界 ⇒ 原始不可行」的情形也保守
    /// 回退（不可行判定仍由 Phase-1 给出，快路径绝不抢答）。
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
            // 新鲜度规范 (c)：外层循环每轮迭代顶部重算检验数——上一轮
            // do_pivot 已改表，下一轮列扫描读取的 rc 必须反映当前表格
            // （扫描“只对过滤列读取”是省工而非新鲜度保证，二者不得混淆）。
            self.compute_rc_batch(2);
            // (1) 选相对违反最重的基行（below = 低于 lb ⇒ 需升到 lb）。
            //     过滤阈值用**绝对** feasibility——与入口 [`Self::basics_violated`]
            //     （绝对 feas）同口径、严于终检的 10·feas：若用相对阈值
            //     feas·(1+|rhs|)，大数值量级（|rhs|≈1e6 ⇒ 阈值≈1e-3）下
            //     小越界被宣布已修复，终检随即降级 IterationLimit（可修复
            //     状态变伪数值失败）。行间**排序**仍按界幅归一（只定优先级，
            //     不定通过与否）。
            let mut worst: Option<(usize, f64, bool)> = None;
            for i in 0..m {
                let (lb, ub) = self.col_bounds(self.basis[i]);
                let v_lb = lb - self.rhs[i];
                let v_ub = self.rhs[i] - ub;
                let (viol, below) = if v_lb > v_ub { (v_lb, true) } else { (v_ub, false) };
                // NaN 不跳过（交由选择与步长校验走 GiveUp 回退，与既有
                // 数值异常口径一致）。
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
            // (2) 对偶比值测试：方向可行（below ⇒ δ·a<0，above ⇒ δ·a>0）的
            //     非基列中选 ratio = 有符号 rc / |a| 最小者；正常期并列取
            //     |主元| 大者（数值稳），退化期 Bland 最小列下标防循环。
            let mut best: Option<(f64, f64, usize, f64)> = None; // (ratio, |a|, q, delta)
            for j in 0..ncs {
                if self.in_basis[j] {
                    continue;
                }
                let (lb, ub) = self.col_bounds(j);
                if ub - lb <= self.tol.feasibility {
                    continue;
                }
                let a = self.matrix[r][j];
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
                let rc = self.rc_scratch[j];
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
                            ratio < br || (ratio == br && (a.abs() > ba || (a.abs() == ba && j < bq)))
                        }
                    }
                };
                if better {
                    best = Some((ratio, a.abs(), j, delta));
                }
            }
            let Some((_, _, q, delta)) = best else {
                // 无方向可行列（对偶无界 ⇒ 原始不可行）：保守回退 Phase-1 判定。
                return DualRepair::GiveUp;
            };
            // (3) 步长 t：出基变量恰好落到被违反的界。
            let a = self.matrix[r][q];
            let (lb_r, ub_r) = self.col_bounds(self.basis[r]);
            let (target, hit_upper) = if below { (lb_r, false) } else { (ub_r, true) };
            let t = (target - self.rhs[r]) / (-delta * a);
            if !t.is_finite() || t < 0.0 {
                return DualRepair::GiveUp;
            }
            self.do_pivot(q, r, delta, t, hit_upper);
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
                return Some(LpStatus::Infeasible); // 契约：lb>ub 直接 Infeasible
            }
        }
        None
    }

    /// 重置本次调用计数器。
    fn start_call(&mut self) {
        self.iter_call = 0;
        self.degenerate_streak = 0;
        self.bland = false;
        self.bland_permanent = false;
        self.bland_dwell = 0;
        self.bland_rearms = 0;
        self.pivots_since_refresh = 0;
    }

    /// 设置单次调用迭代上限提示（R8，割循环专用；None 恢复默认公式）。
    /// 只影响后续 solve/resolve 调用的上限判定，不影响其他语义。
    pub fn set_iter_limit_hint(&mut self, hint: Option<u64>) {
        self.iter_limit_hint = hint;
    }

    /// 构造 LpSolution：提取 x（基本值在容差内夹回界）、按 x 重算目标
    /// （值与目标自洽）、Optimal 时做最终原始可行性/人工变量校验并计算对偶。
    fn build_solution(&mut self, mut status: LpStatus) -> LpSolution {
        // 轴1 回归门：每次 solve/resolve 收口做一次索引全量比对（与每 64
        // 主元的周期断言互补——短调用（< 64 主元）也能捕获簿记漂移；禁用态
        // 索引全空，天然空转）。
        #[cfg(debug_assertions)]
        self.nz_debug_verify();
        let n = self.num_structural;
        let m = self.basis.len();
        let ncs = self.ncols_static();
        let feas = self.tol.feasibility;
        // 提取结构变量值。
        let mut x = vec![0.0; n];
        for j in 0..n {
            if self.in_basis[j] {
                let mut v = f64::NAN;
                for i in 0..m {
                    if self.basis[i] == j {
                        v = self.rhs[i];
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
                if !self.rhs[i].is_finite() {
                    ok = false;
                    break;
                }
                if c >= ncs {
                    if self.rhs[i].abs() > 10.0 * feas {
                        ok = false;
                        break;
                    }
                } else {
                    let (lb, ub) = self.col_bounds(c);
                    if self.rhs[i] < lb - 10.0 * feas || self.rhs[i] > ub + 10.0 * feas {
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
            // 目标：按报告 x 重算（内部 min 系数 → 原始 sense + offset）。
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
            // 对偶：y = c_B·B⁻¹（内部 min 口径），× row_scale 还原到原始 rhs。
            let mut y = vec![0.0; m];
            for i in 0..m {
                let cb = self.cost_of(self.basis[i], 2);
                if cb != 0.0 {
                    for k in 0..m {
                        y[k] += cb * self.binv[i][k];
                    }
                }
            }
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

/// 单阶段主循环结果。
enum PhaseOutcome {
    Optimal,
    Unbounded,
    LimitReached,
    Failed(String),
}

// ===========================================================================
// 单元测试
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::MilpError;
    use crate::model::{Constraint, LinearExpr, VarType, Variable};

    /// 构造默认容差（不可调用 Tolerances::default()——WP1 桩未实现）。
    fn tol() -> Tolerances {
        Tolerances {
            feasibility: 1e-9,
            optimality: 1e-9,
            pivot: 1e-9,
            integrality: 1e-6,
            zero: 1e-12,
        }
    }

    fn var(name: &str, lb: f64, ub: f64) -> Variable {
        Variable {
            name: name.to_string(),
            var_type: VarType::Continuous,
            lb,
            ub,
        }
    }

    /// 直接以结构体字面量构造 Problem（绕过 WP1 尚未实现的构造器方法）。
    /// terms 必须已按变量下标升序、无重复、无零系数（LinearExpr 不变量）。
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
                .map(|(terms, op, rhs)| Constraint {
                    name: None,
                    expr: LinearExpr { terms },
                    op,
                    rhs,
                })
                .collect(),
        }
    }

    /// （M6）add_row 的割松弛列不得与人工列段撞号：1 Eq + 1 Le 问题根解后
    /// add_row，表格必须仍是 B⁻¹A——每个基列在表格中恰为单位向量（行 i
    /// 处 1、他处 0）、基列互异、割松弛列落在静态段内。旧实现把新割松弛列
    /// 分配在 `n + slack.len()`（= 人工列段首），与 Eq 行构建期人工列共址，
    /// 该列在旧 Eq 行残留系数 ⇒ 非单位向量，表格 ≠ B⁻¹A（靠后续整基重建
    /// 偶然兜住）。
    #[test]
    fn test_add_row_slack_column_does_not_alias_artificial_segment() {
        // min 2x + y s.t. x + y = 3（Eq，构建期产人工列）、x − y ≤ 1（Le），
        // x, y ∈ [0, 5] 连续。根解 (0, 3)。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 5.0), var("y", 0.0, 5.0)],
            vec![(0, 2.0), (1, 1.0)],
            vec![
                (vec![(0, 1.0), (1, 1.0)], RelOp::Eq, 3.0),
                (vec![(0, 1.0), (1, -1.0)], RelOp::Le, 1.0),
            ],
        );
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let b = Bounds::from_problem(&p);
        let sol = s.solve(&b);
        assert_eq!(sol.status, LpStatus::Optimal);
        // 前置：Eq 行 ⇒ 构建期存在人工列（Phase-1 清理后列仍占用列空间），
        // 这是撞号成立的前提。
        assert!(s.num_art >= 1, "Eq 行问题构建期必有人工列");
        let ncs_before = s.ncols_static();

        // 根解后加割行 x ≤ 2.5（Le，当前解违反 ⇒ 走人工列路径之外的常规
        // 松弛基路径；再补一条 Ge 割走人工列路径，双向覆盖）。
        let k1 = s.add_row(&[(0, 1.0)], RelOp::Le, 2.5).unwrap();
        let k2 = s.add_row(&[(1, 1.0)], RelOp::Ge, 1.5).unwrap();
        assert_eq!((k1, k2), (0, 1), "割序号 0 起连续递增");

        // 基列互异（撞号会使两行的基列共址）。
        let mut sorted_basis = s.basis.clone();
        sorted_basis.sort_unstable();
        for w in sorted_basis.windows(2) {
            assert_ne!(w[0], w[1], "基列重复（人工列共址）：basis = {:?}", s.basis);
        }
        // 表格 = B⁻¹A 的必要条件：每个基列是单位向量。
        for (i, &c) in s.basis.iter().enumerate() {
            for (r2, row) in s.matrix.iter().enumerate() {
                let expect = if r2 == i { 1.0 } else { 0.0 };
                assert!(
                    (row[c] - expect).abs() <= 1e-9,
                    "基列 {c} 非单位向量：matrix[{r2}][{c}] = {}（basis = {:?}）",
                    row[c],
                    s.basis
                );
            }
        }
        // 割松弛列落在静态段内（避开人工列段）；静态段因新割恰好扩 2。
        assert_eq!(s.ncols_static(), ncs_before + 2);
        for c in s.row_slack.iter().flatten() {
            assert!(*c < s.ncols_static(), "松弛列 {c} 侵入人工列段");
        }
        // in_basis 与 basis 一致（长度 = 全列空间，基位数 = 行数）。
        assert_eq!(
            s.in_basis.iter().filter(|&&v| v).count(),
            s.basis.len(),
            "in_basis 计数须等于行数"
        );
        for &c in &s.basis {
            assert!(s.in_basis[c], "基列 {c} 的 in_basis 未置位");
        }

        // 端到端：加割后 resolve 收敛到受割约束的最优解（目标 min 2x+y 在
        // x+y=3 上 = x+3 ⇒ 唯一最优 x=0、y=3，满足割 x ≤ 2.5、y ≥ 1.5）。
        let after = s.resolve(&b);
        assert_eq!(after.status, LpStatus::Optimal);
        assert_close(after.objective, 3.0);
        assert!((after.x[0] - 0.0).abs() < 1e-7);
        assert!((after.x[1] - 3.0).abs() < 1e-7);
    }

    /// （M7）dual_repair 行违反过滤用绝对 feasibility 口径：大数值量级
    ///（基变量值 ≈ 1.5e6）下，收紧基变量上界 1e-4 产生的小越界必须被对偶
    /// 修复为 Optimal。旧实现的相对阈值 feas·(1+|rhs|) ≈ 1.5e-3 会把该
    /// 1e-4 越界宣布为已修复，终检（10·feas = 1e-8）随即降级 IterationLimit
    ///（可修复状态变伪数值失败）。
    #[test]
    fn test_dual_repair_absolute_threshold_large_magnitude() {
        // max x s.t. x + y ≤ 1.5e6，x ∈ [0, 2e6]、y ∈ [0, 1]：根解 x 基
        //（值 1.5e6，slack 非基 0）；收紧 x.ub 至 1.5e6 − 1e-4 ⇒ 基变量
        // 越上界 1e-4（snap 不触碰基变量）。
        let p = prob(
            Sense::Maximize,
            vec![var("x", 0.0, 2.0e6), var("y", 0.0, 1.0)],
            vec![(0, 1.0)],
            vec![(vec![(0, 1.0), (1, 1.0)], RelOp::Le, 1.5e6)],
        );
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let mut b = Bounds::from_problem(&p);
        let sol = s.solve(&b);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 1.5e6);

        b.ub[0] = 1.5e6 - 1.0e-4;
        let warm = s.resolve(&b);
        assert_eq!(
            warm.status,
            LpStatus::Optimal,
            "1e-4 越界在 |rhs|≈1.5e6 下须被对偶修复（非伪数值失败）"
        );
        assert_close(warm.objective, 1.5e6 - 1.0e-4);
        assert!((warm.x[0] - (1.5e6 - 1.0e-4)).abs() < 1e-6);
        assert!(warm.x[1] >= 0.0 && warm.x[1] <= 1.0);
        // 冷解对照（2026-10 L5 修复后恢复）：旧并列带宽 feas·(1+|t_min|) 在
        // t_min≈1.5e6 时 ≈1.5e-3，把比出基候选小 1e-4 的 x 自身界翻转当作
        // 假并列、按 |主元| 破键输给 Leave（步长越界 1e-4），冷路径 x 越上界
        // 与热解分叉，故当时只钉热路径；带宽与绝对上限取小（封顶
        // feas·1e3=1e-6）后冷路径选真最小候选，热冷一致（y 无目标成本，
        // 最优面在 y 上退化——不钉 y 顶点，只钉终态/目标/x 贴界/行可行）。
        let mut cold = Simplex::new(&p, t).unwrap();
        let cold_sol = cold.solve(&b);
        assert_eq!(
            cold_sol.status,
            LpStatus::Optimal,
            "冷解不得再踩并列带宽扰动"
        );
        assert_close(cold_sol.objective, 1.5e6 - 1.0e-4);
        assert!(
            (cold_sol.x[0] - (1.5e6 - 1.0e-4)).abs() < 1e-6,
            "冷解 x = {}（须贴收紧后上界，不越界）",
            cold_sol.x[0]
        );
        assert!(cold_sol.x[1] >= 0.0 && cold_sol.x[1] <= 1.0);
        assert!(
            cold_sol.x[0] + cold_sol.x[1] <= 1.5e6 + 1e-6,
            "冷解行可行性：x+y = {}",
            cold_sol.x[0] + cold_sol.x[1]
        );
    }

    /// （L5）ratio_test 并列带宽大数量级防护：t_min≈1.5e6 时旧带宽
    /// feas·(1+|t_min|)≈1.5e-3 会把真正更小的界翻转候选（差 1e-4）吞进
    /// 并列集合、按 |主元| 破键输给 Leave（步长越过真最小比值 1e-4）；
    /// 取小后带宽 = feas·RATIO_BAND_CAP = 1e-6 < 1e-4 ⇒ 只有真最小候选
    /// 在集合内 ⇒ Flip 胜出。同时钉 |t_min| ≤ CAP−1 的常规量级行为
    /// 逐值不变（并列仍按 |主元| 破键）。
    #[test]
    fn test_ratio_test_band_large_magnitude() {
        // —— 大数量级：手钉单行表格状态（基 = 行松弛列，rhs = 1.5e6）——
        // max x s.t. x+y ≤ 1.5e6；x ∈ [0, 1.5e6−1e-4]（界翻转距离
        // 1499999.9999 = 真最小 t），y ∈ [0,1]。入基 q = x（delta=+1）。
        let p = prob(
            Sense::Maximize,
            vec![var("x", 0.0, 1.5e6 - 1.0e-4), var("y", 0.0, 1.0)],
            vec![(0, 1.0)],
            vec![(vec![(0, 1.0), (1, 1.0)], RelOp::Le, 1.5e6)],
        );
        let mut s = Simplex::new(&p, tol()).unwrap();
        s.bounds = Bounds::from_problem(&p);
        let slack_col = s.ncols_static() - 1; // 唯一 Le 行的松弛列
        s.matrix = vec![vec![1.0, 1.0, 1.0]]; // B⁻¹A 恒等形态
        s.rhs = vec![1.5e6];
        s.basis = vec![slack_col];
        s.in_basis = vec![false, false, true];
        s.at_upper = vec![false; s.ncols_static()];
        match s.ratio_test(0, 1.0) {
            Step::Flip { t } => assert!(
                (t - 1499999.9999).abs() < 1e-6,
                "真最小比值应为 x 的界翻转（1499999.9999），得到 Flip t = {t}"
            ),
            Step::Leave { row, t, .. } => {
                panic!("t_min≈1.5e6 下带宽不得吞掉更小的界翻转候选：误选出基行 {row}（t = {t}）")
            }
            _ => panic!("预期 Flip（NoBlock/FilteredBlock 均非预期终态）"),
        }

        // —— 常规量级（t_min = 2.0 ≪ CAP−1 = 999）：带宽与旧公式逐值一致，
        // 并列仍按 |主元| 最大破键（回归钉：防护不改变既有取解次序）——
        // 行 0：t = 2.0、|主元| 0.5；行 1：t = 2.0000000005（差 5e-10，在
        // 带宽 feas·3 = 3e-9 内）、|主元| 1.0 ⇒ 旧行为选行 1，新行为相同。
        let p2 = prob(
            Sense::Maximize,
            vec![var("x", 0.0, 10.0), var("y", 0.0, 10.0)],
            vec![(0, 1.0)],
            vec![
                (vec![(0, 0.5)], RelOp::Le, 1.0),
                (vec![(0, 1.0)], RelOp::Le, 2.0000000005),
            ],
        );
        let mut s2 = Simplex::new(&p2, tol()).unwrap();
        s2.bounds = Bounds::from_problem(&p2);
        let slack0 = s2.ncols_static() - 2;
        let slack1 = s2.ncols_static() - 1;
        s2.matrix = vec![vec![0.5, 0.0, 1.0, 0.0], vec![1.0, 0.0, 0.0, 1.0]];
        s2.rhs = vec![1.0, 2.0000000005];
        s2.basis = vec![slack0, slack1];
        s2.in_basis = vec![false, false, true, true];
        s2.at_upper = vec![false; s2.ncols_static()];
        match s2.ratio_test(0, 1.0) {
            Step::Leave { row, t, .. } => {
                assert_eq!(row, 1, "并列（差 5e-10 ≤ 带宽 3e-9）须按 |主元| 最大选行 1");
                assert!((t - 2.0000000005).abs() < 1e-12, "步长 = {t}");
            }
            _ => panic!("预期 Leave（Flip/NoBlock/FilteredBlock 均非预期终态）"),
        }
    }

    fn assert_close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 1e-7 * (1.0 + b.abs()),
            "数值不符: {a} vs {b}"
        );
    }

    // ---- Bounds ----

    #[test]
    fn test_bounds_from_problem_and_tighten() {
        let p = prob(
            Sense::Minimize,
            vec![var("x", -5.0, 8.0), var("y", 0.0, f64::INFINITY)],
            vec![],
            vec![],
        );
        let mut b = Bounds::from_problem(&p);
        assert_eq!(b.lb, vec![-5.0, 0.0]);
        assert!(b.ub[1].is_infinite());
        // 收紧下界。
        assert!(b.tighten_lb(0, -2.0));
        assert_eq!(b.lb[0], -2.0);
        // 只收紧：传入更松的界不放宽。
        assert!(b.tighten_lb(0, -9.0));
        assert_eq!(b.lb[0], -2.0);
        // 新区间为空。
        assert!(!b.tighten_lb(0, 9.0));
        assert_eq!(b.lb[0], -2.0); // 失败不修改
        assert!(b.tighten_ub(0, 8.0));
        assert!(!b.tighten_ub(0, -3.0));
        // 恰好相等允许。
        assert!(b.tighten_lb(0, 8.0));
        assert_eq!((b.lb[0], b.ub[0]), (8.0, 8.0));
        // 越界下标不 panic。
        assert!(!b.tighten_lb(7, 0.0));
        assert!(!b.tighten_ub(7, 0.0));
    }

    // ---- 基础状态 ----

    #[test]
    fn test_small_lp_optimal_and_duals() {
        // min 2x+3y s.t. x+y>=4, x+2y>=6, [0,10]^2 → 最优 (2,2)，目标 10。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0), var("y", 0.0, 10.0)],
            vec![(0, 2.0), (1, 3.0)],
            vec![
                (vec![(0, 1.0), (1, 1.0)], RelOp::Ge, 4.0),
                (vec![(0, 1.0), (1, 2.0)], RelOp::Ge, 6.0),
            ],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 10.0);
        assert_close(sol.x[0], 2.0);
        assert_close(sol.x[1], 2.0);
        // 对偶 = ∂(内部 min 目标)/∂rhs（行 max|a|=1，scale=1）。
        assert_close(sol.duals[0], 1.0);
        assert_close(sol.duals[1], 1.0);
        assert!(s.iterations() > 0);
    }

    #[test]
    fn test_maximize_sense_and_offset() {
        // max 3x+2y s.t. x+y<=4, x+3y<=6, [0,3]^2 → 最优 (3,1)，目标 11。
        let mut p = prob(
            Sense::Maximize,
            vec![var("x", 0.0, 3.0), var("y", 0.0, 3.0)],
            vec![(0, 3.0), (1, 2.0)],
            vec![
                (vec![(0, 1.0), (1, 1.0)], RelOp::Le, 4.0),
                (vec![(0, 1.0), (1, 3.0)], RelOp::Le, 6.0),
            ],
        );
        p.obj_offset = 100.0;
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 111.0);
        assert_close(sol.x[0], 3.0);
        assert_close(sol.x[1], 1.0);
    }

    #[test]
    fn test_unbounded() {
        // min -x, x∈[0,+∞)，无约束 → 无界。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, f64::INFINITY)],
            vec![(0, -1.0)],
            vec![],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Unbounded);
        assert!(sol.objective.is_nan());
    }

    #[test]
    fn test_infeasible() {
        // x<=1 且 x>=2 → 不可行。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0)],
            vec![(0, 1.0)],
            vec![
                (vec![(0, 1.0)], RelOp::Le, 1.0),
                (vec![(0, 1.0)], RelOp::Ge, 2.0),
            ],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Infeasible);
        assert!(sol.objective.is_nan());
    }

    #[test]
    fn test_equality_infeasible() {
        // x+y=5, x,y∈[0,2] → Phase-1 最优 > 0 → 不可行。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 2.0), var("y", 0.0, 2.0)],
            vec![(0, 1.0), (1, 1.0)],
            vec![(vec![(0, 1.0), (1, 1.0)], RelOp::Eq, 5.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        assert_eq!(
            s.solve(&Bounds::from_problem(&p)).status,
            LpStatus::Infeasible
        );
    }

    // ---- 人工变量路径 ----

    #[test]
    fn test_ge_positive_rhs_artificial() {
        // min x s.t. x>=3（x 停在 lb=0，剩余量恒负 → σ=+1 人工路径）。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0)],
            vec![(0, 1.0)],
            vec![(vec![(0, 1.0)], RelOp::Ge, 3.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 3.0);
        assert_close(sol.x[0], 3.0);
    }

    #[test]
    fn test_le_negative_rhs_artificial() {
        // min x s.t. -x<=-3（Le 行负右端，r=-3 → σ=-1 人工路径，行不翻转）。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0)],
            vec![(0, 1.0)],
            vec![(vec![(0, -1.0)], RelOp::Le, -3.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 3.0);
        assert_close(sol.x[0], 3.0);
    }

    #[test]
    fn test_equality_row() {
        // 2x+y=10, min x+y, [0,10]^2 → 最优 x=5,y=0，目标 5。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0), var("y", 0.0, 10.0)],
            vec![(0, 1.0), (1, 1.0)],
            vec![(vec![(0, 2.0), (1, 1.0)], RelOp::Eq, 10.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 5.0);
        assert_close(sol.x[0], 5.0);
        assert_close(sol.x[1], 0.0);
    }

    #[test]
    fn test_ge_residual_zero_degenerate() {
        // x∈[4,10], min y s.t. x+y>=4：x 停在 lb=4，残差 0 → 剩余变量退化入基。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 4.0, 10.0), var("y", 0.0, 10.0)],
            vec![(1, 1.0)],
            vec![(vec![(0, 1.0), (1, 1.0)], RelOp::Ge, 4.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 0.0);
        assert_close(sol.x[0], 4.0);
        assert_close(sol.x[1], 0.0);
    }

    // ---- 退化与 Bland ----

    #[test]
    fn test_degenerate_assignment_3x3() {
        // 3×3 指派（全部等价成本 2）：高度退化，靠完整 Bland 规则保证终止。
        let mut vars = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                vars.push(var(&format!("x{i}{j}"), 0.0, 1.0));
            }
        }
        let mut obj = Vec::new();
        for k in 0..9 {
            obj.push((k, 2.0));
        }
        let mut cons = Vec::new();
        for i in 0..3 {
            cons.push((
                (0..3).map(|j| (i * 3 + j, 1.0)).collect::<Vec<_>>(),
                RelOp::Eq,
                1.0,
            ));
        }
        for j in 0..3 {
            cons.push((
                (0..3).map(|i| (i * 3 + j, 1.0)).collect::<Vec<_>>(),
                RelOp::Eq,
                1.0,
            ));
        }
        let p = prob(Sense::Minimize, vars, obj, cons);
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 6.0);
        let total: f64 = sol.x.iter().sum();
        assert_close(total, 3.0);
        for i in 0..3 {
            let row: f64 = (0..3).map(|j| sol.x[i * 3 + j]).sum();
            assert_close(row, 1.0);
        }
    }

    #[test]
    fn test_assignment_3x3_optimal_value() {
        // 非对称成本：最优指派 (0,1),(1,0),(2,2)，值 1+2+2=5。
        let cost = [[4.0, 1.0, 3.0], [2.0, 0.0, 5.0], [3.0, 2.0, 2.0]];
        let mut vars = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                vars.push(var(&format!("x{i}{j}"), 0.0, 1.0));
            }
        }
        let mut obj = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                obj.push((i * 3 + j, cost[i][j]));
            }
        }
        let mut cons = Vec::new();
        for i in 0..3 {
            cons.push((
                (0..3).map(|j| (i * 3 + j, 1.0)).collect::<Vec<_>>(),
                RelOp::Eq,
                1.0,
            ));
        }
        for j in 0..3 {
            cons.push((
                (0..3).map(|i| (i * 3 + j, 1.0)).collect::<Vec<_>>(),
                RelOp::Eq,
                1.0,
            ));
        }
        let p = prob(Sense::Minimize, vars, obj, cons);
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 5.0);
        assert_close(sol.x[0 * 3 + 1], 1.0);
        assert_close(sol.x[1 * 3 + 0], 1.0);
        assert_close(sol.x[2 * 3 + 2], 1.0);
    }

    // ---- 双侧定价合格性 ----

    #[test]
    fn test_upper_bound_decrease_entering() {
        // min 2x+y s.t. x+y>=4, x∈[0,3], y∈[0,10] → 最优 (0,4)，目标 4。
        // 求解路径：x 先翻到上界 3，随后必须以“贴上界减向入基”回到下界；
        // 只写 rc<-tol 单侧合格性会在 (3,1) 目标 7 处误报 Optimal。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 3.0), var("y", 0.0, 10.0)],
            vec![(0, 2.0), (1, 1.0)],
            vec![(vec![(0, 1.0), (1, 1.0)], RelOp::Ge, 4.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 4.0);
        assert_close(sol.x[0], 0.0);
        assert_close(sol.x[1], 4.0);
    }

    // ---- 热启动 ----

    #[test]
    fn test_warm_tighten_and_relax_ub_cycle() {
        // 契约反例：非基贴 ub=8 → 收紧至 0 snap 到 0 → 回溯恢复 ub=8（lb=-5）
        // 后值 0 成为严格内点，必须补 snap，否则 at_upper 语义失效漏改进。
        let p = prob(
            Sense::Minimize,
            vec![var("x", -5.0, 8.0)],
            vec![(0, -1.0)],
            vec![(vec![(0, 2.0)], RelOp::Le, 30.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let b0 = Bounds::from_problem(&p);
        let sol = s.solve(&b0);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -8.0);

        // 收紧 ub→0：非基贴上界 snap 到 0。
        let mut b1 = b0.clone();
        b1.ub[0] = 0.0;
        let sol = s.resolve(&b1);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 0.0);
        // 与冷启动一致。
        let mut s_cold = Simplex::new(&p, tol()).ok().unwrap();
        let cold = s_cold.solve(&b1);
        assert_eq!(cold.status, LpStatus::Optimal);
        assert_close(cold.objective, 0.0);

        // 回溯恢复 ub=8：值 0 严格内点 → snap 后继续优化到 -8（与冷启动一致）。
        let sol = s.resolve(&b0);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -8.0);
        let mut s_cold = Simplex::new(&p, tol()).ok().unwrap();
        let cold = s_cold.solve(&b0);
        assert_close(cold.objective, -8.0);
    }

    #[test]
    fn test_warm_tighten_lb_negative() {
        // 负下界非基列的收紧/回退一致性：min 2x+y s.t. x+y>=2, x∈[-5,5], y∈[0,10]。
        // 冷启动最优 x=-5（非基贴下界）, y=7，目标 -3。
        let p = prob(
            Sense::Minimize,
            vec![var("x", -5.0, 5.0), var("y", 0.0, 10.0)],
            vec![(0, 2.0), (1, 1.0)],
            vec![(vec![(0, 1.0), (1, 1.0)], RelOp::Ge, 2.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let b0 = Bounds::from_problem(&p);
        let sol = s.solve(&b0);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -3.0);
        assert_close(sol.x[0], -5.0);

        // 收紧 lb→-2：非基贴下界 snap 到 -2，目标 0。
        let mut b1 = b0.clone();
        b1.lb[0] = -2.0;
        let sol = s.resolve(&b1);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 0.0);
        let mut s_cold = Simplex::new(&p, tol()).ok().unwrap();
        let cold = s_cold.solve(&b1);
        assert_close(cold.objective, 0.0);

        // 回退 lb→-5：snap 回下界，目标恢复 -3（与冷启动一致）。
        let sol = s.resolve(&b0);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -3.0);
    }

    #[test]
    fn test_warm_branch_children_consistency() {
        // 分支定界典型流：min -(x+y) s.t. x+y<=3.5, x,y∈[0,3]，LP 最优 -3.5。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 3.0), var("y", 0.0, 3.0)],
            vec![(0, -1.0), (1, -1.0)],
            vec![(vec![(0, 1.0), (1, 1.0)], RelOp::Le, 3.5)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let root = Bounds::from_problem(&p);
        let sol = s.solve(&root);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -3.5);

        // 子节点 A：x<=2（非基贴上界 snap）。
        let mut ba = root.clone();
        ba.ub[0] = 2.0;
        let sol = s.resolve(&ba);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -3.5);
        assert_close(sol.x[0], 2.0);

        // 回溯恢复后子节点 B：x>=3（贴上界值 3 恰为新 lb，无需移动）。
        let mut bb = root.clone();
        bb.lb[0] = 3.0;
        let sol = s.resolve(&bb);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -3.5);
        assert_close(sol.x[0], 3.0);

        // 再回溯到根界：目标仍 -3.5（与冷启动一致）。
        let sol = s.resolve(&root);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -3.5);
        let mut s_cold = Simplex::new(&p, tol()).ok().unwrap();
        let cold = s_cold.solve(&root);
        assert_close(cold.objective, -3.5);
    }

    #[test]
    fn test_warm_repair_after_infeasible_phase1() {
        // 回归（修复判定漏洞）：上一调用以 Infeasible 收场（Phase-1 最优、
        // 人工变量带正值留基）后，snap 传播可把该人工变量基值恰好推到
        // ≤ tol.feasibility，仅凭 basics_violated() 会误判无需修复，Phase-2
        // 在等式行被零成本人工变量吸收的错误 LP 上求解（终检降级
        // IterationLimit，对可解问题谎报 LimitReached）。
        // min x  s.t. x + y = 5（Eq），y 恒固定为 2。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 2.0), var("y", 2.0, 2.0)],
            vec![(0, 1.0)],
            vec![(vec![(0, 1.0), (1, 1.0)], RelOp::Eq, 5.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let b1 = Bounds::from_problem(&p);
        // B1：x∈[0,2], y=2 ⇒ x+y ≤ 4 < 5 ⇒ Infeasible；Phase-1 把 x 翻到
        // 上界 2 后人工变量留基值 1。
        let sol = s.solve(&b1);
        assert_eq!(sol.status, LpStatus::Infeasible);

        // B2：x∈[0,3] 可行，唯一最优 x=3。非基 x 贴上界 2 → 新界内严格
        // 内点，snap 到较近上界 3（Δ=+1）恰把人工变量基值 1 推到 0。
        let mut b2 = b1.clone();
        b2.ub[0] = 3.0;
        let sol = s.resolve(&b2);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 3.0);
        assert_close(sol.x[0], 3.0);
        // 与冷启动一致。
        let mut s_cold = Simplex::new(&p, tol()).ok().unwrap();
        let cold = s_cold.solve(&b2);
        assert_eq!(cold.status, LpStatus::Optimal);
        assert_close(cold.objective, 3.0);
        assert_close(cold.x[0], 3.0);

        // 修复后的状态不再欠账：再次热启动（回到 B1）应正确报 Infeasible。
        let sol = s.resolve(&b1);
        assert_eq!(sol.status, LpStatus::Infeasible);
    }

    // ---- 近零行 ----

    #[test]
    fn test_near_zero_rows() {
        // 行 0（空 Le rhs=5）满足 → 剔除；行 1 正常。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0)],
            vec![(0, 1.0)],
            vec![(vec![], RelOp::Le, 5.0), (vec![(0, 1.0)], RelOp::Le, 3.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        assert_eq!(s.row_of_constraint(0), None);
        assert_eq!(s.row_of_constraint(1), Some(0));
        assert_eq!(s.num_rows(), 1);
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 0.0);
        assert_eq!(sol.duals.len(), 2);
        assert_close(sol.duals[0], 0.0); // 被剔除行恒 0
                                         // 空 Eq rhs=0 剔除、空 Ge rhs=0.5 违反 → Infeasible。
        let p2 = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 1.0)],
            vec![(0, 1.0)],
            vec![(vec![], RelOp::Eq, 0.0), (vec![], RelOp::Ge, 0.5)],
        );
        let mut s2 = Simplex::new(&p2, tol()).ok().unwrap();
        assert_eq!(
            s2.solve(&Bounds::from_problem(&p2)).status,
            LpStatus::Infeasible
        );
        // 空 Le rhs=-1 违反 → Infeasible。
        let p3 = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 1.0)],
            vec![(0, 1.0)],
            vec![(vec![], RelOp::Le, -1.0)],
        );
        let mut s3 = Simplex::new(&p3, tol()).ok().unwrap();
        assert_eq!(
            s3.solve(&Bounds::from_problem(&p3)).status,
            LpStatus::Infeasible
        );
    }

    // ---- 割行 add_row ----

    #[test]
    fn test_add_row_and_resolve() {
        // min -x s.t. x<=4, x∈[0,10] → 最优 x=4，目标 -4。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0)],
            vec![(0, -1.0)],
            vec![(vec![(0, 1.0)], RelOp::Le, 4.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let b = Bounds::from_problem(&p);
        let sol = s.solve(&b);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -4.0);

        // 追加被当前解违反的割 x<=2：resolve 修复后最优 -2。
        let k0 = s.add_row(&[(0, 1.0)], RelOp::Le, 2.0).ok().unwrap();
        assert_eq!(k0, 0);
        assert_eq!(s.row_of_cut(0), Some(1));
        assert_eq!(s.num_rows(), 2);
        assert_eq!(s.row_of_constraint(0), Some(0));
        let sol = s.resolve(&b);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -2.0);
        assert_close(sol.x[0], 2.0);

        // 追加满足的割 x<=10：最优不变。
        let k1 = s.add_row(&[(0, 1.0)], RelOp::Le, 10.0).ok().unwrap();
        assert_eq!(k1, 1);
        assert_eq!(s.row_of_cut(1), Some(2));
        let (row, val) = s.tableau_row(2);
        assert_eq!(row.len(), 1 + 3); // 结构 1 + 松弛 3
        assert!(val >= 0.0);
        let sol = s.resolve(&b);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -2.0);

        // 基/非基访问器一致性。
        assert_eq!(s.basis_columns().len(), s.num_rows());
        let nb = s.nonbasic_columns();
        assert_eq!(nb.len(), 4 - s.num_rows()); // 结构+松弛共 4 列
        assert_eq!(s.column_layout().num_structural, 1);
        assert_eq!(s.column_layout().slack.len(), 3);
    }

    #[test]
    fn test_add_row_ge_and_eq_defensive() {
        // 防御路径：Ge/Eq 语义的 add_row 也能正确修复。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0)],
            vec![(0, 1.0)],
            vec![],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let b = Bounds::from_problem(&p);
        let sol = s.solve(&b);
        assert_close(sol.objective, 0.0);
        // 割 x>=3（Ge 语义，违反当前 x=0）。
        s.add_row(&[(0, 1.0)], RelOp::Ge, 3.0).ok().unwrap();
        let sol = s.resolve(&b);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 3.0);
        // 割 x=5（Eq 语义）。
        s.add_row(&[(0, 1.0)], RelOp::Eq, 5.0).ok().unwrap();
        let sol = s.resolve(&b);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 5.0);
        assert_close(sol.x[0], 5.0);
    }

    // ---- 构建校验 ----

    #[test]
    fn test_new_rejects_invalid_inputs() {
        // lb = -∞。
        let p = prob(
            Sense::Minimize,
            vec![var("x", f64::NEG_INFINITY, 1.0)],
            vec![(0, 1.0)],
            vec![],
        );
        assert!(matches!(
            Simplex::new(&p, tol()),
            Err(MilpError::InvalidModel(_))
        ));
        // 目标系数 NaN。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 1.0)],
            vec![(0, f64::NAN)],
            vec![],
        );
        assert!(matches!(
            Simplex::new(&p, tol()),
            Err(MilpError::InvalidModel(_))
        ));
        // 右端 +∞。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 1.0)],
            vec![(0, 1.0)],
            vec![(vec![(0, 1.0)], RelOp::Le, f64::INFINITY)],
        );
        assert!(matches!(
            Simplex::new(&p, tol()),
            Err(MilpError::InvalidModel(_))
        ));
        // ub = -∞。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, f64::NEG_INFINITY)],
            vec![(0, 1.0)],
            vec![],
        );
        assert!(matches!(
            Simplex::new(&p, tol()),
            Err(MilpError::InvalidModel(_))
        ));
    }

    #[test]
    fn test_solve_lb_greater_than_ub() {
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 1.0)],
            vec![(0, 1.0)],
            vec![],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let bad = Bounds {
            lb: vec![2.0],
            ub: vec![1.0],
        };
        assert_eq!(s.solve(&bad).status, LpStatus::Infeasible);
        assert_eq!(s.resolve(&bad).status, LpStatus::Infeasible);
    }

    // ---- 其他边界 ----

    #[test]
    fn test_fixed_variable() {
        // 固定列 x∈[1,1] 不入基，y 补齐行 → 最优 y=2，目标 2。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 1.0, 1.0), var("y", 0.0, 10.0)],
            vec![(1, 1.0)],
            vec![(vec![(0, 1.0), (1, 1.0)], RelOp::Ge, 3.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 2.0);
        assert_close(sol.x[0], 1.0);
        assert_close(sol.x[1], 2.0);
    }

    #[test]
    fn test_duals_scaling_restoration() {
        // min -x s.t. 2x<=8（行除以 g=2 归一）：duals[0] = ∂z/∂rhs = -0.5。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0)],
            vec![(0, -1.0)],
            vec![(vec![(0, 2.0)], RelOp::Le, 8.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -4.0);
        assert_close(sol.x[0], 4.0);
        assert_close(sol.duals[0], -0.5);
    }

    #[test]
    fn test_solve_twice_and_resolve_stability() {
        // 重复 solve 与同界 resolve 结果稳定。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0), var("y", 0.0, 10.0)],
            vec![(0, 2.0), (1, 3.0)],
            vec![
                (vec![(0, 1.0), (1, 1.0)], RelOp::Ge, 4.0),
                (vec![(0, 1.0), (1, 2.0)], RelOp::Ge, 6.0),
            ],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let b = Bounds::from_problem(&p);
        let s1 = s.solve(&b);
        let s2 = s.solve(&b);
        assert_eq!((s1.status, s1.objective), (s2.status, s2.objective));
        let s3 = s.resolve(&b);
        assert_eq!(s3.status, LpStatus::Optimal);
        assert_close(s3.objective, 10.0);
    }

    #[test]
    fn test_resolve_before_solve_cold_path() {
        // 未 solve 直接 resolve：按冷启动处理。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0)],
            vec![(0, 1.0)],
            vec![(vec![(0, 1.0)], RelOp::Ge, 3.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.resolve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, 3.0);
    }

    #[test]
    fn test_beale_cycling_example() {
        // Beale 经典循环算例：Dantzig 规则会死循环，靠 Bland 模式兜底终止。
        // min -0.75x4 + 150x5 - 0.02x6 + 6x7
        // s.t. 0.25x4-60x5-0.04x6+9x7<=0; 0.5x4-90x5-0.02x6+3x7<=0; x6<=1。
        // 最优 x4=0.04, x6=1，目标 -0.05。
        let p = prob(
            Sense::Minimize,
            vec![
                var("x4", 0.0, f64::INFINITY),
                var("x5", 0.0, f64::INFINITY),
                var("x6", 0.0, f64::INFINITY),
                var("x7", 0.0, f64::INFINITY),
            ],
            vec![(0, -0.75), (1, 150.0), (2, -0.02), (3, 6.0)],
            vec![
                (
                    vec![(0, 0.25), (1, -60.0), (2, -0.04), (3, 9.0)],
                    RelOp::Le,
                    0.0,
                ),
                (
                    vec![(0, 0.5), (1, -90.0), (2, -0.02), (3, 3.0)],
                    RelOp::Le,
                    0.0,
                ),
                (vec![(2, 1.0)], RelOp::Le, 1.0),
            ],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let sol = s.solve(&Bounds::from_problem(&p));
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -0.05);
        assert_close(sol.x[0], 0.04);
        assert_close(sol.x[2], 1.0);
        assert_close(sol.x[1], 0.0);
        assert_close(sol.x[3], 0.0);
    }

    #[test]
    fn test_warm_recover_after_infeasible_node() {
        // B&B 关键路径：不可行子节点 resolve 后回溯恢复界，热启动必须仍能
        // 求得最优（不可行状态不得污染后续节点）。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0)],
            vec![(0, -1.0)],
            vec![(vec![(0, 1.0)], RelOp::Le, 4.0)],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        let root = Bounds::from_problem(&p);
        let sol = s.solve(&root);
        assert_close(sol.objective, -4.0);
        // 不可行节点：x>=5 与 x<=4 矛盾。
        let bad = Bounds {
            lb: vec![5.0],
            ub: vec![10.0],
        };
        assert_eq!(s.resolve(&bad).status, LpStatus::Infeasible);
        // 回溯到根界：恢复最优。
        let sol = s.resolve(&root);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_close(sol.objective, -4.0);
    }

    // ---- 热路径改造：批量检验数 / 主元零列跳过的值级等价 ----

    /// 逐列参考实现（与 rc_of 的运算顺序逐字一致；支持 phase 1/2，供
    /// phase-1 状态对照——rc_of 本体只算 phase 2）。
    fn rc_reference(s: &Simplex, j: usize, phase: u8) -> f64 {
        let m = s.basis.len();
        let mut rc = s.cost_of(j, phase);
        for i in 0..m {
            let cb = s.cost_of(s.basis[i], phase);
            if cb != 0.0 {
                rc -= cb * s.matrix[i][j];
            }
        }
        rc
    }

    /// compute_rc_batch 与逐列计算在全列空间逐列 == 相等（f64 == 断言：
    /// 同序运算 ⇒ 值级相等；±0 符号位差异对 == 不可分辨，禁止按 bit 断言）。
    fn assert_rc_batch_matches(s: &mut Simplex, phase: u8, ctx: &str) {
        let ncols = s.ncols_static() + s.num_art;
        s.compute_rc_batch(phase);
        for j in 0..ncols {
            let expect = if phase == 2 {
                s.rc_of(j)
            } else {
                rc_reference(s, j, 1)
            };
            assert!(
                s.rc_scratch[j] == expect,
                "{ctx}: 列 {j} phase {phase} 批量 {} != 逐列 {}",
                s.rc_scratch[j],
                expect
            );
        }
    }

    #[test]
    fn test_compute_rc_batch_matches_per_column_multi_state() {
        // min x+2y s.t. x+y=5（Eq→人工列）、x>=3（正右端 Ge→人工列）、
        // x+y<=9（Le→松弛列入基）。初始基含人工列（phase-1 真实起点，
        // c_b≠0 行存在），随后覆盖 solve / add_row（割行+新建人工列）/
        // resolve（热启动修复）/ do_pivot（手工改表）多状态。
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0), var("y", 0.0, 10.0)],
            vec![(0, 1.0), (1, 2.0)],
            vec![
                (vec![(0, 1.0), (1, 1.0)], RelOp::Eq, 5.0),
                (vec![(0, 1.0)], RelOp::Ge, 3.0),
                (vec![(0, 1.0), (1, 1.0)], RelOp::Le, 9.0),
            ],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        // 状态 1：全新初始表格（未 solve；基含人工列）。
        assert!(s.num_art > 0);
        assert!(s.basis.iter().any(|&c| c >= s.ncols_static()));
        assert_rc_batch_matches(&mut s, 1, "初始表格");
        assert_rc_batch_matches(&mut s, 2, "初始表格");
        // 状态 2：solve 后（phase-2 最优基；人工列仍占列空间但已离基）。
        let b = Bounds::from_problem(&p);
        assert_eq!(s.solve(&b).status, LpStatus::Optimal);
        assert_rc_batch_matches(&mut s, 2, "solve 后");
        assert_rc_batch_matches(&mut s, 1, "solve 后");
        // 状态 3：add_row 追加两条割——Le 割（满足，松弛入基）与 Ge 割
        // （违反当前解 y=0 → 新建人工列，phase1_pending）。
        s.add_row(&[(1, 1.0)], RelOp::Le, 4.0).ok().unwrap();
        s.add_row(&[(1, 1.0)], RelOp::Ge, 1.0).ok().unwrap();
        assert!(s.phase1_pending);
        assert_rc_batch_matches(&mut s, 2, "add_row 后（含割行+新人工列）");
        // 状态 4：resolve（热启动 + 修复路径）。
        assert_eq!(s.resolve(&b).status, LpStatus::Optimal);
        assert_rc_batch_matches(&mut s, 2, "resolve 后");
        // 状态 5：do_pivot 后（手工驱动一次退化主元改表；不再处于求解器
        // 不变量下的任意中间态恰是批量 vs 逐列一致性最有力的对照）。
        let ncs = s.ncols_static();
        let m = s.basis.len();
        let mut picked = false;
        'outer: for q in 0..ncs {
            if s.in_basis[q] {
                continue;
            }
            for r in 0..m {
                if s.matrix[r][q].abs() >= 1e-6 {
                    s.do_pivot(q, r, 1.0, 0.0, false);
                    picked = true;
                    break 'outer;
                }
            }
        }
        assert!(picked, "应存在可用主元");
        assert_rc_batch_matches(&mut s, 2, "do_pivot 后");
        assert_rc_batch_matches(&mut s, 1, "do_pivot 后");
    }

    // ---- 轴1：稀疏混合表格行索引簿记（规模门后置） ----

    /// 索引簿记一致性内联校验（debug/release 双口径）：稀疏标记行的存活项
    /// 集合 == 该行精确非零集，且低于密度阈；稠密标记行（空索引）无义务。
    fn assert_row_nz_consistent(s: &Simplex) {
        for i in 0..s.matrix.len() {
            let idx = &s.row_nz[i];
            if idx.is_empty() {
                continue;
            }
            let expect: Vec<u32> = s.matrix[i]
                .iter()
                .enumerate()
                .filter(|&(_, &v)| v != 0.0)
                .map(|(j, _)| j as u32)
                .collect();
            let alive: Vec<u32> = idx.iter().copied().filter(|&e| e & NZ_TOMB == 0).collect();
            assert_eq!(alive, expect, "行 {i} 索引簿记与矩阵非零集不符");
            assert!(
                alive.len() * 4 < s.matrix[i].len(),
                "行 {i} 非零数 {} 未低于密度阈（ncols = {}）",
                alive.len(),
                s.matrix[i].len()
            );
        }
    }

    /// 白盒启用轴1簿记（小模型天然低于规模阈；测试直接置位 + 全量重建，
    /// 模拟大表启用态以快速确定性地覆盖簿记路径）。
    fn force_nz_enabled(s: &mut Simplex) {
        s.nz_enabled = true;
        s.nz_rebuild_all();
        assert_row_nz_consistent(s);
    }

    /// （轴1）规模门：小模型天然禁用（旧路径零簿记、row_nz 全空）；大模型
    /// 构建即启用并完成索引重建。
    #[test]
    fn test_nz_scale_gate_disabled_small_enabled_large() {
        // 小模型（34 单元）：禁用 + 全空。
        let p_small = prob(
            Sense::Minimize,
            (0..5).map(|j| var(&format!("x{j}"), 0.0, 8.0)).collect(),
            vec![(0, -1.0)],
            vec![(vec![(0, 1.0), (1, 1.0)], RelOp::Le, 10.0)],
        );
        let mut s = Simplex::new(&p_small, tol()).unwrap();
        assert!(!s.nz_enabled, "阈内模型须禁用轴1（旧路径）");
        assert!(s.row_nz.iter().all(|idx| idx.is_empty()));
        assert_eq!(s.solve(&Bounds::from_problem(&p_small)).status, LpStatus::Optimal);
        assert!(!s.nz_enabled, "求解/加割不跨阈时保持禁用");
        assert!(s.row_nz.iter().all(|idx| idx.is_empty()));

        // 大模型（500 变量 × 700 等式行 = 840000 单元 > 262144）：构建即启用。
        let n = 500;
        let mut vars = Vec::new();
        for j in 0..n {
            vars.push(var(&format!("u{j}"), 0.0, 3.0));
        }
        let mut cons = Vec::new();
        for i in 0..700 {
            let j0 = (i * 13 + 1) % (n - 2);
            cons.push((
                vec![(j0, 1.0), (j0 + 1, 2.0), ((j0 + 37) % n, 1.0)],
                RelOp::Eq,
                4.0,
            ));
        }
        let p_big = prob(Sense::Minimize, vars, (0..n).map(|j| (j, 1.0)).collect(), cons);
        let sb = Simplex::new(&p_big, tol()).unwrap();
        let mb = sb.basis.len();
        let ncolb = sb.ncols_static() + sb.num_art;
        assert!((mb as u64) * (ncolb as u64) > RC_SPARSE_CELL_THRESHOLD);
        assert!(sb.nz_enabled, "超阈模型须启用轴1");
        assert_row_nz_consistent(&sb);
        // 全稀疏行 + 全部行索引非空（3 项结构 + 1 人工列 = 4 nnz ≪ ncols/4）。
        assert!(sb.row_nz.iter().all(|idx| !idx.is_empty()));
    }

    /// （轴1）do_pivot 索引簿记的手算钉（白盒启用态）：fill-in 插入、消元
    /// 恰零惰性墓碑、墓碑复活（nz_dead 同步减）、密度不变量。4 行 30 列模型
    /// （行 nnz=6，6×4=24 < 34 ⇒ 全部稀疏标记），三连手工主元覆盖全部
    /// 簿记分支。
    #[test]
    fn test_row_nz_fill_in_and_tombstone() {
        // min Σ 0（任意目标即可，只手驱 do_pivot）。行 0/2 同支撑同系数
        // （消元恰零），行 1/3 独立支撑（fill-in 目标）。
        let mut cons = Vec::new();
        cons.push((
            (0..5).map(|j| (j, 1.0)).collect::<Vec<_>>(),
            RelOp::Le,
            50.0,
        ));
        cons.push((
            (5..10).map(|j| (j, 1.0)).collect::<Vec<_>>(),
            RelOp::Le,
            50.0,
        ));
        cons.push((
            (0..5).map(|j| (j, 1.0)).collect::<Vec<_>>(),
            RelOp::Le,
            40.0,
        ));
        cons.push((vec![(15, 1.0), (16, 2.0)], RelOp::Le, 30.0));
        let p = prob(
            Sense::Minimize,
            (0..30).map(|j| var(&format!("x{j}"), 0.0, 10.0)).collect(),
            vec![],
            cons,
        );
        let mut s = Simplex::new(&p, tol()).unwrap();
        force_nz_enabled(&mut s);
        // 前置：全部行稀疏标记（nnz=6、ncols=34），松弛列 = 30..33。
        assert_eq!(s.matrix.len(), 4);
        assert!(s.row_nz.iter().all(|idx| !idx.is_empty()), "4 行均应稀疏标记");
        assert_eq!(s.basis, vec![30, 31, 32, 33]);

        // 主元 1：q=0（x0），r=0。行 0 不变（p=1）；f[2]=1 ⇒ 行 2 消元：
        // j∈{0..4} 恰零（1−1）→ 墓碑；j=30（行 0 松弛列）fill-in（0−1=−1）
        // → 插入。行 2 终态存活 = {30, 32}。
        s.do_pivot(0, 0, 1.0, 0.0, false);
        assert_row_nz_consistent(&s);
        let alive2: Vec<usize> = s.row_nz[2]
            .iter()
            .filter(|e| **e & NZ_TOMB == 0)
            .map(|e| (*e & NZ_MASK) as usize)
            .collect();
        assert_eq!(alive2, vec![30, 32], "行 2 消元后存活集");
        assert_eq!(s.matrix[2][30], -1.0);
        assert_eq!(s.matrix[2][32], 1.0);
        for j in 0..5 {
            assert_eq!(s.matrix[2][j], 0.0, "行 2 列 {j} 应恰零");
        }
        assert_eq!(s.nz_dead[2], 5, "行 2 应有 5 个墓碑（j=0..4）");

        // 主元 2：q=30（s0，主元 1 后已离基），r=2（p=−1）。行 2 缩放取反：
        // [30]→+1、[32]→−1（墓碑位 0..4 跳过、值保持 0）；f[0]=1 ⇒ 行 0
        // 消元：j=30 恰零（1−1）→ 墓碑；j=32 fill-in（0−(−1)=+1）→ 插入。
        s.do_pivot(30, 2, 1.0, 0.0, false);
        assert_row_nz_consistent(&s);
        assert_eq!(s.matrix[0][30], 0.0);
        assert_eq!(s.matrix[0][32], 1.0);
        let alive0: Vec<usize> = s.row_nz[0]
            .iter()
            .filter(|e| **e & NZ_TOMB == 0)
            .map(|e| (*e & NZ_MASK) as usize)
            .collect();
        assert_eq!(alive0, vec![0, 1, 2, 3, 4, 32], "行 0 消元后存活集");

        // 主元 3：q=32（s2），r=0（p=1）。f[2]=−1 ⇒ 行 2 消元：
        // j∈{0..4}（行 2 墓碑位、值 0）复活（−(−1)·1=+1，墓碑计数同步减）；
        // j=32 恰零（−1−(−1)·1）→ 墓碑。行 2 终态存活 = {0,1,2,3,4,30}。
        s.do_pivot(32, 0, 1.0, 0.0, false);
        assert_row_nz_consistent(&s);
        let alive2: Vec<usize> = s.row_nz[2]
            .iter()
            .filter(|e| **e & NZ_TOMB == 0)
            .map(|e| (*e & NZ_MASK) as usize)
            .collect();
        assert_eq!(alive2, vec![0, 1, 2, 3, 4, 30], "行 2 复活+墓碑后存活集");
        for j in 0..5 {
            assert_eq!(s.matrix[2][j], 1.0, "行 2 列 {j} 复活值");
        }
        assert_eq!(s.matrix[2][30], 1.0);
        assert_eq!(s.matrix[2][32], 0.0, "行 2 列 32 应恰零（墓碑）");
        assert_eq!(
            s.row_nz[2].iter().filter(|e| **e & NZ_TOMB != 0).count(),
            1,
            "行 2 此时应只有 j=32 一个墓碑"
        );
        // 压缩：墓碑清除后索引即精确非零集（压缩前复活分支已被上方真实
        // 触达——若把压缩提前到主元 3 之前，墓碑已清空、复活变插入，
        // 覆盖就漏了）。
        s.nz_compact();
        assert_eq!(s.row_nz[2], vec![0, 1, 2, 3, 4, 30]);
        assert_eq!(s.nz_dead[2], 0);
        assert_row_nz_consistent(&s);
    }

    /// （轴1，白盒启用态）求解/热启动/加割全流程的索引不变量（累计 > 64
    /// 主元触发周期压缩，add_row 的列插入右移 + 新行索引也被覆盖）。
    #[test]
    fn test_row_nz_invariant_across_solve_cycles() {
        // 50 变量 × 30 行随机稀疏 LP（LCG 确定性），行 nnz≈5+1 松弛 ≪ ncols/4。
        let n = 50;
        let mut lcg = 12345u64;
        let mut rnd = || {
            lcg = lcg.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((lcg >> 33) % 9 + 2) as f64
        };
        let mut vars = Vec::new();
        for j in 0..n {
            vars.push(var(&format!("x{j}"), 0.0, 8.0));
        }
        let obj: Vec<(usize, f64)> = (0..n).map(|j| (j, -rnd())).collect();
        let mut cons = Vec::new();
        for i in 0..30 {
            // 每行 5 个伪随机列（跨行覆盖全部变量空间，保证 fill-in 活跃）。
            let terms: Vec<(usize, f64)> = (0..5)
                .map(|k| (((i * 7 + k * 11) % n), rnd()))
                .collect();
            let rhs: f64 = terms.iter().map(|&(_, a)| a).sum::<f64>() * 0.5;
            cons.push((terms, RelOp::Le, rhs));
        }
        let p = prob(Sense::Minimize, vars, obj, cons);
        let mut s = Simplex::new(&p, tol()).unwrap();
        force_nz_enabled(&mut s);
        let mut b = Bounds::from_problem(&p);
        let sol = s.solve(&b);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert_row_nz_consistent(&s);
        // 注意：solve 内部 rebuild_initial_basis 会按真实规模重判启用位
        // （小模型翻回禁用），此处重新白盒启用以继续覆盖启用态路径。
        force_nz_enabled(&mut s);

        // 热启动序列（界收紧/放松交替，触发 snap + 对偶修复/重建路径）；
        // 中途穿插两次整基重解（冷路径），每次重解后重新启用。
        for step in 0..6 {
            let j = (step * 5) % n;
            if step % 2 == 0 {
                b.ub[j] = 3.0;
            } else {
                b.ub[j] = 8.0;
                b.lb[j] = 1.0;
            }
            assert_eq!(s.resolve(&b).status, LpStatus::Optimal);
            assert_row_nz_consistent(&s);
            if step == 2 || step == 4 {
                assert_eq!(s.solve(&b).status, LpStatus::Optimal);
                force_nz_enabled(&mut s);
            }
        }

        // 加割（Le 违反割 + Ge 割）：割松弛列插入静态段 ⇒ 既有行索引右移，
        // Ge 割再建人工列；resolve 修复后核对。
        let sol = s.resolve(&b);
        let cut0 = sol
            .x
            .iter()
            .enumerate()
            .max_by(|a, c| a.1.partial_cmp(c.1).unwrap())
            .unwrap()
            .0;
        s.add_row(&[(cut0, 1.0)], RelOp::Le, 2.5).unwrap();
        assert_row_nz_consistent(&s);
        assert_eq!(s.resolve(&b).status, LpStatus::Optimal);
        assert_row_nz_consistent(&s);
        s.add_row(&[(cut0, 1.0)], RelOp::Ge, 0.5).unwrap();
        assert_row_nz_consistent(&s);
        assert_eq!(s.resolve(&b).status, LpStatus::Optimal);
        assert_row_nz_consistent(&s);
        assert!(
            s.iterations() > NZ_COMPACT_INTERVAL as u64,
            "累计主元数须超压缩间隔（实际 {}），否则周期压缩路径未被触达",
            s.iterations()
        );
    }

    /// （轴1）tableau_row_terms 与稠密 tableau_row 语义一致：项集 = 稠密行
    /// 非零子集、值逐项相等、升序；越界行同口径 (空, NaN)。禁用态（值扫描
    /// 分支）与白盒启用态（索引分支）双覆盖。
    #[test]
    fn test_tableau_row_terms_matches_dense() {
        let p = prob(
            Sense::Minimize,
            vec![var("x", 0.0, 10.0), var("y", 0.0, 10.0)],
            vec![(0, 1.0), (1, 1.0)],
            vec![
                (vec![(0, 1.0), (1, 1.0)], RelOp::Ge, 2.0),
                (vec![(0, 1.0)], RelOp::Le, 6.0),
            ],
        );
        let mut s = Simplex::new(&p, tol()).unwrap();
        let b = Bounds::from_problem(&p);
        assert!(!s.nz_enabled);
        assert_eq!(s.solve(&b).status, LpStatus::Optimal);
        s.add_row(&[(0, 1.0)], RelOp::Le, 4.0).unwrap();
        assert_eq!(s.resolve(&b).status, LpStatus::Optimal);
        for phase in 0..2 {
            if phase == 1 {
                force_nz_enabled(&mut s);
            }
            for row in 0..s.num_rows() {
                let (dense, dv) = s.tableau_row(row);
                let (terms, tv) = s.tableau_row_terms(row);
                assert_eq!(dv, tv, "行 {row} 基变量值两口接口不一致（phase {phase}）");
                let expect: Vec<(usize, f64)> = dense
                    .iter()
                    .enumerate()
                    .filter(|&(_, &v)| v != 0.0)
                    .map(|(j, &v)| (j, v))
                    .collect();
                assert_eq!(terms, expect, "行 {row} 稀疏项集与稠密非零集不符（phase {phase}）");
            }
        }
        // 越界行：两口接口同口径。
        let (d, dv) = s.tableau_row(999);
        let (t, tv) = s.tableau_row_terms(999);
        assert!(d.is_empty() && t.is_empty());
        assert!(dv.is_nan() && tv.is_nan());
    }

    /// （轴1，自然启用）compute_rc_batch 稀疏散射路径（m×ncols >
    /// RC_SPARSE_CELL_THRESHOLD）与逐列参考 == 相等。状态 A：大等式模型
    /// 初始基（phase-1，人工列 cb=1，全部行稀疏标记）；状态 B：中号 Le
    /// 模型解毕（phase-2，结构基行 cb≠0 + 稠密标记行回退 zip 分支）。
    #[test]
    fn test_rc_sparse_scatter_matches_per_column_large() {
        // —— 状态 A：500 变量 × 700 等式行，m×ncols = 700×1200 = 840000 ——
        let n_a = 500;
        let mut vars = Vec::new();
        for j in 0..n_a {
            vars.push(var(&format!("u{j}"), 0.0, 3.0));
        }
        let mut cons = Vec::new();
        for i in 0..700 {
            let j0 = (i * 13 + 1) % (n_a - 2);
            cons.push((
                vec![(j0, 1.0), (j0 + 1, 2.0), ((j0 + 37) % n_a, 1.0)],
                RelOp::Eq,
                4.0,
            ));
        }
        let pa = prob(
            Sense::Minimize,
            vars,
            (0..n_a).map(|j| (j, 1.0)).collect(),
            cons,
        );
        let mut sa = Simplex::new(&pa, tol()).unwrap();
        let ma = sa.basis.len();
        let ncola = sa.ncols_static() + sa.num_art;
        assert!(
            (ma as u64) * (ncola as u64) > RC_SPARSE_CELL_THRESHOLD,
            "前置：状态 A 须超稀疏散射阈值（{ma}×{ncola}）"
        );
        assert!(sa.nz_enabled, "前置：状态 A 须自然启用轴1");
        assert!(sa.row_nz.iter().all(|idx| !idx.is_empty()));
        assert_rc_batch_matches(&mut sa, 1, "大表初始 phase-1（稀疏散射）");
        assert_rc_batch_matches(&mut sa, 2, "大表初始 phase-2（稀疏散射）");

        // —— 状态 B：600 变量 × 320 行（含 2 条全支撑稠密行），短目标解毕后
        //    phase-2（结构基行 cb≠0：稀疏散射 + 稠密标记行回退 zip 双分支）——
        let n_b = 600;
        let mut vars = Vec::new();
        for j in 0..n_b {
            vars.push(var(&format!("v{j}"), 0.0, 2.0));
        }
        let mut cons = Vec::new();
        for i in 0..320 {
            let terms: Vec<(usize, f64)> = (0..4)
                .map(|k| (((i * 17 + k * 29) % n_b), 1.0 + (k as f64)))
                .collect();
            let rhs: f64 = terms.iter().map(|&(_, a)| a).sum::<f64>() * 3.0;
            cons.push((terms, RelOp::Le, rhs));
        }
        // 2 条全支撑行（nnz = 601 ≥ ncols/4 ⇒ 稠密标记，散射路径回退 zip）。
        for _ in 0..2 {
            cons.push((
                (0..n_b).map(|j| (j, 1.0)).collect::<Vec<_>>(),
                RelOp::Le,
                (n_b as f64) * 1.5,
            ));
        }
        let pb = prob(
            Sense::Minimize,
            vars,
            vec![(0, -1.0), (1, -1.0), (2, -1.0)], // 短目标：少量主元即最优
            cons,
        );
        let mut sb = Simplex::new(&pb, tol()).unwrap();
        assert!(
            sb.row_nz.iter().filter(|idx| idx.is_empty()).count() >= 2,
            "前置：全支撑行须落稠密标记"
        );
        let sol = sb.solve(&Bounds::from_problem(&pb));
        assert_eq!(sol.status, LpStatus::Optimal);
        let mb = sb.basis.len();
        let ncolb = sb.ncols_static() + sb.num_art;
        assert!(
            (mb as u64) * (ncolb as u64) > RC_SPARSE_CELL_THRESHOLD,
            "前置：状态 B 须超稀疏散射阈值（{mb}×{ncolb}）"
        );
        assert!(sb.nz_enabled, "前置：状态 B 须自然启用轴1");
        // 短目标下最优可由非基界翻转达成（全松弛基、cb 全 0）；手工把一个
        // 带成本的结构列主元入基，制造 phase-2 cb≠0 行（与既有
        // test_compute_rc_batch_matches_per_column_multi_state 的状态 5 同法）。
        if !sb.basis.iter().any(|&c| c < sb.num_structural) {
            let mut picked = false;
            'outer: for q in 0..3 {
                if sb.in_basis[q] {
                    continue;
                }
                for i in 0..mb {
                    if !sb.row_nz[i].is_empty() && sb.matrix[i][q].abs() >= 1e-6 {
                        sb.do_pivot(q, i, 1.0, 0.0, false);
                        picked = true;
                        break 'outer;
                    }
                }
            }
            assert!(picked, "应存在含 v0..v2 的稀疏行可供主元");
        }
        assert!(
            sb.basis.iter().any(|&c| c < sb.num_structural),
            "前置：基中须有结构列（phase-2 cb≠0 行存在）"
        );
        assert_row_nz_consistent(&sb);
        assert_rc_batch_matches(&mut sb, 2, "大表解毕 phase-2（稀疏散射+稠密回退）");
    }

    #[test]
    fn test_do_pivot_zero_column_skip_hand_tableau() {
        // 手算小表格：min x+2y+3z s.t. 2x+z<=8（Le）、x+y>=2（Ge 正右端）。
        // 构建期行 0 归一 g=2 ⇒ [1, 0, 0.5 | s1=1]，rhs 4，松弛列 3 入基；
        // 行 1（Ge，残差 2>feas）σ=+1 ⇒ [1, 1, 0 | 0 | −1 | 1]，rhs 2，
        // 人工列 5 入基。列空间 = [x,y,z,s1,s2,a1]（宽 6）。
        let p = prob(
            Sense::Minimize,
            vec![
                var("x", 0.0, 10.0),
                var("y", 0.0, 10.0),
                var("z", 0.0, 10.0),
            ],
            vec![(0, 1.0), (1, 2.0), (2, 3.0)],
            vec![
                (vec![(0, 2.0), (2, 1.0)], RelOp::Le, 8.0),
                (vec![(0, 1.0), (1, 1.0)], RelOp::Ge, 2.0),
            ],
        );
        let mut s = Simplex::new(&p, tol()).ok().unwrap();
        // 初始表格（手算断言；f64 == 口径）。
        assert_eq!(s.matrix[0], vec![1.0, 0.0, 0.5, 1.0, 0.0, 0.0]);
        assert_eq!(s.matrix[1], vec![1.0, 1.0, 0.0, 0.0, -1.0, 1.0]);
        assert_eq!(s.rhs, vec![4.0, 2.0]);
        assert_eq!(s.binv[0], vec![1.0, 0.0]);
        assert_eq!(s.binv[1], vec![0.0, 1.0]);
        assert_eq!(s.basis, vec![3, 5]);

        // 主元 1：q=0（x），r=0，delta=+1，t=2，hit_upper=false。
        // p=1：行 0 不变（[1,0,0.5,1,0,0]）；f[1]=1 ⇒ 行 1 消元
        // [1−1, 1−0, 0−0.5, 0−1, −1−0, 1−0] = [0,1,−0.5,−1,−1,1]，
        // binv[1] = [−1,1]，rhs[1] = 2−2 = 0，rhs[0] = lb[x]+2 = 2。
        // 主元行零列 {y, s2}：跳过是精确无操作（稠密同样得 1、−1）。
        s.do_pivot(0, 0, 1.0, 2.0, false);
        assert_eq!(s.matrix[0], vec![1.0, 0.0, 0.5, 1.0, 0.0, 0.0]);
        assert_eq!(s.matrix[1], vec![0.0, 1.0, -0.5, -1.0, -1.0, 1.0]);
        assert_eq!(s.binv[0], vec![1.0, 0.0]);
        assert_eq!(s.binv[1], vec![-1.0, 1.0]);
        assert_eq!(s.rhs, vec![2.0, 0.0]);
        assert_eq!(s.basis, vec![0, 5]);
        assert!(s.in_basis[0]);
        assert!(!s.in_basis[3]);
        assert!(!s.at_upper[3]);

        // 主元 2：q=2（z），r=1，delta=−1，t=3，hit_upper=true。
        // p=−0.5：行 1 缩放 ×(−2) ⇒ [0,−2,1,2,2,−2]（零列 {x}：稠密缩放
        // 得 −0.0，稀疏跳过保留 +0.0——== 相等）；binv[1] = [2,−2]；
        // f[0]=0.5 ⇒ 行 0 消元 [1,0+1,0.5−0.5,1−1,0−1,0+1] = [1,1,0,0,−1,1]，
        // binv[0] = [1−1, 0+1] = [0,1]，rhs[0] = 2+1.5 = 3.5，
        // rhs[1] = ub[z]−3 = 7。
        s.do_pivot(2, 1, -1.0, 3.0, true);
        assert_eq!(s.matrix[0], vec![1.0, 1.0, 0.0, 0.0, -1.0, 1.0]);
        assert_eq!(s.matrix[1], vec![0.0, -2.0, 1.0, 2.0, 2.0, -2.0]);
        assert_eq!(s.binv[0], vec![0.0, 1.0]);
        assert_eq!(s.binv[1], vec![2.0, -2.0]);
        assert_eq!(s.rhs, vec![3.5, 7.0]);
        assert_eq!(s.basis, vec![0, 2]);
        assert!(!s.in_basis[5]);
        assert!(s.in_basis[2]);
        assert!(s.at_upper[5]);

        // 改表后批量检验数仍与逐列一致（do_pivot 后状态）。
        assert_rc_batch_matches(&mut s, 2, "手算主元后");
        assert_rc_batch_matches(&mut s, 1, "手算主元后");
    }
}
