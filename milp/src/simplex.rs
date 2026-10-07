//! 有界变量两阶段表格单纯形 LP 求解器（WP2 实现）。
//!
//! # 算法结构
//! - 显式 dense 表格 + 有界变量机制（非基变量位于 lb 或 ub，记 `at_upper`），
//!   分支定界只改 [`Bounds`] 不改约束矩阵，热启动代价最低。
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
//!   才允许；入基列无任何阻挡 ⇒ Unbounded。检验数每次定价都由当前表格
//!   重算（不增量维护，无定价漂移）；每 200 次主元用原始数据重算 B⁻¹b；
//!   NaN/Inf 立即按数值失败处理。迭代上限 max(10_000, 20×(m+n)) ⇒
//!   IterationLimit（诚实上报，不冒充 Optimal）。
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
/// 主元计数刷新间隔：每累计本次主元次数重算一次 B⁻¹b 抑制漂移。
const REFRESH_INTERVAL: u64 = 200;
/// Bland 触发的总迭代因子：本次调用总迭代 > BLAND_TOTAL_FACTOR×(m+n) 触发。
const BLAND_TOTAL_FACTOR: u64 = 3;

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
    /// Bland 模式（本次调用内永久生效）。
    bland: bool,
    /// 距上次 B⁻¹b 重算的主元数。
    pivots_since_refresh: u64,
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
            pivots_since_refresh: 0,
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
        self.bounds = bounds.clone();
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
            self.bounds = bounds.clone();
            self.rebuild_initial_basis();
            let st = self.run_two_phase(true);
            self.has_run = true;
            return self.build_solution(st);
        }
        // (1)(2) 非基贴界双向检查 + 传播。
        let prev = self.bounds.clone();
        self.bounds = bounds.clone();
        self.snap_nonbasics(&prev);
        // (3) 基变量越界检查 + Phase-1 欠账检查：任一成立则重建人工列跑
        //     Phase-1 修复，否则直接 Phase-2。欠账 = 上次 Phase-1 未清理即
        //     收场（不可行/上限/数值失败）或 add_row 新建了人工列：此时残基
        //     人工变量行系数活跃，其基值即便被 (1)(2) 的 snap 传播恰好推到
        //     ≤ tol.feasibility，Phase-2 仍会解到该行被零成本人工变量吸收的
        //     错误 LP 上（终检降级 IterationLimit 或谎报 Optimal）。
        let need_repair = self.basics_violated() || self.phase1_pending;
        if need_repair {
            self.rebuild_initial_basis();
        }
        let st = self.run_two_phase(need_repair);
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
        // κ_j = ν_j − Σ_i w_i·M[i][j]（对新行的全旧列消基）。
        let mut kappa = vec![0.0; ncols_old];
        for &(j, c) in &merged {
            kappa[j] = c;
        }
        for i in 0..m {
            if w[i] != 0.0 {
                let wi = w[i];
                let mrow = &self.matrix[i];
                for j in 0..ncols_old {
                    kappa[j] -= wi * mrow[j];
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

        // 新列：松弛列（Eq 无）与人工列（非基时才建）。
        let has_slack = op != RelOp::Eq;
        let has_art = !use_slack_basic;
        let slack_col: Option<usize> = if has_slack {
            let c = n + self.layout.slack.len();
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
        // 旧行补零。
        for row in self.matrix.iter_mut() {
            row.resize(ncols_new, 0.0);
        }
        self.in_basis.resize(ncols_new, false);
        self.at_upper.resize(ncols_new, false);
        // 新表格行 = κ/σ（旧列）| 松弛列系数 slack_sign/σ | 人工列系数 1。
        let mut newrow = vec![0.0; ncols_new];
        for j in 0..ncols_old {
            newrow[j] = kappa[j] / sigma_basic;
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
        self.matrix.push(newrow);
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
    }

    // ===================== 内部：定价 =====================

    /// 定价：双侧合格性扫描。Bland 模式取最小下标合格列；Dantzig 模式取
    /// 得分（贴下界 −rc / 贴上界 +rc）最大者、并列取最小列下标。
    /// 人工列与固定列（ub−lb ≤ feasibility）永不入基。
    fn price(&mut self, phase: u8) -> Pricing {
        let m = self.basis.len();
        let ncs = self.ncols_static();
        let mut c_b = vec![0.0; m];
        for i in 0..m {
            c_b[i] = self.cost_of(self.basis[i], phase);
        }
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
            let mut rc = self.cost_of(j, phase);
            for i in 0..m {
                if c_b[i] != 0.0 {
                    rc -= c_b[i] * self.matrix[i][j];
                }
            }
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
    fn ratio_test(&self, q: usize, delta: f64) -> Step {
        let m = self.basis.len();
        // (t, 行, hit_upper, |主元|, 破键列下标)；行 = usize::MAX 表示界翻转。
        let mut cands: Vec<(f64, usize, bool, f64, usize)> = Vec::new();
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
            cands.push((t_i, i, hit_upper, a.abs(), self.basis[i]));
        }
        // 入基变量自身界翻转（从所贴界移动到对端界）。
        let (lb, ub) = self.col_bounds(q);
        if ub.is_finite() {
            cands.push(((ub - lb).max(0.0), usize::MAX, false, 0.0, q));
        }
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
        for c in &cands {
            if c.0 < t_min {
                t_min = c.0;
            }
        }
        let band = self.tol.feasibility * (1.0 + t_min.abs());
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
    fn do_pivot(&mut self, q: usize, r: usize, delta: f64, t: f64, hit_upper: bool) {
        let m = self.basis.len();
        let ncols = self.matrix[r].len();
        let p = self.matrix[r][q];
        // 行 r 缩放。
        for j in 0..ncols {
            self.matrix[r][j] /= p;
        }
        self.matrix[r][q] = 1.0;
        for k in 0..m {
            self.binv[r][k] /= p;
        }
        // 出基列与入基旧值。
        let leaving = self.basis[r];
        let entering_old = if delta > 0.0 {
            self.col_bounds(q).0
        } else {
            self.col_bounds(q).1
        };
        // 消元（f_k 用未更新的 matrix[k][q]）。
        let f: Vec<f64> = (0..m)
            .map(|k| if k == r { 0.0 } else { self.matrix[k][q] })
            .collect();
        for k in 0..m {
            if k == r {
                continue;
            }
            let fk = f[k];
            if fk != 0.0 {
                for j in 0..ncols {
                    self.matrix[k][j] -= fk * self.matrix[r][j];
                }
                for c in 0..m {
                    self.binv[k][c] -= fk * self.binv[r][c];
                }
                self.rhs[k] -= delta * t * fk;
            }
            self.matrix[k][q] = 0.0; // 数值卫生
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
        let limit: u64 = 10_000.max(20u64.saturating_mul((m + n) as u64));
        let bland_total: u64 = BLAND_TOTAL_FACTOR.saturating_mul((m + n) as u64);
        loop {
            if self.iter_call >= limit {
                return PhaseOutcome::LimitReached;
            }
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
            if self.degenerate_streak > DEGENERATE_STREAK_LIMIT || self.iter_call > bland_total {
                self.bland = true;
            }
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

    /// 非基列贴界双向检查（相对上次有效界 `prev`）：
    /// * 收紧越界：贴上界且新 ub < 当前值 → snap 到新 ub；贴下界且新 lb >
    ///   当前值 → snap 到新 lb。
    /// * 松弛失贴：当前值严格位于新 (lb,ub) 内部 → snap 到较近界；等距并列
    ///   按 Phase-2 检验数符号选侧（rc<0 选下界、rc≥0 选上界）。
    /// 每个 snap 的 Δx_j 按 B⁻¹A_j 传播更新全部基变量值。
    fn snap_nonbasics(&mut self, prev: &Bounds) {
        let n = self.num_structural;
        let feas = self.tol.feasibility;
        for j in 0..n {
            if self.in_basis[j] {
                continue;
            }
            let old_v = if self.at_upper[j] {
                prev.ub[j]
            } else {
                prev.lb[j]
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
        self.pivots_since_refresh = 0;
    }

    /// 构造 LpSolution：提取 x（基本值在容差内夹回界）、按 x 重算目标
    /// （值与目标自洽）、Optimal 时做最终原始可行性/人工变量校验并计算对偶。
    fn build_solution(&mut self, mut status: LpStatus) -> LpSolution {
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
}
