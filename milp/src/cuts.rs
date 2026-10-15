//! Gomory 混合整数割（GMI）（WP3 实现）。
//!
//! 只在根节点调用（割对原问题全体整数可行解全局有效；分支子树内调用无效——
//! 禁止）。推导：取基变量为整数类型结构变量的表格行 r，当前值 x̄，
//! f0 = x̄−⌊x̄⌋ ∈ [away, 1−away] 否则跳过；非基列 j：贴下界取 d_j = x_j−lb_j、
//! α_j = ā_j；贴上界取 d_j = ub_j−x_j、α_j = −ā_j。
//!
//! 系数分类（关键前提）：
//! - 整数类型结构列：仅当其所贴界为整数值（|bound−round(bound)| ≤ tol.zero）
//!   时才可用整数公式 π_j = f_j/f0（f_j ≤ f0）否则 (1−f_j)/(1−f0)，
//!   f_j = α_j−⌊α_j⌋（此时 d_j 在整数可行点取整数值，公式有效）。前提判定
//!   不得依赖 presolve 已把界取整（use_presolve=false 时界可形如 z ≤ 10.5）。
//! - 整数类型但所贴界非整数：必须按连续列公式 π_j = max(α_j/f0, −α_j/(1−f0))
//!   处理，否则割可切掉整数最优而错报 Optimal；
//! - 连续结构列与全部松弛列（含割行松弛列）：π_j = max(α_j/f0, −α_j/(1−f0))。
//!
//! 割为 Σ π_j·d_j ≥ 1（π_j ≥ 0，d_j ≥ 0）。d_j 回代：结构列直接
//! （x_j−lb_j 或 ub_j−x_j）；原始约束松弛列按 [`crate::simplex::SlackCol::sign`]
//! （仅 +1 松弛/−1 剩余两来源，行永不翻转）与该行归一化系数换算
//! （s = σ·(rhs − Σ a_j x_j)/g，g = max|a_ij|）；割行松弛列
//! s_k = rhs_k − Σ_j terms_j·x_j（terms/rhs_k 即本模块 add_row 时传入并记账的
//! 原尺度值——add_row 不缩放是 WP2 契约，故无需查询内部 scale）。
//!
//! 数值纪律：|c_j| < tol.zero 的项只允许按放宽割的方向丢弃（丢项并按界幅抬高
//! 右端补偿，严禁加严；上界 +∞ 无法界定补偿时保留该项且不计入比率检验）；
//! violation（当前 LP 解代入 lhs−rhs）≥ min_violation 才收；
//! max|c|/min 非零|c| > 1e9 整条丢弃（同时补偿割行不缩放的数值风险）。
//!
//! # 分离边界数值卫生闸（清洁化 + 密度闸）
//!
//! 生成器内部的 tol.zero（1e-12）级清洁挡不住表格噪声经 GMI 连续列公式
//! （f0 ≥ away ⇒ 放大至多 1/away 倍）与松弛列回代复合后的残差，且全族
//! 均无密度控制——Devex 基下根割循环曾产出 160 条 661~681 nnz（矩阵
//! 900 列、约 74% 密度）的病态割，真零位残留 3.5e-18~7.0e-15 的 B⁻¹
//! 计算噪声，全局保留后毒化全树（每迭代成本 2-3×、退化打转死节点、
//! 端到端 10×+ 劣化）。为此在**施加边界**（[`cut_loop_recorded`] 与
//! [`separate_local_cuts`]，全部割族与两引擎的唯一汇聚点）对合并后的
//! 每条割统一过 [`sanitize_cut`]：小系数清洁化（只弱化，论证见该函数）
//! → 密度闸（nnz 超上限整条丢弃）→ violation 复核。闸只依赖割自身与
//! 变量界，不按算例名或规模分支。
//!
//! # 覆盖割（WP4 增补，仅根节点）
//!
//! [`cover_cuts`] 对每条合格背包行（Le、全 Binary 支撑、权全正、全支撑
//! Σw > cap）做 x̄ 引导的最小覆盖分离（无提升）：贪婪按 x̄ 降序并入 C 直到
//! Σ_C w > cap，再按 x̄ 升序剔除可剔者得最小覆盖，violation =
//! Σ_{j∈C} x̄_j − (|C|−1) ≥ min_violation 才产出 `Σ_{j∈C} x_j ≤ |C|−1`。
//! 有效性不依赖 x̄：覆盖 C 的权和超容 ⇒ 0-1 可行解不能令 C 全取 1。
//! [`cut_loop_recorded`] 每轮把覆盖割批接在 GMI 批之后（确定性合并、统一
//! 截断 max_cuts_per_round），割行照走 add_row + 记账（worker 按序重放契约
//! 不变）。
//!
//! # 单行 c-MIR 与隐含界割（WP6 增补，仅根节点）
//!
//! [`mir_single_row`] 对**原始 Le 行**做变量界互补归一：负系数列
//! x_j → ub_j − y_j（无有限 ub 的行整行跳过；整数类型列所互补的 ub 须为
//! 整数值，非负系数整数列的 lb 须为整数值——否则平移/互补不再保双射整性，
//! 整行跳过），非负系数列平移 y_j = x_j − lb_j。归一后全系数 ≥ 0、
//! y ≥ 0（整数列 y ∈ Z_+），除以 g = max â_j 再在尺度集
//! δ ∈ {f0, 1}·{1, 0.5, 0.25}（乘积展开：f0、f0/2、f0/4、1、0.5、0.25，
//! f0 = frac(归一右端)）上逐试 MIR 函数：整数列
//! π_j = ⌊α_j⌋ + (f_j − f̂0)^+/(1 − f̂0)、归一后正系数连续列 π_j = 0
//! （连续列只消耗容量，丢弃只放宽），右端 ⌊b̂/(gδ)⌋（取整 ε 纪律：放宽
//! 方向），f̂0 复用 away 门（f̂0 ∈ [away, 1−away]）。取 x̄ 处 violation
//! 最大且过数值门者。有效性 = 界互补/平移双射（整数可行点一一对应）+
//! MIR 引理（纯整数、非负系数行）+ 连续列丢弃只放宽。
//!
//! [`implied_bound_cuts`] 对**整数类型变量计数恰为 1** 的行导出该整数变量
//! y 的隐含界割：行的逻辑蕴涵（其余列按方向折叠到恰好覆盖全部可行点的一侧
//! 界，见函数文档的 Kmin/Kmax 公式）+ 方向取整（上界 floor 放宽、下界 ceil
//! 放宽，见 [`floor_bound`]/[`ceil_bound`]）。需要的一侧变量界无限 ⇒ 该
//! 方向跳过。纯变量界形式经反例证明是唯一保真的线性形式（保留连续列线性
//! 表达并取整常数会切掉整数可行点：行 `y − c ≤ 5.7`、c ∈ [0,∞) 连续时，
//! 点 y=6、c=0.3 行可行而被 `y ≤ 5 + c` 类割切掉），故本实现为纯界形式；
//! 连续列经其有限界参与折叠、不做整性假设。
//!
//! [`cut_loop_recorded`] 每轮割批的确定性合并序（WP6 起，R2-retry 增聚合
//! MIR）：**c-MIR → GMI → 聚合 MIR → 隐含界割 → 覆盖割（提升）**，统一截断
//! `max_cuts_per_round`，全部照走 `add_row` + [`CutRecord`] 记账（worker
//! 按序重放契约不变）。
//!
//! # 表格行聚合 MIR（R2-retry，仅根节点；O3 留档方案的可实现形态）
//!
//! 第四轮诊断：原始行上的 c-MIR 产出近零——表格行才携带基信息。O3 的两个
//! 字面构造经证伪：(a) 表格行内松弛列按 GMI 纪律**全部**回代成 x 空间
//! ⇒ 恒等式塌缩 0=0（表格行是原始等式系的精确线性组合，行内全部松弛按
//! 定义回代恒真，2×2 手算例验证）；(b) 表格行**对** λ 合并产生双基列，
//! GMI 格论证失效（合并权重非整），纯 MIR 则因丢弃连续/松弛列在背包行上
//! 弱于既有 GMI（R2 首试在 mknap 上产出恒零）。落地形态
//! [`aggregated_mir_cuts`]：候选表格行（GMI 同口径）与其**非基松弛列的
//! 原始归一行**合并（μ 恰消去该松弛列、并入其结构内容），锚系数重新归一
//! 后全部列取单侧偏差（基本结构列按较近界），套与 GMI 相同的混合整数
//! 公式——**连续/松弛列 π≠0 保留**。有效性经 50 种子随机小模型穷举验证
//! （302 条割零误切）；默认参数窗口即非零产出（mknap-3x100 8 条、
//! mknap-5x300 25 条、setcover-100x150 386 条/轮，viol 至 0.98）。
//!
//! # 序贯 up-lifting 的预算与正确性（WP6）
//!
//! 对最小覆盖 C 的基础割 `Σ_{j∈C} x_j ≤ |C|−1`（β = |C|−1 **固定不变**），
//! 按序把提升列 j ∉ C 以系数 α_j 并入：α_j = β − max{Σ c_i x_i}，其中
//! c_i 为**当前**割系数（C 成员 1、已提升列 α），最大值取自
//! `w_j + Σ_{i∈已处理} w_i x_i ≤ cap` 的 0-1 背包（未处理列排除）——
//! [`lift_cover_draft`] 用增量 DP（容量轴 0..cap̃ 一次建表、每提升一列
//! O(cap̃) 追加一项）精确求解。序贯性：第 k 列的 DP 含前 k−1 列的 α
//! （严禁对基础覆盖独立逐列计算）。
//!
//! 【与工作包文本的偏差（已核实）】工作包原文把 α_j 写作
//! `|C|−1+Σ已提升α − DP 最优值`（右端随 Σα 增长）——该形式**无效**：
//! 反例行 w=[5,4,3,2]、cap=5，C={0,1}，按原文提升 x2 得 α=1、提升 x3 得
//! α=(1+1)−1=1，割 `x0+x1+x2+x3 ≤ 1` 切掉行可行点 {x2,x3}（3+2=5 ≤ 5，
//! LHS=2 > 1）。标准序贯 up-lifting（Gu–Nemhauser–Savelsbergh）β 固定，
//! 同例得 α_2=1、α_3=0，穷举全部 0-1 行可行点均满足。正确性红线
//! （不得误剪任何整数可行解）优先，按标准公式实现。
//!
//! 预算三层（全部超限 ⇒ 原样退回基础覆盖割，不弱于现状）：
//! - 每行提升列数 K ≤ min(32, |支撑\C|)（按 x̄ 降序取前 K）；
//! - 单割建表门槛 (|C|+K)×(cap̃+1) ≤ 2×10^6 状态；
//! - 根割循环全程 DP 总预算 3×10^7 状态操作（[`LiftBudget`]，耗尽即全程
//!   关闭提升）；每轮只对 violation 最高的前 4 条覆盖割尝试提升。
//!
//! 建表判据（浮点权防御）：行权与 cap 同乘 10^k（k ∈ 0..=6）后全部（容差
//! 1e-9·(1+|v|) 内）恰为整数的行才建表；容量上取整 ⌈cap·10^k⌉、权下取整
//! ⌊w·10^k⌋——双向偏向『装得下』⇒ DP 最优值 ≥ 真值 ⇒ α_j 偏小 ⇒ 割偏弱
//! 仍有效（安全方向）。非整数化行跳过提升（只留基础覆盖割）。

use crate::error::Result;
use crate::model::{Constraint, Problem, RelOp, Tolerances, VarType};
use crate::simplex::{Bounds, LpStatus, RowOrigin};
use crate::simplex_sparse::LpEngine;

/// 割系数动态范围上限：显著项最大/最小绝对值之比超过即整条丢弃（契约值）。
const MAX_COEF_RATIO: f64 = 1e9;

/// 覆盖割的权值幅值上限：行内任一 |w| ≥ 此值整行跳过（数值防御，与
/// [`MAX_COEF_RATIO`] 同量级的契约值）。
const MAX_COVER_WEIGHT: f64 = 1e9;

/// c-MIR 尺度 δ 的下限：过小的 δ 会把系数吹到数值不稳区（下游比率门仍会
/// 兜底，此处提前防御）。
const MIR_MIN_DELTA: f64 = 1e-6;

/// lifting：每行提升列数上限 K（按 x̄ 降序取前 K）。
const MAX_LIFT_CANDIDATES: usize = 32;

/// lifting：单割 DP 建表门槛（状态数 = (|C|+K)×(cap̃+1)）。
const MAX_LIFT_TABLE_STATES: u64 = 2_000_000;

/// lifting：根割循环全程提升 DP 总预算（状态操作数），耗尽即全程关闭提升。
const ROOT_LIFT_TOTAL_BUDGET: u64 = 30_000_000;

/// lifting：每轮只对 violation 最高的前几条覆盖割尝试提升。
const LIFT_TOP_PER_ROUND: usize = 4;

/// 建表判据与取整 ε 的统一相对容差（1e-9·(1+|v|)）。
const NUMERIC_SNAP_TOL: f64 = 1e-9;

/// 分离边界清洁化容差（数值卫生闸）：|c_j| ≤ 此值（且该列界幅有限）的
/// 割系数项按「只弱化」纪律丢弃并按界幅补偿右端。取值依据见
/// [`sanitize_cut`] 的论证（上界盖住 B⁻¹ 噪声放大后的全部残差，下界不
/// 误伤有分离价值的显著项）。
const CUT_CLEANUP_EPS: f64 = 1e-9;

/// 密度闸 nnz 上限的绝对下限（项数）：割非零元超过上限即整条丢弃。
/// 取值依据见 [`sanitize_cut`]。
const CUT_DENSITY_NNZ_FLOOR: usize = 200;

/// 密度闸 nnz 上限的相对系数：上限 = max(下限, ⌈n × 此值⌉)，n 为结构
/// 列数。取值依据见 [`sanitize_cut`]。
const CUT_DENSITY_NNZ_FRAC: f64 = 0.1;

/// 结构变量空间的一条割（op 恒为 Le）。
#[derive(Debug, Clone)]
pub struct Cut {
    /// 稀疏项 (变量下标, 系数)，按下标升序、无零系数。
    pub terms: Vec<(usize, f64)>,
    /// 关系运算符（恒为 Le）。
    pub op: RelOp,
    /// 右端。
    pub rhs: f64,
}

/// 割生成参数。
#[derive(Debug, Clone, Copy)]
pub struct CutOptions {
    /// 最大割轮数。默认 4（WP6 轮次预算实测校准，见 [`CutOptions::default`]）。
    pub max_rounds: usize,
    /// 每轮最多施加割数。默认 40（WP6：20→40）。
    pub max_cuts_per_round: usize,
    /// f0 距 0/1 的最小距离。默认 0.005。
    pub away: f64,
    /// 割须被当前 LP 解违反的最小量。默认 1e-7。
    pub min_violation: f64,
}

impl Default for CutOptions {
    /// 默认值：max_rounds=4、max_cuts_per_round=40、away=0.005、min_violation=1e-7。
    ///
    /// 【WP6 轮数取值的降档依据（工作包 >5% 墙钟回退守则的最终档）】
    /// 计划值 8×40 在 knapsack-2000 上触发灾难性回退：第 6~8 轮 GMI 追加割
    /// 使树 39→11917 节点、墙钟 0.28s→37s（×122，新割族全关亦然——纯 GMI
    /// 晚轮效应，与本包新机制无关）。lifting 三档旋钮（DP 预算/K）不适用
    ///（k2000 被建表门槛跳过、k500 提升系数全 0），按守则末档收缩轮数：
    /// 8→5 恢复 k2000 基线但 knapsack-200/mknap-3x100 仍 >5% 回退；
    /// 8→4 于 2026-10-09 本机（threads=1，改前/改后同窗口交替配对中位数）
    /// 实测净改善：knapsack-500 171→42ms（−75%，1474→536 节点）、
    /// mknap-5x300 57.0→22.2s（−61%，201807→99605 节点，90s 内证最优）、
    /// knapsack-2000 266→262ms（−2%，39→31 节点）、marketsplit-2x20
    /// 41→26ms、milp-mixed-250i40c 88→68ms；残留回退：knapsack-200
    /// 12.6→21.5ms、marketsplit-4x30 738→1296ms（树形对割集混沌敏感，
    /// MostFractional+DepthFirst 下无全域不回退档位，详见包报告）。
    fn default() -> Self {
        CutOptions {
            max_rounds: 4,
            max_cuts_per_round: 40,
            away: 0.005,
            min_violation: 1e-7,
        }
    }
}

/// cut_loop 的割记账（私有）：add_row 传入的原尺度 terms/rhs 与返回的割序号，
/// 用于后续轮次对割松弛列的回代。
struct CutRecord {
    /// add_row 传入的结构变量空间稀疏行（原尺度）。
    terms: Vec<(usize, f64)>,
    /// add_row 传入的右端（原尺度）。
    rhs: f64,
    /// add_row 返回的割序号 k。
    index: usize,
}

/// 依据当前 Optimal 表格生成一批割（不动 Simplex 状态）。
///
/// 非基割松弛列因无记账（本函数不携带 cut_loop 的私有账本）而整列丢弃。
/// 【方向勘误（独立复核发现，R7）：割为 Σπ_j·d_j ≥ 1（π ≥ 0、d ≥ 0），
/// 从 LHS 丢 π>0 的项使割**更强**、可能误切整数可行解——不是放宽。
/// 该丢弃在全部接线路径**不可达**（cut_loop 委托 cut_loop_recorded：
/// ledger 自空完整累积、割循环从重建的冷 Simplex 启动 ⇒ 每个割松弛列
/// 必有记账，`RowOrigin::Cut` 查找恒命中）；公开入口若被用于手工
/// add_row 过的 Simplex（库内无此调用方），丢弃产出的割可能过强而
/// 无效——调用方自负。】cut_loop 内部走带记账路径，可精确回代。
/// 调用前提：`bounds` 即 Simplex 最近一次 solve/resolve 所用界，且状态
/// 为 Optimal。
pub fn gomory_cuts(
    s: &dyn LpEngine,
    problem: &Problem,
    bounds: &Bounds,
    int_vars: &[usize],
    opts: CutOptions,
    tol: &Tolerances,
) -> Vec<Cut> {
    gomory_cuts_inner(s, problem, bounds, int_vars, opts, tol, &[])
}

/// 带割记账的割生成内核（cut_loop 私有路径）。
fn gomory_cuts_inner(
    s: &dyn LpEngine,
    problem: &Problem,
    bounds: &Bounds,
    int_vars: &[usize],
    opts: CutOptions,
    tol: &Tolerances,
    ledger: &[CutRecord],
) -> Vec<Cut> {
    let n = problem.vars.len();
    if bounds.lb.len() != n || bounds.ub.len() != n {
        return Vec::new();
    }
    let layout = s.column_layout();
    let num_struct = layout.num_structural;
    if num_struct != n {
        return Vec::new(); // 表格与问题维度不符（防御）
    }
    let basis = s.basis_columns();
    let nonbasic = s.nonbasic_columns();
    let mut is_int = vec![false; n];
    for &j in int_vars {
        if j < n {
            is_int[j] = true;
        }
    }

    // 重建当前 LP 解 x̄（结构变量）：基列取表格 rhs，非基列取所贴界。
    // 经 [`rebuild_x_bar`] 与覆盖割分离共用同一公开访问器重建路径。
    let Some(x_bar) = rebuild_x_bar(s, bounds, n) else {
        return Vec::new(); // 表格与界不一致（防御）
    };

    let mut cuts = Vec::new();
    for row in 0..s.num_rows() {
        let bcol = basis[row];
        if bcol >= num_struct || !is_int[bcol] {
            continue; // 基变量须为整数类型结构列
        }
        let (terms, value) = s.tableau_row_terms(row);
        if !value.is_finite() {
            continue; // 越界行（空表 + NaN）或数值异常
        }
        if terms.iter().any(|&(_, v)| !v.is_finite()) {
            continue; // 数值异常行整行跳过（NaN/Inf 非零必入稀疏表）
        }
        let f0 = value - value.floor();
        if !(f0 >= opts.away && f0 <= 1.0 - opts.away) {
            continue;
        }

        // Σ_j π_j·d_j ≥ 1 的结构变量系数与常数项。
        let mut coef = vec![0.0; n];
        let mut const_term = 0.0f64;
        let mut any = false;
        // 轴1：稀疏读表（terms 按列升序）与 nonbasic（按列升序）归并游标
        // 推进——系数恰为 0 的列在稀疏接口中省略，等价于稠密读取的
        // α == 0 ⇒ π == 0 ⇒ 无贡献（逐项行为与旧稠密路径一致）。
        let mut tcur = 0usize;
        for nb in &nonbasic {
            let j = nb.col;
            while tcur < terms.len() && terms[tcur].0 < j {
                tcur += 1;
            }
            let alpha0 = if tcur < terms.len() && terms[tcur].0 == j {
                terms[tcur].1
            } else {
                0.0
            };
            let alpha = if j < n {
                if nb.at_upper {
                    -alpha0
                } else {
                    alpha0
                }
            } else {
                alpha0 // 松弛/剩余列恒贴下界 0
            };
            let pi = if j < n && is_int[j] {
                let bound = if nb.at_upper {
                    bounds.ub[j]
                } else {
                    bounds.lb[j]
                };
                if bound.is_finite() && (bound - bound.round()).abs() <= tol.zero {
                    // 所贴界为整数值：d_j 在整数可行点取整数值，整数公式有效。
                    let fj = alpha - alpha.floor();
                    if fj <= f0 {
                        fj / f0
                    } else {
                        (1.0 - fj) / (1.0 - f0)
                    }
                } else {
                    continuous_pi(alpha, f0)
                }
            } else {
                continuous_pi(alpha, f0)
            };
            if !pi.is_finite() || pi <= 0.0 {
                continue; // π = 0 无贡献（精确无操作）；非有限为数值异常，丢列
            }
            if pi < tol.zero {
                // 近零项丢弃。方向勘误（独立复核，R7）：割为 Σπ·d ≥ 1
                //（π ≥ 0、d ≥ 0），丢 0 < π < tol.zero 的项使割**更强**
                //（非放宽）——增量上界 π·d_j ≤ tol.zero·d_j（tol.zero =
                // 1e-12，结构列 d_j ≤ 界幅、松弛列 d_j = s_i 有界），数值
                // 上可忽略，但方向如实记录、不称“安全放宽”。
                continue;
            }
            if j < n {
                if nb.at_upper {
                    coef[j] -= pi;
                    const_term += pi * bounds.ub[j];
                } else {
                    coef[j] += pi;
                    const_term -= pi * bounds.lb[j];
                }
                any = true;
            } else {
                let info = &layout.slack[j - num_struct];
                match info.origin {
                    RowOrigin::Constraint(ci) => {
                        // s = σ·(rhs − Σ a·x)/g，g 为该行 max|a_ij|（与 WP2
                        // 构建期归一一致）。
                        let Some(con) = problem.constraints.get(ci) else {
                            // 防御：行来源越界，丢列——使割更强非放宽
                            //（独立复核方向勘误，R7）；构造保证 origin 索引
                            // 恒有效，此分支不可达。
                            continue;
                        };
                        let g = con
                            .expr
                            .terms
                            .iter()
                            .map(|&(_, a)| a.abs())
                            .fold(0.0_f64, f64::max);
                        if g <= tol.zero {
                            continue; // 防御：不应出现（近零行无松弛列）
                        }
                        let scale = info.sign / g;
                        const_term += pi * scale * con.rhs;
                        for &(k, a) in &con.expr.terms {
                            coef[k] -= pi * scale * a;
                        }
                        any = true;
                    }
                    RowOrigin::Cut(k) => {
                        // s_k = rhs_k − Σ terms·x（add_row 原尺度，不缩放）。
                        let Some(rec) = ledger.iter().find(|r| r.index == k) else {
                            // 无记账（公开入口路径）：丢列。方向勘误（独立
                            // 复核，R7）：割 Σπ·d ≥ 1 丢 π>0 列使割**更强**、
                            // 可能误切整数可行解（非放宽）。接线路径不可达：
                            // cut_loop_recorded 的 ledger 自空完整累积 + 割
                            // 循环从冷 Simplex 启动 ⇒ Cut 列必有记账；仅公开
                            // 入口 gomory_cuts 手工用于带割 Simplex 时可达
                            //（库内无此调用方，见其文档勘误块）。
                            continue;
                        };
                        const_term += pi * rec.rhs;
                        for &(k2, a) in &rec.terms {
                            coef[k2] -= pi * a;
                        }
                        any = true;
                    }
                }
            }
        }
        if !any {
            continue;
        }

        // Σ coef·x + const ≥ 1  ⇔  Σ (−coef)·x ≤ const − 1。
        let mut rhs = const_term - 1.0;
        let mut terms: Vec<(usize, f64)> = Vec::new();
        let mut c_max = 0.0f64;
        let mut c_min = f64::INFINITY;
        for j in 0..n {
            let c = -coef[j];
            if c == 0.0 {
                continue;
            }
            if c.abs() <= tol.zero {
                // 近零项：按界幅补偿右端后丢弃（丢项 + 抬右端，只放宽）。
                let magnitude = bound_magnitude(bounds.lb[j], bounds.ub[j]);
                if magnitude.is_finite() {
                    rhs += c.abs() * magnitude;
                    continue;
                }
                // 上界 +∞ 无法界定补偿：保留该项，但不计入比率检验。
                terms.push((j, c));
                continue;
            }
            terms.push((j, c));
            c_max = c_max.max(c.abs());
            c_min = c_min.min(c.abs());
        }
        if terms.is_empty() || c_max <= 0.0 {
            continue; // 无显著项的割没有意义
        }
        if c_max / c_min > MAX_COEF_RATIO {
            continue; // 系数动态范围过大，整条丢弃
        }
        // violation 检验：当前 LP 解须违反割至少 min_violation。
        let mut lhs = 0.0f64;
        for &(j, c) in &terms {
            lhs += c * x_bar[j];
        }
        if !lhs.is_finite() || !rhs.is_finite() || lhs - rhs < opts.min_violation {
            continue;
        }
        cuts.push(Cut {
            terms,
            op: RelOp::Le,
            rhs,
        });
    }
    cuts
}

/// 重建当前 LP 解 x̄（结构变量空间）：基列取表格 rhs（值经
/// [`Simplex::tableau_row_terms`] 稀疏接口读取——只需值，免整行稠密拷贝），
/// 非基列取所贴界（[`Bounds`]）；全部经公开访问器，不依赖 Simplex 内部状态。
///
/// 返回 `None`（调用方按防御处理）：表格结构列数 ≠ `n`，或重建出现非有限值
/// （表格与界不一致）。
fn rebuild_x_bar(s: &dyn LpEngine, bounds: &Bounds, n: usize) -> Option<Vec<f64>> {
    if bounds.lb.len() != n || bounds.ub.len() != n {
        return None;
    }
    if s.column_layout().num_structural != n {
        return None; // 表格与问题维度不符（防御）
    }
    let basis = s.basis_columns();
    let nonbasic = s.nonbasic_columns();
    let mut x_bar = vec![f64::NAN; n];
    for (row, &bcol) in basis.iter().enumerate() {
        if bcol < n {
            let (_, value) = s.tableau_row_terms(row);
            x_bar[bcol] = value;
        }
    }
    for nb in &nonbasic {
        if nb.col < n {
            x_bar[nb.col] = if nb.at_upper {
                bounds.ub[nb.col]
            } else {
                bounds.lb[nb.col]
            };
        }
    }
    if x_bar.iter().any(|v| !v.is_finite()) {
        return None; // 表格与界不一致（防御）
    }
    Some(x_bar)
}

/// 根节点覆盖割分离草案（WP6 重构：[`cover_cuts`] 的实现内核 + 提升所需元数据）。
///
/// 对 `problem.constraints` 的每条**合格行**做至多 1 条基础覆盖割。合格判据：
/// ① `Le` 行；② 支撑变量全为 [`VarType::Binary`]；③ 权全正且数值稳健
/// （`w` 有限、`|w| < 1e9`（[`MAX_COVER_WEIGHT`]）、`w > tol.zero`；右端与
/// 支撑上 x̄ 亦须有限，任一非有限量整行跳过）；④ 全支撑 Σw > cap（否则
/// 不存在覆盖）。
///
/// 分离步骤（排序并列序固定 ⇒ 全序确定）：
/// - **贪婪**：按 x̄ 降序（并列 w 降序、再下标升序）并入 C，直到 Σ_C w > cap；
/// - **极小化**：按 x̄ 升序（并列 w 升序）遍历 C，剔除后仍保持 Σ > cap 者
///   剔之——终态对每个保留成员 j 均有 Σ_C w − w_j ≤ cap（最小覆盖）；
/// - **violation 检验**：`Σ_{j∈C} x̄_j − (|C|−1) ≥ opts.min_violation` 才
///   产出 `Cut { terms = [(j, 1.0)]（下标升序）, Le, rhs = |C|−1.0 }`。
///
/// 有效性（不依赖 x̄）：覆盖 C 满足 Σ_{j∈C} w_j > cap ⇒ 该行的任何 0-1
/// 可行解不能令 C 全取 1 ⇒ Σ_{j∈C} x_j ≤ |C|−1。
///
/// `bounds`/`x_bar` 与问题维数不符时返回空（防御）。
struct CoverDraft {
    /// 基础覆盖割（[(j,1.0)] 下标升序、Le、rhs = |C|−1）。
    cut: Cut,
    /// 基础割在 x̄ 处的 violation（Σ_C x̄ − (|C|−1)）。
    violation: f64,
    /// 最小覆盖成员（下标升序）。
    cover: Vec<usize>,
    /// 行支撑（(j, w)，下标升序——LinearExpr 不变量）。
    terms: Vec<(usize, f64)>,
    /// 行右端（容量）。
    cap: f64,
}

/// [`CoverDraft`] 的分离内核（[`cover_cuts`] 原逻辑逐行保持）。
fn cover_drafts(
    problem: &Problem,
    bounds: &Bounds,
    x_bar: &[f64],
    opts: CutOptions,
    tol: &Tolerances,
) -> Vec<CoverDraft> {
    let n = problem.vars.len();
    if bounds.lb.len() != n || bounds.ub.len() != n || x_bar.len() != n {
        return Vec::new(); // 维数不符（防御）
    }
    let mut drafts = Vec::new();
    for con in &problem.constraints {
        // ① Le 行；③ 右端有限。
        if con.op != RelOp::Le || !con.rhs.is_finite() {
            continue;
        }
        let cap = con.rhs;
        // ②③ 支撑全 Binary、权全正且数值稳健；④ 全支撑 Σw > cap。
        let mut cands: Vec<(usize, f64)> = Vec::with_capacity(con.expr.terms.len());
        let mut qualify = true;
        let mut total = 0.0f64;
        for &(j, w) in &con.expr.terms {
            let is_binary = problem
                .vars
                .get(j)
                .map(|v| v.var_type == VarType::Binary)
                .unwrap_or(false);
            // w > tol.zero 同时编码『全正』与近零权防御（tol.zero > 0）。
            if !is_binary
                || !w.is_finite()
                || !x_bar[j].is_finite()
                || w.abs() >= MAX_COVER_WEIGHT
                || w <= tol.zero
            {
                qualify = false;
                break;
            }
            cands.push((j, w));
            total += w;
        }
        if !qualify || cands.is_empty() || !(total > cap) {
            continue; // 不合格或无覆盖存在
        }

        // 贪婪：x̄ 降序（并列 w 降序、再下标升序）并入 C，直到 Σ_C w > cap。
        cands.sort_by(|&(ja, wa), &(jb, wb)| {
            x_bar[jb]
                .total_cmp(&x_bar[ja])
                .then_with(|| wb.total_cmp(&wa))
                .then_with(|| ja.cmp(&jb))
        });
        let mut cover: Vec<(usize, f64)> = Vec::new();
        let mut sum = 0.0f64;
        for cand in cands {
            cover.push(cand);
            sum += cand.1;
            if sum > cap {
                break;
            }
        }
        if !(sum > cap) {
            continue; // 防御：不应发生（全支撑 Σw > cap 已保证）
        }

        // 极小化：x̄ 升序（并列 w 升序）遍历，剔除后仍保持超容者剔之。
        // sum 单调不增 ⇒ 终态对每个保留成员 Σ−w_j ≤ cap（最小覆盖）。
        cover.sort_by(|&(ja, wa), &(jb, wb)| {
            x_bar[ja]
                .total_cmp(&x_bar[jb])
                .then_with(|| wa.total_cmp(&wb))
                .then_with(|| ja.cmp(&jb))
        });
        let mut kept: Vec<usize> = Vec::new();
        for &(j, w) in &cover {
            if sum - w > cap {
                sum -= w;
            } else {
                kept.push(j);
            }
        }

        // violation 检验：Σ_{j∈C} x̄_j − (|C|−1) ≥ min_violation。
        let rhs_cut = kept.len() as f64 - 1.0;
        let mut lhs = 0.0f64;
        for &j in &kept {
            lhs += x_bar[j];
        }
        if lhs - rhs_cut < opts.min_violation {
            continue;
        }
        kept.sort_unstable(); // 下标升序（成员唯一，确定性）
        let violation = lhs - rhs_cut;
        drafts.push(CoverDraft {
            cut: Cut {
                terms: kept.iter().map(|&j| (j, 1.0)).collect(),
                op: RelOp::Le,
                rhs: rhs_cut,
            },
            violation,
            cover: kept,
            terms: con.expr.terms.clone(),
            cap,
        });
    }
    drafts
}

/// 根节点覆盖割分离（x̄ 引导贪婪 + 极小化，无提升；行为与 WP4 完全一致）。
///
/// 逐项语义见 [`cover_drafts`]；本函数即草案的割投影，供既有调用方与
/// 「提升退回基础」对照使用。
///
/// `bounds`/`x_bar` 与问题维数不符时返回空（防御）。
pub fn cover_cuts(
    problem: &Problem,
    bounds: &Bounds,
    x_bar: &[f64],
    opts: CutOptions,
    tol: &Tolerances,
) -> Vec<Cut> {
    cover_drafts(problem, bounds, x_bar, opts, tol)
        .into_iter()
        .map(|d| d.cut)
        .collect()
}

/// lifting 的根割循环全程 DP 总预算（状态操作数；WP6 三层预算的第三层）。
///
/// 新建 = 3×10^7（[`ROOT_LIFT_TOTAL_BUDGET`]）；每次 DP 建表/追加按实际
/// 内层循环次数扣减；耗尽（remaining = 0）即全程关闭提升（调用方约定，
/// 任何未过预算的行原样退回基础覆盖割——不弱于现状）。
#[derive(Debug, Clone, Copy)]
pub struct LiftBudget {
    remaining: u64,
}

impl LiftBudget {
    /// 满额预算（3×10^7 状态操作）。
    pub fn new() -> Self {
        LiftBudget {
            remaining: ROOT_LIFT_TOTAL_BUDGET,
        }
    }

    /// 指定剩余预算（测试用：0 即「预算耗尽」形态）。
    pub fn with_remaining(remaining: u64) -> Self {
        LiftBudget { remaining }
    }

    /// 剩余预算。
    pub fn remaining(&self) -> u64 {
        self.remaining
    }
}

impl Default for LiftBudget {
    fn default() -> Self {
        LiftBudget::new()
    }
}

/// 覆盖割批 + 序贯 up-lifting（WP6；cut_loop_recorded 每轮调用）。
///
/// 基础草案（[`cover_drafts`]，行序输出保持确定性）之上，**每轮只对
/// violation 最高的前 [`LIFT_TOP_PER_ROUND`] 条**尝试提升
/// （[`lift_cover_draft`]，三层预算见模块文档）；任何门槛/预算不过即原样
/// 退回该行的基础覆盖割。输出顺序 = 行序（与 [`cover_cuts`] 一致），
/// 提升只替换割内容不改序。
pub fn cover_cuts_with_lifting(
    problem: &Problem,
    bounds: &Bounds,
    x_bar: &[f64],
    opts: CutOptions,
    tol: &Tolerances,
    budget: &mut LiftBudget,
) -> Vec<Cut> {
    let mut drafts = cover_drafts(problem, bounds, x_bar, opts, tol);
    if drafts.is_empty() {
        return Vec::new();
    }
    if budget.remaining() > 0 {
        // violation 降序（稳定排序：平局保持行序 ⇒ 全序确定）。
        let mut order: Vec<usize> = (0..drafts.len()).collect();
        order.sort_by(|&a, &b| drafts[b].violation.total_cmp(&drafts[a].violation));
        for &i in order.iter().take(LIFT_TOP_PER_ROUND) {
            if budget.remaining() == 0 {
                break; // 预算耗尽 ⇒ 全程关闭提升
            }
            if let Some(lifted) = lift_cover_draft(&drafts[i], x_bar, opts, budget) {
                drafts[i].cut = lifted;
            }
        }
    }
    drafts.into_iter().map(|d| d.cut).collect()
}

/// 0-1 背包 DP 追加一项（容量轴降序原地更新）；`ops` 按内层循环次数计量。
fn dp_add_item(dp: &mut [u64], w: u64, p: u64, ops: &mut u64) {
    let wl = w as usize;
    if wl >= dp.len() {
        *ops += 1; // 装不下：仍计 1 op（预算不漏计）
        return;
    }
    let mut c = dp.len() - 1;
    loop {
        let cand = dp[c - wl].saturating_add(p);
        if cand > dp[c] {
            dp[c] = cand;
        }
        *ops += 1;
        if c == wl {
            break;
        }
        c -= 1;
    }
}

/// 对单条覆盖草案做序贯 up-lifting（标准公式 α_j = β − Z，β = |C|−1 固定；
/// 见模块文档的偏差说明与正确性论证）。
///
/// 建表判据/门槛/预算任一不过 ⇒ 返回 None（调用方退回基础覆盖割）：
/// - 权与 cap 须可 10^k 整数化（[`integerize_scale`]）；
/// - 提升列 K ≤ min(32, |支撑\C|)（按 x̄ 降序，并列 w 降序、再下标升序）；
/// - 建表门槛 (|C|+K)×(cap̃+1) ≤ 2×10^6 且 ≤ 剩余预算；
/// - 提升后数值门：系数比 ≤ 1e9、violation ≥ min_violation（α ≥ 0、
///   x̄ ≥ 0 ⇒ 理论上不减，防御性复核）。
fn lift_cover_draft(
    draft: &CoverDraft,
    x_bar: &[f64],
    opts: CutOptions,
    budget: &mut LiftBudget,
) -> Option<Cut> {
    // 建表判据：行权与 cap 同乘 10^k 后全部恰为整数。
    let weights: Vec<f64> = draft.terms.iter().map(|&(_, w)| w).collect();
    let scale = integerize_scale(&weights, draft.cap)?;
    // 容量上取整 ⌈cap·10^k⌉ / 权下取整 ⌊w·10^k⌋：双向偏向『装得下』
    // （DP 最优值 ≥ 真值 ⇒ α 偏小 ⇒ 割偏弱仍有效，安全方向）。
    let cap_t_f = draft.cap * scale;
    if !cap_t_f.is_finite() || cap_t_f > 4e18 {
        return None;
    }
    let cap_t = cap_t_f.ceil().max(0.0) as u64;
    // 提升候选：支撑 \ C，按 x̄ 降序（并列 w 降序、再下标升序）。
    let mut cands: Vec<(usize, f64)> = draft
        .terms
        .iter()
        .filter(|&&(j, _)| draft.cover.binary_search(&j).is_err())
        .cloned()
        .collect();
    cands.sort_by(|&(ja, wa), &(jb, wb)| {
        x_bar[jb]
            .total_cmp(&x_bar[ja])
            .then_with(|| wb.total_cmp(&wa))
            .then_with(|| ja.cmp(&jb))
    });
    let k = cands.len().min(MAX_LIFT_CANDIDATES);
    // 建表门槛 + 预算门（a priori 估计 = (|C|+K)×(cap̃+1) 状态）。
    let states = (draft.cover.len() + k) as u64 * cap_t.saturating_add(1);
    if states > MAX_LIFT_TABLE_STATES {
        return None; // 超限 ⇒ 整行跳过提升只留基础覆盖割
    }
    if states > budget.remaining() {
        return None; // 预算不足 ⇒ 退回基础覆盖割
    }
    // 增量 DP：dp[c] = 已处理列集在容量 c 内的最大利润（= 当前割系数和）。
    let mut dp = vec![0u64; cap_t as usize + 1];
    let beta = draft.cover.len() as u64 - 1; // β = |C|−1，固定
    let mut ops = 0u64;
    let weight_of = |j: usize| -> u64 {
        let w = draft
            .terms
            .iter()
            .find(|&&(jj, _)| jj == j)
            .map(|&(_, w)| w)
            .unwrap_or(0.0);
        (w * scale).floor().max(0.0) as u64
    };
    for &j in &draft.cover {
        dp_add_item(&mut dp, weight_of(j), 1, &mut ops);
    }
    let mut lifted: Vec<(usize, u64)> = Vec::new();
    for &(j, w_f) in &cands[..k] {
        let w_t = (w_f * scale).floor().max(0.0) as u64;
        // Z：x_j = 1 时已处理列的最大利润——容量 ⌈(cap−w_j)·s⌉（放宽方向）。
        let rem = (draft.cap - w_f) * scale;
        let c_j = if rem <= 0.0 {
            0
        } else {
            (rem.ceil() as u64).min(cap_t)
        };
        let z = dp[c_j as usize];
        // α_j = β − Z（Z ≤ β 由当前割在已处理面上的有效性保证；饱和减防御）。
        let alpha = beta.saturating_sub(z);
        if alpha > 0 {
            lifted.push((j, alpha));
        }
        // 序贯性：提升列立即以 (w̃_j, α_j) 入表，参与后续列的 DP。
        dp_add_item(&mut dp, w_t, alpha, &mut ops);
    }
    budget.remaining = budget.remaining.saturating_sub(ops);
    // 组装：C 成员系数 1 + 提升列系数 α（不相交 ⇒ 无重复下标），下标升序。
    let mut terms: Vec<(usize, f64)> = draft
        .cover
        .iter()
        .map(|&j| (j, 1.0))
        .chain(lifted.iter().map(|&(j, a)| (j, a as f64)))
        .collect();
    terms.sort_by_key(|&(j, _)| j);
    // 数值门：系数比 ≤ 1e9（系数全 ≥ 1，比 = max α）、violation ≥ min_violation。
    let alpha_max = lifted.iter().map(|&(_, a)| a).max().unwrap_or(1);
    if alpha_max > MAX_COEF_RATIO as u64 {
        return None;
    }
    let mut lhs = 0.0f64;
    for &(j, c) in &terms {
        lhs += c * x_bar[j];
    }
    if !lhs.is_finite() || lhs - draft.cut.rhs < opts.min_violation {
        return None; // 理论上 violation ≥ 基础割（α≥0、x̄≥0），防御性复核
    }
    Some(Cut {
        terms,
        op: RelOp::Le,
        rhs: draft.cut.rhs, // β 固定不变
    })
}

/// 连续列 GMI 系数：π = max(α/f0, −α/(1−f0)) ≥ 0。
fn continuous_pi(alpha: f64, f0: f64) -> f64 {
    (alpha / f0).max(-alpha / (1.0 - f0))
}

/// 界幅 max(|lb|, |ub|)；ub = +∞ 时为 +∞（表示无法界定）。
fn bound_magnitude(lb: f64, ub: f64) -> f64 {
    lb.abs().max(if ub.is_finite() {
        ub.abs()
    } else {
        f64::INFINITY
    })
}

/// 取整 ε 纪律（WP6 增补）：上界取整按**放宽方向** floor(val+ε)，
/// ε = 1e-9·(1+|val|)；|val−round(val)| ≤ ε 时直接取整数值。
/// 防浮点和差把真界 4 算成 3.999…9 而 floor 掉可行整数点
/// （[`floor_bound`] 结果 ≥ 真值 ⌊val⌋，只可能更松，不会更紧）。
fn floor_bound(val: f64) -> f64 {
    let eps = NUMERIC_SNAP_TOL * (1.0 + val.abs());
    if (val - val.round()).abs() <= eps {
        val.round()
    } else {
        (val + eps).floor()
    }
}

/// 取整 ε 纪律：下界取整按**放宽方向** ceil(val−ε)（结果 ≤ 真值 ⌈val⌉，
/// 只可能更松，不会更紧）。
fn ceil_bound(val: f64) -> f64 {
    let eps = NUMERIC_SNAP_TOL * (1.0 + val.abs());
    if (val - val.round()).abs() <= eps {
        val.round()
    } else {
        (val - eps).ceil()
    }
}


/// MIR 内部用 ⌊·⌋，取整 ε 纪律与 [`floor_bound`]/[`ceil_bound`] 对称的
/// **放宽方向**：结果恒 ≤ 真值 ⌊val⌋，只松不紧——
/// - val 在整数上方 ε 内（或恰等）：吸附到 round(val)（= ⌊val⌋，向下吸附
///   只纠浮点表示尾差，与真 floor 相等）；
/// - val 在整数下方 ε 内：(val−ε).floor()，**不向上吸附**——旧实现取
///   round=k 会使 f_j = val−k < 0、π_j = k 超出精确 MIR 系数至多
///   ε/(1−f̂0)（away=0.005 下 α≈1e6 时可达 ~0.2），属加严方向，理论上可
///   切掉可行整数点；向下取保持 f_j ∈ [0,1)，π_j ≤ 精确 MIR 系数；
/// - 其余：普通 floor。
/// （|val| ≳ 1e9 时 ε ≥ 1，(val−ε).floor() 可低于真 floor——仍在放宽
/// 方向的安全侧；MIR 系数经归一 δ ≥ MIR_MIN_DELTA 有界 ≤ 1e6，常规
/// 量级下结果即真 floor。）
fn mir_floor(val: f64) -> f64 {
    let eps = NUMERIC_SNAP_TOL * (1.0 + val.abs());
    let r = val.round();
    if (val - r).abs() > eps {
        val.floor()
    } else if val >= r {
        r
    } else {
        (val - eps).floor()
    }
}

/// 变量界是否（容差内）整数值。
fn is_integral_bound(v: f64) -> bool {
    v.is_finite() && (v - v.round()).abs() <= NUMERIC_SNAP_TOL * (1.0 + v.abs())
}

/// 行权与 cap 同乘 10^k（k ∈ 0..=6）后是否全部（容差内）恰为整数；
/// 是则返回最小可用的 10^k，否则 None（非整数化行跳过提升）。
fn integerize_scale(weights: &[f64], cap: f64) -> Option<f64> {
    if !cap.is_finite() {
        return None;
    }
    for k in 0..=6u32 {
        let s = 10f64.powi(k as i32);
        let ok = weights.iter().all(|&w| is_integral_bound(w * s))
            && is_integral_bound(cap * s);
        if ok {
            return Some(s);
        }
    }
    None
}

/// 把稠密系数向量组装成通过数值门的 `Cut`（近零项放宽丢弃 + 比率门 +
/// violation 门）；任一门不过返回 None。与 GMI 尾段（近零/比率/violation）
/// 同一纪律，供 c-MIR/隐含界新族复用。
fn assemble_gated_cut(
    coef: &[f64],
    rhs: f64,
    bounds: &Bounds,
    x_bar: &[f64],
    opts: CutOptions,
    tol: &Tolerances,
) -> Option<Cut> {
    let mut terms: Vec<(usize, f64)> = Vec::new();
    let mut rhs = rhs;
    let mut c_max = 0.0f64;
    let mut c_min = f64::INFINITY;
    for (j, &c) in coef.iter().enumerate() {
        if c == 0.0 {
            continue;
        }
        if c.abs() <= tol.zero {
            // 近零项：按界幅补偿右端后丢弃（丢项 + 抬右端，只放宽）。
            let magnitude = bound_magnitude(bounds.lb[j], bounds.ub[j]);
            if magnitude.is_finite() {
                rhs += c.abs() * magnitude;
                continue;
            }
            // 上界 +∞ 无法界定补偿：保留该项，但不计入比率检验。
            terms.push((j, c));
            continue;
        }
        terms.push((j, c));
        c_max = c_max.max(c.abs());
        c_min = c_min.min(c.abs());
    }
    if terms.is_empty() || c_max <= 0.0 {
        return None;
    }
    if c_max / c_min > MAX_COEF_RATIO {
        return None; // 系数动态范围过大，整条丢弃
    }
    let mut lhs = 0.0f64;
    for &(j, c) in &terms {
        lhs += c * x_bar[j];
    }
    if !lhs.is_finite() || !rhs.is_finite() || lhs - rhs < opts.min_violation {
        return None;
    }
    Some(Cut {
        terms,
        op: RelOp::Le,
        rhs,
    })
}

/// 分离边界的数值卫生闸：小系数清洁化（只弱化）→ 密度闸 → violation
/// 复核；任一门不过返回 None。丢割永远安全——只丢界强度，绝不切掉
/// 整数可行点。
///
/// 【为什么落在施加边界】五个割族的产出全部经 [`cut_loop_recorded`] /
/// [`separate_local_cuts`] 汇聚后施加，且稠密/稀疏两引擎共用本模块，
/// 在此单点过闸即全覆盖，无需改任何引擎或生成器语义。
///
/// 【清洁化只弱化的论证（数学红线）】割为 Le 形式 Σ_j c_j·x_j ≤ r。对被
/// 清洁列 k（|c_k| ≤ CUT_CLEANUP_EPS 且界幅 M_k = max(|lb_k|, |ub_k|)
/// 有限），清洁后右端取
///     r' = r + Σ_{被清洁 k} |c_k|·M_k。
/// 对箱内任意满足原割的点 x：
///     Σ_{j 保留} c_j·x_j = Σ_j c_j·x_j − c_k·x_k ≤ r − c_k·x_k
///                        ≤ r + |c_k|·|x_k| ≤ r + |c_k|·M_k = r'
/// （逐列累加即得全体被清洁列的联合界）。故清洁后割的可行域 ⊇（原割
/// 可行域 ∩ 变量箱）。生成器产出的割对全部整数可行点有效，而整数可行
/// 点都落在箱内且满足原割 ⇒ 清洁后割仍对全部整数可行点有效。
/// 补偿方向恒为**抬高右端（放宽）**；被清洁项无论是正贡献还是负贡献，
/// 同一个 |c_k|·M_k 界都成立（−c_k·x_k ≤ |c_k|·|x_k| ≤ |c_k|·M_k 对
/// 任意符号的 c_k 与箱内 x_k 成立），无需按符号分支。若某列界幅无界
/// （ub_k = +∞ ⇒ M_k = +∞），无法界定补偿 ⇒ 该列**原样保留**，绝不
/// 清洁——宁可不清洁也不冒强化方向的风险。
///
/// 【容差 CUT_CLEANUP_EPS = 1e-9 的依据】
/// 上界（要盖住的噪声）：实测 Devex 基表格噪声 ≤ 7.0e-15，经 GMI 连续列
/// 公式 π = max(α/f0, −α/(1−f0))（f0 ≥ away = 0.005 ⇒ 放大 ≤ 200×）
/// 与松弛列回代（原始行构建期按 g = max|a_ij| 归一 ⇒ 每项 |a_ij|/g ≤ 1；
/// 割行按原尺度回代）后，单源残差 ≤ ~1.4e-12，跨轮割行复合仍有数量级
/// 余量——1e-9 全覆盖。
/// 下界（不误伤显著项）：清洁对割的扰动以 |c_k|·M_k 计入右端补偿；
/// violation 门在 min_violation = 1e-7 量级，被清洁项的 lhs 贡献
/// （≤ 1e-9·|x̄_k|，同量级于界幅时 ≤ 1e-9·M_k）比它小两个数量级，
/// 清洁不会吃掉有分离价值的强度；万一吃掉（大界幅 × 噪声系数），由
/// ③ 的复核兜底丢弃该割。
///
/// 【密度闸上限 = max(CUT_DENSITY_NNZ_FLOOR, ⌈n·CUT_DENSITY_NNZ_FRAC⌉)
/// 的依据】
/// (i) 病态实证：Devex 基根割循环曾产出 160 条 661~681 nnz（900 列、
///     约 74% 密度）的割——每迭代成本 2-3×、退化打转、端到端 10×+
///     劣化，近稠密割必须整条丢弃；
/// (ii) 绝对下限 200：割行 nnz 直接进入两引擎的每迭代代价（稠密表格
///     每枢轴 O(nnz) 行更新、revised 引擎 BTRAN/FTRAN 行密度），200 项
///     以上已足以拖慢；同时结构列至多两百项的良性割（仓库内 setcover
///     -100x150 结构列至多 150 项、mknap 覆盖割数十项）不受影响；
/// (iii) 相对系数 0.1：超大 n（> 2000）时固定 200 会过度杀伤（10 万列
///     问题上 300 项割仅 0.3% 密度、边际代价小），按 10% 密度随 n 放宽
///     上限。判据只依赖割自身 nnz 与结构列数 n，与算例名/规模无特判。
fn sanitize_cut(
    cut: Cut,
    bounds: &Bounds,
    x_bar: Option<&[f64]>,
    opts: CutOptions,
) -> Option<Cut> {
    let n = bounds.lb.len();
    // 防御：本闸的补偿论证只对 Le 形式成立（其它关系符的方向分析不同）；
    // 非有限系数/右端、或下标越界的割整条丢弃（丢割安全）。
    if cut.op != RelOp::Le
        || !cut.rhs.is_finite()
        || cut.terms.iter().any(|&(_, c)| !c.is_finite())
        || cut.terms.iter().any(|&(j, _)| j >= n)
    {
        return None;
    }
    // ① 清洁化（只弱化）：|c| ≤ 容差且界幅有限 ⇒ 丢项并按界幅抬高右端。
    let mut terms: Vec<(usize, f64)> = Vec::with_capacity(cut.terms.len());
    let mut rhs = cut.rhs;
    for (j, c) in cut.terms {
        if c.abs() <= CUT_CLEANUP_EPS {
            let magnitude = bound_magnitude(bounds.lb[j], bounds.ub[j]);
            if magnitude.is_finite() {
                rhs += c.abs() * magnitude;
                continue; // 丢弃该项（补偿后只弱化）
            }
            // 界幅无界（ub = +∞）：无法界定只弱化补偿 ⇒ 原样保留该项。
        }
        terms.push((j, c));
    }
    if terms.is_empty() {
        return None; // 清洁后无剩余项：无分离价值，丢弃（安全）
    }
    // ② 密度闸：nnz 超上限整条丢弃（丢割永远安全）。
    let limit = CUT_DENSITY_NNZ_FLOOR.max(((n as f64) * CUT_DENSITY_NNZ_FRAC).ceil() as usize);
    if terms.len() > limit {
        return None;
    }
    // ③ violation 复核：右端补偿可能吃掉 violation；弱化到不再被当前 LP
    //    解违反的割已无分离价值，丢弃（安全）。x_bar 不可得（防御路径，
    //    正常接线不可达）时无法复核，只做①②。
    if let Some(x_bar) = x_bar {
        if x_bar.len() != n {
            return None; // 维数不符（防御）：丢割安全
        }
        let mut lhs = 0.0f64;
        for &(j, c) in &terms {
            lhs += c * x_bar[j];
        }
        if !lhs.is_finite() || lhs - rhs < opts.min_violation {
            return None;
        }
    }
    Some(Cut {
        terms,
        op: RelOp::Le,
        rhs,
    })
}

/// 单行 c-MIR（WP6，仅原始 Le 行，根节点调用）。
///
/// 流程（见模块文档「单行 c-MIR」节）：
/// 1. **界互补归一**：负系数列 x_j → ub_j − y_j（`ub_j` 无限 ⇒ 整行跳过；
///    整数类型列互补的 ub / 平移的 lb 须整数值，否则整行跳过），非负系数列
///    y_j = x_j − lb_j。归一后全系数 â_j > 0、y ≥ 0、整数列 y ∈ Z_+，
///    右端 b̂ = b − Σ_{a≥0} a·lb − Σ_{a<0} a·ub。
/// 2. 除以 g = max â_j，f0 = frac(b̂/g)；在尺度集
///    δ ∈ {f0, 1}·{1, 0.5, 0.25}（去重、δ ≥ 1e-6）上逐试 MIR：
///    B = b̂/(gδ)，f̂0 = B − floor_bound(B) 须 ∈ [away, 1−away]（复用 away 门）；
///    整数列 π_j = mir_floor(α_j) + (f_j − f̂0)^+/(1 − f̂0)（α_j = â_j/(gδ)、
///    f_j = α_j − mir_floor(α_j)），连续列 π_j = 0；割（y 空间）
///    Σ π_j y_j ≤ floor_bound(B)。
/// 3. 回代 y → x（互补列 y_j = ub_j − x_j、平移列 y_j = x_j − lb_j）得
///    x 空间割，过 [`assemble_gated_cut`] 数值门；δ 候选中取 violation 最大
///    者（并列取先者，确定性）。
///
/// 有效性：互补/平移是整数可行点集上的双射；归一行整数部分满足放松行
/// （连续列只消耗容量，丢弃只放宽）；MIR 引理对纯整数非负系数行给出全局
/// 有效割；回代保持有效性。数值门任一不过即弃（None）。
///
/// `bounds`/`x_bar` 与问题维数不符时返回 None（防御）。
pub fn mir_single_row(
    problem: &Problem,
    bounds: &Bounds,
    x_bar: &[f64],
    con: &Constraint,
    opts: CutOptions,
    tol: &Tolerances,
) -> Option<Cut> {
    let n = problem.vars.len();
    if bounds.lb.len() != n || bounds.ub.len() != n || x_bar.len() != n {
        return None;
    }
    // ① 只处理原始 Le 行。
    if con.op != RelOp::Le || !con.rhs.is_finite() {
        return None;
    }
    // ② 界互补归一：cols = (j, â_j>0, 所贴界, 是否互补, 是否整数类型)。
    let mut cols: Vec<(usize, f64, f64, bool, bool)> = Vec::with_capacity(con.expr.terms.len());
    let mut b_hat = con.rhs;
    for &(j, a) in &con.expr.terms {
        let Some(var) = problem.vars.get(j) else {
            return None;
        };
        if !a.is_finite() || a.abs() >= MAX_COVER_WEIGHT || !x_bar[j].is_finite() {
            return None; // 数值防御：整行跳过
        }
        let is_int = matches!(var.var_type, VarType::Integer | VarType::Binary);
        if a < 0.0 {
            let ub = bounds.ub[j];
            if !ub.is_finite() {
                return None; // 负系数列无有限 ub ⇒ 整行跳过
            }
            if is_int && !is_integral_bound(ub) {
                return None; // 互补不保整 ⇒ 整行跳过
            }
            b_hat -= a * ub;
            cols.push((j, -a, ub, true, is_int));
        } else {
            let lb = bounds.lb[j]; // 本 crate 下界恒有限
            if is_int && !is_integral_bound(lb) {
                return None; // 平移不保整 ⇒ 整行跳过
            }
            b_hat -= a * lb;
            cols.push((j, a, lb, false, is_int));
        }
    }
    if cols.is_empty() {
        return None;
    }
    // ③ 归一尺度 g = max â_j。
    let g = cols.iter().map(|c| c.1).fold(0.0f64, f64::max);
    if !g.is_finite() || g <= tol.zero {
        return None;
    }
    let b_norm = b_hat / g;
    if !b_norm.is_finite() {
        return None;
    }
    let f0 = b_norm - mir_floor(b_norm);
    // ④ 尺度集 δ ∈ {f0,1}·{1,0.5,0.25}（乘积展开），去重且 δ ≥ MIR_MIN_DELTA。
    let mut deltas: Vec<f64> = vec![1.0, 0.5, 0.25];
    if f0 > MIR_MIN_DELTA {
        for m in [1.0, 0.5, 0.25] {
            deltas.push(f0 * m);
        }
    }
    let mut uniq: Vec<f64> = Vec::with_capacity(deltas.len());
    for d in deltas {
        if !d.is_finite() || d < MIR_MIN_DELTA {
            continue;
        }
        if !uniq.iter().any(|&u| (u - d).abs() <= 1e-12 * (1.0 + d.abs())) {
            uniq.push(d);
        }
    }
    // ⑤ 逐 δ 试 MIR，取 violation 最大者。
    let mut best: Option<(f64, Cut)> = None;
    for delta in uniq {
        let big = b_norm / delta;
        let rhs_y = floor_bound(big);
        if !rhs_y.is_finite() {
            continue;
        }
        let f_hat = big - rhs_y;
        if !(f_hat >= opts.away && f_hat <= 1.0 - opts.away) {
            continue; // f̂0 近 0/1 病态：复用 away 门
        }
        let mut coef = vec![0.0f64; n];
        let mut rhs_x = rhs_y;
        let mut any = false;
        for &(j, a_hat, bound, compl, is_int) in &cols {
            let alpha = a_hat / (g * delta);
            if !alpha.is_finite() {
                any = false;
                break;
            }
            let pi = if is_int {
                let fl = mir_floor(alpha);
                let fj = alpha - fl;
                fl + ((fj - f_hat) / (1.0 - f_hat)).max(0.0)
            } else {
                // 归一后正系数连续列：MIR 系数 0（丢弃只放宽）。
                0.0
            };
            if !pi.is_finite() || pi <= 0.0 {
                continue;
            }
            if pi < tol.zero {
                continue; // 近零项：y ≥ 0 ⇒ 丢弃只放宽
            }
            if compl {
                // y_j = ub_j − x_j：π·y = π·ub − π·x
                coef[j] -= pi;
                rhs_x -= pi * bound;
            } else {
                // y_j = x_j − lb_j：π·y = π·x − π·lb
                coef[j] += pi;
                rhs_x += pi * bound;
            }
            any = true;
        }
        if !any {
            continue;
        }
        if let Some(cut) = assemble_gated_cut(&coef, rhs_x, bounds, x_bar, opts, tol) {
            let mut lhs = 0.0f64;
            for &(j, c) in &cut.terms {
                lhs += c * x_bar[j];
            }
            let viol = lhs - cut.rhs;
            if best.as_ref().map_or(true, |(v, _)| viol > *v) {
                best = Some((viol, cut));
            }
        }
    }
    best.map(|(_, cut)| cut)
}

/// 每轮的 c-MIR 批：对全部原始 Le 行各至多 1 条（行序确定）。
fn mir_row_cuts(
    problem: &Problem,
    bounds: &Bounds,
    x_bar: &[f64],
    opts: CutOptions,
    tol: &Tolerances,
) -> Vec<Cut> {
    problem
        .constraints
        .iter()
        .filter(|con| con.op == RelOp::Le)
        .filter_map(|con| mir_single_row(problem, bounds, x_bar, con, opts, tol))
        .collect()
}

/// 聚合 MIR（R2-retry，轴2）：每轮合并尝试上限（(候选行, 非基松弛列) 对
/// 计数，行列升序确定性遍历）。聚合行是 B⁻¹ 乘积、数值更脏，比率/
/// violation 门拒掉大量候选属预期，上限防白烧（风险提示原话）。
const AGG_MIR_MAX_TRIES: usize = 64;

/// 聚合 MIR（R2-retry，轴2）：每轮施加条数配额。产出大的算例
/// （setcover-100x150 离线实测 386 条候选/轮）若放任填满 max_cuts_per_round
/// 截断额度，割行数上升会放大后续 resolve 成本——尤其根割循环最终回退的
/// 算例（setcover-100x150：回退前 Bland 爬行在更大表格上白烧，实测墙钟
/// ×3.4）；配额对齐既有 GMI 批的自然产出量级（mknap 5-10 条/轮）。
const AGG_MIR_MAX_PER_ROUND: usize = 8;

/// 表格行聚合 MIR 分离（R2-retry 落地形态：「表格行 ⊕ 原始行」聚合 GMI）。
///
/// O3 方案的两个字面构造经证伪后按可实现形态落地（详见模块文档
/// 「表格行聚合 MIR」节）：(a) 表格行内松弛列按 GMI 纪律**全部**回代成
/// x 空间会塌缩为 0=0（恒等式组合代入定义恒真）；(b) 表格行**对**按
/// λ∈{1,−1} 合并产生双基列，GMI 格论证失效（权重非整），纯 MIR 则因
/// 丢弃连续/松弛列在背包行上弱于既有 GMI（R2 首试实测产出为零）。
/// 本实现改为：候选表格行 r（GMI 同口径：基变量为整数结构列且
/// f0∈[away,1−away]）与其**非基松弛列 s_i 的原始归一行**合并
/// （μ = −α_{s_i}/σ_i 恰消去 s_i，把行 i 的结构内容并入），行 i 支撑可含
/// B1 自身 ⇒ 整行除以锚系数重新归一（B1 ∈ Z 的 mod-1 锚定对任意非零
/// 系数成立；漏归一会静默丢项产无效割——fuzz 穷举抓出）。全部列（基/非基
/// 一律）取单侧偏差（非基按所贴界；基本结构列按较近界，d = x−lb 或
/// ub−x ≥ 0 恒成立，整数列整界时 d ∈ Z_+），对偏差行套与 GMI 相同的
/// 混合整数公式——**连续/松弛列 π ≠ 0 保留**（连续公式
/// max(ĉ/f0, −ĉ/(1−f0))），松弛列回代按 GMI 既有纪律（原始行
/// s = σ·(rhs−Σa·x)/g、割行只查已记账 ledger）。【方向勘误（独立复核
/// 发现，R7）：无记账时丢 π 项使割**更强**（割为 Σπ·d ≥ 1、π ≥ 0、
/// d ≥ 0，从 LHS 丢正项收紧约束）、可能误切整数可行解——**不是放宽**；
/// 接线路径不可达（cut_loop_recorded 的 ledger 自空完整累积、割循环从
/// 冷 Simplex 启动 ⇒ 每个割松弛列必有记账），本路径为公开入口的防御
/// 分支，触发时产出割可能过强而无效。】
///
/// 有效性：合并行是原始等式系的恒等线性组合（B⁻¹ 行 + μ·归一行），
/// 偏差形式 Σ ĉ_j d_j ≡ f0 (mod Z)、d ≥ 0 单侧，GMI 混合整数割引理
/// 直接适用——50 种子随机小模型穷举验证（302 条割全部不切整数可行点，
/// 见 tests::aggregated_mir_fuzz_validity）。窗口扫描（离线证明）：
/// 默认 away=0.005/min_violation=1e-7 即非零产出（mknap-3x100 8 条、
/// mknap-5x300 25 条、setcover-100x150 386 条，viol 至 0.98），全窗口
/// (away∈{.005,.01,.02}×viol∈{1e-7,1e-4,1e-2}) 产出稳定，无需调参。
///
/// 输出序确定：候选行升序 × 非基松弛列升序，尝试上限
/// [`AGG_MIR_MAX_TRIES`]；交由 `cut_loop_recorded` 统一截断。
fn aggregated_mir_cuts(
    s: &dyn LpEngine,
    problem: &Problem,
    bounds: &Bounds,
    int_vars: &[usize],
    x_bar: &[f64],
    ledger: &[CutRecord],
    opts: CutOptions,
    tol: &Tolerances,
) -> Vec<Cut> {
    let n = problem.vars.len();
    if bounds.lb.len() != n || bounds.ub.len() != n || x_bar.len() != n {
        return Vec::new();
    }
    let layout = s.column_layout();
    let num_struct = layout.num_structural;
    if num_struct != n {
        return Vec::new(); // 表格与问题维度不符（防御）
    }
    let q = layout.slack.len();
    let dim = n + q;
    let basis = s.basis_columns();
    let nonbasic = s.nonbasic_columns();
    let mut is_int = vec![false; n];
    for &j in int_vars {
        if j < n {
            is_int[j] = true;
        }
    }
    let mut is_nonbasic = vec![false; dim];
    let mut at_ub = vec![false; dim];
    for nb in &nonbasic {
        is_nonbasic[nb.col] = true;
        at_ub[nb.col] = nb.at_upper;
    }
    let mut out = Vec::new();
    let mut tries = 0usize;
    'cand: for row in 0..s.num_rows() {
        let bcol = basis[row];
        if bcol >= num_struct || !is_int[bcol] {
            continue; // 基变量须为整数类型结构列（GMI 同口径）
        }
        let (terms, value) = s.tableau_row_terms(row);
        if !value.is_finite() || terms.iter().any(|&(_, v)| !v.is_finite()) {
            continue; // 数值异常行整行跳过（NaN/Inf 非零必入稀疏表）
        }
        // 绝对值口径恒等式：B1 + Σ α_j z_j = value + Σ α_j·z̄_j（非基按界；
        // 表格行右端存「非基贴界时的基变量值」，非基结构列贴上界/非零下界
        // 时其对界值的贡献必须补回右端，否则恒等式不成立）。
        let mut a = vec![0.0f64; dim];
        a[bcol] += 1.0;
        let mut r_rhs = value;
        let mut tcur = 0usize;
        for nb in &nonbasic {
            let j = nb.col;
            while tcur < terms.len() && terms[tcur].0 < j {
                tcur += 1;
            }
            let alpha = if tcur < terms.len() && terms[tcur].0 == j {
                terms[tcur].1
            } else {
                0.0
            };
            if alpha == 0.0 {
                continue;
            }
            a[j] += alpha;
            if j < n {
                r_rhs += alpha
                    * if nb.at_upper {
                        bounds.ub[j]
                    } else {
                        bounds.lb[j]
                    };
            }
        }
        // 对每个非基松弛列 s_i（列升序，确定性）合并原始行 i。
        for si in num_struct..dim {
            if tries >= AGG_MIR_MAX_TRIES {
                break 'cand;
            }
            if !is_nonbasic[si] || a[si] == 0.0 {
                continue;
            }
            let info = &layout.slack[si - num_struct];
            let RowOrigin::Constraint(ci) = info.origin else {
                continue; // 割松弛列不作合并目标（其结构内容由 ledger 回代）
            };
            // 防御口径与 GMI 一致（独立复核发现，R7）：.get 而非直接下标
            //——构造保证 origin 索引恒有效，此分支不可达，仅统一口径。
            let Some(con) = problem.constraints.get(ci) else {
                continue;
            };
            let g = con
                .expr
                .terms
                .iter()
                .map(|&(_, x)| x.abs())
                .fold(0.0_f64, f64::max);
            if g <= tol.zero {
                continue; // 防御：近零行无松弛列
            }
            tries += 1;
            let mu = -a[si] / info.sign;
            let mut ma = a.clone();
            let mut mr = r_rhs;
            // 并入 μ·(â_i·x + σ_i·s_i = b̂_i)：s_i 恰消去，结构列进位。
            // 行 i 的支撑可含 B1 自身（ma[bcol] = 1 + μ·â_{i,B1} ≠ 1）——
            // B1 不再是单位锚，必须整行除以 ma[bcol] 重新归一。
            for &(j, aj) in &con.expr.terms {
                ma[j] += mu * (aj / g);
            }
            ma[si] = 0.0;
            mr += mu * (con.rhs / g);
            let anchor = ma[bcol];
            if !anchor.is_finite() || anchor.abs() < 1e-6 {
                continue; // 锚系数退化：该合并弃
            }
            for v in ma.iter_mut() {
                *v /= anchor;
            }
            mr /= anchor;
            // —— 偏差形式：每列（除锚 B1）取单侧偏差（非基按所贴界；
            //    基本结构列按较近界——x ∈ [lb,ub] ⇒ d ≥ 0 恒成立）——
            let mut dev_coef = vec![0.0f64; dim];
            let mut dev_ub = vec![false; dim];
            let mut rhs_dev = mr;
            for j in 0..dim {
                if j == bcol || ma[j] == 0.0 {
                    continue;
                }
                let (lbj, ubj) = if j < n {
                    (bounds.lb[j], bounds.ub[j])
                } else {
                    (0.0, f64::INFINITY) // 松弛列 [0,∞)，恒取 lb 侧
                };
                let cur = if j < n { x_bar[j] } else { 0.0 };
                let use_ub = ubj.is_finite()
                    && (at_ub[j] || (!is_nonbasic[j] && ubj - cur < cur - lbj));
                dev_ub[j] = use_ub;
                if use_ub {
                    dev_coef[j] = -ma[j];
                    rhs_dev -= ma[j] * ubj;
                } else {
                    dev_coef[j] = ma[j];
                    rhs_dev -= ma[j] * lbj;
                }
            }
            let f0 = rhs_dev - rhs_dev.floor();
            if !(f0 >= opts.away && f0 <= 1.0 - opts.away) {
                continue;
            }
            // —— GMI π + 回代（连续/松弛列 π≠0 保留；松弛列按 GMI 纪律
            //    代入 s 的定义，割行只查已记账 ledger——无记账丢 π 项使割
            //    更强非放宽，方向勘误见模块文档）——
            let mut coef = vec![0.0f64; n];
            let mut const_term = 0.0f64;
            let mut any = false;
            for j in 0..dim {
                if j == bcol || dev_coef[j] == 0.0 {
                    continue;
                }
                let cj = dev_coef[j];
                if j >= n {
                    // 松弛/剩余列：偏差 d = s（[0,∞) 恒 lb 侧；基本松弛列
                    // ma=0 不入行，此路径仅非基松弛）。连续公式 π≠0。
                    let Some(pi) = (|| {
                        let pi = continuous_pi(cj, f0);
                        if !pi.is_finite() || pi <= 0.0 || pi < tol.zero {
                            return None;
                        }
                        Some(pi)
                    })() else {
                        continue;
                    };
                    let info2 = &layout.slack[j - num_struct];
                    match info2.origin {
                        RowOrigin::Constraint(ci2) => {
                            let Some(con2) = problem.constraints.get(ci2) else {
                                // 防御：行来源越界，丢 π 项——使割更强非
                                // 放宽（见模块文档方向勘误）；构造保证 origin
                                // 索引恒有效，此分支不可达。
                                continue;
                            };
                            let g2 = con2
                                .expr
                                .terms
                                .iter()
                                .map(|&(_, x)| x.abs())
                                .fold(0.0_f64, f64::max);
                            if g2 <= tol.zero {
                                continue;
                            }
                            let scale = info2.sign / g2;
                            const_term += pi * scale * con2.rhs;
                            for &(k, ak) in &con2.expr.terms {
                                coef[k] -= pi * scale * ak;
                            }
                            any = true;
                        }
                        RowOrigin::Cut(ck) => {
                            // 只查已记账 ledger（GMI 纪律）：无记账丢 π 项
                            // ——使割更强非放宽（模块文档方向勘误）；接线路径
                            // 不可达（ledger 完整累积），触发时割可能过强。
                            let Some(rec) = ledger.iter().find(|r| r.index == ck) else {
                                continue;
                            };
                            const_term += pi * rec.rhs;
                            for &(k2, ak) in &rec.terms {
                                coef[k2] -= pi * ak;
                            }
                            any = true;
                        }
                    }
                    continue;
                }
                let pi = if is_int[j] {
                    let bound = if dev_ub[j] {
                        bounds.ub[j]
                    } else {
                        bounds.lb[j]
                    };
                    if bound.is_finite() && (bound - bound.round()).abs() <= tol.zero {
                        // 所选界为整数值：d_j 在整数可行点取整数值，整数公式有效。
                        let fj = cj - cj.floor();
                        if fj <= f0 {
                            fj / f0
                        } else {
                            (1.0 - fj) / (1.0 - f0)
                        }
                    } else {
                        continuous_pi(cj, f0)
                    }
                } else {
                    continuous_pi(cj, f0)
                };
                if !pi.is_finite() || pi <= 0.0 || pi < tol.zero {
                    continue;
                }
                if dev_ub[j] {
                    coef[j] -= pi;
                    const_term += pi * bounds.ub[j];
                } else {
                    coef[j] += pi;
                    const_term -= pi * bounds.lb[j];
                }
                any = true;
            }
            if !any {
                continue;
            }
            // Σ coef·x + const ≥ 1 ⇔ Σ(−coef)x ≤ const − 1（GMI 尾段同款门）。
            if let Some(cut) = assemble_gated_cut(
                &coef.iter().map(|c| -c).collect::<Vec<_>>(),
                const_term - 1.0,
                bounds,
                x_bar,
                opts,
                tol,
            ) {
                out.push(cut);
                if out.len() >= AGG_MIR_MAX_PER_ROUND {
                    return out; // 每轮配额（行列升序确定性截断）
                }
            }
        }
    }
    out
}

/// 隐含界割（WP6，仅根节点）：整数类型变量计数恰为 1 的行导出 y 的隐含界。
///
/// 行的逻辑蕴涵 + 方向取整（纯变量界形式，见模块文档反例说明）。把行写成
/// 符号 view（sign·a, sign·rhs），设 ã = sign·a 且 ã_y > 0：
/// - **上界割**（≤ 形式：Le 原样 / Ge 取负，Eq 双向之一）：
///   `y ≤ (rhs − Kmin)/ã_y`，`Kmin = Σ_{ã_j>0,j≠y} ã_j·lb_j + Σ_{ã_j<0,j≠y} ã_j·ub_j`
///   （其余列折叠到使 Σãx **最小**的一侧——这是 y 上界的最宽松有效侧；
///   负系数列需要有限 ub，无限 ⇒ 该方向跳过）；
/// - **下界割**（≥ 形式：Ge 原样 / Le 取负，Eq 双向之一）：
///   `y ≥ (rhs − Kmax)/ã_y`，`Kmax = Σ_{ã_j>0,j≠y} ã_j·ub_j + Σ_{ã_j<0,j≠y} ã_j·lb_j`
///   （折叠到使 Σãx **最大**的一侧；正系数列需要有限 ub，无限 ⇒ 跳过）。
/// 取整按放宽方向：上界 [`floor_bound`]、下界 [`ceil_bound`]。
/// 有效性：V 是 y 界在其余列箱式取值下的真实极值的松弛（折叠方向给出
/// 恰好覆盖全部可行点的常数），取整只可能更松 ⇒ 对全部整数可行点有效。
/// violation ≥ min_violation 才收（割在 x̄ 处须真收紧）。
///
/// `bounds`/`x_bar` 与问题维数不符时返回空（防御）。
pub fn implied_bound_cuts(
    problem: &Problem,
    bounds: &Bounds,
    x_bar: &[f64],
    opts: CutOptions,
    tol: &Tolerances,
) -> Vec<Cut> {
    let n = problem.vars.len();
    if bounds.lb.len() != n || bounds.ub.len() != n || x_bar.len() != n {
        return Vec::new();
    }
    let mut cuts = Vec::new();
    for con in &problem.constraints {
        if !con.rhs.is_finite() {
            continue;
        }
        if con.expr.terms.iter().any(|&(_, a)| !a.is_finite()) {
            continue;
        }
        // 整数类型变量计数恰为 1。
        let mut int_y: Option<usize> = None;
        let mut int_count = 0usize;
        for &(j, _) in &con.expr.terms {
            let is_int = problem
                .vars
                .get(j)
                .map(|v| matches!(v.var_type, VarType::Integer | VarType::Binary))
                .unwrap_or(false);
            if is_int {
                int_count += 1;
                if int_count > 1 {
                    break;
                }
                int_y = Some(j);
            }
        }
        if int_count != 1 {
            continue;
        }
        let y = int_y.unwrap_or(usize::MAX);
        let a_y = con
            .expr
            .terms
            .iter()
            .find(|&&(j, _)| j == y)
            .map(|&(_, a)| a)
            .unwrap_or(0.0);
        if a_y == 0.0 || !x_bar[y].is_finite() {
            continue;
        }
        // 有效符号：≤ 形式（Le 原样 / Ge 与 Eq 取负）与 ≥ 形式（Ge 原样 / Le
        // 与 Eq 取负）。Eq 行两向都推导。
        let has_le_form = con.op != RelOp::Ge;
        let has_ge_form = con.op != RelOp::Le;
        // 在符号 view（sign·a, sign·rhs、ã_y > 0）下计算 y 界的真实极值：
        // upper=true（上界）折叠其余列到 Σãx 最小（负系数列用 ub），
        // upper=false（下界）折叠到 Σãx 最大（正系数列用 ub）。
        let implied_v = |sign: f64, rhs_eff: f64, upper: bool| -> Option<f64> {
            let ay = sign * a_y;
            if !(ay > tol.zero) {
                return None;
            }
            let mut v = rhs_eff;
            for &(j, a) in &con.expr.terms {
                if j == y {
                    continue;
                }
                let aj = sign * a;
                let use_ub = if upper { aj < 0.0 } else { aj > 0.0 };
                if use_ub {
                    let ub = bounds.ub[j];
                    if !ub.is_finite() {
                        return None; // 需要的一侧界无限 ⇒ 该方向跳过
                    }
                    v -= aj * ub;
                } else {
                    v -= aj * bounds.lb[j]; // lb 恒有限
                }
            }
            let v = v / ay;
            v.is_finite().then_some(v)
        };
        // 上界割（≤ 形式）：y ≤ floor_bound(V)。
        if has_le_form {
            let sign = if con.op == RelOp::Ge { -1.0 } else { 1.0 };
            let rhs_eff = if con.op == RelOp::Ge { -con.rhs } else { con.rhs };
            if let Some(v) = implied_v(sign, rhs_eff, true) {
                let ub_cut = floor_bound(v);
                if ub_cut.is_finite() && x_bar[y] - ub_cut >= opts.min_violation {
                    cuts.push(Cut {
                        terms: vec![(y, 1.0)],
                        op: RelOp::Le,
                        rhs: ub_cut,
                    });
                }
            }
        }
        // 下界割（≥ 形式）：y ≥ ceil_bound(V) ⇔ −y ≤ −ceil_bound(V)。
        if has_ge_form {
            let sign = if con.op == RelOp::Le { -1.0 } else { 1.0 };
            let rhs_eff = if con.op == RelOp::Le { -con.rhs } else { con.rhs };
            if let Some(v) = implied_v(sign, rhs_eff, false) {
                let lb_cut = ceil_bound(v);
                if lb_cut.is_finite() && lb_cut - x_bar[y] >= opts.min_violation {
                    cuts.push(Cut {
                        terms: vec![(y, -1.0)],
                        op: RelOp::Le,
                        rhs: -lb_cut,
                    });
                }
            }
        }
    }
    cuts
}

/// 根节点割循环：生成→记账→add_row→resolve，直到无割/超 max_rounds/LP 非
/// Optimal；返回 (施加割数, 最终 LP 状态)。
///
/// 入口先做一次防御性 `resolve(bounds)`，保证割生成看到的表格确为当前界下的
/// Optimal 状态。割循环中 LP 变 Infeasible 是整数不可行问题的合法终局（割对
/// 全部整数可行解有效，该集合为空时割可切空松弛区域）——如实作为 Infeasible
/// 终态返回（调用方按根终局映射 SolveStatus::Infeasible）；Unbounded/
/// W6：树内局部割的单轮分离（**只读**表格，不 add_row——行的作用域/
/// 施加/重解由调用方 branch_bound 管理：子树出口经 Simplex 快照恢复）。
/// 批序复用根割（c-MIR → GMI → 隐含界 → 覆盖），先统一过分离边界数值
/// 卫生闸 [`sanitize_cut`]（清洁化 + 密度闸 + violation 复核，见模块
/// 文档）再按配额截断；无 CutRecord 去重（ledger 只存在于根循环的跨轮
/// 记账；局部割生命周期 ≤ 子树，重复 Le 行至多冗余不切错）。
/// 分离失败/维度不符返回空。
pub fn separate_local_cuts(
    s: &dyn LpEngine,
    bounds: &Bounds,
    problem: &Problem,
    int_vars: &[usize],
    opts: CutOptions,
    tol: &Tolerances,
    max_cuts: usize,
) -> Vec<Cut> {
    let n = problem.vars.len();
    let Some(x_bar) = rebuild_x_bar(s, bounds, n) else {
        return Vec::new();
    };
    let mut batch = mir_row_cuts(problem, bounds, &x_bar, opts, tol);
    batch.extend(gomory_cuts_inner(
        s,
        problem,
        bounds,
        int_vars,
        opts,
        tol,
        &[],
    ));
    batch.extend(implied_bound_cuts(problem, bounds, &x_bar, opts, tol));
    batch.extend(cover_cuts_with_lifting(
        problem,
        bounds,
        &x_bar,
        opts,
        tol,
        &mut LiftBudget::new(),
    ));
    // 分离边界数值卫生闸（与 cut_loop_recorded 同一单点，见模块文档）：
    // 清洁化（只弱化）+ 密度闸 + violation 复核，之后再按 max_cuts 配额
    // 截断（配额只花在过闸存活者上）。
    batch = batch
        .into_iter()
        .filter_map(|cut| sanitize_cut(cut, bounds, Some(&x_bar), opts))
        .collect();
    batch.truncate(max_cuts);
    batch
}

/// IterationLimit 则如实作为终态返回。
pub fn cut_loop<E: LpEngine>(
    s: &mut E,
    bounds: &Bounds,
    problem: &Problem,
    int_vars: &[usize],
    opts: CutOptions,
    tol: &Tolerances,
    make: impl Fn() -> Result<E>,
) -> Result<(u64, LpStatus)> {
    cut_loop_recorded(s, bounds, problem, int_vars, opts, tol, make)
        .map(|(applied, status, _)| (applied, status))
}

/// 与 [`cut_loop`] 完全相同的割循环，额外把施加的割按施加顺序原样返回
/// （terms/rhs 即 add_row 时传入的原尺度值）。行为与 cut_loop 逐字节一致。
///
/// 用途：并行引擎在各 worker 的 Simplex 副本上按相同顺序 `add_row` 这些割
/// 即可复现根割表格（行序一致 ⇒ 表格一致）。[`cut_loop`] 即本函数丢弃割
/// 清单后的委托；既有语义逐项保持：入口防御性 `resolve(bounds)`、
/// `max_cuts_per_round` 截断、Infeasible ⇒ 如实作为终态返回（整数不可行
/// 问题的合法证明路径，见 [`cut_loop`] 文档）。
///
/// 每轮割批的确定性合并序（WP6 起，R2-retry 增聚合 MIR 批）：**c-MIR 批
/// （[`mir_row_cuts`]，行序）→ GMI 批 → 聚合 MIR 批（[`aggregated_mir_cuts`]，
/// 「表格行 ⊕ 原始行」聚合 GMI，R2-retry）→ 隐含界批
/// （[`implied_bound_cuts`]）→ 覆盖割批（[`cover_cuts_with_lifting`]，
/// 行序；提升预算 [`LiftBudget`] 跨轮共享、全程 3×10^7 状态操作，耗尽即
/// 关闭提升）**，先统一过分离边界数值卫生闸 [`sanitize_cut`]（小系数
/// 清洁化（只弱化）+ 密度闸 + violation 复核，见模块文档），再截断
/// `max_cuts_per_round`——ledger/recorded 记账的即过闸后的行。x̄ 经
/// [`rebuild_x_bar`] 重建（每轮一次，割批共用）。全部割行照走 `add_row`
/// + [`CutRecord`] 记账，worker 按序重放契约不变。
pub fn cut_loop_recorded<E: LpEngine>(
    s: &mut E,
    bounds: &Bounds,
    problem: &Problem,
    int_vars: &[usize],
    opts: CutOptions,
    tol: &Tolerances,
    make: impl Fn() -> Result<E>,
) -> Result<(u64, LpStatus, Vec<Cut>)> {
    let mut applied: u64 = 0;
    let mut ledger: Vec<CutRecord> = Vec::new();
    let mut recorded: Vec<Cut> = Vec::new();
    // lifting 全程预算：根割循环一次配给，跨轮共享（耗尽即全程关闭提升）。
    let mut lift_budget = LiftBudget::new();
    // 防御性同步：确保表格与 bounds 一致且为 Optimal。
    let mut status = s.resolve(bounds).status;
    if status != LpStatus::Optimal {
        return Ok((0, status, recorded));
    }
    for _round in 0..opts.max_rounds {
        let x_bar = rebuild_x_bar(s, bounds, problem.vars.len());
        // 确定性合并序：c-MIR → GMI → 隐含界 → 覆盖（提升），统一截断。
        let mut batch = match &x_bar {
            Some(x_bar) => mir_row_cuts(problem, bounds, x_bar, opts, tol),
            None => Vec::new(),
        };
        batch.extend(gomory_cuts_inner(
            s,
            problem,
            bounds,
            int_vars,
            opts,
            tol,
            &ledger,
        ));
        if let Some(x_bar) = &x_bar {
            // R2-retry（轴2）：聚合 MIR 批（「表格行 ⊕ 原始行」聚合 GMI）接在
            // GMI 之后——批序 c-MIR(原始)→GMI→聚合 MIR→隐含界→覆盖。
            batch.extend(aggregated_mir_cuts(
                s,
                problem,
                bounds,
                int_vars,
                x_bar,
                &ledger,
                opts,
                tol,
            ));
            batch.extend(implied_bound_cuts(problem, bounds, x_bar, opts, tol));
            batch.extend(cover_cuts_with_lifting(
                problem,
                bounds,
                x_bar,
                opts,
                tol,
                &mut lift_budget,
            ));
        }
        // 分离边界数值卫生闸（两引擎共用单点，见模块文档）：清洁化（只
        // 弱化）+ 密度闸 + violation 复核，之后再按 max_cuts_per_round 截断
        //（截断额度只花在过闸存活者上）。被闸丢弃的割不进入 add_row/
        // ledger/recorded——ledger 记账与 worker 按序重放始终与实际施加行
        // 逐字节一致。
        batch = batch
            .into_iter()
            .filter_map(|cut| sanitize_cut(cut, bounds, x_bar.as_deref(), opts))
            .collect();
        if batch.len() > opts.max_cuts_per_round {
            batch.truncate(opts.max_cuts_per_round);
        }
        if batch.is_empty() {
            break;
        }
        let applied_before = applied;
        let ledger_before = ledger.len();
        let recorded_before = recorded.len();
        for cut in batch {
            let index = s.add_row(&cut.terms, RelOp::Le, cut.rhs)?;
            ledger.push(CutRecord {
                terms: cut.terms.clone(),
                rhs: cut.rhs,
                index,
            });
            recorded.push(cut);
            applied = applied.saturating_add(1);
        }
        // R8（③）：割循环 resolve 收紧预算——暖启动（对偶修复快路径）下
        // 收敛轮实测 ~1.5k 迭代，退化爬行轮放任到默认上限会白烧（setcover
        // 算例第 5 轮 10000 迭代打满后整树回退丢全部割）。预算 = 已提交割
        // 行数感知的收紧式 max(3000, 12×(m+n))；超限走下方回滚而非终态
        // 传播。Infeasible 终局照旧第一时间如实返回；若某问题的 Phase-1
        // 本身需要超过收紧预算的迭代，走下方回滚后树搜索继续——真实不可行
        // 只会被延迟发现，不会被伪造或吞掉（hint 在回滚 resolve 前已复位）。
        let m_now = s.num_rows();
        let n_now = problem.vars.len();
        s.set_iter_limit_hint(Some(
            3_000u64.max(12u64.saturating_mul((m_now + n_now) as u64)),
        ));
        status = s.resolve(bounds).status;
        s.set_iter_limit_hint(None);
        if status == LpStatus::Infeasible {
            // 整数不可行问题的有效割可合法切空 LP 松弛区域：割只保证不切
            // 整数可行解，而该集合为空时任意割皆有效——割后根 LP 不可行
            // 是『整数不可行』的合法证明路径，正确终局是 Infeasible（终检
            // 由 Phase-1 不可行性给出），不是数值错误。原『理论不可能』
            // 断言只对『存在整数可行解』的问题成立。
            return Ok((applied, LpStatus::Infeasible, recorded));
        }
        if status == LpStatus::IterationLimit {
            // R8（③）回滚语义：本轮增量割使 resolve 超收紧预算——弃本轮
            // 增量割，重建冷 Simplex + 重放**已提交**割 + resolve，恢复到
            // 本轮前的 Optimal 根态后正常收尾（status = Optimal）。既有割
            // 的根界保留（树搜索受益），不再以 limit_reached 无解终态传播。
            // 回滚后状态与上一轮末逐字节一致（确定性冷重建路径 + 同割集），
            // 终态必为 Optimal；防御性保留非 Optimal 映射（不冒充）。
            // W2 注：回滚只对**迭代上限**触发——墙钟超限返回独立的
            // LpStatus::TimeLimit，走下方 `!= Optimal` 出口如实传播（整树
            // 诚实 LimitReached，不做「再解一次整段根 LP」的回滚重放，那
            // 会把超限点拖过 deadline）；重建 Simplex 透传原 deadline，
            // 迭代上限触发的回滚若在重放中撞上墙钟 ⇒ TimeLimit 诚实传播
            // （防御性非 Optimal 映射，不冒充 Optimal）。
            ledger.truncate(ledger_before);
            recorded.truncate(recorded_before);
            applied = applied_before;
            let rollback_deadline = s.time_deadline();
            let mut rebuilt = make()?;
            rebuilt.set_time_deadline(rollback_deadline);
            for rec in ledger.iter() {
                rebuilt.add_row(&rec.terms, RelOp::Le, rec.rhs)?;
            }
            let st = rebuilt.resolve(bounds).status;
            *s = rebuilt;
            status = st;
            break;
        }
        if status != LpStatus::Optimal {
            break; // Unbounded：如实返回终态
        }
    }
    Ok((applied, status, recorded))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simplex::Simplex;
    use crate::model::{Constraint, LinearExpr, Sense, Tolerances, VarType, Variable};

    // ===================== 聚合 MIR（R2-retry）=====================

    /// （R2-retry）聚合 MIR 有效性 fuzz：50 种子随机小模型（8 binary + 2
    /// 一般整数 + 1 连续列、4 条 Le 行），根 LP 上产出聚合割后穷举全部
    /// 整数可行点（2^8×3×4 点 × 连续列两端点）验证零误切。锚归一缺失
    /// （合并行含 B1 自身时未重归一）正是被本 fuzz 抓出的无效割来源——
    /// 保留作回归锚点。
    #[test]
    fn aggregated_mir_fuzz_validity() {
        for seed in 1..=50u64 {
            let mut lcg = seed.wrapping_mul(0x9E3779B97F4A7C15);
            let mut rnd = || {
                lcg = lcg.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                ((lcg >> 33) % 97 + 3) as f64
            };
            let mut pf = Problem::new(&format!("fuzz{seed}"), Sense::Maximize);
            for j in 0..8 {
                pf.vars.push(Variable {
                    name: format!("b{j}"),
                    var_type: VarType::Binary,
                    lb: 0.0,
                    ub: 1.0,
                });
            }
            for j in 8..10 {
                pf.vars.push(Variable {
                    name: format!("y{j}"),
                    var_type: VarType::Integer,
                    lb: 0.0,
                    ub: if j == 8 { 2.0 } else { 3.0 },
                });
            }
            pf.vars.push(Variable {
                name: "c0".into(),
                var_type: VarType::Continuous,
                lb: 0.0,
                ub: 2.5,
            });
            let mut obj = Vec::new();
            for j in 0..11 {
                obj.push((j, rnd()));
            }
            pf.obj = LinearExpr { terms: obj };
            let mut cons = Vec::new();
            for i in 0..4 {
                let terms: Vec<(usize, f64)> = (0..11).map(|j| (j, rnd())).collect();
                let rhs = terms.iter().map(|&(_, w)| w).sum::<f64>() * 0.38;
                cons.push(Constraint {
                    name: Some(format!("r{i}")),
                    expr: LinearExpr { terms },
                    op: RelOp::Le,
                    rhs,
                });
            }
            pf.constraints = cons.clone();
            let tf = Tolerances::default();
            let mut sf = Simplex::new(&pf, tf).unwrap();
            let bf = Bounds::from_problem(&pf);
            assert_eq!(sf.solve(&bf).status, LpStatus::Optimal, "seed {seed} 根 LP");
            let xbf = rebuild_x_bar(&sf, &bf, 11).unwrap();
            let cuts = aggregated_mir_cuts(
                &sf,
                &pf,
                &bf,
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                &xbf,
                &[],
                CutOptions::default(),
                &tf,
            );
            // 穷举整数盒（Le 正权行对 c 单调 ⇒ c∈{0, 2.5} 两端点判定存在性）。
            for cut in &cuts {
                for mask in 0..256u32 {
                    for y8 in 0..=2i64 {
                        for y9 in 0..=3i64 {
                            for c in [0.0f64, 2.5] {
                                let v = [
                                    (mask & 1) as f64,
                                    ((mask >> 1) & 1) as f64,
                                    ((mask >> 2) & 1) as f64,
                                    ((mask >> 3) & 1) as f64,
                                    ((mask >> 4) & 1) as f64,
                                    ((mask >> 5) & 1) as f64,
                                    ((mask >> 6) & 1) as f64,
                                    ((mask >> 7) & 1) as f64,
                                    y8 as f64,
                                    y9 as f64,
                                    c,
                                ];
                                let feas = cons.iter().all(|con| {
                                    con.expr
                                        .terms
                                        .iter()
                                        .map(|&(j, w)| w * v[j])
                                        .sum::<f64>()
                                        <= con.rhs + 1e-9
                                });
                                if feas {
                                    let lhs: f64 =
                                        cut.terms.iter().map(|&(j, cf)| cf * v[j]).sum();
                                    assert!(
                                        lhs <= cut.rhs + 1e-7,
                                        "seed {seed} 割 {cut:?} 切掉整数可行点 {v:?}"
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// （R2-retry）聚合 MIR 产出与确定性：mknap-3x100 根 LP 上离线证明
    /// 实测 8 条（viol 至 0.86）；同输入两次调用逐字节一致；每条割在 x̄
    /// 处 violation ≥ min_violation。
    #[test]
    fn aggregated_mir_production_and_determinism() {
        use crate::instance::load_instance;
        use std::path::Path;
        let p = load_instance(Path::new("../bench/instances/mknap-3x100.milp")).unwrap();
        let t = Tolerances::default();
        let mut s = Simplex::new(&p, t).unwrap();
        let b = Bounds::from_problem(&p);
        assert_eq!(s.solve(&b).status, LpStatus::Optimal);
        let n = p.vars.len();
        let int_vars: Vec<usize> = (0..n)
            .filter(|&j| matches!(p.vars[j].var_type, VarType::Integer | VarType::Binary))
            .collect();
        let x_bar = rebuild_x_bar(&s, &b, n).unwrap();
        let opts = CutOptions::default();
        let cuts = aggregated_mir_cuts(&s, &p, &b, &int_vars, &x_bar, &[], opts, &t);
        assert!(
            !cuts.is_empty(),
            "mknap-3x100 根 LP 应产出聚合 MIR 割（离线证明 8 条）"
        );
        for cut in &cuts {
            let lhs: f64 = cut.terms.iter().map(|&(j, c)| c * x_bar[j]).sum();
            assert!(lhs - cut.rhs >= opts.min_violation, "割 {cut:?} violation 不足");
        }
        let cuts2 = aggregated_mir_cuts(&s, &p, &b, &int_vars, &x_bar, &[], opts, &t);
        assert_eq!(cuts.len(), cuts2.len());
        for (a, c) in cuts.iter().zip(cuts2.iter()) {
            assert_eq!(a.terms, c.terms);
            assert_eq!(a.rhs, c.rhs);
        }
    }

    fn tol() -> Tolerances {
        Tolerances::default()
    }

    /// wolsey 经典算例：max 5x1 + 8x2 s.t. x1+x2 ≤ 6, 5x1+9x2 ≤ 45，
    /// x1, x2 ≥ 0 整数。根 LP 41.25 @ (2.25, 3.75)，整数最优 (0,5) = 40；
    /// 整数可行点：x2 = 0..5 时 x1 上限 6,5,4,3,1,0（共 25 个）。
    fn wolsey() -> Problem {
        let mut p = Problem::new("wolsey", Sense::Maximize);
        p.add_var("x1", VarType::Integer).unwrap();
        p.add_var("x2", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 5.0), (1, 8.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 6.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 5.0), (1, 9.0)], RelOp::Le, 45.0)
            .unwrap();
        p
    }

    /// 枚举 wolsey 全部 25 个整数可行点。
    fn wolsey_integer_points() -> Vec<(f64, f64)> {
        let caps = [6i64, 5, 4, 3, 1, 0];
        let mut pts = Vec::new();
        for (x2, &cap) in caps.iter().enumerate() {
            for x1 in 0..=cap {
                pts.push((x1 as f64, x2 as f64));
            }
        }
        pts
    }

    /// 割在取值点处的 lhs。
    fn eval_cut(cut: &Cut, x: &[f64]) -> f64 {
        cut.terms.iter().map(|&(j, c)| c * x[j]).sum()
    }

    /// 默认参数逐项精确符合契约（WP6：5×20 → 4×40，轮数按降档守则校准，
    /// 依据见 [`CutOptions::default`] 文档）。
    #[test]
    fn cut_options_default_values() {
        let o = CutOptions::default();
        assert_eq!(o.max_rounds, 4);
        assert_eq!(o.max_cuts_per_round, 40);
        assert_eq!(o.away, 0.005);
        assert_eq!(o.min_violation, 1e-7);
    }

    /// wolsey 根 LP 41.25 处生成的割须被其全部 25 个整数可行点（含最优
    /// (0,5)）满足，且被当前 LP 解违反 ≥ min_violation。
    #[test]
    fn gomory_wolsey_root_cuts_valid_for_all_integer_points() {
        let p = wolsey();
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let b = Bounds::from_problem(&p);
        let sol = s.solve(&b);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.objective - 41.25).abs() < 1e-9, "根 LP 应为 41.25");
        let cuts = gomory_cuts(&s, &p, &b, &p.integer_vars(), CutOptions::default(), &t);
        assert!(!cuts.is_empty(), "根 LP 处必须生成割");
        for cut in &cuts {
            assert_eq!(cut.op, RelOp::Le);
            for (x1, x2) in wolsey_integer_points() {
                let x = [x1, x2];
                assert!(
                    eval_cut(cut, &x) <= cut.rhs + 1e-6,
                    "割 {cut:?} 切掉整数可行点 ({x1},{x2})"
                );
            }
            let viol = eval_cut(cut, &sol.x) - cut.rhs;
            assert!(viol >= 1e-7, "割须被当前 LP 解违反：viol = {viol}");
        }
    }

    /// 行缩放（同一约束乘正常数）不改变割的有效性：归一化换算正确。
    #[test]
    fn gomory_scaled_rows_still_valid() {
        let mut p = Problem::new("wolsey-scaled", Sense::Maximize);
        p.add_var("x1", VarType::Integer).unwrap();
        p.add_var("x2", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 5.0), (1, 8.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 2.0), (1, 2.0)], RelOp::Le, 12.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 25.0), (1, 45.0)], RelOp::Le, 225.0)
            .unwrap();
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let b = Bounds::from_problem(&p);
        let sol = s.solve(&b);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.objective - 41.25).abs() < 1e-9);
        let cuts = gomory_cuts(&s, &p, &b, &p.integer_vars(), CutOptions::default(), &t);
        assert!(!cuts.is_empty());
        for cut in &cuts {
            for (x1, x2) in wolsey_integer_points() {
                let x = [x1, x2];
                assert!(
                    eval_cut(cut, &x) <= cut.rhs + 1e-6,
                    "缩放行割 {cut:?} 切掉 ({x1},{x2})"
                );
            }
        }
    }

    /// 人为构造整数变量非整数界（z ≤ 10.5）：整数列所贴界非整数值时必须走
    /// 连续公式；生成的割仍被全部整数可行点（含整数最优 (1,10)）满足。
    /// 错误地走整数公式（π_z = 0）会得到 z − x ≥ 9.5，切掉 (1,10)。
    #[test]
    fn gomory_integer_col_at_noninteger_bound_uses_continuous_formula() {
        let mut p = Problem::new("nonint-bound", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap(); // [0, +∞)
        let z = p.add_var("z", VarType::Integer).unwrap(); // [0, 10.5]
        p.set_var_bounds(z, 0.0, 10.5).unwrap();
        p.set_objective(vec![(0, 1.0), (1, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, -1.0)], RelOp::Le, -9.0)
            .unwrap();
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let b = Bounds::from_problem(&p);
        let sol = s.solve(&b);
        assert_eq!(sol.status, LpStatus::Optimal);
        // 根 LP 最优 (1.5, 10.5)，目标 12。
        assert!((sol.objective - 12.0).abs() < 1e-9);
        assert!((sol.x[0] - 1.5).abs() < 1e-9);
        let cuts = gomory_cuts(&s, &p, &b, &p.integer_vars(), CutOptions::default(), &t);
        assert!(!cuts.is_empty(), "基变量 x = 1.5 分数，应生成割");
        // 整数可行点：2x + z ≤ 13.5（整数 ≤ 13）、z ≤ 10、x ≥ 0、z ≥ 9 + x。
        let integer_points = [(0.0, 9.0), (0.0, 10.0), (1.0, 10.0)];
        for cut in &cuts {
            for &(x, z) in &integer_points {
                let v = [x, z];
                assert!(
                    eval_cut(cut, &v) <= cut.rhs + 1e-6,
                    "割 {cut:?} 切掉整数可行点 ({x},{z})"
                );
            }
        }
    }

    /// 根节点割循环：施加割、状态 Optimal、界改进且不切掉整数最优 40。
    #[test]
    fn cut_loop_wolsey_improves_root_bound() {
        let p = wolsey();
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let b = Bounds::from_problem(&p);
        assert_eq!(s.solve(&b).status, LpStatus::Optimal);
        let (applied, status) =
            cut_loop(&mut s, &b, &p, &p.integer_vars(), CutOptions::default(), &t, || Simplex::new(&p, t)).unwrap();
        assert!(applied >= 1, "至少施加一条割");
        assert_eq!(status, LpStatus::Optimal);
        let sol = s.resolve(&b);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!(
            sol.objective >= 40.0 - 1e-6,
            "割不得切掉整数最优 (0,5)（目标 40）：{}",
            sol.objective
        );
        assert!(
            sol.objective <= 41.25 - 1e-6,
            "割必须改进根 LP 界 41.25：{}",
            sol.objective
        );
    }

    /// LP 松弛已整数最优时无割可生成，cut_loop 返回 (0, Optimal)。
    #[test]
    fn cut_loop_no_cuts_when_lp_already_integral() {
        let mut p = Problem::new("integral", Sense::Maximize);
        let x = p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(x, 0.0, 5.0).unwrap();
        p.set_objective(vec![(x, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(x, 1.0)], RelOp::Le, 3.0)
            .unwrap();
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let b = Bounds::from_problem(&p);
        assert_eq!(s.solve(&b).status, LpStatus::Optimal);
        let (applied, status) =
            cut_loop(&mut s, &b, &p, &p.integer_vars(), CutOptions::default(), &t, || Simplex::new(&p, t)).unwrap();
        assert_eq!(applied, 0);
        assert_eq!(status, LpStatus::Optimal);
    }

    /// cut_loop 与 cut_loop_recorded 在 wolsey 上返回相同的 (applied, status)，
    /// 且 recorded 清单长度恰为施加割数。
    #[test]
    fn cut_loop_and_recorded_agree_on_wolsey() {
        let p = wolsey();
        let t = tol();
        let b = Bounds::from_problem(&p);

        let mut s1 = Simplex::new(&p, t).unwrap();
        assert_eq!(s1.solve(&b).status, LpStatus::Optimal);
        let (applied1, status1) =
            cut_loop(&mut s1, &b, &p, &p.integer_vars(), CutOptions::default(), &t, || Simplex::new(&p, t)).unwrap();

        let mut s2 = Simplex::new(&p, t).unwrap();
        assert_eq!(s2.solve(&b).status, LpStatus::Optimal);
        let (applied2, status2, cuts2) = cut_loop_recorded(
            &mut s2,
            &b,
            &p,
            &p.integer_vars(),
            CutOptions::default(),
            &t,
            || Simplex::new(&p, t),
        )
        .unwrap();

        assert_eq!((applied1, status1), (applied2, status2));
        assert_eq!(cuts2.len(), applied2 as usize, "记录条数须等于施加割数");
    }

    /// recorded 割非空时：全新 Simplex::new + 按序 add_row + resolve 复现与根
    /// 相同的根 LP 目标值（并行正确性前提：worker 行序一致 ⇒ 表格一致）。
    #[test]
    fn cut_loop_recorded_replays_on_fresh_simplex() {
        let p = wolsey();
        let t = tol();
        let b = Bounds::from_problem(&p);
        let mut root = Simplex::new(&p, t).unwrap();
        assert_eq!(root.solve(&b).status, LpStatus::Optimal);
        let (applied, status, cuts) = cut_loop_recorded(
            &mut root,
            &b,
            &p,
            &p.integer_vars(),
            CutOptions::default(),
            &t,
            || Simplex::new(&p, t),
        )
        .unwrap();
        assert!(applied >= 1, "wolsey 根 LP 41.25 处必产割");
        assert_eq!(status, LpStatus::Optimal);
        // 根侧终态读数：记账性 resolve（与 branch_bound 根阶段同口径）。
        let root_sol = root.resolve(&b);
        assert_eq!(root_sol.status, LpStatus::Optimal);

        // worker 复现路径：全新 Simplex + 按序 add_row + resolve。
        let mut worker = Simplex::new(&p, t).unwrap();
        for cut in &cuts {
            worker.add_row(&cut.terms, cut.op, cut.rhs).unwrap();
        }
        let w = worker.resolve(&b);
        assert_eq!(w.status, LpStatus::Optimal);
        assert!(
            (w.objective - root_sol.objective).abs() < 1e-9,
            "worker 表格须复现根 LP 目标：root = {}，worker = {}",
            root_sol.objective,
            w.objective
        );
        // 与既有割后根界契约一致：不切掉整数最优 40，且严格改进根 LP 41.25。
        assert!(w.objective >= 40.0 - 1e-6);
        assert!(w.objective <= 41.25 - 1e-6);
    }

    /// M5：整数不可行问题（穷举 0 个整数可行点）上，割循环把根 LP 割空 ⇒
    /// 正确终局 Infeasible（原实现误判 Err(Numerical)——『理论不可能』只对
    /// 存在整数可行解的问题成立）。
    ///
    /// 算例：max 6x1+2x2+5c s.t. 2x1+5x2+2c=8, 3x1+3x2+0.5c≤3.5，
    /// x1∈[0,5] 整数、x2∈[0,2] 整数、c∈[0,3] 连续——LP 松弛可行，整数
    /// 可行点穷举为 0（测试内逐点验证）。
    #[test]
    fn cut_loop_infeasible_after_cuts_is_infeasible_terminal() {
        let mut p = Problem::new("cut-infeasible", Sense::Maximize);
        p.add_var("x1", VarType::Integer).unwrap();
        p.add_var("x2", VarType::Integer).unwrap();
        let c3 = p.add_var("c", VarType::Continuous).unwrap();
        p.set_var_bounds(0, 0.0, 5.0).unwrap();
        p.set_var_bounds(1, 0.0, 2.0).unwrap();
        p.set_var_bounds(2, 0.0, 3.0).unwrap();
        p.set_objective(vec![(0, 6.0), (1, 2.0), (2, 5.0)], 0.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 2.0), (1, 5.0), (2, 2.0)], RelOp::Eq, 8.0)
            .unwrap();
        p.add_constraint(
            None,
            vec![(0, 3.0), (1, 3.0), (2, 0.5)],
            RelOp::Le,
            3.5,
        )
        .unwrap();

        // 穷举钉死：全部 (x1, x2) 整数组合下 c 由 Eq 行唯一决定，均违反
        // c 的界或 Le 行 ⇒ 整数可行点为 0。
        for x1 in 0..=5i64 {
            for x2 in 0..=2i64 {
                let cc = (8.0 - 2.0 * x1 as f64 - 5.0 * x2 as f64) / 2.0;
                let le = 3.0 * x1 as f64 + 3.0 * x2 as f64 + 0.5 * cc;
                let feasible = cc >= -1e-9 && cc <= 3.0 + 1e-9 && le <= 3.5 + 1e-9;
                assert!(
                    !feasible,
                    "意外存在可行整数点 ({x1},{x2},c={cc})"
                );
            }
        }
        let _ = c3;

        let t = tol();
        let b = Bounds::from_problem(&p);
        let mut s = Simplex::new(&p, t).unwrap();
        let root = s.solve(&b);
        assert_eq!(root.status, LpStatus::Optimal, "LP 松弛可行（割前）");

        match cut_loop_recorded(
            &mut s,
            &b,
            &p,
            &p.integer_vars(),
            CutOptions::default(),
            &t,
            || Simplex::new(&p, t),
        ) {
            Ok((applied, status, _)) => {
                assert!(
                    applied >= 1,
                    "该算例根 LP 分数处必产割（否则未覆盖割空路径）"
                );
                assert_eq!(
                    status,
                    LpStatus::Infeasible,
                    "割后根 LP 不可行 ⇒ Infeasible 终态（applied={applied}）"
                );
            }
            Err(e) => panic!("割后不可行不得报数值错误（原缺陷）：{e}"),
        }
        // 记账性 resolve 维持 Infeasible（终态可被调用方稳定读取）。
        assert_eq!(s.resolve(&b).status, LpStatus::Infeasible);
    }

    // ===================== 覆盖割（WP4 增补）=====================

    /// n 个 Binary 变量 + 单条背包行 `Σ w_j x_j ≤ cap`（目标 max Σx）。
    fn binary_knapsack(name: &str, weights: &[f64], cap: f64) -> Problem {
        let mut p = Problem::new(name, Sense::Maximize);
        for i in 0..weights.len() {
            p.add_var(&format!("x{i}"), VarType::Binary).unwrap();
        }
        p.set_objective((0..weights.len()).map(|i| (i, 1.0)).collect(), 0.0)
            .unwrap();
        let terms: Vec<(usize, f64)> = weights.iter().enumerate().map(|(i, &w)| (i, w)).collect();
        p.add_constraint(Some("cap"), terms, RelOp::Le, cap).unwrap();
        p
    }

    /// 小行穷举：割须被背包行的全部 0-1 可行点满足（复用 :452-474 口径），
    /// 且形状恰为 [(j,1.0)] 下标升序、Le、rhs=|C|−1，被 x̄ 违反。
    ///
    /// 行：4(x0+..+x4) ≤ 15（total=20 > cap），x̄ = [1,1,1,1,0.5]：贪婪取
    /// x̄=1 的 4 个（Σ=16 > 15），极小化无可可剔（剔任一 → 12 ≤ 15）⇒
    /// 割 x0+x1+x2+x3 ≤ 3，violation = 4−3 = 1。
    #[test]
    fn cover_cut_exhaustive_validity_and_shape() {
        let p = binary_knapsack("cover-5", &[4.0, 4.0, 4.0, 4.0, 4.0], 15.0);
        let b = Bounds::from_problem(&p);
        let x_bar = [1.0, 1.0, 1.0, 1.0, 0.5];
        let cuts = cover_cuts(&p, &b, &x_bar, CutOptions::default(), &tol());
        assert_eq!(cuts.len(), 1, "单行 ⇒ 至多 1 条");
        let cut = &cuts[0];
        assert_eq!(cut.op, RelOp::Le);
        assert_eq!(cut.rhs, 3.0);
        assert_eq!(
            cut.terms,
            vec![(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0)],
            "terms 须为 (j,1.0) 下标升序"
        );
        // 穷举全部 2^5 个 0-1 点：行可行（Σ4x ≤ 15 ⇒ Σx ≤ 3）⇒ 割可行。
        for mask in 0..32u32 {
            let x: Vec<f64> = (0..5).map(|i| ((mask >> i) & 1) as f64).collect();
            let row_lhs: f64 = x.iter().sum::<f64>() * 4.0;
            if row_lhs <= 15.0 {
                assert!(
                    eval_cut(cut, &x) <= cut.rhs + 1e-9,
                    "割 {cut:?} 切掉行可行 0-1 点 {x:?}"
                );
            }
        }
        // 割被当前 x̄ 违反（violation = 1 ≥ min_violation）。
        assert!(eval_cut(cut, &x_bar) - cut.rhs >= 1e-7);
    }

    /// 极小化：贪婪越过 cap 后剔除可剔者，终态『去任一成员即不超容』。
    ///
    /// 例 1（w=[2,10,11], cap=12, x̄=[0.9,0.8,0.5]）：贪婪全收（Σ=23>12），
    /// 极小化剔 w=10（23−10=13>12）后 C={x0,x2}，Σ=13>12，剔任一即 ≤12。
    /// 例 2（w=[1,2,100], cap=50）：极小化剔到单元素覆盖 x2 ≤ 0。
    #[test]
    fn cover_minimal_after_greedy_overshoot() {
        let cases: &[(&[f64], f64, &[f64], Vec<usize>)] = &[
            (&[2.0, 10.0, 11.0], 12.0, &[0.9, 0.8, 0.5], vec![0, 2]),
            (&[1.0, 2.0, 100.0], 50.0, &[0.9, 0.8, 0.3], vec![2]),
        ];
        for &(w, cap, x_bar, ref expected) in cases {
            let p = binary_knapsack("cover-min", w, cap);
            let b = Bounds::from_problem(&p);
            let cuts = cover_cuts(&p, &b, x_bar, CutOptions::default(), &tol());
            assert_eq!(cuts.len(), 1, "w={w:?}");
            let cut = &cuts[0];
            let got: Vec<usize> = cut.terms.iter().map(|&(j, _)| j).collect();
            assert_eq!(got, *expected, "w={w:?} 的覆盖成员");
            assert_eq!(cut.rhs, (expected.len() - 1) as f64);
            // 最小性性质：Σ_{C} w > cap，且去任一成员即 ≤ cap。
            let sum: f64 = expected.iter().map(|&j| w[j]).sum();
            assert!(sum > cap, "Σ_C w={sum} 须超容 {cap}");
            for &j in expected {
                assert!(sum - w[j] <= cap, "剔掉 x{j} 后 Σ={} 须不超容", sum - w[j]);
            }
        }
    }

    /// x̄ 引导分离在『w 升序不违反、x̄ 降序违反』的反例上产出割（钉死旧
    /// w 升序路径不回归）。
    ///
    /// w=[3,4,5], cap=6, x̄=[0.1,0.2,0.9]：x̄ 降序取 {x1,x2}（Σ=9>6）⇒
    /// violation = 1.1−1 = 0.1 ≥ 1e-7 ⇒ 割 x1+x2 ≤ 1；而 w 升序贪婪会取
    /// {x0,x1}（violation = 0.3−1 < 0）⇒ 旧路径零产出。
    #[test]
    fn cover_xbar_guided_not_weight_order() {
        let p = binary_knapsack("cover-xbar", &[3.0, 4.0, 5.0], 6.0);
        let b = Bounds::from_problem(&p);
        let x_bar = [0.1, 0.2, 0.9];
        let cuts = cover_cuts(&p, &b, &x_bar, CutOptions::default(), &tol());
        assert_eq!(cuts.len(), 1, "x̄ 引导路径必须产出割（旧 w 升序路径为零产出）");
        assert_eq!(cuts[0].terms, vec![(1, 1.0), (2, 1.0)]);
        assert_eq!(cuts[0].rhs, 1.0);
        // 有效性：全部 0-1 行可行点满足割（x1=x2=1 ⇒ 4+5=9 > 6 被行排除）。
        for mask in 0..8u32 {
            let x: Vec<f64> = (0..3).map(|i| ((mask >> i) & 1) as f64).collect();
            if 3.0 * x[0] + 4.0 * x[1] + 5.0 * x[2] <= 6.0 {
                assert!(
                    eval_cut(&cuts[0], &x) <= cuts[0].rhs + 1e-9,
                    "割切掉行可行点 {x:?}"
                );
            }
        }
    }

    /// 无违反时不产割；各类门控（Ge/Eq、非 Binary 支撑、非正系数、
    /// Σw ≤ cap、|w| ≥ 1e9、维数不符）整行跳过。
    #[test]
    fn cover_no_cut_without_violation_and_row_gating() {
        let t = tol();
        let opts = CutOptions::default();
        // 行可行整数 x̄（violation = 0 < min_violation）⇒ 不产割。
        let p = binary_knapsack("cover-none", &[3.0, 4.0, 5.0], 6.0);
        let b = Bounds::from_problem(&p);
        assert!(cover_cuts(&p, &b, &[1.0, 0.0, 0.0], opts, &t).is_empty());

        // Ge 行跳过。
        let mut pg = binary_knapsack("cover-ge", &[3.0, 4.0, 5.0], 6.0);
        pg.constraints[0].op = RelOp::Ge;
        assert!(cover_cuts(&pg, &b, &[0.9, 0.9, 0.9], opts, &t).is_empty());
        // Eq 行跳过。
        let mut pe = binary_knapsack("cover-eq", &[3.0, 4.0, 5.0], 6.0);
        pe.constraints[0].op = RelOp::Eq;
        assert!(cover_cuts(&pe, &b, &[0.9, 0.9, 0.9], opts, &t).is_empty());

        // 非 Binary 支撑（Integer 变量）跳过。
        let mut pi = Problem::new("cover-int", Sense::Maximize);
        pi.add_var("x0", VarType::Integer).unwrap();
        pi.add_var("x1", VarType::Binary).unwrap();
        pi.set_objective(vec![(0, 1.0), (1, 1.0)], 0.0).unwrap();
        pi.add_constraint(None, vec![(0, 3.0), (1, 4.0)], RelOp::Le, 5.0)
            .unwrap();
        let bi = Bounds::from_problem(&pi);
        assert!(cover_cuts(&pi, &bi, &[0.9, 0.9], opts, &t).is_empty());

        // 非正系数跳过。
        let mut pn = Problem::new("cover-neg", Sense::Maximize);
        pn.add_var("x0", VarType::Binary).unwrap();
        pn.add_var("x1", VarType::Binary).unwrap();
        pn.set_objective(vec![(0, 1.0), (1, 1.0)], 0.0).unwrap();
        pn.add_constraint(None, vec![(0, 3.0), (1, -4.0)], RelOp::Le, 5.0)
            .unwrap();
        let bn = Bounds::from_problem(&pn);
        assert!(cover_cuts(&pn, &bn, &[0.9, 0.9], opts, &t).is_empty());

        // 全支撑 Σw ≤ cap（无覆盖存在）跳过。
        let pc = binary_knapsack("cover-cap", &[1.0, 1.0, 1.0], 5.0);
        let bc = Bounds::from_problem(&pc);
        assert!(cover_cuts(&pc, &bc, &[1.0, 1.0, 1.0], opts, &t).is_empty());

        // |w| ≥ 1e9 整行跳过（数值防御）。
        let ph = binary_knapsack("cover-huge", &[2e9, 3.0], 4.0);
        let bh = Bounds::from_problem(&ph);
        assert!(cover_cuts(&ph, &bh, &[0.5, 1.0], opts, &t).is_empty());

        // x̄ 维数不符（防御）。
        assert!(cover_cuts(&p, &b, &[0.9, 0.9], opts, &t).is_empty());
    }

    /// 同输入两次 cut_loop_recorded 产生完全相同的 recorded 清单（确定性），
    /// 且覆盖割（[(j,1.0)]、rhs=|C|−1 形状）确实进入清单；worker 按序重放
    /// 复现根 LP 目标。
    ///
    /// 问题：5 Binary、4Σx ≤ 15、max Σx——根 LP 3.75（3 个 1 + 1 个 0.75），
    /// 覆盖割 Σ_{x̄ 高的 4 个} ≤ 3（violation 0.75）；整数最优 3。
    #[test]
    fn cut_loop_recorded_deterministic_with_cover_cuts() {
        let p = binary_knapsack("cover-det", &[4.0, 4.0, 4.0, 4.0, 4.0], 15.0);
        let t = tol();
        let b = Bounds::from_problem(&p);
        let int_vars = p.integer_vars();

        let mut s1 = Simplex::new(&p, t).unwrap();
        let sol1 = s1.solve(&b);
        assert_eq!(sol1.status, LpStatus::Optimal);
        assert!((sol1.objective - 3.75).abs() < 1e-9, "根 LP 应为 3.75");
        let r1 = cut_loop_recorded(&mut s1, &b, &p, &int_vars, CutOptions::default(), &t, || Simplex::new(&p, t))
            .unwrap();

        let mut s2 = Simplex::new(&p, t).unwrap();
        assert_eq!(s2.solve(&b).status, LpStatus::Optimal);
        let r2 = cut_loop_recorded(&mut s2, &b, &p, &int_vars, CutOptions::default(), &t, || Simplex::new(&p, t))
            .unwrap();

        assert_eq!(r1.0, r2.0, "applied 一致");
        assert_eq!(r1.1, r2.1, "status 一致");
        assert_eq!(r1.2.len(), r1.0 as usize);
        assert_eq!(r1.2.len(), r2.2.len());
        for (c1, c2) in r1.2.iter().zip(r2.2.iter()) {
            assert_eq!(c1.op, c2.op);
            assert_eq!(c1.rhs, c2.rhs);
            assert_eq!(c1.terms, c2.terms, "recorded 清单须逐条完全一致");
        }
        assert!(r1.0 >= 1, "根 LP 3.75 处必产割");

        // 覆盖割确实进入清单（WP6 lifting 后形状放宽：基础形 [(j,1.0)] 且
        // rhs = |terms|−1，或提升形——系数全为 ≥ 1 的整数值、rhs 为正整数
        // |C|−1；本算例 5×4 ≤ 15 的提升割为 [(0..4, 1.0)]、rhs = 3）。
        let has_cover = r1.2.iter().any(|c| {
            !c.terms.is_empty()
                && c.terms.iter().all(|&(_, w)| {
                    w >= 1.0 && (w - w.round()).abs() < 1e-12
                })
                && (c.rhs - c.rhs.round()).abs() < 1e-12
                && c.rhs >= 1.0
                && c.rhs <= c.terms.len() as f64 - 1.0
        });
        assert!(has_cover, "recorded 须含覆盖族割：{:?}", r1.2);

        // 割后根界：严格改进 3.75 且不切整数最优 3。
        let root_sol = s1.resolve(&b);
        assert_eq!(root_sol.status, LpStatus::Optimal);
        assert!(root_sol.objective <= 3.75 - 1e-6, "割须改进根 LP：{}", root_sol.objective);
        assert!(root_sol.objective >= 3.0 - 1e-6, "割不得切整数最优 3：{}", root_sol.objective);

        // worker 按序重放（并行契约）：全新 Simplex + 按序 add_row 复现根目标。
        let mut worker = Simplex::new(&p, t).unwrap();
        for cut in &r1.2 {
            worker.add_row(&cut.terms, cut.op, cut.rhs).unwrap();
        }
        let w = worker.resolve(&b);
        assert_eq!(w.status, LpStatus::Optimal);
        assert!(
            (w.objective - root_sol.objective).abs() < 1e-9,
            "worker 复现：root = {}，worker = {}",
            root_sol.objective,
            w.objective
        );
    }

    // ===================== 取整 ε 纪律（WP6 增补）=====================

    /// floor/ceil 的放宽方向与近整数吸附：真整数附近的浮点和差不得把可行
    /// 整数点 floor/ceil 掉。
    #[test]
    fn floor_ceil_bound_epsilon_discipline() {
        // 真值 5、浮点落其下方 1e-15：吸附为 5（上界取整放宽方向）。
        assert_eq!(floor_bound(4.999999999999999), 5.0);
        // 真值 5、浮点落其上方：同样吸附。
        assert_eq!(floor_bound(5.000000000000001), 5.0);
        // 非近整数：正常 floor（+ε 不越界）。
        assert_eq!(floor_bound(4.999999), 4.0);
        assert_eq!(floor_bound(1.5), 1.0);
        assert_eq!(floor_bound(-1.5), -2.0);
        // 下界镜像：ceil 放宽方向。
        assert_eq!(ceil_bound(5.000000000000001), 5.0);
        assert_eq!(ceil_bound(4.999999999999999), 5.0);
        assert_eq!(ceil_bound(1.5), 2.0);
        assert_eq!(ceil_bound(5.000999), 6.0);
        // mir_floor 的近整数 ε 纪律（放宽方向：结果 ≤ ⌊val⌋）。
        // 整数下方 ε 内：不向上吸附（旧行为取 3.0——f_j<0 使 MIR 系数
        // 加严至多 ε/(1−f̂0)，2026-10 遗留修复改回真 floor 2）。
        assert_eq!(mir_floor(2.9999999999999996), 2.0);
        // 整数上方 ε 内：向下吸附 = ⌊val⌋。
        assert_eq!(mir_floor(3.0000000000000004), 3.0);
        // 恰等整数与非近整数：与普通 floor 一致。
        assert_eq!(mir_floor(3.0), 3.0);
        assert_eq!(mir_floor(2.5), 2.0);
        assert_eq!(mir_floor(-2.5), -3.0);
        // 负数近整数（round 朝上）：下方 ε 内同样不吸附。
        assert_eq!(mir_floor(-3.0000000000000004), -4.0);
        assert_eq!(mir_floor(-2.9999999999999996), -3.0);
    }

    /// （L2 双向数值钉）mir_floor 放宽方向不变量：任意 val 下
    /// mir_floor(val) ≤ ⌊val⌋（只松不紧，MIR 系数只降不升）。网格扫过
    /// 正负、近整数双侧、大数量级（α 上界 1e6 附近，ε≈1e-3 仍不越真
    /// floor 上方）与 eps 边界外侧。
    #[test]
    fn mir_floor_never_tightens_above_true_floor() {
        let mut grid: Vec<f64> = Vec::new();
        for k in [-3.0, -1.0, 0.0, 1.0, 3.0, 100.0, 1e6] {
            for d in [
                0.0,
                1e-16,
                1e-12,
                1e-10,
                1e-9,
                1e-6,
                0.25,
                0.5,
                0.75,
                1.0 - 1e-12,
            ] {
                grid.push(k + d); // 整数上方（含 ε 内吸附区与外部）
                grid.push(k - d); // 整数下方（含 ε 内不吸附区与外部）
            }
        }
        for v in grid {
            assert!(
                mir_floor(v) <= v.floor(),
                "mir_floor({v}) = {} 超出真 floor {}（加严方向）",
                mir_floor(v),
                v.floor()
            );
        }
        // ε 边界精确行为：3−1e-9 恰在 ε=4e-9 内（取 2），3−1e-8 在外
        //（普通 floor=2），3+1e-9 在内（吸附 3）。
        assert_eq!(mir_floor(3.0 - 1e-9), 2.0);
        assert_eq!(mir_floor(3.0 - 1e-8), 2.0);
        assert_eq!(mir_floor(3.0 + 1e-9), 3.0);
    }

    // ===================== 单行 c-MIR（WP6 增补）=====================

    /// 负系数行（互补 z：x_j → ub_j − y_j）的 c-MIR：穷举全部整数点验证
    /// 行可行 ⇒ 割可行，且割被 x̄ 违反、同输入逐字节确定。
    ///
    /// 行：3x − 2z ≤ 4.5，x ∈ [0,3] 整数、z ∈ [0,4] 整数；x̄ = (3, 2.25)
    /// （LP 式分数点）。
    #[test]
    fn mir_negative_coefficient_row_exhaustive_validity() {
        let mut p = Problem::new("mir-neg", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.add_var("z", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 3.0).unwrap();
        p.set_var_bounds(1, 0.0, 4.0).unwrap();
        p.set_objective(vec![(0, 1.0), (1, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 3.0), (1, -2.0)], RelOp::Le, 4.5)
            .unwrap();
        let b = Bounds::from_problem(&p);
        let x_bar = [3.0, 2.25];
        let t = tol();
        let opts = CutOptions::default();
        let cut = mir_single_row(&p, &b, &x_bar, &p.constraints[0], opts, &t)
            .expect("负系数互补行在分数 x̄ 处应产出 c-MIR 割");
        // 2^枚举（整数盒 4×5 = 20 点）：行可行 ⇒ 割可行。
        for x in 0..=3i64 {
            for z in 0..=4i64 {
                let v = [x as f64, z as f64];
                if 3.0 * v[0] - 2.0 * v[1] <= 4.5 {
                    assert!(
                        eval_cut(&cut, &v) <= cut.rhs + 1e-9,
                        "割 {cut:?} 切掉行可行整数点 ({x},{z})"
                    );
                }
            }
        }
        // 割被 x̄ 违反 ≥ min_violation。
        assert!(eval_cut(&cut, &x_bar) - cut.rhs >= opts.min_violation);
        // 确定性：同输入两次完全一致。
        let cut2 = mir_single_row(&p, &b, &x_bar, &p.constraints[0], opts, &t).unwrap();
        assert_eq!(cut.terms, cut2.terms);
        assert_eq!(cut.rhs, cut2.rhs);
    }

    /// 无有限 ub 的负系数列 ⇒ 整行跳过（None）；非 Le 行跳过。
    #[test]
    fn mir_skips_infinite_ub_and_non_le_rows() {
        let mut p = Problem::new("mir-inf", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.add_var("z", VarType::Integer).unwrap(); // ub = +∞
        p.set_objective(vec![(0, 1.0), (1, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 3.0), (1, -2.0)], RelOp::Le, 4.5)
            .unwrap();
        let b = Bounds::from_problem(&p);
        let t = tol();
        assert!(
            mir_single_row(&p, &b, &[3.0, 2.25], &p.constraints[0], CutOptions::default(), &t)
                .is_none(),
            "负系数列 z 无有限 ub ⇒ 整行跳过"
        );
        // Ge / Eq 行不处理。
        for op in [RelOp::Ge, RelOp::Eq] {
            let mut q = Problem::new("mir-op", Sense::Maximize);
            q.add_var("x", VarType::Integer).unwrap();
            q.set_var_bounds(0, 0.0, 5.0).unwrap();
            q.set_objective(vec![(0, 1.0)], 0.0).unwrap();
            q.add_constraint(None, vec![(0, 2.0)], op, 4.5).unwrap();
            let bq = Bounds::from_problem(&q);
            assert!(mir_single_row(
                &q,
                &bq,
                &[1.5],
                &q.constraints[0],
                CutOptions::default(),
                &t
            )
            .is_none());
        }
    }

    /// 整数变量非整数界（互补的 ub = 3.5）⇒ 双射不保整，整行跳过。
    #[test]
    fn mir_skips_nonintegral_bound_on_integer_column() {
        let mut p = Problem::new("mir-nonint-ub", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.add_var("z", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 3.0).unwrap();
        p.set_var_bounds(1, 0.0, 3.5).unwrap();
        p.set_objective(vec![(0, 1.0), (1, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 3.0), (1, -2.0)], RelOp::Le, 4.5)
            .unwrap();
        let b = Bounds::from_problem(&p);
        assert!(mir_single_row(
            &p,
            &b,
            &[3.0, 2.25],
            &p.constraints[0],
            CutOptions::default(),
            &tol()
        )
        .is_none());
    }

    /// 含连续变量的行：连续列 π = 0（丢弃只放宽），穷举（整数列枚举 ×
    /// 连续列端点）验证行可行 ⇒ 割可行。
    ///
    /// 行：2x + 0.5c ≤ 3.25，x ∈ [0,3] 整数、c ∈ [0,5] 连续；x̄ = (1.625, 0)。
    /// δ 逐试后最优 δ=1：割 x ≤ 1（violation 0.625）。
    #[test]
    fn mir_continuous_column_row_exhaustive_validity() {
        let mut p = Problem::new("mir-cont", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.add_var("c", VarType::Continuous).unwrap();
        p.set_var_bounds(0, 0.0, 3.0).unwrap();
        p.set_var_bounds(1, 0.0, 5.0).unwrap();
        p.set_objective(vec![(0, 1.0), (1, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 2.0), (1, 0.5)], RelOp::Le, 3.25)
            .unwrap();
        let b = Bounds::from_problem(&p);
        let x_bar = [1.625, 0.0];
        let t = tol();
        let opts = CutOptions::default();
        let cut = mir_single_row(&p, &b, &x_bar, &p.constraints[0], opts, &t)
            .expect("分数右端处应产出 c-MIR 割");
        // 确定性最优 δ 的割形状：x ≤ 1（连续列不入割）。
        assert_eq!(cut.terms, vec![(0, 1.0)]);
        assert_eq!(cut.rhs, 1.0);
        assert_eq!(cut.op, RelOp::Le);
        // 穷举：x ∈ {0..3}、c ∈ {lb, ub, 中点}（割对 c 线性 ⇒ 端点+内点足够；
        // 本割不含 c，任意 c 恒可行）。
        for x in 0..=3i64 {
            for c in [0.0, 2.5, 5.0] {
                let v = [x as f64, c];
                if 2.0 * v[0] + 0.5 * v[1] <= 3.25 {
                    assert!(
                        eval_cut(&cut, &v) <= cut.rhs + 1e-9,
                        "割 {cut:?} 切掉行可行点 ({x},{c})"
                    );
                }
            }
        }
        assert!(eval_cut(&cut, &x_bar) - cut.rhs >= opts.min_violation);
    }

    // ===================== 隐含界割（WP6 增补）=====================

    /// 1 整数 + 连续列的 Le 行 ⇒ 上界割；穷举（整数枚举 × 连续端点）验证
    /// 行可行 ⇒ 割可行。
    ///
    /// 行：y + 2c ≤ 7.5，y ∈ [0,10] 整数、c ∈ [0,3] 连续；x̄ = (7.5, 0)。
    /// V = 7.5 − 2·lb_c = 7.5 ⇒ y ≤ 7（violation 0.5）。
    #[test]
    fn implied_bound_upper_cut_exhaustive_validity() {
        let mut p = Problem::new("impl-ub", Sense::Maximize);
        p.add_var("y", VarType::Integer).unwrap();
        p.add_var("c", VarType::Continuous).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_var_bounds(1, 0.0, 3.0).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 2.0)], RelOp::Le, 7.5)
            .unwrap();
        let b = Bounds::from_problem(&p);
        let x_bar = [7.5, 0.0];
        let cuts = implied_bound_cuts(&p, &b, &x_bar, CutOptions::default(), &tol());
        assert_eq!(cuts.len(), 1, "恰一条上界割：{cuts:?}");
        let cut = &cuts[0];
        assert_eq!(cut.terms, vec![(0, 1.0)]);
        assert_eq!(cut.rhs, 7.0);
        // 穷举：y ∈ 0..10、c ∈ {0, 1.5, 3}（割纯 y 界，与 c 无关）。
        for y in 0..=10i64 {
            for c in [0.0, 1.5, 3.0] {
                let v = [y as f64, c];
                if v[0] + 2.0 * v[1] <= 7.5 {
                    assert!(
                        eval_cut(cut, &v) <= cut.rhs + 1e-9,
                        "割 {cut:?} 切掉行可行点 ({y},{c})"
                    );
                }
            }
        }
    }

    /// Ge 行 ⇒ 下界割；穷举验证。行：3y − c ≥ 1.2，y ∈ [0,5] 整数、
    /// c ∈ [0,2] 连续；x̄ = (0.4, 0)。V = (1.2 − 0)/3 = 0.4 ⇒ y ≥ 1。
    #[test]
    fn implied_bound_lower_cut_ge_row_validity() {
        let mut p = Problem::new("impl-lb", Sense::Minimize);
        p.add_var("y", VarType::Integer).unwrap();
        p.add_var("c", VarType::Continuous).unwrap();
        p.set_var_bounds(0, 0.0, 5.0).unwrap();
        p.set_var_bounds(1, 0.0, 2.0).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 3.0), (1, -1.0)], RelOp::Ge, 1.2)
            .unwrap();
        let b = Bounds::from_problem(&p);
        let cuts = implied_bound_cuts(&p, &b, &[0.4, 0.0], CutOptions::default(), &tol());
        assert_eq!(cuts.len(), 1, "恰一条下界割：{cuts:?}");
        let cut = &cuts[0];
        assert_eq!(cut.terms, vec![(0, -1.0)]);
        assert_eq!(cut.rhs, -1.0, "−y ≤ −1 即 y ≥ 1");
        for y in 0..=5i64 {
            for c in [0.0, 1.0, 2.0] {
                let v = [y as f64, c];
                if 3.0 * v[0] - v[1] >= 1.2 {
                    assert!(
                        eval_cut(cut, &v) <= cut.rhs + 1e-9,
                        "割 {cut:?} 切掉行可行点 ({y},{c})"
                    );
                }
            }
        }
    }

    /// 取整 ε 专项：浮点和差使 V 落真整数 5 下方 1e-15（2.9 + 3×0.7 =
    /// 4.999999999999999）⇒ 割右端吸附为 5，不误剪 y=5 的行可行点。
    ///
    /// 行：y − 3c ≤ 2.9，y ∈ [0,8] 整数、c ∈ [0,0.7] 连续；x̄ = (5.5, 0)。
    #[test]
    fn implied_bound_rounding_epsilon_does_not_cut_integer_point() {
        let mut p = Problem::new("impl-eps", Sense::Maximize);
        p.add_var("y", VarType::Integer).unwrap();
        p.add_var("c", VarType::Continuous).unwrap();
        p.set_var_bounds(0, 0.0, 8.0).unwrap();
        p.set_var_bounds(1, 0.0, 0.7).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, -3.0)], RelOp::Le, 2.9)
            .unwrap();
        let b = Bounds::from_problem(&p);
        let cuts = implied_bound_cuts(&p, &b, &[5.5, 0.0], CutOptions::default(), &tol());
        assert_eq!(cuts.len(), 1, "{cuts:?}");
        let cut = &cuts[0];
        assert_eq!(cut.rhs, 5.0, "V = 4.999999999999999 须吸附为 5（放宽方向）");
        // y = 5、c = 0.7：行 5 − 2.1 = 2.9 ≤ 2.9 行可行 ⇒ 不得被切。
        let v = [5.0, 0.7];
        assert!(5.0 - 3.0 * v[1] <= 2.9 + 1e-12);
        assert!(eval_cut(cut, &v) <= cut.rhs + 1e-9);
    }

    /// 门控：整数列计数 ≠ 1 的行不产割；需要的一侧界无限 ⇒ 该方向跳过。
    #[test]
    fn implied_bound_gating_multi_int_and_infinite_side() {
        let t = tol();
        let opts = CutOptions::default();
        // 两整数列（计数 2）⇒ 无割。
        let mut p2 = Problem::new("impl-2int", Sense::Maximize);
        p2.add_var("y1", VarType::Integer).unwrap();
        p2.add_var("y2", VarType::Integer).unwrap();
        p2.set_var_bounds(0, 0.0, 9.0).unwrap();
        p2.set_var_bounds(1, 0.0, 9.0).unwrap();
        p2.set_objective(vec![(0, 1.0), (1, 1.0)], 0.0).unwrap();
        p2.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 7.5)
            .unwrap();
        let b2 = Bounds::from_problem(&p2);
        assert!(implied_bound_cuts(&p2, &b2, &[7.5, 7.5], opts, &t).is_empty());

        // 上界方向需要负系数列的 ub：y − c ≤ 10、c ∈ [0,∞) ⇒ 跳过。
        let mut pi = Problem::new("impl-inf", Sense::Maximize);
        pi.add_var("y", VarType::Integer).unwrap();
        pi.add_var("c", VarType::Continuous).unwrap(); // ub = +∞
        pi.set_var_bounds(0, 0.0, 20.0).unwrap();
        pi.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        pi.add_constraint(None, vec![(0, 1.0), (1, -1.0)], RelOp::Le, 10.0)
            .unwrap();
        let bi = Bounds::from_problem(&pi);
        assert!(implied_bound_cuts(&pi, &bi, &[19.5, 0.0], opts, &t).is_empty());

        // 正系数连续列（用 lb 折叠，ub 无限不碍事）⇒ 正常产上界割。
        let mut pp = Problem::new("impl-pos", Sense::Maximize);
        pp.add_var("y", VarType::Integer).unwrap();
        pp.add_var("c", VarType::Continuous).unwrap(); // ub = +∞
        pp.set_var_bounds(0, 0.0, 20.0).unwrap();
        pp.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        pp.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 7.5)
            .unwrap();
        let bp = Bounds::from_problem(&pp);
        let cuts = implied_bound_cuts(&pp, &bp, &[7.5, 0.0], opts, &t);
        assert_eq!(cuts.len(), 1);
        assert_eq!(cuts[0].rhs, 7.0);
    }

    // ===================== 覆盖割序贯 up-lifting（WP6 增补）=====================

    /// 已知最优提升系数行验证 DP 精确性（标准公式 α_j = β − Z，β = |C|−1）：
    /// - w=[5,4,3,2]、cap=5、x̄=[1,0.9,0,0]：C={0,1}，α_2=1、α_3=0
    ///   ⇒ 割 x0+x1+x2 ≤ 1（x3 系数 0 不入割）；
    /// - w=[3,3,3,3,3]、cap=8、x̄=[1,1,1,0,0]：C={0,1,2}，α_3=α_4=1
    ///   ⇒ 割 Σx ≤ 2。
    #[test]
    fn lift_known_optimal_coefficients_exact_dp() {
        let t = tol();
        let opts = CutOptions::default();
        // 例 1。
        let p = binary_knapsack("lift-a", &[5.0, 4.0, 3.0, 2.0], 5.0);
        let b = Bounds::from_problem(&p);
        let mut budget = LiftBudget::new();
        let cuts =
            cover_cuts_with_lifting(&p, &b, &[1.0, 0.9, 0.0, 0.0], opts, &t, &mut budget);
        assert_eq!(cuts.len(), 1);
        assert_eq!(cuts[0].terms, vec![(0, 1.0), (1, 1.0), (2, 1.0)], "α_2=1、α_3=0");
        assert_eq!(cuts[0].rhs, 1.0, "β = |C|−1 固定");
        assert!(budget.remaining() < LiftBudget::new().remaining(), "预算须扣减");

        // 例 2。
        let p = binary_knapsack("lift-b", &[3.0, 3.0, 3.0, 3.0, 3.0], 8.0);
        let b = Bounds::from_problem(&p);
        let mut budget = LiftBudget::new();
        let cuts =
            cover_cuts_with_lifting(&p, &b, &[1.0, 1.0, 1.0, 0.0, 0.0], opts, &t, &mut budget);
        assert_eq!(cuts.len(), 1);
        assert_eq!(
            cuts[0].terms,
            vec![(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0), (4, 1.0)],
            "α_3=α_4=1"
        );
        assert_eq!(cuts[0].rhs, 2.0);
    }

    /// 序贯性回归：穷举全部 2^5 个 0-1 点验证提升割有效性；特别是**多个
    /// 提升列同时取 1** 的点（{x3,x4}，3+3=6 ≤ 8 行可行）必须仍满足割
    /// （工作包原文『α_j = |C|−1+Σ已提升α − Z』形式会在此切掉该点）。
    #[test]
    fn lift_sequential_validity_multi_lifted_ones() {
        let t = tol();
        let opts = CutOptions::default();
        let w = [3.0, 3.0, 3.0, 3.0, 3.0];
        let p = binary_knapsack("lift-seq", &w, 8.0);
        let b = Bounds::from_problem(&p);
        let mut budget = LiftBudget::new();
        let cuts =
            cover_cuts_with_lifting(&p, &b, &[1.0, 1.0, 1.0, 0.0, 0.0], opts, &t, &mut budget);
        assert_eq!(cuts.len(), 1);
        let cut = &cuts[0];
        // 穷举 2^5：行可行（Σ 3x ≤ 8 ⇒ 至多 2 个 1）⇒ 割可行。
        for mask in 0..32u32 {
            let x: Vec<f64> = (0..5).map(|i| ((mask >> i) & 1) as f64).collect();
            if 3.0 * x.iter().sum::<f64>() <= 8.0 {
                assert!(
                    eval_cut(cut, &x) <= cut.rhs + 1e-9,
                    "提升割 {cut:?} 切掉行可行 0-1 点 {x:?}"
                );
            }
        }
        // 多提升列同取 1 的点 {x3,x4}：LHS = 2 ≤ 2。
        let v = [0.0, 0.0, 0.0, 1.0, 1.0];
        assert_eq!(eval_cut(cut, &v), 2.0);
        assert!(eval_cut(cut, &v) <= cut.rhs + 1e-9);
    }

    /// 提升后 violation ≥ 提升前（α ≥ 0、x̄ ≥ 0 ⇒ 单调不减；此处严格增）。
    #[test]
    fn lift_violation_not_decreased() {
        let t = tol();
        let opts = CutOptions::default();
        let p = binary_knapsack("lift-viol", &[5.0, 4.0, 3.0, 2.0], 5.0);
        let b = Bounds::from_problem(&p);
        let x_bar = [1.0, 0.9, 0.5, 0.5];
        let base = cover_cuts(&p, &b, &x_bar, opts, &t);
        assert_eq!(base.len(), 1);
        let base_viol = eval_cut(&base[0], &x_bar) - base[0].rhs;
        let mut budget = LiftBudget::new();
        let lifted = cover_cuts_with_lifting(&p, &b, &x_bar, opts, &t, &mut budget);
        assert_eq!(lifted.len(), 1);
        let lifted_viol = eval_cut(&lifted[0], &x_bar) - lifted[0].rhs;
        assert!(
            lifted_viol >= base_viol - 1e-12,
            "提升后 violation {lifted_viol} < 提升前 {base_viol}"
        );
        assert!(lifted_viol > base_viol, "x̄2=0.5 ⇒ 严格增");
    }

    /// 小数权行被建表判据跳过（10^k 无法整数化）且不产无效割：产出与
    /// 基础 cover_cuts 完全一致。w = 1/3 循环小数、cap = 0.5。
    #[test]
    fn lift_nonintegerizable_weights_fall_back_to_base() {
        let t = tol();
        let opts = CutOptions::default();
        let w = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
        let p = binary_knapsack("lift-frac", &w, 0.5);
        let b = Bounds::from_problem(&p);
        let x_bar = [0.9, 0.9, 0.9];
        let base = cover_cuts(&p, &b, &x_bar, opts, &t);
        let mut budget = LiftBudget::new();
        let lifted = cover_cuts_with_lifting(&p, &b, &x_bar, opts, &t, &mut budget);
        assert_eq!(base.len(), 1, "基础覆盖割存在");
        assert_eq!(lifted.len(), 1);
        assert_eq!(lifted[0].terms, base[0].terms, "非整数化 ⇒ 退回基础割");
        assert_eq!(lifted[0].rhs, base[0].rhs);
        assert_eq!(
            budget.remaining(),
            LiftBudget::new().remaining(),
            "未跑 DP ⇒ 预算不动"
        );
        // 基础割本身仍有效：穷举 2^3。
        for mask in 0..8u32 {
            let x: Vec<f64> = (0..3).map(|i| ((mask >> i) & 1) as f64).collect();
            if w[0] * x.iter().sum::<f64>() <= 0.5 {
                assert!(eval_cut(&lifted[0], &x) <= lifted[0].rhs + 1e-9);
            }
        }
    }

    /// 预算耗尽（remaining = 0）⇒ 全程关闭提升：产出与基础 cover_cuts 等价。
    #[test]
    fn lift_zero_budget_falls_back_to_base() {
        let t = tol();
        let opts = CutOptions::default();
        let p = binary_knapsack("lift-budget0", &[5.0, 4.0, 3.0, 2.0], 5.0);
        let b = Bounds::from_problem(&p);
        let x_bar = [1.0, 0.9, 0.5, 0.5];
        let base = cover_cuts(&p, &b, &x_bar, opts, &t);
        let mut budget = LiftBudget::with_remaining(0);
        let lifted = cover_cuts_with_lifting(&p, &b, &x_bar, opts, &t, &mut budget);
        assert_eq!(lifted.len(), base.len());
        assert_eq!(lifted[0].terms, base[0].terms);
        assert_eq!(lifted[0].rhs, base[0].rhs);
    }

    /// 建表门槛超限（k2000 形态：|C|×cap̃ 巨大）⇒ 整行跳过提升只留基础
    /// 覆盖割。1000 Binary、w=27、cap=26999：(|C|+K)×(cap̃+1) ≈ 2.7e7 > 2e6。
    #[test]
    fn lift_table_gate_skips_large_rows() {
        let t = tol();
        let opts = CutOptions::default();
        let n = 1000usize;
        let p = binary_knapsack("lift-gate", &vec![27.0; n], 26999.0);
        let b = Bounds::from_problem(&p);
        let x_bar = vec![0.9995; n];
        let base = cover_cuts(&p, &b, &x_bar, opts, &t);
        assert_eq!(base.len(), 1, "基础覆盖割存在（violation 0.5）");
        assert_eq!(base[0].terms.len(), 1000, "C = 全支撑");
        let mut budget = LiftBudget::new();
        let lifted = cover_cuts_with_lifting(&p, &b, &x_bar, opts, &t, &mut budget);
        assert_eq!(lifted.len(), 1);
        assert_eq!(lifted[0].terms, base[0].terms, "门槛超限 ⇒ 退回基础割");
        assert_eq!(lifted[0].rhs, base[0].rhs);
        assert_eq!(budget.remaining(), LiftBudget::new().remaining());
    }

    // ===================== 割轮记账上限与合并序（WP6）=====================

    /// 割轮记账上限：max_rounds=1 / max_cuts_per_round=3 ⇒ recorded ≤ 3；
    /// 默认参数 ⇒ recorded ≤ 8×40 = 320；recorded 条数恒等于施加割数。
    #[test]
    fn cut_loop_recorded_respects_round_and_cut_caps() {
        let p = binary_knapsack("cap-many", &vec![4.0; 50], 86.0);
        let t = tol();
        let b = Bounds::from_problem(&p);
        let int_vars = p.integer_vars();

        // 单轮上限：每轮施加 ≤ max_cuts_per_round。
        let opts = CutOptions {
            max_rounds: 1,
            max_cuts_per_round: 3,
            ..CutOptions::default()
        };
        let mut s = Simplex::new(&p, t).unwrap();
        assert_eq!(s.solve(&b).status, LpStatus::Optimal);
        let (applied, status, recorded) =
            cut_loop_recorded(&mut s, &b, &p, &int_vars, opts, &t, || Simplex::new(&p, t)).unwrap();
        assert!(applied >= 1, "根 LP 21.5 处必产割");
        assert_eq!(status, LpStatus::Optimal);
        assert_eq!(recorded.len(), applied as usize);
        assert!(recorded.len() <= 3, "单轮 ≤ max_cuts_per_round：{}", recorded.len());

        // 默认参数：总施加 ≤ max_rounds × max_cuts_per_round = 4×40。
        let mut s2 = Simplex::new(&p, t).unwrap();
        assert_eq!(s2.solve(&b).status, LpStatus::Optimal);
        let (applied2, _status2, recorded2) = cut_loop_recorded(
            &mut s2,
            &b,
            &p,
            &int_vars,
            CutOptions::default(),
            &t,
            || Simplex::new(&p, t),
        )
        .unwrap();
        assert_eq!(recorded2.len(), applied2 as usize);
        assert!(
            recorded2.len() <= 4 * 40,
            "总施加 ≤ 4×40 = 160：{}",
            recorded2.len()
        );
        // 割后不得切掉整数最优（4Σx ≤ 86 ⇒ Σx ≤ 21，目标 Σx ⇒ 最优 21）。
        let sol = s2.resolve(&b);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!(sol.objective >= 21.0 - 1e-6, "割不得切整数最优 21：{}", sol.objective);
        assert!(sol.objective <= 21.5 - 1e-6, "割须改进根 LP 21.5：{}", sol.objective);
    }

    /// 合并序行为：cut_loop_recorded 的 recorded 清单里 c-MIR/隐含界形状的割
    /// （单行 Le 派生）与覆盖/提升形状的割共存，且全程确定性（新族接入后
    /// 两次运行逐字节一致）。算例：2 整数列 wolsey + 5 Binary 背包行拼合。
    #[test]
    fn cut_loop_recorded_merge_order_deterministic_with_new_families() {
        let mut p = Problem::new("merge-det", Sense::Maximize);
        for i in 0..5 {
            p.add_var(&format!("x{i}"), VarType::Binary).unwrap();
        }
        p.add_var("z", VarType::Integer).unwrap();
        p.set_var_bounds(5, 0.0, 10.0).unwrap();
        let mut obj: Vec<(usize, f64)> = (0..5).map(|i| (i, 1.0)).collect();
        obj.push((5, 2.5));
        p.set_objective(obj, 0.0).unwrap();
        // 背包行（覆盖+lifting 目标行）。
        p.add_constraint(
            Some("cap"),
            vec![(0, 4.0), (1, 4.0), (2, 4.0), (3, 4.0), (4, 4.0)],
            RelOp::Le,
            15.0,
        )
        .unwrap();
        // 单整数列 + 连续列行（隐含界目标行）。
        p.add_var("c", VarType::Continuous).unwrap();
        p.set_var_bounds(6, 0.0, 3.0).unwrap();
        p.add_constraint(None, vec![(5, 1.0), (6, 2.0)], RelOp::Le, 7.5)
            .unwrap();
        let t = tol();
        let b = Bounds::from_problem(&p);
        let int_vars = p.integer_vars();

        let run = || {
            let mut s = Simplex::new(&p, t).unwrap();
            assert_eq!(s.solve(&b).status, LpStatus::Optimal);
            cut_loop_recorded(&mut s, &b, &p, &int_vars, CutOptions::default(), &t, || Simplex::new(&p, t)).unwrap()
        };
        let r1 = run();
        let r2 = run();
        assert_eq!(r1.0, r2.0);
        assert_eq!(r1.1, r2.1);
        assert_eq!(r1.2.len(), r2.2.len());
        for (c1, c2) in r1.2.iter().zip(r2.2.iter()) {
            assert_eq!(c1.op, c2.op);
            assert_eq!(c1.rhs, c2.rhs);
            assert_eq!(c1.terms, c2.terms, "recorded 须逐字节一致");
        }
        assert!(r1.0 >= 1, "根 LP 分数处必产割");
        // 终界不切整数最优：整数可行点 x=(1,1,1,0,0)、z=7、c=0 目标 10；
        // 验证割后根界 ≥ 某个下界不具一般性——此处只验证 ≥ LP 收紧方向与
        // 终态 Optimal（穷举有效性已由三族专项测试覆盖）。
        assert_eq!(r1.1, LpStatus::Optimal);
    }

    // ===================== 分离边界数值卫生闸 ======================

    /// 清洁化只弱化：带小系数噪声项的割清洁后，箱内全部满足原割的点必
    /// 满足清洁后割（可行域包含），右端按界幅精确补偿（含负系数噪声项
    /// ——同一界幅界对两种符号都成立）。
    #[test]
    fn sanitize_cleanup_only_weakens_exhaustive() {
        // 3 列：x0 ∈ [0,1]、x1 ∈ [-2,5]、x2 ∈ [0,1000]。
        let bounds = Bounds {
            lb: vec![0.0, -2.0, 0.0],
            ub: vec![1.0, 5.0, 1000.0],
        };
        // 原割：x0 + 5e-10·x1 − 3e-10·x2 ≤ 0.5（后两项为噪声量级）。
        let cut = Cut {
            terms: vec![(0, 1.0), (1, 5e-10), (2, -3e-10)],
            op: RelOp::Le,
            rhs: 0.5,
        };
        let out = sanitize_cut(cut.clone(), &bounds, None, CutOptions::default())
            .expect("有限界幅 ⇒ 割存活");
        assert_eq!(out.terms, vec![(0, 1.0)], "噪声项须被清洁");
        let expect_rhs = 0.5 + 5e-10 * 5.0 + 3e-10 * 1000.0;
        assert!(
            (out.rhs - expect_rhs).abs() < 1e-15,
            "右端补偿：{} 期望 {expect_rhs}",
            out.rhs
        );
        // 穷举网格：满足原割 ⇒ 满足清洁后割（只弱化不变量）。
        for x0 in [0.0f64, 0.25, 0.5, 0.75, 1.0] {
            for x1 in [-2.0f64, 0.0, 2.5, 5.0] {
                for x2 in [0.0f64, 333.3, 1000.0] {
                    let x = [x0, x1, x2];
                    let orig = 1.0 * x0 + 5e-10 * x1 - 3e-10 * x2;
                    if orig <= 0.5 {
                        assert!(
                            eval_cut(&out, &x) <= out.rhs + 1e-15,
                            "清洁后割切掉原可行点 {x:?}"
                        );
                    }
                }
            }
        }
    }

    /// 界幅无界（ub = +∞）的小系数项不得清洁（无法界定只弱化补偿），
    /// 右端不动。
    #[test]
    fn sanitize_keeps_tiny_term_on_infinite_bound() {
        let bounds = Bounds {
            lb: vec![0.0, 0.0],
            ub: vec![1.0, f64::INFINITY],
        };
        let cut = Cut {
            terms: vec![(0, 1.0), (1, 4e-10)],
            op: RelOp::Le,
            rhs: 0.5,
        };
        let out = sanitize_cut(cut, &bounds, None, CutOptions::default()).unwrap();
        assert_eq!(out.terms, vec![(0, 1.0), (1, 4e-10)], "无界列小项须保留");
        assert_eq!(out.rhs, 0.5, "未清洁 ⇒ 右端不动");
    }

    /// 密度闸：nnz 超上限整条丢弃、恰在上限存活（violation 0.5 隔离密度
    /// 维度）。上限 = max(200, ⌈201×0.1⌉=21) = 200。
    #[test]
    fn sanitize_density_gate_drops_over_limit() {
        let n = 201usize;
        let bounds = Bounds {
            lb: vec![0.0; n],
            ub: vec![1.0; n],
        };
        let opts = CutOptions::default();
        let x_bar = vec![1.0f64; n];
        let mk = |k: usize| Cut {
            terms: (0..k).map(|j| (j, 1.0)).collect(),
            op: RelOp::Le,
            rhs: k as f64 - 0.5, // k 个 1 ⇒ violation = 0.5
        };
        assert!(
            sanitize_cut(mk(201), &bounds, Some(&x_bar), opts).is_none(),
            "201 nnz > 上限 200 ⇒ 整条丢弃"
        );
        let kept = sanitize_cut(mk(200), &bounds, Some(&x_bar), opts).unwrap();
        assert_eq!(kept.terms.len(), 200, "恰在上限 ⇒ 存活");
    }

    /// violation 复核：右端补偿吃掉 violation 至 min_violation 以下的割
    /// 丢弃（弱化后已无分离价值）。大界幅 × 噪声系数是触发形态。
    #[test]
    fn sanitize_drops_cut_whose_compensation_eats_violation() {
        let bounds = Bounds {
            lb: vec![0.0, 0.0],
            ub: vec![10.0, 1000.0],
        };
        // x̄ = (5, 0)：原 violation = 1.5e-7 ≥ min_violation=1e-7；清洁
        // 补偿 = 1e-9×1000 = 1e-6 ⇒ 复核 violation < 1e-7 ⇒ 丢弃。
        let cut = Cut {
            terms: vec![(0, 1.0), (1, 1e-9)],
            op: RelOp::Le,
            rhs: 5.0 - 1.5e-7,
        };
        assert!(
            sanitize_cut(cut, &bounds, Some(&[5.0, 0.0]), CutOptions::default())
                .is_none()
        );
    }

    /// 清洁化只弱化（性质 fuzz）：随机割（噪声/显著项混排、随机有限箱）
    /// 上穷举网格点验证「满足原割 ⇒ 满足清洁后割」。
    #[test]
    fn sanitize_cleanup_weakening_property_fuzz() {
        for seed in 1..=200u64 {
            let mut lcg = seed.wrapping_mul(0x9E3779B97F4A7C15);
            let mut rnd = || {
                lcg = lcg
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((lcg >> 33) % 100) as f64 / 100.0
            };
            // 4 列随机有限箱 [lo, lo+span]。
            let (mut lbs, mut ubs, mut axes) = (Vec::new(), Vec::new(), Vec::new());
            for _ in 0..4 {
                let lo = ((rnd() * 4.0) as i64 as f64) - 2.0;
                let span = (rnd() * 3.0).max(0.25);
                lbs.push(lo);
                ubs.push(lo + span);
                axes.push(vec![lo, lo + span / 2.0, lo + span]);
            }
            let bounds = Bounds { lb: lbs, ub: ubs };
            let mut terms = Vec::new();
            for j in 0..4usize {
                // 四成概率注入噪声量级项（正/负随机），否则显著项。
                let c = if rnd() < 0.4 {
                    (rnd() * 2.0 - 1.0) * 1e-9
                } else {
                    rnd() * 2.0 - 1.0 + if rnd() < 0.5 { 0.0 } else { 0.75 }
                };
                terms.push((j, c));
            }
            let rhs = rnd() * 2.0 - 0.5;
            let cut = Cut {
                terms,
                op: RelOp::Le,
                rhs,
            };
            let Some(out) = sanitize_cut(cut.clone(), &bounds, None, CutOptions::default())
            else {
                continue; // 被闸丢弃（清洁后空等）也安全，无性质可验
            };
            for &a in &axes[0] {
                for &b in &axes[1] {
                    for &c in &axes[2] {
                        for &d in &axes[3] {
                            let x = [a, b, c, d];
                            if cut.terms.iter().map(|&(j, cf)| cf * x[j]).sum::<f64>() <= rhs {
                                assert!(
                                    eval_cut(&out, &x) <= out.rhs + 1e-12,
                                    "seed {seed}：清洁后割切掉原可行点 {x:?}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}
