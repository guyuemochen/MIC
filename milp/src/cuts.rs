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

use crate::error::{MilpError, Result};
use crate::model::{Problem, RelOp, Tolerances};
use crate::simplex::{Bounds, LpStatus, RowOrigin, Simplex};

/// 割系数动态范围上限：显著项最大/最小绝对值之比超过即整条丢弃（契约值）。
const MAX_COEF_RATIO: f64 = 1e9;

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
    /// 最大割轮数。默认 5。
    pub max_rounds: usize,
    /// 每轮最多施加割数。默认 20。
    pub max_cuts_per_round: usize,
    /// f0 距 0/1 的最小距离。默认 0.005。
    pub away: f64,
    /// 割须被当前 LP 解违反的最小量。默认 1e-7。
    pub min_violation: f64,
}

impl Default for CutOptions {
    /// 默认值：max_rounds=5、max_cuts_per_round=20、away=0.005、min_violation=1e-7。
    fn default() -> Self {
        CutOptions {
            max_rounds: 5,
            max_cuts_per_round: 20,
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
/// 非基割松弛列因无记账（本函数不携带 cut_loop 的私有账本）而整列丢弃——
/// π ≥ 0、d ≥ 0，丢弃只会放宽割，安全；cut_loop 内部走带记账路径，可精确
/// 回代。调用前提：`bounds` 即 Simplex 最近一次 solve/resolve 所用界，且状态
/// 为 Optimal。
pub fn gomory_cuts(
    s: &Simplex,
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
    s: &Simplex,
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
    let num_slack = layout.slack.len();
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
    let mut x_bar = vec![f64::NAN; n];
    for (row, &bcol) in basis.iter().enumerate() {
        if bcol < n {
            let (_, value) = s.tableau_row(row);
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
        return Vec::new(); // 表格与界不一致（防御）
    }

    let mut cuts = Vec::new();
    for row in 0..s.num_rows() {
        let bcol = basis[row];
        if bcol >= num_struct || !is_int[bcol] {
            continue; // 基变量须为整数类型结构列
        }
        let (table, value) = s.tableau_row(row);
        if table.len() < num_struct + num_slack || !value.is_finite() {
            continue;
        }
        if table.iter().any(|v| !v.is_finite()) {
            continue; // 数值异常行整行跳过
        }
        let f0 = value - value.floor();
        if !(f0 >= opts.away && f0 <= 1.0 - opts.away) {
            continue;
        }

        // Σ_j π_j·d_j ≥ 1 的结构变量系数与常数项。
        let mut coef = vec![0.0; n];
        let mut const_term = 0.0f64;
        let mut any = false;
        for nb in &nonbasic {
            let j = nb.col;
            let alpha = if j < n {
                if nb.at_upper {
                    -table[j]
                } else {
                    table[j]
                }
            } else {
                table[j] // 松弛/剩余列恒贴下界 0
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
                continue; // π = 0 无贡献；非有限为数值异常，丢列（放宽方向）
            }
            if pi < tol.zero {
                continue; // 近零项丢弃（π ≥ 0、d ≥ 0，丢弃只放宽割）
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
                            continue; // 防御：行来源越界，丢列（放宽）
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
                            continue; // 无记账（公开入口路径）：丢列（放宽）
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

/// 根节点割循环：生成→记账→add_row→resolve，直到无割/超 max_rounds/LP 非
/// Optimal；返回 (施加割数, 最终 LP 状态)。
///
/// 入口先做一次防御性 `resolve(bounds)`，保证割生成看到的表格确为当前界下的
/// Optimal 状态。割循环中 LP 变 Infeasible 理论不可能（割不切整数可行解），
/// 若发生返回 `Err(MilpError::Numerical)` 而非报 Infeasible；Unbounded/
/// IterationLimit 则如实作为终态返回。
pub fn cut_loop(
    s: &mut Simplex,
    bounds: &Bounds,
    problem: &Problem,
    int_vars: &[usize],
    opts: CutOptions,
    tol: &Tolerances,
) -> Result<(u64, LpStatus)> {
    cut_loop_recorded(s, bounds, problem, int_vars, opts, tol).map(|(applied, status, _)| {
        (applied, status)
    })
}

/// 与 [`cut_loop`] 完全相同的割循环，额外把施加的割按施加顺序原样返回
/// （terms/rhs 即 add_row 时传入的原尺度值）。行为与 cut_loop 逐字节一致。
///
/// 用途：并行引擎在各 worker 的 Simplex 副本上按相同顺序 `add_row` 这些割
/// 即可复现根割表格（行序一致 ⇒ 表格一致）。[`cut_loop`] 即本函数丢弃割
/// 清单后的委托；既有语义逐项保持：入口防御性 `resolve(bounds)`、
/// `max_cuts_per_round` 截断、Infeasible → `Err(Numerical)`。
pub fn cut_loop_recorded(
    s: &mut Simplex,
    bounds: &Bounds,
    problem: &Problem,
    int_vars: &[usize],
    opts: CutOptions,
    tol: &Tolerances,
) -> Result<(u64, LpStatus, Vec<Cut>)> {
    let mut applied: u64 = 0;
    let mut ledger: Vec<CutRecord> = Vec::new();
    let mut recorded: Vec<Cut> = Vec::new();
    // 防御性同步：确保表格与 bounds 一致且为 Optimal。
    let mut status = s.resolve(bounds).status;
    if status != LpStatus::Optimal {
        return Ok((0, status, recorded));
    }
    for _round in 0..opts.max_rounds {
        let mut batch = gomory_cuts_inner(s, problem, bounds, int_vars, opts, tol, &ledger);
        if batch.len() > opts.max_cuts_per_round {
            batch.truncate(opts.max_cuts_per_round);
        }
        if batch.is_empty() {
            break;
        }
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
        status = s.resolve(bounds).status;
        if status == LpStatus::Infeasible {
            return Err(MilpError::Numerical(
                "割循环后根 LP 不可行：割错误地切掉了全部整数可行解（理论不可能，判定为数值异常）"
                    .to_string(),
            ));
        }
        if status != LpStatus::Optimal {
            break; // Unbounded/IterationLimit：如实返回终态
        }
    }
    Ok((applied, status, recorded))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Sense, Tolerances, VarType};

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

    /// 默认参数逐项精确符合契约。
    #[test]
    fn cut_options_default_values() {
        let o = CutOptions::default();
        assert_eq!(o.max_rounds, 5);
        assert_eq!(o.max_cuts_per_round, 20);
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
            cut_loop(&mut s, &b, &p, &p.integer_vars(), CutOptions::default(), &t).unwrap();
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
            cut_loop(&mut s, &b, &p, &p.integer_vars(), CutOptions::default(), &t).unwrap();
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
            cut_loop(&mut s1, &b, &p, &p.integer_vars(), CutOptions::default(), &t).unwrap();

        let mut s2 = Simplex::new(&p, t).unwrap();
        assert_eq!(s2.solve(&b).status, LpStatus::Optimal);
        let (applied2, status2, cuts2) = cut_loop_recorded(
            &mut s2,
            &b,
            &p,
            &p.integer_vars(),
            CutOptions::default(),
            &t,
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
}
