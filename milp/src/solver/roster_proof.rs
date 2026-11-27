//! 非负容量约束的互斥团割，以及完整整数模型上的有限前沿下界证明。
//! 根点和检验数仅用于挑选割；互斥性由全局行和变量域独立证明。

mod branching;
mod cliques;
pub(super) mod group_cover;
pub(super) mod reduction;

use crate::cuts::Cut;
use crate::heuristics::{round_and_check, Incumbent};
use crate::model::{Problem, RelOp, Sense, Tolerances, VarType};
use crate::simplex::{Bounds, LpStatus};
use crate::simplex_sparse::SparseSimplex;
use std::collections::BTreeMap;
use std::time::Instant;

pub(super) use cliques::separate;

fn usable_tolerances(tol: &Tolerances) -> bool {
    [
        tol.feasibility,
        tol.optimality,
        tol.pivot,
        tol.integrality,
        tol.zero,
    ]
    .into_iter()
    .all(|t| t.is_finite() && t >= 0.0)
}

/// 只采用可靠的正系数上界行。全体列非负，因此任何两项都可与容量单独比较。
/// 门槛留足原行终检、整数容差、单纯形行归一化和浮点加法的余量。
fn guarded_capacity(row: &crate::model::Constraint, p: &Problem, tol: &Tolerances) -> Option<f64> {
    if !matches!(row.op, RelOp::Le | RelOp::Eq) || row.rhs < 0.0 || !row.rhs.is_finite() {
        return None;
    }
    let mut norm = 0.0;
    let mut largest = 0.0_f64;
    for &(j, a) in &row.expr.terms {
        if !(a > 0.0 && a.is_finite() && p.vars[j].lb >= 0.0) {
            return None;
        }
        norm += a;
        largest = largest.max(a);
    }
    if largest <= tol.zero {
        return None;
    }
    let scale = 1.0 + row.rhs + norm;
    let margin = scale * (super::FINAL_VIOLATION_REL + tol.integrality + 8.0 * f64::EPSILON)
        + (10.0 * tol.feasibility + 4.0 * tol.zero) * (1.0 + largest);
    let cap = row.rhs + margin;
    cap.is_finite().then_some(cap)
}

/// 在已证原始根界和调用者提供的可行费用之间，返回可安全合并的证明进展。
pub(super) struct ProofOutcome {
    pub(super) bound: f64,
    pub(super) candidate: Option<Incumbent>,
    pub(super) nodes: u64,
    pub(super) lp_iterations: u64,
}

struct Cell {
    lp: SparseSimplex,
    bounds: Bounds,
    x: Vec<f64>,
    lower: f64,
}

/// 补入由非负整数列构成的冗余合计。只用创建后取得的列号识别，名称不参与分支。
fn add_totals(p: &mut Problem) -> Vec<(usize, f64)> {
    let mut occurrences = vec![0u8; p.vars.len()];
    for row in &p.constraints {
        for &(j, _) in &row.expr.terms {
            occurrences[j] = occurrences[j].saturating_add(1);
        }
    }
    let mut groups = BTreeMap::<u64, Vec<usize>>::new();
    for &(j, cost) in &p.obj.terms {
        let v = &p.vars[j];
        if v.var_type == VarType::Integer
            && v.lb >= 0.0
            && occurrences[j] == 1
            && cost > 0.0
            && cost.is_finite()
            && cost.fract() == 0.0
        {
            groups.entry(cost.to_bits()).or_default().push(j);
        }
    }
    let mut totals = Vec::new();
    for (bits, members) in groups.into_iter().filter(|(_, js)| js.len() >= 2).take(8) {
        let mut name = format!("_mic_sum_{}", p.vars.len());
        while p.var_index(&name).is_some() {
            name.push('_');
        }
        let Ok(j) = p.add_var(&name, VarType::Integer) else {
            break;
        };
        let mut terms: Vec<_> = members.into_iter().map(|i| (i, 1.0)).collect();
        terms.push((j, -1.0));
        p.add_constraint(None, terms, RelOp::Eq, 0.0)
            .expect("有效的非负整数合计等式");
        totals.push((j, f64::from_bits(bits)));
    }
    totals
}

/// 只由列符号证明全行活动非负或非正。避免乘法、求和及容差放宽。
fn implied_by_column_signs(
    row: &crate::model::Constraint,
    vars: &[crate::model::Variable],
) -> bool {
    let nonnegative =
        |j: usize, a: f64| (a > 0. && vars[j].lb >= 0.) || (a < 0. && vars[j].ub <= 0.) || a == 0.;
    let lower = || row.rhs <= 0. && row.expr.terms.iter().all(|&(j, a)| nonnegative(j, a));
    let upper = || row.rhs >= 0. && row.expr.terms.iter().all(|&(j, a)| nonnegative(j, -a));
    match row.op {
        RelOp::Ge => lower(),
        RelOp::Le => upper(),
        RelOp::Eq => lower() && upper(),
    }
}

/// 等式中补一个上下界均为零的无费用列。使用原行最大系数保持单纯形行归一尺度，
/// 近零行原本不入表，新增列也不能让它入表。原工作模型与人工列规则均保持原义。
#[cfg(feature = "sparse_phase1_recovery")]
fn add_proof_zero_columns(p: &mut Problem, tol: &Tolerances) -> usize {
    let first = p.vars.len();
    for r in p.constraints.iter_mut().filter(|r| r.op == RelOp::Eq) {
        let largest = r.expr.terms.iter().map(|t| t.1.abs()).fold(0., f64::max);
        if largest > tol.zero {
            let j = p.vars.len();
            p.vars.push(crate::model::Variable {
                name: String::new(),
                var_type: VarType::Continuous,
                lb: 0.,
                ub: 0.,
            });
            r.expr.terms.push((j, largest));
        }
    }
    first
}

/// 探测用同一绝对时限，根、未采用的备选以及中断子 LP 都消耗调用者的节点预算。
/// 调用者仍须在原始模型复核经还原的 candidate，并按自己的限额与 gap 规则定终态。
#[cfg(test)]
pub(super) fn prove(
    work: &Problem,
    cuts: &[Cut],
    root_bound: f64,
    incumbent: f64,
    gap_percent: f64,
    tol: &Tolerances,
    deadline: Instant,
    node_limit: u64,
) -> ProofOutcome {
    prove_domain(
        work,
        cuts,
        root_bound,
        incumbent,
        gap_percent,
        tol,
        deadline,
        node_limit,
        ProofMode {
            paired: false,
            can_reduce: true,
        },
    )
}

#[derive(Clone, Copy)]
struct ProofMode {
    #[cfg_attr(not(feature = "sparse_phase1_recovery"), allow(dead_code))]
    paired: bool,
    can_reduce: bool,
}

fn prove_domain(
    work: &Problem,
    cuts: &[Cut],
    root_bound: f64,
    incumbent: f64,
    gap_percent: f64,
    tol: &Tolerances,
    deadline: Instant,
    node_limit: u64,
    mode: ProofMode,
) -> ProofOutcome {
    let mut out = ProofOutcome {
        bound: root_bound,
        candidate: None,
        nodes: 0,
        lp_iterations: 0,
    };
    if !usable_tolerances(tol)
        || work.sense != Sense::Minimize
        || !root_bound.is_finite()
        || !incumbent.is_finite()
        || !(gap_percent > 0.0 && gap_percent.is_finite())
        || root_bound >= incumbent - incumbent.abs() * (gap_percent / 100.0)
        || node_limit == 0
        || Instant::now() >= deadline
    {
        return out;
    }
    let mut p = work.clone();
    for cut in cuts {
        if p.add_constraint(None, cut.terms.clone(), cut.op, cut.rhs)
            .is_err()
        {
            return out;
        }
    }
    let totals = add_totals(&mut p);
    p.constraints
        .retain(|r| !implied_by_column_signs(r, &p.vars));
    if Instant::now() >= deadline {
        return out;
    }
    #[cfg(feature = "sparse_phase1_recovery")]
    // 首次定价压缩可能把宽模型降到数千列，仍有足够多的等式从独占零列初基获益。
    let first_zero = ((p.vars.len() >= 4096
        || (mode.paired && (!mode.can_reduce || p.vars.len() >= 1024)))
        && p.obj.terms.iter().any(|&(_, cost)| cost < 0.))
    .then(|| add_proof_zero_columns(&mut p, tol));
    let Ok(mut lp) = SparseSimplex::new(&p, *tol) else {
        return out;
    };
    lp.set_time_deadline(Some(deadline));
    let bounds = Bounds::from_problem(&p);
    #[cfg(feature = "sparse_phase1_recovery")]
    let sol = match first_zero {
        Some(first) => lp.solve_proof(&bounds, first),
        // 宽等式证明根的列若全为非负费用，原目标在各列下界已对偶可行。
        // 用等式的固定零人工列起步，复用对偶修复及原行、原目标核验。
        #[cfg(feature = "choice_dual_root")]
        None if p.vars.len() >= 4096
            && p.vars.len() / p.constraints.len().max(1) >= 32
            && p.constraints.iter().filter(|r| r.op == RelOp::Eq).count() * 2 >= p.constraints.len()
            && p.obj.terms.iter().all(|&(_, c)| c >= 0.) =>
        {
            crate::simplex_sparse::LpEngine::solve_choice_root(&mut lp, &bounds)
        }
        None => lp.solve(&bounds),
    };
    #[cfg(not(feature = "sparse_phase1_recovery"))]
    let sol = lp.solve(&bounds);
    out.nodes = 1;
    out.lp_iterations = lp.iterations();
    if sol.status != LpStatus::Optimal || !sol.objective.is_finite() {
        return out;
    }
    let original_bounds = Bounds::from_problem(work);
    let candidate = |x: &[f64]| {
        let mut inc = round_and_check(work, &original_bounds, x.get(..work.vars.len())?, tol)?;
        super::final_check(work, &mut inc.values).ok()?;
        if !inc
            .values
            .iter()
            .zip(&work.vars)
            .all(|(&v, col)| col.lb <= v && v <= col.ub)
        {
            return None;
        }
        inc.objective = work.obj.eval(&inc.values) + work.obj_offset;
        inc.objective.is_finite().then_some(inc)
    };
    out.candidate = candidate(&sol.x);
    let mut upper = out
        .candidate
        .as_ref()
        .map_or(incumbent, |inc| inc.objective.min(incumbent));
    let target = |u: f64| u - u.abs() * (gap_percent / 100.0);
    // 完整模型已知存在可行解。若新界越过可行费用很多，保留旧证据。
    let consistent = |bound: f64, u: f64| bound <= u + 10.0 * tol.feasibility * (1.0 + u.abs());
    if !consistent(sol.objective, upper) {
        return out;
    }
    let root = root_bound.max(sol.objective.min(upper));
    out.bound = root;
    if root >= target(upper) || out.nodes >= node_limit || Instant::now() >= deadline {
        return out;
    }

    if mode.can_reduce
        && node_limit.saturating_sub(out.nodes) >= 3
        && work.vars.len() >= 4096
        && !cuts.is_empty()
    {
        // 只需证到用户的目标差距。超过门槛的整数列可在这次证明中删去；
        // 最终下界仍取两个子域各自已有证明的较小者。
        let cutoff = target(upper).next_up().min(upper);
        if let Some(screen) = reduction::screen_columns(&p, &sol, cutoff, tol, deadline) {
            if let Some(sub) = reduction::prove_subset(
                work,
                cuts,
                &screen,
                root,
                upper,
                gap_percent,
                tol,
                deadline,
                node_limit - out.nodes,
                ProofMode {
                    can_reduce: false,
                    ..mode
                },
            ) {
                out.nodes += sub.nodes;
                out.lp_iterations += sub.lp_iterations;
                out.bound = out.bound.max(sub.bound);
                if let Some(inc) = sub.candidate.filter(|inc| inc.objective < upper) {
                    out.candidate = Some(inc);
                }
                return out;
            }
        }
    }
    let Some(branches) = branching::Branching::new(&p, totals, tol, deadline) else {
        return out;
    };
    let initial = Cell {
        lp,
        bounds,
        x: sol.x,
        lower: root,
    };
    let mut frontier = vec![initial];
    let mut sealed = f64::INFINITY;
    let int_vars = work.integer_vars();
    // 一个快照包含稀疏矩阵、LU、工作缓冲及变量界；按问题尺寸限制同时保留的快照数。
    let nnz: usize = p.constraints.iter().map(|r| r.expr.terms.len()).sum();
    let bytes = nnz
        .saturating_mul(48)
        .saturating_add(p.vars.len().saturating_mul(192))
        .saturating_add(p.constraints.len().saturating_pow(2).saturating_mul(32))
        .max(1);
    let max_open = (768 * 1024 * 1024 / bytes).clamp(1, 192);
    while !frontier.is_empty()
        && frontier.len() <= max_open
        && out.nodes < node_limit
        && Instant::now() < deadline
    {
        let next = (0..frontier.len())
            .min_by(|&a, &b| frontier[a].lower.total_cmp(&frontier[b].lower))
            .expect("非空前沿");
        let parent = frontier.swap_remove(next);
        let mut best: Option<Vec<Cell>> = None;
        let mut best_low = parent.lower;
        let mut best_high = parent.lower;
        let options = branches.options(&p, &parent.x, tol, deadline);
        for branch in options {
            if best_low > parent.lower + 0.05
                || out.nodes >= node_limit
                || Instant::now() >= deadline
            {
                break;
            }
            let mut children = Vec::new();
            let mut lower = [f64::INFINITY; 2];
            let mut complete = true;
            for (side, value) in lower.iter_mut().enumerate() {
                if out.nodes >= node_limit || Instant::now() >= deadline {
                    complete = false;
                    break;
                }
                let mut b = parent.bounds.clone();
                branches.restrict(&branch, side, &mut b);
                if b.lb.iter().zip(&b.ub).any(|(l, u)| l > u) {
                    continue;
                }
                let mut child = parent.lp.clone();
                child.set_time_deadline(Some(deadline));
                let before = child.iterations();
                let res = child.resolve(&b);
                out.nodes += 1;
                out.lp_iterations += child.iterations() - before;
                match res.status {
                    LpStatus::Optimal if res.objective.is_finite() => {
                        // 即使另一侧随后超限，这个子域的可行点也属于原问题。先收解，
                        // 只有两个子域都得到证明后，才会在下方更新全局下界。
                        if res.objective < upper
                            && crate::branch_bound::is_integral(&int_vars, &res.x, tol.integrality)
                        {
                            if let Some(inc) = candidate(&res.x).filter(|inc| inc.objective < upper)
                            {
                                upper = inc.objective;
                                out.candidate = Some(inc);
                                out.bound = out.bound.min(upper);
                                if out.bound >= target(upper) {
                                    return out;
                                }
                            }
                        }
                        // 子界向下留余量，并继承父域已有界。
                        *value = parent.lower.max(res.objective - 100.0 * tol.optimality);
                        children.push(Cell {
                            lp: child,
                            bounds: b,
                            x: res.x,
                            lower: *value,
                        });
                    }
                    LpStatus::Infeasible => {}
                    _ => {
                        complete = false;
                        break;
                    }
                }
            }
            if !complete {
                continue;
            }
            let lo = lower[0].min(lower[1]);
            let hi = lower[0].max(lower[1]);
            if best.is_none()
                || lo > best_low + 100.0 * tol.optimality
                || (lo >= best_low && hi > best_high + tol.integrality)
            {
                best = Some(children);
                best_low = lo;
                best_high = hi;
            }
        }
        let Some(children) = best else {
            // 未得完备替换时保留父界，全局前沿不丢任何未证子树。
            break;
        };
        for c in children {
            if c.lower >= target(upper) {
                sealed = sealed.min(c.lower);
            } else {
                frontier.push(c);
            }
        }
        let global = frontier.iter().fold(sealed, |b, node| b.min(node.lower));
        if !global.is_finite() || !consistent(global, upper) {
            break;
        }
        out.bound = out.bound.max(global.min(upper));
        if out.bound >= target(upper) {
            break;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simplex::LpSolution;
    use std::time::Duration;

    fn triangle() -> Problem {
        let mut p = Problem::new("互斥三角形", Sense::Minimize);
        for name in ["a", "b", "c"] {
            p.add_var(name, VarType::Binary).unwrap();
        }
        for (a, b) in [(0, 1), (0, 2), (1, 2)] {
            p.add_constraint(None, vec![(a, 1.), (b, 1.)], RelOp::Le, 1.)
                .unwrap();
        }
        p.set_objective(vec![(0, -1.), (1, -1.), (2, -1.)], 0.)
            .unwrap();
        p
    }

    fn end() -> Instant {
        Instant::now() + Duration::from_secs(3)
    }

    #[cfg(feature = "sparse_phase1_recovery")]
    #[test]
    fn proof_zero_column_preserves_row_scale_and_near_zero_row_handling() {
        let mut p = Problem::new("固定零辅助列", Sense::Minimize);
        p.add_var("x", VarType::Continuous).unwrap();
        let tol = Tolerances::default();
        p.add_constraint(None, vec![(0, -8.)], RelOp::Eq, -8.)
            .unwrap();
        p.add_constraint(None, vec![(0, 0.5 * tol.zero)], RelOp::Eq, 0.)
            .unwrap();
        p.add_constraint(None, vec![(0, 2.)], RelOp::Ge, 0.)
            .unwrap();
        assert_eq!(add_proof_zero_columns(&mut p, &tol), 1);
        assert_eq!(p.vars.len(), 2);
        assert_eq!((p.vars[1].lb, p.vars[1].ub), (0., 0.));
        assert_eq!(p.constraints[0].expr.terms, [(0, -8.), (1, 8.)]);
        assert_eq!(p.constraints[1].expr.terms.len(), 1);
        assert_eq!(p.constraints[2].expr.terms.len(), 1);
        p.validate().unwrap();
    }

    #[test]
    fn sign_implied_rows_require_every_term_and_both_sides_of_an_equality() {
        let mut p = Problem::new("符号可删行", Sense::Minimize);
        for j in 0..3 {
            p.add_var(&format!("x{j}"), VarType::Continuous).unwrap();
        }
        p.vars[1].lb = -1.;
        p.vars[1].ub = 0.;
        p.vars[2].ub = 0.;
        let holds = |terms, op, rhs| {
            implied_by_column_signs(
                &crate::model::Constraint {
                    name: None,
                    expr: crate::model::LinearExpr { terms },
                    op,
                    rhs,
                },
                &p.vars,
            )
        };
        assert!(holds(vec![(0, 2.), (1, -3.)], RelOp::Ge, 0.));
        assert!(holds(vec![(0, -2.), (1, 3.)], RelOp::Le, 0.));
        assert!(holds(vec![(2, 1.)], RelOp::Eq, 0.));
        assert!(!holds(vec![(0, 2.), (1, -3.)], RelOp::Eq, 0.));
        assert!(!holds(vec![(0, 2.), (1, 3.)], RelOp::Ge, 0.));
        assert!(!holds(vec![(0, 1.)], RelOp::Ge, f64::MIN_POSITIVE));
        assert!(!holds(vec![(1, 1.)], RelOp::Le, -f64::MIN_POSITIVE));
    }

    #[test]
    fn clique_is_valid_for_every_integer_point_and_restores_root_candidate() {
        let p = triangle();
        let tol = Tolerances::default();
        let root = LpSolution {
            status: LpStatus::Optimal,
            objective: -1.5,
            x: vec![0.5; 3],
            duals: vec![],
        };
        let cuts = separate(&p, &root, &[0.; 3], &tol, end());
        assert_eq!(cuts.len(), 1);
        assert_eq!(cuts[0].terms, vec![(0, 1.), (1, 1.), (2, 1.)]);
        for mask in 0..8 {
            let x: Vec<_> = (0..3).map(|i| ((mask >> i) & 1) as f64).collect();
            if p.constraints.iter().all(|row| row.expr.eval(&x) <= row.rhs) {
                assert!(cuts[0].terms.iter().map(|&(j, w)| w * x[j]).sum::<f64>() <= cuts[0].rhs);
            }
        }
        let proof = prove(&p, &cuts, -1.5, -1., 1., &tol, end(), 1);
        assert_eq!(proof.nodes, 1);
        assert_eq!(proof.bound, -1.);
        assert_eq!(proof.candidate.unwrap().objective, -1.);
    }

    #[test]
    fn tiny_overflow_negative_and_continuous_capacities_do_not_create_edges() {
        let tol = Tolerances::default();
        let root = LpSolution {
            status: LpStatus::Optimal,
            objective: 0.,
            x: vec![0.6; 3],
            duals: vec![],
        };
        let mut p = triangle();
        for r in &mut p.constraints {
            r.rhs = 2. - 1e-12;
        }
        assert!(separate(&p, &root, &[0.; 3], &tol, end()).is_empty());
        for r in &mut p.constraints {
            r.rhs = 1e-15;
            for (_, w) in &mut r.expr.terms {
                *w = 1e-15;
            }
        }
        assert!(separate(&p, &root, &[0.; 3], &tol, end()).is_empty());
        for r in &mut p.constraints {
            r.rhs = f64::MAX;
            for (_, w) in &mut r.expr.terms {
                *w = f64::MAX;
            }
        }
        assert!(separate(&p, &root, &[0.; 3], &tol, end()).is_empty());
        let mut p = triangle();
        p.vars[0].var_type = VarType::Continuous;
        p.vars[1].lb = -1.;
        assert!(separate(&p, &root, &[0.; 3], &tol, end()).is_empty());
    }

    #[test]
    fn interrupted_pair_and_zero_budget_preserve_parent_bound() {
        let p = triangle();
        let tol = Tolerances::default();
        for (deadline, limit) in [(Instant::now(), 50), (end(), 0), (end(), 1), (end(), 2)] {
            let proof = prove(&p, &[], -1.5, -1., 1., &tol, deadline, limit);
            assert!(proof.nodes <= limit);
            assert_eq!(proof.bound, -1.5); // 仅完成一个孩子时不能报告它的 -1。
        }
        let proof = prove(&p, &[], -1.5, -1., 1., &tol, end(), 20);
        assert!(proof.bound >= -1.000001 && proof.bound <= -1.);
        assert!(proof.nodes >= 3);
    }

    #[test]
    fn first_child_keeps_checked_incumbent_when_sibling_exceeds_node_budget() {
        let mut p = triangle();
        p.obj_offset = 7.;
        let tol = Tolerances::default();
        // 全零解费用 7；LP 根三个 0.5 均不能直接取整。第一个子域就有费用 6 的整数点。
        let proof = prove(&p, &[], 5.5, 7., 1., &tol, end(), 2);
        assert_eq!(proof.nodes, 2);
        assert_eq!(proof.bound, 5.5);
        let mut candidate = proof.candidate.unwrap();
        super::super::final_check(&p, &mut candidate.values).unwrap();
        assert_eq!(candidate.objective, 6.);
        assert_eq!(
            p.obj.eval(&candidate.values) + p.obj_offset,
            candidate.objective
        );

        let complete = prove(&p, &[], 5.5, 7., 1., &tol, end(), 10);
        assert_eq!(complete.candidate.unwrap().objective, 6.);
        assert!(complete.bound >= 5.999999 && complete.bound <= 6.);
    }

    #[cfg(feature = "choice_dual_root")]
    #[test]
    fn wide_nonnegative_proof_with_degenerate_equalities_matches_enumeration() {
        let mut p = Problem::new("宽等式覆盖", Sense::Minimize);
        for name in ["甲", "乙", "丙"] {
            p.add_var(name, VarType::Binary).unwrap();
        }
        let sum = p.add_var("两项之和", VarType::Continuous).unwrap();
        // 费用为零的固定列增加宽度，但不增加枚举的整数维数。
        for j in 4..4096 {
            p.vars.push(crate::model::Variable {
                name: format!("固定列{j}"),
                var_type: VarType::Continuous,
                lb: 0.,
                ub: 0.,
            });
        }
        for (a, b) in [(0, 1), (1, 2), (0, 2)] {
            p.add_constraint(None, vec![(a, 1.), (b, 1.)], RelOp::Ge, 1.).unwrap();
            p.add_constraint(None, vec![(0, 1.), (1, 1.), (sum, -1.)], RelOp::Eq, 0.).unwrap();
        }
        p.set_objective((0..3).map(|j| (j, 1.)).collect(), 0.).unwrap();
        let optimum = (0u32..8)
            .filter(|mask| [(0, 1), (1, 2), (0, 2)].iter()
                .all(|&(a, b)| mask & ((1 << a) | (1 << b)) != 0))
            .map(u32::count_ones).min().unwrap() as f64;
        let tol = Tolerances::default();
        let one = prove(&p, &[], 0., 3., 1., &tol, end(), 1);
        assert_eq!(one.nodes, 1);
        assert!((one.bound - 1.5).abs() < 1e-7);
        let full = prove(&p, &[], 0., 3., 1., &tol, end(), 32);
        assert!(full.bound <= optimum && full.bound >= optimum - 1e-5);
        let mut solution = full.candidate.unwrap();
        assert_eq!(solution.objective, optimum);
        super::super::final_check(&p, &mut solution.values).unwrap();
        let invalid = prove(&p, &[], 0., f64::NAN, 1., &tol, end(), 32);
        assert_eq!((invalid.nodes, invalid.bound), (0, 0.));
    }

    #[test]
    fn sum_columns_are_identified_by_index_and_preserve_feasible_assignments() {
        let mut p = Problem::new("合计名称重叠", Sense::Minimize);
        for name in ["proof_sum_17", "_mic_sum_3", "a"] {
            p.add_var(name, VarType::Integer).unwrap();
        }
        for (j, rhs) in [(0, 2.), (1, 3.), (2, 4.)] {
            p.add_constraint(None, vec![(j, 1.)], RelOp::Le, rhs)
                .unwrap();
        }
        p.set_objective(vec![(0, 9.), (1, 9.), (2, 1.)], 0.)
            .unwrap();
        let totals = add_totals(&mut p);
        assert_eq!(totals, vec![(3, 9.)]);
        p.validate().unwrap();
        for a in 0..=2 {
            for b in 0..=3 {
                let x = [a as f64, b as f64, 0., (a + b) as f64];
                assert_eq!(p.constraints.last().unwrap().expr.eval(&x), 0.);
            }
        }
    }
}
