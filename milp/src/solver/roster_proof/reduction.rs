//! 用一份原始根的行乘子证明哪些取正的整数列费用已经超过可行解。
use super::{ProofMode, ProofOutcome};
use crate::cuts::Cut;
use crate::model::{Problem, RelOp, Sense, Tolerances, VarType};
use crate::simplex::{Bounds, LpSolution, LpStatus};
use std::time::Instant;

pub(super) struct Screen {
    keep: Vec<bool>,
    excluded_lower: f64,
}

pub(super) fn screen_columns(
    p: &Problem,
    root: &LpSolution,
    upper: f64,
    tol: &Tolerances,
    end: Instant,
) -> Option<Screen> {
    if p.sense != Sense::Minimize
        || root.status != LpStatus::Optimal
        || !upper.is_finite()
        || root.x.len() != p.vars.len()
        || root.duals.len() != p.constraints.len()
        || p.vars.iter().any(|v| !(v.lb >= 0. && v.lb.is_finite()))
        || !super::usable_tolerances(tol)
    {
        return None;
    }

    // c·x = y·Ax + (c-Aᵀy)·x。乘子只需方向正确，毋须相信求解器给出的检验数或原始目标。
    let mut rc = vec![0.; p.vars.len()];
    for &(j, c) in &p.obj.terms {
        rc[j] = c;
    }
    let mut lower = p.obj_offset;
    for (row, &y) in p.constraints.iter().zip(&root.duals) {
        if !y.is_finite() || Instant::now() >= end {
            return None;
        }
        // 轻微向零收缩，让有正罚价的无上界余量为对偶舍入误差预留空间。
        let y = match row.op {
            RelOp::Le => y.min(0.),
            RelOp::Ge => y.max(0.),
            RelOp::Eq => y,
        } * (1. - 2f64.powi(-24));
        if y == 0. {
            continue;
        }
        let mut norm = 0.;
        for &(j, a) in &row.expr.terms {
            rc[j] = (-a).mul_add(y, rc[j]).next_down();
            norm = (norm + a.abs()).next_up();
        }
        // 计入全行终检允许的违反量。所有乘积之和只向不利于固定的方向舍入。
        let scale = ((1. + row.rhs.abs()).next_up() + norm).next_up();
        let margin = (scale * super::super::FINAL_VIOLATION_REL).next_up().max(
            (tol.integrality
                .mul_add(norm, (10. * tol.feasibility).next_up()))
            .next_up(),
        );
        lower = y.mul_add(row.rhs, lower).next_down();
        lower = (-y.abs()).mul_add(margin, lower).next_down();
    }
    // 无限上界列的检验数即便仅向下舍入一位也不能直接忽略。用全局正系数行
    // 推出有限上界，按最坏端点计入费用；推不出有限上界仍放弃该证书。
    let mut caps: Vec<_> = p.vars.iter().map(|v| v.ub).collect();
    for row in &p.constraints {
        if Instant::now() >= end {
            return None;
        }
        if let Some(cap) = super::guarded_capacity(row, p, tol) {
            for &(j, a) in &row.expr.terms {
                caps[j] = caps[j].min((cap / a).next_up());
            }
        }
    }
    for ((v, &r), &cap) in p.vars.iter().zip(&rc).zip(&caps) {
        if !r.is_finite() {
            return None;
        }
        let b = if r < 0. { cap } else { v.lb };
        if !b.is_finite() {
            return None;
        }
        if b != 0. && r != 0. {
            lower = r.mul_add(b, lower).next_down();
        }
    }
    if !lower.is_finite() || lower > upper {
        return None;
    }
    // 数值余量只会减少固定。限定零位列；对当前在上界的列以及连续列不做替换。
    let guard = (upper + (1e-9 * (1. + upper.abs())).next_up()).next_up();
    let mut excluded_lower = f64::INFINITY;
    let keep = (0..p.vars.len())
        .map(|j| {
            let v = &p.vars[j];
            let positive_lower = (lower + rc[j]).next_down();
            if v.lb == 0.
                && v.ub >= 1.
                && root.x[j] == 0.
                && v.var_type != VarType::Continuous
                && rc[j] > 0.
                && positive_lower > guard
            {
                excluded_lower = excluded_lower.min(positive_lower);
                false
            } else {
                true
            }
        })
        .collect();
    Some(Screen {
        keep,
        excluded_lower,
    })
}

pub(in crate::solver) fn prove_compact(
    work: &Problem,
    cuts: &[Cut],
    root: &LpSolution,
    upper: f64,
    gap: f64,
    tol: &Tolerances,
    end: Instant,
    nodes: u64,
    paired: bool,
) -> ProofOutcome {
    let mode = ProofMode {
        paired,
        can_reduce: true,
    };
    let unchanged = || {
        super::prove_domain(
            work,
            cuts,
            root.objective,
            upper,
            gap,
            tol,
            end,
            nodes,
            mode,
        )
    };
    if work.vars.len() < 4096 || nodes == 0 {
        return unchanged();
    }
    screen_columns(work, root, upper, tol, end)
        .and_then(|screen| {
            prove_subset(
                work,
                cuts,
                &screen,
                root.objective,
                upper,
                gap,
                tol,
                end,
                nodes,
                mode,
            )
        })
        .unwrap_or_else(unchanged)
}

pub(super) fn prove_subset(
    work: &Problem,
    cuts: &[Cut],
    screen: &Screen,
    root_bound: f64,
    upper: f64,
    gap: f64,
    tol: &Tolerances,
    end: Instant,
    nodes: u64,
    mode: ProofMode,
) -> Option<ProofOutcome> {
    // 强化根可能含额外合计和固定零列。只从工作模型删列，辅助列的下界参与取最小仍保守。
    let keep = screen.keep.get(..work.vars.len())?;
    let n = keep.iter().filter(|&&b| b).count();
    if n > keep.len() / 2 || nodes == 0 {
        return None;
    }
    let reduced = crate::solver::restricted_primal::project_zero(work, keep, end)?;
    let mut col = vec![0; keep.len()];
    for (k, j) in (0..keep.len()).filter(|&j| keep[j]).enumerate() {
        col[j] = k;
    }
    let cuts = cuts
        .iter()
        .map(|cut| Cut {
            terms: cut
                .terms
                .iter()
                .filter(|(j, _)| keep[*j])
                .map(|&(j, a)| (col[j], a))
                .collect(),
            op: cut.op,
            rhs: cut.rhs,
        })
        .collect::<Vec<_>>();
    let mut result = super::prove_domain(
        &reduced.model,
        &cuts,
        root_bound,
        upper,
        gap,
        tol,
        end,
        nodes,
        mode,
    );
    // 被删域用定价证书，保留域用子搜索，二者都可继承原父界。
    result.bound = root_bound.max(result.bound.min(screen.excluded_lower).min(upper));
    result.candidate = result.candidate.and_then(|inc| {
        let mut lifted = reduced.lift(work, &Bounds::from_problem(work), &inc.values, tol)?;
        crate::solver::final_check(work, &mut lifted.values).ok()?;
        lifted.objective = work.obj.eval(&lifted.values) + work.obj_offset;
        (lifted.objective.is_finite()
            && lifted
                .values
                .iter()
                .zip(&work.vars)
                .all(|(&v, col)| col.lb <= v && v <= col.ub))
        .then_some(lifted)
    });
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn projected_model_lifts_an_improving_point_and_preserves_discarded_domain_bound() {
        let mut p = Problem::new("压缩后还原证据", Sense::Minimize);
        p.add_vars(
            (0..4096)
                .map(|i| crate::model::Variable {
                    name: format!("x{i}"),
                    var_type: VarType::Binary,
                    lb: 0.,
                    ub: 1.,
                })
                .collect(),
        )
        .unwrap();
        let mut costs = (0..4096).map(|j| (j, 100.)).collect::<Vec<_>>();
        for j in [0, 2, 4] {
            costs[j].1 = -1.;
        }
        p.set_objective(costs, 7.).unwrap();
        for (a, b) in [(0, 2), (2, 4), (0, 4)] {
            p.add_constraint(None, vec![(a, 1.), (b, 1.)], RelOp::Le, 1.)
                .unwrap();
        }
        let mut x = vec![0.; 4096];
        for j in [0, 2, 4] {
            x[j] = 0.5;
        }
        let root = LpSolution {
            status: LpStatus::Optimal,
            objective: 5.5,
            x,
            duals: vec![-0.5; 3],
        };
        let cut = Cut {
            terms: vec![(0, 1.), (2, 1.), (4, 1.)],
            op: RelOp::Le,
            rhs: 1.,
        };
        let t = Tolerances::default();
        let result = prove_compact(
            &p,
            &[cut],
            &root,
            7.,
            1.,
            &t,
            Instant::now() + Duration::from_secs(2),
            1,
            false,
        );
        assert_eq!(result.bound, 6.);
        assert_eq!(result.nodes, 1);
        let inc = result.candidate.unwrap();
        assert_eq!(inc.values.len(), 4096);
        assert_eq!(inc.objective, p.obj.eval(&inc.values) + p.obj_offset);
        assert_eq!(inc.objective, 6.);
        assert!(inc
            .values
            .iter()
            .enumerate()
            .all(|(j, v)| [0, 2, 4].contains(&j) || *v == 0.));
    }

    #[test]
    fn weaker_certificate_for_the_excluded_domain_caps_the_global_bound() {
        let mut p = Problem::new("两个子域分别计界", Sense::Minimize);
        let x = p.add_var("选项", VarType::Binary).unwrap();
        let y = p.add_var("备选", VarType::Binary).unwrap();
        p.set_objective(vec![(x, -0.25), (y, 1.)], 6.).unwrap();
        // 已知 x 取1时的最小费用是5.75。固定 x=0 的模型只能证明自己的费用为6。
        let screen = Screen {
            keep: vec![false, true],
            excluded_lower: 5.75,
        };
        let result = prove_subset(
            &p,
            &[],
            &screen,
            5.5,
            7.,
            1.,
            &Tolerances::default(),
            Instant::now() + Duration::from_secs(2),
            1,
            ProofMode {
                paired: false,
                can_reduce: false,
            },
        )
        .unwrap();
        assert_eq!(result.bound, 5.75);
        assert_eq!(result.candidate.unwrap().objective, 6.);
    }

    #[test]
    fn exclusion_is_certified_with_offset_arbitrary_duals_and_ties() {
        let mut p = Problem::new("费用剔除", Sense::Minimize);
        for name in ["昂贵", "持平", "负费用", "连续", "上界"] {
            p.add_var(
                name,
                if name == "连续" {
                    VarType::Continuous
                } else {
                    VarType::Binary
                },
            )
            .unwrap();
        }
        p.set_var_bounds(3, 0., 1.).unwrap();
        p.set_objective(vec![(0, 15.), (1, 4.), (2, -4.), (3, 10.), (4, 18.)], 37.)
            .unwrap();
        p.add_constraint(None, vec![(1, 1.), (2, 1.)], RelOp::Le, 1.)
            .unwrap();
        let mut root = LpSolution {
            status: LpStatus::Optimal,
            objective: 33.,
            x: vec![0., 0., 1., 0., 1.],
            duals: vec![10.],
        };
        // 故意给上界行错误正号乘子；投影截为零。负目标列须计入箱下界。
        let t = Tolerances::default();
        let end = Instant::now() + Duration::from_secs(1);
        let keep = screen_columns(&p, &root, 37., &t, end).unwrap().keep;
        assert_eq!(keep, vec![false, true, true, true, true]);
        for mask in 0..32 {
            let x = (0..5).map(|j| ((mask >> j) & 1) as f64).collect::<Vec<_>>();
            if x[1] + x[2] <= 1. && !keep[0] && x[0] == 1. {
                assert!(p.obj.eval(&x) + p.obj_offset > 37.);
            }
        }
        root.duals[0] = f64::NAN;
        assert!(screen_columns(&p, &root, 37., &t, end).is_none());
        root.duals[0] = 0.;
        p.vars[2].ub = f64::INFINITY;
        // 第2列虽然没有显式上界，共用容量行仍证明它至多取1。
        assert!(screen_columns(&p, &root, 37., &t, end).is_some());
        p.constraints.clear();
        root.duals.clear();
        assert!(screen_columns(&p, &root, 37., &t, end).is_none());
    }
}
