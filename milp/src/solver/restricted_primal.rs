//! 只在完整根 LP 为零的有界整数列中缩小候选集。小模型仅用于寻找原问题的可行点。

use crate::branch_bound::{
    branch_and_bound_with_preserved, BbOptions, BbStats, BinaryPropagationMode, SearchControl,
    SearchPrelude,
};
use crate::heuristics::{is_choice_row, round_and_check, Incumbent};
use crate::model::{Constraint, LinearExpr, Problem, Tolerances, VarType};
use crate::presolve::{postsolve, presolve_for_solve};
use crate::simplex::{Bounds, LpSolution, LpStatus};
use std::cell::Cell;
use std::time::Instant;

const PER_GROUP: usize = 96;

pub(super) struct Restriction {
    pub(super) model: Problem,
    // 小模型下标到完整工作模型下标；缺席的列都可以精确取零。
    pub(super) columns: Vec<usize>,
}

fn restrict(
    p: &Problem,
    x: &[f64],
    costs: &[f64],
    per_group: usize,
    end: Instant,
) -> Option<Restriction> {
    let n = p.vars.len();
    if x.len() != n || costs.len() != n || x.iter().chain(costs).any(|v| !v.is_finite()) {
        return None;
    }
    let mut degree = vec![0_u8; n];
    let mut group = vec![0; n];
    let mut buckets = vec![Vec::new()];
    for row in &p.constraints {
        if Instant::now() >= end {
            return None;
        }
        let terms = &row.expr.terms;
        for &(j, _) in terms {
            degree[j] = degree[j].saturating_add(1);
        }
        if row.op != crate::model::RelOp::Eq && is_choice_row(p, row) {
            let label = buckets.len();
            for &(j, _) in terms {
                if group[j] == 0 {
                    group[j] = label;
                }
            }
            buckets.push(Vec::new());
        }
    }
    for (j, v) in p.vars.iter().enumerate() {
        if v.var_type != VarType::Continuous
            && v.lb == 0.
            && v.ub.is_finite()
            && v.ub >= 1.
            && degree[j] > 1
            && x[j] == 0.
        {
            buckets[group[j]].push(j);
        }
    }
    let mut keep = vec![true; n];
    let mut removed = 0;
    for (group, bucket) in buckets.iter_mut().enumerate() {
        if Instant::now() >= end {
            return None;
        }
        // 检验数采用 LP 引擎的内部最小化口径，未分组的列共用较宽的候选集。
        let cap = if group == 0 { per_group * 8 } else { per_group };
        if bucket.len() > cap {
            bucket.select_nth_unstable_by(cap, |&i, &j| {
                costs[i].partial_cmp(&costs[j]).unwrap().then(i.cmp(&j))
            });
            for &j in &bucket[cap..] {
                keep[j] = false;
            }
            removed += bucket.len() - cap;
        }
    }
    if removed == 0 {
        return None;
    }
    project_zero(p, &keep, end)
}

/// 删除掩码指定的零列，保持原行和保留列的顺序。
pub(super) fn project_zero(p: &Problem, keep: &[bool], end: Instant) -> Option<Restriction> {
    let n = p.vars.len();
    if keep.len() != n {
        return None;
    }
    let columns: Vec<_> = (0..n).filter(|&j| keep[j]).collect();
    let mut local = vec![usize::MAX; n];
    let mut model = Problem::new("候选列子问题", p.sense);
    for (k, &j) in columns.iter().enumerate() {
        local[j] = k;
        model.vars.push(p.vars[j].clone());
    }
    let project = |e: &LinearExpr| LinearExpr {
        terms: e
            .terms
            .iter()
            .filter(|&&(j, _)| keep[j])
            .map(|&(j, a)| (local[j], a))
            .collect(),
    };
    model.obj = project(&p.obj);
    model.obj_offset = p.obj_offset;
    for row in &p.constraints {
        if Instant::now() >= end {
            return None;
        }
        model.constraints.push(Constraint {
            name: None,
            expr: project(&row.expr),
            op: row.op,
            rhs: row.rhs,
        });
    }
    (Instant::now() < end).then_some(Restriction { model, columns })
}

impl Restriction {
    pub(super) fn lift(
        &self,
        p: &Problem,
        bounds: &Bounds,
        x: &[f64],
        tol: &Tolerances,
    ) -> Option<Incumbent> {
        if x.len() != self.columns.len() {
            return None;
        }
        let mut full = vec![0.; p.vars.len()];
        for (&j, &value) in self.columns.iter().zip(x) {
            full[j] = value;
        }
        let inc = round_and_check(p, bounds, &full, tol)?;
        inc.values
            .iter()
            .zip(&p.vars)
            .all(|(&v, col)| {
                col.lb <= v
                    && v <= col.ub
                    && (col.var_type == VarType::Continuous || v == v.round())
            })
            .then_some(inc)
    }
}

/// 用缩小模型寻找原工作模型的可行解。回调目标重算并包含目标常数；返回值只记工作量。
pub(super) fn search(
    work: &Problem,
    root: &LpSolution,
    reduced_costs: &[f64],
    options: &BbOptions,
    tol: Tolerances,
    deadline: Instant,
    accept: &dyn Fn(Incumbent),
    stop_after_node: Option<&dyn Fn(u64) -> bool>,
) -> BbStats {
    let start = Instant::now();
    // 毫秒向下取整，共用入口时刻；零在 B&B 表示不限时，不能把最后不足一毫秒传进去。
    let millis = deadline
        .saturating_duration_since(start)
        .as_millis()
        .min(u64::MAX as u128) as u64;
    if millis == 0
        || options.node_limit == 0
        || root.status != LpStatus::Optimal
        || !root.objective.is_finite()
    {
        return BbStats::default();
    }
    let Some(restriction) = restrict(work, &root.x, reduced_costs, PER_GROUP, deadline) else {
        return BbStats::default();
    };
    let capture = cfg!(feature = "certified_row_fallback")
        && options.binary_row_propagation.resolve(1) == BinaryPropagationMode::On;
    let Ok((reduced, map, preserved)) = presolve_for_solve(&restriction.model, &tol, capture)
    else {
        return BbStats::default();
    };
    if Instant::now() >= deadline {
        return BbStats::default();
    }
    let bounds = Bounds::from_problem(work);
    let best = Cell::new(None);
    let submit = |y: &[f64]| {
        if Instant::now() >= deadline || y.len() != reduced.vars.len() {
            return;
        }
        let x = postsolve(&map, y, restriction.columns.len(), &tol);
        if let Some(inc) = restriction.lift(work, &bounds, &x, &tol) {
            if best
                .get()
                .is_none_or(|old| super::better(work.sense, inc.objective, old))
                && Instant::now() < deadline
            {
                best.set(Some(inc.objective));
                accept(inc);
            }
        }
    };
    // 子搜索直接进入单线程 B&B，不再经过外层策略选择。只使用其解和工作量。
    let local = BbOptions {
        time_limit_ms: millis,
        ..*options
    };
    let ctrl = SearchControl {
        threads: 1,
        mipgap_percent: Some(0.),
        ..Default::default()
    };
    let prelude = SearchPrelude {
        start: Some(start),
        stop_after_node,
        ..Default::default()
    };
    match branch_and_bound_with_preserved(
        &reduced,
        &local,
        &ctrl,
        tol,
        None,
        preserved,
        Some(&submit),
        prelude,
    ) {
        Ok(result) => {
            if let Some(x) = &result.best_values {
                submit(x);
            }
            result.stats
        }
        Err(_) => BbStats::default(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RelOp, Sense};
    use crate::simplex::Simplex;
    use std::time::Duration;

    fn deadline() -> Instant {
        Instant::now() + Duration::from_secs(5)
    }

    #[test]
    fn keeps_live_columns_auxiliaries_bounds_and_objective_offset() {
        let mut p = Problem::new("筛选", Sense::Maximize);
        for j in 0..6 {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        let continuous = p.add_var("c", VarType::Continuous).unwrap();
        let free = p.add_var("z", VarType::Integer).unwrap();
        let positive = p.add_var("u", VarType::Integer).unwrap();
        p.set_var_bounds(positive, 1., 2.).unwrap();
        let private = p.add_var("v", VarType::Integer).unwrap();
        p.set_var_bounds(private, 0., 2.).unwrap();
        p.add_constraint(
            None,
            vec![(0, -2.), (1, -2.), (2, -2.), (3, -2.)],
            RelOp::Ge,
            -2.,
        )
        .unwrap();
        p.add_constraint(None, vec![(4, 1.), (5, 1.)], RelOp::Le, 1.)
            .unwrap();
        p.add_constraint(
            None,
            (0..p.vars.len()).map(|j| (j, 1.)).collect(),
            RelOp::Le,
            8.,
        )
        .unwrap();
        p.add_constraint(
            None,
            vec![(continuous, 1.), (free, -1.), (positive, 1.)],
            RelOp::Ge,
            -3.,
        )
        .unwrap();
        p.set_objective(
            vec![(1, 10.), (2, 20.), (3, 30.), (5, 40.), (private, 50.)],
            7.,
        )
        .unwrap();
        let mut x = vec![0.; p.vars.len()];
        x[0] = 1e-14; // 根 LP 中非零的列也要保留。
        let mut rc = vec![0.; p.vars.len()];
        rc[1] = 3.;
        rc[2] = -0.;
        rc[4] = 1.;
        let r = restrict(&p, &x, &rc, 1, deadline()).unwrap();
        assert_eq!(
            r.columns,
            vec![0, 2, 5, continuous, free, positive, private]
        );
        assert_eq!(r.model.obj_offset, 7.);
        assert_eq!(r.model.sense, Sense::Maximize);
        assert_eq!(r.model.obj.terms, vec![(1, 20.), (2, 40.), (6, 50.)]);
        assert_eq!(r.model.constraints.len(), p.constraints.len());
        assert_eq!(r.model.constraints[0].expr.terms, vec![(0, -2.), (1, -2.)]);
        assert_eq!(r.model.constraints[0].op, RelOp::Ge);
        assert_eq!(r.model.constraints[0].rhs, -2.);
        let bounds = Bounds::from_problem(&p);
        let mut y = vec![0.; r.columns.len()];
        y[1] = 1.;
        y[5] = 1.;
        let candidate = r.lift(&p, &bounds, &y, &Tolerances::default()).unwrap();
        assert_eq!(candidate.objective, 27.);
        assert_eq!(candidate.values[2], 1.);
        assert_eq!(candidate.values[1], 0.);
        y[5] = 0.; // 完整变量界不满足。
        assert!(r.lift(&p, &bounds, &y, &Tolerances::default()).is_none());
        y[5] = 1.;
        y[0] = 1.; // 违反完整模型的互斥行。
        assert!(r.lift(&p, &bounds, &y, &Tolerances::default()).is_none());
    }

    #[test]
    fn ungrouped_candidates_have_a_separate_quota_and_need_finite_root_data() {
        let mut p = Problem::new("无组", Sense::Minimize);
        for j in 0..12 {
            p.add_var(&format!("u{j}"), VarType::Binary).unwrap();
        }
        for rhs in [11., 10.] {
            p.add_constraint(None, (0..12).map(|j| (j, 1.)).collect(), RelOp::Le, rhs)
                .unwrap();
        }
        let mut costs: Vec<_> = (0..12).map(|i| -(i as f64)).collect();
        let x = vec![0.; 12];
        let fit = restrict(&p, &x, &costs, 1, deadline()).unwrap();
        assert_eq!(fit.columns, (4..12).collect::<Vec<_>>());
        assert!(restrict(&p, &x, &costs, 2, deadline()).is_none());
        assert!(restrict(&p, &x, &costs, 1, Instant::now()).is_none());
        assert!(restrict(&p, &x[1..], &costs, 1, deadline()).is_none());
        costs[0] = f64::NAN;
        assert!(restrict(&p, &x, &costs, 1, deadline()).is_none());
    }

    #[test]
    fn infeasible_restricted_model_emits_no_global_claim() {
        let mut p = Problem::new("子问题无解", Sense::Minimize);
        for j in 0..=PER_GROUP {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        let last = PER_GROUP;
        let t = p.add_var("t", VarType::Binary).unwrap();
        for rhs in [1., 2.] {
            p.add_constraint(None, (0..=last).map(|j| (j, 1.)).collect(), RelOp::Le, rhs)
                .unwrap();
        }
        // LP 取 t=1/2、最后一列为零；完整整数解必须选最后一列。
        p.add_constraint(None, vec![(last, 1.), (t, 2.)], RelOp::Eq, 1.)
            .unwrap();
        p.set_objective(vec![(last, 1.)], 0.).unwrap();
        let tol = Tolerances::default();
        let bounds = Bounds::from_problem(&p);
        let mut engine = Simplex::new(&p, tol).unwrap();
        let root = engine.solve(&bounds);
        let rc = engine.structural_reduced_costs();
        assert_eq!((root.x[t], root.x[last]), (0.5, 0.));
        let selected = restrict(&p, &root.x, &rc, PER_GROUP, deadline()).unwrap();
        assert!(!selected.columns.contains(&last));
        search(
            &p,
            &root,
            &rc,
            &BbOptions::default(),
            tol,
            deadline(),
            &|_| panic!("缩小模型没有整数解"),
            None,
        );
        let mut feasible = vec![0.; p.vars.len()];
        feasible[last] = 1.;
        assert_eq!(
            round_and_check(&p, &bounds, &feasible, &tol)
                .unwrap()
                .objective,
            1.
        );
    }

    #[test]
    fn publishes_only_checked_work_points_and_honors_empty_budgets() {
        let mut p = Problem::new("候选回填", Sense::Maximize);
        let n = PER_GROUP + 4;
        for j in 0..n {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        for rhs in [1., 3.] {
            p.add_constraint(None, (0..n).map(|j| (j, 1.)).collect(), RelOp::Le, rhs)
                .unwrap();
        }
        p.set_objective((0..n).map(|j| (j, 5.)).collect(), 9.)
            .unwrap();
        let t = Tolerances::default();
        let mut lp = Simplex::new(&p, t).unwrap();
        let mut root = lp.solve(&Bounds::from_problem(&p));
        let rc = lp.structural_reduced_costs();
        let best = Cell::new(0.);
        let publish = |inc: Incumbent| {
            assert_eq!(inc.values.len(), n);
            assert_eq!(p.obj.eval(&inc.values) + 9., inc.objective);
            assert!(inc.values.iter().all(|&v| v == 0. || v == 1.));
            assert!(inc.values.iter().sum::<f64>() <= 1.);
            best.set(inc.objective);
        };
        let opts = BbOptions::default();
        assert_eq!(
            search(&p, &root, &rc, &opts, t, Instant::now(), &publish, None).nodes,
            0
        );
        assert_eq!(
            search(
                &p,
                &root,
                &rc,
                &BbOptions {
                    node_limit: 0,
                    ..opts
                },
                t,
                deadline(),
                &publish,
                None
            )
            .nodes,
            0
        );
        assert_eq!(best.get(), 0.);
        assert!(search(&p, &root, &rc, &opts, t, deadline(), &publish, None).nodes > 0);
        assert_eq!(best.get(), 14.);
        root.status = LpStatus::IterationLimit;
        assert_eq!(
            search(
                &p,
                &root,
                &rc,
                &opts,
                t,
                deadline(),
                &|_| panic!("不完整根 LP 不应启动"),
                None
            )
            .nodes,
            0
        );
    }
}
