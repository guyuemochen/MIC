//! 利用等式消去共享余量，先找可行点，再用已求得的完整根解和检验数搜索小候选集。
//! 所有还原点在完整工作模型重新验收。受限子问题的界不传给主树。
mod formulation;
mod packing_lns;
use crate::{
    branch_bound::{
        branch_and_bound_with_preserved, BbOptions, BbStats, BinaryPropagationMode, SearchControl,
        SearchPrelude,
    },
    heuristics::{is_choice_row, private_slack_start, round_and_check, Incumbent},
    model::{Problem, RelOp, Tolerances, VarType},
    presolve::{postsolve, presolve_for_solve},
    search_ctrl::gap_percent,
    simplex::{Bounds, LpSolution, LpStatus},
};
use std::{
    cell::RefCell,
    time::{Duration, Instant},
};

struct Pool<'a> {
    original: &'a Problem,
    model: Problem,
    offset: Vec<f64>,
    count: Vec<usize>,
    fill: Vec<(usize, usize, f64)>,
    group: Vec<usize>,
    num_groups: usize,
}

impl Pool<'_> {
    fn project(
        &self,
        keep: &[bool],
        end: Instant,
    ) -> Option<super::restricted_primal::Restriction> {
        let active: Vec<_> = keep
            .iter()
            .enumerate()
            .map(|(j, &yes)| yes && self.model.vars[j].ub != 0. && self.count[j] > 0)
            .collect();
        super::restricted_primal::project_zero(&self.model, &active, end)
    }

    fn candidates(
        &self,
        root: &[f64],
        reduced_costs: &[f64],
        present: &[bool],
        per_group: usize,
    ) -> Vec<bool> {
        // 消掉的是可由等式恢复的余量。剩余班次列保留原下标；平移下界之后比较支撑。
        let mut keep: Vec<_> = (0..self.model.vars.len())
            .map(|j| {
                self.model.vars[j].var_type != VarType::Binary
                    || present[j]
                    || root[j] > self.offset[j] + 1e-7
            })
            .collect();
        let mut buckets = vec![Vec::new(); self.num_groups + 1];
        for (j, &g) in self.group.iter().enumerate() {
            if !keep[j] {
                buckets[g].push(j);
            }
        }
        for (g, b) in buckets.iter_mut().enumerate() {
            let k = if g == self.num_groups {
                8 * per_group
            } else {
                per_group
            };
            if b.len() > k {
                b.select_nth_unstable_by(k, |&j, &h| {
                    reduced_costs[j]
                        .total_cmp(&reduced_costs[h])
                        .then_with(|| j.cmp(&h))
                });
            }
            for &j in b.iter().take(k) {
                keep[j] = true;
            }
        }
        keep
    }

    fn restore(&self, columns: &[usize], point: &[f64], tol: &Tolerances) -> Option<Incumbent> {
        if columns.len() != point.len() {
            return None;
        }
        let mut x = vec![0.; self.offset.len()];
        for (&j, &v) in columns.iter().zip(point) {
            x[j] = v;
        }
        for &(i, j, a) in &self.fill {
            let r = &self.model.constraints[i];
            x[j] = (r.rhs - r.expr.eval(&x)) / a;
        }
        for (v, base) in x.iter_mut().zip(&self.offset) {
            *v += base;
        }
        let inc = round_and_check(self.original, &Bounds::from_problem(self.original), &x, tol)?;
        super::checked_incumbent(self.original, inc.values)
    }
}

fn packing_items(model: &Problem, end: Instant) -> Option<Vec<packing_lns::Item>> {
    let n = model.vars.len();
    let mut coefficients = vec![Vec::new(); n];
    let mut blocked = vec![false; n];
    for (i, row) in model.constraints.iter().enumerate() {
        if Instant::now() >= end {
            return None;
        }
        let sign = match row.op {
            RelOp::Le => 1.,
            RelOp::Ge => -1.,
            RelOp::Eq => return None,
        };
        let capacity = sign * row.rhs;
        if capacity < 0. {
            return None;
        }
        for &(j, a) in &row.expr.terms {
            let a = sign * a;
            if a < 0. {
                return None;
            }
            if a == 0. {
                continue;
            }
            if a > capacity {
                blocked[j] = true;
            } else {
                let share = a / capacity;
                if share <= 0. || !share.is_finite() {
                    return None;
                }
                coefficients[j].push((i, share));
            }
        }
    }
    let mut costs = vec![0.; n];
    for &(j, a) in &model.obj.terms {
        costs[j] = a;
    }
    let items = coefficients
        .into_iter()
        .enumerate()
        .filter_map(|(j, rows)| {
            let col = &model.vars[j];
            (!blocked[j]
                && col.lb == 0.
                && col.ub >= 1.
                && col.var_type != VarType::Continuous
                && costs[j] < 0.
                && !rows.is_empty())
            .then_some(packing_lns::Item {
                column: j,
                profit: -costs[j],
                rows,
            })
        })
        .collect();
    (Instant::now() < end).then_some(items)
}

fn packing_incumbent(
    pool: &Pool<'_>,
    root: &LpSolution,
    reduced_costs: &[f64],
    tol: &Tolerances,
    deadline: Instant,
    nodes: u64,
) -> (Option<Incumbent>, u64) {
    let end = deadline.min(Instant::now() + Duration::from_secs(2));
    let root: Vec<_> = root
        .x
        .iter()
        .zip(&pool.offset)
        .map(|(a, b)| a - b)
        .collect();
    let Some(items) = packing_items(&pool.model, end) else {
        return (None, 0);
    };
    let (cols, nodes) = packing_lns::selected(
        pool.model.constraints.len(),
        &items,
        &root,
        reduced_costs,
        end,
        nodes,
    );
    if Instant::now() >= deadline || cols.is_empty() {
        return (None, nodes);
    }
    let mut x = vec![0.; pool.model.vars.len()];
    for j in cols {
        x[j] = 1.;
    }
    let columns: Vec<_> = (0..x.len()).collect();
    (pool.restore(&columns, &x, tol), nodes)
}

fn collect_choice_groups(p: &Problem, equalities: bool) -> (Vec<usize>, usize) {
    let mut groups = vec![usize::MAX; p.vars.len()];
    let mut count = 0;
    for row in &p.constraints {
        if (row.op == RelOp::Eq) != equalities
            || row.expr.terms.len() <= 3
            || !is_choice_row(p, row)
        {
            continue;
        }
        let mut used = false;
        for &(j, _) in &row.expr.terms {
            if groups[j] == usize::MAX {
                groups[j] = count;
                used = true;
            }
        }
        if used {
            count += 1;
        }
    }
    for group in &mut groups {
        if *group == usize::MAX {
            *group = count;
        }
    }
    (groups, count)
}

fn choice_groups(p: &Problem) -> Option<(Vec<usize>, usize, bool)> {
    // 单侧互斥行优先。只有等式表达互斥时，才按恰选一行分组。
    for equalities in [false, true] {
        let (groups, count) = collect_choice_groups(p, equalities);
        if count > 0 {
            return Some((groups, count, equalities));
        }
    }
    None
}

/// 共用调用方截止时间与节点预算；正文仅发布已经原工作模型检查的可行费用和向量。
pub(crate) fn search(
    p: &Problem,
    root: &LpSolution,
    reduced_costs: &[f64],
    options: &BbOptions,
    tol: Tolerances,
    deadline: Instant,
    node_budget: u64,
    target: f64,
    accept: &dyn Fn(Incumbent),
) -> BbStats {
    let mut stats = BbStats::default();
    if node_budget == 0
        || Instant::now() >= deadline
        || root.status != LpStatus::Optimal
        || root.x.len() != p.vars.len()
        || reduced_costs.len() != p.vars.len()
        || root.x.iter().chain(reduced_costs).any(|v| !v.is_finite())
    {
        return stats;
    }
    let Some((mut model, offset)) = formulation::projected(p, &tol, Some(deadline)) else {
        return stats;
    };
    let Some(seed) =
        private_slack_start(&model, &Bounds::from_problem(&model), &tol, Some(deadline))
    else {
        return stats;
    };
    let present: Vec<_> = seed.values.iter().map(|&v| v != 0.).collect();
    let x = seed
        .values
        .iter()
        .zip(&offset)
        .map(|(v, b)| v + b)
        .collect();
    let Some(seed) = super::checked_incumbent(p, x) else {
        return stats;
    };
    if Instant::now() >= deadline {
        return stats;
    }
    let best = RefCell::new(seed.clone());
    accept(seed);
    let Some((count, fill)) = formulation::eliminate(&mut model, deadline) else {
        return stats;
    };
    let Some((groups, ng, equalities)) = choice_groups(p) else {
        return stats;
    };
    let pool = Pool {
        original: p,
        model,
        offset,
        count,
        fill,
        group: groups,
        num_groups: ng,
    };
    // 保留根支撑和已选列。宽等式组较多时，每组只补少量低检验数列。
    let candidate_end = deadline.min(Instant::now() + Duration::from_secs(2));
    let keep = pool.candidates(
        &root.x,
        reduced_costs,
        &present,
        if equalities { 2 } else { 16 },
    );
    let Some(spent) = search_pool(
        &pool,
        &keep,
        options,
        tol,
        candidate_end,
        node_budget,
        &best,
        accept,
    ) else {
        // 子树异常返回时，无法确认实际已用的节点数。
        return stats;
    };
    stats = spent;
    // 小候选模型耗尽局部时限后，仍可用父级余下预算在非负装箱上快速补选。
    // 局部枚举只上报原模型可行解和已用节点，局部界不参与全局剪枝。
    if !equalities
        && pool.group.iter().filter(|&&g| g == ng).count() > pool.group.len() / 2
        && gap_percent(p.sense, best.borrow().objective, root.objective)
            .is_some_and(|gap| gap > 2. * target)
        && stats.nodes < node_budget
        && Instant::now() < deadline
    {
        let (inc, nodes) = packing_incumbent(
            &pool,
            root,
            reduced_costs,
            &tol,
            deadline,
            node_budget - stats.nodes,
        );
        stats.nodes += nodes;
        if let Some(inc) = inc.filter(|inc| inc.objective < best.borrow().objective) {
            accept(inc);
        }
    }
    stats
}

fn search_pool(
    pool: &Pool<'_>,
    keep: &[bool],
    options: &BbOptions,
    tol: Tolerances,
    end: Instant,
    node_budget: u64,
    best: &RefCell<Incumbent>,
    accept: &dyn Fn(Incumbent),
) -> Option<BbStats> {
    let Some(restriction) = pool.project(keep, end) else {
        return Some(BbStats::default());
    };
    let small = restriction.model;
    let columns = restriction.columns;
    if columns.len() >= 4096 || Instant::now() >= end {
        return Some(BbStats::default());
    }
    let start = Instant::now();
    let millis = end
        .saturating_duration_since(start)
        .as_millis()
        .min(u64::MAX as u128) as u64;
    if millis == 0 {
        return Some(BbStats::default());
    }
    let capture = cfg!(feature = "certified_row_fallback")
        && options.binary_row_propagation.resolve(1) == BinaryPropagationMode::On;
    let Ok((reduced, map, preserved)) = presolve_for_solve(&small, &tol, capture) else {
        return Some(BbStats::default());
    };
    if Instant::now() >= end {
        return Some(BbStats::default());
    }
    let warm_seed = {
        let point = best.borrow();
        let x: Vec<_> = columns
            .iter()
            .map(|&j| point.values[j] - pool.offset[j])
            .collect();
        let seed = Incumbent {
            objective: small.obj.eval(&x) + small.obj_offset,
            values: x,
        };
        super::work_seed(&small, &reduced, Some(&map), None, &seed, tol, &|x| {
            postsolve(&map, &x, columns.len(), &tol)
        })
    };
    let submit = |y: &[f64]| {
        if Instant::now() >= end || y.len() != reduced.vars.len() {
            return;
        }
        let x = postsolve(&map, y, columns.len(), &tol);
        if let Some(inc) = pool.restore(&columns, &x, &tol) {
            if inc.objective < best.borrow().objective && Instant::now() < end {
                best.replace(inc.clone());
                accept(inc);
            }
        }
    };
    let sub = BbOptions {
        time_limit_ms: millis,
        node_limit: node_budget,
        ..*options
    };
    let ctrl = SearchControl {
        threads: 1,
        mipgap_percent: Some(1.),
        ..Default::default()
    };
    let done = branch_and_bound_with_preserved(
        &reduced,
        &sub,
        &ctrl,
        tol,
        None,
        preserved,
        Some(&submit),
        SearchPrelude {
            start: Some(start),
            incumbent: warm_seed,
            ..Default::default()
        },
    )
    .ok()?;
    if let Some(x) = &done.best_values {
        submit(x);
    }
    Some(done.stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RelOp, Sense};
    use std::time::Duration;

    #[test]
    fn packing_normalizes_signed_capacity_and_excludes_blocked_columns() {
        let mut p = Problem::new("装箱容量", Sense::Minimize);
        for j in 0..3 {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        p.set_objective(vec![(0, -3.), (1, -2.), (2, -1.)], 0.)
            .unwrap();
        p.add_constraint(None, vec![(0, -1.), (1, -2.), (2, -3.)], RelOp::Ge, -2.)
            .unwrap();
        p.add_constraint(None, vec![(1, 1.)], RelOp::Le, 0.)
            .unwrap();
        let end = Instant::now() + Duration::from_secs(1);
        let items = packing_items(&p, end).unwrap();
        assert_eq!(items.len(), 1);
        assert_eq!(items[0].column, 0);
        assert_eq!(items[0].profit, 3.);
        assert_eq!(items[0].rows, vec![(0, 0.5)]);
        assert!(packing_items(&p, Instant::now()).is_none());
        p.add_constraint(None, vec![(0, -1.), (1, 2.)], RelOp::Le, 0.)
            .unwrap();
        assert!(packing_items(&p, end).is_none());
    }

    #[test]
    fn choice_groups_use_equalities_only_without_wide_inequalities() {
        let mut p = Problem::new("互斥行分组", Sense::Minimize);
        for j in 0..8 {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        p.add_constraint(None, (0..4).map(|j| (j, 1.)).collect(), RelOp::Eq, 1.)
            .unwrap();
        p.add_constraint(None, vec![(4, 1.), (5, 1.)], RelOp::Le, 1.)
            .unwrap();
        let (groups, count, equalities) = choice_groups(&p).unwrap();
        assert_eq!(groups, vec![0, 0, 0, 0, 1, 1, 1, 1]);
        assert_eq!(count, 1);
        assert!(equalities);

        p.add_constraint(None, (4..8).map(|j| (j, -2.)).collect(), RelOp::Ge, -2.)
            .unwrap();
        let (groups, count, equalities) = choice_groups(&p).unwrap();
        assert_eq!(groups, vec![1, 1, 1, 1, 0, 0, 0, 0]);
        assert_eq!(count, 1);
        assert!(!equalities);
    }

    #[test]
    fn fractional_eliminated_integer_cannot_be_published() {
        let mut p = Problem::new("整数余量", Sense::Minimize);
        let x = p.add_var("x", VarType::Binary).unwrap();
        let slack = p.add_var("余量", VarType::Integer).unwrap();
        p.add_constraint(None, vec![(x, 1.), (slack, 2.)], RelOp::Eq, 1.)
            .unwrap();
        p.set_objective(vec![(x, -1.), (slack, 2.)], 5.).unwrap();
        let mut model = p.clone();
        let (count, fill) =
            formulation::eliminate(&mut model, Instant::now() + Duration::from_secs(1)).unwrap();
        let pool = Pool {
            original: &p,
            model,
            offset: vec![0.; 2],
            count,
            fill,
            group: vec![0; 2],
            num_groups: 1,
        };
        let tol = Tolerances::default();
        assert!(pool.restore(&[x], &[0.], &tol).is_none());
        let inc = pool.restore(&[x], &[1.], &tol).unwrap();
        assert_eq!(inc.values, vec![1., 0.]);
        assert_eq!(inc.objective, 4.);
    }

    #[test]
    fn root_support_uses_shifted_bounds_and_keeps_existing_choices() {
        let mut p = Problem::new("根支撑", Sense::Minimize);
        for j in 0..5 {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        let pool = Pool {
            original: &p,
            model: p.clone(),
            offset: vec![0., 0., 0., 0., 1.],
            count: vec![1; 5],
            fill: vec![],
            group: vec![0; 5],
            num_groups: 1,
        };
        // 第0列在根解，第1列在原排班。第2及3列检验数并列，配额1取较小下标。
        // 第4列的原值1恰好等于平移下界，在候选模型中为0。
        let kept = pool.candidates(
            &[0.25, 0., 0., 0., 1.],
            &[9., 8., 1., 1., 7.],
            &[false, true, false, false, false],
            1,
        );
        assert_eq!(kept, vec![true, true, true, false, false]);
    }

    #[test]
    fn expired_or_zero_budget_does_not_publish() {
        let p = Problem::new("空", Sense::Minimize);
        let opts = BbOptions::default();
        let root = LpSolution {
            status: LpStatus::Optimal,
            objective: 0.,
            x: vec![],
            duals: vec![],
        };
        let t = Tolerances::default();
        for (deadline, nodes) in [
            (Instant::now(), 1),
            (Instant::now() + Duration::from_secs(1), 0),
        ] {
            assert_eq!(
                search(&p, &root, &[], &opts, t, deadline, nodes, 1., &|_| panic!(
                    "预算不足"
                ))
                .nodes,
                0
            );
        }
    }
}
