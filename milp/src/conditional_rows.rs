//! 将重复使用同一组二元列的短约束合并成有限选项，再求解组间的联结约束。
//! 每个选项是一份满足全部组内约束的原变量赋值。若另一选项目标更小或相等，
//! 且在每条跨组上界中的用量都不大于它，就可以在任意整数解中替换它。
//! 每组用全零选择表示第一个选项，其余选项至多选一；系数记录相对首项的差。
//! 所有取舍用精确整数完成。选项过多或缺少这种结构时交回通用求解路径。

use crate::model::{Constraint, LinearExpr, Problem, RelOp, Sense, Tolerances, VarType, Variable};
use std::collections::BTreeMap;
use std::time::Instant;

const MAX_WIDTH: usize = 16;
const MAX_SUBSETS: usize = 1 << 20;
const MAX_LINKS: usize = 24;
const MAX_ROW_CHECKS: usize = 1 << 25;
const MAX_CHOICES: usize = 256;

fn expired(deadline: Option<Instant>) -> bool {
    deadline.is_some_and(|end| Instant::now() >= end)
}

/// 所有原始和改写后的二元活动度及目标都可逐项精确累加。行容差须小于一个整数单位。
fn integral(expr: &LinearExpr, constant: f64) -> Option<f64> {
    let mut norm = 0.0;
    for &(_, a) in &expr.terms {
        if a != a.trunc() {
            return None;
        }
        norm += a.abs();
    }
    (constant == constant.trunc() && norm + constant.abs() < (1u64 << 40) as f64).then_some(norm)
}

fn separated(p: &Problem, tol: Tolerances) -> bool {
    integral(&p.obj, p.obj_offset).is_some()
        && p.constraints.iter().all(|r| {
            integral(&r.expr, r.rhs).is_some_and(|norm| {
                10.0 * tol.feasibility + tol.integrality * norm < 1.0
                    && crate::solver::FINAL_VIOLATION_REL * (1.0 + r.rhs.abs() + norm) < 1.0
            })
        })
}

#[derive(Debug)]
struct GroupMap {
    columns: Vec<usize>,
    start: usize,
    masks: Vec<usize>,
}

#[derive(Debug)]
pub(crate) struct Reformulation {
    pub(crate) problem: Problem,
    groups: Vec<GroupMap>,
    external: Vec<usize>,
    original_size: usize,
    integer_tol: f64,
}

impl Reformulation {
    pub(crate) fn restore(&self, point: &[f64]) -> Option<Vec<f64>> {
        if point.len() != self.problem.vars.len()
            || point.iter().any(|&v| {
                !(v.round() == 0.0 || v.round() == 1.0) || (v - v.round()).abs() > self.integer_tol
            })
        {
            return None;
        }
        let mut x = vec![0.0; self.original_size];
        for group in &self.groups {
            let mut choice = 0;
            for k in 1..group.masks.len() {
                if point[group.start + k - 1].round() == 1.0 {
                    if choice != 0 {
                        return None;
                    }
                    choice = k;
                }
            }
            let mask = group.masks[choice];
            for (bit, &j) in group.columns.iter().enumerate() {
                x[j] = ((mask >> bit) & 1) as f64;
            }
        }
        for (&j, &value) in self
            .external
            .iter()
            .zip(&point[point.len() - self.external.len()..])
        {
            x[j] = value.round();
        }
        Some(x)
    }
}

struct Partition {
    groups: Vec<Vec<usize>>,
    local: Vec<Vec<usize>>,
    links: Vec<usize>,
    external: Vec<usize>,
    owner: Vec<usize>,
}

impl Partition {
    fn detect(p: &Problem) -> Option<Self> {
        let mut supports = BTreeMap::new();
        for row in &p.constraints {
            if (2..=MAX_WIDTH).contains(&row.expr.terms.len()) {
                let support: Vec<_> = row.expr.terms.iter().map(|&(j, _)| j).collect();
                *supports.entry(support).or_insert(0usize) += 1;
            }
        }
        let mut candidates: Vec<_> = supports
            .into_iter()
            .filter_map(|(s, count)| (count >= 2).then_some(s))
            .collect();
        candidates.sort_by(|a, b| b.len().cmp(&a.len()).then_with(|| a.cmp(b)));
        let mut owner = vec![usize::MAX; p.vars.len()];
        let mut groups = Vec::new();
        let mut work = 0;
        for cols in candidates {
            if cols.iter().all(|&j| owner[j] == usize::MAX) {
                work += 1 << cols.len();
                if work > MAX_SUBSETS {
                    return None;
                }
                for &j in &cols {
                    owner[j] = groups.len();
                }
                groups.push(cols);
            }
        }
        if groups.len() < 2 {
            return None;
        }
        groups.sort_by_key(|cols| cols[0]);
        for (g, cols) in groups.iter().enumerate() {
            for &j in cols {
                owner[j] = g;
            }
        }
        let external: Vec<_> = (0..p.vars.len())
            .filter(|&j| owner[j] == usize::MAX)
            .collect();
        if external.len() > MAX_LINKS {
            return None;
        }
        let mut local = vec![Vec::new(); groups.len()];
        let mut links = Vec::new();
        for (i, r) in p.constraints.iter().enumerate() {
            let &(j, _) = r.expr.terms.first()?;
            let group = owner[j];
            if group != usize::MAX && r.expr.terms.iter().all(|&(k, _)| owner[k] == group) {
                local[group].push(i);
            } else {
                links.push(i);
            }
        }
        let checks: u64 = groups
            .iter()
            .zip(&local)
            .map(|(g, r)| (1u64 << g.len()) * (r.len() + links.len()) as u64)
            .sum();
        if links.is_empty() || links.len() > MAX_LINKS || checks > MAX_ROW_CHECKS as u64 {
            return None;
        }
        Some(Self {
            groups,
            local,
            links,
            external,
            owner,
        })
    }
}

/// 对 Σaⱼxⱼ ≤ M y，若 aⱼ>0、M≥Σaⱼ，则每组选用任一正项都必须启用 y，
/// 而 y=1 时整行必定成立。因此可等价地逐组写成「该组使用」≤y。
fn switch_column(row: &Constraint, owner: &[usize]) -> Option<usize> {
    if row.rhs != 0.0 {
        return None;
    }
    let mut switch = None;
    let mut sum = 0;
    for &(j, a) in &row.expr.terms {
        if a < 0.0 && owner[j] == usize::MAX && switch.is_none() {
            switch = Some((j, -a as i64));
        } else if a > 0.0 && owner[j] != usize::MAX {
            sum += a as i64;
        } else {
            return None;
        }
    }
    switch
        .filter(|&(_, cap)| sum > 0 && sum <= cap)
        .map(|(j, _)| j)
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
struct Choice {
    cost: i64,
    used: Vec<i64>,
    mask: usize,
}

/// 支配关系同时比较目标与每条资源上界。相同签名保留最小原掩码。
fn consider(kept: &mut Vec<Choice>, cost: i64, used: &[i64], mask: usize) {
    if let Some(old) = kept
        .iter_mut()
        .find(|old| old.cost <= cost && old.used.iter().zip(used).all(|(a, b)| a <= b))
    {
        if old.cost == cost && old.used == used {
            old.mask = old.mask.min(mask);
        }
        return;
    }
    kept.retain(|old| !(cost <= old.cost && used.iter().zip(&old.used).all(|(a, b)| a <= b)));
    let position = kept.partition_point(|old| old.cost <= cost);
    kept.insert(
        position,
        Choice {
            cost,
            used: used.to_vec(),
            mask,
        },
    );
}

/// 按格雷码枚举，每步仅翻转一个变量，并增量更新本组有贡献的行。
fn choices(
    problem: &Problem,
    partition: &Partition,
    group: usize,
    switches: &[Option<usize>],
    costs: &[i64],
    deadline: Option<Instant>,
) -> Option<Vec<Choice>> {
    let cols = &partition.groups[group];
    let local = &partition.local[group];
    let mut changes = vec![Vec::new(); cols.len()];
    for (k, &i) in local.iter().chain(&partition.links).enumerate() {
        for &(j, a) in &problem.constraints[i].expr.terms {
            if let Ok(bit) = cols.binary_search(&j) {
                changes[bit].push((k, a as i64));
            }
        }
    }
    let rhs: Vec<_> = local
        .iter()
        .map(|&i| problem.constraints[i].rhs as i64)
        .collect();
    let mut activity = vec![0i64; local.len() + switches.len()];
    let mut kept = Vec::new();
    let mut used = vec![0; switches.len()];
    let mut cost = 0;
    let mut mask = 0;
    for k in 0usize..1 << cols.len() {
        if k % 256 == 0 && expired(deadline) {
            return None;
        }
        if k > 0 {
            let bit = k.trailing_zeros() as usize;
            mask ^= 1 << bit;
            let sign = if mask & (1 << bit) == 0 { -1 } else { 1 };
            cost += sign * costs[cols[bit]];
            for &(row, a) in &changes[bit] {
                activity[row] += sign * a;
            }
        }
        if rhs.iter().zip(&activity).any(|(rhs, a)| a > rhs) {
            continue;
        }
        for (k, switch) in switches.iter().enumerate() {
            let v = activity[local.len() + k];
            used[k] = if switch.is_some() {
                i64::from(v > 0)
            } else {
                v
            };
        }
        consider(&mut kept, cost, &used, mask);
        if kept.len() > MAX_CHOICES {
            return None;
        }
    }
    kept.sort_unstable();
    (!kept.is_empty()).then_some(kept)
}

/// 每条普通联结行按独立组及外部二元列求活动区间。最大值已满足的行可删；
/// 选定一项后，即使其他组都取最小值仍超过右端，该项也不可能出现在可行解中。
fn reduce_links(
    original: &Problem,
    parts: &Partition,
    switches: &[Option<usize>],
    options: &mut [Vec<Choice>],
    deadline: Option<Instant>,
) -> Option<Vec<bool>> {
    let mut redundant = vec![false; parts.links.len()];
    for (k, &i) in parts.links.iter().enumerate() {
        if expired(deadline) {
            return None;
        }
        if switches[k].is_some() {
            continue;
        }
        let row = &original.constraints[i];
        let rhs = row.rhs as i64;
        let (mut lower, mut upper) = (0, 0);
        for &(j, coeff) in &row.expr.terms {
            if parts.owner[j] == usize::MAX {
                lower += (coeff as i64).min(0);
                upper += (coeff as i64).max(0);
            }
        }
        let mut minima = Vec::with_capacity(options.len());
        for group in options.iter() {
            let min = group.iter().map(|p| p.used[k]).min()?;
            lower += min;
            upper += group.iter().map(|p| p.used[k]).max()?;
            minima.push(min);
        }
        if upper <= rhs {
            redundant[k] = true;
            for choice in options.iter_mut().flatten() {
                choice.used[k] = 0;
            }
        } else {
            for (group, min) in options.iter_mut().zip(minima) {
                group.retain(|p| lower - min + p.used[k] <= rhs);
                if group.is_empty() {
                    return None;
                }
            }
        }
    }
    // 删去冗余的资源维之后，同组内还可能有更多选项被支配。
    if !redundant.iter().any(|&dropped| dropped) {
        return Some(redundant);
    }
    for group in options {
        if expired(deadline) {
            return None;
        }
        let mut kept = Vec::new();
        for p in group.drain(..) {
            consider(&mut kept, p.cost, &p.used, p.mask);
        }
        kept.sort_unstable();
        *group = kept;
    }
    Some(redundant)
}

fn row(terms: Vec<(usize, f64)>, rhs: f64) -> Constraint {
    Constraint {
        name: None,
        expr: LinearExpr { terms },
        op: RelOp::Le,
        rhs,
    }
}

fn variable(p: &mut Problem, cost: i64) -> usize {
    let j = p.vars.len();
    p.vars.push(Variable {
        name: format!("x{j}"),
        var_type: VarType::Binary,
        lb: 0.0,
        ub: 1.0,
    });
    if cost != 0 {
        p.obj.terms.push((j, cost as f64));
    }
    j
}

pub(crate) fn reformulate(
    original: &Problem,
    tol: Tolerances,
    deadline: Option<Instant>,
) -> Option<Reformulation> {
    if ![tol.feasibility, tol.integrality, tol.zero]
        .iter()
        .all(|&v| v.is_finite() && (0.0..1.0).contains(&v))
        || original.sense != Sense::Minimize
        || original.vars.len() > 4096
        || original.constraints.len() > 4096
        || !original
            .vars
            .iter()
            .all(|v| v.var_type == VarType::Binary && v.lb == 0.0 && v.ub == 1.0)
        || !original.constraints.iter().all(|r| r.op == RelOp::Le)
        || expired(deadline)
    {
        return None;
    }
    let parts = Partition::detect(original)?;
    if !separated(original, tol) {
        return None;
    }
    let switches: Vec<_> = parts
        .links
        .iter()
        .map(|&i| switch_column(&original.constraints[i], &parts.owner))
        .collect();
    if switches.iter().all(Option::is_none) {
        return None;
    }
    let mut costs = vec![0i64; original.vars.len()];
    for &(j, c) in &original.obj.terms {
        costs[j] = c as i64;
    }
    let mut master = Problem::new(&original.name, Sense::Minimize);
    master.obj_offset = original.obj_offset;
    let mut links: Vec<_> = parts
        .links
        .iter()
        .map(|&i| row(Vec::new(), original.constraints[i].rhs))
        .collect();
    let mut switch_rows: Vec<Vec<Constraint>> = vec![Vec::new(); links.len()];
    let mut options_by_group = Vec::new();
    for g in 0..parts.groups.len() {
        options_by_group.push(choices(original, &parts, g, &switches, &costs, deadline)?);
    }
    let redundant = reduce_links(original, &parts, &switches, &mut options_by_group, deadline)?;
    let mut groups = Vec::new();
    for (cols, options) in parts.groups.iter().zip(options_by_group) {
        let start = master.vars.len();
        let base = &options[0];
        master.obj_offset += base.cost as f64;
        for option in &options[1..] {
            variable(&mut master, option.cost - base.cost);
        }
        if master.vars.len() + parts.external.len() > original.vars.len() {
            return None;
        }
        if options.len() > 1 {
            master.constraints.push(row(
                (start..master.vars.len()).map(|j| (j, 1.0)).collect(),
                1.0,
            ));
        }
        for (k, combined) in links.iter_mut().enumerate() {
            let mut activation = row(Vec::new(), 0.0);
            let target = if switches[k].is_some() {
                &mut activation
            } else {
                combined
            };
            target.rhs -= base.used[k] as f64;
            for (offset, option) in options[1..].iter().enumerate() {
                let change = option.used[k] - base.used[k];
                if change != 0 {
                    target.expr.terms.push((start + offset, change as f64));
                }
            }
            if switches[k].is_some() && (!activation.expr.terms.is_empty() || activation.rhs != 0.0)
            {
                switch_rows[k].push(activation);
            }
        }
        groups.push(GroupMap {
            columns: cols.clone(),
            start,
            masks: options.iter().map(|p| p.mask).collect(),
        });
    }
    let mut external_index = vec![usize::MAX; original.vars.len()];
    for &j in &parts.external {
        external_index[j] = variable(&mut master, costs[j]);
    }
    for (k, mut combined) in links.into_iter().enumerate() {
        if redundant[k] {
            continue;
        }
        if let Some(switch) = switches[k] {
            for mut constraint in switch_rows[k].drain(..) {
                constraint.expr.terms.push((external_index[switch], -1.0));
                master.constraints.push(constraint);
            }
        } else {
            for &(j, a) in &original.constraints[parts.links[k]].expr.terms {
                if parts.owner[j] == usize::MAX {
                    combined.expr.terms.push((external_index[j], a));
                }
            }
            master.constraints.push(combined);
        }
    }
    if !separated(&master, tol) || expired(deadline) {
        return None;
    }
    Some(Reformulation {
        problem: master,
        groups,
        external: parts.external,
        original_size: original.vars.len(),
        integer_tol: tol.integrality,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::branch_bound::SearchControl;
    use crate::model::SolveStatus;
    use crate::solver::{solve_controlled_with_incumbent, SolveOptions};
    use std::cell::RefCell;

    fn fixture() -> Problem {
        let mut p = Problem::new("分组开关", Sense::Minimize);
        for j in 0..6 {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        p.set_objective(
            vec![(0, -4.0), (1, 5.0), (2, -2.0), (3, 6.0), (4, 3.0), (5, 4.0)],
            -7.0,
        )
        .unwrap();
        for j in [0, 2] {
            p.constraints.push(row(vec![(j, 1.0), (j + 1, 1.0)], 1.0));
            p.constraints
                .push(row(vec![(j, -1.0), (j + 1, -1.0)], -1.0));
        }
        // 第一行要求启用外部开关，第二行的容量不足以分拆，必须原样计量。
        p.constraints
            .push(row(vec![(0, 2.0), (2, 3.0), (4, -5.0)], 0.0));
        p.constraints
            .push(row(vec![(0, 1.0), (2, 1.0), (5, -1.0)], 0.0));
        // 组内各恰选一，跨组最大活动度为 -2，仍可删除负右端的冗余行。
        p.constraints
            .push(row(vec![(0, -1.0), (1, -1.0), (2, -1.0), (3, -1.0)], -1.0));
        p
    }

    fn feasible(p: &Problem, x: &[f64]) -> bool {
        x.len() == p.vars.len()
            && x.iter().all(|&v| v == 0.0 || v == 1.0)
            && p.constraints.iter().all(|r| r.expr.eval(x) <= r.rhs)
    }

    fn assignments(n: usize) -> impl Iterator<Item = Vec<f64>> {
        (0..1usize << n).map(move |mask| (0..n).map(|j| ((mask >> j) & 1) as f64).collect())
    }

    #[test]
    fn projection_and_incumbent_match_original_exhaustion() {
        let tol = Tolerances::default();
        let mut p = fixture();
        for reverse in [false, true] {
            if reverse {
                p.constraints.reverse();
                p.vars.reverse();
                let n = p.vars.len();
                for expr in p
                    .constraints
                    .iter_mut()
                    .map(|r| &mut r.expr)
                    .chain([&mut p.obj])
                {
                    for (j, _) in &mut expr.terms {
                        *j = n - 1 - *j;
                    }
                    expr.terms.reverse();
                }
            }
            p.validate().unwrap();
            let work = reformulate(&p, tol, None).unwrap();
            work.problem.validate().unwrap();
            let truth = assignments(p.vars.len())
                .filter(|x| feasible(&p, x))
                .map(|x| (p.obj.eval(&x) + p.obj_offset) as i64)
                .min()
                .unwrap();
            let mut projected = None;
            for y in assignments(work.problem.vars.len()).filter(|y| feasible(&work.problem, y)) {
                let x = work.restore(&y).unwrap();
                assert!(feasible(&p, &x));
                let value = p.obj.eval(&x) + p.obj_offset;
                assert_eq!(value, work.problem.obj.eval(&y) + work.problem.obj_offset);
                projected = Some(projected.map_or(value as i64, |old: i64| old.min(value as i64)));
            }
            assert_eq!(projected, Some(truth));
            let seen = RefCell::new(Vec::new());
            let reply = solve_controlled_with_incumbent(
                &p,
                &SolveOptions::default(),
                &SearchControl::default(),
                None,
                &|_, value, x| {
                    assert!(feasible(&p, x));
                    assert_eq!(value, p.obj.eval(x) + p.obj_offset);
                    seen.borrow_mut().push(value);
                },
            )
            .unwrap();
            assert_eq!(reply.report.status, SolveStatus::Optimal);
            assert_eq!(reply.report.objective, Some(truth as f64));
            assert_eq!(reply.bound, reply.report.objective);
            assert!(feasible(&p, reply.report.values.as_ref().unwrap()));
            assert_eq!(seen.into_inner().last().copied(), reply.report.objective);
        }
    }

    #[test]
    fn uncertain_grid_tight_switch_or_elapsed_deadline_falls_back() {
        let mut p = fixture();
        let tol = Tolerances::default();
        assert!(reformulate(&p, tol, Some(Instant::now())).is_none());
        p.constraints[4].expr.terms[2].1 = -4.0;
        assert!(reformulate(&p, tol, None).is_none());
        p.constraints[4].expr.terms[2].1 = -5.0;
        p.constraints[0].rhs = 1.5;
        assert!(reformulate(&p, tol, None).is_none());
    }

    #[test]
    fn restore_rejects_double_selection_and_rounds_small_residuals() {
        let mut master = Problem::new("选择", Sense::Minimize);
        variable(&mut master, 1);
        variable(&mut master, 2);
        let work = Reformulation {
            problem: master,
            groups: vec![GroupMap {
                columns: vec![0, 1, 2],
                start: 0,
                masks: vec![1, 2, 4],
            }],
            external: vec![],
            original_size: 3,
            integer_tol: 1e-6,
        };
        assert_eq!(work.restore(&[0.0, 0.0]), Some(vec![1.0, 0.0, 0.0]));
        assert_eq!(
            work.restore(&[1.0 - 1e-8, -1e-8]),
            Some(vec![0.0, 1.0, 0.0])
        );
        assert!(work.restore(&[1.0, 1.0]).is_none());
        assert!(work.restore(&[0.9, 0.0]).is_none());
        assert!(work.restore(&[f64::NAN, 0.0]).is_none());
    }
}
