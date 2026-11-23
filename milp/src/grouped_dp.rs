//! 在互不相交的单选组之间传递剩余行界。处理一组选项后，两条行若在所有
//! 未处理组选项上的贡献完全相同，往后只需保留更紧的右端。于是同一组剩余行界
//! 的所有前缀具有相同的可行后缀，可仅保留其中费用最低的前缀及其回溯位置。

use crate::branch_bound::BbStats;
use crate::error::{MilpError, Result};
use crate::model::{Problem, RelOp, Sense, SolveStatus, Tolerances, VarType};
use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

const WIDTH: usize = 12;
const MAX_STATES: usize = 131_072;
const MAX_BACK: usize = 1_000_000;
type Key = [i16; WIDTH];

fn expired(deadline: Option<Instant>) -> bool {
    deadline.is_some_and(|end| Instant::now() >= end)
}

struct Group {
    // 普通二元列有“不选”和“选”两个选项。单选组的每个选项指定一列。
    options: Vec<Option<usize>>,
}

struct Row {
    // 每一组减去所有选项中最小的贡献，故任何剩余右端变负即可剪枝。
    coeff: Vec<Vec<i16>>,
    rhs: i16,
    first: usize,
}

struct Model {
    groups: Vec<Group>,
    rows: Vec<Row>,
    costs: Vec<i64>,
    sign: i64,
}

impl Model {
    fn detect(p: &Problem, tol: Tolerances, deadline: Option<Instant>) -> Option<Self> {
        let n = p.vars.len();
        if !(4..=512).contains(&n)
            || p.constraints.len() > 512
            || p.constraints.len() < 4
            || !tol.feasibility.is_finite()
            || !tol.integrality.is_finite()
            || tol.feasibility < 0.
            || tol.integrality < 0.
            || p.vars
                .iter()
                .any(|v| v.var_type != VarType::Binary || v.lb != 0. || v.ub != 1.)
        {
            return None;
        }
        // 对偶的单选上下限（也接受直接等式）。所有行随后仍接受整数格检查。
        let mut candidates = BTreeMap::<Vec<usize>, u8>::new();
        for r in &p.constraints {
            if !(2..=24).contains(&r.expr.terms.len())
                || r.rhs.abs() != 1.
                || !r.expr.terms.iter().all(|&(_, a)| a == r.rhs)
            {
                continue;
            }
            let flags = match (r.op, r.rhs > 0.) {
                (RelOp::Eq, _) => 3,
                (RelOp::Le, true) | (RelOp::Ge, false) => 1,
                _ => 2,
            };
            let cols = r.expr.terms.iter().map(|&(j, _)| j).collect::<Vec<_>>();
            *candidates.entry(cols).or_default() |= flags;
        }
        let mut sets: Vec<_> = candidates
            .into_iter()
            .filter_map(|(cols, flags)| (flags == 3).then_some(cols))
            .collect();
        if sets.len() < 2 {
            return None;
        }
        let mut grouped = vec![false; n];
        for g in &sets {
            for &j in g {
                if std::mem::replace(&mut grouped[j], true) {
                    return None;
                }
            }
        }
        let mut free: Vec<_> = (0..n).filter(|&j| !grouped[j]).collect();
        if free.is_empty() || sets.len() + free.len() > 100 {
            return None;
        }
        let sign = if p.sense == Sense::Minimize { 1 } else { -1 };
        let mut costs = vec![0; n];
        let mut total = p.obj_offset.abs();
        for &(j, a) in &p.obj.terms {
            if !(a.abs() < (1u64 << 48) as f64 && a == a.trunc()) {
                return None;
            }
            costs[j] = sign * (a as i64);
            total += a.abs();
        }
        if p.obj_offset != p.obj_offset.trunc() || !(total < (1u64 << 48) as f64) {
            return None;
        }

        let mut mixed = vec![Vec::new(); n];
        let mut support = vec![Vec::new(); n];
        for (i, r) in p.constraints.iter().enumerate() {
            if expired(deadline) {
                return None;
            }
            let norm: f64 = r.expr.terms.iter().map(|&(_, a)| a.abs()).sum();
            if !(r.rhs.abs() + norm <= 16_000.
                && r.rhs == r.rhs.trunc()
                && r.expr.terms.iter().all(|&(_, a)| a == a.trunc())
                && 10. * tol.feasibility + tol.integrality * norm < 1.
                && crate::solver::FINAL_VIOLATION_REL * (1. + r.rhs.abs() + norm) < 1.)
            {
                return None;
            }
            let includes_free = r.expr.terms.iter().any(|&(j, _)| !grouped[j]);
            for &(j, _) in &r.expr.terms {
                support[j].push(i);
                if includes_free && grouped[j] {
                    mixed[j].push(i);
                }
            }
        }
        // 先消去出现于较少混合行的组。逐期累积需求的较短后缀先消失，长行随即合并。
        // 相同支撑的组外列最后相邻处理，避免多个独立的整数编码同时展开。
        let weight = |g: &Vec<usize>| {
            let mut ids: Vec<_> = g.iter().flat_map(|&j| mixed[j].iter().copied()).collect();
            ids.sort_unstable();
            ids.dedup();
            ids.len()
        };
        sets.sort_by_cached_key(|g| (weight(g), g[0]));
        if weight(&sets[0]) == 0 || weight(&sets[0]) == weight(sets.last()?) {
            return None;
        }
        free.sort_by(|&a, &b| {
            support[a]
                .cmp(&support[b])
                .then_with(|| costs[b].abs().cmp(&costs[a].abs()))
                .then(a.cmp(&b))
        });
        let mut groups: Vec<_> = sets
            .into_iter()
            .map(|cols| Group {
                options: cols.into_iter().map(Some).collect(),
            })
            .collect();
        groups.extend(free.into_iter().map(|j| Group {
            options: vec![None, Some(j)],
        }));

        let mut rows = Vec::new();
        let mut dense = vec![0i16; n];
        for r in &p.constraints {
            if expired(deadline) {
                return None;
            }
            for &(j, a) in &r.expr.terms {
                dense[j] = a as i16;
            }
            // 等式分成两个上界，单选定义在组内自动成为恒真行。
            for mul in [1i16, -1] {
                if (mul == 1 && r.op == RelOp::Ge) || (mul == -1 && r.op == RelOp::Le) {
                    continue;
                }
                let mut rhs = r.rhs as i16 * mul;
                let mut coeff = Vec::with_capacity(groups.len());
                let mut first = groups.len();
                let mut max = 0;
                for (i, g) in groups.iter().enumerate() {
                    let mut a: Vec<_> = g
                        .options
                        .iter()
                        .map(|col| col.map_or(0, |j| mul * dense[j]))
                        .collect();
                    let min = *a.iter().min()?;
                    rhs -= min;
                    for v in &mut a {
                        *v -= min;
                    }
                    let span = *a.iter().max()?;
                    if span > 0 {
                        first = first.min(i);
                    }
                    max += span;
                    coeff.push(a);
                }
                if rhs < max {
                    rows.push(Row {
                        coeff,
                        rhs,
                        first: first.min(groups.len() - 1),
                    });
                }
            }
            for &(j, _) in &r.expr.terms {
                dense[j] = 0;
            }
        }
        Some(Self {
            groups,
            rows,
            costs,
            sign,
        })
    }
}

struct OptionData {
    feasible: bool,
    initial: Key,
    subtract: Key,
}
struct Step {
    dest: Vec<Option<usize>>,
    data: Vec<OptionData>,
}

// sig[t][r] 是去掉前 t 组后的完整行系数编号，0 表示该行已经没有变量。
// 从尾到头按（本组每个选项的系数，后缀编号）分配编号，等价性不依赖哈希碰撞。
fn compile(model: &Model, deadline: Option<Instant>) -> Option<Vec<Step>> {
    let stages = model.groups.len();
    let nr = model.rows.len();
    let mut sig = vec![vec![0usize; nr]; stages + 1];
    let mut span = vec![vec![0i16; nr]; stages + 1];
    for t in (0..stages).rev() {
        if expired(deadline) {
            return None;
        }
        let mut ids = BTreeMap::new();
        for (i, r) in model.rows.iter().enumerate() {
            let coeff = &r.coeff[t];
            span[t][i] = span[t + 1][i] + coeff.iter().max()?;
            if span[t][i] != 0 {
                let next = ids.len() + 1;
                sig[t][i] = *ids.entry((sig[t + 1][i], coeff)).or_insert(next);
            }
        }
    }
    let mut live = Vec::with_capacity(stages + 1);
    for (t, keys) in sig.iter().enumerate() {
        let mut reps = BTreeMap::new();
        for (i, r) in model.rows.iter().enumerate() {
            if r.first < t && keys[i] != 0 {
                reps.entry(keys[i]).or_insert(i);
            }
        }
        if reps.len() > WIDTH {
            return None;
        }
        live.push(reps.into_values().collect::<Vec<_>>());
    }
    let mut steps = Vec::with_capacity(stages);
    for (t, g) in model.groups.iter().enumerate() {
        if expired(deadline) {
            return None;
        }
        let child = &live[t + 1];
        let position = |i: usize| child.iter().position(|&j| sig[t + 1][i] == sig[t + 1][j]);
        let dest = live[t].iter().map(|&i| position(i)).collect();
        let mut data = Vec::new();
        for h in 0..g.options.len() {
            let mut initial = Key::default();
            for (b, &i) in child.iter().enumerate() {
                initial[b] = span[t + 1][i];
            }
            let mut feasible = true;
            for (i, r) in model.rows.iter().enumerate().filter(|(_, r)| r.first == t) {
                let left = r.rhs as i32 - r.coeff[t][h] as i32;
                if left < 0 {
                    feasible = false;
                    break;
                }
                if let Some(b) = position(i) {
                    initial[b] = initial[b].min(left as i16);
                }
            }
            let mut subtract = Key::default();
            for (a, &i) in live[t].iter().enumerate() {
                subtract[a] = model.rows[i].coeff[t][h];
            }
            data.push(OptionData {
                feasible,
                initial,
                subtract,
            });
        }
        steps.push(Step { dest, data });
    }
    Some(steps)
}

pub(crate) fn solve(
    problem: &Problem,
    tol: Tolerances,
    deadline: Option<Instant>,
    node_limit: u64,
    accept: &mut dyn FnMut(&[f64]),
) -> Result<Option<(SolveStatus, Option<f64>, BbStats)>> {
    let Some(m) = Model::detect(problem, tol, deadline) else {
        return Ok(None);
    };
    let Some(steps) = compile(&m, deadline) else {
        return Ok(None);
    };
    let option_cost = |j: Option<usize>| j.map_or(0, |j| m.costs[j]);
    let lower: i64 = m
        .groups
        .iter()
        .map(|g| g.options.iter().copied().map(option_cost).min().unwrap())
        .sum();
    let objective = |cost| problem.obj_offset + (m.sign * cost) as f64;
    let mut stats = BbStats::default();
    let stopped = |s| Ok(Some((SolveStatus::LimitReached, Some(objective(lower)), s)));
    let mut current = vec![(Key::default(), 0i64)];
    let mut history = Vec::new();
    let mut saved = 0;
    for (t, step) in steps.iter().enumerate() {
        let mut next = Vec::<(Key, i64)>::new();
        let mut back = Vec::new();
        // 不遍历哈希表；相同费用保留第一个来源。处理次序只由原列号和层序确定。
        let mut positions = HashMap::<Key, usize>::new();
        for (parent, (key, cost)) in current.iter().enumerate() {
            if stats.nodes >= node_limit || (parent % 64 == 0 && expired(deadline)) {
                return stopped(stats);
            }
            stats.nodes += 1;
            'options: for (h, data) in step.data.iter().enumerate() {
                if !data.feasible {
                    continue;
                }
                let mut target = data.initial;
                for (a, &dest) in step.dest.iter().enumerate() {
                    let left = key[a] - data.subtract[a];
                    if left < 0 {
                        continue 'options;
                    }
                    if let Some(b) = dest {
                        target[b] = target[b].min(left);
                    }
                }
                let cost = cost + option_cost(m.groups[t].options[h]);
                let from = (parent as u32, h as u8);
                match positions.entry(target) {
                    std::collections::hash_map::Entry::Occupied(slot) => {
                        let i = *slot.get();
                        if cost < next[i].1 {
                            next[i].1 = cost;
                            back[i] = from;
                        }
                    }
                    std::collections::hash_map::Entry::Vacant(slot) => {
                        if next.len() >= MAX_STATES || saved + next.len() >= MAX_BACK {
                            return stopped(stats);
                        }
                        slot.insert(next.len());
                        next.push((target, cost));
                        back.push(from);
                    }
                }
            }
        }
        if expired(deadline) {
            return stopped(stats);
        }
        if next.is_empty() {
            return Ok(Some((SolveStatus::Infeasible, None, stats)));
        }
        saved += back.len();
        history.push(back);
        current = next;
    }
    // 最后一层没有任何未处理变量，所有可行前缀的键必相同。
    let mut x = vec![0.; problem.vars.len()];
    let mut index = 0;
    for (t, back) in history.iter().enumerate().rev() {
        let (parent, h) = back[index];
        if let Some(j) = m.groups[t].options[h as usize] {
            x[j] = 1.;
        }
        index = parent as usize;
    }
    crate::solver::final_check(problem, &mut x)?;
    let expected = objective(current[0].1);
    if current.len() != 1 || expected != problem.obj_offset + problem.obj.eval(&x) {
        return Err(MilpError::Numerical(
            "单选组动态规划费用与原模型不一致".into(),
        ));
    }
    accept(&x);
    if expired(deadline) {
        return stopped(stats);
    }
    Ok(Some((SolveStatus::Optimal, Some(expected), stats)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instance::parse_instance;

    fn small() -> Problem {
        parse_instance("sense min\nvars\na binary\nb binary\nc binary\nd binary\nz binary\nobj\n3 b 6 c -2 d 1 z\ncons\n1 a 1 b = 1\n1 c 1 d >= 1\n-1 c -1 d >= -1\n1 a 1 z >= 1\n1 a 1 c 1 z >= 2\nend\n").unwrap()
    }

    #[test]
    fn suffix_merges_and_keeps_the_cheapest_prefix() {
        let mut p = small();
        p.obj_offset = 7.;
        for (sense, goal) in [(Sense::Minimize, 6.), (Sense::Maximize, 17.)] {
            p.sense = sense;
            let mut calls = 0;
            let (status, bound, _) = solve(&p, Default::default(), None, 1_000, &mut |x| {
                calls += 1;
                assert_eq!(p.obj.eval(x) + p.obj_offset, goal);
            })
            .unwrap()
            .unwrap();
            assert_eq!(
                (calls, status, bound),
                (1, SolveStatus::Optimal, Some(goal))
            );
        }
        let (status, bound, stats) = solve(&p, Default::default(), None, 1, &mut |_| {
            panic!("触限后没有完整赋值")
        })
        .unwrap()
        .unwrap();
        assert_eq!(status, SolveStatus::LimitReached);
        assert!(bound.unwrap() >= 17.);
        assert_eq!(stats.nodes, 1);
    }

    #[test]
    fn incompatible_suffixes_prove_infeasibility() {
        let mut p = small();
        // 两个新增单列等式与原有 a+z>=1 矛盾。
        p.add_constraint(Some("fixed_a"), vec![(0, 1.)], RelOp::Eq, 0.)
            .unwrap();
        p.add_constraint(Some("fixed_z"), vec![(4, 1.)], RelOp::Eq, 0.)
            .unwrap();
        let (status, bound, _) = solve(&p, Default::default(), None, 1_000, &mut |_| panic!())
            .unwrap()
            .unwrap();
        assert_eq!((status, bound), (SolveStatus::Infeasible, None));
    }

    #[test]
    fn unsafe_integer_grids_are_declined() {
        let mut p = small();
        p.sense = Sense::Maximize;
        p.obj.terms[0].1 = -1.0e100;
        assert!(Model::detect(&p, Default::default(), None).is_none());
        let mut q = small();
        q.constraints[3].rhs = 0.999999999;
        assert!(Model::detect(&q, Default::default(), None).is_none());
    }
}
