//! 用等式把一列从其唯一另一行消去，使现有私有余量贪心能给出候选。
//! 此副本仅用于找解，所有候选均返回原工作模型验收。

use crate::model::{Problem, RelOp, Sense, Tolerances, VarType};

use std::{collections::BTreeMap, time::Instant};

pub(super) fn projected(
    p: &Problem,
    t: &Tolerances,
    end: Option<Instant>,
) -> Option<(Problem, Vec<f64>)> {
    let (n, m) = (p.vars.len(), p.constraints.len());
    let expired = || end.is_some_and(|e| Instant::now() >= e);
    if n < 4096 || n / 8 < m || p.sense != Sense::Minimize || expired() {
        return None;
    }
    let mut cols = vec![Vec::new(); n];
    let mut costs = vec![0.; n];
    for &(j, c) in &p.obj.terms {
        costs[j] = c;
    }
    let mut nnz = 0usize;
    for (i, r) in p.constraints.iter().enumerate() {
        if expired() {
            return None;
        }
        nnz += r.expr.terms.len();
        if nnz > 2_000_000 {
            return None;
        }
        for &(j, a) in &r.expr.terms {
            if cols[j].len() < 3 {
                cols[j].push((i, a));
            }
        }
    }
    // 每个等式选择至多一个代价非负的无界余量；用于消元的另一行也只用一次。
    let mut choices = vec![None; m];
    for (j, v) in p.vars.iter().enumerate() {
        if v.ub != f64::INFINITY || costs[j] <= 0. || cols[j].len() != 2 {
            continue;
        }
        let [(mut i, mut a), (mut k, mut b)] = cols[j][..] else {
            continue;
        };
        if p.constraints[k].op == RelOp::Eq {
            std::mem::swap(&mut i, &mut k);
            std::mem::swap(&mut a, &mut b);
        }
        if p.constraints[i].op != RelOp::Eq
            || p.constraints[k].op == RelOp::Eq
            || a == 0.
            || b == 0.
        {
            continue;
        }
        let ratio = b / a;
        if !ratio.is_finite() || ratio.abs() > 1e3 || ratio.abs() < 1e-3 {
            continue;
        }
        let score = costs[j] / a.abs();
        if choices[i].is_none_or(|(s, _, _, _)| score > s) {
            choices[i] = Some((score, k, ratio, j));
        }
    }
    if choices.iter().flatten().count() < m.div_ceil(8) {
        return None;
    }
    let mut q = p.clone();
    let mut used = vec![false; m];
    for (i, choice) in choices.into_iter().enumerate() {
        if expired() {
            return None;
        }
        let Some((_, k, f, j)) = choice else {
            continue;
        };
        if used[k] {
            continue;
        }
        used[k] = true;
        let mut row: BTreeMap<usize, f64> = p.constraints[k].expr.terms.iter().copied().collect();
        for &(h, a) in &p.constraints[i].expr.terms {
            *row.entry(h).or_default() -= f * a;
        }
        // 仅用于生成候选；计算中的舍入余项在原行复核。
        row.remove(&j);
        q.constraints[k].expr.terms = row.into_iter().filter(|(_, a)| a.abs() > t.zero).collect();
        q.constraints[k].rhs = p.constraints[k].rhs - f * p.constraints[i].rhs;
    }
    let mut offset = vec![0.; n];
    for (j, v) in q.vars.iter_mut().enumerate() {
        if !v.lb.is_finite() || (v.var_type != VarType::Continuous && v.lb != v.lb.round()) {
            return None;
        }
        offset[j] = v.lb;
        v.ub -= v.lb;
        v.lb = 0.;
    }
    if offset.iter().any(|&a| a != 0.) {
        for r in &mut q.constraints {
            r.rhs -= r.expr.eval(&offset);
        }
        q.obj_offset += q.obj.eval(&offset);
    }
    Some((q, offset))
}

/// 保留非负变量界对应的半行，把独占余量的费用代入其他列。此变换只用来找候选，
/// 被消去的整数余量还原后若有小数，完整模型的可行性检查会拒绝该候选。
pub(super) fn eliminate(
    p: &mut Problem,
    deadline: Instant,
) -> Option<(Vec<usize>, Vec<(usize, usize, f64)>)> {
    let mut count = vec![0; p.vars.len()];
    for r in &p.constraints {
        for &(j, _) in &r.expr.terms {
            count[j] += 1;
        }
    }
    let mut cost = vec![0.; p.vars.len()];
    for &(j, c) in &p.obj.terms {
        cost[j] = c;
    }
    let mut fill = Vec::new();
    for (i, r) in p.constraints.iter_mut().enumerate() {
        if Instant::now() >= deadline {
            return None;
        }
        if r.op != RelOp::Eq {
            continue;
        }
        let Some(&(j, a)) = r
            .expr
            .terms
            .iter()
            .filter(|&&(j, _)| p.vars[j].ub == f64::INFINITY && count[j] == 1 && cost[j] >= 0.)
            .max_by(|&&(j, a), &&(h, b)| {
                (cost[j] / a.abs())
                    .total_cmp(&(cost[h] / b.abs()))
                    .then_with(|| h.cmp(&j))
            })
        else {
            continue;
        };
        let rate = cost[j] / a;
        if !rate.is_finite() {
            return None;
        }
        p.obj_offset += rate * r.rhs;
        cost[j] = 0.;
        p.vars[j].ub = 0.;
        r.expr.terms.retain(|&(h, _)| h != j);
        for &(h, b) in &r.expr.terms {
            cost[h] -= rate * b;
        }
        r.op = if a > 0. { RelOp::Le } else { RelOp::Ge };
        fill.push((i, j, a));
    }
    if cost.iter().any(|c| !c.is_finite()) || !p.obj_offset.is_finite() {
        return None;
    }
    p.obj.terms = cost
        .into_iter()
        .enumerate()
        .filter(|(_, v)| *v != 0.)
        .collect();
    Some((count, fill))
}
