//! 一个长覆盖行要求至少一个二元变量取一，其余局部约束可以穷举成小表。
//! 逐个尝试只有一个覆盖列取一的起点，以低费用选项补全原解。这是启发式；
//! 下界只假设某一列取一，保留其他覆盖列的全部取值，故也适用于多列同时取一。

use crate::model::{Problem, RelOp, Sense, Tolerances, VarType};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::time::Instant;

const WIDTH: usize = 10;

fn expired(deadline: Option<Instant>) -> bool {
    deadline.is_some_and(|end| Instant::now() >= end)
}

struct Table {
    cols: Vec<usize>,
    // 依次按本表费用、取一列数和掩码排序，第一项匹配就是本步的贪心选择。
    options: Vec<(i64, u16)>,
}

impl Table {
    fn enumerate(
        &mut self,
        p: &Problem,
        ids: &[usize],
        cost: &[i64],
        work: &mut usize,
        deadline: Option<Instant>,
    ) -> Option<()> {
        let mut rows = Vec::with_capacity(ids.len());
        for &i in ids {
            let r = &p.constraints[i];
            let terms: Vec<_> = r
                .expr
                .terms
                .iter()
                .map(|&(j, a)| (1u16 << self.cols.binary_search(&j).unwrap(), a as i32))
                .collect();
            *work = work.checked_sub(terms.len() << self.cols.len())?;
            rows.push((terms, r.op, r.rhs as i32));
        }
        for mask in 0u16..1u16 << self.cols.len() {
            if mask % 64 == 0 && expired(deadline) {
                return None;
            }
            let ok = rows.iter().all(|(terms, op, rhs)| {
                let lhs: i32 = terms
                    .iter()
                    .filter_map(|&(b, a)| (mask & b != 0).then_some(a))
                    .sum();
                match op {
                    RelOp::Le => lhs <= *rhs,
                    RelOp::Ge => lhs >= *rhs,
                    RelOp::Eq => lhs == *rhs,
                }
            });
            if ok {
                let price = self
                    .cols
                    .iter()
                    .enumerate()
                    .filter_map(|(a, &j)| (mask & (1 << a) != 0).then_some(cost[j]))
                    .sum();
                self.options.push((price, mask));
            }
        }
        self.options
            .sort_unstable_by_key(|&(price, mask)| (price, mask.count_ones(), mask));
        Some(())
    }
}

struct Tables {
    anchor: Vec<usize>,
    tables: Vec<Table>,
    order: Vec<usize>,
    lower: i64,
}

impl Tables {
    fn compile(p: &Problem, tol: Tolerances, deadline: Option<Instant>) -> Option<Self> {
        let n = p.vars.len();
        if !(WIDTH + 1..=8192).contains(&n)
            || p.constraints.len() > 16_384
            || p.sense != Sense::Minimize
            || !(tol.feasibility >= 0. && tol.integrality >= 0.)
            || !p
                .vars
                .iter()
                .all(|v| v.var_type == VarType::Binary && v.lb == 0. && v.ub == 1.)
            || !p
                .obj
                .terms
                .iter()
                .all(|&(_, c)| (0.0..=1_000_000.0).contains(&c) && c == c.trunc())
            || p.obj_offset.abs() + n as f64 * 1e6 >= (1u64 << 48) as f64
        {
            return None;
        }
        let mut wide = p.constraints.iter().filter(|r| r.expr.terms.len() > WIDTH);
        let anchor_row = wide.next()?;
        if wide.next().is_some()
            || anchor_row.expr.terms.len() > 512
            || !matches!(
                (anchor_row.op, anchor_row.rhs),
                (RelOp::Ge, 1.) | (RelOp::Le, -1.)
            )
            || !anchor_row
                .expr
                .terms
                .iter()
                .all(|&(_, a)| a == anchor_row.rhs)
        {
            return None;
        }
        let mut local_rows = Vec::new();
        for (i, r) in p.constraints.iter().enumerate() {
            let norm = r.expr.terms.iter().map(|&(_, a)| a.abs()).sum::<f64>();
            if r.expr.terms.is_empty()
                || !(norm + r.rhs.abs() <= 32_000.
                    && r.rhs == r.rhs.trunc()
                    && r.expr.terms.iter().all(|&(_, a)| a == a.trunc())
                    && 10. * tol.feasibility + tol.integrality * norm < 1.
                    && crate::solver::FINAL_VIOLATION_REL * (1. + norm + r.rhs.abs()) < 1.)
            {
                return None;
            }
            if r.expr.terms.len() <= WIDTH {
                local_rows.push(i);
            }
        }
        let mut cost = vec![0i64; n];
        for &(j, c) in &p.obj.terms {
            cost[j] = c as i64;
        }
        // 含于同一列集的短行合成一张表。每行只分配一次。
        local_rows.sort_unstable_by_key(|&i| (Reverse(p.constraints[i].expr.terms.len()), i));
        let mut posts: Vec<Vec<usize>> = vec![vec![]; n];
        let mut tables: Vec<Table> = Vec::new();
        let mut owners: Vec<Vec<usize>> = Vec::new();
        for i in local_rows {
            if expired(deadline) {
                return None;
            }
            let terms = &p.constraints[i].expr.terms;
            let &(rare, _) = terms.iter().min_by_key(|&&(j, _)| posts[j].len())?;
            let owner = posts[rare].iter().copied().find(|&k| {
                terms
                    .iter()
                    .all(|&(j, _)| tables[k].cols.binary_search(&j).is_ok())
            });
            let k = owner.unwrap_or_else(|| {
                let k = tables.len();
                for &(j, _) in terms {
                    posts[j].push(k);
                }
                tables.push(Table {
                    cols: terms.iter().map(|&(j, _)| j).collect(),
                    options: vec![],
                });
                owners.push(Vec::new());
                k
            });
            owners[k].push(i);
        }
        let mut work = 64_000_000usize;
        let mut conditional = vec![0i64; n];
        for (table, ids) in tables.iter_mut().zip(owners) {
            if expired(deadline) {
                return None;
            }
            table.enumerate(p, &ids, &cost, &mut work, deadline)?;
            for (a, &j) in table.cols.iter().enumerate() {
                // 同一列可能出现在多张表，只取最紧的单表下界，避免重复计费。
                let min = table
                    .options
                    .iter()
                    .find(|&&(_, mask)| mask & (1 << a) != 0)
                    .map_or(i64::MAX, |&(price, _)| price);
                conditional[j] = conditional[j].max(min);
            }
        }
        let anchor: Vec<_> = anchor_row.expr.terms.iter().map(|&(j, _)| j).collect();
        let lower = anchor.iter().map(|&j| conditional[j].max(cost[j])).min()?;
        if lower == i64::MAX {
            return None;
        }
        let order = completion_order(&anchor, &tables, &posts, deadline)?;
        Some(Self {
            anchor,
            tables,
            order,
            lower,
        })
    }

    fn complete(
        &self,
        n: usize,
        selected: usize,
        deadline: Option<Instant>,
        work: &mut usize,
    ) -> Option<Vec<f64>> {
        let mut bits = vec![-1i8; n];
        for &j in &self.anchor {
            bits[j] = i8::from(j == selected);
        }
        for (step, &t) in self.order.iter().enumerate() {
            if step % 64 == 0 && expired(deadline) {
                return None;
            }
            let tab = &self.tables[t];
            let mut used = 0;
            let mut ones = 0;
            for (a, &j) in tab.cols.iter().enumerate() {
                if bits[j] >= 0 {
                    used |= 1 << a;
                }
                if bits[j] == 1 {
                    ones |= 1 << a;
                }
            }
            let count = tab.options.len().min(*work);
            let hit = tab.options[..count]
                .iter()
                .position(|&(_, mask)| mask & used == ones);
            *work -= hit.map_or(count, |i| i + 1);
            let (_, mask) = tab.options[hit?];
            for (a, &j) in tab.cols.iter().enumerate() {
                bits[j] = ((mask >> a) & 1) as i8;
            }
        }
        Some(bits.into_iter().map(|v| f64::from(v.max(0))).collect())
    }
}

fn completion_order(
    anchor: &[usize],
    tables: &[Table],
    posts: &[Vec<usize>],
    deadline: Option<Instant>,
) -> Option<Vec<usize>> {
    // 先补已知列占比高的表，避免过早猜测尚未与起点连通的变量。
    let mut known = vec![false; posts.len()];
    for &j in anchor {
        known[j] = true;
    }
    let mut touched: Vec<_> = tables
        .iter()
        .map(|t| t.cols.iter().filter(|&&j| known[j]).count())
        .collect();
    let mut done = vec![false; tables.len()];
    // 2520 是行宽 1 至 10 的公倍数，可以精确比较已知列的比例。
    let rank = |i: usize, count: usize| (count * 2520 / tables[i].cols.len(), Reverse(i));
    let mut queue: BinaryHeap<_> = touched
        .iter()
        .enumerate()
        .map(|(i, &v)| rank(i, v))
        .collect();
    let mut order = Vec::with_capacity(tables.len());
    while let Some((_, Reverse(i))) = queue.pop() {
        if expired(deadline) {
            return None;
        }
        if done[i] {
            continue;
        }
        done[i] = true;
        order.push(i);
        for &j in &tables[i].cols {
            if std::mem::replace(&mut known[j], true) {
                continue;
            }
            for &k in &posts[j] {
                if !done[k] {
                    touched[k] += 1;
                    queue.push(rank(k, touched[k]));
                }
            }
        }
    }
    Some(order)
}

pub(crate) fn seed(
    p: &Problem,
    tol: Tolerances,
    deadline: Option<Instant>,
    limit: u64,
    nodes: &mut u64,
    mut publish: impl FnMut(Vec<f64>),
) -> Option<f64> {
    let m = Tables::compile(p, tol, deadline)?;
    // 单列起点只是寻找可行解，整批尝试共用查表额度，余下预算交回通用求解器。
    let mut work = 16_000_000;
    let mut best = f64::INFINITY;
    for &j in &m.anchor {
        if *nodes >= limit || work == 0 || expired(deadline) {
            break;
        }
        *nodes += 1;
        if let Some(mut x) = m.complete(p.vars.len(), j, deadline, &mut work) {
            let cost = p.obj.eval(&x) + p.obj_offset;
            if cost < best && crate::solver::final_check(p, &mut x).is_ok() {
                publish(x);
                best = cost;
                if best == p.obj_offset + m.lower as f64 {
                    break;
                }
            }
        }
    }
    Some(p.obj_offset + m.lower as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn several_active_columns_can_beat_every_singleton() {
        let mut p = Problem::new("tables", Sense::Minimize);
        for j in 0..12 {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        p.set_objective(
            (0..12)
                .map(|j| (j, if j == 11 { 5. } else { 1. }))
                .collect(),
            -0.5,
        )
        .unwrap();
        p.add_constraint(None, (0..11).map(|j| (j, 1.)).collect(), RelOp::Ge, 1.)
            .unwrap();
        for (j, k) in [(0, 1), (1, 0)] {
            p.add_constraint(None, vec![(j, 1.), (k, -1.), (11, -1.)], RelOp::Le, 0.)
                .unwrap();
        }
        for j in 2..11 {
            p.add_constraint(None, vec![(j, 1.), (11, -1.)], RelOp::Le, 0.)
                .unwrap();
        }
        let mut nodes = 0;
        let mut published = vec![];
        let bound = seed(&p, Default::default(), None, 100, &mut nodes, |x| {
            published.push(p.obj.eval(&x) + p.obj_offset);
        });
        // 两列同时取一时，费用为 1.5，单列起点至少需要 5.5。
        assert_eq!(bound, Some(1.5));
        assert_eq!(published, vec![5.5]);
        assert_eq!(nodes, 11);
        let mut nodes = 0;
        assert_eq!(
            seed(&p, Default::default(), None, 1, &mut nodes, |_| {}),
            bound
        );
        assert_eq!(nodes, 1);
        p.constraints[1].rhs = 0.25;
        assert!(Tables::compile(&p, Default::default(), None).is_none());
    }
}
