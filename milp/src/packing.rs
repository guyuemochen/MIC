//! 少量整数容量约束的二元装载。把每件物品的收益分给使用的资源，每行独立做一次
//! 一维背包。任一原可行解在这些背包中均可行，且分摊的收益之和等于原收益，
//! 所以各行最优值的和是可证明的上界。改变分摊只影响界的紧度。
//! 搜索时按剩余容量查后缀背包表，尾部少量物品以位集求交，直接还原最优组合。

use crate::branch_bound::BbStats;
use crate::error::{MilpError, Result};
use crate::model::{Problem, RelOp, Sense, SolveStatus, Tolerances, VarType};
use crate::simplex::{Bounds, LpStatus, Simplex};
use std::cmp::Reverse;
use std::time::Instant;

#[path = "packing_core.rs"]
mod packing_core;

const ROWS: usize = 8;
const SCALE: i64 = 10_000;
const MAX_COLUMNS: usize = 512;
const MAX_CELLS: usize = 1 << 22;
// 收益、单行表、各行边际之和及中间差值的绝对值均需能用 i64 表示。
const MAX_SHARE: i64 = i64::MAX / (16 * ROWS * MAX_COLUMNS) as i64;

fn expired(deadline: Option<Instant>) -> bool {
    deadline.is_some_and(|end| Instant::now() >= end)
}

#[derive(Clone)]
struct Item {
    column: usize,
    value: i64,
    weight: [usize; ROWS],
    share: [i64; ROWS],
    bias: i64,
}

struct Packing {
    items: Vec<Item>,
    capacity: Vec<usize>,
    fixed: Vec<f64>,
    fixed_value: i64,
}

impl Packing {
    fn detect(p: &Problem, tol: Tolerances) -> Option<Self> {
        let n = p.vars.len();
        let m = p.constraints.len();
        if !(2..=ROWS).contains(&m)
            || !(64..=MAX_COLUMNS).contains(&n)
            || p.sense != Sense::Maximize
            || p.obj_offset.abs() > (1_u64 << 40) as f64
            || p.obj_offset.fract() != 0.0
            || p.vars
                .iter()
                .any(|v| v.var_type != VarType::Binary || v.lb != 0.0 || v.ub != 1.0)
        {
            return None;
        }
        let mut items: Vec<_> = (0..n)
            .map(|column| Item {
                column,
                value: 0,
                weight: [0; ROWS],
                share: [0; ROWS],
                bias: 0,
            })
            .collect();
        for &(j, c) in &p.obj.terms {
            if c.fract() != 0.0 || !(0.0..=32767.0).contains(&c) {
                return None;
            }
            items[j].value = c as i64;
        }
        let mut capacity = Vec::with_capacity(m);
        for (k, row) in p.constraints.iter().enumerate() {
            if row.op != RelOp::Le || row.rhs.fract() != 0.0 || !(0.0..=4096.0).contains(&row.rhs) {
                return None;
            }
            let mut norm = 0.;
            for &(j, a) in &row.expr.terms {
                if a.fract() != 0. || !(1.0..=32767.0).contains(&a) {
                    return None;
                }
                items[j].weight[k] = a as usize;
                norm += a;
            }
            if !(10.0 * tol.feasibility + tol.integrality * norm < 1.0
                && crate::solver::FINAL_VIOLATION_REL * (1.0 + row.rhs + norm) < 1.0)
            {
                return None;
            }
            capacity.push(row.rhs as usize);
        }
        if capacity.iter().map(|b| b + 1).sum::<usize>() * (n + 1) > MAX_CELLS {
            return None;
        }
        let mut fixed = vec![0.; n];
        let mut fixed_value = 0;
        items.retain(|item| {
            if item.value == 0
                || capacity
                    .iter()
                    .enumerate()
                    .any(|(k, &b)| item.weight[k] > b)
            {
                return false;
            }
            if item.weight.iter().all(|&a| a == 0) {
                fixed[item.column] = 1.;
                fixed_value += item.value;
                return false;
            }
            true
        });
        Some(Self {
            items,
            capacity,
            fixed,
            fixed_value,
        })
    }

    // 沿当前分摊偏向的次序装入物品，每次按原容量检查。
    fn greedy(&self) -> (i64, Vec<f64>) {
        let mut base = self.fixed.clone();
        let mut remaining = self.capacity.clone();
        let mut best = self.fixed_value;
        let mut order: Vec<_> = self.items.iter().collect();
        order.sort_unstable_by_key(|it| (Reverse(it.bias), it.column));
        for it in order {
            if remaining
                .iter()
                .enumerate()
                .all(|(k, &b)| b >= it.weight[k])
            {
                for (k, b) in remaining.iter_mut().enumerate() {
                    *b -= it.weight[k];
                }
                best += it.value;
                base[it.column] = 1.;
            }
        }
        (best, base)
    }

    fn initial_shares(&mut self, prices: &[i64]) {
        for item in &mut self.items {
            let mut remaining = item.value * SCALE;
            let mut uses = 0;
            for (k, &price) in prices.iter().enumerate() {
                item.share[k] = price * item.weight[k] as i64;
                remaining -= item.share[k];
                if item.weight[k] > 0 {
                    uses += 1;
                }
            }
            item.bias = remaining;
            for k in 0..prices.len() {
                if item.weight[k] > 0 {
                    let part = remaining / uses;
                    item.share[k] += part;
                    remaining -= part;
                    uses -= 1;
                }
            }
        }
        self.items
            .sort_unstable_by_key(|it| (Reverse(it.bias.abs()), it.column));
    }
}

struct RowTable {
    stride: usize,
    values: Vec<i64>,
}

impl RowTable {
    fn layer(&self, i: usize) -> &[i64] {
        &self.values[i * self.stride..(i + 1) * self.stride]
    }
    fn fill(&mut self, items: &[Item], k: usize, deadline: Option<Instant>) -> bool {
        for (i, it) in items.iter().enumerate().rev() {
            if i % 32 == 0 && expired(deadline) {
                return false;
            }
            let (left, right) = self.values.split_at_mut((i + 1) * self.stride);
            let row = &mut left[i * self.stride..];
            let next = &right[..self.stride];
            row.copy_from_slice(next);
            let w = it.weight[k];
            if it.share[k] > 0 && w > 0 {
                for (v, &tail) in row[w..].iter_mut().zip(next) {
                    *v = (*v).max(it.share[k] + tail);
                }
            }
        }
        true
    }
}

fn fill_tables(tables: &mut [RowTable], items: &[Item], deadline: Option<Instant>) -> bool {
    tables
        .iter_mut()
        .enumerate()
        .all(|(k, d)| d.fill(items, k, deadline))
}

// 每行计算该物品取一和取零的最优值差。把这些差的平均值分回各行，
// 余数按行下标分配，逐件逐次保持整数收益总和。前缀使用更新后的分摊。
fn balance(items: &mut [Item], tables: &[RowTable], deadline: Option<Instant>) -> bool {
    let mut prefix: Vec<_> = tables.iter().map(|t| vec![0_i64; t.stride]).collect();
    for (i, item) in items.iter_mut().enumerate() {
        if i % 16 == 0 && expired(deadline) {
            return false;
        }
        let mut delta = [0; ROWS];
        let mut sum = 0;
        let mut count = 0;
        for (k, table) in tables.iter().enumerate() {
            let w = item.weight[k];
            if w == 0 {
                continue;
            }
            let pre = &prefix[k];
            let suffix = table.layer(i + 1);
            let best = |limit: usize| {
                pre[..=limit]
                    .iter()
                    .zip(suffix[..=limit].iter().rev())
                    .map(|(a, b)| a + b)
                    .max()
                    .unwrap_or(0)
            };
            delta[k] = item.share[k] + best(table.stride - 1 - w) - best(table.stride - 1);
            sum += delta[k];
            count += 1;
        }
        item.bias = sum;
        let avg = sum / count;
        let mut extra = sum % count;
        for (k, pre) in prefix.iter_mut().enumerate() {
            let w = item.weight[k];
            if w == 0 {
                continue;
            }
            let adj = extra.signum();
            extra -= adj;
            item.share[k] += avg + adj - delta[k];
            if item.share[k].abs() > MAX_SHARE {
                return false;
            }
            if item.share[k] > 0 {
                for b in (w..pre.len()).rev() {
                    pre[b] = pre[b].max(item.share[k] + pre[b - w]);
                }
            }
        }
    }
    true
}

struct Tail {
    start: usize,
    // 按原收益降序，收益相同时按组合编号。
    options: Vec<(i64, u16)>,
    // 第 k 行按容量保存可行组合的位集，组合顺序与 options 相同。
    bits: Vec<Vec<u64>>,
    limits: Vec<usize>,
    words: usize,
}

impl Tail {
    fn new(items: &[Item], capacity: &[usize], deadline: Option<Instant>) -> Option<Self> {
        // 限制位表内存；若行容量较大则把更少的尾部列交给查表。
        let mut len = items.len().min(15);
        while ((1_usize << len).div_ceil(64)) * capacity.iter().map(|b| b + 1).sum::<usize>()
            > MAX_CELLS
        {
            len -= 1;
        }
        let start = items.len() - len;
        let tail = &items[start..];
        let count = 1_usize << len;
        let words = count.div_ceil(64);
        let mut profit = vec![0; count];
        for mask in 1..count {
            profit[mask] = profit[mask & (mask - 1)] + tail[mask.trailing_zeros() as usize].value;
        }
        let mut options: Vec<_> = profit
            .into_iter()
            .enumerate()
            .map(|(mask, v)| (v, mask as u16))
            .collect();
        options.sort_unstable_by_key(|&(v, mask)| (Reverse(v), mask));
        let mut bits = Vec::new();
        let mut limits = Vec::new();
        let mut loads = vec![0; count];
        for (k, &b) in capacity.iter().enumerate() {
            if expired(deadline) {
                return None;
            }
            let limit = b.min(tail.iter().map(|it| it.weight[k]).sum());
            limits.push(limit);
            let mut row = vec![0u64; (limit + 1) * words];
            for mask in 1..count {
                loads[mask] =
                    loads[mask & (mask - 1)] + tail[mask.trailing_zeros() as usize].weight[k];
            }
            for (i, &(_, mask)) in options.iter().enumerate() {
                let w = loads[mask as usize];
                if w <= limit {
                    row[w * words + i / 64] |= 1u64 << (i % 64);
                }
            }
            for b in 1..=limit {
                let (lo, hi) = row.split_at_mut(b * words);
                for (v, &prev) in hi[..words].iter_mut().zip(&lo[(b - 1) * words..]) {
                    *v |= prev;
                }
            }
            bits.push(row);
        }
        Some(Self {
            start,
            options,
            bits,
            limits,
            words,
        })
    }
    fn best(&self, rem: &[usize], threshold: i64) -> Option<(i64, u16)> {
        let end = self.options.partition_point(|&(v, _)| v >= threshold);
        let mut base = [0; ROWS];
        for (k, &b) in rem.iter().enumerate() {
            base[k] = b.min(self.limits[k]) * self.words;
        }
        for word in 0..end.div_ceil(64) {
            let mut candidates = u64::MAX;
            for (k, row) in self.bits.iter().enumerate() {
                candidates &= row[base[k] + word];
                if candidates == 0 {
                    break;
                }
            }
            if candidates != 0 {
                let idx = word * 64 + candidates.trailing_zeros() as usize;
                return (idx < end).then(|| self.options[idx]);
            }
        }
        None
    }
}

struct Search<'a> {
    items: &'a [Item],
    tables: &'a [RowTable],
    tail: Tail,
    deadline: Option<Instant>,
    limit: u64,
    nodes: u64,
    stopped: bool,
    rem: Vec<usize>,
    values: Vec<f64>,
    best: i64,
    accept: &'a mut dyn FnMut(&[f64]),
}

impl Search<'_> {
    fn bound(&self, i: usize, profit: i64) -> i64 {
        let sum: i64 = self
            .tables
            .iter()
            .zip(&self.rem)
            .map(|(table, &b)| table.values[i * table.stride + b])
            .sum();
        profit + sum / SCALE
    }
    fn visit(&mut self, i: usize, profit: i64) {
        if self.nodes >= self.limit || (self.nodes % 256 == 0 && expired(self.deadline)) {
            self.stopped = true;
            return;
        }
        self.nodes += 1;
        if profit > self.best {
            self.best = profit;
            (self.accept)(&self.values);
        }
        if self.bound(i, profit) <= self.best {
            return;
        }
        if i == self.tail.start {
            if let Some((extra, mask)) = self.tail.best(&self.rem, self.best + 1 - profit) {
                let tail = &self.items[i..];
                for (j, it) in tail.iter().enumerate() {
                    self.values[it.column] = f64::from((mask >> j) & 1);
                }
                self.best = profit + extra;
                (self.accept)(&self.values);
                for it in tail {
                    self.values[it.column] = 0.;
                }
            }
            return;
        }
        let it = &self.items[i];
        let drop = self.bound(i + 1, profit);
        let take = if self.rem.iter().enumerate().all(|(k, &b)| it.weight[k] <= b) {
            let sum: i64 = self
                .tables
                .iter()
                .enumerate()
                .map(|(k, t)| t.layer(i + 1)[self.rem[k] - it.weight[k]])
                .sum();
            profit + it.value + sum / SCALE
        } else {
            -1
        };
        let prefer_take = take > drop || (take == drop && it.bias > 0);
        for take_first in [prefer_take, !prefer_take] {
            if self.stopped {
                break;
            }
            if take_first && take > self.best {
                self.values[it.column] = 1.;
                for (k, r) in self.rem.iter_mut().enumerate() {
                    *r -= it.weight[k];
                }
                self.visit(i + 1, profit + it.value);
                for (k, r) in self.rem.iter_mut().enumerate() {
                    *r += it.weight[k];
                }
                self.values[it.column] = 0.;
            } else if !take_first && drop > self.best {
                self.visit(i + 1, profit);
            }
        }
    }
}

pub(crate) fn solve(
    p: &Problem,
    tol: Tolerances,
    deadline: Option<Instant>,
    limit: u64,
    accept: &mut dyn FnMut(&[f64]),
) -> Result<Option<(SolveStatus, Option<f64>, BbStats)>> {
    let Some(mut model) = Packing::detect(p, tol) else {
        return packing_core::solve(p, tol, deadline, limit, accept);
    };
    let mut stats = BbStats::default();
    let limited = |stats| Ok(Some((SolveStatus::LimitReached, None, stats)));
    if expired(deadline) {
        return limited(stats);
    }
    let mut lp = Simplex::new(p, tol)?;
    lp.set_time_deadline(deadline);
    let sol = lp.solve(&Bounds::from_problem(p));
    stats.lp_iterations = lp.iterations();
    if expired(deadline) {
        return limited(stats);
    }
    if matches!(sol.status, LpStatus::IterationLimit | LpStatus::TimeLimit) {
        return limited(stats);
    }
    if sol.status != LpStatus::Optimal {
        return Ok(None);
    }
    // LP 仅提供起始资源价格。取整后的分摊及其上界完全用整数重新计算。
    let prices: Vec<_> = sol
        .duals
        .iter()
        .map(|&d| ((-d * SCALE as f64).round().clamp(0., 32767. * SCALE as f64)) as i64)
        .collect();
    model.initial_shares(&prices);
    let (mut best, first) = model.greedy();
    if p.constraints
        .iter()
        .any(|row| row.expr.eval(&first) > row.rhs)
    {
        return Err(MilpError::Numerical("装载初解未通过原约束核验".into()));
    }
    accept(&first);
    // 任意非负整数资源价格都给出此界，价格来自 LP 的近似程度只影响界的强弱。
    let priced = model.fixed_value
        + (model
            .capacity
            .iter()
            .zip(&prices)
            .map(|(&b, &y)| b as i64 * y)
            .sum::<i64>()
            + model.items.iter().map(|it| it.bias.max(0)).sum::<i64>())
            / SCALE;
    let limited = |stats| {
        Ok(Some((
            SolveStatus::LimitReached,
            Some(priced as f64 + p.obj_offset),
            stats,
        )))
    };
    let mut tables: Vec<_> = model
        .capacity
        .iter()
        .map(|b| RowTable {
            stride: b + 1,
            values: vec![0; (model.items.len() + 1) * (b + 1)],
        })
        .collect();
    for _ in 0..12 {
        if !fill_tables(&mut tables, &model.items, deadline) {
            return limited(stats);
        }
        if !balance(&mut model.items, &tables, deadline) {
            if expired(deadline) {
                return limited(stats);
            }
            // 恢复每件物品的完整分摊，舍去这一次超出整数保护范围的调整。
            model.initial_shares(&prices);
            break;
        }
        model.items.reverse();
    }
    model
        .items
        .sort_unstable_by_key(|it| (Reverse(it.bias.abs()), it.column));
    if !fill_tables(&mut tables, &model.items, deadline) {
        return limited(stats);
    }
    let Some(tail) = Tail::new(&model.items, &model.capacity, deadline) else {
        return limited(stats);
    };
    let bad = std::cell::Cell::new(false);
    let mut checked = |x: &[f64]| {
        if p.constraints.iter().any(|r| r.expr.eval(x) > r.rhs) {
            bad.set(true);
        } else {
            accept(x);
        }
    };
    let (next, base) = model.greedy();
    if next > best {
        best = next;
        checked(&base);
    }
    let mut search = Search {
        items: &model.items,
        tables: &tables,
        tail,
        deadline,
        limit,
        nodes: 0,
        stopped: false,
        rem: model.capacity,
        values: model.fixed,
        best,
        accept: &mut checked,
    };
    let root_bound = priced.min(search.bound(0, model.fixed_value));
    search.visit(0, model.fixed_value);
    if bad.get() {
        return Err(MilpError::Numerical("多容量背包解未通过原约束核验".into()));
    }
    stats.nodes = search.nodes;
    let completed = !search.stopped && !expired(deadline);
    let status = if completed {
        SolveStatus::Optimal
    } else {
        SolveStatus::LimitReached
    };
    let bound = if completed { search.best } else { root_bound };
    Ok(Some((status, Some(bound as f64 + p.obj_offset), stats)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shares_bounds_and_tail_match_joint_capacity_table() {
        // 用两行联合容量的完整表作独立真值，覆盖免费物品、零收益、超重和相同收益。
        for seed in 1..=6 {
            let mut p = Problem::new("两种容量", Sense::Maximize);
            let mut profit = Vec::new();
            let mut a = Vec::new();
            let mut b = Vec::new();
            let mut joint = [[0; 12]; 14];
            for j in 0..68 {
                p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
                let w = (j * 7 + seed * (j % 3)) % 17;
                let h = (j * 5 + seed * (j % 7)) % 13;
                let c = ((j * seed + 3) % 11) as i64;
                profit.push((j, c as f64));
                a.push((j, w as f64));
                b.push((j, h as f64));
                for high in (w..14).rev() {
                    for low in (h..12).rev() {
                        joint[high][low] = joint[high][low].max(joint[high - w][low - h] + c);
                    }
                }
            }
            p.set_objective(profit, -11.).unwrap();
            p.add_constraint(None, a, RelOp::Le, 13.).unwrap();
            p.add_constraint(None, b, RelOp::Le, 11.).unwrap();
            let optimum = joint[13][11] as f64 - 11.;
            for limit in [1, 1_000_000] {
                let mut last = None;
                let (status, bound, stats) = solve(&p, Default::default(), None, limit, &mut |x| {
                    assert!(x.iter().all(|&v| v == 0. || v == 1.));
                    assert!(p.constraints.iter().all(|r| r.expr.eval(x) <= r.rhs));
                    let value = p.obj.eval(x) + p.obj_offset;
                    assert!(last.is_none_or(|old| value > old));
                    last = Some(value);
                })
                .unwrap()
                .unwrap();
                assert!(stats.nodes <= limit && last.is_some() && bound.unwrap() >= optimum);
                if status == SolveStatus::Optimal {
                    assert_eq!(bound, Some(optimum));
                    assert_eq!(last, bound);
                } else {
                    assert_eq!(limit, 1);
                    assert_eq!(status, SolveStatus::LimitReached);
                }
            }
            let (status, _, _) =
                solve(&p, Default::default(), Some(Instant::now()), 1, &mut |_| {
                    panic!("到时后不应接纳解")
                })
                .unwrap()
                .unwrap();
            assert_eq!(status, SolveStatus::LimitReached);
            p.constraints[0].rhs += 0.25;
            assert!(solve(&p, Default::default(), None, 1, &mut |_| panic!(
                "不应处理非整数容量"
            ))
            .unwrap()
            .is_none());
        }
    }
}
