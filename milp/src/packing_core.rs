//! 稠密多容量背包。资源价格给每列一个整数化的净收益：偏离其符号一定消耗相应的上界。
//! 先在净收益接近零的列上找可行解，再逐步降低待证门槛。证界时保留所有可能达到门槛
//! 的列，对尾部子集按原收益分桶，其余列用代理容量的后缀背包表剪枝。

use super::{
    expired, BbStats, Bounds, Instant, LpStatus, MilpError, Problem, RelOp, Result, Reverse, Sense,
    Simplex, SolveStatus, Tolerances, VarType,
};

const WIDTH: usize = 20;
const SCALE: i64 = 10_000;
const CELLS: usize = 12_000_000;
type Weights = [u16; WIDTH];

#[derive(Clone)]
struct Item {
    col: usize,
    profit: usize,
    weights: Weights,
    net: i64,
}

struct Model {
    items: Vec<Item>,
    caps: Vec<i32>,
}

impl Model {
    fn detect(p: &Problem, tol: Tolerances) -> Option<Self> {
        let m = p.constraints.len();
        let n = p.vars.len();
        if !(9..=WIDTH).contains(&m)
            || !(192..=768).contains(&n)
            || p.sense != Sense::Maximize
            || p.obj_offset.fract() != 0.0
            || p.obj_offset.abs() > (1u64 << 40) as f64
            || p.vars
                .iter()
                .any(|v| v.var_type != VarType::Binary || v.lb != 0.0 || v.ub != 1.0)
        {
            return None;
        }
        let mut items: Vec<_> = (0..n)
            .map(|col| Item {
                col,
                profit: 0,
                weights: [0; WIDTH],
                net: 0,
            })
            .collect();
        for &(j, c) in &p.obj.terms {
            if c.fract() != 0.0 || !(1.0..=127.0).contains(&c) {
                return None;
            }
            items[j].profit = c as usize;
        }
        if items.iter().any(|it| it.profit == 0) {
            return None;
        }
        let mut caps = Vec::new();
        for (k, row) in p.constraints.iter().enumerate() {
            if row.op != RelOp::Le
                || row.rhs.fract() != 0.0
                || !(1.0..=65535.0).contains(&row.rhs)
                || row.expr.terms.len() != n
            {
                return None;
            }
            let mut norm = 0.;
            for &(j, a) in &row.expr.terms {
                if a.fract() != 0.0 || !(1.0..=127.0).contains(&a) {
                    return None;
                }
                norm += a;
                items[j].weights[k] = a as u16;
            }
            if !(10.0 * tol.feasibility + tol.integrality * norm < 1.0
                && crate::solver::FINAL_VIOLATION_REL * (1.0 + row.rhs + norm) < 1.0)
            {
                return None;
            }
            caps.push(row.rhs as i32);
        }
        items.retain(|it| {
            caps.iter()
                .zip(&it.weights)
                .all(|(&b, &w)| i32::from(w) <= b)
        });
        Some(Self { items, caps })
    }

    // 净收益相反的决定至少损失其绝对值。达到门槛的解只可能改变这批列。
    fn improving_columns(&self, upper_scaled: i64, floor: usize) -> Vec<bool> {
        let gap = upper_scaled - floor as i64 * SCALE;
        self.items.iter().map(|it| it.net.abs() <= gap).collect()
    }
}

#[derive(Clone, Copy)]
struct HalfSubset {
    first: u16,
    second: u16,
    mask: u32,
}

struct Tail {
    split: usize,
    first: usize,
    second: usize,
    left: Vec<Weights>,
    right: Vec<Weights>,
    left_by_profit: Vec<Vec<HalfSubset>>,
    right_by_profit: Vec<Vec<HalfSubset>>,
    bins: Vec<Option<Vec<(u32, u32)>>>,
}
impl Tail {
    fn sums(items: &[&Item], m: usize) -> (Vec<Weights>, Vec<usize>) {
        let mut w = vec![[0; WIDTH]; 1 << items.len()];
        let mut p = vec![0; w.len()];
        for z in 1..w.len() {
            let i = &items[z.trailing_zeros() as usize];
            let prev = z & (z - 1);
            p[z] = p[prev] + i.profit;
            for k in 0..m {
                w[z][k] = w[prev][k] + i.weights[k];
            }
        }
        (w, p)
    }
    fn new(items: &[&Item], price: &[i64], size: usize, deadline: Option<Instant>) -> Option<Self> {
        let start = items.len().saturating_sub(size);
        let tail = &items[start..];
        let split = tail.len() / 2;
        let (left, a) = Self::sums(&tail[..split], price.len());
        let (right, b) = Self::sums(&tail[split..], price.len());
        let mut dims: Vec<_> = (0..price.len()).collect();
        dims.sort_unstable_by_key(|&k| (Reverse(price[k]), k));
        let first = dims[0];
        let second = dims[1];
        let bins = vec![None; a.last()? + b.last()? + 1];
        let by_profit = |weights: &[Weights], profits: &[usize], offset: usize| {
            let mut groups = vec![Vec::new(); profits.last().copied().unwrap_or(0) + 1];
            for (mask, (&gain, w)) in profits.iter().zip(weights).enumerate() {
                groups[gain].push(HalfSubset {
                    first: w[first],
                    second: w[second],
                    mask: (mask << offset) as u32,
                });
            }
            groups
        };
        let left_by_profit = by_profit(&left, &a, 0);
        let right_by_profit = by_profit(&right, &b, split);
        if expired(deadline) {
            return None;
        }
        Some(Self {
            split,
            first,
            second,
            left,
            right,
            left_by_profit,
            right_by_profit,
            bins,
        })
    }
    fn find(&mut self, value: usize, rem: &[i32], price: &[i64], slack: i64) -> Option<u32> {
        // 最多 21 列、单项最多 127，任何真实尾部用量均小于 4096。可行时剩余
        // 容量的加权和不超过 slack，故各正价资源的剩余量分别不超过 slack/price[k]。
        let bounds = |k: usize| {
            let hi = rem[k] as u32;
            let lo = if price[k] > 0 {
                (i64::from(hi) - slack / price[k]).max(0) as u32
            } else {
                0
            };
            (lo, hi.min(4095))
        };
        let (a, b) = bounds(self.first);
        let (c, d) = bounds(self.second);
        if c > d || a / 8 > b / 8 {
            return None;
        }
        // 分支实际只会询问少数尾部收益。第一次查询该收益时才连接左右子集，
        // 同一收益的后续查询复用按资源占用排序的桶。每对左右子集只属于一个桶。
        let bin = self.bins[value].get_or_insert_with(|| {
            let first = value.saturating_sub(self.right_by_profit.len() - 1);
            let last = value.min(self.left_by_profit.len() - 1);
            let groups = (first..=last).map(|gain| {
                (
                    &self.left_by_profit[gain],
                    &self.right_by_profit[value - gain],
                )
            });
            let capacity = groups.clone().map(|(l, r)| l.len() * r.len()).sum();
            let mut pairs = Vec::with_capacity(capacity);
            for (left, right) in groups {
                for l in left {
                    for r in right {
                        // 第一种资源按八单位分箱，箱内按第二种资源和掩码排序。
                        let primary = u32::from(l.first + r.first);
                        let secondary = u32::from(l.second + r.second);
                        pairs.push(((primary / 8) * 4096 + secondary, l.mask | r.mask));
                    }
                }
            }
            pairs.sort_unstable();
            pairs
        });
        if bin.is_empty() {
            return None;
        }
        for block in a / 8..=b / 8 {
            let key = block * 4096;
            let pos = bin.partition_point(|&(k, _)| k < key + c);
            for &(k, mask) in &bin[pos..] {
                if k > key + d {
                    break;
                }
                let l = &self.left[(mask as usize) & (self.left.len() - 1)];
                let r = &self.right[mask as usize >> self.split];
                if rem
                    .iter()
                    .enumerate()
                    .all(|(j, &cap)| i32::from(l[j] + r[j]) <= cap)
                {
                    return Some(mask);
                }
            }
        }
        None
    }
}

struct Table {
    prices: Vec<u32>,
    layers: Vec<Vec<u32>>,
}

fn price_weight(weights: &Weights, prices: &[u32]) -> u32 {
    prices
        .iter()
        .zip(weights)
        .map(|(&y, &w)| y * u32::from(w))
        .sum()
}

struct Window<'a> {
    items: Vec<&'a Item>,
    base: usize,
    rem: Vec<i32>,
    x: Vec<f64>,
    tail: usize,
    leaf: Tail,
    first_table: usize,
    middle: Vec<(Weights, usize, i64)>,
    inner: usize,
    tables: Vec<Table>,
    price: &'a [i64],
    bonus: Vec<i64>,
}

impl<'a> Window<'a> {
    fn new(
        model: &'a Model,
        price: &'a [i64],
        selected: &[bool],
        outside: &[f64],
        deadline: Option<Instant>,
        tail_size: usize,
        join: usize,
    ) -> Option<Self> {
        let m = model.caps.len();
        let mut items = Vec::new();
        let mut x = outside.to_vec();
        let mut rem = model.caps.clone();
        let mut base = 0;
        for (it, &inside) in model.items.iter().zip(selected) {
            if inside {
                items.push(it);
                x[it.col] = 0.;
            } else if x[it.col] == 1. {
                base += it.profit;
                for (cap, &w) in rem.iter_mut().zip(&it.weights) {
                    *cap -= i32::from(w);
                }
            }
        }
        items.sort_unstable_by_key(|it| (Reverse(it.net.abs()), it.col));
        let inner = items.len().saturating_sub(tail_size);
        let tail = inner.saturating_sub(join);
        let mut sums = vec![0; items.len() + 1];
        let mut bonus = vec![0; items.len() + 1];
        for i in (0..items.len()).rev() {
            sums[i] = sums[i + 1] + items[i].profit;
            bonus[i] = bonus[i + 1] + items[i].net.max(0);
        }
        let mut first_table = tail;
        let mut cells = sums[tail] + 1;
        while first_table > 0 && cells + sums[first_table - 1] + 1 <= CELLS {
            first_table -= 1;
            cells += sums[first_table] + 1;
        }
        if expired(deadline) {
            return None;
        }
        let leaf = Tail::new(&items, price, tail_size, deadline)?;
        let (weights, p) = Tail::sums(&items[tail..inner], m);
        let middle = weights
            .into_iter()
            .zip(p)
            .map(|(w, v)| {
                let cost = price.iter().zip(&w).map(|(&y, &a)| y * i64::from(a)).sum();
                (w, v, cost)
            })
            .collect();
        let mut tables = Vec::new();
        for h in 0..(1 + 2 * m).min(CELLS / cells) {
            if expired(deadline) {
                return None;
            }
            let mut prices: Vec<_> = price.iter().map(|&v| v as u32).collect();
            if h > 0 {
                let k = (h - 1) % m;
                if h <= m {
                    prices[k] += 900;
                } else if prices[k] > 0 {
                    prices[k] = prices[k].saturating_sub(900);
                } else {
                    continue;
                }
            }
            // 非负价格加权仍是原容量的不等式。表中存达到至少给定收益所需的最小加权资源。
            // u32 饱和只会把最小资源需求往下压，仍是有效下界；每件物品自身的
            // 加权费用最多 20*127*(127*SCALE+900)，不会在 price_weight 中溢出。
            let mut layers = vec![Vec::new(); tail + 1];
            let mut last = vec![u32::MAX; sums[tail] + 1];
            last[0] = 0;
            let mut len = 0;
            for it in &items[tail..] {
                let cost = price_weight(&it.weights, &prices);
                len += it.profit;
                for v in (it.profit..=len).rev() {
                    last[v] = last[v].min(last[v - it.profit].saturating_add(cost));
                }
            }
            for v in (1..last.len()).rev() {
                last[v - 1] = last[v - 1].min(last[v]);
            }
            layers[tail] = last;
            for i in (first_table..tail).rev() {
                if i % 16 == 0 && expired(deadline) {
                    return None;
                }
                let it = items[i];
                let cost = price_weight(&it.weights, &prices);
                let next = &layers[i + 1];
                let mut row = vec![u32::MAX; sums[i] + 1];
                row[..next.len()].copy_from_slice(next);
                row[..it.profit]
                    .iter_mut()
                    .for_each(|x| *x = (*x).min(cost));
                for (v, &a) in row[it.profit..].iter_mut().zip(next) {
                    *v = (*v).min(a.saturating_add(cost));
                }
                layers[i] = row;
            }
            tables.push(Table { prices, layers });
        }
        Some(Self {
            items,
            base,
            rem,
            x,
            tail,
            leaf,
            first_table,
            tables,
            price,
            bonus,
            middle,
            inner,
        })
    }
}

struct Progress<'a> {
    best: usize,
    x: Vec<f64>,
    nodes: u64,
    stop: u64,
    deadline: Option<Instant>,
    accept: &'a mut dyn FnMut(&[f64]),
}

impl Window<'_> {
    fn search(&mut self, i: usize, profit: usize, floor: usize, run: &mut Progress<'_>) -> bool {
        if run.nodes >= run.stop || (run.nodes % 128 == 0 && expired(run.deadline)) {
            return false;
        }
        run.nodes += 1;
        // 净收益的绝对值大于剩余预算时，走反方向连目标上界也够不到。
        // 这些列已按绝对值降序排列，一次把整个被迫前缀固定在当前节点。
        let cap: i64 = self
            .price
            .iter()
            .zip(&self.rem)
            .map(|(&y, &b)| y * i64::from(b))
            .sum();
        let gap =
            profit as i64 * SCALE + cap + self.bonus[i] - (run.best + 1).max(floor) as i64 * SCALE;
        if gap < 0 {
            return true;
        }
        let end = i + self.items[i..self.tail].partition_point(|it| it.net.abs() > gap);
        let mut fixed = profit;
        for &it in &self.items[i..end] {
            self.x[it.col] = f64::from(it.net > 0);
            if it.net > 0 {
                fixed += it.profit;
                for (b, &w) in self.rem.iter_mut().zip(&it.weights) {
                    *b -= i32::from(w);
                }
            }
        }
        let done = self.branches(end, fixed, floor, run);
        for &it in &self.items[i..end] {
            if it.net > 0 {
                for (b, &w) in self.rem.iter_mut().zip(&it.weights) {
                    *b += i32::from(w);
                }
            }
        }
        done
    }
    fn branches(&mut self, i: usize, profit: usize, floor: usize, run: &mut Progress<'_>) -> bool {
        let threshold = (run.best + 1).max(floor);
        let need = threshold.saturating_sub(profit);
        if self.rem.iter().any(|&b| b < 0) {
            return true;
        }
        let cap: i64 = self
            .price
            .iter()
            .zip(&self.rem)
            .map(|(&y, &b)| y * i64::from(b))
            .sum();
        let priced = cap + self.bonus[i];
        if priced < need as i64 * SCALE {
            return true;
        }
        if i >= self.first_table {
            if need >= self.tables[0].layers[i].len() {
                return true;
            }
            for t in &self.tables {
                let cap: i64 = t
                    .prices
                    .iter()
                    .zip(&self.rem)
                    .map(|(&y, &b)| i64::from(y) * i64::from(b))
                    .sum();
                if i64::from(t.layers[i][need]) > cap {
                    return true;
                }
            }
        }
        if i == self.tail {
            // 尾部左右两张子集表按所需收益连接。左表的每种取法只需在右表中查对应收益桶。
            for (left_mask, (weights, left_profit, cost)) in self.middle.iter().enumerate() {
                let mut residual = [0; WIDTH];
                for (k, &b) in self.rem.iter().enumerate() {
                    residual[k] = b - i32::from(weights[k]);
                }
                let rem = &residual[..self.rem.len()];
                if rem.iter().any(|&b| b < 0) {
                    continue;
                }
                let leaf_cap = cap - cost + self.bonus[self.inner];
                let used = profit + left_profit;
                let target = (run.best + 1).max(floor).saturating_sub(used);
                let highest = (leaf_cap / SCALE) as usize;
                for value in (target..=highest.min(self.leaf.bins.len() - 1)).rev() {
                    if expired(run.deadline) {
                        return false;
                    }
                    if let Some(mask) =
                        self.leaf
                            .find(value, rem, self.price, leaf_cap - value as i64 * SCALE)
                    {
                        for (k, it) in self.items[self.tail..self.inner].iter().enumerate() {
                            self.x[it.col] = ((left_mask >> k) & 1) as f64;
                        }
                        for (k, it) in self.items[self.inner..].iter().enumerate() {
                            self.x[it.col] = f64::from((mask >> k) & 1);
                        }
                        run.best = used + value;
                        run.x.clone_from(&self.x);
                        (run.accept)(&run.x);
                        break;
                    }
                }
            }
            return true;
        }
        let it = self.items[i];
        for take in [it.net >= 0, it.net < 0] {
            self.x[it.col] = f64::from(take);
            if take {
                for (b, &w) in self.rem.iter_mut().zip(&it.weights) {
                    *b -= i32::from(w);
                }
            }
            let done = self.search(i + 1, profit + usize::from(take) * it.profit, floor, run);
            if take {
                for (b, &w) in self.rem.iter_mut().zip(&it.weights) {
                    *b += i32::from(w);
                }
            }
            if !done {
                return false;
            }
        }
        true
    }
}

fn neighborhood(items: &[Item], best: &[f64], rep: usize, rand: &mut u32) -> Vec<bool> {
    let n = items.len();
    let mut rank = Vec::with_capacity(n);
    for (j, it) in items.iter().enumerate() {
        *rand = rand.wrapping_mul(1664525).wrapping_add(1013904223);
        let key = if rep == 0 {
            j as u32
        } else if j < 14 {
            0
        } else if j >= 50 + (rep % 4) * 8 {
            u32::MAX
        } else {
            (*rand >> 13) * (5 + j as u32)
                / if (it.net > 0) != (best[it.col] == 1.) {
                    8
                } else {
                    1
                }
        };
        rank.push((key, j));
    }
    rank.sort_unstable();
    let mut selected = vec![false; n];
    for &(_, j) in &rank[..36.min(n)] {
        selected[j] = true;
    }
    selected
}

pub(super) fn solve(
    p: &Problem,
    tol: Tolerances,
    deadline: Option<Instant>,
    limit: u64,
    accept: &mut dyn FnMut(&[f64]),
) -> Result<Option<(SolveStatus, Option<f64>, BbStats)>> {
    let Some(mut model) = Model::detect(p, tol) else {
        return Ok(None);
    };
    let mut stats = BbStats::default();
    if expired(deadline) {
        return Ok(Some((SolveStatus::LimitReached, None, stats)));
    }
    let mut lp = Simplex::new(p, tol)?;
    lp.set_time_deadline(deadline);
    let sol = lp.solve(&Bounds::from_problem(p));
    stats.lp_iterations = lp.iterations();
    if matches!(sol.status, LpStatus::TimeLimit | LpStatus::IterationLimit) || expired(deadline) {
        return Ok(Some((SolveStatus::LimitReached, None, stats)));
    }
    if sol.status != LpStatus::Optimal {
        return Ok(None);
    }
    // 价格只作选列与剪枝。取整后重新计算整数上界，不依赖浮点 LP 的最优性。
    let prices: Vec<i64> = sol
        .duals
        .iter()
        .map(|&v| (-v * SCALE as f64).round().clamp(0., 127. * SCALE as f64) as i64)
        .collect();
    let mut ub: i64 = prices
        .iter()
        .zip(&model.caps)
        .map(|(&y, &b)| y * i64::from(b))
        .sum();
    let mut signs = vec![0.; p.vars.len()];
    for it in &mut model.items {
        it.net = it.profit as i64 * SCALE
            - prices
                .iter()
                .zip(&it.weights)
                .map(|(&y, &w)| y * i64::from(w))
                .sum::<i64>();
        ub += it.net.max(0);
        signs[it.col] = f64::from(it.net > 0);
    }
    let mut order: Vec<_> = model.items.iter().collect();
    order.sort_unstable_by_key(|it| (Reverse(it.net), it.col));
    let mut rem = model.caps.clone();
    let mut best = 0;
    let mut x = vec![0.; p.vars.len()];
    for it in order {
        if rem
            .iter()
            .zip(&it.weights)
            .all(|(&b, &w)| b >= i32::from(w))
        {
            for (b, &w) in rem.iter_mut().zip(&it.weights) {
                *b -= i32::from(w);
            }
            best += it.profit;
            x[it.col] = 1.;
        }
    }
    let invalid = std::cell::Cell::new(false);
    let mut checked = |x: &[f64]| {
        if p.constraints.iter().any(|r| r.expr.eval(x) > r.rhs) {
            invalid.set(true);
        } else {
            accept(x);
        }
    };
    checked(&x);
    let mut run = Progress {
        best,
        x,
        nodes: 0,
        stop: limit,
        deadline,
        accept: &mut checked,
    };
    model
        .items
        .sort_unstable_by_key(|it| (it.net.abs(), it.col));
    let mut rand = 101u32;
    let mut upper = (ub / SCALE) as usize;
    // 固定邻域外的列只用来找解，未解完的邻域不会产生全局界。
    let quota = limit.min(700_000) / 2;
    let mut first_complete = false;
    for rep in 0..36 {
        if expired(deadline) || run.nodes >= quota || run.best >= upper {
            break;
        }
        // 首个局部窗口若在配额内搜完，用下一份相同配额搜索全部可能改进的列。
        // 该窗口完整结束就已排除了更好解；中断只保留已复核的可行解。
        let try_proof = rep == 1 && first_complete;
        let floor = if try_proof { run.best + 1 } else { 0 };
        let selected = if try_proof {
            model.improving_columns(ub, floor)
        } else {
            neighborhood(&model.items, &run.x, rep, &mut rand)
        };
        let outside = if try_proof || rep % 3 == 0 {
            &signs
        } else {
            &run.x
        };
        let Some(mut w) = Window::new(&model, &prices, &selected, outside, deadline, 20, 0) else {
            break;
        };
        run.stop = (run.nodes + 40_000).min(quota);
        let complete = w.search(0, w.base, floor, &mut run);
        if try_proof && complete {
            upper = run.best;
            break;
        }
        first_complete = rep == 0 && complete;
    }
    run.stop = limit;
    // 若损失超过上界与目标之间的差，则该列在达到门槛的任何解中都必须遵从净收益符号。
    // 完整排除一段收益后才降低全局上界。时间或节点用尽时保留最后已证明的上界。
    while run.best < upper && run.nodes < limit && !expired(deadline) {
        let floor = upper.saturating_sub(1).max(run.best + 1);
        let selected = model.improving_columns(ub, floor);
        let Some(mut w) = Window::new(&model, &prices, &selected, &signs, deadline, 21, 8) else {
            break;
        };
        if !w.search(0, w.base, floor, &mut run) {
            break;
        }
        upper = run.best.max(floor - 1);
    }
    if invalid.get() {
        return Err(MilpError::Numerical("多容量装载解未通过原行复核".into()));
    }
    stats.nodes = run.nodes;
    let status = if run.best == upper && !expired(deadline) {
        SolveStatus::Optimal
    } else {
        SolveStatus::LimitReached
    };
    Ok(Some((status, Some(p.obj_offset + upper as f64), stats)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tail_builds_only_queried_profits_and_matches_exhaustion() {
        let items: Vec<_> = (0..8)
            .map(|j| {
                let mut weights = [0; WIDTH];
                weights[..3].copy_from_slice(&[
                    (j % 3 + 1) as u16,
                    (2 * j % 5 + 1) as u16,
                    (j % 4 + 2) as u16,
                ]);
                Item {
                    col: j,
                    profit: (j % 3 + 1) * 2,
                    weights,
                    net: 0,
                }
            })
            .collect();
        let refs: Vec<_> = items.iter().collect();
        let prices = [6, 3, 0];
        let mut tail = Tail::new(&refs, &prices, refs.len(), None).unwrap();
        assert!(tail.bins.iter().all(Option::is_none));
        // 相同收益在不同的剩余容量下应复用同一个完整索引；奇数收益的桶为空。
        for profit in [12, 11, 16] {
            for cap in [[3, 7, 8], [11, 10, 14], [8, 19, 24]] {
                let expected = (0u32..1 << items.len())
                    .filter_map(|mask| {
                        let mut weights = [0; 3];
                        let mut gain = 0;
                        for item in &items {
                            if mask & (1 << item.col) != 0 {
                                gain += item.profit;
                                for (w, &a) in weights.iter_mut().zip(&item.weights) {
                                    *w += u32::from(a);
                                }
                            }
                        }
                        (gain == profit && weights.iter().zip(cap).all(|(&w, b)| w <= b as u32))
                            .then_some((weights[0] / 8, weights[1], mask))
                    })
                    .min()
                    .map(|(_, _, mask)| mask);
                assert_eq!(tail.find(profit, &cap, &prices, i64::MAX), expected);
            }
        }
        assert_eq!(tail.bins.iter().filter(|b| b.is_some()).count(), 3);
    }

    #[test]
    fn dense_core_matches_exhaustion_and_keeps_valid_interrupted_bounds() {
        for seed in 1..=4 {
            let n = 384;
            let mut p = Problem::new("多资源", Sense::Maximize);
            let gains: Vec<_> = (0..n).map(|j| 1 + (j * 37 + seed * 17) % 70).collect();
            for j in 0..n {
                p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
            }
            p.set_objective(
                gains
                    .iter()
                    .enumerate()
                    .map(|(j, &v)| (j, v as f64))
                    .collect(),
                -11.,
            )
            .unwrap();
            for k in 0..9 {
                let terms = (0..n)
                    .map(|j| {
                        (
                            j,
                            if j < 16 {
                                (1 + (j * (k + 3) + seed * (j + 2)) % 9) as f64
                            } else {
                                127.
                            },
                        )
                    })
                    .collect();
                p.add_constraint(None, terms, RelOp::Le, (32 + k % 3) as f64)
                    .unwrap();
            }
            let mut truth = 0usize;
            for mask in 0..(1u32 << 16) {
                let value = (0..16)
                    .filter(|&j| mask >> j & 1 == 1)
                    .map(|j| gains[j])
                    .sum();
                if value > truth
                    && p.constraints.iter().all(|r| {
                        r.expr.terms[..16]
                            .iter()
                            .filter(|&&(j, _)| mask >> j & 1 == 1)
                            .map(|&(_, a)| a)
                            .sum::<f64>()
                            <= r.rhs
                    })
                {
                    truth = value;
                }
            }
            let optimum = truth as f64 - 11.;
            for limit in [1, 100_000] {
                let mut best = f64::NEG_INFINITY;
                let (status, bound, stats) = solve(&p, Default::default(), None, limit, &mut |x| {
                    assert!(x.iter().all(|&v| v == 0. || v == 1.));
                    assert!(p.constraints.iter().all(|r| r.expr.eval(x) <= r.rhs));
                    let z = p.obj.eval(x) + p.obj_offset;
                    assert!(z > best && z <= optimum);
                    best = z;
                })
                .unwrap()
                .unwrap();
                assert!(stats.nodes <= limit && bound.unwrap() >= optimum);
                if status == SolveStatus::Optimal {
                    assert_eq!(bound, Some(optimum));
                    assert_eq!(best, optimum);
                } else {
                    assert_eq!(limit, 1);
                    assert_eq!(status, SolveStatus::LimitReached);
                }
            }
            assert_eq!(
                solve(
                    &p,
                    Default::default(),
                    Some(Instant::now()),
                    5,
                    &mut |_| panic!("逾期不应回调")
                )
                .unwrap()
                .unwrap()
                .0,
                SolveStatus::LimitReached
            );
        }
    }
}
