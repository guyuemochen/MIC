//! 消去每行私有的补差列后，共享连续列的活动度必须落在若干等距区间之一。
//! 该投影的根松弛供应全局界，固定部分共享列的邻域只供应可行点。

use super::{expired, small_integer, Repair};
use crate::branch_bound::BbStats;
use crate::model::{LinearExpr, Problem, RelOp, Sense, VarType};
use crate::simplex::{Bounds, LpSolution, LpStatus};
use crate::solver::SolveOptions;
use std::collections::BTreeSet;
use std::time::Instant;

#[derive(Clone)]
struct Band {
    terms: Vec<(usize, i32)>,
    low: i32,
    high: i32,
    period: i32,
    count: i32,
}

impl Band {
    // 用盒界剔除够不到的两端区间，再取剩余区间的凸包。私有周数上下界
    // 已折算为从第零个到第 count 个区间，包括周期系数为负的情形。
    fn hull(&self, bounds: &[(i32, i32)]) -> Option<(i32, i32, i32, i32)> {
        let mut lo = 0;
        let mut hi = 0;
        for &(j, a) in &self.terms {
            let (l, u) = bounds[j];
            lo += (a * l).min(a * u);
            hi += (a * l).max(a * u);
        }
        let (first, last) = if self.period == 0 {
            (0, 0)
        } else {
            (
                -((self.high - lo).div_euclid(self.period)).min(0),
                ((hi - self.low).div_euclid(self.period)).min(self.count),
            )
        };
        (first <= last && hi >= self.low && lo <= self.high + self.period * self.count).then_some((
            lo,
            hi,
            self.low + first * self.period,
            self.high + last * self.period,
        ))
    }

    fn split(&self, x: &[f64]) -> Option<(f64, i32, i32)> {
        let gap = self.period - (self.high - self.low);
        if self.period == 0 || gap <= 0 {
            return None;
        }
        let sum = self
            .terms
            .iter()
            .map(|&(j, a)| x[j] * f64::from(a))
            .sum::<f64>();
        let k = ((sum - f64::from(self.low)) / f64::from(self.period)).floor() as i32;
        let left = self.high + k * self.period;
        let right = self.low + (k + 1) * self.period;
        let dist = (sum - f64::from(left)).min(f64::from(right) - sum);
        (dist > 1e-7).then_some((dist / f64::from(gap), left, right))
    }

    // 角点处活动度在周期缺口时，以每个变量到该角点的非负距离重写活动度
    // 为 正距离和−负距离和。要越过缺口，正距离和至少为 r 或负距离和至少为 s，
    // 故 s·正距离和+r·负距离和 ≥ rs。任一端点选择都产生全局有效不等式。
    fn cut(&self, x: &[f64], bounds: &[(i32, i32)]) -> Option<(Vec<(usize, i32)>, i64)> {
        let gap = self.period - (self.high - self.low);
        if self.period == 0 || gap <= 0 {
            return None;
        }
        let mut base = 0;
        let mut pos = 0.;
        let mut neg = 0.;
        let mut ends = Vec::with_capacity(self.terms.len());
        for &(j, a) in &self.terms {
            let (l, u) = bounds[j];
            let reversed = x[j] - f64::from(l) > f64::from(u) - x[j];
            let end = if reversed { u } else { l };
            let sign = if reversed { -1 } else { 1 };
            let distance = (x[j] - f64::from(end)) * f64::from(sign);
            if a * sign > 0 {
                pos += distance;
            } else {
                neg += distance;
            }
            base += a * end;
            ends.push((sign, end, distance));
        }
        let violation = |b: i32, p: f64, n: f64| {
            let r = (self.low - b).rem_euclid(self.period);
            if r == 0 || r >= gap {
                return 0.;
            }
            f64::from(r * (gap - r)) - f64::from(gap - r) * p - f64::from(r) * n
        };
        let mut best = violation(base, pos, neg);
        let mut flip = None;
        for (i, &(j, a)) in self.terms.iter().enumerate() {
            let (sign, _, distance) = ends[i];
            let w = bounds[j].1 - bounds[j].0;
            let opposite = f64::from(w) - distance;
            let (dp, dn) = if a * sign > 0 {
                (-distance, opposite)
            } else {
                (opposite, -distance)
            };
            let v = violation(base + a * sign * w, pos + dp, neg + dn);
            if v > best {
                best = v;
                flip = Some(i);
            }
        }
        if best <= 1e-6 {
            return None;
        }
        if let Some(i) = flip {
            let (j, a) = self.terms[i];
            let (s, end, _) = &mut ends[i];
            let delta = *s * (bounds[j].1 - bounds[j].0);
            *s = -*s;
            *end += delta;
            base += a * delta;
        }
        let r = (self.low - base).rem_euclid(self.period);
        let mut rhs = i64::from(r * (gap - r));
        let mut cut = Vec::with_capacity(self.terms.len());
        for (&(j, a), &(s, end, _)) in self.terms.iter().zip(&ends) {
            let c = s * if a * s > 0 { gap - r } else { r };
            rhs += i64::from(c) * i64::from(end);
            cut.push((j, c));
        }
        Some((cut, rhs))
    }
}

struct Projection {
    columns: Vec<usize>,
    bounds: Vec<(i32, i32)>,
    costs: Vec<i64>,
    offset: i64,
    bands: Vec<Band>,
}

impl Projection {
    fn build(p: &Problem, r: &Repair, deadline: Option<Instant>) -> Option<Self> {
        let sign = if p.sense == Sense::Minimize { 1 } else { -1 };
        let mut coeff = vec![0i64; p.vars.len()];
        let mut constant = i64::from(small_integer(p.obj_offset)?) * sign;
        for &(j, a) in &p.obj.terms {
            coeff[j] = i64::from(small_integer(a)?) * sign;
        }
        for row in &r.rows {
            if expired(deadline) {
                return None;
            }
            if let Some((j, c)) = row.local {
                let price = coeff[j] * i64::from(c);
                coeff[j] = 0;
                constant = constant.checked_add(price * i64::from(row.rhs))?;
                for &(j, a) in &row.terms {
                    coeff[j] -= price * i64::from(a);
                }
                if let Some(q) = row.cycles {
                    coeff[q.column] -= price * i64::from(q.coefficient);
                }
            }
            if row.cycles.is_some_and(|q| coeff[q.column] != 0) {
                return None;
            }
        }
        let mut map = vec![0; p.vars.len()];
        let mut columns = Vec::new();
        let mut bounds = Vec::new();
        let mut costs = Vec::new();
        for (j, &price) in coeff.iter().enumerate() {
            let (l, u) = r.bounds[j];
            if r.incidence[j].is_empty() {
                constant =
                    constant.checked_add(price * i64::from(if price < 0 { u } else { l }))?;
            } else {
                map[j] = columns.len();
                columns.push(j);
                bounds.push((l, u));
                costs.push(price);
            }
        }
        // 该常数及投影行的整数系数在转给浮点 LP 时必须能精确表示。
        if columns.is_empty() || constant.unsigned_abs() > 1 << 50 {
            return None;
        }
        let mut bands = Vec::new();
        for &i in &r.constraints {
            let row = &r.rows[i];
            let (period, origin, count) = row.cycles.map_or((0, 0, 0), |q| {
                (
                    q.coefficient.abs(),
                    (q.coefficient * q.low).max(q.coefficient * q.high),
                    q.high - q.low,
                )
            });
            bands.push(Band {
                terms: row.terms.iter().map(|&(j, a)| (map[j], a)).collect(),
                low: row.lower - origin,
                high: row.upper - origin,
                period,
                count,
            });
        }
        Some(Self {
            columns,
            bounds,
            costs,
            offset: constant,
            bands,
        })
    }

    fn linear_model(&self) -> Option<Problem> {
        let mut p = Problem::new("周期区间", Sense::Minimize);
        for (j, &(l, u)) in self.bounds.iter().enumerate() {
            p.vars.push(crate::model::Variable {
                name: format!("x{j}"),
                var_type: VarType::Continuous,
                lb: l.into(),
                ub: u.into(),
            });
        }
        p.obj = LinearExpr {
            terms: self
                .costs
                .iter()
                .enumerate()
                .filter_map(|(j, &c)| (c != 0).then_some((j, c as f64)))
                .collect(),
        };
        p.obj_offset = self.offset as f64;
        for row in &self.bands {
            let (min, max, lo, hi) = row.hull(&self.bounds)?;
            if lo > min {
                add_row(&mut p, &row.terms, RelOp::Ge, lo.into());
            }
            if hi < max {
                add_row(&mut p, &row.terms, RelOp::Le, hi.into());
            }
        }
        Some(p)
    }
}

fn add_row(p: &mut Problem, terms: &[(usize, i32)], op: RelOp, rhs: i64) {
    let expr = LinearExpr {
        terms: terms.iter().map(|&(j, a)| (j, f64::from(a))).collect(),
    };
    p.constraints.push(crate::model::Constraint {
        name: None,
        expr,
        op,
        rhs: rhs as f64,
    });
}

// 量化后的对偶乘子不必最优。对剩余目标的每个系数独立取盒界下界，
// 就给出严格的弱对偶证书。输入仅来自上述整数投影，由 i128 完成全部求和。
fn lower_bound(p: &Problem, duals: &[f64]) -> Option<f64> {
    const SCALE: i128 = 1 << 24;
    if duals.len() != p.constraints.len() {
        return None;
    }
    let mut cost = vec![0i128; p.vars.len()];
    for &(j, a) in &p.obj.terms {
        cost[j] = a as i128 * SCALE;
    }
    let mut value = (p.obj_offset as i128) * SCALE;
    for (r, &d) in p.constraints.iter().zip(duals) {
        if !d.is_finite() || d.abs() > 1e9 {
            return None;
        }
        let mut v = (d * SCALE as f64).round() as i128;
        match r.op {
            RelOp::Ge => v = v.max(0),
            RelOp::Le => v = v.min(0),
            RelOp::Eq => {}
        }
        value += v * (r.rhs as i128);
        for &(j, a) in &r.expr.terms {
            cost[j] -= v * (a as i128);
        }
    }
    for (&a, var) in cost.iter().zip(&p.vars) {
        value += a * (if a < 0 { var.ub } else { var.lb } as i128);
    }
    let f = value as f64;
    Some((if f as i128 > value { f.next_down() } else { f }) / (SCALE as f64))
}

fn relaxation(
    p: &Problem,
    opts: &SolveOptions,
    deadline: Option<Instant>,
    stats: &mut BbStats,
) -> Option<LpSolution> {
    if expired(deadline) {
        return None;
    }
    let mut engine = crate::simplex_sparse::new_engine(p, opts.tol, opts.bb.lp_engine).ok()?;
    engine.set_time_deadline(deadline);
    let sol = engine.solve(&Bounds::from_problem(p));
    stats.lp_iterations += engine.iterations();
    (sol.status == LpStatus::Optimal
        && sol.x.len() == p.vars.len()
        && sol
            .x
            .iter()
            .zip(&p.vars)
            .all(|(&x, v)| x >= v.lb - 1e-6 && x <= v.ub + 1e-6))
    .then_some(sol)
}

// 从一个仍可改善目标的列出发，每次展开一条相邻周期行，得到共享列的局部邻域。
// 未选列保持当前值。每行独立重算活动度上下界，可直接裁去不可达区间。
fn restrict(g: &Projection, indices: &[usize], incumbent: &[i32]) -> Projection {
    let mut map = vec![None; incumbent.len()];
    for (k, &j) in indices.iter().enumerate() {
        map[j] = Some(k);
    }
    let mut bands = Vec::new();
    let bounds: Vec<_> = indices.iter().map(|&j| g.bounds[j]).collect();
    for b in &g.bands {
        let mut terms = Vec::new();
        let mut fixed = 0;
        for &(j, c) in &b.terms {
            if let Some(k) = map[j] {
                terms.push((k, c));
            } else {
                fixed += c * incumbent[j];
            }
        }
        if !terms.is_empty() {
            let row = Band {
                terms,
                low: b.low - fixed,
                high: b.high - fixed,
                period: b.period,
                count: b.count,
            };
            bands.push(row);
        }
    }
    let offset = g.offset
        + g.costs
            .iter()
            .enumerate()
            .filter(|(j, _)| map[*j].is_none())
            .map(|(j, &c)| c * i64::from(incumbent[j]))
            .sum::<i64>();
    Projection {
        columns: indices.to_vec(),
        bounds,
        costs: indices.iter().map(|&j| g.costs[j]).collect(),
        offset,
        bands,
    }
}

struct Neighborhood<'a> {
    original: &'a Problem,
    repair: &'a Repair,
    global: &'a Projection,
    opts: &'a SolveOptions,
    deadline: Option<Instant>,
    remaining: u64,
    stats: &'a mut BbStats,
    accept: &'a mut dyn FnMut(&[f64]),
    values: Vec<i32>,
    best: f64,
    global_bound: Option<f64>,
}

impl Neighborhood<'_> {
    fn search(&mut self, lp: &mut Problem, local: &Projection, window_nodes: &mut u64) -> bool {
        if self.remaining == 0 || *window_nodes == 0 || expired(self.deadline) {
            return false;
        }
        self.remaining -= 1;
        *window_nodes -= 1;
        self.stats.nodes += 1;
        let Some(sol) = relaxation(lp, self.opts, self.deadline, self.stats) else {
            return false;
        };
        if sol.objective >= self.best - 1e-7 {
            return false;
        }
        let mut split = None;
        for b in &local.bands {
            if let Some((fraction, left, right)) = b.split(&sol.x) {
                if split.is_none_or(|(old, _, _, _)| fraction > old) {
                    split = Some((fraction, b, left, right));
                }
            }
        }
        if let Some((_, b, left, right)) = split {
            let incumbent_sum: i32 = b
                .terms
                .iter()
                .map(|&(j, c)| c * self.values[self.global.columns[local.columns[j]]])
                .sum();
            let mut branches = [(RelOp::Le, left), (RelOp::Ge, right)];
            if incumbent_sum >= right {
                branches.swap(0, 1);
            }
            for (op, rhs) in branches {
                add_row(lp, &b.terms, op, rhs.into());
                let found = self.search(lp, local, window_nodes);
                lp.constraints.pop();
                if found {
                    return true;
                }
            }
            return false;
        }
        let mut point = self.values.clone();
        for (&j, &v) in local.columns.iter().zip(&sol.x) {
            if small_integer(v.round()).is_none() {
                return false;
            }
            point[self.global.columns[j]] = v.round() as i32;
        }
        let Some(x) = self.repair.restore(self.original, point, self.deadline) else {
            return false;
        };
        let raw = self.original.obj.eval(&x) + self.original.obj_offset;
        let cost = if self.original.sense == Sense::Minimize {
            raw
        } else {
            -raw
        };
        if cost >= self.best {
            return false;
        }
        self.best = cost;
        self.values = x.iter().map(|&v| v as i32).collect();
        (self.accept)(&x);
        true
    }

    fn run(&mut self) {
        let g = self.global;
        let mut adj = vec![Vec::new(); g.columns.len()];
        for (i, b) in g.bands.iter().enumerate() {
            if b.period > 0 {
                for &(j, _) in &b.terms {
                    adj[j].push(i);
                }
            }
        }
        let mut state = 0x6a09e667f3bcc909u64;
        let mut pick = |count: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % count as u64) as usize
        };
        for turn in 0..768 {
            if self.remaining == 0
                || expired(self.deadline)
                || self.global_bound.is_some_and(|b| self.best <= b)
            {
                break;
            }
            let incumbent: Vec<_> = g.columns.iter().map(|&j| self.values[j]).collect();
            let unfinished: Vec<_> = (0..incumbent.len())
                .filter(|&j| {
                    g.costs[j] > 0 && incumbent[j] > g.bounds[j].0
                        || g.costs[j] < 0 && incumbent[j] < g.bounds[j].1
                })
                .collect();
            if unfinished.is_empty() {
                break;
            }
            let mut indices = vec![unfinished[pick(unfinished.len())]];
            for _ in 0..1 + turn % 7 {
                let root = indices[pick(indices.len())];
                let incident = &adj[root];
                if incident.is_empty() {
                    continue;
                }
                let row = &g.bands[incident[pick(incident.len())]];
                if row.terms.len() > 48 {
                    continue;
                }
                indices.extend(row.terms.iter().map(|&(j, _)| j));
                indices.sort_unstable();
                indices.dedup();
                if indices.len() > 32 {
                    break;
                }
            }
            let local = restrict(g, &indices, &incumbent);
            if let Some(mut lp) = local.linear_model() {
                self.search(&mut lp, &local, &mut 36);
            }
        }
    }
}

// 两个分支覆盖当前位置所在的整个缺口。分别用局部 LP 的定点对偶下界，
// 二者取小即可提升父域下界。任何一侧算不出证书都沿用父界，不能凭数值
// 不可行或启发式失败删除一侧。每层只深入使二者较小值最大的分支对。
struct BoundSearch<'a> {
    bands: &'a [Band],
    opts: &'a SolveOptions,
    deadline: Option<Instant>,
    node_limit: u64,
    stats: &'a mut BbStats,
}

impl BoundSearch<'_> {
    fn strengthen(&mut self, lp: &mut Problem, x: &[f64], inherited: f64, depth: usize) -> f64 {
        if depth == 0 || self.stats.nodes >= self.node_limit || expired(self.deadline) {
            return inherited;
        }
        let mut splits: Vec<_> = self
            .bands
            .iter()
            .enumerate()
            .filter_map(|(i, b)| {
                b.split(x)
                    .map(|(fraction, left, right)| (fraction, i, left, right))
            })
            .collect();
        splits.sort_unstable_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        let mut best = None;
        let mut proved = inherited;
        for (_, i, left, right) in splits.into_iter().take(if depth == 3 { 6 } else { 2 }) {
            let mut children = Vec::new();
            let mut both = f64::INFINITY;
            for (op, rhs) in [(RelOp::Le, left), (RelOp::Ge, right)] {
                if self.stats.nodes >= self.node_limit || expired(self.deadline) {
                    return proved;
                }
                self.stats.nodes += 1;
                add_row(lp, &self.bands[i].terms, op, rhs.into());
                if let Some(sol) = relaxation(lp, self.opts, self.deadline, self.stats) {
                    let bound = lower_bound(lp, &sol.duals)
                        .unwrap_or(inherited)
                        .max(inherited);
                    both = both.min(bound);
                    children.push((op, rhs, sol.x, bound));
                } else {
                    both = both.min(inherited);
                }
                lp.constraints.pop();
            }
            if both > proved || best.is_none() {
                proved = both;
                best = Some((i, children));
            }
        }
        if let Some((i, children)) = best {
            if depth > 1 && children.len() == 2 {
                let mut both = f64::INFINITY;
                for (op, rhs, point, bound) in children {
                    add_row(lp, &self.bands[i].terms, op, rhs.into());
                    both = both.min(self.strengthen(lp, &point, bound, depth - 1));
                    lp.constraints.pop();
                }
                proved = proved.max(both);
            }
        }
        proved
    }
}

pub(crate) fn improve(
    p: &Problem,
    current: &[f64],
    opts: &SolveOptions,
    deadline: Option<Instant>,
    node_limit: u64,
    stats: &mut BbStats,
    accept: &mut dyn FnMut(&[f64]),
) -> Option<f64> {
    if stats.nodes >= node_limit || expired(deadline) {
        return None;
    }
    let repair = Repair::detect(p, deadline)?;
    let projection = Projection::build(p, &repair, deadline)?;
    let mut lp = projection.linear_model()?;
    if lp.constraints.len() > 2048 || projection.columns.len() > 1024 {
        return None;
    }
    let mut bound = None::<f64>;
    let mut known = BTreeSet::new();
    let mut root_point = Vec::new();
    stats.nodes += 1;
    for _ in 0..4 {
        let Some(sol) = relaxation(&lp, opts, deadline, stats) else {
            break;
        };
        if let Some(b) = lower_bound(&lp, &sol.duals) {
            bound = Some(bound.map_or(b, |old| old.max(b)));
        }
        root_point = sol.x.clone();
        if !opts.bb.root_cuts {
            break;
        }
        let old = lp.constraints.len();
        for row in &projection.bands {
            if expired(deadline) {
                break;
            }
            if let Some(cut) = row.cut(&sol.x, &projection.bounds) {
                if known.insert(cut.clone()) {
                    add_row(&mut lp, &cut.0, RelOp::Ge, cut.1);
                }
            }
            if lp.constraints.len() >= old + 256 {
                break;
            }
        }
        if lp.constraints.len() == old {
            break;
        }
    }
    if let Some(b) = bound {
        bound = Some(
            BoundSearch {
                bands: &projection.bands,
                opts,
                deadline,
                node_limit,
                stats,
            }
            .strengthen(&mut lp, &root_point, b, 3),
        );
    }
    if current.len() == p.vars.len() && current.iter().all(|&v| small_integer(v).is_some()) {
        let objective = p.obj.eval(current) + p.obj_offset;
        Neighborhood {
            original: p,
            repair: &repair,
            global: &projection,
            opts,
            deadline,
            remaining: node_limit - stats.nodes,
            global_bound: bound,
            stats,
            accept,
            values: current.iter().map(|&v| v as i32).collect(),
            best: if p.sense == Sense::Minimize {
                objective
            } else {
                -objective
            },
        }
        .run();
    }
    bound.map(|b| if p.sense == Sense::Minimize { b } else { -b })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn corner_cuts_preserve_all_periodic_slabs_and_finite_cycle_bounds() {
        let bounds = [(-2, 4), (0, 6), (-1, 3)];
        for l in -3..=2 {
            for width in 0..=5 {
                let b = Band {
                    terms: vec![(0, 1), (1, -1), (2, 1)],
                    low: l,
                    high: l + width,
                    period: 5,
                    count: 2,
                };
                let (_, _, lo, hi) = b.hull(&bounds).unwrap();
                let mut cuts = BTreeSet::new();
                for x in -2..=4 {
                    for y in 0..=6 {
                        for z in -1..=3 {
                            let point = [x as f64, y as f64, z as f64];
                            if let Some(c) = b.cut(&point, &bounds) {
                                cuts.insert(c);
                            }
                        }
                    }
                }
                for x in -4..=8 {
                    for y in 0..=12 {
                        for z in -2..=6 {
                            let activity = x - y + z;
                            if !(0..=2).any(|k| {
                                2 * (l + 5 * k) <= activity && activity <= 2 * (l + width + 5 * k)
                            }) {
                                continue;
                            }
                            assert!(2 * lo <= activity && activity <= 2 * hi);
                            assert!(b
                                .split(&[x as f64 / 2., y as f64 / 2., z as f64 / 2.])
                                .is_none());
                            let point = [x, y, z];
                            for (a, rhs) in &cuts {
                                assert!(
                                    a.iter().map(|&(j, t)| i64::from(t * point[j])).sum::<i64>()
                                        >= 2 * rhs
                                );
                            }
                        }
                    }
                }
            }
        }
        let exact = Band {
            terms: vec![(0, 1)],
            low: 1,
            high: 1,
            period: 4,
            count: 1,
        };
        assert!(exact.hull(&[(6, 12)]).is_none());
        assert_eq!(exact.hull(&[(2, 7)]), Some((2, 7, 5, 5)));
    }

    fn sample(sense: Sense, cycle: f64) -> Problem {
        let mut p = Problem::new("投影", sense);
        for (name, v, lo, hi) in [
            ("x", VarType::Continuous, 0., 4.),
            ("y", VarType::Continuous, 0., 2.),
            ("z", VarType::Continuous, 0., 4.),
            ("q", VarType::Integer, -1., 1.),
            ("w", VarType::Continuous, 2., 2.),
        ] {
            let j = p.add_var(name, v).unwrap();
            p.set_var_bounds(j, lo, hi).unwrap();
        }
        p.add_constraint(
            None,
            vec![(0, 1.), (1, -1.), (3, cycle), (4, 1.)],
            RelOp::Eq,
            5.,
        )
        .unwrap();
        p.add_constraint(None, vec![(0, 1.), (2, 1.)], RelOp::Eq, 4.)
            .unwrap();
        let sign = if sense == Sense::Minimize { 1. } else { -1. };
        p.set_objective(vec![(0, 2. * sign), (2, sign), (4, 2. * sign)], 3. * sign)
            .unwrap();
        p
    }

    #[test]
    fn global_certificate_and_verified_neighbors_follow_sense_and_budgets() {
        for sense in [Sense::Minimize, Sense::Maximize] {
            for cycle in [6., -6.] {
                let p = sample(sense, cycle);
                let seed = vec![4., 1., 0., 0., 2.];
                let sign = if sense == Sense::Minimize { 1. } else { -1. };
                // y=x−3，整数周期取 ±1 要么 y 超界，要么 x 超界。x=3 时目标最优。
                let optimum = 14. * sign;
                let mut s = BbStats::default();
                let mut found = Vec::new();
                let mut record = |x: &[f64]| {
                    assert!(p.constraints.iter().all(|r| r.expr.eval(x) == r.rhs));
                    assert!(x
                        .iter()
                        .zip(&p.vars)
                        .all(|(&v, var)| var.lb <= v && v <= var.ub));
                    found.push(p.obj.eval(x) + p.obj_offset);
                };
                let opts = SolveOptions::default();
                assert_eq!(
                    improve(
                        &p,
                        &seed,
                        &opts,
                        Some(Instant::now()),
                        2,
                        &mut s,
                        &mut record
                    ),
                    None
                );
                assert_eq!(s.nodes, 0);
                let b = improve(&p, &seed, &opts, None, 1, &mut s, &mut record).unwrap();
                assert_eq!(s.nodes, 1);
                assert!((b - optimum).abs() < 1e-7);
                assert_eq!(
                    improve(&p, &seed, &opts, None, 1, &mut s, &mut record),
                    None
                );
                let _ = improve(&p, &seed, &opts, None, 4, &mut s, &mut record);
                assert!(s.nodes <= 4);
                assert_eq!(found, vec![optimum]);
            }
        }
    }

    #[test]
    fn both_halves_require_a_certificate_before_raising_the_parent_bound() {
        let bands = [Band {
            terms: vec![(0, 1)],
            low: 0,
            high: 1,
            period: 3,
            count: 1,
        }];
        let mut lp = Problem::new("缺口", Sense::Minimize);
        for name in ["x", "z"] {
            let j = lp.add_var(name, VarType::Continuous).unwrap();
            lp.set_var_bounds(j, 0., 4.).unwrap();
        }
        lp.set_objective(vec![(1, 1.)], 0.).unwrap();
        lp.add_constraint(None, vec![(0, 1.), (1, 1.)], RelOp::Ge, 2.)
            .unwrap();
        lp.add_constraint(None, vec![(0, -1.), (1, 1.)], RelOp::Ge, -2.)
            .unwrap();
        // x=2 在线性松弛中让 |x−2|=0，两个周期区间上该目标至少为 1。
        for (limit, expected) in [(0, 0.), (1, 0.), (2, 1.)] {
            let opts = SolveOptions::default();
            let mut stats = BbStats::default();
            let b = BoundSearch {
                bands: &bands,
                opts: &opts,
                deadline: None,
                node_limit: limit,
                stats: &mut stats,
            }
            .strengthen(&mut lp, &[2., 0.], 0., 1);
            assert_eq!(b, expected);
            assert_eq!(stats.nodes, limit);
            assert_eq!(lp.constraints.len(), 2);
        }
    }

    #[test]
    fn wrong_signed_duals_fractional_costs_and_private_cycle_costs_are_safe() {
        let mut lp = Problem::new("证书", Sense::Minimize);
        let j = lp.add_var("x", VarType::Continuous).unwrap();
        lp.set_var_bounds(j, -2., 3.).unwrap();
        lp.set_objective(vec![(j, 2.)], -1.).unwrap();
        lp.add_constraint(None, vec![(j, 1.)], RelOp::Ge, 1.)
            .unwrap();
        assert_eq!(lower_bound(&lp, &[-2.]), Some(-5.));
        assert!(lower_bound(&lp, &[2.00000003]).unwrap() <= 1.);
        assert_eq!(lower_bound(&lp, &[f64::NAN]), None);
        let mut p = sample(Sense::Minimize, 6.);
        let r = Repair::detect(&p, None).unwrap();
        p.obj.add_term(3, 1.);
        assert!(Projection::build(&p, &r, None).is_none());
        p.obj.add_term(3, -1.);
        p.obj_offset = 0.5;
        assert!(Projection::build(&p, &r, None).is_none());
    }
}
