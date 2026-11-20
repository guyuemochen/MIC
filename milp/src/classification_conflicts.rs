//! 非负特征的线性分类模型。相同特征上，较小样本判正且较大样本判负时，
//! 至少支付一笔二元误差费用。用最小割计算这些冲突的费用下界，再用已核验
//! 的可行费用固定误差开关。所有删行、合并都在固定后进行，固定前检查整张原模型。
use crate::model::{Problem, RelOp, Sense, Tolerances, VarType};
use crate::simplex::{Bounds, LpStatus, Simplex};
use std::{
    collections::{BTreeMap, VecDeque},
    time::Instant,
};

const SCALE: f64 = (1u64 << 36) as f64;
const GRID: f64 = (1u64 << 42) as f64;

fn expired(end: Option<Instant>) -> bool {
    end.is_some_and(|t| Instant::now() >= t)
}
fn subset(a: u64, b: u64) -> bool {
    a & !b == 0
}

#[derive(Clone)]
struct Observation {
    row: usize,
    private: usize,
    features: u64,
    positive: bool,
    binary: bool,
    cost: f64,
}
struct Data {
    cols: Vec<usize>,
    minus: usize,
    plus: usize,
    samples: Vec<Observation>,
    margin: f64,
    release: f64,
}

impl Data {
    fn detect(p: &Problem, tol: Tolerances) -> Option<Self> {
        if p.sense != Sense::Minimize || p.constraints.len() > 640 || p.obj_offset != 0.0 {
            return None;
        }
        let mut eqs = p.constraints.iter().filter(|r| r.op == RelOp::Eq);
        let norm = eqs.next()?;
        if eqs.next().is_some() || norm.rhs != 1.0 || norm.expr.terms.len() > 120 {
            return None;
        }
        let mut shared = vec![false; p.vars.len()];
        let mut cols = Vec::new();
        for &(j, a) in &norm.expr.terms {
            let v = &p.vars[j];
            if a != 1.0
                || v.var_type != VarType::Continuous
                || v.lb != 0.0
                || (v.ub < 1.0 && v.ub != 0.0)
            {
                return None;
            }
            shared[j] = true;
            if v.ub != 0.0 {
                cols.push(j);
            }
        }
        if cols.len() < 2 || cols.len() > 60 {
            return None;
        }
        let mut index = vec![None; p.vars.len()];
        for (k, &j) in cols.iter().enumerate() {
            index[j] = Some(k);
        }
        let mut cost = vec![0.; p.vars.len()];
        for &(j, c) in &p.obj.terms {
            if shared[j] || c <= 0. || c > 1e4 {
                return None;
            }
            cost[j] = c;
        }
        let mut visited = shared.clone();
        let mut samples = Vec::new();
        let mut margin = None;
        let mut release = None;
        let mut neg_col = None;
        let mut every_positive = (1u64 << cols.len()) - 1;
        for (r, row) in p
            .constraints
            .iter()
            .enumerate()
            .filter(|(_, r)| r.op != RelOp::Eq)
        {
            let pos = row.op == RelOp::Ge;
            let eps = if pos { row.rhs } else { -row.rhs };
            if !(eps > 0. && eps < 0.1 && eps > 100. * (tol.feasibility + tol.integrality))
                || margin.replace(eps).is_some_and(|e| e != eps)
            {
                return None;
            }
            let mut private = None;
            let mut neg = None;
            let mut features = 0;
            for &(j, a) in &row.expr.terms {
                if shared[j] {
                    if let Some(k) = index[j] {
                        if a == 1.0 {
                            features |= 1u64 << k;
                        } else if a != -1.0 || neg.replace(j).is_some() {
                            return None;
                        }
                    }
                } else {
                    let m = if pos { a } else { -a };
                    let v = &p.vars[j];
                    if std::mem::replace(&mut visited[j], true)
                        || private.replace(j).is_some()
                        || !(m >= 1. && m <= 2.)
                        || release.replace(m).is_some_and(|z| z != m)
                        || !matches!(v.var_type, VarType::Continuous | VarType::Binary)
                        || v.lb != 0.0
                        || v.ub != 1.0
                        || cost[j] <= 0.0
                    {
                        return None;
                    }
                }
            }
            let j = private?;
            let n = neg?;
            if neg_col.replace(n).is_some_and(|k| k != n) {
                return None;
            }
            every_positive &= features;
            samples.push(Observation {
                row: r,
                private: j,
                features,
                positive: pos,
                binary: p.vars[j].var_type == VarType::Binary,
                cost: cost[j],
            });
        }
        let count = samples.iter().filter(|s| s.positive && s.binary).count();
        if !visited.iter().all(|&v| v) || every_positive == 0 || !(2..=32).contains(&count) {
            return None;
        }
        let plus = cols[every_positive.trailing_zeros() as usize];
        Some(Self {
            cols,
            minus: neg_col?,
            plus,
            samples,
            margin: margin?,
            release: release?,
        })
    }

    fn residual(&self, s: &Observation, score: f64) -> f64 {
        self.margin + if s.positive { -score } else { score }
    }

    // 放松后的 LP 仅产生候选方向。将方向投到二进制定点格，增大一点正确分类
    // 的间隔，再由两个相反的截距列补齐总和。此时原行中的特征和可精确计算。
    fn candidate(&self, p: &Problem, guess: &[f64]) -> Option<Vec<f64>> {
        let mut x = vec![0.; p.vars.len()];
        let factor = 1.000001;
        let mut remaining = GRID as i64;
        for (k, &j) in self.cols.iter().enumerate() {
            if j == self.plus || j == self.minus {
                continue;
            }
            let n = (guess[k].max(0.) * factor * GRID).round() as i64;
            remaining = remaining.checked_sub(n)?;
            x[j] = n as f64 / GRID;
        }
        let find = |j| self.cols.binary_search(&j).ok();
        let mut diff =
            ((guess[find(self.plus)?] - guess[find(self.minus)?]) * factor * GRID).round() as i64;
        if (diff ^ remaining) & 1 != 0 {
            diff += 1;
        }
        if remaining < diff.abs() {
            return None;
        }
        x[self.plus] = ((remaining + diff) / 2) as f64 / GRID;
        x[self.minus] = ((remaining - diff) / 2) as f64 / GRID;
        for s in &self.samples {
            let score = p.constraints[s.row].expr.eval(&x);
            let d = self.residual(s, score);
            // 向不利方向外扩，选出的松弛量在实数原行上也足够。
            let hi = if d <= 0.
                && (if s.positive {
                    score >= self.margin
                } else {
                    score <= -self.margin
                }) {
                0.
            } else {
                (d.next_up() / self.release).next_up()
            };
            x[s.private] = if s.binary {
                f64::from(hi > 0.)
            } else {
                hi.max(0.)
            };
            if x[s.private] > 1. || self.release * x[s.private] < d.next_up() && hi > 0. {
                return None;
            }
        }
        for r in &p.constraints {
            let a = r.expr.eval(&x);
            if !match r.op {
                RelOp::Le => a <= r.rhs,
                RelOp::Ge => a >= r.rhs,
                RelOp::Eq => a == r.rhs,
            } {
                return None;
            }
        }
        Some(x)
    }

    fn seed(
        &self,
        p: &Problem,
        bits: &[bool],
        end: Option<Instant>,
        tol: Tolerances,
        iterations: &mut u64,
    ) -> Option<Vec<f64>> {
        let mut lp = Problem::new("", Sense::Minimize);
        for &j in &self.cols {
            lp.add_var(&p.vars[j].name, VarType::Continuous).ok()?;
        }
        let mut obj = Vec::new();
        for (i, s) in self.samples.iter().enumerate() {
            if expired(end) {
                return None;
            }
            if s.binary && bits[i] {
                continue;
            }
            let mut terms = Vec::new();
            for (k, &j) in self.cols.iter().enumerate() {
                if s.features & (1 << k) != 0 {
                    terms.push((k, 1.));
                }
                if j == self.minus {
                    terms.push((k, -1.));
                }
            }
            if !s.binary || !s.positive {
                let j = lp
                    .add_var(&p.vars[s.private].name, VarType::Continuous)
                    .ok()?;
                lp.vars[j].ub = 1.;
                obj.push((j, s.cost * if s.binary { 0.5 / self.margin } else { 1. }));
                terms.push((
                    j,
                    if s.positive {
                        self.release
                    } else {
                        -self.release
                    },
                ));
            }
            let row = &p.constraints[s.row];
            lp.add_constraint(None, terms, row.op, row.rhs).ok()?;
        }
        lp.add_constraint(
            None,
            (0..self.cols.len()).map(|j| (j, 1.)).collect(),
            RelOp::Eq,
            1.,
        )
        .ok()?;
        lp.set_objective(obj, 0.).ok()?;
        let mut simplex = Simplex::new(&lp, tol).ok()?;
        simplex.set_time_deadline(end);
        let sol = simplex.solve(&Bounds::from_problem(&lp));
        *iterations += simplex.iterations();
        if sol.status != LpStatus::Optimal {
            return None;
        }
        self.candidate(p, &sol.x)
    }
}

#[derive(Clone, Copy)]
struct Edge {
    dest: usize,
    back: usize,
    cap: i64,
}
#[derive(Clone)]
struct Network {
    arcs: Vec<Vec<Edge>>,
    source: usize,
    sink: usize,
    inf: i64,
}

impl Network {
    fn arc(&mut self, a: usize, b: usize, n: i64) {
        let (u, v) = (self.arcs[a].len(), self.arcs[b].len());
        self.arcs[a].push(Edge {
            dest: b,
            back: v,
            cap: n,
        });
        self.arcs[b].push(Edge {
            dest: a,
            back: u,
            cap: 0,
        });
    }
    fn new(data: &Data) -> Self {
        let n = data.samples.len();
        let mut f = Self {
            arcs: vec![Vec::new(); n + 2],
            source: n,
            sink: n + 1,
            inf: 1,
        };
        for (i, s) in data.samples.iter().enumerate().filter(|(_, s)| s.binary) {
            // 乘二的幂精确，逐列向下取整后，即使割不唯一也给出保守费用。
            let c = (s.cost * SCALE).floor() as i64;
            f.inf += c;
            if s.positive {
                f.arc(n, i, c)
            } else {
                f.arc(i, n + 1, c)
            }
        }
        for (i, s) in data
            .samples
            .iter()
            .enumerate()
            .filter(|(_, s)| s.positive && s.binary)
        {
            for (j, t) in data
                .samples
                .iter()
                .enumerate()
                .filter(|(_, t)| !t.positive && t.binary)
            {
                if subset(s.features, t.features) {
                    f.arc(i, j, f.inf);
                }
            }
        }
        f
    }
    fn force(&mut self, i: usize, source_side: bool) {
        if source_side {
            self.arc(self.source, i, self.inf)
        } else {
            self.arc(i, self.sink, self.inf)
        }
    }
    fn flow(&mut self) -> (i64, Vec<bool>) {
        let n = self.arcs.len();
        let mut total = 0;
        loop {
            let mut level = vec![usize::MAX; n];
            level[self.source] = 0;
            let mut q = VecDeque::from([self.source]);
            while let Some(i) = q.pop_front() {
                for e in &self.arcs[i] {
                    if e.cap > 0 && level[e.dest] == usize::MAX {
                        level[e.dest] = level[i] + 1;
                        q.push_back(e.dest);
                    }
                }
            }
            if level[self.sink] == usize::MAX {
                return (total, level.iter().map(|&d| d != usize::MAX).collect());
            }
            let mut it = vec![0; n];
            loop {
                let a = self.push(self.source, self.inf, &level, &mut it);
                if a == 0 {
                    break;
                }
                total += a;
            }
        }
    }
    fn push(&mut self, i: usize, cap: i64, level: &[usize], it: &mut [usize]) -> i64 {
        if i == self.sink {
            return cap;
        }
        while it[i] < self.arcs[i].len() {
            let e = self.arcs[i][it[i]];
            if e.cap > 0 && level[e.dest] == level[i] + 1 {
                let f = self.push(e.dest, cap.min(e.cap), level, it);
                if f > 0 {
                    self.arcs[i][it[i]].cap -= f;
                    self.arcs[e.dest][e.back].cap += f;
                    return f;
                }
            }
            it[i] += 1;
        }
        0
    }
}

pub(crate) struct Reformulation {
    pub(crate) problem: Problem,
    pub(crate) seed: Vec<f64>,
    pub(crate) lp_iterations: u64,
    affine: Vec<(usize, usize, f64, usize)>,
}

impl Reformulation {
    pub(crate) fn project_seed(&self) -> Option<Vec<f64>> {
        let mut work = self.seed.clone();
        work.resize(self.problem.vars.len(), 0.);
        for &(_, j, _, g) in &self.affine {
            work[g] += work[j];
            work[j] = 0.;
        }
        work.iter()
            .zip(&self.problem.vars)
            .all(|(&v, b)| b.lb <= v && v <= b.ub)
            .then_some(work)
    }

    pub(crate) fn restore(&self, orig: &Problem, point: &[f64]) -> Option<Vec<f64>> {
        let mut x = point.get(..orig.vars.len())?.to_vec();
        for &(row, j, m, _) in &self.affine {
            x[j] = 0.;
            let r = &orig.constraints[row];
            let v = (r.rhs - r.expr.eval(&x)) / m;
            if !v.is_finite() || !(-1e-7..=1.0000001).contains(&v) {
                return None;
            }
            x[j] = v.clamp(0., 1.);
        }
        Some(x)
    }
}

fn add_cut(w: &mut Problem, terms: Vec<(usize, f64)>, op: RelOp, rhs: f64) -> Option<()> {
    // 多个同型样本会产生二次规模的传递行；过密时直接交回通用搜索。
    if w.constraints.len() >= 4096 {
        return None;
    }
    w.add_constraint(None, terms, op, rhs).ok()?;
    Some(())
}

pub(crate) fn reformulate(
    p: &Problem,
    tol: Tolerances,
    end: Option<Instant>,
) -> Option<Reformulation> {
    let data = Data::detect(p, tol)?;
    let obs = &data.samples;
    let mut net = Network::new(&data);
    if expired(end) {
        return None;
    }
    let (lower, side) = net.flow();
    let flags = |cut: &[bool]| {
        obs.iter()
            .enumerate()
            .map(|(i, s)| s.binary && (s.positive != cut[i]))
            .collect::<Vec<_>>()
    };
    let mut second = None;
    for (i, _) in obs
        .iter()
        .enumerate()
        .filter(|(_, s)| s.positive && s.binary)
    {
        if expired(end) {
            return None;
        }
        let mut f = net.clone();
        f.force(i, !side[i]);
        let (extra, cut) = f.flow();
        if second.as_ref().is_none_or(|(best, _)| extra < *best) {
            second = Some((extra, cut));
        }
    }
    let mut seed = None;
    let mut upper = f64::INFINITY;
    let mut iterations = 0;
    for cut in [side, second?.1] {
        if let Some(x) = data.seed(p, &flags(&cut), end, tol, &mut iterations) {
            let ub = p
                .obj
                .terms
                .iter()
                .fold(0f64, |sum, &(j, c)| (sum + (x[j] * c).next_up()).next_up());
            if ub < upper {
                upper = ub;
                seed = Some(x);
            }
        }
    }
    let seed = seed?;
    let budget = (upper * SCALE).ceil() as i64;
    let mut fixed = vec![None; obs.len()];
    for (i, s) in obs.iter().enumerate().filter(|(_, s)| s.binary) {
        if expired(end) {
            return None;
        }
        for penalty in [false, true] {
            let side = s.positive != penalty;
            let mut f = net.clone();
            // 同标签的误差可取逐点最小值。假定此开关的状态时，一并计入
            // 单调性所迫使的同标签开关；比较的两行共享同一间隔和松弛系数。
            for (j, t) in obs
                .iter()
                .enumerate()
                .filter(|(_, t)| t.binary && t.positive == s.positive)
            {
                if if side {
                    subset(s.features, t.features)
                } else {
                    subset(t.features, s.features)
                } {
                    f.force(j, side);
                }
            }
            if lower + f.flow().0 > budget {
                fixed[i] = Some(!penalty);
                break;
            }
        }
    }
    if expired(end) {
        return None;
    }
    build(p, data, fixed, seed, iterations, end)
}

fn build(
    p: &Problem,
    data: Data,
    fixed: Vec<Option<bool>>,
    seed: Vec<f64>,
    iterations: u64,
    end: Option<Instant>,
) -> Option<Reformulation> {
    let obs = &data.samples;
    let right: Vec<_> = (0..obs.len())
        .filter(|&i| fixed[i] == Some(false))
        .collect();
    let positive: Vec<_> = right.iter().copied().filter(|&i| obs[i].positive).collect();
    if positive.is_empty() || right.len() == positive.len() {
        return None;
    }
    // 至少一条正确负样本迫使负截距 >= 间隔，因此任一得分 <= 1-2*间隔。
    // 至少一条正确正样本又迫使负截距 <= (1-间隔)/2。因此两侧任一违反量
    // 均小于 1，开关取 1 的原行可安全删除，连续补位也一定在上界之内。
    let unknown: Vec<_> = (0..obs.len())
        .filter(|&i| obs[i].binary && fixed[i].is_none())
        .collect();
    let mut w = p.clone();
    let mut keep = vec![true; p.constraints.len()];
    let mut aggregate = BTreeMap::<(bool, u64), Vec<usize>>::new();
    for (i, s) in obs.iter().enumerate() {
        if let Some(v) = fixed[i] {
            w.vars[s.private].lb = f64::from(v);
            w.vars[s.private].ub = f64::from(v);
        }
        if fixed[i] == Some(true) {
            keep[s.row] = false;
            continue;
        }
        let implies = |k: usize| {
            let t = &obs[k];
            if t.positive {
                subset(t.features, s.features)
            } else {
                subset(s.features, t.features)
            }
        };
        if fixed[i] == Some(false) {
            if right.iter().any(|&k| {
                obs[k].positive == s.positive
                    && implies(k)
                    && (obs[k].features != s.features || k < i)
            }) {
                keep[s.row] = false;
            }
        } else if !s.binary {
            if right
                .iter()
                .any(|&k| obs[k].positive == s.positive && implies(k))
            {
                w.vars[s.private].ub = 0.;
                keep[s.row] = false;
            } else if right
                .iter()
                .any(|&k| obs[k].positive != s.positive && implies(k))
            {
                keep[s.row] = false;
                aggregate
                    .entry((s.positive, s.cost.to_bits()))
                    .or_default()
                    .push(i);
            }
        }
    }
    w.constraints = std::mem::take(&mut w.constraints)
        .into_iter()
        .zip(keep.iter())
        .filter_map(|(r, &k)| k.then_some(r))
        .collect();
    for &i in &unknown {
        if expired(end) {
            return None;
        }
        let a = &obs[i];
        for &j in &unknown {
            let b = &obs[j];
            if i == j || !subset(a.features, b.features) {
                continue;
            }
            if a.positive && !b.positive {
                add_cut(
                    &mut w,
                    vec![(a.private, 1.), (b.private, 1.)],
                    RelOp::Ge,
                    1.,
                )?;
            } else if a.positive == b.positive {
                if unknown.iter().any(|&k| {
                    k != i
                        && k != j
                        && obs[k].positive == a.positive
                        && obs[k].features != a.features
                        && obs[k].features != b.features
                        && subset(a.features, obs[k].features)
                        && subset(obs[k].features, b.features)
                }) {
                    continue;
                }
                add_cut(
                    &mut w,
                    vec![(a.private, 1.), (b.private, -1.)],
                    if a.positive { RelOp::Ge } else { RelOp::Le },
                    0.,
                )?;
            }
        }
    }
    // 两条正确正样本的特征和逐维不大于两条负样本的和时，两条负样本不可能都正确。
    // 对二值特征，两数之和可由并集和交集精确表示。
    let pairs: Vec<_> = positive
        .iter()
        .enumerate()
        .flat_map(|(k, &i)| {
            positive[k..].iter().map(move |&j| {
                (
                    obs[i].features | obs[j].features,
                    obs[i].features & obs[j].features,
                )
            })
        })
        .collect();
    let negatives: Vec<_> = (0..obs.len())
        .filter(|&i| obs[i].binary && !obs[i].positive && fixed[i] != Some(true))
        .collect();
    if negatives.len() * negatives.len() * pairs.len() > 60_000_000 {
        return None;
    }
    for (k, &i) in negatives.iter().enumerate() {
        if expired(end) {
            return None;
        }
        for &j in &negatives[k + 1..] {
            let a = obs[i].features;
            let b = obs[j].features;
            if pairs
                .iter()
                .any(|&(any, both)| subset(any, a | b) && subset(both, a & b))
            {
                add_cut(
                    &mut w,
                    vec![(obs[i].private, 1.), (obs[j].private, 1.)],
                    RelOp::Ge,
                    1.,
                )?;
            }
        }
    }
    let mut affine = Vec::new();
    for ((pos, cost), list) in aggregate {
        let mut tail = &list[..];
        while !tail.is_empty() {
            // 组长取二的幂，常数项放大精确；同价列相加，目标无需做浮点代入。
            let n = 1usize << tail.len().ilog2();
            let (part, rest) = tail.split_at(n);
            tail = rest;
            if n == 1 {
                let s = &obs[part[0]];
                // 原行的误差在该方向始终非负，因此最优修复恰取等号。
                add_cut(
                    &mut w,
                    p.constraints[s.row].expr.terms.clone(),
                    RelOp::Eq,
                    p.constraints[s.row].rhs,
                )?;
            } else {
                let g = w
                    .add_var(&format!("errsum{}", w.vars.len()), VarType::Continuous)
                    .ok()?;
                w.vars[g].ub = n as f64;
                w.obj.add_term(g, f64::from_bits(cost));
                let sign = if pos { 1. } else { -1. };
                let coef = sign * data.release;
                let mut sum = BTreeMap::<usize, f64>::new();
                for &i in part {
                    let s = &obs[i];
                    w.vars[s.private].ub = 0.;
                    affine.push((s.row, s.private, coef, g));
                    for &(j, a) in &p.constraints[s.row].expr.terms {
                        if j != s.private && p.vars[j].ub != 0. {
                            *sum.entry(j).or_default() += a;
                        }
                    }
                }
                sum.insert(g, coef);
                add_cut(
                    &mut w,
                    sum.into_iter().collect(),
                    RelOp::Eq,
                    sign * data.margin * n as f64,
                )?;
            }
        }
    }
    Some(Reformulation {
        problem: w,
        seed,
        lp_iterations: iterations,
        affine,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn example() -> Problem {
        let mut p = Problem::new("误差分类", Sense::Minimize);
        for j in 0..4 {
            p.add_var(&format!("w{j}"), VarType::Continuous).unwrap();
        }
        let mut obj = Vec::new();
        for (i, (mask, pos, c, bt)) in [
            (1, true, 0.07, true),
            (3, true, 0.1, true),
            (0, false, 0.2, true),
            (2, false, 0.25, true),
            (0, true, 0.01, true),
            (3, false, 0.02, true),
            (1, false, 0.001, false),
            (1, false, 0.001, false),
        ]
        .into_iter()
        .enumerate()
        {
            let j = p
                .add_var(
                    &format!("e{i}"),
                    if bt {
                        VarType::Binary
                    } else {
                        VarType::Continuous
                    },
                )
                .unwrap();
            p.vars[j].ub = 1.;
            let mut terms = vec![(2, 1.), (3, -1.)];
            for k in 0..2 {
                if mask & (1 << k) != 0 {
                    terms.push((k, 1.));
                }
            }
            let sign = if pos { 1.0 } else { -1.0 };
            terms.push((j, sign * 1.25));
            p.add_constraint(
                None,
                terms,
                if pos { RelOp::Ge } else { RelOp::Le },
                sign / 1024.,
            )
            .unwrap();
            obj.push((j, c));
        }
        p.set_objective(obj, 0.).unwrap();
        p.add_constraint(None, (0..4).map(|j| (j, 1.)).collect(), RelOp::Eq, 1.)
            .unwrap();
        p
    }

    fn feasible(p: &Problem, x: &[f64]) -> bool {
        x.len() == p.vars.len()
            && x.iter()
                .zip(&p.vars)
                .all(|(&v, b)| v.is_finite() && b.lb <= v && v <= b.ub)
            && p.constraints.iter().all(|r| {
                let d = r.expr.eval(x) - r.rhs;
                match r.op {
                    RelOp::Le => d <= 1e-12,
                    RelOp::Ge => d >= -1e-12,
                    RelOp::Eq => d.abs() <= 1e-12,
                }
            })
    }

    #[test]
    fn combines_active_equal_price_slacks_and_keeps_the_optimum() {
        let p = example();
        let f = reformulate(&p, Default::default(), None).unwrap();
        assert!(feasible(&p, &f.seed));
        let seed = f.project_seed().unwrap();
        assert!(feasible(&f.problem, &seed));
        assert_eq!(f.affine.len(), 2);
        let mut solver = Simplex::new(&f.problem, Default::default()).unwrap();
        let sol = solver.solve(&Bounds::from_problem(&f.problem));
        assert_eq!(sol.status, LpStatus::Optimal);
        let x = f.restore(&p, &sol.x).unwrap();
        assert!(feasible(&p, &x));
        // x0=2ε、截距为 -ε 时，前四个二元样本全对；剩余两个相反标签及
        // 两个连续误差不得不付费。任意错分前四个样本的费用均大于该解。
        let optimum = 0.03 + 2. * 0.001 * (2. / 1024.) / 1.25;
        assert!((p.obj.eval(&x) - optimum).abs() < 1e-12);
        assert!((sol.objective - optimum).abs() < 1e-12);
    }

    #[test]
    fn unmatched_bounds_coefficients_and_interruption_fall_back() {
        let p = example();
        assert!(reformulate(&p, Default::default(), Some(Instant::now())).is_none());
        let mut changed = p.clone();
        changed.constraints[0].expr.add_term(0, 0.5);
        assert!(reformulate(&changed, Default::default(), None).is_none());
        let mut changed = p.clone();
        changed.vars[0].ub = 0.5;
        assert!(reformulate(&changed, Default::default(), None).is_none());
        let mut changed = p.clone();
        changed.obj.add_term(0, 1.);
        assert!(reformulate(&changed, Default::default(), None).is_none());
    }
    #[test]
    fn small_error_prices_keep_the_required_slack() {
        // 两个高价开关都应为零。于是同一行必须花一笔二元费用，连续误差至少为 2*间隔/1.25。
        // 把价格直接加成约束时，预处理容易放大这些相差九个数量级的系数。
        let p = crate::instance::parse_instance(
            "name small-prices
sense min
vars
w continuous
m continuous
p continuous
a binary
b binary
c binary
e continuous
obj
100 a 99 b 0.0000023 c 0.00000007 e
cons
1 w 1 m 1 p = 1
1 w -1 m 1 p 1.25 a >= 0.00048828125
-1 m 1 p -1.25 b <= -0.00048828125
-1 m 1 p 1.25 c >= 0.00048828125
-1 m 1 p 1.25 e >= 0.00048828125
bounds
e <= 1
end",
        )
        .unwrap();
        assert!(reformulate(&p, Default::default(), None).is_some());
        let result = crate::solver::solve_default(&p).unwrap();
        assert_eq!(result.status, crate::model::SolveStatus::Optimal);
        assert!(feasible(&p, &result.values.unwrap()));
        let expected = 2.3e-6 + 7e-8 * 2.0 * 0.00048828125 / 1.25;
        assert!((result.objective.unwrap() - expected).abs() < 1e-15);
    }
}
