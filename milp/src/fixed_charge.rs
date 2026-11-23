//! 带固定启用费用的连续分配。用等式回代消去库存、欠交等辅助列后，
//! 在各项余量允许的范围内增加流量，同时为所涉及的开关付费。
//! 候选解只给分支定界提供上界或下界的可行一侧，不能用来证明最优。
use crate::capacity_gate::{at_one, is_binary};
use crate::model::{LinearExpr, Problem, RelOp, Sense, VarType};
use std::collections::{BTreeMap, VecDeque};

type Terms = Vec<(usize, f64)>;

struct Allocation {
    // 每个原连续变量的常数与流量系数。
    affine: Vec<(f64, Terms)>,
    capacity: Vec<f64>,
    flow_limits: Vec<f64>,
    resources: Vec<Terms>,
    gates: Vec<Vec<usize>>,
    conflicts: Vec<Vec<usize>>,
    gain: Vec<f64>,
    fee: Vec<f64>,
}

fn charge_projection(expr: &LinearExpr, affine: &[(f64, Terms)], work: &mut usize) -> Option<()> {
    let total = expr
        .terms
        .iter()
        .try_fold(0usize, |s, &(j, _)| s.checked_add(affine[j].1.len() + 1))?;
    if total > 65_536 {
        return None;
    }
    *work = work.checked_sub(total)?;
    Some(())
}

fn project(
    expr: &LinearExpr,
    affine: &[(f64, Terms)],
    factor: f64,
    skip: Option<usize>,
    scratch: &mut [f64],
    touched: &mut Vec<usize>,
) -> (f64, Terms) {
    let mut offset = 0.;
    for &(j, a) in &expr.terms {
        if Some(j) == skip {
            continue;
        }
        let a = factor * a;
        offset += a * affine[j].0;
        for &(k, b) in &affine[j].1 {
            if scratch[k] == 0. {
                touched.push(k);
            }
            scratch[k] += a * b;
        }
    }
    touched.sort_unstable();
    touched.dedup();
    let terms = touched
        .drain(..)
        .filter_map(|k| {
            let a = std::mem::take(&mut scratch[k]);
            (a != 0.).then_some((k, a))
        })
        .collect();
    (offset, terms)
}

impl Allocation {
    fn extract(p: &Problem, stop: &mut dyn FnMut() -> bool) -> Option<Self> {
        let n = p.vars.len();
        // 限制矩阵、邻接表、等式消元的总开销，也限制单次不能打断的行展开。
        let mut work = 4_000_000usize;
        work = work.checked_sub(n.checked_add(p.constraints.len())?.checked_mul(8)?)?;
        if stop() {
            return None;
        }
        let signed = if p.sense == Sense::Minimize { 1. } else { -1. };
        if !p.vars.iter().any(is_binary)
            || !p.vars.iter().any(|v| v.var_type == VarType::Continuous)
            || p.vars.iter().any(|v| {
                v.var_type != VarType::Continuous && (!is_binary(v) || v.lb != 0. || v.ub != 1.)
            })
        {
            return None;
        }
        let mut fee = vec![0.; n];
        for &(j, a) in &p.obj.terms {
            if is_binary(&p.vars[j]) {
                fee[j] = signed * a;
                if fee[j] < 0. {
                    return None;
                }
            }
        }
        if !fee.iter().any(|&a| a > 0.) {
            return None;
        }
        let mut flow = vec![false; n];
        let mut has_paid_gate = false;
        for (i, r) in p.constraints.iter().enumerate() {
            if i % 128 == 0 && stop() {
                return None;
            }
            if r.expr.terms.len() > 65_536 {
                return None;
            }
            work = work.checked_sub(r.expr.terms.len())?;
            if let Some((j, _)) = at_one(p, r) {
                has_paid_gate |= fee[j] > 0.;
                for &(k, _) in &r.expr.terms {
                    if k != j {
                        flow[k] = true;
                    }
                }
            }
        }
        if !has_paid_gate {
            return None;
        }
        let cols: Vec<_> = (0..n).filter(|&j| flow[j]).collect();
        let m = cols.len();
        let mut affine = vec![(0., vec![]); n];
        let mut known = vec![false; n];
        let mut map = vec![usize::MAX; n];
        for (k, &j) in cols.iter().enumerate() {
            if p.vars[j].lb != 0. {
                return None;
            }
            affine[j].1.push((k, 1.));
            known[j] = true;
            map[j] = k;
        }
        let mut scratch = vec![0.; m];
        let mut touched = vec![];
        // 只在一个未知列已经可以回代时访问该行。即使等式按反拓扑序给出，
        // 每个输入非零元也只登记一次，常数与系数的实际展开另计工作量。
        let mut unknown = vec![0; p.constraints.len()];
        let mut dependents = vec![Vec::new(); n];
        let mut ready = VecDeque::new();
        let mut use_eq = vec![false; p.constraints.len()];
        for (i, r) in p
            .constraints
            .iter()
            .enumerate()
            .filter(|(_, r)| r.op == RelOp::Eq)
        {
            if i % 128 == 0 && stop() {
                return None;
            }
            for &(j, _) in &r.expr.terms {
                if is_binary(&p.vars[j]) {
                    return None;
                }
                if !known[j] {
                    unknown[i] += 1;
                    dependents[j].push(i);
                }
            }
            if unknown[i] == 1 {
                ready.push_back(i);
            }
        }
        while let Some(i) = ready.pop_front() {
            if stop() {
                return None;
            }
            if unknown[i] != 1 {
                continue;
            }
            let r = &p.constraints[i];
            let &(j, a) = r.expr.terms.iter().find(|&&(k, _)| !known[k])?;
            charge_projection(&r.expr, &affine, &mut work)?;
            let (off, terms) = project(
                &r.expr,
                &affine,
                -1. / a,
                Some(j),
                &mut scratch,
                &mut touched,
            );
            affine[j] = (r.rhs / a + off, terms);
            known[j] = true;
            use_eq[i] = true;
            for &other in &dependents[j] {
                unknown[other] -= 1;
                if unknown[other] == 1 {
                    ready.push_back(other);
                }
            }
        }
        if p.vars
            .iter()
            .enumerate()
            .any(|(j, v)| v.var_type == VarType::Continuous && !known[j])
        {
            return None;
        }

        let mut capacity = Vec::new();
        let mut resources = vec![vec![]; m];
        let mut gates = vec![vec![]; m];
        let mut conflicts = vec![vec![]; n];
        // 比例等价的装载行共用一个余量，保留其中最紧者。
        let mut by_terms = BTreeMap::new();
        let mut add = |terms: &[(usize, f64)], cap: f64, sign: f64| -> Option<()> {
            if !cap.is_finite() || cap < 0. || terms.iter().any(|&(_, a)| !a.is_finite()) {
                return None;
            }
            if terms.iter().all(|&(_, a)| a * sign <= 0.) {
                return Some(());
            }
            if terms.iter().any(|&(_, a)| a * sign < 0.) {
                return None;
            }
            let norm = terms[0].1.abs();
            let key: Vec<_> = terms
                .iter()
                .map(|&(k, a)| (k, (a * sign / norm).to_bits()))
                .collect();
            let cap = cap / norm;
            if !cap.is_finite() {
                return None;
            }
            if let Some(&i) = by_terms.get(&key) {
                let old: &mut f64 = &mut capacity[i];
                *old = old.min(cap);
            } else {
                let i = capacity.len();
                by_terms.insert(key, i);
                capacity.push(cap);
                for &(k, a) in terms {
                    resources[k].push((i, a * sign / norm));
                }
            }
            Some(())
        };
        for (j, v) in p
            .vars
            .iter()
            .enumerate()
            .filter(|(_, v)| v.var_type == VarType::Continuous)
        {
            let (off, terms) = &affine[j];
            add(terms, off - v.lb, -1.)?;
            if v.ub.is_finite() {
                add(terms, v.ub - off, 1.)?;
            }
        }
        for (i, r) in p.constraints.iter().enumerate() {
            if use_eq[i] {
                continue;
            }
            if i % 128 == 0 && stop() {
                return None;
            }
            let s = if r.op == RelOp::Ge { -1. } else { 1. };
            if let Some((gate, cap)) = at_one(p, r) {
                let terms: Terms = r
                    .expr
                    .terms
                    .iter()
                    .filter(|&&(j, _)| j != gate)
                    .map(|&(j, a)| {
                        let k = map[j];
                        gates[k].push(gate);
                        (k, s * a)
                    })
                    .collect();
                add(&terms, cap, 1.)?;
            } else if r.expr.terms.iter().any(|&(j, _)| is_binary(&p.vars[j])) {
                if r.op == RelOp::Eq
                    || r.rhs * s <= 0.
                    || !r
                        .expr
                        .terms
                        .iter()
                        .all(|&(j, a)| is_binary(&p.vars[j]) && a == r.rhs)
                {
                    return None;
                }
                work = work.checked_sub(r.expr.terms.len().checked_mul(r.expr.terms.len())?)?;
                for (pos, &(j, _)) in r.expr.terms.iter().enumerate() {
                    if pos % 32 == 0 && stop() {
                        return None;
                    }
                    conflicts[j].extend(r.expr.terms.iter().map(|&(k, _)| k).filter(|&k| k != j));
                }
            } else {
                charge_projection(&r.expr, &affine, &mut work)?;
                let (off, terms) = project(&r.expr, &affine, s, None, &mut scratch, &mut touched);
                add(&terms, s * r.rhs - off, 1.)?;
                if r.op == RelOp::Eq {
                    add(&terms, off - s * r.rhs, -1.)?;
                }
            }
        }
        for list in gates.iter_mut().chain(&mut conflicts) {
            list.sort_unstable();
            list.dedup();
        }
        charge_projection(&p.obj, &affine, &mut work)?;
        let (_, mut terms) = project(&p.obj, &affine, -signed, None, &mut scratch, &mut touched);
        let mut gain = vec![0.; m];
        for (k, a) in terms.drain(..) {
            gain[k] = a;
        }
        for k in 0..m {
            if gates[k].iter().any(|&j| {
                conflicts[j]
                    .iter()
                    .any(|g| gates[k].binary_search(g).is_ok())
            }) {
                gain[k] = 0.;
            }
        }
        // 非负流量使每条装载行也给单列一个上界。显式保留这些上界后，
        // 即使同时把所有列推到上界也不会满载的行，无需在每次邻域扫描时访问。
        let flow_limits: Vec<_> = resources
            .iter()
            .map(|rows| {
                rows.iter()
                    .map(|&(i, c)| capacity[i] / c)
                    .fold(f64::INFINITY, f64::min)
            })
            .collect();
        let mut load = vec![0.; capacity.len()];
        for (k, rows) in resources.iter().enumerate() {
            for &(i, c) in rows {
                load[i] += flow_limits[k] * c;
            }
        }
        let mut compact = Vec::new();
        let mut index = vec![None; capacity.len()];
        for (i, cap) in capacity.into_iter().enumerate() {
            if load[i] > cap {
                index[i] = Some(compact.len());
                compact.push(cap);
            }
        }
        for rows in &mut resources {
            rows.retain_mut(|(i, _)| {
                if let Some(j) = index[*i] {
                    *i = j;
                    true
                } else {
                    false
                }
            });
        }
        Some(Self {
            affine,
            capacity: compact,
            flow_limits,
            resources,
            gates,
            conflicts,
            gain,
            fee,
        })
    }

    fn expand(&self, p: &Problem, x: &[f64]) -> Vec<f64> {
        let mut values: Vec<_> = self
            .affine
            .iter()
            .map(|(b, terms)| b + terms.iter().map(|&(k, a)| x[k] * a).sum::<f64>())
            .collect();
        for (k, &quantity) in x.iter().enumerate() {
            if quantity > 0. {
                for &j in &self.gates[k] {
                    values[j] = 1.;
                }
            }
        }
        // 没有互斥或启用费的开关开足容量，使可行解也符合单调固定后的变量界。
        for (j, v) in p.vars.iter().enumerate() {
            if is_binary(v) && self.fee[j] == 0. && self.conflicts[j].is_empty() {
                values[j] = 1.;
            }
        }
        values
    }
}

// 每次重启使用固定序列，不读取墙钟生成扰动。
struct Noise(u64);
impl Noise {
    fn sample(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / ((1u64 << 53) as f64)
    }
}

#[derive(Clone)]
struct State {
    x: Vec<f64>,
    left: Vec<f64>,
    active: Vec<bool>,
    forbidden: Vec<bool>,
    profit: f64,
}

impl State {
    fn new(a: &Allocation) -> Self {
        Self {
            x: vec![0.; a.gain.len()],
            left: a.capacity.clone(),
            active: vec![false; a.fee.len()],
            forbidden: vec![false; a.fee.len()],
            profit: 0.,
        }
    }

    fn put(&mut self, a: &Allocation, k: usize, amount: f64) {
        self.x[k] += amount;
        self.profit += amount * a.gain[k];
        for &(i, c) in &a.resources[k] {
            self.left[i] -= c * amount;
        }
        for &j in &a.gates[k] {
            if !self.active[j] {
                self.active[j] = true;
                self.profit -= a.fee[j];
                for &other in &a.conflicts[j] {
                    self.forbidden[other] = true;
                }
            }
        }
    }

    fn fill(
        &mut self,
        a: &Allocation,
        rng: &mut Noise,
        iteration: usize,
        work: &mut usize,
        stop: &mut dyn FnMut() -> bool,
    ) {
        let mut candidates = Vec::new();
        for (k, &gain) in a.gain.iter().enumerate() {
            if gain > 0. {
                // 部分重启按单位收益排序，部分按一次能实现的总收益排序。
                let noise = if iteration == 0 {
                    1.
                } else {
                    0.45 + rng.sample()
                };
                candidates.push((k, noise));
            }
        }
        let mut steps = 0;
        while !candidates.is_empty() {
            if *work == 0 || (steps % 8 == 0 && stop()) {
                return;
            }
            steps += 1;
            let mut best: Option<(usize, f64, f64)> = None;
            candidates.retain(|&(k, noise)| {
                *work = work.saturating_sub(1 + a.gates[k].len() * 2 + a.resources[k].len());
                if *work == 0 {
                    return false;
                }
                if a.gates[k].iter().any(|&j| self.forbidden[j]) {
                    return false;
                }
                let quantity = a.resources[k]
                    .iter()
                    .map(|&(i, c)| self.left[i] / c)
                    .fold(a.flow_limits[k] - self.x[k], f64::min);
                if quantity <= 1e-8 || !quantity.is_finite() {
                    return false;
                }
                let fee: f64 = a.gates[k]
                    .iter()
                    .filter(|&&j| !self.active[j])
                    .map(|&j| a.fee[j])
                    .sum();
                // 共享启用费可由多列摊薄；改变构造时估算的费用，评价解始终照原价。
                let estimated_fee = fee * [1.0, 1.0, 0.7, 1.4, 0.4][iteration % 5];
                let value = quantity * a.gain[k] - estimated_fee;
                if value > 0. {
                    let score = match iteration % 3 {
                        0 => value,
                        1 => value / quantity.sqrt(),
                        _ => value / quantity,
                    };
                    let score = score * noise;
                    if best.is_none_or(|(_, _, old)| score > old) {
                        best = Some((k, quantity, score));
                    }
                }
                true
            });
            let Some((k, amount, _)) = best else { return };
            self.put(a, k, amount);
        }
    }
}

pub(crate) fn primal(
    p: &Problem,
    stop: &mut dyn FnMut() -> bool,
    report: &mut dyn FnMut(&[f64]),
) -> Option<Vec<f64>> {
    let a = Allocation::extract(p, stop)?;
    let mut best = State::new(&a);
    if !a.gain.iter().any(|&g| g > 0.) {
        return Some(a.expand(p, &best.x));
    }
    let mut current = best.clone();
    let mut last_improvement = 0;
    let mut rng = Noise(0x9216_d5d9_8979_fb1b);
    let n = a.gain.len();
    // 初次按纯收益构造，随后重新分配部分已启用开关上的流量。
    // 除截止时间外，总量也受列数与非零元数约束。
    let mut work = 750_000_000usize;
    // 小模型很快耗尽有意义的邻域；大模型在扫描候选列时按访问的非零元扣费。
    let state_width = n + a.fee.len() + a.capacity.len();
    let trials = n
        .saturating_mul(64)
        .min(32_768)
        .min((100_000_000 / state_width).max(1));
    for iteration in 0..trials {
        if work == 0 || stop() {
            break;
        }
        let mut next = State::new(&a);
        if iteration % 20 != 0 {
            let incumbent = if iteration % 4 == 0 { &best } else { &current };
            let target = (rng.sample() * n as f64) as usize;
            let mut remove_gate = vec![false; a.fee.len()];
            for &g in &a.gates[target] {
                remove_gate[g] = true;
                for &j in &a.conflicts[g] {
                    remove_gate[j] = true;
                }
            }
            let mut remove_row = vec![false; a.capacity.len()];
            let proportion = rng.sample() * 0.8;
            for &(i, _) in &a.resources[target] {
                if incumbent.left[i] < 1e-8 && rng.sample() < proportion {
                    remove_row[i] = true;
                }
            }
            let rate = 0.01 + rng.sample() * 0.25;
            for (k, &x) in incumbent.x.iter().enumerate() {
                if x > 0.
                    && !a.gates[k].iter().any(|&j| remove_gate[j])
                    && !a.resources[k].iter().any(|&(i, _)| remove_row[i])
                    && rng.sample() >= rate
                {
                    next.put(&a, k, x);
                }
            }
        }
        next.fill(&a, &mut rng, iteration, &mut work, stop);
        if next.profit > best.profit {
            last_improvement = iteration;
            best = next.clone();
            report(&a.expand(p, &best.x));
        }
        if next.profit >= current.profit
            || (iteration - last_improvement > 50 && iteration % 51 == 0)
        {
            current = next;
        }
        if iteration % 250 == 0 {
            current = best.clone();
        }
    }
    Some(a.expand(p, &best.x))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instance::parse_instance;

    #[test]
    fn shared_fee_and_sign_reversal() {
        // 每票单独都不能回本，两票共用一次启用费则可以。
        let source = "name shared\nsense max\nvars\nx continuous\nz continuous\ny integer\nobj\n2 x 2 z -3 y\ncons\n-5 x 5 y >= 0\n3 z -3 y <= 0\nbounds\ny <= 1\nend";
        let p = parse_instance(source).unwrap();
        let mut reported = Vec::new();
        let x = primal(&p, &mut || false, &mut |x| reported.push(p.obj.eval(x))).unwrap();
        assert_eq!(x, [1., 1., 1.]);
        assert_eq!(reported.last(), Some(&1.));
    }

    #[test]
    fn reconstructs_equalities_and_respects_mutex() {
        let source = "name inventory\nsense min\nvars\nstock continuous\nx continuous\nz continuous\nleft binary\nright binary\nfree binary\nobj\n3 stock 2 left 1 right\ncons\n-3 stock -3 x -3 z = -12\n1 x -4 left <= 0\n1 z -4 right <= 0\n1 x 1 z -4 free <= 0\n-2 left -2 right >= -2\nend";
        let p = parse_instance(source).unwrap();
        let x = primal(&p, &mut || false, &mut |_| {}).unwrap();
        assert_eq!(x, [0., 0., 4., 0., 1., 1.]);
        assert_eq!(p.obj.eval(&x), 1.);
    }

    #[test]
    fn reverse_order_chain_is_processed_once() {
        let mut p = Problem::new("reverse", Sense::Minimize);
        let x = p.add_var("flow", VarType::Continuous).unwrap();
        let gate = p.add_var("gate", VarType::Binary).unwrap();
        let chain: Vec<_> = (0..2200)
            .map(|j| {
                p.add_var(&format!("stock{j}"), VarType::Continuous)
                    .unwrap()
            })
            .collect();
        for pair in chain.windows(2).rev() {
            p.add_constraint(None, vec![(pair[0], 1.), (pair[1], -1.)], RelOp::Eq, 0.)
                .unwrap();
        }
        p.add_constraint(None, vec![(chain[0], 1.), (x, 1.)], RelOp::Eq, 4.)
            .unwrap();
        p.add_constraint(None, vec![(x, 1.), (gate, -4.)], RelOp::Le, 0.)
            .unwrap();
        p.set_objective(vec![(*chain.last().unwrap(), 2.), (gate, 1.)], 0.)
            .unwrap();
        let mut polls = 0;
        let a = Allocation::extract(&p, &mut || {
            polls += 1;
            polls > 3000
        })
        .unwrap();
        assert_eq!(a.gain, [2.]);
    }

    #[test]
    fn mandatory_flow_and_interrupted_extraction_fall_back() {
        let source = "name mandatory\nsense max\nvars\nx continuous\ny binary\nobj\n2 x -1 y\ncons\n1 x -4 y <= 0\n1 x >= 1\nend";
        let p = parse_instance(source).unwrap();
        assert!(primal(&p, &mut || false, &mut |_| {}).is_none());
        assert!(primal(&p, &mut || true, &mut |_| {}).is_none());
    }
}
