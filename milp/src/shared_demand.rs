//! 局部二元变量只参与一条需求，其余约束仅决定它何时可用。固定共享变量后，
//! 各需求可独立求最小成本。按影响同一批需求的共享变量分组，删去被支配的赋值，
//! 再按最小联结表依次消去各组。回溯恢复全部原变量，并对原模型复核。
//! 各表在整数上精确计算。先估算表宽及总工作量，规模不合适时交回通用求解器。

use crate::error::{MilpError, Result};
use crate::model::{LinearExpr, Problem, RelOp, Sense, SolveStatus, Tolerances, VarType};
use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

const INF: u32 = u32::MAX;
const TABLE_CELLS: usize = 1 << 20;
const TOTAL_CELLS: usize = 1 << 25;
const GROUP_BITS: usize = 16;
const LOCAL_BITS: usize = 8;

fn expired(end: Option<Instant>) -> bool {
    end.is_some_and(|t| Instant::now() >= t)
}

fn integer_norm(expr: &LinearExpr, rhs: f64) -> Option<f64> {
    let mut norm = 0.0;
    for &(_, a) in &expr.terms {
        if a != a.trunc() {
            return None;
        }
        norm += a.abs();
    }
    (rhs == rhs.trunc() && norm + rhs.abs() < (1 << 28) as f64).then_some(norm)
}

#[derive(Default, Clone, Copy)]
struct Effect {
    weight: i32,
    blocked: [u8; 2],
}

struct Demand {
    rhs: i32,
    blocked: u8,
    locals: Vec<usize>,
    // 按原列下标顺序的位掩码枚举，记录活动度和成本，供回溯使用。
    options: Vec<(i32, u32)>,
    minimum: i32,
    span: usize,
    lookup: Vec<u32>,
}

impl Demand {
    fn new(rhs: i32, terms: &[(usize, i32)], costs: &[u32], capacity: &mut usize) -> Option<Self> {
        if terms.len() > LOCAL_BITS {
            return None;
        }
        let minimum: i32 = terms.iter().map(|&(_, a)| a.min(0)).sum();
        let maximum: i32 = terms.iter().map(|&(_, a)| a.max(0)).sum();
        let span = (maximum - minimum + 1) as usize;
        let count = 1usize << terms.len();
        let cells = count.checked_mul(span)?;
        *capacity = capacity.checked_sub(cells)?;
        let mut options = vec![(0, 0); count];
        for mask in 1..count {
            let (j, a) = terms[mask.trailing_zeros() as usize];
            let (w, c) = options[mask & (mask - 1)];
            options[mask] = (w + a, c + costs[j]);
        }
        let mut lookup = vec![INF; cells];
        for (blocked, row) in lookup.chunks_exact_mut(span).enumerate() {
            let free = count - 1 - blocked;
            let mut mask = free;
            loop {
                let (w, c) = options[mask];
                let best = &mut row[(w - minimum) as usize];
                *best = (*best).min(c);
                if mask == 0 {
                    break;
                }
                mask = (mask - 1) & free;
            }
            for i in 1..span {
                row[i] = row[i].min(row[i - 1]);
            }
        }
        Some(Self {
            rhs,
            blocked: 0,
            locals: terms.iter().map(|&(j, _)| j).collect(),
            options,
            minimum,
            span,
            lookup,
        })
    }

    #[inline]
    fn cost(&self, rhs: i32, blocked: u8) -> u32 {
        if rhs < self.minimum {
            INF
        } else {
            self.lookup
                [blocked as usize * self.span + ((rhs - self.minimum) as usize).min(self.span - 1)]
        }
    }

    fn restore(&self, rhs: i32, blocked: u8, x: &mut [f64]) -> Option<u32> {
        let cost = self.cost(rhs, blocked);
        let mask = self
            .options
            .iter()
            .enumerate()
            .position(|(mask, &(w, c))| mask & blocked as usize == 0 && w <= rhs && c == cost)?;
        for (k, &j) in self.locals.iter().enumerate() {
            x[j] = ((mask >> k) & 1) as f64;
        }
        Some(cost)
    }
}

struct Structure {
    costs: Vec<u32>,
    demands: Vec<Demand>,
    // 每条需求的直接系数，以及某共享变量取 0 或 1 时禁止的局部选择。
    effects: Vec<Vec<Effect>>,
    shared: Vec<usize>,
    selected: Vec<usize>,
    fixed_cost: u32,
}

/// 对 x+y≥1，若从 (1,0) 换为 (0,1) 在其他各行的活动度和成本均不增，
/// 且从 (1,1) 删掉 x 也如此，则一定存在 x=0、y=1 的最优解。检查包含门控在内的原行。
fn dominated_covers(p: &Problem, gates: &[Option<usize>], costs: &[u32]) -> Vec<Option<bool>> {
    let mut fixed = vec![None; costs.len()];
    for (i, r) in p.constraints.iter().enumerate() {
        let &[(j, a), (k, b)] = r.expr.terms.as_slice() else {
            continue;
        };
        if a >= 0.0
            || a != b
            || a != r.rhs
            || gates[j].is_some()
            || gates[k].is_some()
            || fixed[j].is_some()
            || fixed[k].is_some()
        {
            continue;
        }
        for (keep, drop) in [(j, k), (k, j)] {
            if costs[keep] > costs[drop] {
                continue;
            }
            let valid = p.constraints.iter().enumerate().all(|(h, row)| {
                if h == i {
                    return true;
                }
                let terms = &row.expr.terms;
                let coeff = |j| {
                    terms
                        .binary_search_by_key(&j, |t| t.0)
                        .map_or(0.0, |v| terms[v].1)
                };
                let old = coeff(drop);
                old >= 0.0 && coeff(keep) <= old
            });
            if valid {
                fixed[keep] = Some(true);
                fixed[drop] = Some(false);
                break;
            }
        }
    }
    fixed
}

impl Structure {
    fn detect(p: &Problem, tol: Tolerances, end: Option<Instant>) -> Option<Self> {
        let n = p.vars.len();
        if p.sense != Sense::Minimize
            || n > 512
            || p.constraints.len() > 512
            || p.vars
                .iter()
                .any(|v| v.var_type != VarType::Binary || v.lb != 0.0 || v.ub != 1.0)
        {
            return None;
        }
        integer_norm(&p.obj, p.obj_offset)?;
        let mut costs = vec![0; n];
        for &(j, c) in &p.obj.terms {
            if c < 0.0 {
                return None;
            }
            costs[j] = c as u32;
        }
        let mut local_gate = vec![None; n];
        let mut gate_rows = BTreeSet::new();
        for (i, r) in p.constraints.iter().enumerate() {
            let norm = integer_norm(&r.expr, r.rhs)?;
            if r.op != RelOp::Le
                || !(10.0 * tol.feasibility + tol.integrality * norm < 1.0
                    && crate::solver::FINAL_VIOLATION_REL * (1.0 + r.rhs.abs() + norm) < 1.0)
            {
                return None;
            }
            if r.rhs != 0.0 || r.expr.terms.len() < 3 {
                continue;
            }
            let mut negatives = r.expr.terms.iter().filter(|&&(_, a)| a < 0.0);
            let Some(&(g, a)) = negatives.next() else {
                continue;
            };
            // g=0 时所有正项必须为零，g=1 时这行必须容得下全部正项。
            if negatives.next().is_some() || norm + a > -a {
                continue;
            }
            gate_rows.insert(i);
            for &(j, w) in &r.expr.terms {
                if w > 0.0 && local_gate[j].replace(g).is_some() {
                    return None;
                }
            }
        }
        if gate_rows.is_empty()
            || local_gate
                .iter()
                .flatten()
                .any(|&g| local_gate[g].is_some())
        {
            return None;
        }
        let fixed = dominated_covers(p, &local_gate, &costs);
        let mut owned = vec![false; n];
        let mut conflicts = vec![Vec::new(); n];
        let mut locations = vec![None; n];
        let mut demands = Vec::new();
        let mut effects = vec![Vec::new(); n];
        let mut lut_capacity = TABLE_CELLS;
        for (i, r) in p.constraints.iter().enumerate() {
            if expired(end) {
                return None;
            }
            if gate_rows.contains(&i) {
                continue;
            }
            if let &[(j, a), (k, b)] = r.expr.terms.as_slice() {
                if a > 0.0
                    && a == b
                    && r.rhs == a
                    && local_gate[j].is_some() != local_gate[k].is_some()
                {
                    let (local, shared) = if local_gate[j].is_some() {
                        (j, k)
                    } else {
                        (k, j)
                    };
                    conflicts[local].push(shared);
                    continue;
                }
            }
            let mut terms = Vec::new();
            let mut rhs = r.rhs as i32;
            let mut shared_max = 0;
            for &(j, a) in &r.expr.terms {
                if local_gate[j].is_some() {
                    if std::mem::replace(&mut owned[j], true) {
                        return None;
                    }
                    terms.push((j, a as i32));
                } else {
                    match fixed[j] {
                        Some(true) => rhs -= a as i32,
                        Some(false) => (),
                        None => shared_max += (a as i32).max(0),
                    }
                }
            }
            // 任何共享赋值下，局部列取零都可行且成本为零，无须建立这张表。
            if shared_max <= rhs {
                continue;
            }
            let f = demands.len();
            for column in &mut effects {
                column.push(Effect::default());
            }
            for &(j, a) in &r.expr.terms {
                if local_gate[j].is_none() && fixed[j].is_none() {
                    effects[j][f].weight = a as i32;
                }
            }
            for (k, &(j, _)) in terms.iter().enumerate() {
                locations[j] = Some((f, k));
            }
            demands.push(Demand::new(rhs, &terms, &costs, &mut lut_capacity)?);
        }
        let mut shared = Vec::new();
        for (j, gate) in local_gate.into_iter().enumerate() {
            if let Some(g) = gate {
                let Some((f, k)) = locations[j] else { continue };
                effects[g][f].blocked[0] |= 1 << k;
                for &other in &conflicts[j] {
                    effects[other][f].blocked[1] |= 1 << k;
                }
            } else if fixed[j].is_none() {
                shared.push(j);
            }
        }
        let mut selected = Vec::new();
        let mut fixed_cost = 0;
        for (j, value) in fixed.into_iter().enumerate() {
            if let Some(value) = value {
                if value {
                    selected.push(j);
                    fixed_cost += costs[j];
                }
                for (f, d) in demands.iter_mut().enumerate() {
                    d.blocked |= effects[j][f].blocked[value as usize];
                }
            }
        }
        Some(Self {
            costs,
            demands,
            effects,
            shared,
            selected,
            fixed_cost,
        })
    }

    fn restore(&self, mut x: Vec<f64>) -> Option<(u32, Vec<f64>)> {
        let mut cost: u32 = self
            .shared
            .iter()
            .filter(|&&j| x[j] == 1.0)
            .map(|&j| self.costs[j])
            .sum();
        for &j in &self.selected {
            x[j] = 1.0;
        }
        cost += self.fixed_cost;
        for (f, d) in self.demands.iter().enumerate() {
            let mut rhs = d.rhs;
            let mut blocked = d.blocked;
            for &j in &self.shared {
                let e = self.effects[j][f];
                let value = x[j] as usize;
                if value == 1 {
                    rhs -= e.weight;
                }
                blocked |= e.blocked[value];
            }
            cost += d.restore(rhs, blocked, &mut x)?;
        }
        Some((cost, x))
    }
}

struct Group {
    cols: Vec<usize>,
    // 选项的原变量掩码与成本。
    states: Vec<(u16, u32)>,
    // 每个受影响需求中，各选项的直接活动度和局部禁用位。
    effects: Vec<Vec<(i32, u8)>>,
}

fn groups(model: &Structure, end: Option<Instant>) -> Option<Vec<Group>> {
    let mut buckets = BTreeMap::<(Vec<usize>, Option<usize>), Vec<usize>>::new();
    for &j in &model.shared {
        let col = &model.effects[j];
        let scope: Vec<_> = col
            .iter()
            .enumerate()
            .filter(|(_, e)| e.weight != 0 || e.blocked != [0; 2])
            .map(|(f, _)| f)
            .collect();
        let Some(&first) = scope.first() else {
            continue;
        };
        let varying = scope.iter().any(|&f| col[f].weight != col[first].weight);
        buckets
            .entry((scope, varying.then_some(j)))
            .or_default()
            .push(j);
    }
    let mut result = Vec::new();
    let mut work = TOTAL_CELLS;
    for ((scope, varying), cols) in buckets {
        if expired(end) || cols.len() > GROUP_BITS {
            return None;
        }
        let controlled: Vec<_> = cols
            .iter()
            .enumerate()
            .filter(|&(_, &j)| scope.iter().any(|&f| model.effects[j][f].blocked != [0; 2]))
            .map(|(k, _)| k)
            .collect();
        // 当控制位过多，判定两个禁用集的包含关系本身可能超过求解预算。
        if controlled.len() > 6 {
            return None;
        }
        let classes = 1usize << controlled.len();
        let count = 1usize << cols.len();
        work = work.checked_sub(classes * count + classes * classes * scope.len())?;
        let mut blocked = vec![vec![0u8; scope.len()]; classes];
        let mut ctrl_bit = vec![0usize; cols.len()];
        for (b, &k) in controlled.iter().enumerate() {
            ctrl_bit[k] = 1 << b;
            for (class, masks) in blocked.iter_mut().enumerate() {
                for (u, &f) in scope.iter().enumerate() {
                    masks[u] |= model.effects[cols[k]][f].blocked[(class >> b) & 1];
                }
            }
        }
        let subset: Vec<Vec<usize>> = blocked
            .iter()
            .map(|row| {
                blocked
                    .iter()
                    .enumerate()
                    .filter(|(_, other)| other.iter().zip(row).all(|(&a, &b)| a & !b == 0))
                    .map(|(k, _)| k)
                    .collect()
            })
            .collect();
        let bits: Vec<u32> = blocked
            .iter()
            .map(|row| row.iter().map(|b| b.count_ones()).sum())
            .collect();
        // 排序保证已保留的选项成本不更高。仅当每一行的活动度和禁用集合均不增才可替换。
        let mut options = Vec::with_capacity(count);
        options.push((0u32, 0i32, 0usize, 0u16));
        for mask in 1..count {
            let k = mask.trailing_zeros() as usize;
            let (c, w, class, _) = options[mask & (mask - 1)];
            options.push((
                c + model.costs[cols[k]],
                w + model.effects[cols[k]][scope[0]].weight,
                class | ctrl_bit[k],
                mask as u16,
            ));
        }
        options.sort_unstable_by_key(|&(c, w, class, mask)| (c, w, bits[class], mask));
        let mut least = vec![i32::MAX; classes];
        let mut states = Vec::new();
        let mut effects = vec![Vec::new(); model.demands.len()];
        for (c, w, class, mask) in options {
            if varying.is_none() && subset[class].iter().any(|&other| least[other] <= w) {
                continue;
            }
            if states.len() >= 1024 {
                return None;
            }
            least[class] = least[class].min(w);
            states.push((mask, c));
            for (k, &f) in scope.iter().enumerate() {
                let weight = varying.map_or(w, |j| i32::from(mask) * model.effects[j][f].weight);
                effects[f].push((weight, blocked[class][k]));
            }
        }
        result.push(Group {
            cols,
            states,
            effects,
        });
    }
    // 较大的维放在内层，连续查询成本表。相同大小按首列下标破键。
    result.sort_by_key(|g| (g.states.len(), g.cols[0]));
    Some(result)
}

type Scope = Vec<usize>;

fn size(s: &[usize], groups: &[Group]) -> Option<usize> {
    s.iter().try_fold(1usize, |n, &j| {
        n.checked_mul(groups[j].states.len())
            .filter(|&v| v <= TABLE_CELLS)
    })
}

fn index(s: &[usize], assignment: &[usize], groups: &[Group]) -> usize {
    s.iter()
        .fold(0, |i, &j| i * groups[j].states.len() + assignment[j])
}

struct Step {
    group: usize,
    union: Scope,
    inputs: Vec<Scope>,
    output: Scope,
}

struct Plan {
    rows: Vec<Scope>,
    steps: Vec<Step>,
}

impl Plan {
    fn new(groups: &[Group], demands: &[Demand]) -> Option<Self> {
        if groups.len() > 32 {
            return None;
        }
        let mut budget = TOTAL_CELLS;
        let mut rows = Vec::new();
        let mut scopes = BTreeSet::new();
        for f in 0..demands.len() {
            let s: Scope = (0..groups.len())
                .filter(|&j| !groups[j].effects[f].is_empty())
                .collect();
            budget = budget.checked_sub(size(&s, groups)?)?;
            scopes.insert(s.clone());
            rows.push(s);
        }
        scopes.extend((0..groups.len()).map(|j| vec![j]));
        let mut steps = Vec::new();
        for _ in 0..groups.len() {
            let mut best = None;
            for (j, g) in groups.iter().enumerate() {
                let union: BTreeSet<_> = scopes
                    .iter()
                    .filter(|s| s.contains(&j))
                    .flatten()
                    .copied()
                    .collect();
                if union.is_empty() {
                    continue;
                }
                let union: Scope = union.into_iter().collect();
                if let Some(cells) = size(&union, groups) {
                    let key = (cells, g.cols[0]);
                    if best.as_ref().is_none_or(|(old, _, _)| key < *old) {
                        best = Some((key, j, union));
                    }
                }
            }
            let ((cells, _), group, union) = best?;
            budget = budget.checked_sub(cells)?;
            let inputs: Vec<_> = scopes
                .iter()
                .filter(|s| s.contains(&group))
                .cloned()
                .collect();
            for s in &inputs {
                scopes.remove(s);
            }
            let output: Scope = union.iter().copied().filter(|&j| j != group).collect();
            scopes.insert(output.clone());
            steps.push(Step {
                group,
                union,
                inputs,
                output,
            });
        }
        Some(Self { rows, steps })
    }
}

fn fill(
    row: &Demand,
    axes: &[&[(i32, u8)]],
    rhs: i32,
    blocked: u8,
    table: &mut [u32],
    end: Option<Instant>,
) -> bool {
    if table.len() >= 4096 && expired(end) {
        return false;
    }
    match axes.split_first() {
        None => table[0] = table[0].saturating_add(row.cost(rhs, blocked)),
        Some((last, [])) => {
            for (v, &(w, b)) in table.iter_mut().zip(*last) {
                *v = v.saturating_add(row.cost(rhs - w, blocked | b));
            }
        }
        Some((axis, tail)) => {
            let stride = table.len() / axis.len();
            for (chunk, &(w, b)) in table.chunks_exact_mut(stride).zip(*axis) {
                if !fill(row, tail, rhs - w, blocked | b, chunk, end) {
                    return false;
                }
            }
        }
    }
    true
}

struct Trace {
    group: usize,
    output: Scope,
    choice: Vec<u16>,
}

fn eliminate(
    step: Step,
    tables: &mut BTreeMap<Scope, Vec<u32>>,
    groups: &[Group],
    end: Option<Instant>,
) -> Option<Trace> {
    let inputs: Vec<_> = step
        .inputs
        .into_iter()
        .map(|s| {
            let values = tables.remove(&s).expect("消元源表必须存在");
            (s, values)
        })
        .collect();
    let cells = size(&step.output, groups).expect("消元宽度已验证");
    let mut minimum = vec![INF; cells];
    let mut choice = vec![0u16; cells];
    let mut x = vec![0; groups.len()];
    // 最内层是最长的选项表。固定外层取值后，各张输入表在这一维要么恒定，
    // 要么是一段连续切片；一次定位整段即可，避免每个格反复计算多维下标。
    let (&inner, outer) = step.union.split_last().expect("消元至少涉及本组");
    let stride = groups[inner].states.len();
    let mut totals = vec![0u32; stride];
    let outer_size = size(outer, groups).expect("联结宽度已验证");
    // 每个外层块最多 1024 项。
    for i in 0..outer_size {
        if i % (4096 / stride).max(1) == 0 && expired(end) {
            return None;
        }
        let mut constant = 0u32;
        totals.fill(0);
        for (scope, values) in &inputs {
            let base = index(scope, &x, groups);
            if scope.last() == Some(&inner) {
                for (total, &value) in totals.iter_mut().zip(&values[base..base + stride]) {
                    *total = total.saturating_add(value);
                }
            } else {
                constant = constant.saturating_add(values[base]);
            }
        }
        let cell = index(&step.output, &x, groups);
        if inner == step.group {
            for (i, &value) in totals.iter().enumerate() {
                let cost = value.saturating_add(constant);
                if cost < minimum[cell] {
                    minimum[cell] = cost;
                    choice[cell] = i as u16;
                }
            }
        } else {
            let state = x[step.group] as u16;
            for ((best, arg), &value) in minimum[cell..cell + stride]
                .iter_mut()
                .zip(&mut choice[cell..cell + stride])
                .zip(&totals)
            {
                let cost = value.saturating_add(constant);
                if cost < *best {
                    *best = cost;
                    *arg = state;
                }
            }
        }
        for &j in outer.iter().rev() {
            x[j] += 1;
            if x[j] < groups[j].states.len() {
                break;
            }
            x[j] = 0;
        }
    }
    if let Some(t) = tables.get_mut(&step.output) {
        for (a, b) in t.iter_mut().zip(minimum) {
            *a = a.saturating_add(b);
        }
    } else {
        tables.insert(step.output.clone(), minimum);
    }
    Some(Trace {
        group: step.group,
        output: step.output,
        choice,
    })
}

fn verified(p: &Problem, cost: u32, x: &[f64]) -> bool {
    p.obj.eval(x) == f64::from(cost) && p.constraints.iter().all(|r| r.expr.eval(x) <= r.rhs)
}

/// 单个根节点内的整数动态规划。预算预检不通过时，调用方可无副作用地继续原搜索。
pub(crate) fn solve(
    p: &Problem,
    tol: Tolerances,
    end: Option<Instant>,
    accept: &mut dyn FnMut(&[f64]),
) -> Result<Option<(SolveStatus, Option<f64>)>> {
    let limited = Ok(Some((SolveStatus::LimitReached, None)));
    let Some(model) = Structure::detect(p, tol, end) else {
        return Ok(None);
    };
    let Some(groups) = groups(&model, end) else {
        return if expired(end) { limited } else { Ok(None) };
    };
    let Some(plan) = Plan::new(&groups, &model.demands) else {
        return Ok(None);
    };
    if expired(end) {
        return limited;
    }
    let mut seed = vec![0.0; p.vars.len()];
    for &j in &model.shared {
        seed[j] = 1.0;
    }
    if let Some((mut cost, mut x)) = model.restore(seed) {
        if !verified(p, cost, &x) {
            return Err(MilpError::Numerical("局部补位的初解未通过原行核验".into()));
        }
        accept(&x);
        // 默认先启用共享列以取得首解，再优先关掉昂贵的列；每次由各需求重新选局部补位。
        let mut order = model.shared.clone();
        order.sort_unstable_by_key(|&j| (std::cmp::Reverse(model.costs[j]), j));
        for j in order {
            if expired(end) {
                return limited;
            }
            let mut trial = x.clone();
            trial[j] = 0.0;
            if let Some((next, candidate)) = model.restore(trial) {
                if next < cost {
                    if !verified(p, next, &candidate) {
                        return Err(MilpError::Numerical(
                            "局部补位的改进解未通过原行核验".into(),
                        ));
                    }
                    accept(&candidate);
                    cost = next;
                    x = candidate;
                }
            }
        }
    }
    let mut tables = BTreeMap::from([(Vec::new(), vec![model.fixed_cost])]);
    for (f, scope) in plan.rows.into_iter().enumerate() {
        if expired(end) {
            return limited;
        }
        let axes: Vec<_> = scope
            .iter()
            .map(|&j| groups[j].effects[f].as_slice())
            .collect();
        let d = &model.demands[f];
        let cells = size(&scope, &groups).expect("因子宽度已验证");
        let table = tables.entry(scope).or_insert_with(|| vec![0; cells]);
        if !fill(d, &axes, d.rhs, d.blocked, table, end) {
            return limited;
        }
    }
    for (j, g) in groups.iter().enumerate() {
        let t = tables
            .entry(vec![j])
            .or_insert_with(|| vec![0; g.states.len()]);
        for (v, &(_, c)) in t.iter_mut().zip(&g.states) {
            *v = v.saturating_add(c);
        }
    }
    let mut history = Vec::new();
    for step in plan.steps {
        let Some(trace) = eliminate(step, &mut tables, &groups, end) else {
            return limited;
        };
        history.push(trace);
    }
    let best = tables.get(&vec![]).expect("消元完成后必须有常数因子")[0];
    if expired(end) {
        return limited;
    }
    if best == INF {
        return Ok(Some((SolveStatus::Infeasible, None)));
    }
    let mut choices = vec![0; groups.len()];
    for trace in history.into_iter().rev() {
        choices[trace.group] = trace.choice[index(&trace.output, &choices, &groups)] as usize;
    }
    let mut x = vec![0.0; p.vars.len()];
    for (g, &i) in groups.iter().zip(&choices) {
        let mask = g.states[i].0;
        for (k, &j) in g.cols.iter().enumerate() {
            x[j] = f64::from((mask >> k) & 1);
        }
    }
    let (cost, x) = model
        .restore(x)
        .ok_or_else(|| MilpError::Numerical("需求消元未能还原局部变量".into()))?;
    if cost != best || !verified(p, cost, &x) {
        return Err(MilpError::Numerical("需求消元的界与原模型解不一致".into()));
    }
    accept(&x);
    if expired(end) {
        return limited;
    }
    Ok(Some((
        SolveStatus::Optimal,
        Some(f64::from(best) + p.obj_offset),
    )))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(p: &mut Problem, terms: &[(usize, i32)], rhs: i32) {
        p.add_constraint(
            None,
            terms.iter().map(|&(j, a)| (j, f64::from(a))).collect(),
            RelOp::Le,
            f64::from(rhs),
        )
        .unwrap();
    }

    fn instance(seed: u32) -> Problem {
        let mut p = Problem::new("", Sense::Minimize);
        for j in 0..9 {
            p.add_var(&format!("v{j}"), VarType::Binary).unwrap();
        }
        let mut state = seed;
        let mut next = |n: u32| {
            state = state.wrapping_mul(1664525).wrapping_add(1013904223);
            (state >> 12) % n
        };
        p.set_objective((0..9).map(|j| (j, f64::from(next(8)))).collect(), -7.0)
            .unwrap();
        row(&mut p, &[(0, 2), (1, 3), (4, -5)], 0);
        row(&mut p, &[(2, 1), (3, 1), (5, -3)], 0);
        row(&mut p, &[(1, 2), (6, 2)], 2);
        row(&mut p, &[(3, 1), (4, 1)], 1);
        for (i, locals) in [[0, 2], [1, 3]].iter().enumerate() {
            row(
                &mut p,
                &[
                    (locals[0], -1 - next(4) as i32),
                    (locals[1], next(5) as i32 - 3),
                    (6, next(7) as i32 - 3),
                    (7, next(7) as i32 - 3),
                    (8, next(7) as i32 - 3),
                ],
                -(i as i32) - 1,
            );
        }
        row(&mut p, &[(7, -1), (8, -1)], -1);
        p
    }

    #[test]
    fn restores_optimum_of_gated_demands_with_shared_only_rows() {
        let mut feasible = 0;
        let mut infeasible = 0;
        for seed in 0..80 {
            let p = instance(seed);
            let brute = (0..1u32 << p.vars.len())
                .filter_map(|mask| {
                    let x: Vec<_> = (0..p.vars.len())
                        .map(|j| f64::from((mask >> j) & 1))
                        .collect();
                    p.constraints
                        .iter()
                        .all(|r| r.expr.eval(&x) <= r.rhs)
                        .then(|| (p.obj.eval(&x) + p.obj_offset) as i64)
                })
                .min();
            let mut points = Vec::new();
            let result = solve(&p, Tolerances::default(), None, &mut |x| {
                points.push(x.to_vec())
            })
            .unwrap()
            .unwrap();
            for x in &points {
                assert!(x.iter().all(|&v| v == 0.0 || v == 1.0));
                assert!(p.constraints.iter().all(|r| r.expr.eval(x) <= r.rhs));
            }
            if let Some(best) = brute {
                feasible += 1;
                assert_eq!(result, (SolveStatus::Optimal, Some(best as f64)));
                assert_eq!(
                    p.obj.eval(points.last().unwrap()) + p.obj_offset,
                    best as f64
                );
            } else {
                infeasible += 1;
                assert_eq!(result, (SolveStatus::Infeasible, None));
                assert!(points.is_empty());
            }
        }
        assert!(feasible > 0 && infeasible > 0);
    }

    #[test]
    fn refuses_insufficient_gate_capacity_and_overlapping_demands() {
        let mut p = instance(0);
        p.constraints[0].expr.add_term(4, 1.0);
        row(&mut p, &[(0, -1), (1, -1), (2, -1)], -1);
        assert_eq!(
            solve(&p, Tolerances::default(), None, &mut |_| panic!(
                "拒绝前产生了解"
            ))
            .unwrap(),
            None
        );
    }

    #[test]
    fn all_zero_local_choices_can_remove_every_factor() {
        let mut p = Problem::new("", Sense::Minimize);
        for j in 0..5 {
            p.add_var(&format!("v{j}"), VarType::Binary).unwrap();
        }
        p.set_objective(vec![(3, 2.0), (4, 1.0)], 3.0).unwrap();
        row(&mut p, &[(0, 1), (1, 1), (2, -2)], 0);
        row(&mut p, &[(0, 9), (1, -3), (3, 4)], 4);
        row(&mut p, &[(3, -1), (4, -1)], -1);
        let mut x = Vec::new();
        assert_eq!(
            solve(&p, Tolerances::default(), None, &mut |v| x = v.to_vec()).unwrap(),
            Some((SolveStatus::Optimal, Some(4.0)))
        );
        assert_eq!(x, vec![0.0, 0.0, 0.0, 0.0, 1.0]);
    }
}
