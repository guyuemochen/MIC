//! 固定费用流量模型的下界。需求余额沿时间累加，付费开关关闭时相应流量必为零。
//! 将每次开通放宽为立即清空本组全部余额，按最近一次开通的时刻做动态规划。
//! 两组只有在同一付费开关下的流量两两互斥时才能分别收取费用。供给不足时，
//! 把同一供给链下各余额共有的一部分单位费用提出来，对累计缺口单独计价。

use crate::capacity_gate::is_binary;
use crate::model::{Constraint, Problem, RelOp, Sense, VarType};
use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

// 对整数需求和费用用整数运算。下列限制也约束探测与动态规划的最坏开销。
const MAX_PERIODS: usize = 64;
const MAX_CHAINS: usize = 256;
const MAX_COMPARISONS: usize = 250_000;
const WORK: usize = 2_000_000;
const BIG: i128 = i128::MAX / 4;

struct Budget {
    left: usize,
    deadline: Option<Instant>,
}

impl Budget {
    fn spend(&mut self, units: usize) -> Option<()> {
        self.left = self.left.checked_sub(units)?;
        (!self.deadline.is_some_and(|t| Instant::now() >= t)).then_some(())
    }
}

struct Balance {
    previous: Option<usize>,
    flows: Vec<usize>,
    demand: i64,
    price: i64,
}

struct Chain {
    periods: Vec<Balance>,
}

fn whole(x: f64) -> Option<i64> {
    (x.is_finite() && x.fract() == 0.0 && x.abs() <= (1u64 << 30) as f64).then_some(x as i64)
}

// 等式中流量与余额系数同号、绝对值相等。只接受能精确还原右端的整单位需求。
fn balance(
    row: &Constraint,
    is_flow: impl Fn(usize) -> bool,
    price: &[i64],
) -> Option<(usize, Balance)> {
    if row.op != RelOp::Eq {
        return None;
    }
    let terms = &row.expr.terms;
    let (_, scale) = terms.iter().copied().find(|&(j, _)| is_flow(j))?;
    if scale.abs() < 1.0 / (1u64 << 30) as f64 || terms.iter().any(|&(_, a)| a.abs() != scale.abs())
    {
        return None;
    }
    let demand = whole(row.rhs / scale)?;
    if demand < 0 || scale.mul_add(demand as f64, -row.rhs) != 0.0 {
        return None;
    }
    let mut flows = Vec::new();
    let mut end = None;
    let mut previous = None;
    for &(j, a) in terms {
        if is_flow(j) {
            if a != scale {
                return None;
            }
            flows.push(j);
        } else if a == scale {
            if end.replace(j).is_some() {
                return None;
            }
        } else if previous.replace(j).is_some() {
            return None;
        }
    }
    let end = end?;
    Some((
        end,
        Balance {
            previous,
            flows,
            demand,
            price: price[end],
        },
    ))
}

fn chains(mut rows: BTreeMap<usize, Balance>) -> Vec<Chain> {
    let mut next = BTreeMap::new();
    let mut starts = Vec::new();
    let mut split = BTreeSet::new();
    for (&j, r) in &rows {
        match r.previous {
            None => starts.push(j),
            Some(prev) if next.insert(prev, j).is_some() => {
                split.insert(prev);
            }
            Some(_) => {}
        }
    }
    let mut chains = Vec::new();
    for mut j in starts {
        let mut periods = Vec::new();
        while let Some(r) = rows.remove(&j) {
            periods.push(r);
            if periods.len() == MAX_PERIODS || split.contains(&j) {
                break;
            }
            let Some(&k) = next.get(&j) else { break };
            j = k;
        }
        if !periods.is_empty() {
            chains.push(Chain { periods });
        }
    }
    chains
}

fn mutex(row: &Constraint) -> Option<(usize, usize)> {
    let &[(a, x), (b, y)] = row.expr.terms.as_slice() else {
        return None;
    };
    // 即使约束是等式，其上界方向也足以排除双一。
    let direction = match row.op {
        RelOp::Le => 1.0,
        RelOp::Ge => -1.0,
        RelOp::Eq if x < 0.0 => -1.0,
        RelOp::Eq => 1.0,
    };
    let (x, y, rhs) = (x * direction, y * direction, row.rhs * direction);
    // 右端可精确表示；只有双一的活动度舍入后仍严格超过它才排斥。
    (x > 0.0 && y > 0.0 && x + y > rhs && (x + y).is_finite()).then_some((a, b))
}

struct Switches {
    required: Vec<Vec<usize>>,
    paid: Vec<Option<usize>>,
    exclusive: BTreeSet<(usize, usize)>,
}

impl Switches {
    fn new(p: &Problem, prices: &[i64], budget: &mut Budget) -> Option<Self> {
        let mut required = vec![Vec::new(); p.vars.len()];
        let mut exclusive = BTreeSet::new();
        for row in &p.constraints {
            budget.spend(1 + row.expr.terms.len())?;
            if let Some((gate, _)) = crate::capacity_gate::at_one(p, row) {
                for &(flow, _) in &row.expr.terms {
                    if flow != gate {
                        required[flow].push(gate);
                    }
                }
            } else if let Some((a, b)) = mutex(row) {
                if [a, b].iter().all(|&j| is_binary(&p.vars[j])) {
                    exclusive.insert((a, b));
                }
            }
        }
        let mut paid = vec![None; p.vars.len()];
        for (j, gates) in required.iter_mut().enumerate() {
            gates.sort_unstable();
            gates.dedup();
            paid[j] = gates
                .iter()
                .copied()
                .filter(|&g| prices[g] > 0)
                .max_by_key(|&g| (prices[g], std::cmp::Reverse(g)));
        }
        Some(Self {
            required,
            paid,
            exclusive,
        })
    }

    fn cannot_share(&self, a: usize, b: usize) -> bool {
        self.required[a].iter().any(|&x| {
            self.required[b]
                .iter()
                .any(|&y| self.exclusive.contains(&(x.min(y), x.max(y))))
        })
    }
}

fn group_chains(
    demand: &[Chain],
    switches: &Switches,
    budget: &mut Budget,
) -> Option<Vec<Vec<usize>>> {
    let mut uses: BTreeMap<usize, Vec<(usize, usize, usize)>> = BTreeMap::new();
    let mut components: Vec<_> = (0..demand.len()).collect();
    for (i, chain) in demand.iter().enumerate() {
        for (t, r) in chain.periods.iter().enumerate() {
            for &j in &r.flows {
                let entries = uses.entry(switches.paid[j]?).or_default();
                budget.spend(1 + entries.len())?;
                for &(k, other_time, flow) in entries.iter() {
                    if other_time != t {
                        return None;
                    } // 同一费用不能跨期重复计算。
                    if components[i] == components[k] {
                        continue;
                    }
                    budget.spend(
                        switches.required[flow]
                            .len()
                            .saturating_mul(switches.required[j].len()),
                    )?;
                    if !switches.cannot_share(j, flow) {
                        let a = components[k];
                        let b = components[i];
                        for component in &mut components {
                            if *component == a {
                                *component = b;
                            }
                        }
                    }
                }
                entries.push((i, t, j));
            }
        }
    }
    let mut groups = BTreeMap::new();
    for (i, g) in components.into_iter().enumerate() {
        groups.entry(g).or_insert_with(Vec::new).push(i);
    }
    Some(groups.into_values().collect())
}

// 同一商品的多个客户共用供给时，余额总和至少是累计需求减去累计到货量。
// 只对各期流量均包含在供给等式内、且彼此没有重复流量的需求链合并计费。
struct Deficit {
    period: usize,
    gap: i128,
    // 按原始单位费用递增排列的（需求链、累计最大余额、原费用）。
    clients: Vec<(usize, i128, i64)>,
}

impl Deficit {
    // 对确知无法满足的数量从最便宜的余额开始分配；各余额不可能大于累计需求。
    // 把计入缺口的价格从动态规划中逐项移除，避免对同一余额重复计费。
    fn assign(&self, rate: i64, demand: &mut [Chain]) -> i128 {
        let mut left = self.gap;
        let mut fee = 0;
        for &(i, cap, price) in &self.clients {
            let share = rate.min(price);
            demand[i].periods[self.period].price = price - share;
            let take = left.min(cap);
            fee += take * i128::from(share);
            left -= take;
        }
        fee
    }
}

fn shortages(
    p: &Problem,
    demand: &[Chain],
    flow: &[bool],
    prices: &[i64],
    budget: &mut Budget,
) -> Option<Vec<Deficit>> {
    let mut supplies = BTreeMap::new();
    for row in &p.constraints {
        budget.spend(1 + row.expr.terms.len())?;
        if let Some((end, r)) = balance(row, |j| flow[j], prices) {
            if prices[end] == 0 && r.previous.is_none_or(|j| prices[j] == 0) {
                supplies.entry(end).or_insert(r);
            }
        }
    }
    let mut out = Vec::new();
    let widths: Vec<usize> = demand
        .iter()
        .map(|c| c.periods.iter().map(|r| 1 + r.flows.len()).sum())
        .collect();
    let mut used = vec![false; demand.len()];
    for supply in chains(supplies) {
        let sets: Vec<BTreeSet<_>> = supply
            .periods
            .iter()
            .map(|r| r.flows.iter().copied().collect())
            .collect();
        let mut seen = BTreeSet::new();
        let mut clients = Vec::new();
        for (i, c) in demand.iter().enumerate() {
            budget.spend(widths[i] + 1)?;
            if used[i]
                || c.periods.len() != sets.len()
                || !c
                    .periods
                    .iter()
                    .zip(&sets)
                    .all(|(r, s)| r.flows.iter().all(|j| s.contains(j) && !seen.contains(j)))
            {
                continue;
            }
            // 同一需求链内各期也要拥有独立的流量列。
            let count: usize = c.periods.iter().map(|r| r.flows.len()).sum();
            let ids: BTreeSet<_> = c
                .periods
                .iter()
                .flat_map(|r| r.flows.iter().copied())
                .collect();
            if count != ids.len() {
                continue;
            }
            seen.extend(ids);
            clients.push(i);
        }
        if clients.is_empty() {
            continue;
        }
        let mut gap = 0i128;
        let mut caps = vec![0i128; clients.len()];
        for (t, r) in supply.periods.iter().enumerate() {
            gap -= i128::from(r.demand);
            for (&i, cap) in clients.iter().zip(&mut caps) {
                let d = i128::from(demand[i].periods[t].demand);
                gap += d;
                *cap += d;
            }
            if gap > 0 {
                let mut entries: Vec<_> = clients
                    .iter()
                    .zip(&caps)
                    .map(|(&i, &c)| (i, c, demand[i].periods[t].price))
                    .collect();
                entries.sort_by_key(|&(i, _, price)| (price, i));
                out.push(Deficit {
                    period: t,
                    gap,
                    clients: entries,
                });
            }
        }
        for i in clients {
            used[i] = true;
        }
    }
    Some(out)
}

struct Group {
    flows: Vec<Vec<usize>>,
    // carry[t][s] 是从 s 到 t 均未开通时，第 t 期必付的余额费用。
    carry: Vec<Vec<i128>>,
}

impl Group {
    fn new(chains: &[Chain], ids: &[usize]) -> Self {
        let len = ids.iter().map(|&i| chains[i].periods.len()).max().unwrap();
        let mut flows = vec![Vec::new(); len];
        let mut carry: Vec<Vec<i128>> = (0..len).map(|t| vec![0; t + 1]).collect();
        for &i in ids {
            let chain = &chains[i].periods;
            for (t, r) in chain.iter().enumerate() {
                flows[t].extend_from_slice(&r.flows);
                let mut amount = 0;
                for s in (0..=t).rev() {
                    amount += i128::from(chain[s].demand);
                    carry[t][s] += amount * i128::from(r.price);
                }
            }
        }
        Self { flows, carry }
    }

    fn fees(&self, switches: &Switches, price: &[i64], off: &[bool]) -> Vec<i128> {
        self.flows
            .iter()
            .map(|row| {
                row.iter()
                    .filter(|&&j| switches.required[j].iter().all(|&g| !off[g]))
                    .map(|&j| i128::from(price[switches.paid[j].unwrap()]))
                    .min()
                    .unwrap_or(BIG)
            })
            .collect()
    }

    fn minimum(&self, fees: &[i128]) -> i128 {
        let mut values = vec![0i128];
        for (r, &fee) in self.carry.iter().zip(fees) {
            let serve = values.iter().copied().min().unwrap() + fee;
            for (v, &holding) in values.iter_mut().zip(r) {
                *v += holding;
            }
            values.push(serve.min(BIG));
        }
        values.into_iter().min().unwrap()
    }
}

/// 返回原始目标方向下的下界；不使用求解器当前解或外部数据。
pub(crate) fn objective_bound(p: &Problem, deadline: Option<Instant>) -> Option<f64> {
    let expired = || deadline.is_some_and(|end| Instant::now() >= end);
    let sign = if p.sense == Sense::Minimize {
        1.0
    } else {
        -1.0
    };
    if expired()
        || p.vars.len() > MAX_COMPARISONS
        || p.constraints.len() > MAX_COMPARISONS
        || !p
            .obj
            .terms
            .iter()
            .any(|&(j, c)| c * sign > 0.0 && p.vars[j].var_type == VarType::Continuous)
        || !p
            .obj
            .terms
            .iter()
            .any(|&(j, c)| c * sign > 0.0 && is_binary(&p.vars[j]))
        || p.vars.iter().any(|v| v.lb < 0.0)
    {
        return None;
    }
    let mut budget = Budget {
        left: MAX_COMPARISONS,
        deadline,
    };
    let mut price = vec![0; p.vars.len()];
    for &(j, c) in &p.obj.terms {
        price[j] = whole(c * sign)?;
        if price[j] < 0 {
            return None;
        }
    }
    let offset = i128::from(whole(p.obj_offset * sign)?);
    let switches = Switches::new(p, &price, &mut budget)?;
    let mut rows = BTreeMap::new();
    for r in &p.constraints {
        budget.spend(1 + r.expr.terms.len())?;
        if let Some((end, r)) = balance(r, |j| switches.paid[j].is_some(), &price) {
            if price[end] > 0
                && p.vars[end].var_type == VarType::Continuous
                && r.previous
                    .is_none_or(|j| price[j] > 0 && p.vars[j].var_type == VarType::Continuous)
            {
                rows.entry(end).or_insert(r);
            }
        }
    }
    let mut demand = chains(rows);
    if demand.is_empty() || demand.len() > MAX_CHAINS || expired() {
        return None;
    }
    let ids = group_chains(&demand, &switches, &mut budget)?;
    let mut flows = vec![false; p.vars.len()];
    for chain in &demand {
        for r in &chain.periods {
            for &j in &r.flows {
                flows[j] = true;
            }
        }
    }
    let deficits = shortages(p, &demand, &flows, &price, &mut budget)?;
    let mut rates: Vec<_> = deficits.iter().map(|d| d.clients[0].2).collect();
    let mut terms: Vec<_> = deficits
        .iter()
        .zip(&rates)
        .map(|(d, &r)| d.assign(r, &mut demand))
        .collect();
    let rebuild: usize = demand
        .iter()
        .map(|c| {
            c.periods
                .iter()
                .enumerate()
                .map(|(t, r)| r.flows.len() + t + 1)
                .sum::<usize>()
        })
        .sum();
    budget.spend(rebuild)?;
    let groups: Vec<_> = ids.iter().map(|g| Group::new(&demand, g)).collect();
    let work_per_mask: usize = groups
        .iter()
        .map(|g| {
            g.flows
                .iter()
                .flatten()
                .map(|&j| switches.required[j].len())
                .sum::<usize>()
                + g.carry.iter().map(Vec::len).sum::<usize>()
        })
        .sum();
    budget.spend(work_per_mask)?;
    let mut off = vec![false; p.vars.len()];
    let base: Vec<_> = groups
        .iter()
        .map(|g| g.fees(&switches, &price, &off))
        .collect();
    let compute = |groups: &[Group], fees: &[Vec<i128>]| {
        groups
            .iter()
            .zip(fees)
            .map(|(g, f)| g.minimum(f))
            .sum::<i128>()
    };
    let constant = offset + terms.iter().sum::<i128>();
    let mut bound = constant + compute(&groups, &base);
    let convert = |v: i128| (v.unsigned_abs() <= (1u128 << 52)).then_some((v as f64) * sign);
    // 两个互斥开关至少有一个关闭。只枚举两个方向都增加最低开通费的互斥对；
    // 任何中断均退回已完成的整棵枚举树所给出的下界。
    let scenarios = WORK.checked_div(work_per_mask.max(1)).unwrap_or(0);
    if scenarios < 2 || expired() {
        return convert(bound);
    }
    let mut pairs = Vec::new();
    for &(a, b) in switches.exclusive.iter().take(scenarios / 2) {
        if expired() {
            return convert(bound);
        }
        let mut scores = [0i128; 2];
        for (s, &bit) in [a, b].iter().enumerate() {
            off[bit] = true;
            for (g, original) in groups.iter().zip(&base) {
                let fs = g.fees(&switches, &price, &off);
                scores[s] += original
                    .iter()
                    .zip(fs)
                    .map(|(&x, y)| (y - x).min(1i128 << 40))
                    .sum::<i128>();
            }
            off[bit] = false;
        }
        if scores[0].min(scores[1]) > 0 {
            pairs.push((std::cmp::Reverse(scores[0].min(scores[1])), a, b));
        }
    }
    pairs.sort_unstable();
    let bits = (scenarios.ilog2() as usize).min(pairs.len()).min(16);
    let mut cases = Vec::new();
    let mut worst = BIG;
    for mask in 0..1usize << bits {
        if expired() {
            return convert(bound);
        }
        for (s, &(_, a, b)) in pairs.iter().take(bits).enumerate() {
            off[if mask & (1 << s) == 0 { a } else { b }] = true;
        }
        let fees: Vec<_> = groups
            .iter()
            .map(|g| g.fees(&switches, &price, &off))
            .collect();
        worst = worst.min(compute(&groups, &fees));
        cases.push(fees);
        for &(_, a, b) in pairs.iter().take(bits) {
            off[a] = false;
            off[b] = false;
        }
    }
    bound = bound.max(constant + worst);
    // 分配给短缺的每一种非负费率都给出独立有效的界；在单位费用的折点上比较。
    let visits = groups
        .iter()
        .map(|g| g.carry.iter().map(Vec::len).sum::<usize>())
        .sum::<usize>();
    let mut tune = Budget {
        left: WORK * 2,
        deadline,
    };
    for (k, deficit) in deficits.iter().enumerate() {
        let candidates: BTreeSet<_> = std::iter::once(0)
            .chain(deficit.clients.iter().map(|c| c.2))
            .collect();
        for rate in candidates {
            if rate == rates[k] {
                continue;
            }
            if tune
                .spend(rebuild + visits.saturating_mul(cases.len()))
                .is_none()
            {
                return convert(bound);
            }
            let cost = deficit.assign(rate, &mut demand);
            let changed: Vec<_> = ids.iter().map(|g| Group::new(&demand, g)).collect();
            let mut best = BIG;
            for case in &cases {
                if expired() {
                    return convert(bound);
                }
                best = best.min(compute(&changed, case));
            }
            let candidate = offset + terms.iter().sum::<i128>() - terms[k] + cost + best;
            if candidate > bound {
                bound = candidate;
                rates[k] = rate;
                terms[k] = cost;
            }
        }
        deficit.assign(rates[k], &mut demand);
    }
    convert(bound)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instance::parse_instance;

    // 两个需求分别为 3、7，单位费用不同，共用同一个付费开关和有限供给。
    fn small(supply: usize, conflict: bool, maximize: bool) -> Problem {
        let direction = if maximize { "max" } else { "min" };
        let sign = if maximize { -1 } else { 1 };
        let exclusion = if conflict { "m: -2 a -2 b >= -2" } else { "" };
        let data = format!("name example\nsense {direction}\nvars\nf continuous\ng continuous\nh continuous\nk continuous\ns continuous\ny integer\na binary\nb binary\nobj\n{} h {} k {} y\ncons\n\
d0: -2 f -2 h = -6\nd1: 0.5 g 0.5 k = 3.5\n\
p: 4 f 4 g 4 s = {}\ncap: 1 f 1 g -10 y <= 0\n\
ga: -1 f 10 a >= 0\ngb: 1 g -10 b <= 0\n{exclusion}\nbounds\ny <= 1\nend\n",
            90*sign, 50*sign,100*sign,supply*4);
        let mut p = parse_instance(&data).unwrap();
        p.obj_offset = 17.0;
        p
    }

    #[test]
    fn shared_charge_mutex_shortage_and_objective_direction() {
        for supply in 0..=10 {
            for conflict in [false, true] {
                // 固定开关时可行流量是一个矩形与整数容量半平面的交，
                // 顶点都是整数，因此逐个列举流量也涵盖了连续最优解。
                let mut exact = i64::MAX;
                for a in 0..=3 {
                    for b in 0..=7 {
                        if a + b <= supply && !(conflict && a > 0 && b > 0) {
                            exact = exact.min(
                                90 * (3 - a as i64)
                                    + 50 * (7 - b as i64)
                                    + 100 * i64::from(a + b > 0),
                            );
                        }
                    }
                }
                let min = objective_bound(&small(supply, conflict, false), None).unwrap();
                let max = objective_bound(&small(supply, conflict, true), None).unwrap();
                assert!(min <= 17.0 + exact as f64);
                assert_eq!(17.0 - min, max - 17.0);
                if supply == 10 || supply == 0 {
                    assert_eq!(min, 17.0 + exact as f64);
                }
            }
        }
    }

    #[test]
    fn one_fee_spanning_two_periods_is_never_charged_twice() {
        let p = parse_instance("name horizon\nsense min\nvars\nu continuous\nv continuous\nh continuous\nk continuous\ny binary\nobj\n100 h 100 k 5 y\ncons\na: 1 u 1 h = 1\nb: -1 h 1 v 1 k = 1\nc: 1 u 1 v -2 y <= 0\nend\n").unwrap();
        assert_eq!(objective_bound(&p, None), None);
        assert_eq!(
            objective_bound(&small(5, false, false), Some(Instant::now())),
            None
        );
    }

    #[test]
    fn interrupted_grouping_discards_partial_fee_partition() {
        let chains = vec![Chain {
            periods: vec![Balance {
                previous: None,
                flows: (0..2000).collect(),
                demand: 1,
                price: 1,
            }],
        }];
        let switches = Switches {
            required: vec![vec![2000]; 2000],
            paid: vec![Some(2000); 2000],
            exclusive: BTreeSet::new(),
        };
        let mut budget = Budget {
            left: 5000,
            deadline: None,
        };
        assert!(group_chains(&chains, &switches, &mut budget).is_none());
    }
}
