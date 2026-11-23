//! 带时间窗的路线覆盖。完整枚举后用路线 LP 找对偶，精确复核对偶和剩余搜索。
//
//! 任何原整数解都是恰选一图上的置换。带正时延的任务弧不能组成纯任务环；
//! 删除站点后，每个置换只留下跳过的任务以及从站点出发并回站的任务链。
//! 对称且完整的站点可按路线分配，每条链用一个站点，其余站点选零成本自环。
//! 因此站点容量和逐任务恰覆盖描述了必要的组合约束；构造的上界另用原模型检查。

use crate::model::{dyadic, Constraint, LinearExpr, Problem, RelOp, Sense, Tolerances, VarType, Variable};
use crate::simplex::{Bounds, LpStatus};
use crate::simplex_sparse::{new_engine, LpEngineKind};
use crate::time_window_routes::RouteModel;
use std::collections::{HashMap, HashSet};
use std::time::Instant;

const MAX_JOBS: usize = 96;
const MAX_ROUTES: usize = 8192;
const MAX_MEMO: usize = 250_000;

type Exact = i128;

fn expired(deadline: Option<Instant>) -> bool {
    deadline.is_some_and(|end| Instant::now() >= end)
}

fn grid_exponent(mut values: impl Iterator<Item = f64>) -> Option<i32> {
    let exponent = values.try_fold(0, |acc, value| {
        let (_, exponent) = dyadic(value)?;
        Some(acc.min(exponent))
    })?;
    (-900..=0).contains(&exponent).then_some(exponent)
}

fn grid_value(value: f64, exponent: i32) -> Option<Exact> {
    let (mantissa, power) = dyadic(value)?;
    if mantissa == 0 {
        return Some(0);
    }
    let factor = 1_i128.checked_shl(u32::try_from(power - exponent).ok()?)?;
    // 给时间链和路径目标的多项加和留出空间。
    let ticks = mantissa.checked_mul(factor)?;
    (factor > 0 && ticks.unsigned_abs() < (1 << 90)).then_some(ticks)
}

/// 将有理格点值向下取整为浮点，供只认原方向有限对偶界的调用方使用。
fn bound_as_float(value: Exact, denominator: Exact, exponent: i32) -> Option<f64> {
    let integer = (value as f64).next_down();
    let scale = 2.0_f64.powi(exponent);
    let rounded_divisor = denominator as f64;
    let divisor = if integer < 0.0 {
        rounded_divisor.next_down()
    } else {
        rounded_divisor.next_up()
    };
    // 在除法和缩放合并时也向下保守取界，包括下溢。
    let result = ((integer / divisor).next_down() * scale).next_down();
    result.is_finite().then_some(result)
}

#[derive(Default)]
pub(crate) struct RouteOutcome {
    pub(crate) nodes: u64,
    pub(crate) lp_iterations: u64,
    pub(crate) proven_optimal: bool,
    /// 最小化方向，未加目标常数。只有精确复核的界会写入。
    pub(crate) bound: Option<f64>,
    pub(crate) interrupted: bool,
    /// 经过原模型检查的解。证明可能沿用已知的 seed。
    pub(crate) best: Option<Vec<f64>>,
}

struct Route {
    jobs: u128,
    first: usize,
    last: usize,
    inside: Vec<usize>,
    uses_station: bool,
    cost: Exact,
}

struct RouteNetwork<'a> {
    model: &'a RouteModel,
    jobs: Vec<usize>,
    stations: Vec<usize>,
    costs: Vec<Exact>,
    windows: Vec<(Exact, Exact)>,
    duration: Vec<Vec<Option<Exact>>>,
    cost_exp: i32,
}

impl<'a> RouteNetwork<'a> {
    fn new(model: &'a RouteModel) -> Option<Self> {
        let (jobs, stations): (Vec<_>, Vec<_>) =
            (0..model.timing.len()).partition(|&j| model.timing[j].is_some());
        if jobs.is_empty() || jobs.len() > MAX_JOBS || stations.is_empty() {
            return None;
        }

        for &station in &stations {
            for &other in &stations {
                if model.costs[model.arcs[station][other]?] != 0.0 {
                    return None;
                }
            }
            let first_station = stations[0];
            for &job in &jobs {
                let from = model.arcs[station][job]?;
                let to = model.arcs[job][station]?;
                if model.costs[from] != model.costs[model.arcs[first_station][job]?]
                    || model.costs[to] != model.costs[model.arcs[job][first_station]?]
                {
                    return None;
                }
            }
        }

        let cost_exp = grid_exponent(model.costs.iter().copied())?;
        let costs = model
            .costs
            .iter()
            .map(|&c| grid_value(c, cost_exp))
            .collect::<Option<_>>()?;
        let values = model
            .timing
            .iter()
            .flatten()
            .flat_map(|w| [w.lower, w.upper]);
        let coefficients = model
            .duration
            .iter()
            .flatten()
            .flatten()
            .flat_map(|&(a, b)| [a, b]);
        let time_exp = grid_exponent(values.chain(coefficients))?;
        let windows = jobs
            .iter()
            .map(|&j| {
                let w = model.timing[j]?;
                Some((
                    grid_value(w.lower, time_exp)?,
                    grid_value(w.upper, time_exp)?,
                ))
            })
            .collect::<Option<_>>()?;
        let mut duration = vec![vec![None; jobs.len()]; jobs.len()];
        for (i, &from) in jobs.iter().enumerate() {
            for (j, &to) in jobs.iter().enumerate() {
                if i != j && model.arcs[from][to].is_some() {
                    let (rhs, coefficient) = model.duration[from][to]?;
                    let delay = grid_value(rhs, time_exp)? - grid_value(coefficient, time_exp)?;
                    if delay <= 0 {
                        return None;
                    }
                    duration[i][j] = Some(delay);
                }
            }
        }
        Some(Self {
            model,
            jobs,
            stations,
            costs,
            windows,
            duration,
            cost_exp,
        })
    }

    fn enumerate(&self, deadline: Option<Instant>) -> Option<Vec<Route>> {
        let mut routes = Vec::new();
        for (i, &job) in self.jobs.iter().enumerate() {
            if let Some(col) = self.model.arcs[job][job] {
                routes.push(Route {
                    jobs: 1 << i,
                    first: i,
                    last: i,
                    inside: vec![col],
                    uses_station: false,
                    cost: self.costs[col],
                });
            }
        }
        let station = self.stations[0];
        for (i, &job) in self.jobs.iter().enumerate() {
            let col = self.model.arcs[station][job]?;
            let mut chain = Route {
                jobs: 1 << i,
                first: i,
                last: i,
                inside: Vec::new(),
                uses_station: true,
                cost: self.costs[col],
            };
            self.extend(&mut chain, self.windows[i].0, &mut routes, deadline)?;
        }
        Some(routes)
    }

    fn extend(
        &self,
        chain: &mut Route,
        earliest: Exact,
        routes: &mut Vec<Route>,
        deadline: Option<Instant>,
    ) -> Option<()> {
        if routes.len() == MAX_ROUTES || expired(deadline) {
            return None;
        }
        let last = chain.last;
        let from = self.jobs[last];
        let to_station = self.model.arcs[from][self.stations[0]]?;
        routes.push(Route {
            jobs: chain.jobs,
            first: chain.first,
            last,
            inside: chain.inside.clone(),
            uses_station: true,
            cost: chain.cost + self.costs[to_station],
        });
        for (next, delay) in self.duration[last].iter().enumerate() {
            let bit = 1 << next;
            if chain.jobs & bit != 0 {
                continue;
            }
            if let Some(delay) = delay {
                // 最早到达时间恰好是这些单向差分约束的可行性判据。
                // 这是原始 f64 系数的精确整数和，临界等号也必须保留。
                let (lower, upper) = self.windows[next];
                let arrival = lower.max(earliest + delay);
                if arrival <= upper {
                    let col = self.model.arcs[from][self.jobs[next]]?;
                    chain.last = next;
                    chain.jobs |= bit;
                    chain.cost += self.costs[col];
                    chain.inside.push(col);
                    self.extend(chain, arrival, routes, deadline)?;
                    chain.inside.pop();
                    chain.cost -= self.costs[col];
                    chain.jobs ^= bit;
                    chain.last = last;
                }
            }
        }
        Some(())
    }

    fn reconstruct(&self, routes: &[Route], selected: &[usize]) -> Option<Vec<usize>> {
        let mut columns = Vec::new();
        let mut stations = self.stations.iter();
        for &id in selected {
            let route = &routes[id];
            if route.uses_station {
                let &station = stations.next()?;
                columns.push(self.model.arcs[station][self.jobs[route.first]]?);
                columns.push(self.model.arcs[self.jobs[route.last]][station]?);
            }
            columns.extend(&route.inside);
        }
        for &station in stations {
            columns.push(self.model.arcs[station][station]?);
        }
        Some(columns)
    }
}

fn route_lp(routes: &[Route], jobs: usize, stations: usize, cost_exp: i32) -> Problem {
    let mut problem = Problem::new("路线覆盖", Sense::Minimize);
    problem.vars = (0..routes.len())
        .map(|i| Variable {
            name: format!("p{i}"),
            var_type: VarType::Continuous,
            lb: 0.0,
            ub: f64::INFINITY,
        })
        .collect();
    let mut rows = vec![LinearExpr::new(); jobs + 1];
    let step = 2.0_f64.powi(cost_exp);
    for (i, path) in routes.iter().enumerate() {
        if path.cost != 0 {
            problem.obj.terms.push((i, path.cost as f64 * step));
        }
        if path.uses_station {
            rows[0].terms.push((i, 1.0));
        }
        for j in Bits(path.jobs) {
            rows[j + 1].terms.push((i, 1.0));
        }
    }
    problem.constraints = rows
        .into_iter()
        .enumerate()
        .map(|(i, expr)| Constraint {
            name: None,
            expr,
            op: if i == 0 { RelOp::Le } else { RelOp::Eq },
            rhs: if i == 0 { stations as f64 } else { 1.0 },
        })
        .collect();
    problem
}

struct Bits(u128);

impl Iterator for Bits {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.0 == 0 {
            return None;
        }
        let bit = self.0.trailing_zeros() as usize;
        self.0 &= self.0 - 1;
        Some(bit)
    }
}

fn exact_div(numerator: Exact, denominator: Exact) -> Option<Exact> {
    let result = numerator.checked_div(denominator)?;
    (result.checked_mul(denominator)? == numerator).then_some(result)
}

/// Bareiss 消元保持所有中间项为整数。返回对偶各项的共同分母和整数分子。
fn solve_exact_system(
    mut rows: Vec<Vec<Exact>>,
    deadline: Option<Instant>,
) -> Option<(Exact, Vec<Exact>)> {
    let n = rows.len();
    let mut previous = 1;
    for k in 0..n {
        if expired(deadline) {
            return None;
        }
        let pivot_row = (k..n).find(|&i| rows[i][k] != 0)?;
        rows.swap(k, pivot_row);
        let pivot = rows[k][k];
        for i in k + 1..n {
            let left = rows[i][k];
            for j in k + 1..=n {
                let value = pivot
                    .checked_mul(rows[i][j])?
                    .checked_sub(left.checked_mul(rows[k][j])?)?;
                rows[i][j] = exact_div(value, previous)?;
            }
            rows[i][k] = 0;
        }
        previous = pivot;
    }
    let denominator = previous.checked_abs()?;
    let mut solution = vec![0; n];
    for i in (0..n).rev() {
        let remainder = rows[i][n].checked_mul(denominator)?;
        let numerator = (i + 1..n).try_fold(remainder, |acc, j| {
            acc.checked_sub(rows[i][j].checked_mul(solution[j])?)
        })?;
        solution[i] = exact_div(numerator, rows[i][i])?;
    }
    Some((denominator, solution))
}

struct Certificate {
    denominator: Exact,
    dual_bound: Exact,
    station_penalty: Exact,
    reduced_costs: Vec<Exact>,
}

impl Certificate {
    fn from_basis(
        basis: &[usize],
        routes: &[Route],
        jobs: usize,
        stations: usize,
        deadline: Option<Instant>,
    ) -> Option<Self> {
        if basis.len() != jobs + 1 {
            return None;
        }
        let mut matrix = vec![vec![0; jobs + 2]; jobs + 1];
        for (row, &col) in matrix.iter_mut().zip(basis) {
            if let Some(route) = routes.get(col) {
                row[0] = Exact::from(route.uses_station);
                for j in Bits(route.jobs) {
                    row[j + 1] = 1;
                }
                row[jobs + 1] = route.cost;
            } else if col == routes.len() {
                row[0] = 1; // 唯一的松弛列为站点容量。
            } else {
                return None; // 残基人工列只可回退。
            }
        }
        let (denominator, dual) = solve_exact_system(matrix, deadline)?;
        if dual[0] > 0 {
            return None;
        }
        let dual_bound = dual[1..]
            .iter()
            .try_fold(dual[0].checked_mul(stations as Exact)?, |acc, &d| {
                acc.checked_add(d)
            })?;
        let mut reduced_costs = Vec::with_capacity(routes.len());
        for route in routes {
            let cost = route.cost.checked_mul(denominator)?;
            let base = if route.uses_station {
                cost.checked_sub(dual[0])?
            } else {
                cost
            };
            let reduced = Bits(route.jobs).try_fold(base, |acc, j| acc.checked_sub(dual[j + 1]))?;
            if reduced < 0 {
                return None;
            }
            reduced_costs.push(reduced);
        }
        Some(Self {
            denominator,
            dual_bound,
            station_penalty: dual[0].checked_neg()?,
            reduced_costs,
        })
    }

    fn bound(&self, extra: Exact, exponent: i32) -> Option<f64> {
        bound_as_float(
            self.dual_bound.checked_add(extra)?,
            self.denominator,
            exponent,
        )
    }
}

#[derive(Clone, Copy)]
struct SearchRoute {
    jobs: u128,
    station: usize,
    cost: Exact,
    share: Exact,
}

struct JobClass {
    jobs: u128,
    // 下标为剩余任务数，对应保留该组最小的若干任务编号。
    prefixes: Vec<u128>,
}

struct CoverSearch {
    routes: Vec<SearchRoute>,
    by_job: Vec<Vec<usize>>,
    stations: usize,
    penalty: Exact,
    classes: Vec<JobClass>,
    seen: HashMap<(u128, usize), Exact>,
    nodes: u64,
    node_limit: u64,
    deadline: Option<Instant>,
    interrupted: bool,
    best: Exact,
    selected: Option<Vec<usize>>,
}

impl CoverSearch {
    fn new(
        jobs: usize,
        stations: usize,
        routes: &[Route],
        proof: &Certificate,
        deadline: Option<Instant>,
        node_limit: u64,
    ) -> Self {
        let routes: Vec<_> = routes
            .iter()
            .zip(&proof.reduced_costs)
            .map(|(r, &cost)| SearchRoute {
                jobs: r.jobs,
                station: usize::from(r.uses_station),
                cost,
                share: cost / Exact::from(r.jobs.count_ones()),
            })
            .collect();
        let mut by_job = vec![Vec::new(); jobs];
        for (i, r) in routes.iter().enumerate() {
            for j in Bits(r.jobs) {
                by_job[j].push(i);
            }
        }
        for paths in &mut by_job {
            paths.sort_unstable_by_key(|&i| (routes[i].cost, i));
        }
        Self {
            routes,
            by_job,
            stations,
            penalty: proof.station_penalty,
            classes: Vec::new(),
            seen: HashMap::new(),
            nodes: 0,
            node_limit,
            deadline,
            interrupted: false,
            best: 0,
            selected: None,
        }
    }

    /// 两个任务互换后，每条可用路线都须有任务集合、站点占用和减价完全
    /// 对应的路线。同组内与代表的两两交换能生成该组的任意排列。
    fn interchangeable_jobs(&self, cutoff: Exact) -> Option<Vec<JobClass>> {
        let active: HashSet<_> = self
            .routes
            .iter()
            .filter(|r| r.cost < cutoff)
            .map(|r| (r.jobs, r.station, r.cost))
            .collect();
        let by_job: Vec<_> = self
            .by_job
            .iter()
            .map(|ids| {
                let end = ids.partition_point(|&i| self.routes[i].cost < cutoff);
                &ids[..end]
            })
            .collect();
        let mut classes = Vec::new();
        let mut assigned = 0u128;
        for a in 0..by_job.len() {
            if expired(self.deadline) {
                return None;
            }
            if assigned & (1 << a) != 0 || by_job[a].is_empty() {
                continue;
            }
            let mut jobs = 1u128 << a;
            let mut prefixes = vec![0, jobs];
            for b in a + 1..by_job.len() {
                // 列数不同可以保守跳过这对，重复列只会减少识别机会。
                if assigned & (1 << b) != 0 || by_job[a].len() != by_job[b].len() {
                    continue;
                }
                let pair = (1u128 << a) | (1u128 << b);
                let interchangeable = by_job[a].iter().chain(by_job[b]).all(|&i| {
                    let r = self.routes[i];
                    r.jobs & pair == pair || active.contains(&(r.jobs ^ pair, r.station, r.cost))
                });
                if interchangeable {
                    jobs |= 1u128 << b;
                    prefixes.push(jobs);
                }
            }
            assigned |= jobs;
            if prefixes.len() > 2 {
                classes.push(JobClass { jobs, prefixes });
            }
        }
        Some(classes)
    }

    fn run(&mut self, cutoff: Exact) {
        if self.nodes >= self.node_limit {
            self.interrupted = true;
            return;
        }
        let Some(classes) = self.interchangeable_jobs(cutoff) else {
            self.interrupted = true;
            return;
        };
        self.classes = classes;
        self.seen.clear();
        self.best = cutoff;
        let remaining = (1 << self.by_job.len()) - 1;
        self.visit(remaining, 0, 0, &mut Vec::new());
    }

    fn visit(&mut self, remaining: u128, stations: usize, used: Exact, choices: &mut Vec<usize>) {
        if used >= self.best {
            return;
        }
        if self.nodes >= self.node_limit || expired(self.deadline) {
            self.interrupted = true;
            return;
        }
        self.nodes += 1;
        if remaining == 0 {
            let total = used + (self.stations - stations) as Exact * self.penalty;
            if total < self.best {
                self.best = total;
                self.selected = Some(choices.clone());
            }
            return;
        }

        // 等价只用于记忆已搜状态。真实路线编号保留在 choices 中供原模型复核。
        let mut canonical = remaining;
        for class in &self.classes {
            let count = (remaining & class.jobs).count_ones() as usize;
            canonical = (canonical & !class.jobs) | class.prefixes[count];
        }
        let state = (canonical, stations);
        if self
            .seen
            .get(&state)
            .is_some_and(|&previous| previous <= used)
        {
            return;
        }
        if self.seen.len() < MAX_MEMO {
            self.seen.insert(state, used);
        }

        // 每个剩余任务至少支付包含它的最小路线份额，同一路线各份额之和
        // 不超过该路线减价。任一任务没有可用路线的状态也可以直接剪去。
        let gap = self.best - used;
        let mut share = 0;
        let mut fewest = (usize::MAX, usize::MAX);
        for job in Bits(remaining) {
            let mut count = 0;
            let mut cheapest = Exact::MAX;
            for &i in &self.by_job[job] {
                let route = self.routes[i];
                if route.cost >= gap {
                    break;
                }
                if route.jobs & remaining == route.jobs && stations + route.station <= self.stations
                {
                    count += 1;
                    cheapest = cheapest.min(route.share);
                }
            }
            if count == 0 {
                return;
            }
            fewest = fewest.min((count, job));
            share += cheapest;
            if share >= gap {
                return;
            }
        }

        // 覆盖选中任务的所有合法路线给出不遗漏的分支。路径列号稳定破并列。
        let job = fewest.1;
        for k in 0..self.by_job[job].len() {
            let i = self.by_job[job][k];
            let route = self.routes[i];
            if route.cost >= self.best - used {
                break;
            }
            if route.jobs & remaining == route.jobs && stations + route.station <= self.stations {
                choices.push(i);
                self.visit(
                    remaining ^ route.jobs,
                    stations + route.station,
                    used + route.cost,
                    choices,
                );
                choices.pop();
                if self.interrupted {
                    return;
                }
            }
        }
    }
}

/// 完整路线上先由浮点 LP 选择基，再以原始系数的格点检查建立可用界。
/// 未识别的结构、精确算术溢出或无法复核的候选，交还通用求解器。
pub(crate) fn optimize(
    problem: &Problem,
    model: &RouteModel,
    tol: Tolerances,
    deadline: Option<Instant>,
    node_limit: u64,
    seed: Option<&[f64]>,
) -> RouteOutcome {
    let mut outcome = RouteOutcome::default();
    if node_limit == 0 || expired(deadline) {
        outcome.interrupted = true;
        return outcome;
    }
    let _ = optimize_inner(
        problem,
        model,
        tol,
        deadline,
        node_limit,
        seed,
        &mut outcome,
    );
    if expired(deadline) {
        outcome.interrupted = true;
        outcome.proven_optimal = false;
    }
    outcome
}

fn optimize_inner(
    problem: &Problem,
    model: &RouteModel,
    tol: Tolerances,
    deadline: Option<Instant>,
    node_limit: u64,
    seed: Option<&[f64]>,
    outcome: &mut RouteOutcome,
) -> Option<()> {
    let network = RouteNetwork::new(model)?;
    let routes = network.enumerate(deadline)?;
    let lp = route_lp(
        &routes,
        network.jobs.len(),
        network.stations.len(),
        network.cost_exp,
    );
    // 路线只经过少数任务。先用稀疏引擎选基，若最优基无法精确复核，
    // 再尝试稠密基。两次计算同一个根松弛，共用墙钟并累计迭代数。
    let bounds = Bounds::from_problem(&lp);
    let mut certificate = None;
    for kind in [LpEngineKind::Revised, LpEngineKind::Tableau] {
        if expired(deadline) {
            return None;
        }
        let Ok(mut simplex) = new_engine(&lp, tol, kind) else {
            continue;
        };
        simplex.set_time_deadline(deadline);
        let relaxed = simplex.solve(&bounds);
        outcome.nodes = 1;
        outcome.lp_iterations += simplex.iterations();
        if relaxed.status != LpStatus::Optimal {
            return None;
        }
        certificate = Certificate::from_basis(
            simplex.basis_columns(),
            &routes,
            network.jobs.len(),
            network.stations.len(),
            deadline,
        );
        if certificate.is_some() {
            break;
        }
    }
    let certificate = certificate?;
    outcome.bound = certificate.bound(0, network.cost_exp);

    let seed = seed.and_then(|seed| {
        if seed.len() != model.costs.len() {
            return None;
        }
        let columns: Vec<_> = seed
            .iter()
            .enumerate()
            .filter(|&(i, x)| *x == 1.0 && problem.vars[i].var_type == VarType::Binary)
            .map(|(i, _)| i)
            .collect();
        let values = model.candidate_from_arcs(problem, &columns)?;
        let cost = columns.iter().map(|&i| network.costs[i]).sum::<Exact>();
        let residual = cost
            .checked_mul(certificate.denominator)?
            .checked_sub(certificate.dual_bound)?;
        (residual >= 0).then_some((residual, values))
    });

    let largest = certificate.reduced_costs.iter().copied().max()?;
    // 这也界住搜索里 used、share 和惩罚的中间和。每条路线至少覆盖一个任务。
    let max_residual = largest
        .checked_mul(network.jobs.len() as Exact)?
        .checked_add(
            certificate
                .station_penalty
                .checked_mul(network.stations.len() as Exact)?,
        )?;
    let final_cutoff = seed
        .as_ref()
        .map(|(cost, _)| *cost)
        .unwrap_or(max_residual.checked_add(1)?);
    let first_positive = certificate
        .reduced_costs
        .iter()
        .copied()
        .chain([certificate.station_penalty])
        .filter(|&cost| cost > 0)
        .min();
    let mut search = CoverSearch::new(
        network.jobs.len(),
        network.stations.len(),
        &routes,
        &certificate,
        deadline,
        node_limit - 1,
    );
    let mut cutoff = 1.min(final_cutoff);
    loop {
        if cutoff > 0 {
            search.run(cutoff);
        }
        outcome.nodes = 1 + search.nodes;
        outcome.interrupted = search.interrupted || expired(deadline);
        if let Some(selected) = search.selected.take() {
            let choices = network.reconstruct(&routes, &selected)?;
            let values = model.candidate_from_arcs(problem, &choices)?;
            let cost = choices.iter().map(|&i| network.costs[i]).sum::<Exact>();
            let exact = cost.checked_mul(certificate.denominator)?;
            if exact != certificate.dual_bound.checked_add(search.best)? {
                return None;
            }
            outcome.best = Some(values);
            if !outcome.interrupted {
                outcome.bound = certificate.bound(search.best, network.cost_exp);
                outcome.proven_optimal = true;
            }
            return Some(());
        }
        if outcome.interrupted {
            return Some(());
        }
        if cutoff == final_cutoff {
            if let Some((_, values)) = seed {
                outcome.best = Some(values);
                outcome.bound = certificate.bound(cutoff, network.cost_exp);
                outcome.proven_optimal = true;
            }
            return Some(());
        }

        outcome.bound = certificate.bound(cutoff, network.cost_exp);
        // 界是严格上界。首轮仅含零减价；下一轮先查小于两倍最小正减价的解，
        // 两项最小正减价之和留到后续层。最后一层仍覆盖所有可能严格改善。
        let step = cutoff.max(first_positive.unwrap_or(1));
        cutoff = step.saturating_mul(2).min(final_cutoff);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simplex::Simplex;

    fn two_jobs(upper: f64) -> Problem {
        let mut p = Problem::new("时间窗临界值", Sense::Minimize);
        for i in 0..16 {
            p.add_var(&format!("a{i}"), VarType::Binary).unwrap();
        }
        let t0 = p.add_var("t0", VarType::Continuous).unwrap();
        let t1 = p.add_var("t1", VarType::Continuous).unwrap();
        p.set_var_bounds(t0, 0.1, 1.0).unwrap();
        p.set_var_bounds(t1, 0.0, upper).unwrap();
        for i in 0..4 {
            p.add_constraint(
                None,
                (0..4).map(|j| (4 * i + j, 1.0)).collect(),
                RelOp::Eq,
                1.0,
            )
            .unwrap();
            p.add_constraint(
                None,
                (0..4).map(|j| (4 * j + i, 1.0)).collect(),
                RelOp::Eq,
                1.0,
            )
            .unwrap();
        }
        p.add_constraint(
            None,
            vec![(t1, 1.0), (t0, -1.0), (11, -10.0), (10, 1.0)],
            RelOp::Ge,
            -9.5,
        )
        .unwrap();
        p.add_constraint(
            None,
            vec![(t0, 1.0), (t1, -1.0), (14, -10.0), (15, 1.0)],
            RelOp::Ge,
            -9.5,
        )
        .unwrap();
        p.set_objective(
            vec![
                (2, 2.0),
                (3, 2.0),
                (6, 2.0),
                (7, 2.0),
                (8, 3.0),
                (9, 3.0),
                (12, 3.0),
                (13, 3.0),
                (10, 10.0),
                (14, 1.0),
                (15, 10.0),
            ],
            7.0,
        )
        .unwrap();
        p
    }

    fn try_cover(p: &Problem, nodes: u64, deadline: Option<Instant>) -> RouteOutcome {
        p.validate().unwrap();
        let model = RouteModel::detect(p).unwrap();
        optimize(p, &model, Tolerances::default(), deadline, nodes, None)
    }

    #[test]
    fn exact_windows_preserve_and_exclude_adjacent_float_cutoffs() {
        // 0.1 + 0.5 的精确值严格高于 f64 字面量 0.6，下一浮点上端则可行。
        for (upper, cost) in [(0.6, 6.0), (0.6_f64.next_up(), 5.0)] {
            let p = two_jobs(upper);
            let result = try_cover(&p, 100, None);
            assert!(result.proven_optimal && !result.interrupted);
            assert_eq!(p.obj.eval(&result.best.unwrap()) + p.obj_offset, cost + 7.0);
            let bound = result.bound.unwrap();
            assert!(bound <= cost && cost - bound < 1.0e-12);
        }
    }

    #[test]
    fn permutation_keeps_proof_and_station_asymmetry_falls_back() {
        let mut p = two_jobs(0.6_f64.next_up());
        p.constraints.reverse();
        let n = p.vars.len();
        p.vars.reverse();
        for expr in p
            .constraints
            .iter_mut()
            .map(|row| &mut row.expr)
            .chain([&mut p.obj])
        {
            expr.terms.reverse();
            for (j, _) in &mut expr.terms {
                *j = n - 1 - *j;
            }
        }
        let a = try_cover(&p, 100, None);
        let b = try_cover(&p, 100, None);
        assert!(a.proven_optimal && b.proven_optimal);
        assert_eq!(a.best, b.best);
        assert_eq!(p.obj.eval(&a.best.unwrap()), 5.0);

        p.obj.add_term(n - 1 - 2, 1.0);
        let result = try_cover(&p, 100, None);
        assert!(result.bound.is_none() && result.best.is_none() && !result.proven_optimal);
        assert_eq!(result.nodes, 0);
    }

    #[test]
    fn node_and_clock_limits_do_not_claim_optimal() {
        let p = two_jobs(0.6_f64.next_up());
        let limited = try_cover(&p, 1, None);
        assert!(limited.interrupted && !limited.proven_optimal);
        assert_eq!(limited.nodes, 1);
        assert!(limited.bound.is_some_and(|b| b <= 5.0));
        let expired = try_cover(&p, 100, Some(Instant::now()));
        assert!(expired.interrupted && !expired.proven_optimal);
        assert_eq!(expired.nodes, 0);
        assert!(expired.bound.is_none());
    }

    #[test]
    fn rational_system_checks_swaps_singularities_and_overflow() {
        let (denominator, solution) =
            solve_exact_system(vec![vec![0, 1, 2], vec![2, 1, 3]], None).unwrap();
        assert_eq!(denominator, 2);
        assert_eq!(solution, vec![1, 4]);
        assert!(solve_exact_system(vec![vec![1, 1, 2], vec![1, 1, 3]], None).is_none());
        assert!(
            solve_exact_system(vec![vec![Exact::MAX, 0, 2], vec![0, Exact::MAX, 3]], None,)
                .is_none()
        );
        assert!(grid_value(1.0, -90).is_none());
    }

    fn route(jobs: u128, uses_station: bool, cost: Exact) -> Route {
        Route {
            jobs,
            first: 0,
            last: 0,
            inside: Vec::new(),
            uses_station,
            cost,
        }
    }

    #[test]
    fn symmetry_respects_prices_and_station_use_in_each_layer() {
        // 前三个任务在便宜路线中可互换；开放更多路线后，站点占用或减价各不相同。
        let mut routes = vec![route(16, false, 0)];
        for bit in [1, 2, 4] {
            routes.push(route(bit | 8, true, 0));
            routes.push(route(bit, false, 2));
        }
        routes.extend([route(17, true, 3), route(18, false, 3), route(20, true, 4)]);
        for stations in 0..=3 {
            for penalty in [0, 2] {
                let proof = Certificate {
                    dual_bound: 0,
                    denominator: 1,
                    station_penalty: penalty,
                    reduced_costs: routes.iter().map(|r| r.cost).collect(),
                };
                let mut search = CoverSearch::new(5, stations, &routes, &proof, None, 10000);
                let cheap = search.interchangeable_jobs(3).unwrap();
                assert!(cheap.iter().any(|class| class.jobs == 7));
                let full = search.interchangeable_jobs(5).unwrap();
                assert!(full.iter().all(|class| (class.jobs & 7).count_ones() <= 1));

                // 直接检查全部路线子集，独立核对严格上界、最优代价和重建的选择。
                let evaluate = |ids: &[usize]| {
                    let (mut jobs, mut used, mut cost) = (0, 0, 0);
                    for &i in ids {
                        let r = &routes[i];
                        if jobs & r.jobs != 0 {
                            return None;
                        }
                        jobs |= r.jobs;
                        used += usize::from(r.uses_station);
                        cost += r.cost;
                    }
                    if jobs != 31 || used > stations {
                        return None;
                    }
                    Some(cost + (stations - used) as Exact * penalty)
                };
                let optimum = (0..(1 << routes.len()))
                    .filter_map(|subset| evaluate(&Bits(subset).collect::<Vec<_>>()))
                    .min();
                for cutoff in 1..=10 {
                    search.run(cutoff);
                    assert!(!search.interrupted);
                    let expected = optimum.filter(|&v| v < cutoff);
                    let selected = search.selected.take();
                    assert_eq!(selected.is_some(), expected.is_some());
                    assert_eq!(selected.as_deref().and_then(evaluate), expected);
                    assert_eq!(search.best, expected.unwrap_or(cutoff));
                }
            }
        }
    }

    #[test]
    fn fractional_cover_needs_an_exact_positive_cost_search() {
        // 三条双任务路线各选一半得 LP 界 6。整数覆盖至少须有一个闲置任务，
        // 选一条双任务路线加一个闲置任务的目标为 13。
        let routes = [
            route(0b011, true, 4),
            route(0b101, true, 4),
            route(0b110, true, 4),
            route(0b001, false, 9),
            route(0b010, false, 9),
            route(0b100, false, 9),
        ];
        let lp = route_lp(&routes, 3, 2, 0);
        let mut simplex = Simplex::new(&lp, Tolerances::default()).unwrap();
        let relaxed = simplex.solve(&Bounds::from_problem(&lp));
        assert_eq!(relaxed.status, LpStatus::Optimal);
        assert_eq!(relaxed.objective, 6.0);
        let proof = Certificate::from_basis(simplex.basis_columns(), &routes, 3, 2, None).unwrap();
        assert_eq!(proof.dual_bound, 6 * proof.denominator);
        assert!(Certificate::from_basis(&[3, 4, 5, 6], &routes, 3, 2, None).is_none());

        let mut search = CoverSearch::new(3, 2, &routes, &proof, None, 100);
        search.run(7 * proof.denominator);
        assert!(search.selected.is_none() && !search.interrupted);
        search.run(7 * proof.denominator + 1);
        assert_eq!(search.selected, Some(vec![0, 5]));
        assert_eq!(search.best, 7 * proof.denominator);
    }
}
