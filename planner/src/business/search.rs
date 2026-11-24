//! 先构造初始解，再对分配和列表进行局部搜索。试探与提交共用同一条增量评分路径。

mod pairs;
use pairs::Pairs;

use super::{
    score::{self, Breakdown, Scorer},
    sequence::{self, Edit},
    Model, State,
};
use crate::{
    localsearch::{annealing_accepts, DEFAULT_SAMPLES_PER_STEP},
    rng::{derive_seed, Rng},
    ConstructionHeuristic, Metaheuristic,
};
use std::{
    collections::{BTreeSet, HashMap},
    thread,
    time::{Duration, Instant},
};

#[derive(Clone, Debug)]
pub struct Config {
    pub seed: u64,
    pub steps: usize,
    pub time_limit: Option<Duration>,
    pub threads: usize,
    pub construction: ConstructionHeuristic,
    pub algorithm: Metaheuristic,
}

impl Default for Config {
    fn default() -> Self {
        let defaults = crate::solver::SolverConfig::default();
        let local = defaults.local_search;
        let crate::localsearch::Termination::Steps(steps) = local.termination else {
            unreachable!()
        };
        Self {
            seed: local.seed,
            steps,
            time_limit: None,
            threads: 1,
            construction: defaults.construction,
            algorithm: local.metaheuristic,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StopReason {
    Steps,
    Time,
    NoMoves,
}

impl StopReason {
    pub fn status(self) -> &'static str {
        match self {
            Self::NoMoves => "Completed",
            Self::Steps | Self::Time => "LimitReached",
        }
    }
}

#[derive(Debug)]
pub struct Outcome {
    pub state: State,
    pub score: Breakdown,
    pub steps: usize,
    pub complete: bool,
    pub reason: StopReason,
}

impl Outcome {
    pub fn feasible(&self) -> bool {
        self.complete && self.score.hard == 0
    }
}

type Rank = (i64, i64);
type PlanResult<T> = Result<T, String>;
#[derive(Clone, Copy)]
struct SearchClock {
    start: Instant,
    limit: Option<Duration>,
}

impl SearchClock {
    fn expired(self) -> bool {
        self.limit
            .is_some_and(|limit| self.start.elapsed() >= limit)
    }

    /// 记录已提交状态的时间，用同一个读数判断它是否在预算内。
    fn commit_time(self) -> Option<Duration> {
        let elapsed = self.start.elapsed();
        self.limit
            .is_none_or(|limit| elapsed <= limit)
            .then_some(elapsed)
    }
}

/// 排序顺序同时用于同分时稳定选择移动。列表插入位置按删除源区段后的列表计算。
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Action {
    Assign(usize, usize),
    Exchange(usize, usize),
    AssignPair(usize, usize, usize, usize),
    ExchangePairs([usize; 2], [usize; 2]),
    Insert {
        group: usize,
        owner: usize,
        pos: usize,
        item: usize,
    },
    Relocate {
        group: usize,
        src: usize,
        start: usize,
        count: usize,
        dst: usize,
        pos: usize,
    },
    Swap {
        group: usize,
        a: usize,
        i: usize,
        b: usize,
        j: usize,
    },
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum TabuKey {
    Var(usize),
    Item(usize, usize),
}

#[derive(Default)]
struct Change {
    values: Vec<(usize, Option<i64>)>,
    edits: Vec<Edit>,
}

impl Change {
    fn keys(&self) -> impl Iterator<Item = TabuKey> + '_ {
        self.values.iter().map(|(i, _)| TabuKey::Var(*i)).chain(
            self.edits
                .iter()
                .flat_map(|e| e.added.iter().map(move |i| TabuKey::Item(e.group, *i))),
        )
    }
    fn restore(&self, state: &mut State) {
        for &(i, v) in &self.values {
            state.values[i] = v;
        }
        sequence::undo(state, &self.edits);
    }
    fn vars(&self) -> Vec<usize> {
        self.values.iter().map(|x| x.0).collect()
    }
}

struct Working<'a> {
    model: &'a Model,
    state: State,
    scorer: Scorer,
}

impl<'a> Working<'a> {
    fn new(model: &'a Model, state: State) -> PlanResult<Self> {
        let scorer = Scorer::new(model, &state)?;
        Ok(Self {
            model,
            state,
            scorer,
        })
    }

    /// 成环的候选可以跳过；模型或数值错误需要交还业务调用方。
    fn apply(&mut self, action: Action) -> PlanResult<Option<Change>> {
        let mut change = Change::default();
        match action {
            Action::Assign(i, idx) => {
                let value = self.model.domains[i][idx];
                change.values.push((i, self.state.values[i].replace(value)));
            }
            Action::Exchange(i, j) => {
                change
                    .values
                    .extend([(i, self.state.values[i]), (j, self.state.values[j])]);
                self.state.values.swap(i, j);
            }
            Action::AssignPair(a, x, b, y) => {
                for (var, idx) in [(a, x), (b, y)] {
                    let value = self.model.domains[var][idx];
                    change
                        .values
                        .push((var, self.state.values[var].replace(value)));
                }
            }
            Action::ExchangePairs([a, b], [c, d]) => {
                for (i, j) in [(a, c), (b, d)] {
                    change
                        .values
                        .extend([(i, self.state.values[i]), (j, self.state.values[j])]);
                    self.state.values.swap(i, j);
                }
            }
            Action::Insert {
                group,
                owner,
                pos,
                item,
            } => {
                change.edits.push(sequence::splice(
                    &mut self.state,
                    group,
                    owner,
                    pos,
                    0,
                    &[item],
                )?);
            }
            Action::Relocate {
                group,
                src,
                start,
                count,
                dst,
                pos,
            } => {
                change.edits =
                    sequence::relocate(&mut self.state, group, src, start, count, dst, pos)?;
            }
            Action::Swap { group, a, i, b, j } => {
                change.edits = sequence::swap(&mut self.state, group, a, i, b, j)?;
            }
        }
        match self
            .scorer
            .update(self.model, &self.state, &change.vars(), &change.edits)
        {
            Ok(()) => Ok(Some(change)),
            Err(e) if e == sequence::CYCLE => {
                change.restore(&mut self.state);
                Ok(None)
            }
            Err(e) => Err(e),
        }
    }

    fn rollback(&mut self, change: Change) -> PlanResult<()> {
        change.restore(&mut self.state);
        let inverse: Vec<_> = change.edits.iter().rev().map(Edit::inverse).collect();
        self.scorer
            .update(self.model, &self.state, &change.vars(), &inverse)
    }

    fn preview(&mut self, action: Action) -> PlanResult<Option<Rank>> {
        let Some(change) = self.apply(action)? else {
            return Ok(None);
        };
        let rank = self.scorer.score();
        self.rollback(change)?;
        Ok(Some(rank))
    }
}

fn choose(
    work: &mut Working<'_>,
    actions: impl Iterator<Item = Action>,
    first_fit: bool,
    clock: SearchClock,
) -> PlanResult<Option<Duration>> {
    let before = work.scorer.score();
    let mut best: Option<(Rank, Action)> = None;
    for action in actions {
        if clock.expired() {
            return Ok(None);
        }
        let Some(rank) = work.preview(action)? else {
            continue;
        };
        if first_fit && rank.0 <= before.0 {
            best = Some((rank, action));
            break;
        }
        if best.is_none_or(|old| (rank, action) < old) {
            best = Some((rank, action));
        }
    }
    let (_, action) = best.ok_or("无可用的构造移动")?;
    let change = work.apply(action)?.ok_or("已验证的构造移动出现环")?;
    let accepted = clock.commit_time();
    if accepted.is_none() {
        work.rollback(change)?;
    }
    Ok(accepted)
}

fn construct(
    model: &Model,
    cfg: &Config,
    pairs: &Pairs,
    clock: SearchClock,
) -> PlanResult<(State, Option<Duration>)> {
    let mut work = Working::new(model, State::empty(model))?;
    let mut last_commit = clock.commit_time();
    for (i, domain) in model.domains.iter().enumerate() {
        if let Some(j) = pairs.partner(i) {
            if j < i {
                continue;
            }
            let Some(at) = pairs::choose_pair(&mut work, i, j, clock)? else {
                return Ok((work.state, None));
            };
            last_commit = Some(at);
            continue;
        }
        let actions = (0..domain.len()).map(|idx| Action::Assign(i, idx));
        let Some(at) = choose(
            &mut work,
            actions,
            cfg.construction == ConstructionHeuristic::FirstFit,
            clock,
        )?
        else {
            return Ok((work.state, None));
        };
        last_commit = Some(at);
    }
    for (group, spec) in model.lists.iter().enumerate() {
        let n = spec.eligible.len();
        let mut indegree = vec![0; n];
        let mut successors = vec![Vec::new(); n];
        for &(a, b) in &spec.precedences {
            indegree[b] += 1;
            successors[a].push(b);
        }
        let mut ready: BTreeSet<_> = (0..n).filter(|&i| indegree[i] == 0).collect();
        while let Some(item) = ready.pop_first() {
            let targets: Vec<_> = spec.eligible[item]
                .iter()
                .map(|&owner| (owner, work.state.lists[group][owner].len()))
                .collect();
            let actions = targets.into_iter().flat_map(|(owner, len)| {
                (0..=len).map(move |pos| Action::Insert {
                    group,
                    owner,
                    pos,
                    item,
                })
            });
            let Some(at) = choose(&mut work, actions, false, clock)? else {
                return Ok((work.state, None));
            };
            last_commit = Some(at);
            for &next in &successors[item] {
                indegree[next] -= 1;
                if indegree[next] == 0 {
                    ready.insert(next);
                }
            }
        }
    }
    // 未完成的构造已从上面的分支返回。这里读取最后一次正式提交的时间，
    // 不把释放评分器、复制初始解或创建局搜线程的开销算入首次可行。
    let first_feasible = last_commit.filter(|_| work.scorer.score().0 == 0);
    Ok((work.state, first_feasible))
}

fn pick(rng: &mut Rng, len: usize) -> usize {
    rng.gen_range(0, len)
}

/// 尽量从可变的维度采样；空列表仍能作为接收区段的目的地。
struct Neighborhood {
    variables: Vec<usize>,
    groups: Vec<usize>,
}

impl Neighborhood {
    fn new(model: &Model) -> Self {
        let variables = (0..model.domains.len())
            .filter(|&i| model.domains[i].len() > 1)
            .collect();
        let groups = (0..model.lists.len())
            .filter(|&i| {
                let spec = &model.lists[i];
                let mut used = BTreeSet::new();
                spec.eligible
                    .iter()
                    .any(|o| o.len() > 1 || !used.insert(o[0]))
            })
            .collect();
        Self { variables, groups }
    }

    fn sample(&self, model: &Model, pairs: &Pairs, state: &State, rng: &mut Rng) -> Option<Action> {
        let choice = pick(rng, self.variables.len() + self.groups.len());
        if choice < self.variables.len() {
            let i = self.variables[choice];
            if let Some(action) = pairs.sample(i, model, state, rng) {
                return Some(action);
            }
            if self.variables.len() > 1 && rng.gen_ratio(0.25) {
                let j = self.variables[pick(rng, self.variables.len())];
                let (x, y) = (state.values[i]?, state.values[j]?);
                if i != j
                    && x != y
                    && model.domains[i].contains(&y)
                    && model.domains[j].contains(&x)
                {
                    return Some(Action::Exchange(i.min(j), i.max(j)));
                }
            }
            let domain = &model.domains[i];
            let idx = pick(rng, domain.len());
            return (Some(domain[idx]) != state.values[i]).then_some(Action::Assign(i, idx));
        }
        let group = self.groups[choice - self.variables.len()];
        let spec = &model.lists[group];
        let owners = &state.lists[group];
        let src = pick(rng, owners.len());
        if owners[src].is_empty() {
            return None;
        }
        let start = pick(rng, owners[src].len());
        let item = owners[src][start];
        let eligible = &spec.eligible[item];
        let dst = eligible[pick(rng, eligible.len())];
        match pick(rng, 3) {
            0 if !owners[dst].is_empty() => {
                let j = pick(rng, owners[dst].len());
                if (src, start) == (dst, j) || !spec.eligible[owners[dst][j]].contains(&src) {
                    return None;
                }
                let ((a, i), (b, j)) = if (src, start) < (dst, j) {
                    ((src, start), (dst, j))
                } else {
                    ((dst, j), (src, start))
                };
                Some(Action::Swap { group, a, i, b, j })
            }
            kind => {
                let count = if kind == 2 {
                    1 + pick(rng, owners[src].len() - start)
                } else {
                    1
                };
                if src != dst
                    && owners[src][start..start + count]
                        .iter()
                        .any(|&id| !spec.eligible[id].contains(&dst))
                {
                    return None;
                }
                let len = owners[dst].len() - if src == dst { count } else { 0 };
                let pos = pick(rng, len + 1);
                if src == dst && start == pos {
                    return None;
                }
                Some(Action::Relocate {
                    group,
                    src,
                    start,
                    count,
                    dst,
                    pos,
                })
            }
        }
    }
}

fn outcome(model: &Model, state: State, steps: usize, reason: StopReason) -> PlanResult<Outcome> {
    let score = score::full(model, &state)?;
    let complete = state.complete(model);
    Ok(Outcome {
        state,
        score,
        steps,
        reason,
        complete,
    })
}

fn improve(
    model: &Model,
    cfg: &Config,
    initial: &State,
    pairs: &Pairs,
    capture_first: bool,
    clock: SearchClock,
    seed: u64,
) -> PlanResult<(Outcome, Option<FirstFeasible>)> {
    let mut work = Working::new(model, initial.clone())?;
    let mut rng = Rng::new(seed);
    let neighborhood = Neighborhood::new(model);
    let mut best_state = initial.clone();
    let mut first_feasible = None;
    let mut best = work.scorer.score();
    let mut tabu = HashMap::<TabuKey, usize>::new();
    let mut history = match cfg.algorithm {
        Metaheuristic::LateAcceptance { history_len } => vec![best; history_len],
        _ => Vec::new(),
    };
    let mut inspect_best = best.0 == 0;
    let mut away_steps = 0;
    let mut steps = 0;
    let reason = loop {
        if clock.expired() {
            break StopReason::Time;
        }
        if steps >= cfg.steps {
            break StopReason::Steps;
        }
        if neighborhood.variables.is_empty() && neighborhood.groups.is_empty() {
            break StopReason::NoMoves;
        }
        let before = work.scorer.score();
        match cfg.algorithm {
            Metaheuristic::TabuSearch { tabu_size } => {
                // 给禁忌搜索一段离开最好解的探索期。持续偏离一轮仍未回到最好分数时，
                // 从保存的方案重新采样，沿用随机数进度，避免一直在较差的区域搜索。
                if pairs.len() > 0 && best.0 == 0 && away_steps >= pairs.len().max(tabu_size) {
                    work = Working::new(model, best_state.clone())?;
                    tabu.clear();
                    away_steps = 0;
                    inspect_best = false;
                }
                // 严格改善历史最好解的移动本来就享有禁忌豁免。优先检查容易漏采的复合键。
                let descent = if inspect_best {
                    inspect_best = false;
                    pairs.best_improvement(&mut work, clock)?
                } else {
                    None
                };
                let selected = if let Some(action) = descent {
                    Some(action)
                } else {
                    let mut winner: Option<(Rank, Action)> = None;
                    for _ in 0..DEFAULT_SAMPLES_PER_STEP {
                        if clock.expired() {
                            break;
                        }
                        let Some(action) = neighborhood.sample(model, pairs, &work.state, &mut rng)
                        else {
                            continue;
                        };
                        let Some(change) = work.apply(action)? else {
                            continue;
                        };
                        let rank = work.scorer.score();
                        let allowed = rank < best
                            || change
                                .keys()
                                .all(|key| tabu.get(&key).copied().unwrap_or(0) <= steps);
                        if allowed && winner.is_none_or(|old| (rank, action) < old) {
                            winner = Some((rank, action));
                        }
                        work.rollback(change)?;
                    }
                    winner.map(|(_, action)| action)
                };
                if let Some(winner) = selected {
                    let change = work.apply(winner)?.ok_or("已验证的搜索移动出现环")?;
                    // 同分的另一种排法可能开放新的改善移动，例如先挪动一间教室。
                    inspect_best = best.0 == 0 && work.scorer.score() == best;
                    for key in change.keys() {
                        tabu.insert(key, steps.saturating_add(tabu_size));
                    }
                }
                if work.scorer.score() > best {
                    away_steps += 1;
                } else {
                    away_steps = 0;
                }
            }
            _ => {
                if let Some(action) = neighborhood.sample(model, pairs, &work.state, &mut rng) {
                    if let Some(change) = work.apply(action)? {
                        let rank = work.scorer.score();
                        let accept = match cfg.algorithm {
                            Metaheuristic::LateAcceptance { .. } => {
                                rank <= history[steps % history.len()]
                            }
                            Metaheuristic::SimulatedAnnealing {
                                start_temp,
                                end_temp,
                            } => {
                                // 只在字典序的首个变化层计算温度差，改进始终接受。
                                let (next, old) = if rank.0 != before.0 {
                                    (rank.0, before.0)
                                } else {
                                    (rank.1, before.1)
                                };
                                let delta = i128::from(next) - i128::from(old);
                                annealing_accepts(delta, steps, start_temp, end_temp, &mut rng)
                            }
                            _ => unreachable!(),
                        };
                        if !accept {
                            work.rollback(change)?;
                        }
                    }
                }
                if !history.is_empty() {
                    let slot = steps % history.len();
                    history[slot] = work.scorer.score();
                }
            }
        }
        if work.scorer.score() < best {
            let Some(at) = clock.commit_time() else {
                break StopReason::Time;
            };
            best = work.scorer.score();
            inspect_best = best.0 == 0;
            if capture_first && best.0 == 0 && first_feasible.is_none() {
                first_feasible = Some(FirstFeasible {
                    at,
                    state: work.state.clone(),
                });
            }
            best_state = work.state.clone();
        } else if clock.expired() {
            break StopReason::Time;
        }
        steps += 1;
    };
    let result = outcome(model, best_state, steps, reason)?;
    if (result.score.hard, result.score.soft) != best {
        return Err("增量评分与最终全量评分不一致".into());
    }
    Ok((result, first_feasible))
}

/// 构造起点只做一次。每个线程有自己的步数预算，共享同一墙钟截止时间；平分取最小线程编号。
pub fn optimize(model: &Model, cfg: &Config) -> PlanResult<Outcome> {
    Ok(optimize_observed(model, cfg, false)?.0)
}

/// 首次提交的完整可行解及从求解入口开始的用时，供独立验分。
#[derive(Debug)]
pub struct FirstFeasible {
    pub at: Duration,
    pub state: State,
}

/// 自求解开始到首次提交完整且硬罚分为零的解。所有线程共用同一时钟；
/// 只统计预算内的历史最好解，未找到时返回 None。
pub fn optimize_timed(model: &Model, cfg: &Config) -> PlanResult<(Outcome, Option<FirstFeasible>)> {
    optimize_observed(model, cfg, true)
}

fn optimize_observed(
    model: &Model,
    cfg: &Config,
    capture: bool,
) -> PlanResult<(Outcome, Option<FirstFeasible>)> {
    let clock = SearchClock {
        start: Instant::now(),
        limit: cfg.time_limit,
    };
    if !(1..=1024).contains(&cfg.threads) {
        return Err("线程数须在 1 至 1024 之间".into());
    }
    match cfg.algorithm {
        Metaheuristic::TabuSearch { tabu_size: 0 } => return Err("禁忌期限须为正数".into()),
        Metaheuristic::LateAcceptance { history_len }
            if !(1..=1_000_000).contains(&history_len) =>
        {
            return Err("历史窗口须在 1 至 1000000 之间".into())
        }
        Metaheuristic::SimulatedAnnealing {
            start_temp,
            end_temp,
        } if !start_temp.is_finite()
            || !end_temp.is_finite()
            || end_temp <= 0.0
            || end_temp > start_temp =>
        {
            return Err("退火温度无效".into())
        }
        _ => (),
    }
    score::validate(model)?;
    let pairs = match cfg.construction {
        ConstructionHeuristic::BestFit => Pairs::discover(model),
        ConstructionHeuristic::FirstFit => Pairs::default(),
    };
    let (initial, at) = construct(model, cfg, &pairs, clock)?;
    let initial_feasible = at.filter(|_| capture).map(|at| FirstFeasible {
        at,
        state: initial.clone(),
    });
    let capture_first = capture && initial_feasible.is_none();
    if clock.expired() {
        return Ok((
            outcome(model, initial, 0, StopReason::Time)?,
            initial_feasible,
        ));
    }
    if !initial.complete(model) {
        return Err("构造结束时仍有未分配的元素".into());
    }
    if cfg.threads == 1 {
        let (result, first) =
            improve(model, cfg, &initial, &pairs, capture_first, clock, cfg.seed)?;
        return Ok((result, initial_feasible.or(first)));
    }
    thread::scope(|scope| {
        let mut workers = Vec::new();
        let initial = &initial;
        let pairs = &pairs;
        for idx in 0..cfg.threads {
            let seed = if idx == 0 {
                cfg.seed
            } else {
                derive_seed(cfg.seed, idx as u64)
            };
            workers.push(
                thread::Builder::new()
                    .spawn_scoped(scope, move || {
                        improve(model, cfg, initial, pairs, capture_first, clock, seed)
                    })
                    .map_err(|e| format!("无法启动搜索线程：{e}"))?,
            );
        }
        let mut winner: Option<Outcome> = None;
        let mut first_feasible = initial_feasible;
        for w in workers {
            let (result, first) = w.join().map_err(|_| "搜索线程异常结束")??;
            if let Some(event) = first {
                if first_feasible.as_ref().is_none_or(|old| event.at < old.at) {
                    first_feasible = Some(event);
                }
            }
            if winner.as_ref().is_none_or(|old| {
                (result.score.hard, result.score.soft) < (old.score.hard, old.score.soft)
            }) {
                winner = Some(result);
            }
        }
        Ok((winner.ok_or("缺少搜索结果")?, first_feasible))
    })
}

#[cfg(test)]
mod timing_tests {
    use super::*;
    use crate::business::{Expr, Rule};

    /// 第二个变量的值固定为 0。构造第一个变量时看不到冲突，局搜可以把它移到 1。
    fn conflict_model(domain: Vec<i64>) -> Model {
        let mut model = Model {
            domains: vec![domain, vec![0]],
            ..Model::default()
        };
        let a = model.expr(Expr::Var(0));
        let b = model.expr(Expr::Var(1));
        let conflict = model.expr(Expr::Eq(a, b));
        model.rules.push(Rule {
            name: "同时占用".into(),
            hard: true,
            weight: 1,
            expr: conflict,
        });
        model
    }

    #[test]
    fn only_complete_committed_states_count_as_first_feasible() {
        let cfg = Config {
            steps: 0,
            ..Config::default()
        };
        let (r, first) = optimize_timed(&conflict_model(vec![1]), &cfg).unwrap();
        assert!(r.feasible() && first.is_some());
        let (r, first) = optimize_timed(&conflict_model(vec![0, 1]), &cfg).unwrap();
        assert!(r.complete && !r.feasible() && first.is_none());
        let (r, first) = optimize_timed(
            &conflict_model(vec![1]),
            &Config {
                time_limit: Some(Duration::ZERO),
                ..cfg
            },
        )
        .unwrap();
        assert!(!r.complete && first.is_none());
        let (r, first) = optimize_timed(&conflict_model(vec![0]), &Config::default()).unwrap();
        assert!(r.complete && !r.feasible() && first.is_none());
    }

    #[test]
    fn first_feasible_may_come_from_a_search_move() {
        for threads in [1, 3] {
            let cfg = Config {
                threads,
                steps: 15,
                ..Config::default()
            };
            let started = Instant::now();
            let (r, first) = optimize_timed(&conflict_model(vec![0, 1]), &cfg).unwrap();
            assert!(r.feasible());
            let event = first.unwrap();
            assert!(event.at <= started.elapsed());
            assert_eq!(event.state.values, vec![Some(1), Some(0)]);
        }
    }
}
