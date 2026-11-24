//! 从硬冲突规则中识别联合决策。例如两个属性组成唯一键时，逐个选值
//! 会在另一个属性尚未赋值时漏看冲突。模型本身仍按原表达式评分。

use super::{pick, Action, Change, Model, PlanResult, Rng, SearchClock, State, Working};
use crate::business::Expr;
use std::collections::BTreeSet;
use std::time::Duration;

/// 限制每次联合构造的评分次数；大值域仍走原有的逐变量构造。
const MAX_COMBINATIONS: usize = 1024;

#[derive(Default)]
pub(super) struct Pairs {
    /// 按键的高位、低位保留顺序，交换时保持属性对应。
    members: Vec<[usize; 2]>,
    /// 仅在确有复合键时分配，支持按变量注册顺序构造。
    indices: Vec<Option<usize>>,
}

/// 只识别常数进位的双属性键，允许交换加法和乘法的操作数。
fn packed(model: &Model, node: usize) -> Option<([usize; 2], i64)> {
    let nodes = &model.nodes;
    let Expr::Add(a, b) = nodes[node] else {
        return None;
    };
    let (x, y, low) = match (&nodes[a], &nodes[b]) {
        (Expr::Mul(x, y), Expr::Var(low)) | (Expr::Var(low), Expr::Mul(x, y)) => (*x, *y, *low),
        _ => return None,
    };
    let (high, radix) = match (&nodes[x], &nodes[y]) {
        (Expr::Var(i), Expr::Const(k)) | (Expr::Const(k), Expr::Var(i)) => (*i, *k),
        _ => return None,
    };
    (high != low && radix != 0).then_some(([high, low], radix))
}

impl Pairs {
    pub(super) fn len(&self) -> usize {
        self.members.len()
    }

    pub(super) fn discover(model: &Model) -> Self {
        let mut keys = BTreeSet::new();
        let mut visit = |id: usize| {
            if let Expr::Eq(a, b) = model.nodes[id] {
                if let (Some((left, k)), Some((right, l))) = (packed(model, a), packed(model, b)) {
                    if k == l && left.iter().all(|i| !right.contains(i)) {
                        keys.insert((left, k));
                        keys.insert((right, k));
                    }
                }
            }
        };
        for rule in model.rules.iter().filter(|r| r.hard && r.weight > 0) {
            if let Expr::Sum(terms) = &model.nodes[rule.expr] {
                for &id in terms {
                    visit(id);
                }
            } else {
                visit(rule.expr);
            }
        }
        if keys.is_empty() {
            return Self::default();
        }

        let mut members = BTreeSet::new();
        for (pair, radix) in keys {
            let high = &model.domains[pair[0]];
            let low = &model.domains[pair[1]];
            if high.len() < 2 || low.len() < 2 || high.len() > MAX_COMBINATIONS / low.len() {
                continue;
            }
            let span =
                i128::from(*low.iter().max().unwrap()) - i128::from(*low.iter().min().unwrap());
            if span < i128::from(radix).abs() {
                members.insert(pair);
            }
        }
        let mut counts = vec![0; model.domains.len()];
        for &[a, b] in &members {
            counts[a] += 1;
            counts[b] += 1;
        }
        // 一个属性同时出现在多个不同的键中时，没有唯一的联合决策单位。
        let members: Vec<_> = members
            .into_iter()
            .filter(|&[a, b]| counts[a] == 1 && counts[b] == 1)
            .collect();
        let mut indices = vec![None; counts.len()];
        for (id, &[a, b]) in members.iter().enumerate() {
            indices[a] = Some(id);
            indices[b] = Some(id);
        }
        Self { members, indices }
    }

    pub(super) fn partner(&self, var: usize) -> Option<usize> {
        let [a, b] = self.members[self.indices.get(var).copied().flatten()?];
        Some(if a == var { b } else { a })
    }

    fn exchange(&self, model: &Model, state: &State, i: usize, j: usize) -> Option<Action> {
        let left @ [a, b] = self.members[i.min(j)];
        let right @ [c, d] = self.members[i.max(j)];
        let compatible = |x: usize, y: usize| {
            model.domains[x].contains(&state.values[y].unwrap())
                && model.domains[y].contains(&state.values[x].unwrap())
        };
        ((state.values[a], state.values[b]) != (state.values[c], state.values[d])
            && compatible(a, c)
            && compatible(b, d))
        .then_some(Action::ExchangePairs(left, right))
    }

    /// 在可行的最好解附近检查完整复合键移动，保留严格改善最大的动作。
    /// 流式枚举不缓存邻域；选出的动作用普通提交路径再次确认时间和分数。
    pub(super) fn best_improvement(
        &self,
        work: &mut Working<'_>,
        clock: SearchClock,
    ) -> PlanResult<Option<Action>> {
        if self.members.is_empty() {
            return Ok(None);
        }
        // 大邻域可能来不及扫完，给提交已发现的候选和普通禁忌搜索留出时间。
        let scan_clock = SearchClock {
            start: clock.start,
            limit: clock.limit.map(|limit| {
                let left = limit.saturating_sub(clock.start.elapsed());
                limit - left / 2
            }),
        };
        let current = work.scorer.score();
        let mut best = None;
        let mut consider = |work: &mut Working<'_>, action| -> PlanResult<()> {
            if let Some(rank) = work.preview(action)? {
                if rank < current && best.is_none_or(|old| (rank, action) < old) {
                    best = Some((rank, action));
                }
            }
            Ok(())
        };
        'rows: for (i, &[a, b]) in self.members.iter().enumerate() {
            for j in i + 1..self.members.len() {
                if scan_clock.expired() {
                    break 'rows;
                }
                if let Some(action) = self.exchange(work.model, &work.state, i, j) {
                    consider(work, action)?;
                }
            }
            for x in 0..work.model.domains[a].len() {
                for y in 0..work.model.domains[b].len() {
                    if scan_clock.expired() {
                        break 'rows;
                    }
                    if (work.state.values[a], work.state.values[b])
                        != (
                            Some(work.model.domains[a][x]),
                            Some(work.model.domains[b][y]),
                        )
                    {
                        consider(work, Action::AssignPair(a, x, b, y))?;
                    }
                }
            }
        }
        Ok(best.map(|(_, action)| action))
    }

    /// 互换整对属性让两项决策直接交换占用的键，无需先经过一次撞键的中间状态。
    pub(super) fn sample(
        &self,
        var: usize,
        model: &Model,
        state: &State,
        rng: &mut Rng,
    ) -> Option<Action> {
        let i = self.indices.get(var).copied().flatten()?;
        if !rng.gen_ratio(0.25) {
            return None;
        }
        let [a, b] = self.members[i];
        if self.members.len() > 1 && rng.gen_ratio(0.75) {
            let other = pick(rng, self.members.len() - 1);
            let j = if other < i { other } else { other + 1 };
            if let Some(action) = self.exchange(model, state, i, j) {
                return Some(action);
            }
        }
        let x = pick(rng, model.domains[a].len());
        let y = pick(rng, model.domains[b].len());
        (state.values[a] != Some(model.domains[a][x])
            || state.values[b] != Some(model.domains[b][y]))
        .then_some(Action::AssignPair(a, x, b, y))
    }
}

/// 相邻试探使用同一个工作解：高位不变时只动低位，评分等价于从原点逐个试探。
/// 只有选定的最佳组合可以正式提交；任何一步耗尽预算都恢复两项原值。
pub(super) fn choose_pair(
    work: &mut Working<'_>,
    a: usize,
    b: usize,
    clock: SearchClock,
) -> PlanResult<Option<Duration>> {
    let original = Change {
        values: vec![(a, work.state.values[a]), (b, work.state.values[b])],
        edits: Vec::new(),
    };
    let mut best = None;
    for x in 0..work.model.domains[a].len() {
        for y in 0..work.model.domains[b].len() {
            if clock.expired() {
                if best.is_some() {
                    work.rollback(original)?;
                }
                return Ok(None);
            }
            let action = if y == 0 {
                Action::AssignPair(a, x, b, y)
            } else {
                Action::Assign(b, y)
            };
            work.apply(action)?.ok_or("变量赋值被评分器拒绝")?;
            let trial = (work.scorer.score(), x, y);
            if best.is_none_or(|old| trial < old) {
                best = Some(trial);
            }
        }
    }
    let (_, x, y) = best.ok_or("无可用的联合构造移动")?;
    work.apply(Action::AssignPair(a, x, b, y))?
        .ok_or("已验证的联合构造移动被拒绝")?;
    let accepted = clock.commit_time();
    if accepted.is_none() {
        work.rollback(original)?;
    }
    Ok(accepted)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::business::{
        optimize_timed, score,
        search::{SearchClock, Working},
        Config, Rule,
    };
    use crate::ConstructionHeuristic;
    use std::time::{Duration, Instant};

    fn add_key(model: &mut Model, a: usize, b: usize, radix: i64) -> usize {
        let a = model.expr(Expr::Var(a));
        let b = model.expr(Expr::Var(b));
        let k = model.expr(Expr::Const(radix));
        let scaled = model.expr(Expr::Mul(k, a));
        model.expr(Expr::Add(b, scaled))
    }

    fn model(radix: i64) -> Model {
        let mut model = Model {
            domains: vec![vec![0, 1]; 6],
            ..Model::default()
        };
        let keys: Vec<_> = (0..3)
            .map(|i| add_key(&mut model, 2 * i, 2 * i + 1, radix))
            .collect();
        let collisions: Vec<_> = [(0, 1), (1, 2), (0, 2)]
            .into_iter()
            .map(|(i, j)| model.expr(Expr::Eq(keys[i], keys[j])))
            .collect();
        let expr = model.expr(Expr::Sum(collisions));
        model.rules.push(Rule {
            name: "重复".into(),
            hard: true,
            weight: 1,
            expr,
        });
        let costs = [0, 2, 4].map(|i| model.expr(Expr::Var(i)));
        model.objective = Some(model.expr(Expr::Sum(costs.to_vec())));
        model
    }

    #[test]
    fn joint_selection_breaks_ties_by_domain_order_and_avoids_conflict() {
        let mut m = model(2);
        for i in [1, 3, 5] {
            m.domains[i].reverse();
        }
        let cfg = Config {
            steps: 0,
            construction: ConstructionHeuristic::BestFit,
            ..Config::default()
        };
        let (result, first) = optimize_timed(&m, &cfg).unwrap();
        assert!(result.feasible());
        assert_eq!(result.score.soft, 1);
        assert_eq!(result.state.values, [0, 1, 0, 0, 1, 1].map(Some));
        assert_eq!(first.unwrap().state, result.state);

        let (scalar, _) = optimize_timed(
            &m,
            &Config {
                construction: ConstructionHeuristic::FirstFit,
                ..cfg
            },
        )
        .unwrap();
        assert_eq!(scalar.score.hard, 1);
    }

    #[test]
    fn tabu_finds_the_best_feasible_exchange_in_one_step() {
        // 两个低代价位置已占满。第三项只有和前两项之一同时交换，才能从 9 降分。
        let mut m = model(2);
        let mut costs = Vec::new();
        for (i, price) in [(0, 2), (2, 4), (4, 9)] {
            let var = m.expr(Expr::Var(i));
            costs.push(m.expr(Expr::Lookup(vec![0, price], var)));
        }
        m.objective = Some(m.expr(Expr::Sum(costs)));
        for seed in [0, 19] {
            let cfg = Config {
                seed,
                steps: 1,
                construction: ConstructionHeuristic::BestFit,
                algorithm: crate::Metaheuristic::TabuSearch { tabu_size: 8 },
                ..Config::default()
            };
            let (result, first) = optimize_timed(&m, &cfg).unwrap();
            let initial = first.unwrap().state;
            let score = score::full(&m, &initial).unwrap();
            assert_eq!((score.hard, score.soft), (0, 9));
            assert_eq!(result.steps, 1);
            assert_eq!((result.score.hard, result.score.soft), (0, 2));
            assert_eq!(result.score, score::full(&m, &result.state).unwrap());
        }
    }

    #[test]
    fn descent_respects_domains_and_restores_its_starting_state() {
        let mut m = model(4);
        m.domains[1] = vec![2, 3];
        let mut state = State::empty(&m);
        state.values = [1, 2, 0, 0, 0, 1].map(Some).to_vec();
        let pairs = Pairs::discover(&m);
        assert!(pairs.exchange(&m, &state, 0, 1).is_none());
        let mut work = Working::new(&m, state.clone()).unwrap();
        let clock = SearchClock {
            start: Instant::now(),
            limit: Some(Duration::ZERO),
        };
        assert_eq!(pairs.best_improvement(&mut work, clock).unwrap(), None);
        assert_eq!(work.state, state);
        let action = pairs
            .best_improvement(
                &mut work,
                SearchClock {
                    limit: None,
                    ..clock
                },
            )
            .unwrap();
        // 原位置 (1,2) 可改成空闲的 (0,2)，不需要和不兼容的变量对互换。
        assert_eq!(action, Some(Action::AssignPair(0, 0, 1, 0)));
        assert_eq!(work.state, state);
        assert_eq!(work.scorer.score(), (0, 1));
    }

    #[test]
    fn ambiguous_nonunique_and_large_keys_fall_back() {
        let m = model(1);
        assert!(Pairs::discover(&m).members.is_empty());
        let mut m = model(2);
        m.domains[0] = (0..1024).collect();
        assert_eq!(Pairs::discover(&m).members, [[2, 3], [4, 5]]);
        let mut m = model(2);
        let lhs = add_key(&mut m, 0, 3, 2);
        let rhs = add_key(&mut m, 4, 5, 2);
        let eq = m.expr(Expr::Eq(lhs, rhs));
        m.rules.push(Rule {
            name: "交叠".into(),
            hard: true,
            weight: 1,
            expr: eq,
        });
        assert_eq!(Pairs::discover(&m).members, [[4, 5]]);
        m.rules.iter_mut().for_each(|r| r.hard = false);
        assert!(Pairs::discover(&m).members.is_empty());
    }

    #[test]
    fn timed_out_pair_leaves_previous_assignments_intact() {
        let m = model(2);
        let mut work = Working::new(&m, State::empty(&m)).unwrap();
        work.apply(Action::AssignPair(0, 0, 1, 0)).unwrap().unwrap();
        let previous = work.state.clone();
        let prev_score = work.scorer.score();
        let clock = SearchClock {
            start: Instant::now(),
            limit: Some(Duration::ZERO),
        };
        assert!(choose_pair(&mut work, 2, 3, clock).unwrap().is_none());
        assert_eq!(work.state, previous);
        assert_eq!(work.scorer.score(), prev_score);
    }

    #[test]
    fn successive_trials_match_separate_previews_including_late_errors() {
        for costs in [vec![7, -2, 6, -2], vec![0, 1, 2]] {
            let mut m = model(2);
            let key = add_key(&mut m, 2, 3, 2);
            let lookup = m.expr(Expr::Lookup(costs, key));
            m.objective = Some(lookup);
            let mut initial = State::empty(&m);
            initial.values[0] = Some(0);
            initial.values[1] = Some(1);
            let mut separate = Working::new(&m, initial.clone()).unwrap();
            let mut successive = Working::new(&m, initial).unwrap();
            let clock = SearchClock {
                start: Instant::now(),
                limit: None,
            };
            let actions = (0..2).flat_map(|x| (0..2).map(move |y| Action::AssignPair(2, x, 3, y)));
            let expected = crate::business::search::choose(&mut separate, actions, false, clock);
            let actual = choose_pair(&mut successive, 2, 3, clock);
            match (expected, actual) {
                (Ok(Some(_)), Ok(Some(_))) => {
                    assert_eq!(separate.state, successive.state);
                    assert_eq!(separate.scorer.score(), successive.scorer.score());
                }
                (Err(a), Err(b)) => {
                    assert_eq!(a, b);
                    assert!(a.contains("查表下标 3"));
                }
                result => panic!("连续试探与独立试探结果不符：{result:?}"),
            }
        }
    }

    #[test]
    fn paired_changes_and_undo_agree_with_full_score() {
        let m = model(2);
        let mut work = Working::new(&m, State::empty(&m)).unwrap();
        for (a, b, x, y) in [(0, 1, 1, 1), (2, 3, 1, 0), (4, 5, 0, 0)] {
            work.apply(Action::AssignPair(a, x, b, y)).unwrap().unwrap();
        }
        for action in [
            Action::AssignPair(0, 0, 1, 0),
            Action::ExchangePairs([0, 1], [4, 5]),
        ] {
            let state = work.state.clone();
            let before = work.scorer.score();
            let change = work.apply(action).unwrap().unwrap();
            let full = score::full(&m, &work.state).unwrap();
            assert_eq!(work.scorer.score(), (full.hard, full.soft));
            work.rollback(change).unwrap();
            assert_eq!(work.state, state);
            assert_eq!(work.scorer.score(), before);
        }
    }
}
