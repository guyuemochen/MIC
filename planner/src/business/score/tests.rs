//! 在两个互不依赖的计算路径上比较每一次移动和撤销后的逐项分数。

use super::*;
use crate::business::Rule;
use crate::rng::Rng;

fn rule(name: &str, hard: bool, weight: i64, expr: usize) -> Rule {
    Rule {
        name: name.to_owned(),
        hard,
        weight,
        expr,
    }
}

fn assert_score(model: &Model, state: &State, scorer: &Scorer) {
    assert_eq!(scorer.scores, full(model, state).unwrap());
}

#[test]
fn unknown_values_duplicate_aggregates_and_signed_objective() {
    let mut m = Model::default();
    m.domains = vec![vec![0, 1, 2]; 3];
    let a = m.expr(Expr::Var(0));
    let b = m.expr(Expr::Var(1));
    let c = m.expr(Expr::Var(2));
    let eq = m.expr(Expr::Eq(a, b));
    let zero = m.expr(Expr::Const(0));
    let cost = m.expr(Expr::Lookup(vec![-9, 8, 17], c));
    let partial = m.expr(Expr::If(eq, cost, zero));
    let max = m.expr(Expr::Max(vec![cost, partial, cost]));
    let sum = m.expr(Expr::Sum(vec![max, cost, cost]));
    let changed = m.expr(Expr::Ne(b, c));
    m.rules.push(rule("同班", true, 3, eq));
    m.rules.push(rule("偏好", false, 2, changed));
    m.objective = Some(sum);
    let mut s = State::empty(&m);
    let mut inc = Scorer::new(&m, &s).unwrap();
    assert_eq!((0, 0), inc.score());
    let mut rng = Rng::new(854630);
    for _ in 0..800 {
        let x = rng.gen_range(0, 3);
        let y = rng.gen_range(0, 3);
        s.values[x] = Some(rng.gen_range(0, 3) as i64);
        s.values[y] = if rng.gen_range(0, 5) == 0 {
            None
        } else {
            Some(rng.gen_range(0, 3) as i64)
        };
        inc.update(&m, &s, &[x, y], &[]).unwrap();
        assert_score(&m, &s, &inc);
    }
    s.values = vec![Some(0); 3];
    inc.update(&m, &s, &[0, 1, 2], &[]).unwrap();
    assert_eq!(inc.score(), (3, -27));
}

#[test]
fn list_moves_cycles_and_undo_match_full_breakdown() {
    let mut model = Model::default();
    let n = 6;
    model.lists.push(ListSpec {
        owners: 3,
        eligible: vec![vec![0, 1, 2]; n],
        durations: vec![1, 2, 3, 7, 11, 0],
        precedences: vec![(0, 3), (1, 2), (0, 4), (2, 4), (3, 5)],
    });
    model.domains.push(vec![0, 1]);
    let mut cost = vec![];
    let mut finish = vec![];
    let matrix = (0..(n + 1) * (n + 1))
        .map(|i| (i % 19) as i64)
        .collect::<Vec<_>>();
    for owner in 0..3 {
        cost.push(model.expr(Expr::Distance {
            group: 0,
            owner,
            matrix: matrix.clone(),
        }));
        let load = model.expr(Expr::Load {
            group: 0,
            owner,
            demands: vec![2, 4, 1, 9, 5, 6],
        });
        let cap = model.expr(Expr::Const(10));
        let excess = model.expr(Expr::Sub(load, cap));
        let zero = model.expr(Expr::Const(0));
        let pen = model.expr(Expr::Max(vec![zero, excess]));
        model
            .rules
            .push(rule(&format!("载重{owner}"), true, 4, pen));
    }
    for owner in 0..3 {
        for (before, after) in [(0, 1), (1, 0), (2, 5)] {
            let adj = model.expr(Expr::Adjacent {
                group: 0,
                owner,
                before,
                after,
            });
            model.rules.push(rule(
                &format!("相邻{owner}_{before}_{after}"),
                false,
                2,
                adj,
            ));
        }
    }
    for item in 0..n {
        finish.push(model.expr(Expr::Finish { group: 0, item }));
    }
    let end = model.expr(Expr::Max(finish));
    let pref = model.expr(Expr::Var(0));
    model.rules.push(rule("偏好", false, 2, pref));
    cost.push(end);
    let obj = model.expr(Expr::Sum(cost));
    model.objective = Some(obj);
    let mut s = State::empty(&model);
    let mut inc = Scorer::new(&model, &s).unwrap();
    for i in 0..n {
        let owner = i % 3;
        let pos = s.lists[0][owner].len();
        let e = sequence::splice(&mut s, 0, owner, pos, 0, &[i]).unwrap();
        inc.update(&model, &s, &[], &[e]).unwrap();
        assert_score(&model, &s, &inc);
    }
    let mut rng = Rng::new(445954);
    let mut cyclic = 0;
    for k in 0..600 {
        let owner = rng.gen_range(0, 3);
        let dest = rng.gen_range(0, 3);
        let len = s.lists[0][owner].len();
        if len == 0 {
            continue;
        }
        let pos = rng.gen_range(0, len);
        let count = rng.gen_range(1, len - pos + 1);
        let max_dest = s.lists[0][dest].len() - if owner == dest { count } else { 0 };
        let to = rng.gen_range(0, max_dest + 1);
        let edits = sequence::relocate(&mut s, 0, owner, pos, count, dest, to).unwrap();
        match inc.update(&model, &s, &[], &edits) {
            Ok(()) => {
                assert_score(&model, &s, &inc);
                if k % 3 == 0 {
                    sequence::undo(&mut s, &edits);
                    let undo = edits.iter().rev().map(Edit::inverse).collect::<Vec<_>>();
                    inc.update(&model, &s, &[], &undo).unwrap();
                }
            }
            Err(e) => {
                assert_eq!(e, sequence::CYCLE);
                assert_eq!(full(&model, &s).unwrap_err(), e);
                sequence::undo(&mut s, &edits);
                cyclic += 1;
            }
        }
        assert_score(&model, &s, &inc);
        s.values[0] = Some((k % 2) as i64);
        inc.update(&model, &s, &[0], &[]).unwrap();
        assert_score(&model, &s, &inc);
    }
    assert!(cyclic > 10);
}

#[test]
fn oversized_owner_count_bad_edges_and_undefined_expression_are_errors() {
    let mut m = Model::default();
    m.lists.push(ListSpec {
        owners: usize::MAX,
        eligible: vec![],
        durations: vec![],
        precedences: vec![],
    });
    assert!(validate(&m).is_err());
    m.lists[0] = ListSpec {
        owners: 1,
        eligible: vec![vec![0]; 2],
        durations: vec![],
        precedences: vec![(0, 1), (1, 0)],
    };
    assert!(validate(&m).is_err());
    m.lists[0].precedences.clear();
    m.expr(Expr::Add(1, 2));
    assert!(validate(&m).is_err());
}

#[test]
fn bad_lookup_negative_penalty_and_overflow_propagate_as_errors() {
    let mut m = Model::default();
    m.domains.push(vec![0, 1]);
    let x = m.expr(Expr::Var(0));
    let c = m.expr(Expr::Lookup(vec![9], x));
    m.objective = Some(c);
    let mut s = State::empty(&m);
    let mut inc = Scorer::new(&m, &s).unwrap();
    s.values[0] = Some(1);
    assert!(inc.update(&m, &s, &[0], &[]).is_err());
    assert!(full(&m, &s).is_err());
    m.nodes[c] = Expr::Const(-1);
    m.rules.push(rule("负数", true, 1, c));
    assert!(full(&m, &s).is_err());
    m.nodes[c] = Expr::Const(i64::MAX);
    m.rules[0].weight = 2;
    assert!(full(&m, &s).is_err());
}

#[test]
fn numeric_error_leaves_cached_rule_aggregate_and_schedule_reusable() {
    let mut m = Model::default();
    m.lists.push(ListSpec {
        owners: 2,
        eligible: vec![vec![0, 1]; 2],
        durations: vec![1, 1],
        precedences: vec![],
    });
    let t = m.expr(Expr::Finish { group: 0, item: 1 });
    let load = m.expr(Expr::Load {
        group: 0,
        owner: 0,
        demands: vec![1, 1],
    });
    let twice = m.expr(Expr::Sum(vec![load, load]));
    let last = m.expr(Expr::Max(vec![t, load]));
    let lookup = m.expr(Expr::Lookup(vec![10, 20, 30], twice));
    m.rules.push(rule("人数", true, 10, load));
    let total = m.expr(Expr::Add(last, lookup));
    m.objective = Some(total);
    let mut s = State::empty(&m);
    let first = sequence::splice(&mut s, 0, 0, 0, 0, &[0]).unwrap();
    let mut inc = Scorer::new(&m, &s).unwrap();
    let e = sequence::splice(&mut s, 0, 0, 1, 0, &[1]).unwrap();
    assert!(inc.update(&m, &s, &[], std::slice::from_ref(&e)).is_err());
    sequence::undo(&mut s, std::slice::from_ref(&e));
    assert_score(&m, &s, &inc);
    let e = sequence::splice(&mut s, 0, 1, 0, 0, &[1]).unwrap();
    inc.update(&m, &s, &[], &[e]).unwrap();
    assert_score(&m, &s, &inc);
    sequence::undo(&mut s, std::slice::from_ref(&first));
    inc.update(&m, &s, &[], &[first.inverse()]).unwrap();
    assert_score(&m, &s, &inc);
}

#[test]
fn guard_masks_lookup_fault_and_unused_overflow_during_optimization() {
    use crate::business::{optimize, Config};
    let mut m = Model::default();
    m.domains.push(vec![-1, 0]);
    let x = m.expr(Expr::Var(0));
    let zero = m.expr(Expr::Const(0));
    let available = m.expr(Expr::Ge(x, zero));
    let price = m.expr(Expr::Lookup(vec![5], x));
    let obj = m.expr(Expr::If(available, price, zero));
    m.objective = Some(obj);
    let limit = m.expr(Expr::Const(i64::MAX));
    m.expr(Expr::Mul(limit, limit));
    let result = optimize(
        &m,
        &Config {
            steps: 30,
            ..Config::default()
        },
    )
    .unwrap();
    assert_eq!(result.state.values, vec![Some(-1)]);
    assert_eq!(result.score.soft, 0);
    let mut s = State::empty(&m);
    let mut inc = Scorer::new(&m, &s).unwrap();
    for (x, score) in [(-1, 0), (0, 5), (-1, 0)] {
        s.values[0] = Some(x);
        inc.update(&m, &s, &[0], &[]).unwrap();
        assert_eq!(inc.score(), (0, score));
        assert_score(&m, &s, &inc);
    }
}

#[test]
fn selected_branch_fault_is_recoverable_when_aggregates_and_condition_change() {
    let mut m = Model::default();
    m.domains = vec![vec![0, 1], vec![-1, 0]];
    let switch = m.expr(Expr::Var(0));
    let idx = m.expr(Expr::Var(1));
    let val = m.expr(Expr::Lookup(vec![3], idx));
    let sum = m.expr(Expr::Sum(vec![val, val]));
    let max = m.expr(Expr::Max(vec![sum, val]));
    let alt = m.expr(Expr::Const(9));
    let obj = m.expr(Expr::If(switch, max, alt));
    m.objective = Some(obj);
    let mut state = State {
        values: vec![Some(0), Some(-1)],
        lists: vec![],
    };
    let mut scorer = Scorer::new(&m, &state).unwrap();
    assert_eq!(scorer.score().1, 9);
    state.values[0] = Some(1);
    assert!(scorer.update(&m, &state, &[0], &[]).is_err());
    assert!(full(&m, &state).is_err());
    state.values[1] = Some(0);
    scorer.update(&m, &state, &[0, 1], &[]).unwrap();
    assert_eq!(scorer.score().1, 6);
    assert_score(&m, &state, &scorer);
}

#[test]
fn hidden_rule_arithmetic_and_list_leaf_overflows_do_not_escape_branch() {
    let mut m = Model::default();
    m.domains = vec![vec![0, 1]];
    m.lists.push(ListSpec {
        owners: 2,
        eligible: vec![vec![0, 1]; 2],
        durations: vec![i64::MAX, 1],
        precedences: vec![],
    });
    let c = m.expr(Expr::Var(0));
    let zero = m.expr(Expr::Const(0));
    let load = m.expr(Expr::Load {
        group: 0,
        owner: 0,
        demands: vec![i64::MAX, 1],
    });
    let finish = m.expr(Expr::Finish { group: 0, item: 1 });
    let one = m.expr(Expr::Const(1));
    let huge = m.expr(Expr::Const(i64::MAX));
    let overflow = m.expr(Expr::Add(huge, one));
    for (i, term) in [load, finish, overflow].into_iter().enumerate() {
        let expr = m.expr(Expr::If(c, term, zero));
        m.rules.push(rule(&format!("罚项{i}"), true, 1, expr));
    }
    let mut s = State {
        values: vec![Some(0)],
        lists: vec![vec![vec![0], vec![1]]],
    };
    let mut inc = Scorer::new(&m, &s).unwrap();
    let edits = sequence::relocate(&mut s, 0, 1, 0, 1, 0, 1).unwrap();
    inc.update(&m, &s, &[], &edits).unwrap();
    assert_score(&m, &s, &inc);
    s.values[0] = Some(1);
    assert!(inc.update(&m, &s, &[0], &[]).is_err());
    assert!(full(&m, &s).is_err());
    s.values[0] = Some(0);
    sequence::undo(&mut s, &edits);
    let inv = edits.iter().rev().map(Edit::inverse).collect::<Vec<_>>();
    inc.update(&m, &s, &[], &inv).unwrap();
    assert_score(&m, &s, &inc);
}

fn pair_terms(model: &mut Model, keys: &[usize]) -> Vec<usize> {
    let mut terms = Vec::new();
    for (i, &a) in keys.iter().enumerate() {
        for &b in &keys[..i] {
            terms.push(model.expr(Expr::Eq(a, b)));
        }
    }
    terms
}

#[test]
fn conflict_groups_match_every_partial_move_swap_and_undo() {
    let mut m = Model::default();
    m.domains = vec![vec![-2, 0, 1, 2]; 9];
    let mut keys: Vec<_> = (0..9).map(|v| m.expr(Expr::Var(v))).collect();
    let two = m.expr(Expr::Const(2));
    keys[0] = m.expr(Expr::Mul(keys[0], two));
    let mut terms = pair_terms(&mut m, &keys[..5]);
    let shared = terms[0];
    let pruned = terms[1];
    terms.extend(pair_terms(&mut m, &keys[5..]));
    let sum = m.expr(Expr::Sum(terms));
    m.rules.push(rule("各组内部冲突", true, 2, sum));
    m.rules.push(rule("同一指标带权", false, 5, sum));
    m.rules.push(rule("另用单项比较", false, 3, shared));
    m.objective = Some(sum);
    let mut state = State::empty(&m);
    let mut scorer = Scorer::new(&m, &state).unwrap();
    assert!(matches!(scorer.aggregates[sum], Some(Aggregate::Pairs(_))));
    assert!(scorer.parents.iter().flatten().any(|&(p, _)| p == shared));
    assert!(!scorer.parents.iter().flatten().any(|&(p, _)| p == pruned));
    let mut rng = Rng::new(51870);
    let domain = [None, Some(-2), Some(0), Some(1), Some(2)];
    for step in 0..400 {
        let (x, y) = (rng.gen_range(0, 9), rng.gen_range(0, 9));
        let (old_x, old_y) = (state.values[x], state.values[y]);
        if step % 3 == 0 {
            state.values.swap(x, y);
        } else {
            state.values[x] = domain[rng.gen_range(0, domain.len())];
            state.values[y] = domain[rng.gen_range(0, domain.len())];
        }
        scorer.update(&m, &state, &[x, y], &[]).unwrap();
        assert_score(&m, &state, &scorer);
        if step % 4 == 0 {
            state.values[x] = old_x;
            state.values[y] = old_y;
            scorer.update(&m, &state, &[y, x], &[]).unwrap();
            assert_score(&m, &state, &scorer);
        }
    }
}

#[test]
fn incomplete_repeated_and_self_comparisons_keep_the_original_sum() {
    let mut m = Model::default();
    m.domains = vec![vec![0, 1]; 4];
    let keys: Vec<_> = (0..4).map(|v| m.expr(Expr::Var(v))).collect();
    let edges = pair_terms(&mut m, &keys);
    let flipped = m.expr(Expr::Eq(keys[0], keys[1]));
    let same = m.expr(Expr::Eq(keys[0], keys[0]));
    for (case, replacement) in [None, Some(edges[0]), Some(flipped), Some(same)]
        .into_iter()
        .enumerate()
    {
        let mut altered = edges.clone();
        altered.pop();
        altered.extend(replacement);
        let sum = m.expr(Expr::Sum(altered));
        m.rules
            .push(rule(&format!("不规则冲突{case}"), true, 1, sum));
    }
    let mut state = State::empty(&m);
    let mut scorer = Scorer::new(&m, &state).unwrap();
    assert!(m
        .rules
        .iter()
        .all(|r| matches!(scorer.aggregates[r.expr], Some(Aggregate::Sum(_)))));
    for i in 0..4 {
        state.values[i] = Some(0);
        scorer.update(&m, &state, &[i], &[]).unwrap();
        assert_score(&m, &state, &scorer);
    }
}

#[test]
fn conflict_bucket_fault_order_and_rollback_match_original_comparisons() {
    let mut m = Model::default();
    m.domains = vec![vec![0, 1, 2]; 3];
    let switch = m.expr(Expr::Var(0));
    let x = m.expr(Expr::Var(1));
    let y = m.expr(Expr::Var(2));
    let a = m.expr(Expr::Lookup(vec![0], x));
    let b = m.expr(Expr::Lookup(vec![1, 0], y));
    let zero = m.expr(Expr::Const(0));
    let terms = pair_terms(&mut m, &[a, b, zero]);
    let sum = m.expr(Expr::Sum(terms));
    let guarded = m.expr(Expr::If(switch, sum, zero));
    m.objective = Some(guarded);
    let mut state = State::empty(&m);
    state.values = vec![Some(0), Some(1), Some(2)];
    let mut scorer = Scorer::new(&m, &state).unwrap();
    assert!(matches!(scorer.aggregates[sum], Some(Aggregate::Pairs(_))));
    // x、y 同时越界，原比较先读 b。触及多个键的失败候选必须恢复原桶和错误。
    for candidate in [
        [1, 1, 2],
        [1, 1, 0],
        [1, 0, 2],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 2],
        [1, 0, 0],
    ] {
        let before = state.values.clone();
        state.values = candidate.map(Some).to_vec();
        match full(&m, &state) {
            Ok(expected) => {
                scorer.update(&m, &state, &[0, 1, 2], &[]).unwrap();
                assert_eq!(scorer.scores, expected);
            }
            Err(expected) => {
                assert_eq!(
                    scorer.update(&m, &state, &[0, 1, 2], &[]).unwrap_err(),
                    expected
                );
                state.values = before;
                assert_score(&m, &state, &scorer);
            }
        }
    }
}
