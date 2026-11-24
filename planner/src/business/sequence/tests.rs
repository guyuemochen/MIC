use super::*;
use crate::rng::Rng;

fn state(lists: &[&[usize]]) -> State {
    State {
        values: vec![],
        lists: vec![lists.iter().map(|x| x.to_vec()).collect()],
    }
}

fn flexible(items: usize, owners: usize) -> ListSpec {
    ListSpec {
        owners,
        eligible: vec![(0..owners).collect(); items],
        durations: (0..items as i64).map(|x| x + 1).collect(),
        precedences: vec![],
    }
}

#[test]
fn directed_distance_and_load_cover_empty_ends_and_multisegment_undo() {
    let mut s = state(&[&[0, 1, 2], &[]]);
    let original = s.clone();
    // 有向、不对称，仓库到仓库的元素故意取非零。
    let matrix = [91, 2, 12, 6, 3, 92, 4, 7, 14, 5, 93, 8, 9, 10, 11, 94];
    let demands = [2, 5, 7];
    let route_costs = |s: &State| {
        s.lists[0]
            .iter()
            .map(|r| distance(r, &matrix, 3).unwrap())
            .collect::<Vec<_>>()
    };
    let route_loads = |s: &State| {
        s.lists[0]
            .iter()
            .map(|r| load(r, &demands).unwrap())
            .collect::<Vec<_>>()
    };
    assert_eq!(route_costs(&s), [23, 0]);
    let mut all = vec![];
    for (src, start, count, dst, pos) in [
        (0, 1, 2, 1, 0),
        (0, 0, 1, 1, 2),
        (1, 1, 2, 1, 0),
        (1, 0, 3, 0, 0),
    ] {
        let before = route_costs(&s);
        let loads = route_loads(&s);
        let edits = relocate(&mut s, 0, src, start, count, dst, pos).unwrap();
        let mut predicted: Vec<_> = before.into_iter().map(i128::from).collect();
        let mut predicted_loads: Vec<_> = loads.into_iter().map(i128::from).collect();
        for edit in &edits {
            predicted[edit.owner] += distance_delta(edit, &matrix, 3).unwrap();
            predicted_loads[edit.owner] += load_delta(edit, &demands).unwrap();
        }
        assert_eq!(
            predicted,
            route_costs(&s)
                .into_iter()
                .map(i128::from)
                .collect::<Vec<_>>()
        );
        assert_eq!(
            predicted_loads,
            route_loads(&s)
                .into_iter()
                .map(i128::from)
                .collect::<Vec<_>>()
        );
        all.extend(edits);
    }
    undo(&mut s, &all);
    assert_eq!(s, original);
}

#[test]
fn both_swap_orders_and_invalid_splices_preserve_state() {
    let mut s = state(&[&[0, 1, 2, 3], &[4]]);
    let before = s.clone();
    for &(a, i, b, j, ref result) in &[
        (0, 0, 0, 2, vec![vec![2, 1, 0, 3], vec![4]]),
        (0, 2, 0, 0, vec![vec![2, 1, 0, 3], vec![4]]),
        (0, 1, 0, 2, vec![vec![0, 2, 1, 3], vec![4]]),
        (0, 2, 1, 0, vec![vec![0, 1, 4, 3], vec![2]]),
    ] {
        let edits = swap(&mut s, 0, a, i, b, j).unwrap();
        assert_eq!(&s.lists[0], result);
        undo(&mut s, &edits);
        assert_eq!(s, before);
    }
    assert!(relocate(&mut s, 0, 0, 0, 3, 0, 2).is_err());
    assert!(splice(&mut s, 0, 0, 3, usize::MAX, &[]).is_err());
    assert!(swap(&mut s, 0, 0, 0, 1, 1).is_err());
    assert_eq!(s, before);
}

#[test]
fn missing_predecessor_and_rejected_cycle_keep_schedule_reusable() {
    let mut spec = flexible(4, 2);
    spec.precedences = vec![(0, 1), (2, 3)];
    let mut s = state(&[&[0, 3], &[2, 1]]);
    let mut cache = Schedule::new(&spec, &s.lists[0]).unwrap();
    assert_eq!(cache.finish, [Some(1), Some(5), Some(3), Some(7)]);
    // 先在一台机器换序；再换另一台会连接成环。
    let first = swap(&mut s, 0, 0, 0, 0, 1).unwrap();
    assert_eq!(cache.update(&spec, &first).unwrap(), vec![0, 1]);
    let intermediate = cache.finish.clone();
    let second = swap(&mut s, 0, 1, 0, 1, 1).unwrap();
    assert_eq!(cache.update(&spec, &second), Err(CYCLE.into()));
    assert!(!acyclic(&spec, &s.lists[0]).unwrap());
    assert_eq!(cache.finish, intermediate);
    undo(&mut s, &second);
    undo(&mut s, &first);
    cache
        .update(
            &spec,
            &first.iter().rev().map(Edit::inverse).collect::<Vec<_>>(),
        )
        .unwrap();
    assert_eq!(cache.finish, [Some(1), Some(5), Some(3), Some(7)]);
    let del = splice(&mut s, 0, 1, 0, 1, &[]).unwrap();
    cache.update(&spec, &[del]).unwrap();
    assert_eq!(cache.finish, [Some(1), Some(3), None, Some(5)]);
}

#[test]
fn schedule_keeps_wide_time_and_rejects_ineligible_items_without_cache_changes() {
    let mut spec = flexible(2, 2);
    spec.durations = vec![i64::MAX, 1];
    spec.eligible[0] = vec![0];
    let mut s = state(&[&[0], &[1]]);
    let mut cache = Schedule::new(&spec, &s.lists[0]).unwrap();
    let before = cache.finish.clone();
    let move_to_first = relocate(&mut s, 0, 1, 0, 1, 0, 1).unwrap();
    cache.update(&spec, &move_to_first).unwrap();
    assert_eq!(cache.finish(1), Some(i128::from(i64::MAX) + 1));
    assert_eq!(
        cache.finish,
        Schedule::new(&spec, &s.lists[0]).unwrap().finish
    );
    undo(&mut s, &move_to_first);
    cache
        .update(
            &spec,
            &move_to_first
                .iter()
                .rev()
                .map(Edit::inverse)
                .collect::<Vec<_>>(),
        )
        .unwrap();
    let invalid = relocate(&mut s, 0, 0, 0, 1, 1, 1).unwrap();
    assert!(cache.update(&spec, &invalid).unwrap_err().contains("候选"));
    undo(&mut s, &invalid);
    assert_eq!(cache.finish, before);
    assert_eq!(
        cache.links,
        Schedule::new(&spec, &s.lists[0]).unwrap().links
    );
}

#[test]
fn adjacent_edges_cover_both_boundaries_order_and_undo() {
    let mut s = state(&[&[0, 1, 2, 3], &[]]);
    for (src, start, count, dst, pos) in [
        (0, 1, 2, 1, 0),
        (1, 1, 1, 0, 1),
        (0, 0, 1, 0, 2),
        (0, 0, 3, 1, 1),
    ] {
        let previous = s.clone();
        let edits = relocate(&mut s, 0, src, start, count, dst, pos).unwrap();
        // 编号 4 对应假想仓库；首尾两侧的缺省邻居不得算作它。
        for before in 0..=4 {
            for after in 0..=4 {
                for owner in 0..2 {
                    let count = |state: &State| {
                        state.lists[0][owner]
                            .windows(2)
                            .filter(|edge| *edge == [before, after])
                            .count() as i128
                    };
                    let delta: i128 = edits
                        .iter()
                        .filter(|e| e.owner == owner)
                        .map(|e| adjacent_delta(e, before, after))
                        .sum();
                    let reversed: i128 = edits
                        .iter()
                        .rev()
                        .filter(|e| e.owner == owner)
                        .map(|e| adjacent_delta(&e.inverse(), before, after))
                        .sum();
                    assert_eq!(delta, count(&s) - count(&previous));
                    assert_eq!(reversed, -delta);
                }
            }
        }
    }
}

#[test]
fn random_walk_cache_matches_fresh_build_on_moves_rollbacks_and_cycle_rejections() {
    let mut spec = flexible(12, 3);
    spec.precedences = vec![
        (0, 1),
        (1, 2),
        (2, 3),
        (4, 5),
        (6, 7),
        (8, 9),
        (9, 10),
        (10, 11),
    ];
    let mut s = state(&[&[0, 3, 6, 9], &[1, 4, 7, 10], &[2, 5, 8, 11]]);
    let mut cache = Schedule::new(&spec, &s.lists[0]).unwrap();
    let mut rng = Rng::new(17);
    let mut rejected = 0;
    let mut accepted = 0;
    for _ in 0..4000 {
        let a = rng.gen_range(0, 3);
        let b = rng.gen_range(0, 3);
        if s.lists[0][a].is_empty() {
            continue;
        }
        let i = rng.gen_range(0, s.lists[0][a].len());
        let old = cache.finish.clone();
        let edits = if rng.gen_ratio(0.5) && !s.lists[0][b].is_empty() {
            let j = rng.gen_range(0, s.lists[0][b].len());
            swap(&mut s, 0, a, i, b, j).unwrap()
        } else {
            let count = rng.gen_range(1, s.lists[0][a].len() - i + 1);
            let len = s.lists[0][b].len() - if a == b { count } else { 0 };
            let j = rng.gen_range(0, len + 1);
            relocate(&mut s, 0, a, i, count, b, j).unwrap()
        };
        let expected = Schedule::new(&spec, &s.lists[0]);
        match (cache.update(&spec, &edits), expected) {
            (Ok(changes), Ok(fresh)) => {
                assert_eq!(cache.links, fresh.links);
                assert_eq!(cache.finish, fresh.finish);
                assert_eq!(
                    changes,
                    (0..12)
                        .filter(|&i| old[i] != fresh.finish[i])
                        .collect::<Vec<_>>()
                );
                if rng.gen_ratio(0.25) {
                    undo(&mut s, &edits);
                    cache
                        .update(
                            &spec,
                            &edits.iter().rev().map(Edit::inverse).collect::<Vec<_>>(),
                        )
                        .unwrap();
                    assert_eq!(cache.finish, old);
                }
                accepted += 1;
            }
            (Err(a), Err(b)) => {
                assert_eq!(a, CYCLE);
                assert_eq!(a, b);
                undo(&mut s, &edits);
                assert_eq!(cache.finish, old);
                rejected += 1;
            }
            other => panic!("增量与重建结果不符：{other:?}"),
        }
    }
    assert!(rejected > 100 && accepted > 100);
}
