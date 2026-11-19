//! 预处理中的短行活动度推断，依赖已认证的二元互斥组。

use super::{next_down, next_up, outward_product, WorkRow};
use crate::binary_propagation::{certifies_binary_mutex, standard_bounds};
use crate::model::{RelOp, Tolerances, VarType};

/// 只检查包含两个自由二元变量的短行。推断都基于本轮同一变量盒，
/// 行序不影响收集结果。证据来自当前存活的工作约束和变量界。
pub(super) fn zero_fixes(
    rows: &[WorkRow],
    types: &[VarType],
    alive: &[bool],
    lb: &[f64],
    ub: &[f64],
    tol: Tolerances,
) -> Vec<bool> {
    // 多于 8 项先留给原规则。该上限限制单条混合行试算成本。
    const MAX_TERMS: usize = 8;
    let standard_binary = |j: usize| {
        types.get(j) == Some(&VarType::Binary) && alive[j] && standard_bounds(lb[j], ub[j])
    };
    let free_binary = |j| standard_binary(j) && lb[j] == 0.0 && ub[j] == 1.0;
    let mut groups = Vec::<Vec<usize>>::new();
    for (i, row) in rows.iter().enumerate().filter(|(_, row)| row.alive) {
        if !certifies_binary_mutex(&row.terms, row.op, row.rhs, tol, standard_binary) {
            continue;
        }
        if groups.is_empty() {
            groups.resize_with(types.len(), Vec::new);
        }
        for &(j, _) in &row.terms {
            if free_binary(j) {
                groups[j].push(i);
            }
        }
    }
    if groups.is_empty() {
        return Vec::new();
    }
    let mut fixed = vec![false; types.len()];
    for row in rows
        .iter()
        .filter(|row| row.alive && (2..=MAX_TERMS).contains(&row.terms.len()))
    {
        let mut candidates = row.terms.iter().copied().filter(|&(j, _)| free_binary(j));
        let (Some(first), Some(second), None) =
            (candidates.next(), candidates.next(), candidates.next())
        else {
            continue;
        };
        let (a, b) = (&groups[first.0], &groups[second.0]);
        let (shorter, longer) = if a.len() <= b.len() { (a, b) } else { (b, a) };
        if !shorter
            .iter()
            .any(|group| longer.binary_search(group).is_ok())
        {
            continue;
        }
        let Some((lower, upper, margin)) = mutex_row_activity(row, lb, ub, [first, second], tol)
        else {
            continue;
        };
        let rhs_lower = next_down(row.rhs - margin);
        let rhs_upper = next_up(row.rhs + margin);
        if !rhs_lower.is_finite() || !rhs_upper.is_finite() {
            continue;
        }
        for (j, coefficient) in [first, second] {
            let assumed_lower = next_down(lower + coefficient);
            let assumed_upper = next_up(upper + coefficient);
            if !assumed_lower.is_finite() || !assumed_upper.is_finite() {
                continue;
            }
            let violates_le = row.op != RelOp::Ge && assumed_lower > rhs_upper;
            let violates_ge = row.op != RelOp::Le && assumed_upper < rhs_lower;
            fixed[j] |= violates_le || violates_ge;
        }
    }
    fixed
}

/// 两个候选列之外的行活动度区间，以及按行系数尺度放宽的数值余量。
fn mutex_row_activity(
    row: &WorkRow,
    lb: &[f64],
    ub: &[f64],
    pair: [(usize, f64); 2],
    tol: Tolerances,
) -> Option<(f64, f64, f64)> {
    if row.terms.windows(2).any(|terms| terms[0].0 >= terms[1].0) {
        return None;
    }
    let (mut lower, mut upper, mut other_norm) = (0.0, 0.0, 0.0);
    for &(j, coefficient) in &row.terms {
        let (&lo, &hi) = (lb.get(j)?, ub.get(j)?);
        if !coefficient.is_finite()
            || coefficient == 0.0
            || !lo.is_finite()
            || !hi.is_finite()
            || lo > hi
        {
            return None;
        }
        if j == pair[0].0 || j == pair[1].0 {
            continue;
        }
        other_norm = next_up(other_norm + coefficient.abs());
        let (min_bound, max_bound) = if coefficient > 0.0 {
            (lo, hi)
        } else {
            (hi, lo)
        };
        lower = next_down(lower + outward_product(coefficient, min_bound)?.0);
        upper = next_up(upper + outward_product(coefficient, max_bound)?.1);
    }
    let norm = next_up(other_norm + next_up(pair[0].1.abs() + pair[1].1.abs()));
    let feasibility = next_up(10.0 * tol.feasibility);
    let normalized_lp = next_up(feasibility.max(tol.zero) * norm);
    // 按工作行尺度，同时覆盖取整复核、终检形式和其它列贴界误差。
    let rounded_row = next_up(feasibility + next_up(tol.integrality * norm));
    let final_scale = next_up(next_up(1.0 + row.rhs.abs()) + norm);
    let final_row = next_up(crate::solver::FINAL_VIOLATION_REL * final_scale);
    let other_bounds = next_up(next_up(feasibility + tol.integrality) * other_norm);
    let margin = next_up(normalized_lp.max(rounded_row).max(final_row) + other_bounds);
    if lower.is_finite() && upper.is_finite() && margin.is_finite() {
        Some((lower, upper, margin))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::super::presolve;
    use super::*;
    use crate::model::{Problem, Sense, SolveStatus};
    use crate::solver::solve_default;

    fn reference_problem(group_op: RelOp, short_op: RelOp) -> Problem {
        let mut p = Problem::new("二元互斥", Sense::Maximize);
        p.add_var("a", VarType::Binary).unwrap();
        p.add_var("b", VarType::Binary).unwrap();
        p.add_var("t", VarType::Continuous).unwrap();
        p.set_var_bounds(2, -2.0, f64::INFINITY).unwrap();
        p.set_objective(vec![(0, 2.0), (1, 7.0)], 1.0).unwrap();
        let group_sign = if group_op == RelOp::Ge { -1.0 } else { 1.0 };
        p.add_constraint(
            None,
            vec![(0, group_sign), (1, group_sign)],
            group_op,
            group_sign,
        )
        .unwrap();
        let sign = if short_op == RelOp::Le { -1.0 } else { 1.0 };
        p.add_constraint(
            None,
            vec![(0, 2.0 * sign), (1, -3.0 * sign), (2, sign)],
            short_op,
            0.0,
        )
        .unwrap();
        p.add_constraint(None, vec![(2, 1.0)], RelOp::Le, 1.0)
            .unwrap();
        p
    }

    #[test]
    fn mutex_short_row_uses_r4_bounds_and_restores_original_values() {
        for group_op in [RelOp::Le, RelOp::Ge, RelOp::Eq] {
            for short_op in [RelOp::Le, RelOp::Ge, RelOp::Eq] {
                let p = reference_problem(group_op, short_op);
                let (reduced, map) = presolve(&p, &Tolerances::default()).unwrap();
                assert!(map.fixed.contains(&(1, 0.0)), "{group_op:?} {short_op:?}");
                assert!(reduced.vars.iter().all(|v| v.name != "b"));
                assert_eq!(
                    reduced.obj_offset,
                    if group_op == RelOp::Eq { 3.0 } else { 1.0 }
                );
                let result = solve_default(&p).unwrap();
                assert_eq!(result.status, SolveStatus::Optimal);
                assert_eq!(result.objective, Some(3.0));
                let values = result.values.unwrap();
                assert_eq!(&values[..2], &[1.0, 0.0]);
                if short_op == RelOp::Eq {
                    assert_eq!(values[2], -2.0);
                }
            }
        }
    }

    #[test]
    fn mutex_short_row_rebuilds_after_batched_removal_and_preserves_other_fixes() {
        let mut p = reference_problem(RelOp::Le, RelOp::Ge);
        let c = p.add_var("c", VarType::Binary).unwrap();
        let d = p.add_var("d", VarType::Binary).unwrap();
        let e = p.add_var("e", VarType::Binary).unwrap();
        p.add_constraint(None, vec![(c, 1.0), (d, 2.0), (e, 1.0)], RelOp::Le, 2.0)
            .unwrap();
        p.add_constraint(None, vec![(e, 1.0)], RelOp::Le, 0.0)
            .unwrap();
        // 第一次短行检查时 b/c/d 仍有三列自由二元。批量移除 b 后，
        // 下一轮才能用 c 与 d 的互斥，以及经 R4 收紧的 t，固定 d。
        p.add_constraint(
            None,
            vec![(1, 2.0), (2, 1.0), (c, 2.0), (d, -3.0)],
            RelOp::Ge,
            0.0,
        )
        .unwrap();
        let (_, map) = presolve(&p, &Tolerances::default()).unwrap();
        assert_eq!(map.fixed, vec![(1, 0.0), (e, 0.0), (d, 0.0)]);
        let result = solve_default(&p).unwrap();
        assert_eq!(result.status, SolveStatus::Optimal);
        assert_eq!(result.objective, Some(3.0));
        for j in [1, d, e] {
            assert_eq!(result.values.as_ref().unwrap()[j], 0.0);
        }
    }

    #[test]
    fn mutex_short_row_leaves_a_contradictory_empty_row_for_the_solver() {
        let mut p = Problem::new("互斥后空行矛盾", Sense::Minimize);
        p.add_var("a", VarType::Binary).unwrap();
        p.add_var("b", VarType::Binary).unwrap();
        let terms = vec![(0, 1.0), (1, 1.0)];
        p.add_constraint(None, terms.clone(), RelOp::Eq, 1.0)
            .unwrap();
        p.add_constraint(None, terms, RelOp::Le, 0.5).unwrap();

        let (reduced, map) = presolve(&p, &Tolerances::default()).unwrap();
        assert_eq!(map.fixed, vec![(0, 0.0), (1, 0.0)]);
        assert!(reduced.vars.is_empty());
        assert_eq!(reduced.constraints.len(), 1);
        assert!(reduced.constraints[0].expr.terms.is_empty());
        assert_eq!(reduced.constraints[0].op, RelOp::Eq);
        assert_eq!(reduced.constraints[0].rhs, 1.0);
        let result = solve_default(&p).unwrap();
        assert_eq!(result.status, SolveStatus::Infeasible);
    }

    #[test]
    fn mutex_short_row_rejects_uncertified_groups_and_uncertain_activity() {
        let row = |terms, op, rhs| WorkRow {
            origin: 0,
            terms,
            op,
            rhs,
            alive: true,
        };
        let detect = |group, short, types: &[VarType], lb: &[f64], ub: &[f64]| {
            zero_fixes(
                &[group, short],
                types,
                &[true; 3],
                lb,
                ub,
                Tolerances::default(),
            )
            .get(1)
            .copied()
            .unwrap_or(false)
        };
        let types = [VarType::Binary, VarType::Binary, VarType::Continuous];
        let group = || row(vec![(0, 1.0), (1, 1.0)], RelOp::Eq, 1.0);
        let short = || row(vec![(0, 1.0), (1, -1.0), (2, 1.0)], RelOp::Ge, 0.0);
        let lb = [0.0; 3];
        let ub = [1.0, 1.0, 0.5];
        assert!(detect(group(), short(), &types, &lb, &ub));
        assert!(!detect(
            group(),
            short(),
            &types,
            &lb,
            &[1.0, 1.0, 1.0 - 1e-8]
        ));
        // 全行系数范数决定取整复核容差。即便只固定二元列，也不能略去
        // 大连续系数带来的余量。假设 b=1 时最大活动度是 99999。
        let big_row = |rhs| row(vec![(0, 1.0), (1, -1.0), (2, 100_000.0)], RelOp::Ge, rhs);
        let big_ub = [1.0; 3];
        assert!(!detect(group(), big_row(99_999.01), &types, &lb, &big_ub));
        assert!(detect(group(), big_row(100_000.0), &types, &lb, &big_ub));
        assert!(!detect(
            group(),
            short(),
            &types,
            &lb,
            &[1.0, 1.0, f64::INFINITY]
        ));
        assert!(!detect(group(), short(), &types, &[0.0, 0.0, 1.0], &ub));
        assert!(!detect(
            group(),
            short(),
            &[VarType::Integer, types[1], types[2]],
            &lb,
            &ub
        ));
        assert!(!detect(
            row(vec![(0, 1.0), (1, 1.0), (2, 1.0)], RelOp::Eq, 1.0),
            short(),
            &types,
            &lb,
            &ub
        ));
        assert!(!detect(
            row(vec![(0, 1.0), (1, 1.0)], RelOp::Le, 2.0),
            short(),
            &types,
            &lb,
            &ub
        ));
        assert!(!detect(
            row(vec![(0, 1.0), (1, 1.0), (1, 1.0)], RelOp::Eq, 1.0),
            short(),
            &types,
            &lb,
            &ub
        ));
        assert!(!detect(
            group(),
            row(vec![(1, 1.0), (1, -1.0), (2, 1.0)], RelOp::Ge, 2.0),
            &types,
            &lb,
            &ub
        ));
        assert!(!detect(
            group(),
            row(vec![(0, 1.0), (1, f64::NAN), (2, 1.0)], RelOp::Ge, 0.0),
            &types,
            &lb,
            &ub
        ));
        assert!(!detect(
            group(),
            row(vec![(0, 1.0), (1, -1.0), (2, f64::MAX)], RelOp::Ge, 0.0),
            &types,
            &lb,
            &[1.0, 1.0, 2.0]
        ));
    }

    #[test]
    fn mutex_short_row_fixes_only_binary_values_excluded_by_enumeration() {
        for seed in 0..60_i32 {
            let mut p = Problem::new("枚举验证", Sense::Maximize);
            for j in 0..4 {
                p.add_var(&format!("b{j}"), VarType::Binary).unwrap();
            }
            p.add_var("z", VarType::Integer).unwrap();
            p.set_var_bounds(4, -2.0, 2.0).unwrap();
            p.set_objective((0..5).map(|j| (j, (2 * j + 1) as f64)).collect(), 0.0)
                .unwrap();
            for j in [0, 2] {
                p.add_constraint(None, vec![(j, 1.0), (j + 1, 1.0)], RelOp::Le, 1.0)
                    .unwrap();
                let sign = if (seed + j as i32) % 2 == 0 {
                    1.0
                } else {
                    -1.0
                };
                let op = [RelOp::Le, RelOp::Eq, RelOp::Ge][(seed as usize + j) % 3];
                let rhs = (seed % 11 - 5) as f64;
                p.add_constraint(
                    None,
                    vec![
                        (j, (seed % 4 + 1) as f64),
                        (j + 1, -(seed % 5 + 1) as f64),
                        (4, sign),
                    ],
                    op,
                    rhs,
                )
                .unwrap();
            }
            let (_, map) = presolve(&p, &Tolerances::default()).unwrap();
            let mut best = None::<f64>;
            for mask in 0..16 {
                for z in -2..=2 {
                    let values: Vec<_> = (0..4)
                        .map(|j| ((mask >> j) & 1) as f64)
                        .chain([z as f64])
                        .collect();
                    let feasible = p.constraints.iter().all(|row| {
                        let lhs: f64 = row.expr.terms.iter().map(|&(j, a)| a * values[j]).sum();
                        match row.op {
                            RelOp::Le => lhs <= row.rhs,
                            RelOp::Ge => lhs >= row.rhs,
                            RelOp::Eq => lhs == row.rhs,
                        }
                    });
                    if feasible {
                        for &(j, value) in &map.fixed {
                            assert_eq!(values[j], value, "种子 {seed} 可行点 {values:?}");
                        }
                        let obj: f64 = p.obj.terms.iter().map(|&(j, a)| a * values[j]).sum();
                        best = Some(best.map_or(obj, |current| current.max(obj)));
                    }
                }
            }
            let result = solve_default(&p).unwrap();
            assert_eq!(result.objective, best, "种子 {seed}");
            let expected = if best.is_some() {
                SolveStatus::Optimal
            } else {
                SolveStatus::Infeasible
            };
            assert_eq!(result.status, expected, "种子 {seed}");
        }
    }
}
