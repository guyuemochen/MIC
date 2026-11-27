//! 同号容量行可排除单列已超额的整数选项；等式与同左端反向不等式可夹住私有余量。

use super::WorkRow;
use crate::model::{dyadic, RelOp, Tolerances, VarType, Variable};
use std::collections::{hash_map::DefaultHasher, BTreeMap, BTreeSet};
use std::hash::{Hash, Hasher};

// 首个系数的符号和二进制指数给整行定向、定标。尾数不做除法，正负零同键。
fn signature(mut terms: impl Iterator<Item = (usize, f64)>, rhs: f64, lower: bool) -> Option<u64> {
    let first = terms.next()?;
    let (mantissa, base) = dyadic(first.1)?;
    if mantissa == 0 {
        return None;
    }
    let signed = mantissa.signum() as i64;
    let key = |v: f64| {
        let (m, e) = dyadic(v)?;
        Some((m as i64 * signed, if m == 0 { 0 } else { e - base }))
    };
    let mut hash = DefaultHasher::new();
    (lower == (signed > 0)).hash(&mut hash);
    key(rhs)?.hash(&mut hash);
    for (j, a) in std::iter::once(first).chain(terms) {
        j.hash(&mut hash);
        key(a)?.hash(&mut hash);
    }
    Some(hash.finish())
}

// 一步整数变化须超过两行各自容许误差之和。在每行自己的单位下各留两倍余量。
// R3 可能视近零行为无约束，单侧行必须用自身的系数计算尺度。
fn separates_step(row: &WorkRow, step: f64, tol: &Tolerances) -> bool {
    let mut norm = 0.;
    // 其他列以后可能代入，留下的公共项仍不能被 R3 当作零行。
    for &(_, a) in &row.terms {
        if a.abs() <= tol.zero {
            return false;
        }
        norm += a.abs();
    }
    let lp = tol.feasibility.max(tol.zero) * norm;
    let rounded = 10. * tol.feasibility + tol.integrality * norm;
    let checked = crate::solver::FINAL_VIOLATION_REL * (1. + row.rhs.abs() + norm);
    norm.is_finite() && step.is_finite() && step.abs() > 2. * lp.max(rounded).max(checked)
}

fn matches_side(eq: &WorkRow, j: usize, a: f64, side: &WorkRow, tol: &Tolerances) -> bool {
    if eq.terms.len() != side.terms.len() + 1 || side.terms.is_empty() || side.op == RelOp::Eq {
        return false;
    }
    let terms = eq.terms.iter().copied().filter(|&(col, _)| col != j);
    let Some((_, first)) = terms.clone().next() else {
        return false;
    };
    let scale = side.terms[0].1 / first;
    let step = a * scale;
    if !matches!(dyadic(scale), Some((1 | -1, _)))
        || step / scale != a
        || (side.op == RelOp::Ge) != (step > 0.)
        || !separates_step(side, step, tol)
    {
        return false;
    }
    // 常见的同尺度匹配不做乘除。其他比例要求正反两向完全吻合，也排除次正规数舍入。
    if scale == 1.0 {
        return eq.rhs == side.rhs && terms.eq(side.terms.iter().copied());
    }
    let same = |x: f64, y: f64| x * scale == y && y / scale == x;
    same(eq.rhs, side.rhs)
        && terms
            .zip(&side.terms)
            .all(|((j, x), &(k, y))| j == k && same(x, y))
}

pub(super) fn tighten(rows: &[WorkRow], vars: &[Variable], ub: &mut [f64], tol: &Tolerances) {
    if !vars
        .iter()
        .any(|v| v.var_type != VarType::Continuous && v.lb == 0.0 && v.ub >= 1.0)
    {
        return;
    }
    // 只有等式和单侧行能够配对时，才统计私有列。
    let pair_possible = rows
        .iter()
        .any(|r| r.alive && r.op == RelOp::Eq && r.terms.len() >= 2)
        && rows
            .iter()
            .any(|r| r.alive && r.op != RelOp::Eq && !r.terms.is_empty());
    let mut occurrences = if pair_possible {
        vec![0usize; vars.len()]
    } else {
        Vec::new()
    };
    for row in rows.iter().filter(|r| r.alive) {
        if pair_possible {
            for &(j, _) in &row.terms {
                occurrences[j] += 1;
            }
        }
        if row.terms.len() < 2 {
            continue;
        }
        // 在所有变量非负的容量行中，一项取 1 就超过总容量时，该整数列只能取零。
        let sign = if row.op == RelOp::Ge || (row.op == RelOp::Eq && row.terms[0].1 < 0.0) {
            -1.0
        } else {
            1.0
        };
        let rhs = sign * row.rhs;
        if rhs < 0.0 {
            continue;
        }
        let largest = row.terms.iter().try_fold(0.0_f64, |largest, &(j, a)| {
            (vars[j].lb >= 0.0 && ub[j] >= vars[j].lb && sign * a > 0.0)
                .then(|| largest.max(a.abs()))
        });
        let Some(largest) = largest else { continue };
        if largest <= rhs {
            continue;
        }
        let margin = tol.zero.max(tol.feasibility * largest.max(rhs));
        for &(j, a) in &row.terms {
            if vars[j].var_type != VarType::Continuous
                && vars[j].lb == 0.0
                && sign * a - rhs > margin
            {
                ub[j] = 0.0;
            }
        }
    }

    // A + a*z = b 与 A ≥ b（a>0）或 A ≤ b（a<0）共同推出 z=0。
    // 仅搜索含唯一私有整数列的等式，每行只计算一次签名。散列相同后逐项复核。
    if !pair_possible {
        return;
    }
    let mut candidates = Vec::new();
    for row in rows
        .iter()
        .filter(|r| r.alive && r.op == RelOp::Eq && r.terms.len() >= 2)
    {
        let mut private = row.terms.iter().copied().filter(|&(j, a)| {
            occurrences[j] == 1
                && vars[j].var_type != VarType::Continuous
                && vars[j].lb == 0.0
                && ub[j] >= 1.0
                && a.abs() > tol.zero
        });
        if let Some((j, a)) = private.next() {
            if private.next().is_none() {
                candidates.push((row, j, a));
            }
        }
    }
    if candidates.is_empty() {
        return;
    }
    // 同支持集的不同系数与右端进入不同签名。即使人为制造散列碰撞，
    // 每个等式最多逐项复核八次；忽略多余候选只会少做一次收紧。
    const MAX_MATCHES: usize = 8;
    let lengths: BTreeSet<_> = candidates
        .iter()
        .map(|(r, _, _)| r.terms.len() - 1)
        .collect();
    let mut sides = BTreeMap::<u64, Vec<&WorkRow>>::new();
    for row in rows
        .iter()
        .filter(|r| r.alive && r.op != RelOp::Eq && lengths.contains(&r.terms.len()))
    {
        if let Some(key) = signature(row.terms.iter().copied(), row.rhs, row.op == RelOp::Ge) {
            let bucket = sides.entry(key).or_default();
            if bucket.len() < MAX_MATCHES {
                bucket.push(row);
            }
        }
    }
    for (eq, j, a) in candidates {
        let remaining = eq.terms.iter().copied().filter(|&(k, _)| k != j);
        let Some(key) = signature(remaining, eq.rhs, a > 0.) else {
            continue;
        };
        if let Some(matches) = sides.get(&key) {
            if separates_step(eq, a, tol) && matches.iter().any(|r| matches_side(eq, j, a, r, tol))
            {
                ub[j] = 0.0;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(op: RelOp, rhs: f64, terms: &[(usize, f64)]) -> WorkRow {
        WorkRow {
            origin: 0,
            terms: terms.to_vec(),
            op,
            rhs,
            alive: true,
        }
    }

    fn bounds(rows: &[WorkRow], kinds: &[VarType], lower: &[f64]) -> Vec<f64> {
        let vars: Vec<_> = kinds
            .iter()
            .zip(lower)
            .enumerate()
            .map(|(j, (&var_type, &lb))| Variable {
                name: format!("x{j}"),
                var_type,
                lb,
                ub: 1.0,
            })
            .collect();
        let mut upper = vec![1.0; vars.len()];
        tighten(rows, &vars, &mut upper, &Tolerances::default());
        upper
    }

    #[test]
    fn capacity_uses_sign_domain_and_integer_step() {
        use VarType::*;
        let terms = [(0, 4.), (1, 2.), (2, 3.), (3, 2. + 1e-11)];
        for op in [RelOp::Le, RelOp::Eq] {
            let r = row(op, 2., &terms);
            assert_eq!(
                bounds(&[r], &[Integer, Binary, Continuous, Binary], &[0.; 4]),
                [0., 1., 1., 1.]
            );
        }
        let negative: Vec<_> = terms.iter().map(|&(j, a)| (j, -a)).collect();
        assert_eq!(
            bounds(
                &[row(RelOp::Ge, -2., &negative)],
                &[Integer, Binary, Continuous, Binary],
                &[0.; 4]
            ),
            [0., 1., 1., 1.]
        );
        let r = row(RelOp::Le, 2., &terms);
        assert_eq!(
            bounds(
                &[r],
                &[Integer, Binary, Continuous, Binary],
                &[0., 0., -1., 0.]
            ),
            [1.; 4]
        );
    }

    #[test]
    fn private_slack_requires_identical_terms_rhs_and_opposite_side() {
        use VarType::*;
        let kinds = [Binary, Integer, Binary];
        for sign in [1., -1.] {
            let op = if sign > 0. { RelOp::Ge } else { RelOp::Le };
            let eq = || row(RelOp::Eq, 1., &[(0, 1.), (1, sign), (2, 1.)]);
            let side = || row(op, 1., &[(0, 1.), (2, 1.)]);
            assert_eq!(bounds(&[eq(), side()], &kinds, &[0.; 3]), [1., 0., 1.]);
            let mut almost = side();
            almost.terms[0].1 += 1e-10;
            assert_eq!(bounds(&[eq(), almost], &kinds, &[0.; 3]), [1.; 3]);
            let mut other = side();
            other.rhs += 1e-10;
            assert_eq!(bounds(&[eq(), other], &kinds, &[0.; 3]), [1.; 3]);
            let mut reversed = side();
            reversed.op = if sign > 0. { RelOp::Le } else { RelOp::Ge };
            assert_eq!(bounds(&[eq(), reversed], &kinds, &[0.; 3]), [1.; 3]);
        }
    }

    #[test]
    fn powers_of_two_may_reverse_the_side_and_leave_zero_rhs() {
        use VarType::*;
        for rhs in [-3., -0., 0., 4.] {
            for a in [-1., 1.] {
                for factor in [-8.0_f64, -0.5, 0.25, 4.] {
                    for (op, upper) in [
                        (RelOp::Ge, u8::from(a * factor < 0.) as f64),
                        (RelOp::Le, u8::from(a * factor > 0.) as f64),
                    ] {
                        let eq = row(RelOp::Eq, rhs, &[(0, 2.), (1, a), (2, -1.)]);
                        let side = row(op, rhs * factor, &[(0, 2. * factor), (2, -factor)]);
                        assert_eq!(
                            bounds(&[eq, side], &[Continuous, Integer, Continuous], &[0.; 3])[1],
                            upper
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn near_powers_rounding_or_nonproportional_terms_do_not_fix() {
        let eq = row(RelOp::Eq, 3., &[(0, 1.), (1, -2.), (2, 1.)]);
        for (factor, second, rhs) in [
            (1.5, -3., 4.5),
            (
                2.0_f64.next_up(),
                -2. * 2.0_f64.next_up(),
                3. * 2.0_f64.next_up(),
            ),
            (2., (-4.0_f64).next_up(), 6.),
            (2., -4., 6.0_f64.next_up()),
        ] {
            let side = row(RelOp::Ge, rhs, &[(0, factor), (1, second)]);
            assert!(!matches_side(&eq, 2, 1., &side, &Tolerances::default()));
        }
        let side = row(RelOp::Eq, 6., &[(0, 2.), (1, -4.)]);
        assert!(!matches_side(&eq, 2, 1., &side, &Tolerances::default()));
    }

    #[test]
    fn refuses_overflow_underflow_and_rounded_subnormals() {
        let tol = Tolerances {
            zero: 0.,
            ..Tolerances::default()
        };
        for (eq_first, side_first, a) in [
            (2.0_f64.powi(-700), 2.0_f64.powi(700), 1.),
            (2.0_f64.powi(700), 2.0_f64.powi(-700), 1.),
            (1., 2., f64::MAX),
        ] {
            let eq = row(RelOp::Eq, 0., &[(0, eq_first), (1, a)]);
            let side = row(RelOp::Ge, 0., &[(0, side_first)]);
            assert!(!matches_side(&eq, 1, a, &side, &tol));
        }
        // 3×最小次正规数乘 1/2 舍入到 2 个单位，但反向得到 4 个。
        let tiny = f64::from_bits(3);
        let eq = row(RelOp::Eq, tiny, &[(0, 2.), (1, 1.)]);
        let side = row(RelOp::Ge, tiny * 0.5, &[(0, 1.)]);
        assert_eq!(side.rhs.to_bits(), 2);
        assert!(!matches_side(&eq, 1, 1., &side, &tol));
    }

    #[test]
    fn side_near_zero_now_or_after_later_substitution_cannot_fix() {
        use VarType::*;
        let tiny = 2.0_f64.powi(-30);
        let eq = row(RelOp::Eq, 1., &[(0, 1e-5), (1, 1.)]);
        let side = row(RelOp::Ge, tiny, &[(0, tiny * 1e-5)]);
        // 侧行按 R3 可忽略，x=0、z=1 是工作模型的可行解。
        assert_eq!(bounds(&[eq, side], &[Binary; 2], &[0.; 2]), [1.; 2]);
        let eq = row(RelOp::Eq, 1e-10, &[(0, 1e-13), (1, 1e-10)]);
        let side = row(RelOp::Ge, 1e-10, &[(0, 1e-13)]);
        assert_eq!(bounds(&[eq, side], &[Binary; 2], &[0.; 2]), [1.; 2]);

        let eq = row(RelOp::Eq, 1., &[(0, 1e-13), (1, -1.), (2, 1.)]);
        let side = row(RelOp::Ge, 1., &[(0, 1e-13), (1, -1.)]);
        assert_eq!(bounds(&[eq, side], &[Binary; 3], &[0.; 3]), [1.; 3]);

        // 侧行有非零项，整数私有列的一个步长仍可能落入该侧的绝对容差。
        let eq = row(RelOp::Eq, 1., &[(0, 0.25), (1, 1.)]);
        let side = row(RelOp::Ge, tiny, &[(0, tiny * 0.25)]);
        assert_eq!(bounds(&[eq, side], &[Binary; 2], &[0.; 2]), [1.; 2]);
    }
}
