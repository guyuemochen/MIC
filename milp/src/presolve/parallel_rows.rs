//! 对同一线性表达式的多条限制，只需要保留最紧的下界和上界。
//! 方向和右端一律按输入 f64 精确代表的有理数比较，避免除法舍入把两行并成一行。

use super::WorkRow;
use crate::model::{dyadic, gcd_u64, RelOp};
use std::cmp::Ordering;
use std::collections::BTreeMap;

type Direction = Vec<(usize, i64, i32)>;

#[derive(Clone, Copy)]
struct Level {
    numerator: i128,
    denominator: i128,
    exponent: i32,
}

// 非零系数都可唯一写成带符号的奇数乘以二的整数次幂。除去共同的奇数因子、
// 首列的符号和公共指数，相同的键当且仅当整行左端精确成比例。
fn canonical(row: &WorkRow, zero: f64) -> Option<(Direction, Level)> {
    let (first, power) = dyadic(row.terms[0].1)?;
    let mut common = 0;
    let mut key = Vec::with_capacity(row.terms.len());
    let mut nonzero_row = false;
    for &(j, a) in &row.terms {
        let (mantissa, exponent) = dyadic(a)?;
        if mantissa == 0 {
            return None;
        }
        nonzero_row |= a.abs() > zero;
        let signed = (mantissa * first.signum()) as i64;
        if common != 1 {
            common = gcd_u64(common, signed.unsigned_abs());
        }
        key.push((j, signed, exponent - power));
    }
    // R3 保留的违反零行须交给 LP，不与放大后有非零列的行合并。
    if !nonzero_row {
        return None;
    }
    if common > 1 {
        for (_, m, _) in &mut key {
            *m /= common as i64;
        }
    }
    let (rhs, rhs_power) = dyadic(row.rhs)?;
    let level = Level {
        numerator: rhs * first.signum(),
        denominator: i128::from(common),
        exponent: rhs_power - power,
    };
    Some((key, level))
}

// 归一化右端为 (numerator/denominator) × 2^exponent。交叉乘积至多 106 位，
// 指数单独处理，即使普通除法上溢、下溢或舍入到同一数也能比较。
fn compare_levels(a: Level, b: Level) -> Ordering {
    let left = a.numerator * b.denominator;
    let right = b.numerator * a.denominator;
    let signs = left.signum().cmp(&right.signum());
    if signs != Ordering::Equal || left == 0 {
        return signs;
    }
    let magnitude = |v: i128, exp: i32| {
        let m = v.unsigned_abs();
        let shift = m.leading_zeros();
        (exp - shift as i32, m << shift)
    };
    let order = magnitude(left, a.exponent).cmp(&magnitude(right, b.exponent));
    if left < 0 {
        order.reverse()
    } else {
        order
    }
}

#[derive(Default)]
struct Bounds {
    lower: Option<(usize, Level)>,
    upper: Option<(usize, Level)>,
}

impl Bounds {
    fn consider(&mut self, i: usize, level: Level, rows: &[WorkRow]) {
        let row = &rows[i];
        let (lower_op, upper_op) = if row.terms[0].1 > 0. {
            (RelOp::Ge, RelOp::Le)
        } else {
            (RelOp::Le, RelOp::Ge)
        };
        for (side, op, tighter) in [
            (&mut self.lower, lower_op, Ordering::Greater),
            (&mut self.upper, upper_op, Ordering::Less),
        ] {
            if row.op != RelOp::Eq && row.op != op {
                continue;
            }
            if side.is_none_or(|(old, previous)| {
                let order = compare_levels(level, previous);
                if order == Ordering::Equal {
                    // 同界先留同时限制另一侧的等式，再留最早的原行。
                    row.op == RelOp::Eq && rows[old].op != RelOp::Eq
                } else {
                    order == tighter
                }
            }) {
                *side = Some((i, level));
            }
        }
    }
}

pub(super) fn reduce(rows: &mut [WorkRow], zero: f64, dropped: &mut Vec<usize>) {
    // 先按容易取得的列号筛选。大部分稀疏行在这一步便是单例，不必为它们建立方向键。
    let mut entries = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let t = &row.terms;
        if row.alive && t.len() >= 2 {
            entries.push(((t.len(), t[0].0, t[t.len() / 2].0, t[t.len() - 1].0), i));
        }
    }
    entries.sort_unstable();
    let mut keep = vec![true; rows.len()];
    for candidates in entries.chunk_by(|a, b| a.0 == b.0).filter(|g| g.len() > 1) {
        let mut groups = BTreeMap::<_, Bounds>::new();
        for &(_, i) in candidates {
            let row = &rows[i];
            if let Some((key, level)) = canonical(row, zero) {
                groups.entry(key).or_default().consider(i, level, rows);
                keep[i] = false;
            }
        }
        for limits in groups.values() {
            for &(i, _) in limits.lower.iter().chain(&limits.upper) {
                keep[i] = true;
            }
        }
    }
    for (row, retain) in rows.iter_mut().zip(keep) {
        if !retain {
            row.alive = false;
            dropped.push(row.origin);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Problem, Sense, Tolerances, VarType};
    use crate::presolve::presolve;

    fn row(origin: usize, a: f64, b: f64, op: RelOp, rhs: f64) -> WorkRow {
        WorkRow {
            origin,
            terms: vec![(0, a), (1, b)],
            op,
            rhs,
            alive: true,
        }
    }

    #[test]
    fn presolve_keeps_original_extremes_for_mixed_variables() {
        let mut p = Problem::new("parallel", Sense::Maximize);
        p.add_var("x", VarType::Continuous).unwrap();
        p.add_var("z", VarType::Integer).unwrap();
        for (i, (a, b, op, rhs)) in [
            (2., -1., RelOp::Le, 24.),   // 上界 12
            (-3., 1.5, RelOp::Ge, -18.), // 上界 6
            (1., -0.5, RelOp::Ge, 1.),   // 下界 1
            (-4., 2., RelOp::Le, -16.),  // 下界 4
            (6., -3., RelOp::Le, 36.),   // 上界 6 的重复
            (2., 1., RelOp::Ge, 2.),     // 独立方向
        ]
        .into_iter()
        .enumerate()
        {
            p.add_constraint(Some(&format!("r{i}")), vec![(0, a), (1, b)], op, rhs)
                .unwrap();
        }
        let (q, map) = presolve(&p, &Tolerances::default()).unwrap();
        assert_eq!(map.dropped_rows, [0, 2, 4]);
        assert_eq!(
            q.constraints
                .iter()
                .map(|r| r.name.as_deref().unwrap())
                .collect::<Vec<_>>(),
            ["r1", "r3", "r5"]
        );
        assert_eq!(map.var_map, [0, 1]);
        for x in -12..=12 {
            for z in -8..=8 {
                let point = [x as f64 / 2.0, z as f64];
                let feasible = |p: &Problem| {
                    p.constraints.iter().all(|r| {
                        let lhs = r.expr.eval(&point);
                        match r.op {
                            RelOp::Le => lhs <= r.rhs,
                            RelOp::Ge => lhs >= r.rhs,
                            RelOp::Eq => lhs == r.rhs,
                        }
                    })
                };
                assert_eq!(feasible(&p), feasible(&q));
            }
        }
    }

    #[test]
    fn equalities_and_conflicting_bounds_stay_sound() {
        let mut rows = vec![
            row(0, 1., 2., RelOp::Le, 1.),
            row(1, -2., -4., RelOp::Le, -2.),
            row(2, 3., 6., RelOp::Eq, 3.),
            row(3, -5., -10., RelOp::Eq, -5.),
        ];
        let mut dropped = vec![];
        reduce(&mut rows, 0., &mut dropped);
        assert_eq!(dropped, [0, 1, 3]);
        rows.push(row(4, -2., -4., RelOp::Le, -3.));
        reduce(&mut rows, 0., &mut dropped);
        assert!(rows[2].alive && rows[4].alive);
        assert_eq!(dropped, [0, 1, 3]);
    }

    #[test]
    fn almost_parallel_rows_and_ambiguous_float_ratios() {
        let next = f64::from_bits(3f64.to_bits() + 1);
        let tiny = f64::from_bits(1);
        let mut rows = vec![
            row(0, 2., 3., RelOp::Le, 1.),
            row(1, 2., next, RelOp::Le, 1.),
            row(2, 2e-14, 3e-14, RelOp::Ge, 1.),
            // 如果用浮点商作为方向键，前一对同时上溢，后一对同时下溢。
            row(3, tiny, 1., RelOp::Ge, 1.),
            row(4, 2. * tiny, 3., RelOp::Ge, 1.),
            row(5, 1e300, tiny, RelOp::Ge, 1.),
            row(6, 1e300, 2. * tiny, RelOp::Ge, 1.),
        ];
        let mut dropped = vec![];
        reduce(&mut rows, 1e-12, &mut dropped);
        assert!(dropped.is_empty());

        // 原浮点除法都会得到无穷，精确比较仍留下更紧的上界。
        let mut big_rhs = vec![
            row(0, tiny, 1., RelOp::Le, 2.),
            row(1, 2. * tiny, 2., RelOp::Le, 3.),
        ];
        reduce(&mut big_rhs, 0., &mut dropped);
        assert_eq!(dropped, [0]);
        // 两个负的右端除以首列都会下溢到 -0。
        let mut small_rhs = vec![
            row(0, 1e300, 1., RelOp::Le, -tiny),
            row(1, 1e300, 1., RelOp::Le, -2. * tiny),
        ];
        dropped.clear();
        reduce(&mut small_rhs, 0., &mut dropped);
        assert_eq!(dropped, [0]);
    }
}
