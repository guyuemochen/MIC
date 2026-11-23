//! 给结构识别器统一的二元行表示。变量与目标保持原序，通用求解仍使用输入模型。
//! 每个有限双精度系数都是奇整数乘以二的幂。包括右端在内，奇数部分的公约数
//! 乘以最小幂就是共同格距，除后可得到同一组互素整数，无须近似取整。

use crate::model::{dyadic, gcd_u64, Constraint, Problem, RelOp, Tolerances, VarType};
use crate::solver::FINAL_VIOLATION_REL;
use std::time::Instant;

fn common_unit(row: &Constraint, tol: Tolerances) -> Option<f64> {
    let mut odd = 0;
    let mut power = i32::MAX;
    for v in std::iter::once(row.rhs).chain(row.expr.terms.iter().map(|&(_, a)| a)) {
        if v == 0. || (odd == 1 && power == 0 && v.fract() == 0.) {
            continue;
        }
        let (mantissa, exponent) = dyadic(v)?;
        odd = gcd_u64(odd, mantissa.unsigned_abs() as u64);
        power = power.min(exponent);
    }
    if odd == 0 || (odd == 1 && power == 0) || !(-1022..=1023).contains(&power) {
        return None;
    }
    let unit = odd as f64 * 2f64.powi(power);
    let norm: f64 = row.expr.terms.iter().map(|&(_, a)| a.abs() / unit).sum();
    let total = norm + row.rhs.abs() / unit;
    // 所有整数活动度及右端都须能精确表示。在原行单位下，加上求和舍入误差、
    // 整数/行容差仍分得开相邻格点，才可用约分后的行证界。
    let roundoff = 4. * f64::EPSILON * (row.expr.terms.len() as f64 + 2.) * total;
    let violation = 10. * tol.feasibility / unit + tol.integrality * norm;
    let final_allowance = FINAL_VIOLATION_REL * (1. / unit + total);
    (unit.is_finite()
        && (unit * total).is_finite()
        && total < (1u64 << 48) as f64
        && violation + roundoff < 1.
        && final_allowance + roundoff < 1.)
        .then_some(unit)
}

pub(crate) fn for_recognition(
    p: &Problem,
    tol: Tolerances,
    deadline: Option<Instant>,
) -> Option<Problem> {
    if !(tol.feasibility >= 0. && tol.integrality >= 0.)
        || !p.vars.iter().all(|v| {
            v.var_type == VarType::Binary
                || (v.var_type == VarType::Integer && v.lb == 0. && v.ub == 1.)
        })
    {
        return None;
    }
    let mut changed = None;
    if p.vars.iter().any(|v| v.var_type == VarType::Integer) {
        let mut form = p.clone();
        for v in &mut form.vars {
            v.var_type = VarType::Binary;
        }
        changed = Some(form);
    }
    for (i, row) in p.constraints.iter().enumerate() {
        if i % 32 == 0 && deadline.is_some_and(|end| Instant::now() >= end) {
            return None;
        }
        let unit = common_unit(row, tol).unwrap_or(1.);
        // 所有不等式用上界表示。等式第一个非零系数取正，列序保持不变。
        let flip = row.op == RelOp::Ge
            || (row.op == RelOp::Eq && row.expr.terms.first().is_some_and(|&(_, a)| a < 0.));
        if unit == 1. && !flip {
            continue;
        }
        let new = &mut changed.get_or_insert_with(|| p.clone()).constraints[i];
        let divisor = if flip { -unit } else { unit };
        for (_, a) in &mut new.expr.terms {
            *a /= divisor;
        }
        new.rhs /= divisor;
        if new.op == RelOp::Ge {
            new.op = RelOp::Le;
        }
    }
    changed
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Sense;

    fn form(coefficients: [f64; 2], rhs: f64, op: RelOp) -> Problem {
        let mut p = Problem::new("规范行", Sense::Minimize);
        p.add_var("a", VarType::Binary).unwrap();
        p.add_var("b", VarType::Binary).unwrap();
        p.add_constraint(
            None,
            vec![(0, coefficients[0]), (1, coefficients[1])],
            op,
            rhs,
        )
        .unwrap();
        p
    }

    #[test]
    fn signed_factors_and_integer_declarations() {
        for factor in [0.25, -0.5, 3., -15., 0.125 * 27.] {
            let op = if factor > 0. { RelOp::Le } else { RelOp::Ge };
            let mut p = form([3. * factor, -2. * factor], 5. * factor, op);
            p.vars[1].var_type = VarType::Integer;
            let q = for_recognition(&p, Tolerances::default(), None).unwrap();
            assert_eq!(q.constraints[0].expr.terms, [(0, 3.), (1, -2.)]);
            assert_eq!(q.constraints[0].op, RelOp::Le);
            assert_eq!(q.constraints[0].rhs, 5.);
            assert_eq!(q.vars[1].var_type, VarType::Binary);
            assert_eq!(p.vars[1].var_type, VarType::Integer);
        }
    }

    #[test]
    fn uncertain_grid_and_continuous_box_are_kept() {
        for p in [
            form([3., -2f64.next_up()], 5., RelOp::Le),
            form([1e-9, 1e-9], 1e-9, RelOp::Le),
        ] {
            assert!(for_recognition(&p, Tolerances::default(), None).is_none());
        }
        let mut p = form([0.25, 0.25], 0.25, RelOp::Le);
        p.vars[0].var_type = VarType::Continuous;
        assert!(for_recognition(&p, Tolerances::default(), None).is_none());
    }

    #[test]
    fn zero_rhs_and_negative_equations_are_canonical() {
        let p = form([-4.5, 7.5], 0., RelOp::Eq);
        let q = for_recognition(&p, Tolerances::default(), None).unwrap();
        assert_eq!(q.constraints[0].expr.terms, [(0, 3.), (1, -5.)]);
        assert_eq!(q.constraints[0].op, RelOp::Eq);
        assert_eq!(q.constraints[0].rhs, 0.);
    }
}
