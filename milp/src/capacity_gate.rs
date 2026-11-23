//! 一条容量约束由单个二元开关控制，两侧同时变号时语义相同。
use crate::model::{Constraint, Problem, RelOp, VarType, Variable};

pub(crate) fn is_binary(v: &Variable) -> bool {
    v.var_type == VarType::Binary || (v.var_type == VarType::Integer && v.lb == 0.0 && v.ub == 1.0)
}

/// 返回开关列、满容量及启用值。允许流量项按任意正系数加权。
fn gate(p: &Problem, row: &Constraint) -> Option<(usize, f64, bool)> {
    if row.op == RelOp::Eq || row.expr.terms.len() < 2 {
        return None;
    }
    let sign = if row.op == RelOp::Le { 1.0 } else { -1.0 };
    let mut binary = None;
    for &(j, a) in &row.expr.terms {
        let v = &p.vars[j];
        if is_binary(v) && binary.is_none() {
            binary = Some((j, a * sign));
        } else if v.var_type != VarType::Continuous || a * sign <= 0.0 {
            return None;
        }
    }
    let (j, a) = binary?;
    if a < 0.0 && row.rhs == 0.0 {
        Some((j, -a, true))
    } else if a > 0.0 && row.rhs * sign == a {
        Some((j, a, false))
    } else {
        None
    }
}

pub(crate) fn at_one(p: &Problem, row: &Constraint) -> Option<(usize, f64)> {
    gate(p, row).and_then(|(j, capacity, enabled)| enabled.then_some((j, capacity)))
}

/// 只固定不计费、在所有行中都只负责放宽容量的自由开关。改变它不影响
/// 其他变量或目标；若还受预算、互斥、等式、反向开关约束，就保留给搜索。
pub(crate) fn free_switches(p: &Problem) -> Vec<(usize, bool)> {
    if !p.vars.iter().any(|v| v.var_type == VarType::Continuous) {
        return Vec::new();
    }
    let mut possible: Vec<_> = p
        .vars
        .iter()
        .map(|v| is_binary(v) && v.lb == 0.0 && v.ub == 1.0)
        .collect();
    for &(j, _) in &p.obj.terms {
        possible[j] = false;
    }
    if !possible.contains(&true) {
        return Vec::new();
    }
    let mut enabled = vec![None; p.vars.len()];
    for row in &p.constraints {
        if !row.expr.terms.iter().any(|&(j, _)| possible[j]) {
            continue;
        }
        if let Some((j, _, at_one)) = gate(p, row) {
            if enabled[j].is_some_and(|old| old != at_one) {
                possible[j] = false;
            }
            enabled[j] = Some(at_one);
        } else {
            for &(j, _) in &row.expr.terms {
                possible[j] = false;
            }
        }
    }
    enabled
        .into_iter()
        .enumerate()
        .filter_map(|(j, value)| value.filter(|_| possible[j]).map(|v| (j, v)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Sense, Tolerances};

    #[test]
    fn only_free_consistent_switches_are_fixed_and_restored() {
        let mut p = Problem::new("容量开关", Sense::Maximize);
        let x = p.add_var("x", VarType::Continuous).unwrap();
        let switches: Vec<_> = (0..6)
            .map(|i| p.add_var(&format!("y{i}"), VarType::Binary).unwrap())
            .collect();
        let [open, closed, mixed, shared, paid, fixed] = switches[..] else {
            unreachable!()
        };
        for j in [open, mixed, shared, paid, fixed] {
            p.add_constraint(None, vec![(x, -0.5), (j, 2.0)], RelOp::Ge, 0.0)
                .unwrap();
        }
        for j in [closed, mixed] {
            p.add_constraint(None, vec![(x, 0.25), (j, 1.0)], RelOp::Le, 1.0)
                .unwrap();
        }
        p.add_constraint(None, vec![(shared, 1.0), (paid, 1.0)], RelOp::Le, 1.0)
            .unwrap();
        p.set_var_bounds(fixed, 0.0, 0.0).unwrap();
        p.set_objective(vec![(x, 2.0), (paid, -1.0)], 1.0).unwrap();
        assert_eq!(free_switches(&p), vec![(open, true), (closed, false)]);
        let (work, map) = crate::presolve::presolve(&p, &Tolerances::default()).unwrap();
        assert!(map.fixed.contains(&(open, 1.0)) && map.fixed.contains(&(closed, 0.0)));
        let lifted = crate::presolve::postsolve(
            &map,
            &vec![0.0; work.vars.len()],
            p.vars.len(),
            &Tolerances::default(),
        );
        assert_eq!((lifted[open], lifted[closed]), (1.0, 0.0));
    }
}
