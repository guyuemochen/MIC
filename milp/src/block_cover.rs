//! 识别带启用费用和共享替代项的二元覆盖模型，供自动分支规则选择。
//! 一个共享项替代多个局部项时，选它会同时改变多个需求行；只看 LP
//! 当前值离整数的距离，很难反映这种影响。伪代价分支用子问题的界变化估价。

use crate::model::{Constraint, Problem, RelOp, Sense, VarType};
use std::collections::HashMap;

/// 统一写成小于等于；等式不参与模式识别。
fn le_sign(row: &Constraint) -> Option<f64> {
    match row.op {
        RelOp::Le => Some(1.0),
        RelOp::Ge => Some(-1.0),
        RelOp::Eq => None,
    }
}

/// Σ aᵢ xᵢ ≤ M y，aᵢ、M 均正。局部项 xᵢ 只能在 y 开启后使用。
fn activation(row: &Constraint) -> Option<usize> {
    if row.rhs != 0.0 || row.expr.terms.len() < 3 {
        return None;
    }
    let sign = le_sign(row)?;
    let mut negative = row.expr.terms.iter().filter(|&&(_, a)| sign * a < 0.0);
    let &(gate, _) = negative.next()?;
    negative.next().is_none().then_some(gate)
}

pub(crate) fn prefers_pseudocost(problem: &Problem) -> bool {
    if problem.sense != Sense::Minimize
        || !problem.constraints.iter().any(|r| activation(r).is_some())
        || problem
            .vars
            .iter()
            .any(|v| v.var_type != VarType::Binary || v.lb != 0.0 || v.ub != 1.0)
        || problem.obj.terms.iter().any(|&(_, c)| c < 0.0)
    {
        return false;
    }
    let n = problem.vars.len();
    let mut price = vec![0.0; n];
    for &(j, c) in &problem.obj.terms {
        let Some(value) = price.get_mut(j) else {
            return false;
        };
        *value = c;
    }
    let mut gates = vec![None; n];
    for row in &problem.constraints {
        let Some(gate) = activation(row) else {
            continue;
        };
        if !price.get(gate).is_some_and(|&c| c > 0.0) {
            continue;
        }
        for &(j, _) in &row.expr.terms {
            if j != gate {
                let Some(entry) = gates.get_mut(j) else {
                    return false;
                };
                if entry.replace(gate).is_some_and(|other| other != gate) {
                    return false; // 本检测器只识别单一启用项。
                }
            }
        }
    }

    // 互斥行 xᵢ + z ≤ 1：z 是同一启用组之外、另付费用的替代项。
    let mut alternatives = vec![Vec::new(); n];
    for row in &problem.constraints {
        let [(i, a), (j, b)] = row.expr.terms.as_slice() else {
            continue;
        };
        if *a != row.rhs || *b != row.rhs || !le_sign(row).is_some_and(|s| s * row.rhs > 0.0) {
            continue;
        }
        for (&local, &shared) in [(i, j), (j, i)] {
            let (Some(&Some(gate)), Some(&None)) = (gates.get(local), gates.get(shared)) else {
                continue;
            };
            if gate != shared && price[shared] > 0.0 {
                alternatives[local].push(shared);
            }
        }
    }

    // 至少两个不同的局部列、在不同需求行中能用同一个共享列等量替代。
    // 此判断只选分支策略，不用作界、割或可行性证明。
    let mut witness = HashMap::new();
    for (r, row) in problem.constraints.iter().enumerate() {
        let Some(sign) = le_sign(row) else { continue };
        if sign * row.rhs >= 0.0 {
            continue;
        }
        let terms = &row.expr.terms;
        for &(j, a) in terms {
            let Some(&Some(gate)) = gates.get(j) else {
                continue;
            };
            if sign * a >= 0.0 {
                continue;
            }
            for &shared in &alternatives[j] {
                let Ok(k) = terms.binary_search_by_key(&shared, |&(col, _)| col) else {
                    continue;
                };
                if terms[k].1 == a {
                    let (first_row, first_local) = *witness.entry((gate, shared)).or_insert((r, j));
                    if first_row != r && first_local != j {
                        return true;
                    }
                }
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_alternative_requires_distinct_demands_and_a_paid_gate() {
        let mut p = Problem::new("共用与局部供应", Sense::Minimize);
        for name in ["局部甲", "局部乙", "开设", "共用"] {
            p.add_var(name, VarType::Binary).unwrap();
        }
        p.set_objective(vec![(2, 10.0), (3, 6.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 2.0), (1, 1.0), (2, -3.0)], RelOp::Le, 0.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 1.0), (3, 1.0)], RelOp::Le, 1.0)
            .unwrap();
        p.add_constraint(None, vec![(1, -2.0), (3, -2.0)], RelOp::Ge, -2.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 2.0), (3, 2.0)], RelOp::Ge, 1.0)
            .unwrap();
        assert!(!prefers_pseudocost(&p));
        p.add_constraint(None, vec![(1, -3.0), (3, -2.0)], RelOp::Le, -1.0)
            .unwrap();
        assert!(!prefers_pseudocost(&p));
        p.constraints.last_mut().unwrap().expr.terms[1].1 = -3.0;
        assert!(prefers_pseudocost(&p));
        p.constraints.reverse();
        assert!(prefers_pseudocost(&p));
        p.set_objective(vec![(3, 6.0)], 0.0).unwrap();
        assert!(!prefers_pseudocost(&p));
    }
}
