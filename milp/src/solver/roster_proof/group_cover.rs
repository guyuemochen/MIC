//! 用一条等式与共享余量的下界行消元，给整数余量构造覆盖不等式。
use super::{guarded_capacity, usable_tolerances};
use crate::cuts::Cut;
use crate::model::{Constraint, Problem, RelOp, Tolerances, VarType};
use std::{
    collections::{BTreeMap, BTreeSet},
    time::Instant,
};

fn row_error(row: &Constraint, tol: &Tolerances) -> f64 {
    let norm = row
        .expr
        .terms
        .iter()
        .fold(0.0_f64, |s, t| (s + t.1.abs()).next_up());
    let scale = ((1. + row.rhs.abs()).next_up() + norm).next_up();
    (scale
        * (crate::solver::FINAL_VIOLATION_REL
            + tol.integrality
            + 10. * tol.feasibility
            + 4. * tol.zero))
        .next_up()
}

pub(in crate::solver) fn separate(p: &Problem, tol: &Tolerances, end: Instant) -> Vec<Cut> {
    covers(p, tol, Some(end))
}

/// 根松弛可直接采用全域有效的整数零上界，即使旧根点恰好也满足这个界。
/// 在常规预处理完成后就地收紧，保留列号与等式，后续仍可据此还原候选。
pub(in crate::solver) fn fix_zero_bounds(p: &mut Problem, tol: &Tolerances, end: Option<Instant>) {
    if p.vars.len() < 4096
        || p.vars.len() / 8 < p.constraints.len()
        || !p
            .constraints
            .iter()
            .any(|r| r.op == RelOp::Le && crate::heuristics::is_choice_row(p, r))
    {
        return;
    }
    for cut in covers(p, tol, end) {
        if let [(j, 1.)] = cut.terms.as_slice() {
            if cut.rhs == 0. && p.vars[*j].lb == 0. && p.vars[*j].ub > 0. {
                p.vars[*j].ub = 0.;
            }
        }
    }
}

fn covers(p: &Problem, tol: &Tolerances, end: Option<Instant>) -> Vec<Cut> {
    let expired = || end.is_some_and(|end| Instant::now() >= end);
    if !usable_tolerances(tol) || p.vars.len() < 4096 || p.vars.iter().any(|v| v.lb < 0.) {
        return Vec::new();
    }
    let mut groups = vec![None; p.vars.len()];
    let mut shared = vec![(None, None, 0u8); p.vars.len()];
    for (i, r) in p.constraints.iter().enumerate() {
        if expired() {
            return Vec::new();
        }
        // 一列只分入一个已证互斥组。重叠的互斥行可以忽略，仍保有每组至多取一的证明。
        if r.op == RelOp::Le
            && crate::heuristics::is_choice_row(p, r)
            && guarded_capacity(r, p, tol).is_some_and(|cap| 2. * r.rhs > cap)
        {
            for &(j, _) in &r.expr.terms {
                groups[j].get_or_insert(i);
            }
        }
        for &(j, a) in &r.expr.terms {
            if p.vars[j].var_type != VarType::Integer || !p.vars[j].ub.is_infinite() {
                continue;
            }
            let item = &mut shared[j];
            item.2 = item.2.saturating_add(1);
            match r.op {
                RelOp::Eq => item.0 = Some((i, a)),
                RelOp::Ge => item.1 = Some((i, a)),
                RelOp::Le => {}
            }
        }
    }
    let mut pairs = BTreeSet::new();
    for (eq, side, count) in shared {
        if let (Some((i, a)), Some((k, b)), 2) = (eq, side, count) {
            // 只缩放二的整数幂；反算核实有限乘积没有下溢舍入，避免伪造消元。
            let scale = b / a;
            let scaled_rhs = p.constraints[i].rhs * scale;
            if scale.is_normal()
                && scale > 0.
                && scale.to_bits() & ((1u64 << 52) - 1) == 0
                && b.is_finite()
                && a * scale == b
                && b / scale == a
                && scaled_rhs.is_finite()
                && scaled_rhs == p.constraints[k].rhs
                && scaled_rhs / scale == p.constraints[i].rhs
            {
                pairs.insert((i, k, scale.to_bits()));
            }
        }
    }
    let mut cuts = Vec::new();
    for (i, k, scale) in pairs {
        let scale = f64::from_bits(scale);
        if cuts.len() >= 512 || expired() {
            break;
        }
        let eq = &p.constraints[i];
        let side = &p.constraints[k];
        // side - scale * eq >= 0。对每个正余项向上舍入；不能证明为负的微小项同样计入上界。
        let mut delta: BTreeMap<_, _> = side.expr.terms.iter().copied().collect();
        if eq.expr.terms.iter().any(|&(_, a)| {
            let scaled = scale * a;
            !scaled.is_finite() || scaled / scale != a
        }) {
            continue;
        }
        for &(j, a) in &eq.expr.terms {
            let a = scale * a;
            let b = delta.entry(j).or_default();
            *b = if *b == a { 0. } else { (*b - a).next_up() };
        }
        let Some((&z, &neg)) = delta
            .iter()
            .filter(|&(&j, &a)| a < -tol.zero && p.vars[j].var_type == VarType::Integer)
            .min_by(|a, b| a.1.total_cmp(b.1))
        else {
            continue;
        };
        let step = -neg;
        let mut rest = (((row_error(eq, tol) * scale).next_up() + row_error(side, tol)).next_up()
            / step)
            .next_up();
        let mut grouped: BTreeMap<usize, Vec<(usize, f64)>> = BTreeMap::new();
        for (&j, &a) in &delta {
            if a <= 0. {
                continue;
            }
            let weight = (a / step).next_up();
            if weight <= tol.zero || p.vars[j].var_type == VarType::Continuous {
                rest = weight.mul_add(p.vars[j].ub, rest).next_up();
            } else if let Some(g) = groups[j] {
                grouped.entry(g).or_default().push((j, weight));
            } else if p.vars[j].ub <= 1. {
                grouped
                    .entry(p.constraints.len() + j)
                    .or_default()
                    .push((j, weight));
            } else {
                rest = weight.mul_add(p.vars[j].ub, rest).next_up();
            }
        }
        if !rest.is_finite() || grouped.len() > 6 {
            continue;
        }
        let groups: Vec<_> = grouped.into_values().collect();
        let maxima: Vec<_> = groups
            .iter()
            .map(|g| g.iter().map(|t| t.1).fold(0., f64::max))
            .collect();
        // 移出若干组后，余下组的最大贡献可向下取整：左侧 z 与移出组的取整权重之和均为整数。
        let rhs: Vec<_> = (0..1usize << groups.len())
            .map(|mask| {
                maxima
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| mask & (1 << i) == 0)
                    .fold(rest, |s, (_, &w)| (s + w).next_up())
                    .floor()
            })
            .collect();
        for (mask, &b) in rhs.iter().enumerate() {
            if !b.is_finite() || b >= (1u64 << 50) as f64 {
                continue;
            }
            if mask == 0 && p.vars[z].ub <= b {
                continue;
            }
            if (0..groups.len()).any(|i| mask & (1 << i) != 0 && rhs[mask ^ (1 << i)] <= b) {
                continue;
            }
            let mut terms = vec![(z, 1.)];
            for (i, g) in groups
                .iter()
                .enumerate()
                .filter(|(i, _)| mask & (1 << i) != 0)
            {
                if !maxima[i].is_finite() || maxima[i] > (1u64 << 40) as f64 {
                    terms.clear();
                    break;
                }
                terms.extend(g.iter().map(|&(j, w)| (j, -w.ceil())));
            }
            if !terms.is_empty() {
                cuts.push(Cut {
                    terms,
                    op: RelOp::Le,
                    rhs: b,
                });
            }
        }
    }
    cuts
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Sense, Variable};
    use std::time::Duration;

    #[test]
    fn single_choice_bounds_integer_residual_before_root() {
        let mut p = Problem::new("共享余量与互斥组选项", Sense::Minimize);
        let a = p.add_var("甲", VarType::Binary).unwrap();
        let b = p.add_var("乙", VarType::Binary).unwrap();
        let s = p.add_var("共同项", VarType::Integer).unwrap();
        let z = p.add_var("余量", VarType::Integer).unwrap();
        p.add_vars(
            (4..4096)
                .map(|i| Variable {
                    name: format!("填充{i}"),
                    var_type: VarType::Binary,
                    lb: 0.,
                    ub: 1.,
                })
                .collect(),
        )
        .unwrap();
        p.add_constraint(None, vec![(a, 1.), (b, 1.)], RelOp::Le, 1.)
            .unwrap();
        p.add_constraint(None, vec![(s, 1.), (z, 1.)], RelOp::Eq, 1.)
            .unwrap();
        p.add_constraint(None, vec![(a, 0.5), (b, 0.5), (s, 1.)], RelOp::Ge, 1.)
            .unwrap();
        let tol = Tolerances::default();
        let mut fixed = p.clone();
        fix_zero_bounds(&mut fixed, &tol, None);
        assert_eq!(fixed.vars[z].ub, 0.);
        assert_eq!(fixed.vars.len(), p.vars.len());
        fix_zero_bounds(&mut p, &tol, Some(Instant::now()));
        assert!(p.vars[z].ub.is_infinite());
        // 少了互斥行时，两项均可取一，整数余量 1 是可行的。
        let mut independent = p.clone();
        independent.constraints.remove(0);
        fix_zero_bounds(&mut independent, &tol, None);
        assert!(independent.vars[z].ub.is_infinite());
        let mut continuous = p.clone();
        continuous.vars[z].var_type = VarType::Continuous;
        fix_zero_bounds(&mut continuous, &tol, None);
        assert!(continuous.vars[z].ub.is_infinite());
        // 接近整数门槛的系数须留足误差，不能把一个本可取一的余量固定。
        p.constraints[2].expr.terms[0].1 = 1. - tol.feasibility;
        fix_zero_bounds(&mut p, &tol, None);
        assert!(p.vars[z].ub.is_infinite());
    }

    #[test]
    fn paired_cover_requires_two_distinct_choices_and_respects_tolerance() {
        let mut p = Problem::new("两组各贡献半份", Sense::Minimize);
        for name in ["甲早", "甲晚", "乙早", "乙晚"] {
            p.add_var(name, VarType::Binary).unwrap();
        }
        let slack = p.add_var("共享余量", VarType::Integer).unwrap();
        let z = p.add_var("整数需求", VarType::Integer).unwrap();
        p.add_vars(
            (6..4096)
                .map(|i| Variable {
                    name: format!("填充{i}"),
                    var_type: VarType::Binary,
                    lb: 0.,
                    ub: 1.,
                })
                .collect(),
        )
        .unwrap();
        for g in [[0, 1], [2, 3]] {
            p.add_constraint(None, g.map(|j| (j, 1.)).to_vec(), RelOp::Le, 1.)
                .unwrap();
        }
        p.add_constraint(None, vec![(slack, 1.), (z, 1.)], RelOp::Eq, 1.)
            .unwrap();
        p.add_constraint(
            None,
            vec![(0, 0.5), (1, 0.5), (2, 0.5), (3, 0.5), (slack, 1.)],
            RelOp::Ge,
            1.,
        )
        .unwrap();
        let deadline = Instant::now() + Duration::from_secs(2);
        let cuts = separate(&p, &Tolerances::default(), deadline);
        assert!(cuts
            .iter()
            .any(|c| c.rhs == 0. && c.terms.contains(&(0, -1.)) && !c.terms.contains(&(2, -1.))));
        for mask in 0..16usize {
            let bits: Vec<_> = (0..4).map(|j| ((mask >> j) & 1) as f64).collect();
            for zz in [0., 1.] {
                if bits[0] + bits[1] > 1.
                    || bits[2] + bits[3] > 1.
                    || 0.5 * bits.iter().sum::<f64>() < zz
                {
                    continue;
                }
                assert!(cuts.iter().all(|c| c
                    .terms
                    .iter()
                    .map(|&(j, a)| a * if j == z { zz } else { bits[j] })
                    .sum::<f64>()
                    <= c.rhs));
            }
        }
        // 把下界行扩大一倍，原有整数可行域及可证覆盖割都应相同。
        let ge = p.constraints.last_mut().unwrap();
        ge.rhs *= 2.;
        for (_, a) in &mut ge.expr.terms {
            *a *= 2.;
        }
        let scaled = separate(&p, &Tolerances::default(), deadline);
        assert_eq!(cuts.len(), scaled.len());
        for (old, new) in cuts.iter().zip(&scaled) {
            assert_eq!(old.rhs, new.rhs);
            assert_eq!(old.terms, new.terms);
        }
        // 若只有右端不同，余量等式无法消去常数。
        p.constraints.last_mut().unwrap().rhs = 1.;
        assert!(separate(&p, &Tolerances::default(), deadline).is_empty());
        p.constraints.last_mut().unwrap().rhs = 2.;
        // 放宽归一化误差会让“至多取一”失去证据，不能继续发布依赖该组的割。
        let mut tol = Tolerances::default();
        tol.feasibility = 1.;
        assert!(!separate(&p, &tol, deadline).iter().any(|c| c.rhs <= 1.));
        p.vars[0].lb = -1.;
        assert!(separate(&p, &Tolerances::default(), deadline).is_empty());
    }
}
