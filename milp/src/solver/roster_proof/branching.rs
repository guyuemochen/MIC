//! 完备分支：整数合计的向下/向上界，或正容量行中互斥选择的两组禁用域。

use super::guarded_capacity;
use crate::model::{Problem, Tolerances, VarType};
use crate::simplex::Bounds;
use std::collections::HashSet;
use std::time::Instant;

pub(super) enum Branch {
    Integer(usize, f64),
    Partition(usize, usize),
}

pub(super) struct Branching {
    totals: Vec<(usize, f64)>,
    // (工作行号，互斥整数列)。正系数支撑按工作行号索引。
    groups: Vec<(usize, Vec<usize>)>,
    positive: Vec<Vec<u64>>,
}

impl Branching {
    pub(super) fn new(
        p: &Problem,
        totals: Vec<(usize, f64)>,
        tol: &Tolerances,
        deadline: Instant,
    ) -> Option<Self> {
        let words = p.vars.len().div_ceil(64);
        // 检查额外索引大小再分配。无法建表仍可返回完整根界。
        if words.checked_mul(p.constraints.len())? > 8_000_000 {
            return None;
        }
        let mut groups = Vec::new();
        let mut positive = Vec::with_capacity(p.constraints.len());
        for (i, r) in p.constraints.iter().enumerate() {
            if Instant::now() >= deadline {
                return None;
            }
            let mut bits = vec![0u64; words];
            for &(j, w) in &r.expr.terms {
                if w > 0.0 {
                    bits[j / 64] |= 1 << (j % 64);
                }
            }
            positive.push(bits);
            if let Some(cap) = guarded_capacity(r, p, tol) {
                let group: Vec<usize> = r
                    .expr
                    .terms
                    .iter()
                    .filter_map(|&(j, w)| {
                        (w > cap * 0.5 && p.vars[j].var_type != VarType::Continuous).then_some(j)
                    })
                    .collect();
                if group.len() >= 2 {
                    groups.push((i, group));
                }
            }
        }
        Some(Self {
            totals,
            groups,
            positive,
        })
    }

    fn contains(&self, row: usize, j: usize) -> bool {
        self.positive[row][j / 64] & (1 << (j % 64)) != 0
    }

    pub(super) fn options(
        &self,
        p: &Problem,
        x: &[f64],
        tol: &Tolerances,
        deadline: Instant,
    ) -> Vec<Branch> {
        let frac = |j: usize| (x[j] - x[j].round()).abs();
        let mut totals: Vec<_> = self
            .totals
            .iter()
            .copied()
            .filter(|&(j, _)| frac(j) > tol.integrality)
            .collect();
        totals.sort_by(|&(a, ca), &(b, cb)| {
            (frac(b) * cb).total_cmp(&(frac(a) * ca)).then(a.cmp(&b))
        });
        let mut options: Vec<_> = totals
            .into_iter()
            .take(6)
            .map(|(j, _)| Branch::Integer(j, x[j].floor()))
            .collect();
        if options.len() == 6 {
            return options;
        }
        // 两个候选对当前分数点作相同的分割时，仅保留最先出现的一对工作行。
        let mut seen = HashSet::new();
        let mut ranked = Vec::new();
        for (index, (row, group)) in self.groups.iter().enumerate() {
            if Instant::now() >= deadline {
                return options;
            }
            let positive: Vec<_> = group.iter().copied().filter(|&j| x[j] > 1e-7).collect();
            if !(2..=256).contains(&positive.len()) {
                continue;
            }
            for h in 0..self.positive.len() {
                if h == *row {
                    continue;
                }
                let (inside, outside): (Vec<_>, Vec<_>) =
                    positive.iter().copied().partition(|&j| self.contains(h, j));
                if inside.is_empty() || outside.is_empty() {
                    continue;
                }
                let mass = |js: &[usize]| js.iter().map(|&j| x[j]).sum::<f64>();
                let score = mass(&inside).min(mass(&outside));
                if score < 0.02 {
                    continue;
                }
                let key = if inside < outside {
                    (inside, outside)
                } else {
                    (outside, inside)
                };
                if seen.insert(key) {
                    ranked.push((score, index, h));
                }
                // 最多用前六个候选探测。此上限也让高度稠密的分数点有界退化。
                if ranked.len() >= 4096 {
                    break;
                }
            }
            if ranked.len() >= 4096 {
                break;
            }
        }
        ranked.sort_unstable_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
        options.extend(
            ranked
                .into_iter()
                .take(6 - options.len())
                .map(|(_, g, h)| Branch::Partition(g, h)),
        );
        if options.is_empty() {
            let mut best = None;
            let mut distance = tol.integrality;
            for (j, v) in p.vars.iter().enumerate() {
                if v.var_type != VarType::Continuous && frac(j) > distance {
                    best = Some(j);
                    distance = frac(j);
                }
            }
            if let Some(j) = best {
                options.push(Branch::Integer(j, x[j].floor()));
            }
        }
        options
    }

    pub(super) fn restrict(&self, branch: &Branch, side: usize, bounds: &mut Bounds) {
        match *branch {
            Branch::Integer(j, floor) => {
                if side == 0 {
                    bounds.ub[j] = bounds.ub[j].min(floor);
                } else {
                    bounds.lb[j] = bounds.lb[j].max(floor + 1.0);
                }
            }
            Branch::Partition(g, h) => {
                // 整数解至多有一列取正，故禁用交集或禁用其补集覆盖所有整数解，两侧可重叠。
                for &j in &self.groups[g].1 {
                    if self.contains(h, j) == (side == 0) {
                        bounds.ub[j] = bounds.ub[j].min(0.0);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RelOp, Sense};
    use std::time::Duration;

    #[test]
    fn both_partition_domains_cover_every_feasible_point_with_general_integers() {
        let mut p = Problem::new("成对分区", Sense::Minimize);
        for name in ["a", "b", "c", "d"] {
            p.add_var(name, VarType::Integer).unwrap();
        }
        p.add_constraint(
            None,
            vec![(0, 3.), (1, 2.), (2, 2.), (3, 1.)],
            RelOp::Le,
            3.,
        )
        .unwrap();
        p.add_constraint(None, vec![(0, 1.), (3, 1.)], RelOp::Ge, 0.)
            .unwrap();
        let branching = Branching::new(
            &p,
            vec![],
            &Tolerances::default(),
            Instant::now() + Duration::from_secs(1),
        )
        .unwrap();
        let mut left = Bounds::from_problem(&p);
        let mut right = left.clone();
        branching.restrict(&Branch::Partition(0, 1), 0, &mut left);
        branching.restrict(&Branch::Partition(0, 1), 1, &mut right);
        assert_eq!(left.ub, [0., f64::INFINITY, f64::INFINITY, f64::INFINITY]);
        assert_eq!(right.ub, [f64::INFINITY, 0., 0., f64::INFINITY]);
        for a in 0..=3 {
            for b in 0..=3 {
                for c in 0..=3 {
                    for d in 0..=3 {
                        let x = [a as f64, b as f64, c as f64, d as f64];
                        if p.constraints[0].expr.eval(&x) <= 3. {
                            assert!([&left, &right]
                                .iter()
                                .any(|bounds| x.iter().zip(&bounds.ub).all(|(x, u)| x <= u)));
                        }
                    }
                }
            }
        }
    }
}
