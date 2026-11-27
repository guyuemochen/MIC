//! 在当前分数支撑上枚举互斥团，再把与团内每列均冲突的整数列加入团。

use super::guarded_capacity;
use crate::cuts::Cut;
use crate::model::{Problem, RelOp, Tolerances, VarType};
use crate::simplex::{LpSolution, LpStatus};
use std::collections::BTreeMap;
use std::time::Instant;

struct Conflicts {
    // 一个容量行的相同权重只用一张冲突位表。列只保存它参加的表号。
    by_column: Vec<Vec<usize>>,
    masks: Vec<Vec<u64>>,
    eligible: Vec<bool>,
}

impl Conflicts {
    fn build(p: &Problem, tol: &Tolerances, deadline: Instant) -> Option<Self> {
        let n = p.vars.len();
        let words = n.div_ceil(64);
        let mut eligible = vec![false; n];
        let mut rows = Vec::new();
        for r in &p.constraints {
            if Instant::now() >= deadline {
                return None;
            }
            if let Some(cap) = guarded_capacity(r, p, tol) {
                for &(j, w) in &r.expr.terms {
                    if p.vars[j].var_type != VarType::Continuous
                        && (p.vars[j].ub <= 1.0 || w > cap * 0.5)
                    {
                        eligible[j] = true;
                    }
                }
                rows.push((r, cap));
            }
        }
        let mut masks = Vec::new();
        let mut by_column = vec![Vec::new(); n];
        let mut memberships = 0usize;
        for (r, cap) in rows {
            if Instant::now() >= deadline {
                return None;
            }
            if memberships.checked_add(r.expr.terms.len())? > 2_000_000 {
                break;
            }
            let mut weights = BTreeMap::<u64, Vec<usize>>::new();
            for &(j, w) in &r.expr.terms {
                if eligible[j] {
                    weights.entry(w.to_bits()).or_default().push(j);
                    if weights.len() > 64 {
                        break;
                    }
                }
            }
            // 独特权重很多的行可略过。所有已保留的边都有完整的原行证据。
            if weights.len() > 64 || (masks.len() + weights.len()).checked_mul(words)? > 4_000_000 {
                continue;
            }
            memberships += r.expr.terms.len();
            for (&w, members) in &weights {
                let mut mask = vec![0u64; words];
                for (&other, cols) in &weights {
                    if f64::from_bits(w) + f64::from_bits(other) > cap {
                        for &j in cols {
                            mask[j / 64] |= 1 << (j % 64);
                        }
                    }
                }
                for &j in members {
                    by_column[j].push(masks.len());
                }
                masks.push(mask);
            }
        }
        Some(Self {
            by_column,
            masks,
            eligible,
        })
    }

    fn adjacent(&self, j: usize, k: usize) -> bool {
        let bit = 1 << (k % 64);
        self.by_column[j]
            .iter()
            .any(|&m| self.masks[m][k / 64] & bit != 0)
    }

    /// 与一列的全部冲突邻居求交。只有剩余非空的机器字才需要查询该列的容量行。
    fn intersect(&self, remaining: &mut [u64], j: usize) {
        for (i, word) in remaining
            .iter_mut()
            .enumerate()
            .filter(|(_, word)| **word != 0)
        {
            let neighbors = self.by_column[j]
                .iter()
                .fold(0, |u, &m| u | self.masks[m][i]);
            *word &= neighbors;
        }
        remaining[j / 64] &= !(1 << (j % 64));
    }
}

struct Seeds<'a> {
    frac: &'a [usize],
    x: &'a [f64],
    neighbors: Vec<u128>,
    found: Vec<(f64, Vec<usize>)>,
    calls: usize,
    deadline: Instant,
}

impl Seeds<'_> {
    /// Bron–Kerbosch 带主元。并列主元与候选按较小工作列号优先。
    fn visit(&mut self, clique: &mut Vec<usize>, mut possible: u128, mut excluded: u128) {
        if self.calls >= 32_768 || self.found.len() >= 4096 || Instant::now() >= self.deadline {
            return;
        }
        self.calls += 1;
        let mut pivots = possible | excluded;
        if pivots == 0 {
            let excess = clique.iter().map(|&j| self.x[j]).sum::<f64>() - 1.0;
            if excess > 1e-6 {
                self.found.push((excess, clique.clone()));
            }
            return;
        }
        let mut pivot = pivots.trailing_zeros() as usize;
        let mut overlap = (possible & self.neighbors[pivot]).count_ones();
        while pivots != 0 {
            let k = pivots.trailing_zeros() as usize;
            let count = (possible & self.neighbors[k]).count_ones();
            if count > overlap {
                pivot = k;
                overlap = count;
            }
            pivots &= pivots - 1;
        }
        let mut candidates = possible & !self.neighbors[pivot];
        while candidates != 0 {
            let k = candidates.trailing_zeros() as usize;
            let bit = 1u128 << k;
            clique.push(self.frac[k]);
            self.visit(
                clique,
                possible & self.neighbors[k],
                excluded & self.neighbors[k],
            );
            clique.pop();
            possible &= !bit;
            excluded |= bit;
            candidates &= !bit;
        }
    }
}

pub(in crate::solver) fn separate(
    p: &Problem,
    root: &LpSolution,
    rc: &[f64],
    tol: &Tolerances,
    deadline: Instant,
) -> Vec<Cut> {
    let n = p.vars.len();
    if !super::usable_tolerances(tol)
        || root.status != LpStatus::Optimal
        || root.x.len() != n
        || rc.len() != n
        || Instant::now() >= deadline
    {
        return Vec::new();
    }
    let Some(graph) = Conflicts::build(p, tol, deadline) else {
        return Vec::new();
    };
    let frac: Vec<_> = (0..n)
        .filter(|&j| graph.eligible[j] && root.x[j] > 1e-7 && root.x[j] < 1.0 - 1e-7)
        .take(128)
        .collect();
    if frac.len() < 2 {
        return Vec::new();
    }
    let mut neighbors = vec![0u128; frac.len()];
    for (i, &j) in frac.iter().enumerate() {
        for (k, &other) in frac[..i].iter().enumerate() {
            if graph.adjacent(j, other) {
                neighbors[i] |= 1 << k;
                neighbors[k] |= 1 << i;
            }
        }
    }
    let mut seeds = Seeds {
        frac: &frac,
        x: &root.x,
        neighbors,
        found: Vec::new(),
        calls: 0,
        deadline,
    };
    let all = u128::MAX >> (128 - frac.len());
    seeds.visit(&mut Vec::new(), all, 0);
    let mut found = seeds.found;
    let violation_key = |value: f64| (value * 1e9).round_ties_even();
    found.sort_unstable_by(|a, b| {
        violation_key(b.0)
            .total_cmp(&violation_key(a.0))
            .then(a.1.cmp(&b.1))
    });
    let mut order: Vec<_> = (0..n).filter(|&j| graph.eligible[j]).collect();
    let reduced_key = |j: usize| (rc[j].max(0.0) * 1e5).round_ties_even();
    order.sort_unstable_by(|&a, &b| reduced_key(a).total_cmp(&reduced_key(b)).then(a.cmp(&b)));
    let mut cuts: Vec<Cut> = Vec::new();
    for (_, mut clique) in found.into_iter().take(24) {
        if Instant::now() >= deadline {
            break;
        }
        let mut remaining = vec![u64::MAX; n.div_ceil(64)];
        for &j in &clique {
            graph.intersect(&mut remaining, j);
        }
        for (i, &j) in order.iter().enumerate() {
            if i % 128 == 0 && Instant::now() >= deadline {
                return cuts;
            }
            if remaining[j / 64] & (1 << (j % 64)) != 0 {
                clique.push(j);
                graph.intersect(&mut remaining, j);
            }
        }
        clique.sort_unstable();
        if cuts.iter().any(|cut| {
            let overlap = cut
                .terms
                .iter()
                .filter(|(j, _)| clique.binary_search(j).is_ok())
                .count();
            overlap as f64 > 0.98 * clique.len() as f64
        }) {
            continue;
        }
        cuts.push(Cut {
            terms: clique.into_iter().map(|j| (j, 1.0)).collect(),
            op: RelOp::Le,
            rhs: 1.0,
        });
        if cuts.len() == 8 {
            break;
        }
    }
    cuts
}
