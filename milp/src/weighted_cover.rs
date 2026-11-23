//! 正整数费用的集合覆盖。每行至少选一列，选中一列即可同时覆盖它出现的所有行。
//! 贪心与局部修补先产生可行解，整数下界搜索负责证明能否找到更便宜的解。

use crate::branch_bound::BbStats;
use crate::error::{MilpError, Result};
use crate::model::{Problem, RelOp, Sense, SolveStatus, Tolerances, VarType};
use std::time::Instant;

mod primal;
mod search;

struct Cover {
    rows: Vec<Vec<usize>>,
    cols: Vec<Vec<usize>>,
    costs: Vec<i64>,
}

fn expired(deadline: Option<Instant>) -> bool {
    deadline.is_some_and(|end| Instant::now() >= end)
}

impl Cover {
    fn detect(p: &Problem, tol: Tolerances) -> Option<Self> {
        let n = p.vars.len();
        if !(64..=1024).contains(&n)
            || !(32..=1024).contains(&p.constraints.len())
            || p.sense != Sense::Minimize
            || p.obj.terms.len() != n
            || !p
                .vars
                .iter()
                .all(|v| v.var_type == VarType::Binary && v.lb == 0. && v.ub == 1.)
            || !p
                .obj
                .terms
                .iter()
                .all(|&(_, c)| (1.0..=1_000_000.0).contains(&c) && c == c.trunc())
            || !tol.feasibility.is_finite()
            || tol.feasibility < 0.
            || !tol.integrality.is_finite()
            || tol.integrality < 0.
            || p.obj_offset.abs() + n as f64 * 1_000_000. >= (1u64 << 48) as f64
        {
            return None;
        }
        let costs: Vec<_> = p.obj.terms.iter().map(|&(_, c)| c as i64).collect();
        if costs.iter().all(|&c| c == costs[0]) {
            return None;
        }
        let mut rows = Vec::with_capacity(p.constraints.len());
        let mut cols = vec![Vec::new(); n];
        let mut nz = 0;
        for (i, r) in p.constraints.iter().enumerate() {
            let len = r.expr.terms.len();
            nz += len;
            if nz > 100_000
                || len == 0
                || 10. * tol.feasibility + len as f64 * tol.integrality >= 1.
            {
                return None;
            }
            let unit = match (r.op, r.rhs) {
                (RelOp::Ge, 1.) => 1.,
                (RelOp::Le, -1.) => -1.,
                _ => return None,
            };
            if r.expr.terms.iter().any(|&(_, a)| a != unit) {
                return None;
            }
            let row: Vec<_> = r.expr.terms.iter().map(|&(j, _)| j).collect();
            for &j in &row {
                cols[j].push(i);
            }
            rows.push(row);
        }
        Some(Self { rows, cols, costs })
    }

    fn value(&self, point: &[bool]) -> i64 {
        point
            .iter()
            .zip(&self.costs)
            .filter_map(|(&on, &c)| on.then_some(c))
            .sum()
    }
}

pub(crate) fn solve(
    p: &Problem,
    tol: Tolerances,
    deadline: Option<Instant>,
    limit: u64,
    accept: &mut dyn FnMut(&[f64]),
) -> Result<Option<(SolveStatus, Option<f64>, BbStats)>> {
    let Some(m) = Cover::detect(p, tol) else {
        return Ok(None);
    };
    let mut stats = BbStats::default();
    if expired(deadline) || limit == 0 {
        return Ok(Some((SolveStatus::LimitReached, Some(p.obj_offset), stats)));
    }
    let mut error = None;
    let mut publish = |bits: &[bool]| {
        let mut x: Vec<_> = bits.iter().map(|&b| if b { 1. } else { 0. }).collect();
        match crate::solver::final_check(p, &mut x) {
            Ok(()) => accept(&x),
            Err(e) => error = Some(e),
        }
    };
    // 时限很短时也要留时间做证明。贪心首解在第一次局部重建之前发布。
    let seed_deadline = deadline.map(|end| {
        let now = Instant::now();
        now + end.saturating_duration_since(now) / 8
    });
    let mut best = m.seed(seed_deadline, &mut publish);
    let (proved, nodes, lower) = search::search(
        &m.rows,
        &m.cols,
        &m.costs,
        &mut best,
        deadline,
        limit,
        &mut publish,
    );
    if let Some(e) = error {
        return Err(e);
    }
    stats.nodes = nodes;
    let value = m.value(&best);
    if lower > value || (proved && lower != value) {
        return Err(MilpError::Numerical("集合覆盖的费用与证明界不符".into()));
    }
    let status = if proved && !expired(deadline) {
        SolveStatus::Optimal
    } else {
        SolveStatus::LimitReached
    };
    Ok(Some((status, Some(p.obj_offset + lower as f64), stats)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signed_cover_with_offset_and_limits() {
        let mut p = Problem::new("cover", Sense::Minimize);
        for j in 0..256 {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        let mut c = vec![(0, 4.), (1, 2.), (2, 3.)];
        c.extend((3..256).map(|j| (j, 20.)));
        p.set_objective(c, -2.5).unwrap();
        for _ in 0..16 {
            p.add_constraint(None, vec![(0, 1.), (1, 1.)], RelOp::Ge, 1.)
                .unwrap();
            p.add_constraint(None, vec![(0, -1.), (2, -1.)], RelOp::Le, -1.)
                .unwrap();
        }
        let mut last = None;
        let (status, bound, _) = solve(&p, Tolerances::default(), None, 1000, &mut |x| {
            assert!(x[0] + x[1] >= 1. && x[0] + x[2] >= 1.);
            last = Some(p.obj_offset + p.obj.eval(x));
        })
        .unwrap()
        .unwrap();
        assert_eq!(
            (status, bound, last),
            (SolveStatus::Optimal, Some(1.5), Some(1.5))
        );
        let (status, bound, stats) = solve(
            &p,
            Tolerances::default(),
            Some(Instant::now()),
            8,
            &mut |_| panic!("已过期限时不应开始构造可行解"),
        )
        .unwrap()
        .unwrap();
        assert_eq!(
            (status, bound, stats.nodes),
            (SolveStatus::LimitReached, Some(-2.5), 0)
        );
        p.obj.terms[1].1 = 2.5;
        assert!(Cover::detect(&p, Tolerances::default()).is_none());
        p.obj.terms[1].1 = 2.;
        p.constraints[1].expr.terms[1].1 = -2.;
        assert!(Cover::detect(&p, Tolerances::default()).is_none());
    }
}
