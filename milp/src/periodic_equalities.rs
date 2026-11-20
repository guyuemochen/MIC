//! 周期等式的初解。共享列给定后，一行的私有连续列可在其区间内补差，
//! 私有整数列可加减若干个周期。只调整共享列来减小各行到可行区间的距离，
//! 最后回填私有列并按原模型逐行核验。搜索没有找到解时交回通用求解器。

use crate::model::{Problem, RelOp, Sense, VarType};
use std::time::Instant;

fn expired(deadline: Option<Instant>) -> bool {
    deadline.is_some_and(|d| Instant::now() >= d)
}

fn small_integer(value: f64) -> Option<i32> {
    (value.abs() <= 32767.0 && value.fract() == 0.0).then_some(value as i32)
}

#[derive(Clone, Copy)]
struct Cycles {
    column: usize,
    coefficient: i32,
    low: i32,
    high: i32,
}

struct Row {
    terms: Vec<(usize, i32)>,
    local: Option<(usize, i32)>,
    cycles: Option<Cycles>,
    rhs: i32,
    lower: i32,
    upper: i32,
    minimum: i32,
    distances: Vec<i32>,
    sum: i32,
    weight: i64,
}

impl Row {
    // 有限范围内只需比较达到区间左端所需周数的向下和向上取整。
    fn nearest(&self, sum: i32) -> (i32, i32) {
        let Some(c) = self.cycles else {
            return ((self.lower - sum).max(sum - self.upper).max(0), 0);
        };
        let sign = c.coefficient.signum();
        let period = c.coefficient.abs();
        let floor = (self.lower - sum).div_euclid(period);
        let candidate = |k: i32| {
            let q = (k * sign).clamp(c.low, c.high);
            let activity = sum + q * c.coefficient;
            ((self.lower - activity).max(activity - self.upper).max(0), q)
        };
        candidate(floor).min(candidate(floor + 1))
    }

    fn distance(&self, sum: i32) -> i32 {
        self.distances[(sum - self.minimum) as usize]
    }
}

struct Repair {
    rows: Vec<Row>,
    constraints: Vec<usize>,
    incidence: Vec<Vec<(usize, i32)>>,
    bounds: Vec<(i32, i32)>,
    values: Vec<i32>,
}

impl Repair {
    fn detect(p: &Problem, deadline: Option<Instant>) -> Option<Self> {
        let n = p.vars.len();
        if n > 4096
            || p.constraints.len() > 4096
            || !p.constraints.iter().all(|r| r.op == RelOp::Eq)
            || !p.vars.iter().any(|v| v.var_type == VarType::Continuous)
            || !p.vars.iter().any(|v| v.var_type != VarType::Continuous)
        {
            return None;
        }
        let mut bounds = Vec::with_capacity(n);
        for var in &p.vars {
            let lo = small_integer(var.lb)?;
            let hi = small_integer(var.ub)?;
            if !(0..=128).contains(&(hi - lo)) {
                return None;
            }
            bounds.push((lo, hi));
        }
        let mut degree = vec![0; n];
        for row in &p.constraints {
            if expired(deadline) {
                return None;
            }
            for &(j, _) in &row.expr.terms {
                degree[j] += 1;
            }
        }
        let mut capacity = 1_000_000usize;
        let mut period = None;
        let mut rows = Vec::new();
        let mut incidence = vec![Vec::new(); n];
        let mut constraints = Vec::new();
        for input in &p.constraints {
            if expired(deadline) {
                return None;
            }
            let mut rhs = small_integer(input.rhs)?;
            let mut local = None;
            let mut cycles = None;
            let mut terms = Vec::new();
            for &(j, a) in &input.expr.terms {
                let c = small_integer(a)?;
                let (lo, hi) = bounds[j];
                if lo == hi {
                    rhs = rhs.checked_sub(c.checked_mul(lo)?)?;
                } else if p.vars[j].var_type != VarType::Continuous {
                    if degree[j] != 1
                        || c.abs() < 2
                        || c.abs() > 1024
                        || cycles.is_some()
                        || period.is_some_and(|k| k != c.abs())
                    {
                        return None;
                    }
                    period = Some(c.abs());
                    cycles = Some(Cycles {
                        column: j,
                        coefficient: c,
                        low: lo,
                        high: hi,
                    });
                } else {
                    if c.abs() != 1 {
                        return None;
                    }
                    // 单行列优先留跨度最大的补差；原列升序使同宽时取最小下标。
                    if degree[j] == 1
                        && local.is_none_or(|(old, _): (usize, i32)| {
                            hi - lo > bounds[old].1 - bounds[old].0
                        })
                    {
                        local = Some((j, c));
                    }
                    terms.push((j, c));
                }
            }
            // 将活动度与右端限制在 i32 的安全计算范围内；超出只回退。
            if rhs.unsigned_abs() > 100_000_000 {
                return None;
            }
            let (lower, upper) = if let Some((j, c)) = local {
                terms.retain(|&(k, _)| k != j);
                let (lo, hi) = bounds[j];
                (rhs - (c * lo).max(c * hi), rhs - (c * lo).min(c * hi))
            } else {
                (rhs, rhs)
            };
            let mut minimum = 0;
            let mut maximum = 0;
            let mut sum = 0;
            for &(j, c) in &terms {
                let (lo, hi) = bounds[j];
                minimum += (lo * c).min(hi * c);
                maximum += (lo * c).max(hi * c);
                sum += lo * c;
            }
            let mut row = Row {
                terms,
                local,
                cycles,
                rhs,
                lower,
                upper,
                minimum,
                distances: Vec::new(),
                sum,
                weight: 1,
            };
            if !(cycles.is_none() && minimum >= lower && maximum <= upper) {
                let span = (maximum - minimum + 1) as usize;
                capacity = capacity.checked_sub(span)?;
                row.distances.reserve(span);
                for v in minimum..=maximum {
                    if v % 4096 == 0 && expired(deadline) {
                        return None;
                    }
                    row.distances.push(row.nearest(v).0);
                }
                for &(j, c) in &row.terms {
                    incidence[j].push((rows.len(), c));
                }
                constraints.push(rows.len());
            }
            rows.push(row);
        }
        period?;
        let values = bounds.iter().map(|&(lo, _)| lo).collect();
        Some(Self {
            rows,
            constraints,
            incidence,
            bounds,
            values,
        })
    }

    fn find(&mut self, p: &Problem, deadline: Option<Instant>) -> Option<()> {
        // 修复同样多的违反量时，优先沿原目标系数较有利的方向移动。
        let mut direction = vec![0i32; self.values.len()];
        let sign = if p.sense == Sense::Minimize {
            1.0
        } else {
            -1.0
        };
        for &(j, c) in &p.obj.terms {
            direction[j] = (sign * c).signum() as i32;
        }
        let mut tabu_until = vec![0; self.values.len()];
        let mut violated = Vec::new();
        let mut record = i64::MAX;
        // 全部评价次数合并计限；困难模型只付有限的启发式代价。
        let mut work = 3_000_000usize;
        for step in 0..1500 {
            if expired(deadline) {
                return None;
            }
            violated.clear();
            let mut total = 0i64;
            for &i in &self.constraints {
                let r = &self.rows[i];
                let d = r.distance(r.sum);
                if d != 0 {
                    violated.push(i);
                    total += i64::from(d);
                }
            }
            if violated.is_empty() {
                return Some(());
            }
            record = record.min(total);
            let target = &self.rows[violated[step % violated.len()]];
            let mut best = None;
            for &(j, a) in &target.terms {
                if expired(deadline) {
                    return None;
                }
                let (lo, hi) = self.bounds[j];
                for v in lo..=hi {
                    let delta = v - self.values[j];
                    if target.distance(target.sum + delta * a) >= target.distance(target.sum) {
                        continue;
                    }
                    work = work.checked_sub(self.incidence[j].len())?;
                    let mut score = 0;
                    let mut change = 0;
                    for &(i, b) in &self.incidence[j] {
                        let r = &self.rows[i];
                        let d = i64::from(r.distance(r.sum + delta * b) - r.distance(r.sum));
                        score += r.weight * d;
                        change += d;
                    }
                    if step < tabu_until[j] && total + change >= record {
                        continue;
                    }
                    let proposed = (score, direction[j] * delta, j, v);
                    if best.is_none_or(|old| proposed < old) {
                        best = Some(proposed);
                    }
                }
            }
            if best.is_none_or(|(score, price, _, _)| (score, price) >= (0, 0)) && step % 3 == 0 {
                for &i in &violated {
                    self.rows[i].weight += 1;
                }
            }
            if let Some((_, _, j, v)) = best {
                let delta = v - self.values[j];
                for &(i, c) in &self.incidence[j] {
                    self.rows[i].sum += delta * c;
                }
                self.values[j] = v;
                tabu_until[j] = step + 1 + step % 9;
            }
        }
        None
    }

    fn restore(
        &self,
        p: &Problem,
        mut values: Vec<i32>,
        deadline: Option<Instant>,
    ) -> Option<Vec<f64>> {
        for row in &self.rows {
            if expired(deadline) {
                return None;
            }
            // 恒满足的行未维护增量活动度，统一按共享列的最终赋值重算。
            let sum = row.terms.iter().map(|&(j, a)| values[j] * a).sum();
            let (d, q) = row.nearest(sum);
            if d != 0 {
                return None;
            }
            let mut residual = row.rhs - sum;
            if let Some(c) = row.cycles {
                values[c.column] = q;
                residual -= q * c.coefficient;
            }
            if let Some((j, c)) = row.local {
                values[j] = c * residual;
            }
        }
        let x: Vec<_> = values.into_iter().map(f64::from).collect();
        if x.iter()
            .zip(&p.vars)
            .any(|(&v, var)| v < var.lb || v > var.ub)
            || p.constraints.iter().any(|r| r.expr.eval(&x) != r.rhs)
            || !(p.obj.eval(&x) + p.obj_offset).is_finite()
            || expired(deadline)
        {
            return None;
        }
        Some(x)
    }
}

/// 返回经过原始等式和变量界复核的整数点，只供应可行解。
pub(crate) fn seed(p: &Problem, deadline: Option<Instant>) -> Option<Vec<f64>> {
    if expired(deadline) {
        return None;
    }
    let mut model = Repair::detect(p, deadline)?;
    model.find(p, deadline)?;
    let values = std::mem::take(&mut model.values);
    model.restore(p, values, deadline)
}

/// 固定刚修复的周期数，再求一次线性规划来改善连续列的费用。这个子问题的
/// 对偶界只适用于固定后的模型，调用方只接收按原始等式重新核验的可行点。
pub(crate) fn polish(
    p: &Problem,
    seed: &[f64],
    opts: &crate::solver::SolveOptions,
    deadline: Option<Instant>,
    iterations: &mut u64,
) -> Option<Vec<f64>> {
    if expired(deadline) {
        return None;
    }
    let mut fixed = p.clone();
    for (var, &x) in fixed.vars.iter_mut().zip(seed) {
        if var.var_type != VarType::Continuous {
            var.lb = x;
            var.ub = x;
            var.var_type = VarType::Continuous;
        }
    }
    let (work, map) = crate::presolve::presolve(&fixed, &opts.tol).ok()?;
    if expired(deadline) {
        return None;
    }
    let mut lp = crate::simplex_sparse::new_engine(&work, opts.tol, opts.bb.lp_engine).ok()?;
    lp.set_time_deadline(deadline);
    let sol = lp.solve(&crate::simplex::Bounds::from_problem(&work));
    *iterations += lp.iterations();
    if sol.status != crate::simplex::LpStatus::Optimal || expired(deadline) {
        return None;
    }
    let x = crate::presolve::postsolve(&map, &sol.x, p.vars.len(), &opts.tol);
    if x.iter().zip(&p.vars).any(|(&v, var)| {
        !v.is_finite()
            || v < var.lb
            || v > var.ub
            || (var.var_type != VarType::Continuous && v != v.round())
    }) || p.constraints.iter().any(|r| r.expr.eval(&x) != r.rhs)
        || !(p.obj.eval(&x) + p.obj_offset).is_finite()
    {
        return None;
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(sign: f64) -> Problem {
        let mut p = Problem::new("周期补差", Sense::Minimize);
        for (name, hi) in [("共享", 2.0), ("首行补差", 1.0), ("次行补差", 2.0)] {
            let j = p.add_var(name, VarType::Continuous).unwrap();
            p.set_var_bounds(j, 0.0, hi).unwrap();
        }
        p.add_var("周数", VarType::Binary).unwrap();
        p.add_constraint(
            None,
            vec![(0, 1.0), (1, 1.0), (3, sign * 5.0)],
            RelOp::Eq,
            3.0 + sign * 5.0,
        )
        .unwrap();
        p.add_constraint(None, vec![(0, 1.0), (2, 1.0)], RelOp::Eq, 2.0)
            .unwrap();
        p.set_objective(vec![(0, 1.0)], 0.5).unwrap();
        p
    }

    #[test]
    fn private_period_of_either_sign_and_local_slack_are_restored() {
        for sign in [1.0, -1.0] {
            let p = sample(sign);
            let x = seed(&p, None).unwrap();
            assert_eq!(x, vec![2.0, 1.0, 0.0, 1.0]);
            assert_eq!(seed(&p, None), Some(x.clone()));
            assert!(p.constraints.iter().all(|r| r.expr.eval(&x) == r.rhs));
            assert_eq!(p.obj.eval(&x) + p.obj_offset, 2.5);
        }
    }

    #[test]
    fn shared_integer_mismatched_row_and_elapsed_time_decline() {
        let mut p = sample(-1.0);
        assert!(seed(&p, Some(Instant::now())).is_none());
        p.constraints[1].expr.add_term(3, -5.0);
        assert!(seed(&p, None).is_none());
        p.constraints[1].expr.add_term(3, 5.0);
        p.constraints[1].rhs = 3.25;
        assert!(seed(&p, None).is_none());
        p.constraints[1].rhs = f64::NAN;
        assert!(seed(&p, None).is_none());
    }
}

#[path = "periodic_projection.rs"]
mod projection;
pub(crate) use projection::improve;
