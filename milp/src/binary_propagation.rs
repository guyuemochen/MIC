//! 默认关闭的二元静态行传播。只使用传入分支定界的工作模型，不读取割。
//!
//! 每个认证行转成正权字面量之和不超过整数右端；已为真的字面量之和
//! 超过右端即矛盾，加上一个未定字面量的权重超过右端即强迫它为假。
//! 位级恢复共同二进制格点，且格距严格大于解验收容差，避免把浮点误差
//! 或验收允许的微小违反变成不可行证明。所有活动度都是 checked i128。

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::time::Instant;

use crate::model::{Problem, RelOp, Tolerances, VarType};
use crate::simplex::Bounds;

// 测试在同一生产传播入口精确触发限时；线程隔离、离开作用域恢复旧状态。
#[cfg(test)]
std::thread_local! {
    static TEST_DEADLINE: std::cell::Cell<Option<(usize, usize)>> = const { std::cell::Cell::new(None) };
}

#[cfg(test)]
pub(crate) struct TestDeadlineGuard(Option<(usize, usize)>);

#[cfg(test)]
impl TestDeadlineGuard {
    pub(crate) fn on_check(target: usize) -> Self {
        Self(TEST_DEADLINE.with(|state| state.replace(Some((target, 0)))))
    }
}

#[cfg(test)]
impl Drop for TestDeadlineGuard {
    fn drop(&mut self) {
        TEST_DEADLINE.with(|state| state.set(self.0));
    }
}

#[cfg(test)]
fn test_deadline_expired() -> bool {
    TEST_DEADLINE.with(|state| match state.get() {
        Some((target, count)) => {
            state.set(Some((target, count + 1)));
            count + 1 >= target
        }
        None => false,
    })
}

/// 整数绝对值总量封顶，使原尺度任意二元部分和保持精确可表示。
const MAX_MAGNITUDE: i128 = 1 << 50;

#[derive(Clone, Copy)]
struct Literal {
    var: usize,
    positive: bool,
    weight: i128,
}

struct Row {
    literals: Vec<Literal>,
    rhs: i128,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct PropagationStats {
    /// 接受的规范不等式数；等式可贡献两条。
    pub accepted_rows: u64,
    /// 未能精确认证的原始行数。
    pub skipped_rows: u64,
    /// 节点传播调用次数。
    pub calls: u64,
    /// 出队检查的行数。
    pub row_visits: u64,
    /// 活动度与固定扫描访问的项数。
    pub term_visits: u64,
    /// 推导的二元固定数。
    pub fixes: u64,
    /// 已证明矛盾的节点数。
    pub conflicts: u64,
}

#[derive(Debug, PartialEq, Eq)]
pub(crate) enum Outcome {
    Complete,
    Conflict,
    Deadline,
}

pub(crate) struct BinaryPropagation {
    rows: Vec<Row>,
    adjacent: Vec<Vec<usize>>,
    queue: BinaryHeap<Reverse<usize>>,
    queued: Vec<bool>,
    pub(crate) stats: PropagationStats,
}

/// 返回奇数尾数与二进制指数，不通过近似取整猜测格点。
fn dyadic(value: f64) -> Option<(i128, i32)> {
    if !value.is_finite() {
        return None;
    }
    if value == 0.0 {
        return Some((0, 0));
    }
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let mut mantissa = bits & ((1_u64 << 52) - 1);
    let mut power = if exponent == 0 {
        -1074
    } else {
        exponent - 1023 - 52
    };
    if exponent != 0 {
        mantissa |= 1_u64 << 52;
    }
    let zeros = mantissa.trailing_zeros();
    mantissa >>= zeros;
    power += zeros as i32;
    let signed = if bits >> 63 == 0 {
        mantissa as i128
    } else {
        -(mantissa as i128)
    };
    Some((signed, power))
}

fn on_grid(value: f64, exponent: i32) -> Option<i128> {
    let (mantissa, power) = dyadic(value)?;
    if mantissa == 0 {
        return Some(0);
    }
    let shift = u32::try_from(power.checked_sub(exponent)?).ok()?;
    if shift >= 127 {
        return None;
    }
    mantissa.checked_mul(1_i128.checked_shl(shift)?)
}

fn standard_bounds(lb: f64, ub: f64) -> bool {
    (lb == 0.0 || lb == 1.0) && (ub == 0.0 || ub == 1.0) && lb <= ub
}

fn certify(problem: &Problem, index: usize, tol: Tolerances) -> Option<Vec<Row>> {
    let row = &problem.constraints[index];
    if !tol.zero.is_finite()
        || tol.zero < 0.0
        || !tol.feasibility.is_finite()
        || tol.feasibility < 0.0
        || !tol.integrality.is_finite()
        || tol.integrality < 0.0
    {
        return None;
    }
    let mut exponent = if row.rhs == 0.0 {
        i32::MAX
    } else {
        dyadic(row.rhs)?.1
    };
    let mut max_coefficient = 0.0_f64;
    let mut previous = None;
    for &(j, coefficient) in &row.expr.terms {
        let var = problem.vars.get(j)?;
        if previous.is_some_and(|old| old >= j)
            || var.var_type != VarType::Binary
            || !standard_bounds(var.lb, var.ub)
            || coefficient == 0.0
        {
            return None;
        }
        previous = Some(j);
        exponent = exponent.min(dyadic(coefficient)?.1);
        max_coefficient = max_coefficient.max(coefficient.abs());
    }
    if max_coefficient <= tol.zero || !(-1022..=1023).contains(&exponent) {
        return None;
    }
    let rhs = on_grid(row.rhs, exponent)?;
    let mut magnitude = rhs.checked_abs()?;
    let mut norm = 0_i128;
    let mut coefficients = Vec::with_capacity(row.expr.terms.len());
    for &(j, coefficient) in &row.expr.terms {
        let weight = on_grid(coefficient, exponent)?;
        let absolute = weight.checked_abs()?;
        magnitude = magnitude.checked_add(absolute)?;
        norm = norm.checked_add(absolute)?;
        if magnitude > MAX_MAGNITUDE {
            return None;
        }
        coefficients.push((j, weight));
    }
    let grid = 2_f64.powi(exponent);
    let norm_original = grid * norm as f64;
    let allowance = 10.0 * tol.feasibility + tol.integrality * norm_original;
    if !(grid * magnitude as f64).is_finite() || !allowance.is_finite() || allowance >= grid {
        return None;
    }
    let make = |sign: i128| -> Option<Row> {
        let mut bound = rhs.checked_mul(sign)?;
        let mut literals = Vec::with_capacity(coefficients.len());
        for &(var, coefficient) in &coefficients {
            let signed = coefficient.checked_mul(sign)?;
            let weight = signed.checked_abs()?;
            if signed < 0 {
                bound = bound.checked_add(weight)?;
            }
            literals.push(Literal {
                var,
                positive: signed > 0,
                weight,
            });
        }
        Some(Row {
            literals,
            rhs: bound,
        })
    };
    match row.op {
        RelOp::Le => Some(vec![make(1)?]),
        RelOp::Ge => Some(vec![make(-1)?]),
        RelOp::Eq => Some(vec![make(1)?, make(-1)?]),
    }
}

impl BinaryPropagation {
    pub(crate) fn new(problem: &Problem, tol: Tolerances) -> Self {
        let mut rows = Vec::new();
        let mut stats = PropagationStats::default();
        for i in 0..problem.constraints.len() {
            if let Some(mut certified) = certify(problem, i, tol) {
                rows.append(&mut certified);
            } else {
                stats.skipped_rows += 1;
            }
        }
        stats.accepted_rows = rows.len() as u64;
        let mut adjacent = vec![Vec::new(); problem.vars.len()];
        for (i, row) in rows.iter().enumerate() {
            for literal in &row.literals {
                adjacent[literal.var].push(i);
            }
        }
        let queued = vec![false; rows.len()];
        Self {
            rows,
            adjacent,
            queue: BinaryHeap::new(),
            queued,
            stats,
        }
    }

    fn seed_var(&mut self, var: usize) {
        for &row in &self.adjacent[var] {
            if !self.queued[row] {
                self.queued[row] = true;
                self.queue.push(Reverse(row));
            }
        }
    }

    /// 队列仅在本次节点中有效；无跨节点活动度或闭包缓存。
    /// depth=1 全播种，其余从分支列与全部 RC 覆盖日志播种。每次固定的旧界
    /// 交给调用方记入 DFS trail，混合搜索则由已有完整域快照捕获。
    pub(crate) fn propagate(
        &mut self,
        bounds: &mut Bounds,
        all_rows: bool,
        seeds: impl IntoIterator<Item = usize>,
        deadline: Option<Instant>,
        record: impl FnMut(usize, bool, f64),
    ) -> Outcome {
        self.propagate_until(
            bounds,
            all_rows,
            seeds,
            || {
                #[cfg(test)]
                if test_deadline_expired() {
                    return true;
                }
                deadline.is_some_and(|end| Instant::now() >= end)
            },
            record,
        )
    }

    // 注入检查函数使测试能精确打断项扫描；正常路径只读已有真实deadline。
    fn propagate_until(
        &mut self,
        bounds: &mut Bounds,
        all_rows: bool,
        seeds: impl IntoIterator<Item = usize>,
        mut expired: impl FnMut() -> bool,
        mut record: impl FnMut(usize, bool, f64),
    ) -> Outcome {
        self.stats.calls += 1;
        if all_rows {
            for row in 0..self.rows.len() {
                self.queue.push(Reverse(row));
                self.queued[row] = true;
            }
        } else {
            for var in seeds {
                self.seed_var(var);
            }
        }
        let result = self.run(bounds, &mut expired, &mut record);
        // 矛盾和限时同样清理去重位，下一兄弟不能继承未处理事件。
        while let Some(Reverse(row)) = self.queue.pop() {
            self.queued[row] = false;
        }
        result
    }

    fn run(
        &mut self,
        bounds: &mut Bounds,
        expired: &mut impl FnMut() -> bool,
        record: &mut impl FnMut(usize, bool, f64),
    ) -> Outcome {
        if expired() {
            return Outcome::Deadline;
        }
        while let Some(Reverse(index)) = self.queue.pop() {
            self.queued[index] = false;
            self.stats.row_visits += 1;
            let row = &self.rows[index];
            let mut sum = 0_i128;
            let mut valid = true;
            for literal in &row.literals {
                self.stats.term_visits += 1;
                if self.stats.term_visits % 256 == 0 && expired() {
                    return Outcome::Deadline;
                }
                let (lb, ub) = (bounds.lb[literal.var], bounds.ub[literal.var]);
                if !standard_bounds(lb, ub) {
                    valid = false;
                    break;
                }
                if lb == ub && ((lb == 1.0) == literal.positive) {
                    let Some(next) = sum.checked_add(literal.weight) else {
                        valid = false;
                        break;
                    };
                    sum = next;
                }
            }
            if !valid {
                continue;
            }
            if sum > row.rhs {
                self.stats.conflicts += 1;
                return Outcome::Conflict;
            }
            // 当前行强迫的字面量均为假，不增加其活动度；可沿用同一 sum。
            for k in 0..self.rows[index].literals.len() {
                self.stats.term_visits += 1;
                if self.stats.term_visits % 256 == 0 && expired() {
                    return Outcome::Deadline;
                }
                let literal = self.rows[index].literals[k];
                let var = literal.var;
                if bounds.lb[var] == bounds.ub[var] {
                    continue;
                }
                let Some(with_literal) = sum.checked_add(literal.weight) else {
                    continue;
                };
                if with_literal <= self.rows[index].rhs {
                    continue;
                }
                let is_lb = !literal.positive;
                if is_lb {
                    record(var, true, bounds.lb[var]);
                    bounds.lb[var] = 1.0;
                } else {
                    record(var, false, bounds.ub[var]);
                    bounds.ub[var] = 0.0;
                }
                self.stats.fixes += 1;
                self.seed_var(var);
            }
        }
        if expired() {
            Outcome::Deadline
        } else {
            Outcome::Complete
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Constraint, LinearExpr, Sense};

    fn model(coefficients: &[f64], op: RelOp, rhs: f64) -> Problem {
        let mut p = Problem::new("精确行传播", Sense::Maximize);
        for j in 0..coefficients.len() {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        p.add_constraint(
            None,
            coefficients.iter().copied().enumerate().collect(),
            op,
            rhs,
        )
        .unwrap();
        p
    }

    fn feasible(p: &Problem, bounds: &Bounds, mask: usize) -> bool {
        let x: Vec<_> = (0..p.vars.len())
            .map(|j| ((mask >> j) & 1) as f64)
            .collect();
        if x.iter()
            .enumerate()
            .any(|(j, &v)| v < bounds.lb[j] || v > bounds.ub[j])
        {
            return false;
        }
        p.constraints.iter().all(|row| {
            let lhs = row.expr.eval(&x);
            match row.op {
                RelOp::Le => lhs <= row.rhs,
                RelOp::Ge => lhs >= row.rhs,
                RelOp::Eq => lhs == row.rhs,
            }
        })
    }

    #[test]
    fn every_small_domain_preserves_all_feasible_completions() {
        // 六变量729种三态域；正负权、负右端、三种行方向全部对二元真值枚举。
        for op in [RelOp::Le, RelOp::Ge, RelOp::Eq] {
            for rhs in -7..=7 {
                let p = model(&[1.0, -2.0, 3.0, -1.0, 2.0, -3.0], op, rhs as f64);
                let mut propagation = BinaryPropagation::new(&p, Tolerances::default());
                assert!(propagation.stats.accepted_rows > 0);
                for domain in 0..729 {
                    let mut bounds = Bounds::from_problem(&p);
                    let mut digits = domain;
                    for j in 0..6 {
                        match digits % 3 {
                            1 => bounds.ub[j] = 0.0,
                            2 => bounds.lb[j] = 1.0,
                            _ => {}
                        }
                        digits /= 3;
                    }
                    let completions: Vec<_> = (0..64)
                        .filter(|&mask| feasible(&p, &bounds, mask))
                        .collect();
                    let result = propagation.propagate(&mut bounds, true, [], None, |_, _, _| {});
                    assert_ne!(result, Outcome::Deadline);
                    if result == Outcome::Conflict {
                        assert!(completions.is_empty());
                    }
                    for mask in completions {
                        assert!(
                            feasible(&p, &bounds, mask),
                            "op={op:?} rhs={rhs} domain={domain} mask={mask}"
                        );
                    }
                    if result == Outcome::Complete {
                        let fixed = propagation.stats.fixes;
                        assert_eq!(
                            propagation.propagate(&mut bounds, true, [], None, |_, _, _| {}),
                            Outcome::Complete
                        );
                        assert_eq!(fixed, propagation.stats.fixes, "返回完成必须已闭包");
                    }
                }
            }
        }
    }

    #[test]
    fn grid_certificate_rejects_uncertain_or_invalid_rows() {
        let tol = Tolerances::default();
        for scale in [0.25, 1.0, 8.0, 1024.0] {
            let p = model(&[scale, -2.0 * scale], RelOp::Le, scale);
            assert_eq!(BinaryPropagation::new(&p, tol).stats.accepted_rows, 1);
        }
        for (coeff, rhs) in [
            (1.0, 1.0 - 1e-10),
            (0.1, 0.3),
            (f64::MIN_POSITIVE, 0.0),
            (f64::from_bits(1), 0.0),
            (1.0, 2_f64.powi(200)),
            (tol.zero, 0.0),
        ] {
            let p = model(&[coeff], RelOp::Le, rhs);
            assert_eq!(
                BinaryPropagation::new(&p, tol).stats.accepted_rows,
                0,
                "{coeff:?}, {rhs:?}"
            );
        }
        let p = model(&[1.0], RelOp::Le, 0.0);
        for feasibility in [0.1, 0.100001, f64::NAN, -1.0] {
            let custom = Tolerances {
                feasibility,
                integrality: 0.0,
                ..tol
            };
            assert_eq!(BinaryPropagation::new(&p, custom).stats.accepted_rows, 0);
        }
        for kind in [VarType::Integer, VarType::Continuous] {
            let mut bad = p.clone();
            bad.vars[0].var_type = kind;
            assert_eq!(BinaryPropagation::new(&bad, tol).stats.accepted_rows, 0);
        }
        for (lb, ub) in [(0.1, 1.0), (0.0, 2.0), (1.0, 0.0), (f64::NAN, 1.0)] {
            let mut bad = p.clone();
            bad.vars[0].lb = lb;
            bad.vars[0].ub = ub;
            assert_eq!(BinaryPropagation::new(&bad, tol).stats.accepted_rows, 0);
        }
        for (terms, rhs) in [
            (vec![(5, 1.0)], 0.0),
            (vec![(0, f64::INFINITY)], 0.0),
            (vec![(0, 1.0)], f64::NAN),
            (vec![(0, 1.0), (0, -1.0)], 0.0),
            (vec![(0, 2_f64.powi(200)), (1, 1.0)], 0.0),
        ] {
            let mut bad = model(&[1.0, 1.0], RelOp::Le, 0.0);
            bad.constraints = vec![Constraint {
                name: None,
                expr: LinearExpr { terms },
                op: RelOp::Le,
                rhs,
            }];
            assert_eq!(BinaryPropagation::new(&bad, tol).stats.accepted_rows, 0);
        }
    }

    #[test]
    fn cascade_rollback_conflict_queue_and_deadline_are_isolated() {
        let mut p = model(&[1.0, 1.0, 0.0, 0.0], RelOp::Le, 1.0);
        p.add_constraint(None, vec![(1, 1.0), (2, -1.0)], RelOp::Eq, 0.0)
            .unwrap();
        p.add_constraint(None, vec![(2, 1.0), (3, 1.0)], RelOp::Ge, 1.0)
            .unwrap();
        let mut propagation = BinaryPropagation::new(&p, Tolerances::default());
        let mut bounds = Bounds::from_problem(&p);
        bounds.lb[0] = 1.0;
        let parent = bounds.clone();
        let mut trail = Vec::new();
        assert_eq!(
            propagation.propagate(&mut bounds, true, [], None, |j, lb, old| trail
                .push((j, lb, old))),
            Outcome::Complete
        );
        assert_eq!(bounds.ub[1], 0.0);
        assert_eq!(bounds.ub[2], 0.0);
        assert_eq!(bounds.lb[3], 1.0);
        for (j, lb, old) in trail.into_iter().rev() {
            if lb {
                bounds.lb[j] = old;
            } else {
                bounds.ub[j] = old;
            }
        }
        assert_eq!(bounds.lb, parent.lb);
        assert_eq!(bounds.ub, parent.ub);
        bounds.lb[1] = 1.0;
        assert_eq!(
            propagation.propagate(&mut bounds, true, [], None, |_, _, _| {}),
            Outcome::Conflict
        );
        bounds = Bounds::from_problem(&p);
        assert_eq!(
            propagation.propagate(&mut bounds, true, [], Some(Instant::now()), |_, _, _| {}),
            Outcome::Deadline
        );
        assert_eq!(
            propagation.propagate(&mut bounds, true, [], None, |_, _, _| {}),
            Outcome::Complete
        );
        // 非标准运行时界只跳过相关行，不产生伪固定/伪冲突。
        let p = model(&[1.0], RelOp::Le, 0.0);
        let mut propagation = BinaryPropagation::new(&p, Tolerances::default());
        let mut bounds = Bounds::from_problem(&p);
        bounds.lb[0] = 0.5;
        assert_eq!(
            propagation.propagate(&mut bounds, true, [], None, |_, _, _| {}),
            Outcome::Complete
        );
        assert_eq!(bounds.lb[0], 0.5);
        assert_eq!(bounds.ub[0], 1.0);
    }
    #[test]
    fn deadline_mid_scan_preserves_only_proved_fixes_and_cleans_queue() {
        let p = model(&vec![1.0; 300], RelOp::Le, 0.0);
        let mut propagation = BinaryPropagation::new(&p, Tolerances::default());
        let mut bounds = Bounds::from_problem(&p);
        let mut checks = 0;
        let mut trail = Vec::new();
        let outcome = propagation.propagate_until(
            &mut bounds,
            true,
            [],
            || {
                checks += 1;
                checks == 3
            },
            |j, lb, old| trail.push((j, lb, old)),
        );
        assert_eq!(outcome, Outcome::Deadline);
        assert!(propagation.stats.fixes > 0 && propagation.stats.fixes < 300);
        assert_eq!(trail.len() as u64, propagation.stats.fixes);
        assert!(propagation.queue.is_empty() && propagation.queued.iter().all(|&v| !v));
        assert_eq!(
            propagation.propagate(&mut bounds, true, [], None, |_, _, _| {}),
            Outcome::Complete
        );
        assert!(bounds.ub.iter().all(|&v| v == 0.0));
    }
}
