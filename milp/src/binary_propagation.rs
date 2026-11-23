//! 单线程自动启用的二元静态行传播。只使用工作模型静态行，不读取割。
//! 默认关闭的 certified_row_fallback 可携带同一缩减模型强化前的认证行；
//! 新行认证失败才回退，不改变精确格点、整数算术或容差门槛。
//!
//! 每个认证行转成正权字面量之和不超过整数右端；已为真的字面量之和
//! 超过右端即矛盾，加上一个未定字面量的权重超过右端即强迫它为假。
//! 位级恢复共同二进制格点，且格距严格大于解验收容差，避免把浮点误差
//! 或验收允许的微小违反变成不可行证明。所有活动度都是 checked i128。

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::time::Instant;

use crate::model::{dyadic, Constraint, Problem, RelOp, Tolerances, VarType, Variable};
use crate::simplex::Bounds;

// 测试在同一生产传播入口精确触发限时；线程隔离、离开作用域恢复旧状态。
#[cfg(test)]
std::thread_local! {
    static TEST_DEADLINE: std::cell::Cell<Option<(usize, usize)>> = const { std::cell::Cell::new(None) };
}

#[cfg(test)]
pub(crate) struct TestDeadlineGuard(Option<(usize, usize)>);

#[cfg(all(test, feature = "certified_row_fallback"))]
std::thread_local! {
    static CAPTURE_ATTEMPTS: std::cell::Cell<Option<usize>> = const { std::cell::Cell::new(None) };
}

/// 只供门控测试观察实际认证尝试；不会使根闭合时的传播统计提前变成非零。
#[cfg(all(test, feature = "certified_row_fallback"))]
pub(crate) struct CaptureAudit(Option<usize>);

#[cfg(all(test, feature = "certified_row_fallback"))]
impl CaptureAudit {
    pub(crate) fn new() -> Self {
        Self(CAPTURE_ATTEMPTS.with(|state| state.replace(Some(0))))
    }
    pub(crate) fn count(&self) -> usize {
        CAPTURE_ATTEMPTS.with(|state| state.get().unwrap_or(0))
    }
}

#[cfg(all(test, feature = "certified_row_fallback"))]
impl Drop for CaptureAudit {
    fn drop(&mut self) {
        CAPTURE_ATTEMPTS.with(|state| state.set(self.0));
    }
}

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

/// 只在本次 presolve→求解链中移动的旧认证行；不实现 Clone，不公开构造字段。
/// 原行号属于最终缩减模型，不能拿原始输入的列号或节点局部行填充。
#[derive(Default)]
pub(crate) struct PreservedBinaryRows {
    context: Option<PreservedContext>,
    rows: Vec<PreservedRow>,
}

struct PreservedRow {
    index: usize,
    certified: Vec<Row>,
    strengthened: Option<RowWitness>,
}

/// 只为保存的旧行记录强化后稀疏见证；逐位核对，不以地址或哈希代替证明来源。
struct RowWitness {
    op: RelOp,
    rhs_bits: u64,
    terms: Vec<(usize, u64)>,
}

impl RowWitness {
    #[cfg(feature = "certified_row_fallback")]
    fn of(row: &Constraint) -> Self {
        Self { op: row.op, rhs_bits: row.rhs.to_bits(),
            terms: row.expr.terms.iter().map(|&(j, a)| (j, a.to_bits())).collect() }
    }

    fn matches(&self, row: &Constraint) -> bool {
        self.op == row.op && self.rhs_bits == row.rhs.to_bits()
            && self.terms.len() == row.expr.terms.len()
            && self.terms.iter().zip(&row.expr.terms)
                .all(|(&(j, bits), &(k, a))| j == k && bits == a.to_bits())
    }
}

#[derive(PartialEq, Eq)]
struct PreservedContext {
    num_vars: usize,
    num_rows: usize,
    tolerance_bits: [u64; 5],
}

impl PreservedContext {
    fn of(problem: &Problem, tol: Tolerances) -> Self {
        Self { num_vars: problem.vars.len(), num_rows: problem.constraints.len(),
            tolerance_bits: tolerance_bits(tol) }
    }
}

fn tolerance_bits(tol: Tolerances) -> [u64; 5] {
    [tol.feasibility.to_bits(), tol.optimality.to_bits(), tol.pivot.to_bits(),
        tol.integrality.to_bits(), tol.zero.to_bits()]
}

impl PreservedBinaryRows {
    #[cfg(feature = "certified_row_fallback")]
    pub(crate) fn for_problem(problem: &Problem, tol: Tolerances) -> Self {
        Self { context: Some(PreservedContext::of(problem, tol)), rows: Vec::new() }
    }

    /// 调用方已确认本行即将强化；保存之前仍须通过唯一的生产认证函数。
    #[cfg(feature = "certified_row_fallback")]
    pub(crate) fn capture(&mut self, index: usize, vars: &[Variable], row: &Constraint, tol: Tolerances) -> bool {
        #[cfg(test)]
        CAPTURE_ATTEMPTS.with(|state| {
            if let Some(count) = state.get() { state.set(Some(count + 1)); }
        });
        if row.op == RelOp::Eq
            || !self.context.as_ref().is_some_and(|context| context.num_vars == vars.len()
            && index < context.num_rows && context.tolerance_bits == tolerance_bits(tol))
            || self.rows.last().is_some_and(|last| last.index >= index)
        {
            return false;
        }
        let Some(certified) = certify(vars, row, tol) else { return false; };
        self.rows.push(PreservedRow { index, certified, strengthened: None });
        true
    }

    /// 同一行完成系数强化后封存见证；已有见证不能被重写以掩盖模型变化。
    #[cfg(feature = "certified_row_fallback")]
    pub(crate) fn seal_last(&mut self, index: usize, row: &Constraint) {
        if let Some(last) = self.rows.last_mut().filter(|last| last.index == index
            && last.strengthened.is_none()) {
            last.strengthened = Some(RowWitness::of(row));
        }
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// 入口只读校验，不构造邻接表；容差与最终行列空间必须与捕获时一致。
    /// 载荷由私有 presolve 入口产生并按值移交，禁止跨模型缓存或复用。
    pub(crate) fn validate(&self, problem: &Problem, tol: Tolerances) -> bool {
        if self.rows.is_empty() {
            return true;
        }
        if self.context != Some(PreservedContext::of(problem, tol)) {
            return false;
        }
        let mut previous = None;
        for saved in &self.rows {
            let Some(source) = problem.constraints.get(saved.index) else {
                return false;
            };
            let certified = &saved.certified;
            if previous.is_some_and(|old| old >= saved.index) || certified.len() != 1
                || !saved.strengthened.as_ref().is_some_and(|witness| witness.matches(source))
                || source.op == RelOp::Eq
                || certified[0].literals.len() != source.expr.terms.len()
                || certified[0].literals.iter().zip(&source.expr.terms).any(|(literal, &(j, _))| {
                    literal.var != j || !problem.vars.get(j).is_some_and(|v| {
                        v.var_type == VarType::Binary && standard_bounds(v.lb, v.ub)
                    })
                })
            {
                return false;
            }
            previous = Some(saved.index);
        }
        true
    }

    #[cfg(all(test, feature = "certified_row_fallback"))]
    pub(crate) fn len(&self) -> usize {
        self.rows.len()
    }
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

pub(crate) fn standard_bounds(lb: f64, ub: f64) -> bool {
    (lb == 0.0 || lb == 1.0) && (ub == 0.0 || ub == 1.0) && lb <= ub
}

fn certify(vars: &[Variable], row: &Constraint, tol: Tolerances) -> Option<Vec<Row>> {
    certify_binary_terms(&row.expr.terms, row.op, row.rhs, tol, |j| {
        vars.get(j).is_some_and(|v| {
            v.var_type == VarType::Binary && standard_bounds(v.lb, v.ub)
        })
    })
}

fn certify_binary_terms(
    terms: &[(usize, f64)],
    op: RelOp,
    rhs: f64,
    tol: Tolerances,
    is_binary: impl Fn(usize) -> bool,
) -> Option<Vec<Row>> {
    if !tol.zero.is_finite()
        || tol.zero < 0.0
        || !tol.feasibility.is_finite()
        || tol.feasibility < 0.0
        || !tol.integrality.is_finite()
        || tol.integrality < 0.0
    {
        return None;
    }
    let mut exponent = if rhs == 0.0 {
        i32::MAX
    } else {
        dyadic(rhs)?.1
    };
    let mut max_coefficient = 0.0_f64;
    let mut previous = None;
    for &(j, coefficient) in terms {
        if previous.is_some_and(|old| old >= j)
            || !is_binary(j)
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
    let rhs = on_grid(rhs, exponent)?;
    let mut magnitude = rhs.checked_abs()?;
    let mut norm = 0_i128;
    let mut coefficients = Vec::with_capacity(terms.len());
    for &(j, coefficient) in terms {
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
    match op {
        RelOp::Le => Some(vec![make(1)?]),
        RelOp::Ge => Some(vec![make(-1)?]),
        RelOp::Eq => Some(vec![make(1)?, make(-1)?]),
    }
}

/// 用原有格点和容差认证一组正字面量。两项最小权重之和已经超过右端时，
/// 这一组内任意两个二元变量都不能同时取 1。每组只记录一份成员索引。
#[cfg(feature = "binary_mutex_row_fixing")]
pub(crate) fn certifies_binary_mutex(
    terms: &[(usize, f64)],
    op: RelOp,
    rhs: f64,
    tol: Tolerances,
    is_binary: impl Fn(usize) -> bool,
) -> bool {
    if terms.len() < 2 {
        return false;
    }
    let Some(certified) = certify_binary_terms(terms, op, rhs, tol, is_binary) else {
        return false;
    };
    certified.iter().any(|row| {
        if row.rhs < 0 || row.literals.iter().any(|literal| !literal.positive) {
            return false;
        }
        let (mut smallest, mut second) = (i128::MAX, i128::MAX);
        for literal in &row.literals {
            if literal.weight < smallest {
                second = smallest;
                smallest = literal.weight;
            } else {
                second = second.min(literal.weight);
            }
        }
        smallest.checked_add(second).is_some_and(|sum| sum > row.rhs)
    })
}

impl BinaryPropagation {
    pub(crate) fn new(problem: &Problem, tol: Tolerances) -> Self {
        Self::build(problem, tol, PreservedBinaryRows::default())
    }

    /// 首次孩子展开时消费载荷；新行先认证，只有失败才在同一位置采用旧行。
    pub(crate) fn with_preserved(problem: &Problem, tol: Tolerances, preserved: PreservedBinaryRows) -> Option<Self> {
        if preserved.is_empty() {
            return Some(Self::new(problem, tol));
        }
        if !preserved.validate(problem, tol) { return None; }
        Some(Self::build(problem, tol, preserved))
    }

    fn build(problem: &Problem, tol: Tolerances, preserved: PreservedBinaryRows) -> Self {
        let mut rows = Vec::new();
        let mut stats = PropagationStats::default();
        let mut old_rows = preserved.rows.into_iter().peekable();
        for i in 0..problem.constraints.len() {
            let old = if old_rows.peek().is_some_and(|saved| saved.index == i) {
                old_rows.next().map(|saved| saved.certified)
            } else { None };
            // 两个来源最多采用一个；被救回的行不再计 skipped，不产生重复行。
            if let Some(mut certified) = certify(&problem.vars, &problem.constraints[i], tol).or(old) {
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

    #[cfg(feature = "certified_row_fallback")]
    fn preserve_first(p: &Problem, tol: Tolerances) -> PreservedBinaryRows {
        let mut saved = PreservedBinaryRows::for_problem(p, tol);
        saved.capture(0, &p.vars, &p.constraints[0], tol);
        saved
    }

    #[cfg(feature = "certified_row_fallback")]
    fn signature(state: &BinaryPropagation) -> Vec<(i128, Vec<(usize, bool, i128)>)> {
        state.rows.iter().map(|row| (row.rhs, row.literals.iter()
            .map(|literal| (literal.var, literal.positive, literal.weight)).collect())).collect()
    }

    #[cfg(feature = "certified_row_fallback")]
    #[test]
    fn certified_fallback_prefers_new_then_old_and_counts_rejection_once() {
        for op in [RelOp::Le, RelOp::Ge] {
            let sign = if op == RelOp::Le { 1.0 } else { -1.0 };
            for tol in [Tolerances::default(), Tolerances { feasibility: 1e-8,
                integrality: 1e-7, ..Tolerances::default() }] {
                // 两个精确行的二元可行集相同；权重不同可检验实际采用了新行。
                let mut p = model(&[-10.0 * sign, 2.0 * sign, sign], op, 0.0);
                let mut saved = preserve_first(&p, tol);
                assert_eq!(saved.len(), 1);
                let old_signature = signature(&BinaryPropagation::new(&p, tol));
                p.constraints[0].expr.terms[0].1 = -3.0 * sign;
                saved.seal_last(0, &p.constraints[0]);
                let new_signature = signature(&BinaryPropagation::new(&p, tol));
                assert_ne!(old_signature, new_signature);
                let state = BinaryPropagation::with_preserved(&p, tol, saved).unwrap();
                assert_eq!(signature(&state), new_signature);
                assert_eq!((state.stats.accepted_rows, state.stats.skipped_rows), (1, 0));

                let mut p = model(&[-10.0 * sign, 2.0 * sign, sign], op, 0.0);
                let mut saved = preserve_first(&p, tol);
                p.constraints[0].expr.terms[0].1 = -f64::from_bits(3.0_f64.to_bits() + 1) * sign;
                saved.seal_last(0, &p.constraints[0]);
                assert_eq!(BinaryPropagation::new(&p, tol).stats.accepted_rows, 0);
                let state = BinaryPropagation::with_preserved(&p, tol, saved).unwrap();
                assert_eq!(signature(&state), old_signature);
                assert_eq!((state.stats.accepted_rows, state.stats.skipped_rows), (1, 0));

                let mut p = model(&[-10.1 * sign, 2.0 * sign, sign], op, 0.0);
                let mut saved = preserve_first(&p, tol);
                assert_eq!(saved.len(), 0, "不能把未认证旧行作为回退来源");
                p.constraints[0].expr.terms[0].1 = -f64::from_bits(3.0_f64.to_bits() + 1) * sign;
                saved.seal_last(0, &p.constraints[0]);
                let state = BinaryPropagation::with_preserved(&p, tol, saved).unwrap();
                assert_eq!((state.stats.accepted_rows, state.stats.skipped_rows), (0, 1));
            }
        }
        // 前面等式贡献两条，回退仍按源行号插入，不能误把数据库序号当源行号。
        let tol = Tolerances::default();
        let mut p = model(&[1.0, -1.0, 0.0], RelOp::Eq, 0.0);
        p.add_constraint(None, vec![(0, -10.0), (1, 2.0), (2, 1.0)], RelOp::Le, 0.0).unwrap();
        let original = signature(&BinaryPropagation::new(&p, tol));
        let mut saved = PreservedBinaryRows::for_problem(&p, tol);
        assert!(saved.capture(1, &p.vars, &p.constraints[1], tol));
        p.constraints[1].expr.terms[0].1 = -f64::from_bits(3.0_f64.to_bits() + 1);
        saved.seal_last(1, &p.constraints[1]);
        let state = BinaryPropagation::with_preserved(&p, tol, saved).unwrap();
        assert_eq!(signature(&state), original);
        assert_eq!((state.stats.accepted_rows, state.stats.skipped_rows), (3, 0));
    }

    #[cfg(feature = "certified_row_fallback")]
    #[test]
    fn certified_fallback_rejects_foreign_model_tolerance_and_column_mapping() {
        let tol = Tolerances::default();
        let mut p = model(&[-10.0, 2.0, 1.0], RelOp::Le, 0.0);
        let mut saved = preserve_first(&p, tol);
        assert!(!saved.validate(&p, tol), "未封存强化后见证不能放行");
        saved.seal_last(0, &p.constraints[0]);
        assert!(saved.validate(&p, tol));
        let mut other = p.clone();
        other.constraints[0].rhs = 1.0;
        assert!(!saved.validate(&other, tol), "同维模型也必须逐位匹配已封存的行");
        p.constraints[0].expr.terms[0].1 = -9.0;
        saved.seal_last(0, &p.constraints[0]);
        assert!(!saved.validate(&p, tol), "原位改系数不能靠重封存洗掉差异");
        p.constraints[0].expr.terms[0].1 = -10.0;
        assert!(!saved.validate(&p, Tolerances { integrality: 1e-5, ..tol }));
        p.constraints[0].expr.terms.swap(0, 1);
        assert!(!saved.validate(&p, tol));
        p.constraints[0].expr.terms.swap(0, 1);
        p.vars[0].var_type = VarType::Integer;
        assert!(!saved.validate(&p, tol));
        p.vars[0].var_type = VarType::Binary;
        p.vars[0].ub = 0.5;
        assert!(!saved.validate(&p, tol));
        p.vars[0].ub = 1.0;
        p.constraints[0].op = RelOp::Ge;
        assert!(!saved.validate(&p, tol));
        p.constraints[0].op = RelOp::Le;
        let moved = p;
        assert!(saved.validate(&moved, tol), "按值移动不应破坏已封存见证");
        let mut invalid_index = preserve_first(&moved, tol);
        invalid_index.seal_last(0, &moved.constraints[0]);
        invalid_index.rows[0].index = 1;
        assert!(!invalid_index.validate(&moved, tol));
        let mut invalid_dimension = moved;
        invalid_dimension.vars.pop();
        assert!(!saved.validate(&invalid_dimension, tol));
    }

    #[cfg(feature = "certified_row_fallback")]
    #[test]
    fn certified_fallback_preserves_all_small_domains_and_sibling_rollback() {
        let tol = Tolerances::default();
        let mut p = model(&[-10.0, 2.0, 1.0], RelOp::Le, 0.0);
        let original = p.clone();
        let mut saved = preserve_first(&p, tol);
        p.constraints[0].expr.terms[0].1 = -f64::from_bits(3.0_f64.to_bits() + 1);
        saved.seal_last(0, &p.constraints[0]);
        let mut state = BinaryPropagation::with_preserved(&p, tol, saved).unwrap();
        let mut fixes = 0;
        let mut conflicts = 0;
        for domain in 0..27 {
            let mut digits = domain;
            let mut bounds = Bounds::from_problem(&p);
            for j in 0..3 {
                match digits % 3 { 1 => bounds.ub[j] = 0.0, 2 => bounds.lb[j] = 1.0, _ => {} }
                digits /= 3;
            }
            let before = bounds.clone();
            let completions: Vec<_> = (0..8).filter(|&mask| feasible(&original, &bounds, mask)).collect();
            let mut trail = Vec::new();
            let outcome = state.propagate(&mut bounds, true, [], None, |j, lb, old| trail.push((j, lb, old)));
            if outcome == Outcome::Conflict { conflicts += 1; assert!(completions.is_empty()); }
            for mask in completions { assert!(feasible(&original, &bounds, mask)); }
            fixes += trail.len();
            for (j, is_lb, old) in trail.into_iter().rev() {
                if is_lb { bounds.lb[j] = old; } else { bounds.ub[j] = old; }
            }
            assert_eq!(bounds.lb, before.lb);
            assert_eq!(bounds.ub, before.ub);
            assert!(state.queue.is_empty() && state.queued.iter().all(|&queued| !queued));
        }
        assert!(fixes > 0 && conflicts > 0);
    }

    #[cfg(feature = "certified_row_fallback")]
    #[test]
    fn certified_fallback_deadline_after_partial_fixes_cleans_its_queue() {
        let tol = Tolerances::default();
        let mut coefficients = vec![1.0; 301];
        coefficients[0] = -1000.0;
        let mut p = model(&coefficients, RelOp::Le, 0.0);
        let mut saved = preserve_first(&p, tol);
        p.constraints[0].expr.terms[0].1 = -f64::from_bits(300.0_f64.to_bits() + 1);
        saved.seal_last(0, &p.constraints[0]);
        assert_eq!(BinaryPropagation::new(&p, tol).stats.accepted_rows, 0);
        let mut state = BinaryPropagation::with_preserved(&p, tol, saved).unwrap();
        let mut bounds = Bounds::from_problem(&p);
        bounds.ub[0] = 0.0;
        let mut checks = 0;
        let mut trail = Vec::new();
        let outcome = state.propagate_until(&mut bounds, true, [], || {
            checks += 1; checks == 3
        }, |j, is_lb, old| trail.push((j, is_lb, old)));
        assert_eq!(outcome, Outcome::Deadline);
        assert!(state.stats.fixes > 0 && state.stats.fixes < 300);
        assert_eq!(state.stats.fixes as usize, trail.len());
        assert_eq!(state.stats.conflicts, 0);
        assert!(state.queue.is_empty() && state.queued.iter().all(|&queued| !queued));
        assert_eq!(state.propagate(&mut bounds, true, [], None, |_, _, _| {}), Outcome::Complete);
        assert!(bounds.ub.iter().all(|&ub| ub == 0.0));
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
