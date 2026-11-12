//! 生产求解路径的只读证明进展诊断；仅 proof_profile 特性编译。
//! 不补跑LP、不改候选排序/学习/预算。调用方显式reset后采集，take后停止。
//! 本轮仅单线程；事件缺失表示未到达该阶段，不代表零成本或零界改善。
use crate::model::Sense;
use crate::simplex::{LpSolution, LpStatus};
use std::cell::RefCell;
use std::io::{self, Write};

#[derive(Clone, Debug)]
pub enum Value {
    Int(u64),
    Number(f64),
    Text(&'static str),
    String(String),
    Bool(bool),
    Null,
}
impl From<u64> for Value {
    fn from(v: u64) -> Self {
        Self::Int(v)
    }
}
impl From<usize> for Value {
    fn from(v: usize) -> Self {
        Self::Int(v as u64)
    }
}
impl From<f64> for Value {
    fn from(v: f64) -> Self {
        Self::Number(v)
    }
}
impl From<bool> for Value {
    fn from(v: bool) -> Self {
        Self::Bool(v)
    }
}
impl From<&'static str> for Value {
    fn from(v: &'static str) -> Self {
        Self::Text(v)
    }
}
pub fn optional(v: Option<f64>) -> Value {
    v.map_or(Value::Null, Value::Number)
}
pub fn status(v: LpStatus) -> &'static str {
    match v {
        LpStatus::Optimal => "optimal",
        LpStatus::Infeasible => "infeasible",
        LpStatus::Unbounded => "unbounded",
        LpStatus::IterationLimit => "iteration_limit",
        LpStatus::TimeLimit => "time_limit",
    }
}
pub fn bound(sol: &LpSolution) -> Value {
    if sol.status == LpStatus::Optimal {
        optional(Some(sol.objective).filter(|v| v.is_finite()))
    } else {
        Value::Null
    }
}
#[derive(Clone, Debug)]
pub struct Event {
    pub kind: &'static str,
    pub fields: Vec<(&'static str, Value)>,
}
#[derive(Clone, Debug, Default)]
pub struct ProbeTotals {
    pub calls: u64,
    pub iterations: u64,
    pub seconds: f64,
    pub positive: u64,
    pub zero: u64,
    pub negative: u64,
    pub nonfinite_optimal_gain: u64,
    pub infeasible: u64,
    pub iteration_limit: u64,
    pub time_limit: u64,
    pub unbounded: u64,
}
#[derive(Clone, Debug, Default)]
pub struct Trace {
    pub available: bool,
    pub requested_threads: usize,
    pub events: Vec<Event>,
    pub probes: [ProbeTotals; 9],
    pub chosen_omitted: u64,
    pub failures_omitted: u64,
    failure_samples: usize,
    pub tree_lp_calls: u64,
    pub tree_lp_iterations: u64,
    pub unresolved_lp_nodes: u64,
    pub unresolved_est: Option<f64>,
    pub processed_min_depth: Option<u32>,
    pub processed_max_depth: Option<u32>,
    probed_vars: Vec<usize>,
    enabled: bool,
    depth: Option<usize>,
    seen_frontiers: u64,
    best_samples: usize,
    chosen_samples: usize,
}
thread_local! { static TRACE:RefCell<Trace> = RefCell::new(Trace::default()); }
pub fn reset(requested_threads: usize) {
    TRACE.with(|trace| {
        *trace.borrow_mut() = Trace {
            requested_threads,
            available: requested_threads == 1,
            enabled: requested_threads == 1,
            ..Trace::default()
        }
    });
}
pub fn take() -> Trace {
    TRACE.with(|trace| std::mem::take(&mut *trace.borrow_mut()))
}
pub fn record(kind: &'static str, fields: Vec<(&'static str, Value)>) -> Option<usize> {
    TRACE.with(|trace| {
        let mut t = trace.borrow_mut();
        if !t.enabled {
            return None;
        }
        let id = t.events.len();
        t.events.push(Event { kind, fields });
        Some(id)
    })
}
pub fn append(id: Option<usize>, fields: Vec<(&'static str, Value)>) {
    if let Some(id) = id {
        TRACE.with(|trace| {
            if let Some(event) = trace.borrow_mut().events.get_mut(id) {
                event.fields.extend(fields);
            }
        });
    }
}
pub fn lp(
    kind: &'static str,
    sol: &LpSolution,
    iterations: u64,
    seconds: Option<f64>,
    retained_cuts: u64,
) {
    record(
        kind,
        vec![
            ("status", status(sol.status).into()),
            ("bound", bound(sol)),
            ("iterations", iterations.into()),
            (
                "iterations_scope",
                if kind == "root_after_cuts" {
                    "bb_cumulative_at_stage"
                } else {
                    "call"
                }
                .into(),
            ),
            // B&B既有计数会随cut_loop内部重建重置一段；不冒充实际全部LP成本。
            (
                "iterations_note",
                if kind == "root_after_cuts" {
                    "production_counter_may_reset_on_rebuild"
                } else {
                    "observed_call_delta"
                }
                .into(),
            ),
            ("seconds", optional(seconds)),
            (
                "seconds_scope",
                if kind == "root_after_cuts" {
                    "root_cut_stage_including_rebuild"
                } else {
                    "call"
                }
                .into(),
            ),
            ("retained_cuts", retained_cuts.into()),
        ],
    );
}

pub struct ProbeScope(Option<usize>);
impl ProbeScope {
    pub fn new(depth: usize) -> Self {
        Self(TRACE.with(|t| {
            let mut t = t.borrow_mut();
            let old = t.depth;
            t.depth = Some(depth);
            t.probed_vars.clear();
            old
        }))
    }
}
impl Drop for ProbeScope {
    fn drop(&mut self) {
        TRACE.with(|t| t.borrow_mut().depth = self.0);
    }
}
pub fn probe(
    sol: &LpSolution,
    sense: Sense,
    parent: f64,
    iterations: u64,
    seconds: f64,
    var: usize,
) {
    TRACE.with(|trace| {
        let mut t = trace.borrow_mut();
        if !t.enabled {
            return;
        }
        let Some(depth) = t.depth.filter(|&d| d < 9) else {
            return;
        };
        if !t.probed_vars.contains(&var) {
            t.probed_vars.push(var);
        }
        let p = &mut t.probes[depth];
        p.calls += 1;
        p.iterations += iterations;
        p.seconds += seconds;
        match sol.status {
            LpStatus::Optimal => {
                let gain = match sense {
                    Sense::Maximize => parent - sol.objective,
                    Sense::Minimize => sol.objective - parent,
                };
                if !gain.is_finite() {
                    p.nonfinite_optimal_gain += 1;
                } else if gain > 0.0 {
                    p.positive += 1;
                } else if gain < 0.0 {
                    p.negative += 1;
                } else {
                    p.zero += 1;
                }
            }
            LpStatus::Infeasible => p.infeasible += 1,
            LpStatus::IterationLimit => {
                p.iteration_limit += 1;
            }
            LpStatus::TimeLimit => p.time_limit += 1,
            LpStatus::Unbounded => p.unbounded += 1,
        }
    });
}
pub fn clear_probed() {
    TRACE.with(|t| t.borrow_mut().probed_vars.clear());
}
pub fn chosen(depth: u32, var: usize, observations: u32) {
    let probed = TRACE.with(|t| t.borrow().probed_vars.contains(&var));
    let sample = TRACE.with(|t| {
        let mut t = t.borrow_mut();
        if !t.enabled {
            return false;
        }
        if t.chosen_samples >= 32 {
            t.chosen_omitted += 1;
            false
        } else {
            t.chosen_samples += 1;
            true
        }
    });
    if sample {
        record(
            "chosen",
            vec![
                ("depth", (depth as u64).into()),
                ("var", var.into()),
                ("observations", (observations as u64).into()),
                ("probed_this_node", probed.into()),
            ],
        );
    }
}
pub fn node(
    depth: u32,
    var: usize,
    up: bool,
    est: f64,
    parent: f64,
    best_bound: bool,
) -> Option<usize> {
    let sample = TRACE.with(|t| {
        let mut t = t.borrow_mut();
        if !t.enabled {
            return false;
        }
        if best_bound && t.best_samples < 32 {
            t.best_samples += 1;
            true
        } else {
            depth == 1
        }
    });
    if !sample {
        return None;
    }
    record(
        "node",
        vec![
            ("depth", (depth as u64).into()),
            ("var", var.into()),
            ("up", up.into()),
            ("inherited_est", est.into()),
            ("parent_lp_bound", parent.into()),
            ("best_bound_selection", best_bound.into()),
        ],
    )
}
pub fn tree_lp(sol: &LpSolution, sense: Sense, est: f64, depth: u32, iterations: u64) {
    TRACE.with(|t| {
        let mut t = t.borrow_mut();
        if !t.enabled {
            return;
        }
        t.tree_lp_calls += 1;
        t.tree_lp_iterations += iterations;
        t.processed_min_depth = Some(t.processed_min_depth.map_or(depth, |old| old.min(depth)));
        t.processed_max_depth = Some(t.processed_max_depth.map_or(depth, |old| old.max(depth)));
        if matches!(sol.status, LpStatus::IterationLimit | LpStatus::TimeLimit) {
            t.unresolved_lp_nodes += 1;
            if est.is_finite() {
                t.unresolved_est = Some(t.unresolved_est.map_or(est, |old| match sense {
                    Sense::Maximize => old.max(est),
                    Sense::Minimize => old.min(est),
                }));
            }
        }
    });
}
pub fn node_outcome(id: Option<usize>, outcome: &'static str) {
    append(id, vec![("outcome", outcome.into())]);
}
pub fn node_lp(id: Option<usize>, sol: &LpSolution, iterations: u64) {
    append(
        id,
        vec![
            ("outcome", "lp".into()),
            ("status", status(sol.status).into()),
            ("actual_bound", bound(sol)),
            ("iterations", iterations.into()),
        ],
    );
}
/// 栈未弹出时采样。并列采用浮点==，只把有限est当界；退出样本总是保留。
pub fn frontier(
    sense: Sense,
    nodes: u64,
    exit: bool,
    incumbent: Option<f64>,
    unresolved_bound: Option<f64>,
    combined_bound: impl FnOnce() -> Option<f64>,
    entries: impl Iterator<Item = (f64, u32)>,
) {
    const POINTS: [u64; 7] = [1, 64, 128, 256, 512, 1024, 4096];
    let sample = TRACE.with(|t| {
        let mut t = t.borrow_mut();
        if !t.enabled {
            return false;
        }
        if exit {
            return true;
        }
        if let Some(i) = POINTS.iter().position(|&n| n == nodes) {
            let mask = 1 << i;
            if t.seen_frontiers & mask == 0 {
                t.seen_frontiers |= mask;
                return true;
            }
        }
        false
    });
    if !sample {
        return;
    }
    // 只有真正采样时才扫描生产前沿，不能让实参求值每个节点都付O(stack)。
    let combined_bound = combined_bound();
    let mut len = 0usize;
    let mut root_children = 0usize;
    let mut minimum = u32::MAX;
    let mut maximum = 0u32;
    let mut weakest: Option<f64> = None;
    let mut ties = 0usize;
    let mut tie_depth = u32::MAX;
    let mut unknown = 0usize;
    for (est, depth) in entries {
        len += 1;
        root_children += usize::from(depth == 1);
        minimum = minimum.min(depth);
        maximum = maximum.max(depth);
        if !est.is_finite() {
            unknown += 1;
            continue;
        }
        let better = weakest.is_none_or(|old| match sense {
            Sense::Maximize => est > old,
            Sense::Minimize => est < old,
        });
        if better {
            weakest = Some(est);
            ties = 1;
            tie_depth = depth;
        } else if weakest == Some(est) {
            ties += 1;
            tie_depth = tie_depth.min(depth);
        }
    }
    record(
        "frontier",
        vec![
            ("nodes", nodes.into()),
            ("exit", exit.into()),
            ("incumbent", optional(incumbent)),
            (
                "stack_bound",
                if unknown == 0 {
                    optional(weakest)
                } else {
                    Value::Null
                },
            ),
            (
                "unresolved_count",
                TRACE.with(|t| t.borrow().unresolved_lp_nodes).into(),
            ),
            ("unresolved_bound", optional(unresolved_bound)),
            ("combined_bound", optional(combined_bound)),
            ("size", len.into()),
            ("unknown_est_count", unknown.into()),
            ("weakest_ties", ties.into()),
            (
                "weakest_min_depth",
                if ties > 0 {
                    (tie_depth as u64).into()
                } else {
                    Value::Null
                },
            ),
            ("depth1_pending", root_children.into()),
            (
                "min_depth",
                if len > 0 {
                    (minimum as u64).into()
                } else {
                    Value::Null
                },
            ),
            (
                "max_depth",
                if len > 0 {
                    (maximum as u64).into()
                } else {
                    Value::Null
                },
            ),
        ],
    );
}
fn write_string(out: &mut impl Write, value: &str) -> io::Result<()> {
    write!(out, "\"")?;
    for ch in value.chars() {
        match ch {
            '"' => write!(out, "\\\"")?,
            '\\' => write!(out, "\\\\")?,
            '\u{0}'..='\u{1f}' => write!(out, "\\u{:04x}", u32::from(ch))?,
            _ => write!(out, "{ch}")?,
        }
    }
    write!(out, "\"")
}

pub fn lp_failure(
    stage: &'static str,
    sol: &LpSolution,
    engine: &dyn crate::simplex_sparse::LpEngine,
) {
    if sol.status == LpStatus::Optimal {
        return;
    }
    let sample = TRACE.with(|trace| {
        let mut t = trace.borrow_mut();
        if !t.enabled {
            return false;
        }
        if t.failure_samples >= 64 {
            t.failures_omitted += 1;
            false
        } else {
            t.failure_samples += 1;
            true
        }
    });
    if !sample {
        return;
    }
    let detail = engine.proof_diagnostics();
    record(
        "lp_failure",
        vec![
            ("stage", stage.into()),
            ("status", status(sol.status).into()),
            (
                "numerical_error",
                detail.numerical_error.map_or(Value::Null, Value::String),
            ),
            (
                "phase1_pending",
                detail.phase1_pending.map_or(Value::Null, Value::Bool),
            ),
        ],
    );
}

fn write_value(out: &mut impl Write, v: &Value) -> io::Result<()> {
    match v {
        Value::Int(n) => write!(out, "{n}"),
        Value::Number(n) if n.is_finite() => write!(out, "{n:?}"),
        Value::Text(s) => write_string(out, s),
        Value::String(s) => write_string(out, s),
        Value::Bool(v) => write!(out, "{v}"),
        _ => write!(out, "null"),
    }
}
impl Trace {
    pub fn write_json(&self, out: &mut impl Write) -> io::Result<()> {
        write!(
            out,
            "{{\"available\":{},\"requested_threads\":{},\"chosen_omitted\":{},\"events\":[",
            self.available, self.requested_threads, self.chosen_omitted
        )?;
        for (i, event) in self.events.iter().enumerate() {
            if i > 0 {
                write!(out, ",")?;
            }
            write!(out, "{{\"kind\":\"{}\"", event.kind)?;
            for (key, value) in &event.fields {
                write!(out, ",\"{key}\":")?;
                write_value(out, value)?;
            }
            write!(out, "}}")?;
        }
        write!(out, "],\"failures_omitted\":{},", self.failures_omitted)?;
        write!(out,"\"tree_lp_calls\":{},\"tree_lp_iterations\":{},\"unresolved_lp_nodes\":{},\"unresolved_est\":",self.tree_lp_calls,self.tree_lp_iterations,self.unresolved_lp_nodes)?;
        write_value(out, &optional(self.unresolved_est))?;
        write!(out, ",\"processed_min_depth\":")?;
        write_value(
            out,
            &self
                .processed_min_depth
                .map_or(Value::Null, |v| Value::Int(v as u64)),
        )?;
        write!(out, ",\"processed_max_depth\":")?;
        write_value(
            out,
            &self
                .processed_max_depth
                .map_or(Value::Null, |v| Value::Int(v as u64)),
        )?;
        write!(out, ",\"probes\":[")?;
        for (depth, p) in self.probes.iter().enumerate() {
            if depth > 0 {
                write!(out, ",")?;
            }
            write!(
                out,
                "{{\"depth\":{depth},\"calls\":{},\"iterations\":{},\"seconds\":",
                p.calls, p.iterations
            )?;
            write_value(out, &Value::Number(p.seconds))?;
            write!(out,",\"positive_gain\":{},\"zero_gain\":{},\"negative_gain\":{},\"nonfinite_optimal_gain\":{},\"infeasible\":{},\"iteration_limit\":{},\"time_limit\":{},\"unbounded\":{}}}",p.positive,p.zero,p.negative,p.nonfinite_optimal_gain,p.infeasible,p.iteration_limit,p.time_limit,p.unbounded)?;
        }
        write!(out, "]}}")
    }
}

// 仅测试：在真实add_row/resolve之后注入一次返回状态，驱动生产回滚控制流。
// bit1=小批失败，bit2=内部重建失败（进一步覆盖B&B外层弃全部根割）。
#[cfg(test)]
thread_local! {static CUT_FAILURES:std::cell::Cell<u8>=const{std::cell::Cell::new(0)};}
#[cfg(test)]
struct CutFailureGuard(u8);
#[cfg(test)]
impl CutFailureGuard {
    fn new(mask: u8) -> Self {
        Self(CUT_FAILURES.with(|v| v.replace(mask)))
    }
}
#[cfg(test)]
impl Drop for CutFailureGuard {
    fn drop(&mut self) {
        CUT_FAILURES.with(|v| v.set(self.0));
    }
}
#[cfg(test)]
pub(crate) fn test_cut_result(mut result: LpSolution, bit: u8) -> LpSolution {
    if CUT_FAILURES.with(|v| {
        let old = v.get();
        v.set(old & !bit);
        old & bit != 0
    }) {
        result.status = LpStatus::IterationLimit;
        result.objective = f64::NAN;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::branch_bound::{
        branch_and_bound, BbOptions, BinaryPropagationMode, BranchRule, SearchControl,
    };
    use crate::model::{Problem, RelOp, Sense, SolveStatus, Tolerances, VarType};
    use crate::simplex_sparse::LpEngineKind;
    use crate::solver::{solve_controlled, SolveOptions};

    fn field<'a>(event: &'a Event, key: &str) -> &'a Value {
        &event
            .fields
            .iter()
            .find(|(name, _)| *name == key)
            .unwrap()
            .1
    }
    fn number(event: &Event, key: &str) -> f64 {
        match field(event, key) {
            Value::Number(n) => *n,
            Value::Int(n) => *n as f64,
            other => panic!("非数字字段 {key}: {other:?}"),
        }
    }
    fn fractional_problem() -> Problem {
        let mut p = Problem::new("真实根阶段", Sense::Maximize);
        p.add_var("a", VarType::Binary).unwrap();
        p.add_var("b", VarType::Binary).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 7.0)], 0.25).unwrap();
        p.add_constraint(None, vec![(0, 5.0), (1, 3.0)], RelOp::Le, 7.0)
            .unwrap();
        p
    }
    fn plain(engine: LpEngineKind) -> BbOptions {
        BbOptions {
            root_cuts: false,
            heuristics: false,
            root_rc_fixing: false,
            local_cuts: false,
            binary_row_propagation: BinaryPropagationMode::Off,
            branch_rule_auto: false,
            branch_rule: BranchRule::MostFractional,
            node_limit: 1,
            lp_engine: engine,
            ..BbOptions::default()
        }
    }
    #[test]
    fn stages_come_from_actual_scaled_work_and_capture_does_not_change_result() {
        for engine in [LpEngineKind::Tableau, LpEngineKind::Revised] {
            for scaling in [false, true] {
                for presolve in [false, true] {
                    let mut p = fractional_problem();
                    p.add_var("fixed", VarType::Binary).unwrap();
                    p.set_var_bounds(2, 0.0, 0.0).unwrap();
                    let opts = SolveOptions {
                        bb: plain(engine),
                        use_scaling: scaling,
                        use_presolve: presolve,
                        ..SolveOptions::default()
                    };
                    reset(1);
                    let observed =
                        solve_controlled(&p, &opts, &SearchControl::default(), None).unwrap();
                    let trace = take();
                    reset(2);
                    let baseline =
                        solve_controlled(&p, &opts, &SearchControl::default(), None).unwrap();
                    let unsupported = take();
                    assert!(!unsupported.available && unsupported.events.is_empty());
                    assert_eq!(observed.report.status, baseline.report.status);
                    assert_eq!(observed.report.objective, baseline.report.objective);
                    assert_eq!(observed.report.values, baseline.report.values);
                    assert_eq!(observed.bound, baseline.bound);
                    assert_eq!(observed.report.stats.nodes, baseline.report.stats.nodes);
                    assert_eq!(
                        observed.report.stats.lp_iterations,
                        baseline.report.stats.lp_iterations
                    );
                    let context = trace.events.iter().find(|e| e.kind == "context").unwrap();
                    assert_eq!(
                        number(context, "work_variables"),
                        if presolve { 2.0 } else { 3.0 }
                    );
                    assert!(
                        matches!(field(context,"effective_engine"),Value::Text(value) if *value==if engine==LpEngineKind::Revised{"revised"}else{"tableau"})
                    );
                    let first = trace.events.iter().find(|e| e.kind == "root_lp").unwrap();
                    assert_eq!(number(first, "bound"), 15.25);
                    let after = trace
                        .events
                        .iter()
                        .find(|e| e.kind == "root_after_cuts")
                        .unwrap();
                    assert_eq!(number(after, "bound"), 15.25);
                    assert!(trace.events.iter().any(|e| e.kind == "root_ready"));
                    assert!(trace
                        .events
                        .iter()
                        .any(|e| e.kind == "frontier"
                            && matches!(field(e, "exit"), Value::Bool(true))));
                }
            }
        }
    }
    #[test]
    fn actual_zero_gain_probe_is_counted_without_learning_it() {
        let mut p = Problem::new("零增益真实探测", Sense::Minimize);
        p.add_var("x", VarType::Binary).unwrap();
        p.add_var("y", VarType::Binary).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Ge, 0.5)
            .unwrap();
        for engine in [LpEngineKind::Tableau, LpEngineKind::Revised] {
            let opts = BbOptions {
                branch_rule: BranchRule::Pseudocost,
                pc_reliability_probes: 2,
                ..plain(engine)
            };
            reset(1);
            let r = branch_and_bound(&p, &opts, Tolerances::default()).unwrap();
            let trace = take();
            assert_eq!(r.status, SolveStatus::LimitReached);
            assert!(trace.probes[0].calls >= 2);
            assert_eq!(trace.probes[0].zero, trace.probes[0].calls);
            let chosen = trace.events.iter().find(|e| e.kind == "chosen").unwrap();
            assert_eq!(
                number(chosen, "observations"),
                0.0,
                "诊断不能把零gain送进学习"
            );
            assert!(matches!(
                field(chosen, "probed_this_node"),
                Value::Bool(true)
            ));
        }
    }
    #[test]
    fn root_cut_rollback_observes_actual_retained_rows() {
        let p = fractional_problem();
        for kind in [LpEngineKind::Tableau, LpEngineKind::Revised] {
            for mask in [1, 3] {
                let _failure = CutFailureGuard::new(mask);
                reset(1);
                let r = branch_and_bound(
                    &p,
                    &BbOptions {
                        root_cuts: true,
                        ..plain(kind)
                    },
                    Tolerances::default(),
                )
                .unwrap();
                let trace = take();
                let rollback = trace
                    .events
                    .iter()
                    .find(|e| e.kind == "cut_rollback")
                    .expect("必须真实执行小批丢弃与重建控制流");
                assert!(number(rollback, "discarded") > 0.0);
                let ready = trace
                    .events
                    .iter()
                    .find(|e| e.kind == "root_after_cuts")
                    .unwrap();
                assert_eq!(number(ready, "retained_cuts"), r.stats.cuts as f64);
                assert_eq!(
                    number(ready, "bound"),
                    15.25,
                    "不能使用失败小批返回的NaN或尝试中界"
                );
                assert_eq!(
                    trace.events.iter().any(|e| e.kind == "cut_outer_rollback"),
                    mask == 3
                );
                if mask == 3 {
                    assert!(matches!(field(rollback, "bound"), Value::Null));
                } else {
                    assert!((number(ready, "bound") - number(rollback, "bound")).abs() < 1e-7);
                }
                assert!(trace.events.iter().any(|e| e.kind == "cut_round"));
            }
        }
    }
    #[test]
    fn nonoptimal_objectives_and_nonfinite_numbers_never_become_json_bounds() {
        reset(1);
        for state in [
            LpStatus::IterationLimit,
            LpStatus::TimeLimit,
            LpStatus::Infeasible,
            LpStatus::Unbounded,
        ] {
            lp(
                "root_lp",
                &LpSolution {
                    status: state,
                    objective: 123.0,
                    x: Vec::new(),
                    duals: Vec::new(),
                },
                0,
                Some(f64::NAN),
                0,
            );
        }
        record(
            "special",
            vec![
                ("infinity", f64::INFINITY.into()),
                ("negative_infinity", f64::NEG_INFINITY.into()),
            ],
        );
        let trace = take();
        let mut bytes = Vec::new();
        trace.write_json(&mut bytes).unwrap();
        let json = String::from_utf8(bytes).unwrap();
        assert!(!json.contains(":NaN") && !json.contains(":inf") && !json.contains(":-inf"));
        for event in trace.events.iter().filter(|e| e.kind == "root_lp") {
            assert!(matches!(field(event, "bound"), Value::Null));
        }
        assert!(json.contains("\"seconds\":null"));
    }
    #[test]
    fn diagnostic_error_is_fresh_per_call_for_both_engines() {
        let p = fractional_problem();
        for kind in [LpEngineKind::Tableau, LpEngineKind::Revised] {
            let mut engine =
                crate::simplex_sparse::new_engine(&p, Tolerances::default(), kind).unwrap();
            let valid = crate::simplex::Bounds::from_problem(&p);
            let mut invalid = valid.clone();
            invalid.lb[0] = f64::NAN;
            let failed = engine.solve(&invalid);
            assert_ne!(failed.status, LpStatus::Optimal);
            assert!(engine.proof_diagnostics().numerical_error.is_some());
            assert_eq!(engine.solve(&valid).status, LpStatus::Optimal);
            assert!(engine.proof_diagnostics().numerical_error.is_none());
            let _ = engine.resolve(&invalid);
            assert!(engine.proof_diagnostics().numerical_error.is_some());
            assert_eq!(engine.resolve(&valid).status, LpStatus::Optimal);
            assert!(engine.proof_diagnostics().numerical_error.is_none());
        }
        let mut bytes = Vec::new();
        write_value(&mut bytes, &Value::String("错误\"\\\n".into())).unwrap();
        assert_eq!(String::from_utf8(bytes).unwrap(), "\"错误\\\"\\\\\\u000a\"");
    }
    #[test]
    fn interrupted_child_is_reported_separately_from_empty_execution_stack() {
        use crate::branch_bound::{test_node_lp_interrupt::Guard, unresolved_bound_tests};
        for sense in [Sense::Minimize, Sense::Maximize] {
            let mut p = Problem::new("诊断未完成子树", sense);
            p.add_var("b", VarType::Binary).unwrap();
            p.set_objective(
                vec![(0, if sense == Sense::Minimize { -1.0 } else { 1.0 })],
                0.25,
            )
            .unwrap();
            p.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 0.5)
                .unwrap();
            for engine in [LpEngineKind::Tableau, LpEngineKind::Revised] {
                for state in [LpStatus::IterationLimit, LpStatus::TimeLimit] {
                    let opts = BbOptions {
                        branch_rule: BranchRule::Pseudocost,
                        reuse_probe_evidence: true,
                        ..unresolved_bound_tests::options(engine)
                    };
                    reset(1);
                    let result = {
                        let _guard =
                            Guard::install(Some(unresolved_bound_tests::interrupt(state, 0, 1)));
                        branch_and_bound(&p, &opts, Tolerances::default()).unwrap()
                    };
                    let trace = take();
                    assert_eq!(result.status, SolveStatus::LimitReached);
                    assert_eq!(trace.unresolved_lp_nodes, 1);
                    assert_eq!(trace.unresolved_est, result.bound);
                    let frontier = trace
                        .events
                        .iter()
                        .rev()
                        .find(|e| e.kind == "frontier")
                        .unwrap();
                    assert_eq!(number(frontier, "size"), 0.0);
                    assert!(matches!(field(frontier, "stack_bound"), Value::Null));
                    assert_eq!(number(frontier, "unresolved_count"), 1.0);
                    assert_eq!(Some(number(frontier, "unresolved_bound")), result.bound);
                    assert_eq!(Some(number(frontier, "combined_bound")), result.bound);
                    let node = trace
                        .events
                        .iter()
                        .find(|e| {
                            e.kind == "node" && e.fields.iter().any(|(k, _)| *k == "actual_bound")
                        })
                        .unwrap();
                    assert!(matches!(field(node, "actual_bound"), Value::Null));
                }
            }
        }
    }
    #[test]
    fn noncheckpoint_frontier_does_not_evaluate_or_scan() {
        reset(1);
        frontier(
            Sense::Minimize,
            2,
            false,
            None,
            None,
            || panic!("非采样节点不能聚合生产前沿"),
            std::iter::from_fn(|| -> Option<(f64, u32)> {
                panic!("非采样节点不能扫描栈")
            }),
        );
        frontier(
            Sense::Minimize,
            64,
            false,
            None,
            None,
            || Some(3.0),
            std::iter::once((3.0, 1)),
        );
        frontier(
            Sense::Minimize,
            64,
            false,
            None,
            None,
            || panic!("重复节点数不能重复采样"),
            std::iter::empty(),
        );
        let trace = take();
        assert_eq!(
            trace.events.iter().filter(|e| e.kind == "frontier").count(),
            1
        );
    }
}
