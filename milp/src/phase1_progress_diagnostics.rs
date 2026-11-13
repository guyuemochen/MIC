//! 显式启用、线程内且有界的第一阶段只读证据；所有浮点量按IEEE754位串输出。
//! reset必须在待观察引擎构造前调用；只认领随后第一个实际第一阶段调用。
use crate::model::{Problem, RelOp, Tolerances};
use crate::simplex::Bounds;
use std::cell::RefCell;
use std::io::{self, Write};
use std::sync::atomic::{AtomicU64, Ordering};

pub const BYTE_LIMIT: usize = 256 * 1024 * 1024;
pub const EVENT_LIMIT: usize = 4096;
const SOURCE_LIMIT: usize = 8;
const SNAPSHOT_LIMIT: usize = 4;
// 预留事件和快照的固定容量，以及线程状态、来源表和两条短错误消息。
const FIXED_BYTES: usize = 16_384
    + EVENT_LIMIT * std::mem::size_of::<Summary>()
    + SNAPSHOT_LIMIT * std::mem::size_of::<Snapshot>();
static NEXT_ID: AtomicU64 = AtomicU64::new(1);
fn id() -> u64 {
    NEXT_ID
        .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |v| v.checked_add(1))
        .unwrap_or(0)
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Pricing {
    pub kind: &'static str,
    pub column: Option<usize>,
    pub direction: Option<f64>,
    pub evaluated_revision: Option<u64>,
}
#[derive(Clone, Debug)]
pub struct Summary {
    pub kind: &'static str,
    pub state_revision: u64,
    pub iterations: u64,
    pub call_iterations: u64,
    pub eta_count: usize,
    pub bland: bool,
    pub artificial_basic: usize,
    pub artificial_nonzero: usize,
    pub nonfinite: usize,
    pub positive_sum: f64,
    pub negative_sum: f64,
    pub minimum: Option<f64>,
    pub maximum: Option<f64>,
    pub basic_bound_violation: f64,
    pub pricing: Pricing,
}
#[derive(Clone, Debug)]
pub struct Snapshot {
    pub kind: &'static str,
    pub summary: Summary,
    pub x_unclipped: Vec<f64>,
    pub basis: Vec<usize>,
    pub basic_values: Vec<f64>,
    pub at_upper: Vec<bool>,
    pub lower: Vec<f64>,
    pub upper: Vec<f64>,
    pub factorization_broken: bool,
    /// 重构失败的局部新值独立命名，不把旧基本值冒充重算结果。
    pub failed_recomputed_values: Option<Vec<f64>>,
    pub failed_refresh_value: Option<(usize, f64)>,
}
#[derive(Clone, Debug)]
pub struct Matrix {
    pub structural_columns: usize,
    pub static_columns: usize,
    pub row_terms: Vec<Vec<(usize, f64)>>,
    pub row_rhs: Vec<f64>,
    pub row_op: Vec<RelOp>,
    pub row_scale: Vec<f64>,
    pub row_slack: Vec<Option<usize>>,
    pub constraint_row: Vec<Option<usize>>,
    pub cut_row: Vec<Option<usize>>,
    pub artificial_columns: Option<Vec<(usize, f64)>>,
}
#[derive(Debug)]
pub struct Capture {
    pub tolerances: Tolerances,
    pub slot_attempt_iterations: [Option<u64>; 4],
    pub slot_omitted: [bool; 4],
    pub engine: &'static str,
    pub source_id: u64,
    pub engine_input: Option<Problem>,
    pub matrix: Option<Matrix>,
    pub events: Vec<Summary>,
    pub snapshots: Vec<Snapshot>,
    pub events_omitted: u64,
    pub objects_omitted: u64,
    pub omitted_estimated_bytes: usize,
    pub source_unavailable: bool,
    pub artificial_mapping_unavailable: bool,
    pub estimated_bytes: usize,
    pub first_failure_detail: Option<String>,
    pub terminal_detail: Option<String>,
    pub terminal: &'static str,
}
struct Source {
    id: u64,
    problem: Problem,
    bytes: usize,
}
#[derive(Default)]
struct State {
    armed: bool,
    claimed: bool,
    owner: usize,
    sources: Vec<Source>,
    bytes: usize,
    omitted: u64,
    omitted_bytes: usize,
    capture: Option<Capture>,
    pricing: Pricing,
    revision: u64,
    failed_refresh_value: Option<(usize, f64)>,
    slots: [bool; 4],
}
thread_local! { static STATE: RefCell<State> = RefCell::new(State::default()); }

/// 仅在求解调用之外启用；来源编号跨reset不复用，引擎克隆只复制编号。
pub fn reset() {
    STATE.with(|s| {
        *s.borrow_mut() = State {
            armed: true,
            bytes: FIXED_BYTES,
            ..State::default()
        }
    });
}
/// 取走即关闭；未发生实际第一阶段时返回None。
pub fn take() -> Option<Capture> {
    STATE.with(|s| {
        let mut s = s.borrow_mut();
        s.armed = false;
        s.owner = 0;
        s.sources.clear();
        s.capture.take()
    })
}
fn sized(mut parts: impl Iterator<Item = usize>) -> usize {
    parts
        .try_fold(0usize, |a, b| a.checked_add(b))
        .unwrap_or(usize::MAX)
}
fn reserve(s: &mut State, bytes: usize) -> bool {
    if s.bytes
        .checked_add(bytes)
        .map_or(false, |n| n <= BYTE_LIMIT)
    {
        s.bytes += bytes;
        true
    } else {
        s.omitted = s.omitted.saturating_add(1);
        s.omitted_bytes = s.omitted_bytes.saturating_add(bytes);
        false
    }
}
fn problem_size(p: &Problem) -> usize {
    sized(
        std::iter::once(std::mem::size_of::<Problem>().saturating_add(p.name.len()))
            .chain(
                p.vars
                    .iter()
                    .map(|v| std::mem::size_of_val(v).saturating_add(v.name.len())),
            )
            .chain(std::iter::once(p.obj.terms.len().saturating_mul(16)))
            .chain(p.constraints.iter().map(|r| {
                std::mem::size_of_val(r)
                    .saturating_add(r.name.as_ref().map_or(0, String::len))
                    .saturating_add(r.expr.terms.len().saturating_mul(16))
            })),
    )
}
pub(crate) fn register_source(p: &Problem) -> u64 {
    STATE.with(|cell| {
        let mut s = cell.borrow_mut();
        if !s.armed || s.claimed {
            return 0;
        }
        let token = id();
        let bytes = problem_size(p);
        if token == 0 {
            s.omitted = s.omitted.saturating_add(1);
            s.omitted_bytes = s.omitted_bytes.saturating_add(bytes);
            return 0;
        }
        if s.sources.len() >= SOURCE_LIMIT {
            s.omitted = s.omitted.saturating_add(1);
            s.omitted_bytes = s.omitted_bytes.saturating_add(bytes);
        } else if reserve(&mut s, bytes) {
            s.sources.push(Source {
                id: token,
                problem: p.clone(),
                bytes,
            });
        }
        token
    })
}

pub(crate) struct MatrixView<'a> {
    pub tolerances: Tolerances,
    pub n: usize,
    pub ncs: usize,
    pub terms: &'a [Vec<(usize, f64)>],
    pub rhs: &'a [f64],
    pub op: &'a [RelOp],
    pub scale: &'a [f64],
    pub slack: &'a [Option<usize>],
    pub constraint_row: &'a [Option<usize>],
    pub cut_row: &'a [Option<usize>],
    pub art: Option<&'a [(usize, f64)]>,
    pub dense_initial: Option<(&'a [usize], &'a [Vec<f64>])>,
    pub art_count: usize,
}
pub(crate) struct View<'a> {
    pub n: usize,
    pub ncs: usize,
    pub basis: &'a [usize],
    pub xb: &'a [f64],
    pub in_basis: &'a [bool],
    pub at_upper: &'a [bool],
    pub bounds: &'a Bounds,
    pub iterations: u64,
    pub call_iterations: u64,
    pub eta: usize,
    pub bland: bool,
    pub feasibility: f64,
    pub broken: bool,
}
pub(crate) struct Scope {
    owner: usize,
}
impl Drop for Scope {
    fn drop(&mut self) {
        finish(self.owner, "unknown_or_unwound");
    }
}
pub(crate) fn active(owner: usize) -> bool {
    STATE.with(|s| {
        let s = s.borrow();
        s.armed && s.claimed && s.owner == owner && owner != 0
    })
}
pub(crate) fn begin(owner: usize, source: u64, engine: &'static str, m: MatrixView<'_>) -> Scope {
    STATE.with(|cell| {
        let mut s = cell.borrow_mut();
        if !s.armed || s.claimed {
            return;
        }
        s.claimed = true;
        s.owner = owner;
        let input = s
            .sources
            .iter()
            .position(|v| source != 0 && v.id == source)
            .map(|i| s.sources.swap_remove(i));
        s.sources.clear();
        s.bytes = FIXED_BYTES.saturating_add(input.as_ref().map_or(0, |v| v.bytes));
        // 稠密生产第一阶段入口刚建±单位基；不对演进后的逆矩阵做反推。
        let dense_valid = m.dense_initial.map_or(false, |(basis, inverse)| {
            basis.len() == m.rhs.len()
                && inverse.len() == basis.len()
                && inverse.iter().enumerate().all(|(i, row)| {
                    row.len() == basis.len()
                        && row.iter().enumerate().all(|(j, &v)| {
                            if i == j {
                                v == 1.0 || v == -1.0
                            } else {
                                v == 0.0
                            }
                        })
                })
                && (0..m.art_count).all(|k| basis.iter().filter(|&&j| j == m.ncs + k).count() == 1)
        });
        let art_available = m.art.is_some() || dense_valid;
        let bytes = sized(
            std::iter::once(std::mem::size_of::<Matrix>())
                .chain(m.terms.iter().map(|r| {
                    std::mem::size_of::<Vec<(usize, f64)>>()
                        .saturating_add(r.len().saturating_mul(16))
                }))
                .chain([
                    m.rhs.len().saturating_mul(8),
                    m.op.len().saturating_mul(std::mem::size_of::<RelOp>()),
                    m.scale.len().saturating_mul(8),
                    m.slack
                        .len()
                        .saturating_add(m.constraint_row.len())
                        .saturating_add(m.cut_row.len())
                        .saturating_mul(std::mem::size_of::<Option<usize>>()),
                    m.art_count.saturating_mul(16),
                ]),
        );
        let matrix = if reserve(&mut s, bytes) {
            Some(Matrix {
                structural_columns: m.n,
                static_columns: m.ncs,
                row_terms: m.terms.to_vec(),
                row_rhs: m.rhs.to_vec(),
                row_op: m.op.to_vec(),
                row_scale: m.scale.to_vec(),
                row_slack: m.slack.to_vec(),
                constraint_row: m.constraint_row.to_vec(),
                cut_row: m.cut_row.to_vec(),
                artificial_columns: if let Some(a) = m.art {
                    Some(a.to_vec())
                } else if dense_valid {
                    let (basis, inverse) = m.dense_initial.unwrap();
                    let mut a = vec![(0, 0.0); m.art_count];
                    for (i, &j) in basis.iter().enumerate() {
                        if j >= m.ncs && j < m.ncs + m.art_count {
                            a[j - m.ncs] = (i, inverse[i][i]);
                        }
                    }
                    Some(a)
                } else {
                    None
                },
            })
        } else {
            None
        };
        s.capture = Some(Capture {
            tolerances: m.tolerances,
            slot_attempt_iterations: [None; 4],
            slot_omitted: [false; 4],
            engine,
            source_id: source,
            source_unavailable: input.is_none(),
            engine_input: input.map(|v| v.problem),
            matrix,
            events: Vec::with_capacity(EVENT_LIMIT),
            snapshots: Vec::with_capacity(SNAPSHOT_LIMIT),
            events_omitted: 0,
            objects_omitted: s.omitted,
            omitted_estimated_bytes: s.omitted_bytes,
            artificial_mapping_unavailable: !art_available,
            estimated_bytes: s.bytes,
            first_failure_detail: None,
            terminal_detail: None,
            terminal: "running",
        });
    });
    Scope { owner }
}
/// 状态变化前只递增诊断世代，旧定价明确不适用于新基/重算值。
pub(crate) fn invalidate(owner: usize) {
    STATE.with(|cell| {
        let mut s = cell.borrow_mut();
        if s.owner == owner && s.armed {
            s.revision = s.revision.saturating_add(1);
        }
    });
}
pub(crate) fn failed_refresh(owner: usize, row: usize, value: f64) {
    STATE.with(|cell| {
        let mut s = cell.borrow_mut();
        if s.owner == owner && s.armed {
            s.failed_refresh_value = Some((row, value));
        }
    });
}
pub(crate) fn pricing(owner: usize, mut value: Pricing, bland: bool) {
    STATE.with(|cell| {
        let mut s = cell.borrow_mut();
        if s.owner != owner || !s.armed {
            return;
        }
        value.evaluated_revision = Some(s.revision);
        s.pricing = value;
        if let Some(c) = s.capture.as_mut() {
            if let Some(candidate) = c.snapshots.iter().find(|v| v.kind == "first_candidate") {
                if candidate.summary.state_revision == value.evaluated_revision.unwrap() {
                    let mut event = candidate.summary.clone();
                    event.kind = "pricing_at_first_candidate";
                    event.pricing = value;
                    event.bland = bland;
                    if c.events.len() < EVENT_LIMIT {
                        c.events.push(event);
                    } else {
                        c.events_omitted = c.events_omitted.saturating_add(1);
                    }
                }
            }
        }
    });
}
/// 第一次故障原因与最终退出原因分开；重试恢复不能改写或冒充最终错误。
pub(crate) fn failure_reason(owner: usize, value: &str) {
    STATE.with(|cell| {
        let mut s = cell.borrow_mut();
        if s.owner == owner && s.armed {
            if let Some(c) = s.capture.as_mut() {
                if c.first_failure_detail.is_none() {
                    c.first_failure_detail = Some(value.chars().take(256).collect());
                }
            }
        }
    });
}
pub(crate) fn detail(owner: usize, value: Option<&str>) {
    if let Some(reason) = value {
        failure_reason(owner, reason);
    }
    STATE.with(|cell| {
        let mut s = cell.borrow_mut();
        if s.owner == owner && s.armed {
            if let Some(c) = s.capture.as_mut() {
                c.terminal_detail = value.map(|v| v.chars().take(256).collect());
            }
        }
    });
}
pub(crate) fn finish(owner: usize, terminal: &'static str) {
    STATE.with(|cell| {
        let mut s = cell.borrow_mut();
        if s.owner == owner && s.armed && s.claimed {
            if let Some(c) = s.capture.as_mut() {
                c.terminal = terminal;
            }
            s.armed = false;
            s.owner = 0;
        }
    });
}

/// 额外扫描只在显式诊断时执行；不调用定价、重构或求解，不修改引擎。
pub(crate) fn observe(owner: usize, kind: &'static str, v: View<'_>, failed: Option<&[f64]>) {
    if !active(owner) {
        return;
    }
    STATE.with(|cell| {
        let mut s = cell.borrow_mut();
        let mut summary = Summary {
            kind,
            state_revision: s.revision,
            iterations: v.iterations,
            call_iterations: v.call_iterations,
            eta_count: v.eta,
            bland: v.bland,
            artificial_basic: 0,
            artificial_nonzero: 0,
            nonfinite: 0,
            positive_sum: 0.0,
            negative_sum: 0.0,
            minimum: None,
            maximum: None,
            basic_bound_violation: 0.0,
            pricing: s.pricing,
        };
        for (i, &j) in v.basis.iter().enumerate() {
            let value = v.xb[i];
            if !value.is_finite() {
                summary.nonfinite += 1;
            }
            let (lb, ub) = if j < v.n {
                (v.bounds.lb[j], v.bounds.ub[j])
            } else {
                (0.0, f64::INFINITY)
            };
            if value.is_finite() {
                summary.basic_bound_violation = summary
                    .basic_bound_violation
                    .max((lb - value).max(value - ub).max(0.0));
            } else {
                summary.basic_bound_violation = f64::NAN;
            }
            if j >= v.ncs {
                summary.artificial_basic += 1;
                if value != 0.0 {
                    summary.artificial_nonzero += 1;
                }
                if value.is_nan() {
                    summary.positive_sum = f64::NAN;
                    summary.negative_sum = f64::NAN;
                } else if value >= 0.0 {
                    summary.positive_sum += value;
                } else {
                    summary.negative_sum += value;
                }
                if value.is_nan() {
                    summary.minimum = Some(value);
                    summary.maximum = Some(value);
                } else {
                    summary.minimum = Some(summary.minimum.map_or(value, |x| {
                        if x.is_nan() {
                            x
                        } else {
                            x.min(value)
                        }
                    }));
                    summary.maximum = Some(summary.maximum.map_or(value, |x| {
                        if x.is_nan() {
                            x
                        } else {
                            x.max(value)
                        }
                    }));
                }
            }
        }
        if summary.nonfinite > 0 {
            summary.basic_bound_violation = f64::NAN;
        }
        let candidate = summary.nonfinite == 0
            && summary.positive_sum.is_finite()
            && summary.positive_sum <= v.feasibility
            && summary.minimum.map_or(true, |x| x >= -v.feasibility);
        let slot = match kind {
            "entry" => Some(0),
            "failure" | "failed_recompute" | "partial_refresh" => Some(2),
            "before_cleanup" => Some(3),
            _ => {
                if candidate && !s.slots[1] {
                    Some(1)
                } else {
                    None
                }
            }
        };
        let emit = kind != "step" || v.call_iterations % 1024 == 0 || slot == Some(1);
        if emit {
            if let Some(c) = s.capture.as_mut() {
                if c.events.len() < EVENT_LIMIT {
                    c.events.push(summary.clone());
                } else {
                    c.events_omitted = c.events_omitted.saturating_add(1);
                }
            }
        }
        if let Some(slot) = slot {
            if !s.slots[slot] {
                s.slots[slot] = true;
                s.capture.as_mut().unwrap().slot_attempt_iterations[slot] = Some(v.call_iterations);
                let bytes =
                    v.n.saturating_mul(24)
                        .saturating_add(v.basis.len().saturating_mul(16))
                        .saturating_add(v.at_upper.len())
                        .saturating_add(failed.map_or(0, |x| x.len().saturating_mul(8)));
                if reserve(&mut s, bytes) {
                    let mut x = vec![0.0; v.n];
                    for j in 0..v.n {
                        if !v.in_basis[j] {
                            x[j] = if v.at_upper[j] {
                                v.bounds.ub[j]
                            } else {
                                v.bounds.lb[j]
                            };
                        }
                    }
                    for (i, &j) in v.basis.iter().enumerate() {
                        if j < v.n {
                            x[j] = v.xb[i];
                        }
                    }
                    let snapshot = Snapshot {
                        kind: if slot == 1 { "first_candidate" } else { kind },
                        summary,
                        x_unclipped: x,
                        basis: v.basis.to_vec(),
                        basic_values: v.xb.to_vec(),
                        at_upper: v.at_upper.to_vec(),
                        lower: v.bounds.lb.clone(),
                        upper: v.bounds.ub.clone(),
                        factorization_broken: v.broken,
                        failed_recomputed_values: failed.map(|x| x.to_vec()),
                        failed_refresh_value: s.failed_refresh_value,
                    };
                    s.capture.as_mut().unwrap().snapshots.push(snapshot);
                } else {
                    s.capture.as_mut().unwrap().slot_omitted[slot] = true;
                }
            }
        }
        let (bytes, omitted, omitted_bytes) = (s.bytes, s.omitted, s.omitted_bytes);
        if let Some(c) = s.capture.as_mut() {
            c.estimated_bytes = bytes;
            c.objects_omitted = omitted;
            c.omitted_estimated_bytes = omitted_bytes;
        }
    });
}

fn text(w: &mut impl Write, s: &str) -> io::Result<()> {
    write!(w, "\"")?;
    for c in s.chars() {
        match c {
            '"' => write!(w, "\\\"")?,
            '\\' => write!(w, "\\\\")?,
            '\u{0}'..='\u{1f}' => write!(w, "\\u{:04x}", c as u32)?,
            _ => write!(w, "{c}")?,
        }
    }
    write!(w, "\"")
}
fn float(w: &mut impl Write, x: f64) -> io::Result<()> {
    write!(w, "\"{:016x}\"", x.to_bits())
}
fn floats(w: &mut impl Write, x: &[f64]) -> io::Result<()> {
    write!(w, "[")?;
    for (i, &v) in x.iter().enumerate() {
        if i > 0 {
            write!(w, ",")?;
        }
        float(w, v)?;
    }
    write!(w, "]")
}
fn indices(w: &mut impl Write, x: &[usize]) -> io::Result<()> {
    write!(w, "{x:?}")
}
fn optional_indices(w: &mut impl Write, x: &[Option<usize>]) -> io::Result<()> {
    write!(w, "[")?;
    for (i, v) in x.iter().enumerate() {
        if i > 0 {
            write!(w, ",")?;
        }
        match v {
            Some(j) => write!(w, "{j}")?,
            None => write!(w, "null")?,
        }
    }
    write!(w, "]")
}
fn terms(w: &mut impl Write, x: &[(usize, f64)]) -> io::Result<()> {
    write!(w, "[")?;
    for (i, &(j, a)) in x.iter().enumerate() {
        if i > 0 {
            write!(w, ",")?;
        }
        write!(w, "[{j},")?;
        float(w, a)?;
        write!(w, "]")?;
    }
    write!(w, "]")
}
impl Summary {
    fn write_json(&self, w: &mut impl Write) -> io::Result<()> {
        write!(w, "{{\"kind\":")?;
        text(w, self.kind)?;
        write!(w,",\"state_revision\":{},\"phase\":1,\"iterations\":{},\"call_iterations\":{},\"eta_count\":{},\"bland\":{},\"artificial_basic\":{},\"artificial_nonzero\":{},\"nonfinite\":{},\"positive_sum\":",self.state_revision,self.iterations,self.call_iterations,self.eta_count,self.bland,self.artificial_basic,self.artificial_nonzero,self.nonfinite)?;
        float(w, self.positive_sum)?;
        write!(w, ",\"negative_sum\":")?;
        float(w, self.negative_sum)?;
        write!(w, ",\"minimum\":")?;
        match self.minimum {
            Some(v) => float(w, v)?,
            None => write!(w, "null")?,
        }
        write!(w, ",\"maximum\":")?;
        match self.maximum {
            Some(v) => float(w, v)?,
            None => write!(w, "null")?,
        }
        write!(w, ",\"basic_bound_violation\":")?;
        float(w, self.basic_bound_violation)?;
        write!(w, ",\"pricing\":{{\"kind\":")?;
        text(w, self.pricing.kind)?;
        write!(w, ",\"column\":")?;
        match self.pricing.column {
            Some(j) => write!(w, "{j}")?,
            None => write!(w, "null")?,
        }
        write!(w, ",\"evaluated_revision\":")?;
        match self.pricing.evaluated_revision {
            Some(v) => write!(w, "{v}")?,
            None => write!(w, "null")?,
        }
        write!(w, ",\"direction\":")?;
        match self.pricing.direction {
            Some(v) => float(w, v)?,
            None => write!(w, "null")?,
        }
        write!(w, "}}}}")
    }
}
/// 引擎构造LP的精确输入；与外层文件/预处理还原坐标严格区分。
pub fn write_problem_json(p: &Problem, w: &mut impl Write) -> io::Result<()> {
    write!(w, "{{\"name\":")?;
    text(w, &p.name)?;
    write!(w, ",\"sense\":\"{:?}\",\"offset\":", p.sense)?;
    float(w, p.obj_offset)?;
    write!(w, ",\"objective\":")?;
    terms(w, &p.obj.terms)?;
    write!(w, ",\"variables\":[")?;
    for (i, v) in p.vars.iter().enumerate() {
        if i > 0 {
            write!(w, ",")?;
        }
        write!(w, "{{\"name\":")?;
        text(w, &v.name)?;
        write!(w, ",\"type\":\"{:?}\",\"lower\":", v.var_type)?;
        float(w, v.lb)?;
        write!(w, ",\"upper\":")?;
        float(w, v.ub)?;
        write!(w, "}}")?;
    }
    write!(w, "],\"rows\":[")?;
    for (i, r) in p.constraints.iter().enumerate() {
        if i > 0 {
            write!(w, ",")?;
        }
        write!(w, "{{\"terms\":")?;
        terms(w, &r.expr.terms)?;
        write!(w, ",\"op\":\"{:?}\",\"rhs\":", r.op)?;
        float(w, r.rhs)?;
        write!(w, "}}")?;
    }
    write!(w, "]}}")
}
impl Capture {
    pub fn write_json(&self, w: &mut impl Write) -> io::Result<()> {
        write!(
            w,
            "{{\"schema\":1,\"float_encoding\":\"ieee754_binary64_hex_bits\",\"engine\":"
        )?;
        text(w, self.engine)?;
        write!(w,",\"source_id\":{},\"estimated_bytes\":{},\"byte_limit\":{},\"events_omitted\":{},\"objects_omitted\":{},\"omitted_estimated_bytes\":{},\"source_unavailable\":{},\"artificial_mapping_unavailable\":{},\"terminal\":",self.source_id,self.estimated_bytes,BYTE_LIMIT,self.events_omitted,self.objects_omitted,self.omitted_estimated_bytes,self.source_unavailable,self.artificial_mapping_unavailable)?;
        text(w, self.terminal)?;
        write!(w, ",\"first_failure_detail\":")?;
        match &self.first_failure_detail {
            Some(v) => text(w, v)?,
            None => write!(w, "null")?,
        }
        write!(w, ",\"terminal_detail\":")?;
        match &self.terminal_detail {
            Some(v) => text(w, v)?,
            None => write!(w, "null")?,
        }
        write!(w, ",\"tolerances\":{{\"feasibility\":")?;
        float(w, self.tolerances.feasibility)?;
        write!(w, ",\"optimality\":")?;
        float(w, self.tolerances.optimality)?;
        write!(w, ",\"pivot\":")?;
        float(w, self.tolerances.pivot)?;
        write!(w, ",\"zero\":")?;
        float(w, self.tolerances.zero)?;
        write!(w, ",\"integrality\":")?;
        float(w, self.tolerances.integrality)?;
        write!(w, "}}")?;
        write!(w,",\"slot_order\":[\"entry\",\"first_candidate\",\"first_failure\",\"before_cleanup\"],\"slot_omitted\":{:?},\"slot_attempt_iterations\":[",self.slot_omitted)?;
        for (i, v) in self.slot_attempt_iterations.iter().enumerate() {
            if i > 0 {
                write!(w, ",")?;
            }
            match v {
                Some(v) => write!(w, "{v}")?,
                None => write!(w, "null")?,
            }
        }
        write!(w, "]")?;
        write!(w, ",\"engine_input_lp\":")?;
        match &self.engine_input {
            Some(p) => write_problem_json(p, w)?,
            None => write!(w, "null")?,
        }
        write!(w, ",\"normalized_matrix\":")?;
        if let Some(m) = &self.matrix {
            write!(
                w,
                "{{\"structural_columns\":{},\"static_columns\":{},\"row_terms\":[",
                m.structural_columns, m.static_columns
            )?;
            for (i, r) in m.row_terms.iter().enumerate() {
                if i > 0 {
                    write!(w, ",")?;
                }
                terms(w, r)?;
            }
            write!(w, "],\"row_rhs\":")?;
            floats(w, &m.row_rhs)?;
            write!(w, ",\"row_op\":[")?;
            for (i, op) in m.row_op.iter().enumerate() {
                if i > 0 {
                    write!(w, ",")?;
                }
                write!(w, "\"{op:?}\"")?;
            }
            write!(w, "],\"row_scale\":")?;
            floats(w, &m.row_scale)?;
            write!(w, ",\"row_slack\":")?;
            optional_indices(w, &m.row_slack)?;
            write!(w, ",\"constraint_row\":")?;
            optional_indices(w, &m.constraint_row)?;
            write!(w, ",\"cut_row\":")?;
            optional_indices(w, &m.cut_row)?;
            write!(w, ",\"artificial_columns\":")?;
            match &m.artificial_columns {
                Some(a) => terms(w, a)?,
                None => write!(w, "null")?,
            }
            write!(w, "}}")?;
        } else {
            write!(w, "null")?;
        }
        write!(w, ",\"events\":[")?;
        for (i, e) in self.events.iter().enumerate() {
            if i > 0 {
                write!(w, ",")?;
            }
            e.write_json(w)?;
        }
        write!(w, "],\"snapshots\":[")?;
        for (i, s) in self.snapshots.iter().enumerate() {
            if i > 0 {
                write!(w, ",")?;
            }
            write!(w, "{{\"kind\":")?;
            text(w, s.kind)?;
            write!(w, ",\"summary\":")?;
            s.summary.write_json(w)?;
            write!(w, ",\"x_unclipped\":")?;
            floats(w, &s.x_unclipped)?;
            write!(w, ",\"basis\":")?;
            indices(w, &s.basis)?;
            write!(w, ",\"basic_values\":")?;
            floats(w, &s.basic_values)?;
            write!(w, ",\"at_upper\":{:?},\"lower\":", s.at_upper)?;
            floats(w, &s.lower)?;
            write!(w, ",\"upper\":")?;
            floats(w, &s.upper)?;
            write!(
                w,
                ",\"factorization_broken\":{},\"failed_recomputed_values\":",
                s.factorization_broken
            )?;
            match &s.failed_recomputed_values {
                Some(v) => floats(w, v)?,
                None => write!(w, "null")?,
            }
            write!(w, ",\"failed_refresh_row_value\":")?;
            match s.failed_refresh_value {
                Some((row, value)) => {
                    write!(w, "[{row},")?;
                    float(w, value)?;
                    write!(w, "]")?;
                }
                None => write!(w, "null")?,
            }
            write!(w, "}}")?;
        }
        write!(w, "]}}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Sense, Tolerances, VarType};
    use crate::simplex::{LpStatus, Simplex};
    use crate::simplex_sparse::SparseSimplex;

    fn problem(negative: bool) -> Problem {
        let mut p = Problem::new("phase1_readonly", Sense::Minimize);
        let x = p.add_var("x", VarType::Continuous).unwrap();
        let y = p.add_var("y", VarType::Continuous).unwrap();
        p.add_constraint(
            None,
            vec![(x, 1.0)],
            RelOp::Eq,
            if negative {
                1000.0 + 2.0_f64.powi(-21)
            } else {
                1000.0
            },
        )
        .unwrap();
        p.add_constraint(None, vec![(x, 1.0), (y, 0.5)], RelOp::Eq, 1000.0)
            .unwrap();
        p
    }

    #[test]
    fn runtime_disabled_and_enabled_keep_public_fixed_limit_results() {
        for negative in [false, true] {
            let p = problem(negative);
            let bounds = Bounds::from_problem(&p);
            take();
            let mut dense = Simplex::new(&p, Tolerances::default()).unwrap();
            dense.set_iter_limit_hint(Some(1000));
            let expected = dense.solve(&bounds);
            reset();
            let mut observed = Simplex::new(&p, Tolerances::default()).unwrap();
            observed.set_iter_limit_hint(Some(1000));
            let actual = observed.solve(&bounds);
            let capture = take().unwrap();
            assert_eq!(actual.status, expected.status);
            assert_eq!(actual.objective.to_bits(), expected.objective.to_bits());
            assert_eq!(observed.iterations(), dense.iterations());
            assert_eq!(observed.basis_columns(), dense.basis_columns());
            assert_eq!(
                actual.x.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                expected.x.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
            );
            assert!(!capture.source_unavailable && !capture.artificial_mapping_unavailable);
            take();
            let mut sparse = SparseSimplex::new(&p, Tolerances::default()).unwrap();
            sparse.set_iter_limit_hint(Some(1000));
            let expected = sparse.solve(&bounds);
            reset();
            let mut observed = SparseSimplex::new(&p, Tolerances::default()).unwrap();
            observed.set_iter_limit_hint(Some(1000));
            let actual = observed.solve(&bounds);
            let capture = take().unwrap();
            assert_eq!(actual.status, expected.status);
            assert_eq!(actual.objective.to_bits(), expected.objective.to_bits());
            assert_eq!(observed.iterations(), sparse.iterations());
            assert_eq!(observed.basis_columns(), sparse.basis_columns());
            assert_eq!(
                actual.x.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                expected.x.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
            );
            assert!(!capture.source_unavailable);
        }
    }

    #[test]
    fn former_negative_counterexample_now_keeps_infeasibility_proof() {
        // 新安全比值修复不会再制造负人工；此处钉住当前公开接口正确状态。
        let p = problem(true);
        reset();
        let mut s = SparseSimplex::new(&p, Tolerances::default()).unwrap();
        let result = s.solve(&Bounds::from_problem(&p));
        assert_eq!(result.status, LpStatus::Infeasible);
        let capture = take().unwrap();
        assert_eq!(capture.terminal, "infeasible");
        assert!(capture
            .events
            .iter()
            .any(|e| e.kind == "terminal" && e.positive_sum > capture.tolerances.feasibility));
        assert!(!capture.snapshots.iter().any(|s| s.kind == "before_cleanup"));
    }

    #[test]
    fn injected_negative_entry_evidence_is_not_cleaned_away() {
        // 这里只验证守卫证据，不把内部注入冒充新核心的自然输入反例。
        let p = problem(false);
        reset();
        let mut s = Simplex::new(&p, Tolerances::default()).unwrap();
        s.rhs[0] = -2.0_f64.powi(-21);
        s.start_call();
        let status = s.run_two_phase(true);
        assert_eq!(status, LpStatus::IterationLimit);
        let capture = take().unwrap();
        let failure = capture
            .snapshots
            .iter()
            .find(|s| s.kind == "failure")
            .unwrap();
        assert_eq!(failure.summary.call_iterations, 0);
        assert!(failure
            .basic_values
            .iter()
            .any(|v| *v < -Tolerances::default().feasibility));
        assert!(!capture.snapshots.iter().any(|s| s.kind == "before_cleanup"));
    }

    #[test]
    fn first_candidate_and_reset_do_not_capture_later_clones() {
        let p = problem(false);
        reset();
        let mut s = Simplex::new(&p, Tolerances::default()).unwrap();
        let mut clone = s.clone();
        s.solve(&Bounds::from_problem(&p));
        clone.solve(&Bounds::from_problem(&p));
        let capture = take().unwrap();
        assert_eq!(
            capture
                .snapshots
                .iter()
                .filter(|s| s.kind == "first_candidate")
                .count(),
            1
        );
        assert!(capture.snapshots.len() <= 4);
        assert_eq!(capture.terminal, "phase1_cleaned");
        reset();
        clone.solve(&Bounds::from_problem(&p));
        let later = take().unwrap();
        assert!(later.source_unavailable, "旧引擎编号不得误认新reset的来源");
    }

    #[test]
    fn nonfinite_and_negative_zero_json_preserve_exact_bits() {
        let mut out = Vec::new();
        floats(
            &mut out,
            &[
                -0.0,
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::from_bits(0x7ff8000000000123),
            ],
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "[\"8000000000000000\",\"7ff0000000000000\",\"fff0000000000000\",\"7ff8000000000123\"]"
        );

        // 非有限人工行在有限违界行之前或之后，违界汇总都必须保留 NaN。
        let p = problem(false);
        for values in [[f64::NAN, -1.0], [-1.0, f64::NAN]] {
            reset();
            let mut solver = Simplex::new(&p, Tolerances::default()).unwrap();
            solver.rhs[..2].copy_from_slice(&values);
            solver.start_call();
            assert_eq!(solver.run_two_phase(true), LpStatus::IterationLimit);
            let capture = take().unwrap();
            let entry = &capture.snapshots[0];
            assert_eq!(entry.kind, "entry");
            assert_eq!(entry.summary.nonfinite, 1);
            assert!(entry.summary.basic_bound_violation.is_nan());
        }
    }

    #[test]
    fn source_budget_omits_before_clone_and_reports_it() {
        reset();
        STATE.with(|s| s.borrow_mut().bytes = BYTE_LIMIT);
        let p = problem(false);
        let mut engine = Simplex::new(&p, Tolerances::default()).unwrap();
        engine.solve(&Bounds::from_problem(&p));
        let capture = take().unwrap();
        assert!(capture.source_unavailable);
        assert!(capture.objects_omitted > 0);
        assert!(capture.estimated_bytes <= BYTE_LIMIT);
    }
    #[test]
    fn event_and_snapshot_truncation_are_visible_without_fabricating_values() {
        let p = problem(false);
        reset();
        let source = register_source(&p);
        let basis = [1usize];
        let xb = [0.0];
        let in_basis = [false, true];
        let at_upper = [false, false];
        let bounds = Bounds {
            lb: vec![0.0],
            ub: vec![f64::INFINITY],
        };
        let rows = vec![vec![(0usize, 1.0)]];
        let rhs = [0.0];
        let ops = [RelOp::Eq];
        let scales = [1.0];
        let slack = [None];
        let mapping = [Some(0usize)];
        let art = [(0usize, 1.0)];
        let _scope = begin(
            17,
            source,
            "test",
            MatrixView {
                tolerances: Tolerances::default(),
                n: 1,
                ncs: 1,
                terms: &rows,
                rhs: &rhs,
                op: &ops,
                scale: &scales,
                slack: &slack,
                constraint_row: &mapping,
                cut_row: &[],
                art: Some(&art),
                dense_initial: None,
                art_count: 1,
            },
        );
        let view = |iteration| View {
            n: 1,
            ncs: 1,
            basis: &basis,
            xb: &xb,
            in_basis: &in_basis,
            at_upper: &at_upper,
            bounds: &bounds,
            iterations: iteration,
            call_iterations: iteration,
            eta: 0,
            bland: false,
            feasibility: 1e-9,
            broken: false,
        };
        observe(17, "entry", view(0), None);
        STATE.with(|s| s.borrow_mut().bytes = BYTE_LIMIT);
        observe(17, "step", view(0), None);
        for k in 1..=EVENT_LIMIT + 2 {
            observe(17, "step", view(k as u64 * 1024), None);
        }
        let capture = take().unwrap();
        assert_eq!(capture.events.len(), EVENT_LIMIT);
        assert_eq!(capture.events.capacity(), EVENT_LIMIT);
        assert_eq!(capture.snapshots.capacity(), SNAPSHOT_LIMIT);
        assert!(
            FIXED_BYTES
                >= capture.events.capacity() * std::mem::size_of::<Summary>()
                    + capture.snapshots.capacity() * std::mem::size_of::<Snapshot>()
        );
        assert!(capture.events_omitted > 0);
        assert_eq!(capture.slot_attempt_iterations[1], Some(0));
        assert!(capture.slot_omitted[1]);
        assert!(!capture
            .snapshots
            .iter()
            .any(|s| s.kind == "first_candidate"));
        assert_eq!(capture.estimated_bytes, BYTE_LIMIT);
        let mut out = Vec::new();
        capture.write_json(&mut out).unwrap();
        let json = String::from_utf8(out).unwrap();
        assert!(json.contains("\"slot_omitted\":[false, true, false, false]"));
        assert!(
            json.contains("7ff0000000000000"),
            "合法无限上界保留原始位串"
        );
    }

    #[test]
    fn first_candidate_pricing_is_explicitly_same_state_and_old_error_is_not_terminal_reason() {
        let p = problem(false);
        reset();
        let mut s = Simplex::new(&p, Tolerances::default()).unwrap();
        s.solve(&Bounds::from_problem(&p));
        let capture = take().unwrap();
        let candidate = capture
            .snapshots
            .iter()
            .find(|s| s.kind == "first_candidate")
            .unwrap();
        let pricing = capture
            .events
            .iter()
            .find(|e| e.kind == "pricing_at_first_candidate")
            .unwrap();
        assert_eq!(pricing.state_revision, candidate.summary.state_revision);
        assert_eq!(
            pricing.pricing.evaluated_revision,
            Some(candidate.summary.state_revision)
        );
        reset();
        let mut s = Simplex::new(&p, Tolerances::default()).unwrap();
        s.numerical_error = Some("上一次调用的旧错误".to_string());
        s.set_time_deadline(Some(
            std::time::Instant::now() - std::time::Duration::from_secs(1),
        ));
        let result = s.solve(&Bounds::from_problem(&p));
        assert_eq!(result.status, LpStatus::TimeLimit);
        let capture = take().unwrap();
        assert_eq!(capture.terminal, "time_limit");
        assert!(capture.terminal_detail.is_none());
    }
}
