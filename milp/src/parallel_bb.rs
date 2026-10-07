//! 并行分支定界引擎（工作包：并行分支定界引擎与搜索控制层）。
//!
//! 架构：调用线程先完成与顺序 [`crate::branch_bound::branch_and_bound`]
//! 同构的根阶段（根 LP → 纯 LP 捷径 → 根割循环（记录割清单）→ 根启发式 →
//! resolve 回取），随后把根子节点（各自携带完整 [`Bounds`] 快照与 est =
//! 父节点目标）压入共享池，再经 [`std::thread::scope`] 派生 N 个 worker
//! 抢占式处理池节点。worker 各自持有独立 [`Simplex`]（Send、不共享），
//! 按序 `add_row` 重放根阶段记录的割清单，首节点 resolve 自动冷启动，
//! 任意界跳转由 snap/Phase-1 重建处理。
//!
//! 正确性骨架（与顺序路径逐字相同的判定）：
//! - 每节点恰处理一次（池互斥弹出）；已处理节点的全部子节点必入池
//!   （空界分支除外——与顺序 tighten_child 语义一致）⇒ 树覆盖无遗漏；
//! - 剪枝判定 (d)/(e)/(f) 与 [`crate::branch_bound`] 逐字移植（本模块私有
//!   副本，测试向量钉死）；incumbent 单调 ⇒ 剪枝永不丢更优解；
//! - stop 粘滞（CAS 一次性置位，镜像顺序 limit_hit 语义）：截断绝不冒充
//!   完备；gap 触发用触发时刻证据快照（粘滞、clamp 后非负有限）。
//!
//! gap 保守性：聚合 db 取『池内 est ∪ 在处理槽 est』的方向极值，est 为父
//! 谱系 LP 目标（对子树有效）；快照后 db 沿谱系只可能变严、incumbent 只
//! 可能更优 ⇒ 计算得到的 gap ≥ 真实 gap（N1 的 clamp 只把负 gap 提到 0，
//! 而 db ≤ inc 蕴含最优 ≤ inc，0% 触发诚实）⇒ 触发只推迟、不虚假。
//!
//! 锁纪律：任一线程任一时刻至多持一把锁，绝不嵌套；gap 检查两段式
//! （池锁内聚合 → 释放 → incumbent 读 → 释放 → 无锁计算 + CAS）；
//! 锁中毒 `into_inner` 恢复并置 stop；join Err（worker panic）映射
//! `MilpError::Internal`，不外泄 panic。
//!
//! 与顺序路径的已知行为差异（不影响正确性与最优值，README 口径）：
//! nodes/迭代数/耗时/取值（多最优时）跨运行可变。根割循环经
//! `cuts::cut_loop_recorded`（带记账、与顺序 `cuts::cut_loop` 逐字节一致），
//! worker 按序 `add_row` 重放记录的割清单 ⇒ 根割行为与顺序路径一致。

use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicU64, AtomicU8, Ordering};
use std::sync::{Condvar, Mutex, MutexGuard, RwLock};
use std::thread;
use std::time::{Duration, Instant};

use crate::branch_bound::{BbOptions, BbStats, BranchRule, NodeRule};
use crate::cuts::{self, CutOptions};
use crate::error::{MilpError, Result};
use crate::heuristics::{self, Incumbent};
use crate::model::{Problem, RelOp, Sense, SolveStatus, Tolerances};
use crate::search_ctrl::{
    aggregate_bound, gap_percent, AutoDegrade, DegradePolicy, DegradeSnapshot, HeuristicController,
    HeuristicDecision, ProgressCallback, ProgressInfo, SearchSnapshot,
};
use crate::simplex::{Bounds, LpSolution, LpStatus, Simplex};

/// stop 槽位：未置位。
const STOP_NONE: u8 = 0;
/// stop 槽位：节点上限（含迭代上限与异常中止的通用终止码；error 优先判定）。
const STOP_NODE: u8 = 1;
/// stop 槽位：时间上限。
const STOP_TIME: u8 = 2;
/// stop 槽位：gap 达到。
const STOP_GAP: u8 = 3;
/// stop 槽位：策略级自动降级触发（auto 模式 + 时间预算 + 双停滞）。
const STOP_DEGRADE: u8 = 4;

/// 并行参数。仅 derive(Clone)——含 `Arc<dyn Fn>` 字段不得 derive(Debug)
/// （rustc 实证 E0277：dyn Fn 不实现 Debug）。
#[derive(Clone)]
pub struct ParallelConfig {
    /// 工作线程数（≥ 2；1 由调用方走顺序路径）。
    pub threads: usize,
    /// mipgap 百分比；None 或 0 = 关闭。
    pub mipgap_percent: Option<f64>,
    /// 是否启用自适应启发式（auto：树节点 dive 由
    /// [`HeuristicController`] 决策；根启发式仍由 `BbOptions::heuristics`
    /// 控制）。
    pub adaptive_heuristics: bool,
    /// 进度行周期秒数；None = 不输出周期行（首行契约仍由 progress 回调
    /// 决定，见 R6）。
    pub log_period_secs: Option<f64>,
    /// 策略级自动降级参数（仅在 `adaptive_heuristics = true` 且
    /// `BbOptions::time_limit_ms != 0` 时评估；触发 ⇒ 与 Gap 同径提前终止，
    /// 终态映射为 [`ParallelStop::AutoDegrade`]）。
    pub degrade: DegradePolicy,
    /// 进度回调（首行在根阶段结束汇合点无条件发射）。
    pub progress: Option<ProgressCallback>,
}

/// 停止原因（粘滞：一旦置位终态按此判定，即使池随后排空也不得报完备——
/// 镜像 branch_bound 的 limit_hit 粘滞语义）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParallelStop {
    /// 节点上限（含 LP 迭代上限）。
    NodeLimit,
    /// 时间上限。
    TimeLimit,
    /// mipgap 达到。
    Gap,
    /// 策略级自动降级触发（auto 模式 + 时间预算 + 双停滞；证据口径同 Gap）。
    AutoDegrade,
}

/// 并行引擎结果（自有类型；终态映射由 branch_bound 侧完成）。
#[derive(Debug, Clone)]
pub struct ParallelOutcome {
    /// 根终局捷径直返；None = 进入过树搜索。
    pub root_terminal: Option<SolveStatus>,
    /// 停止原因；None = 树探索完备（池空 ∧ active == 0 且从未置 stop）。
    pub stop: Option<ParallelStop>,
    /// 最优目标（原始 sense + obj_offset 口径）。
    pub best_objective: Option<f64>,
    /// 最优解（与入参 problem 同维同序）。
    pub best_values: Option<Vec<f64>>,
    /// 统计：nodes = 根 + performed；lp_iterations = 根 + Σ worker；
    /// time_ms = 墙钟。
    pub stats: BbStats,
    /// 终局对偶界（原始 sense 口径）：
    /// - Some(Gap)：触发时刻证据 db（粘滞、有限）；
    /// - Some(NodeLimit|TimeLimit)：剩余池 est 聚合（可 None：截断不可知的
    ///   诚实 None）；
    /// - None（完备）：有 incumbent ⇒ best_objective、Infeasible ⇒ None；
    /// - root_terminal：纯 LP Optimal ⇒ Some(objective)。
    pub bound: Option<f64>,
    /// 终局 gap：
    /// - Some(Gap)：触发时刻证据值（clamp 后非负有限）；
    /// - 完备且有 incumbent ⇒ Some(0.0)；
    /// - 其余按 [`gap_percent`]（缺任一侧 ⇒ None）。
    pub gap_percent: Option<f64>,
}

// ===========================================================================
// 共享状态（模块私有）
// ===========================================================================

/// 池节点：完整界快照 + est（父 resolve 目标，原始 sense 口径，Optimal
/// 路径必有限）。
struct PoolNode {
    bounds: Bounds,
    est: f64,
}

/// 池状态。不变量：`active` == inflight 中 `Some` 槽位数。
struct PoolState {
    /// 待处理节点（尾弹出 = LIFO，下行分支先探索）。
    nodes: Vec<PoolNode>,
    /// 在处理节点数。
    active: usize,
    /// 搜索结束标志（完备或 stop 后置位）。
    done: bool,
    /// Some(est) = worker w 在处理节点（est 为该节点父谱系目标）；
    /// None = 空闲（聚合时剔除，无哨兵毒化）。
    inflight: Vec<Option<f64>>,
}

/// gap 触发证据（粘滞）：仅 stop CAS 0→Gap 的胜者写入一次（db、inc 均有限，
/// gap 经 clamp 非负）。`incumbent` 为触发时刻的 incumbent 溯源字段（契约
/// 草图成员；终局组装只消费 db/gap_percent，保留作证据记录）。
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
struct GapEvidence {
    db: f64,
    incumbent: f64,
    gap_percent: Option<f64>,
}

/// 共享状态。问题与割清单不入 Shared——worker 闭包经 thread::scope 直接
/// 借用调用栈上的 `&Problem` 与 `&[(terms, rhs)]`（无 'static 要求）。
struct Shared {
    pool: Mutex<PoolState>,
    /// worker 等节点：push → notify_one；done/stop → notify_all。
    cv_workers: Condvar,
    /// 仅协调线程 wait_timeout：仅 done/stop → notify_all。
    cv_progress: Condvar,
    incumbent: RwLock<Option<(f64, Vec<f64>)>>,
    /// 预约计数：根阶段后初始化 = 根节点数（= 1）——对齐顺序
    /// `stats.nodes >= node_limit`（根计 1）。
    reserved: AtomicU64,
    /// 已处理节点数（含根，同步初始化 = 1）。
    performed: AtomicU64,
    /// 最近一次 incumbent 改进所在节点。
    last_improve_node: AtomicU64,
    /// 最近一次 incumbent 改进的时刻（毫秒；降级停滞信号）。
    last_improve_millis: AtomicU64,
    /// 0=None 1=NodeLimit 2=TimeLimit 3=Gap 4=Degrade（CAS 一次性置位，粘滞）。
    stop: AtomicU8,
    gap_evidence: Mutex<Option<GapEvidence>>,
    controller: Mutex<HeuristicController>,
    /// 策略级自动降级判定器（窗口界基准有状态；锁内串行判定）。
    degrade: Mutex<AutoDegrade>,
    /// 降级参数副本（资格预检无锁先行）。
    degrade_cfg: DegradePolicy,
    /// 降级生效的时间预算；None = 关闭（无预算或非 auto 模式）。
    degrade_budget_secs: Option<f64>,
    /// 求解起点（降级判定的耗时基准）。
    start: Instant,
    error: Mutex<Option<MilpError>>,
}

impl Shared {
    #[allow(clippy::too_many_arguments)]
    fn new(
        threads: usize,
        degrade_policy: DegradePolicy,
        degrade_budget_secs: Option<f64>,
        start: Instant,
    ) -> Self {
        Shared {
            pool: Mutex::new(PoolState {
                nodes: Vec::new(),
                active: 0,
                done: false,
                inflight: vec![None; threads],
            }),
            cv_workers: Condvar::new(),
            cv_progress: Condvar::new(),
            incumbent: RwLock::new(None),
            reserved: AtomicU64::new(0),
            performed: AtomicU64::new(0),
            last_improve_node: AtomicU64::new(0),
            last_improve_millis: AtomicU64::new(0),
            stop: AtomicU8::new(STOP_NONE),
            gap_evidence: Mutex::new(None),
            controller: Mutex::new(HeuristicController::new(Default::default())),
            degrade: Mutex::new(AutoDegrade::new(degrade_policy)),
            degrade_cfg: degrade_policy,
            degrade_budget_secs,
            start,
            error: Mutex::new(None),
        }
    }

    // —— 锁获取（中毒 into_inner 恢复 + 置 stop）——

    fn lock_pool(&self) -> MutexGuard<'_, PoolState> {
        self.pool
            .lock()
            .unwrap_or_else(|poisoned| {
                self.abort();
                poisoned.into_inner()
            })
    }

    fn lock_controller(&self) -> MutexGuard<'_, HeuristicController> {
        self.controller
            .lock()
            .unwrap_or_else(|poisoned| {
                self.abort();
                poisoned.into_inner()
            })
    }

    fn lock_gap_evidence(&self) -> MutexGuard<'_, Option<GapEvidence>> {
        self.gap_evidence
            .lock()
            .unwrap_or_else(|poisoned| {
                self.abort();
                poisoned.into_inner()
            })
    }

    fn lock_degrade(&self) -> MutexGuard<'_, AutoDegrade> {
        self.degrade
            .lock()
            .unwrap_or_else(|poisoned| {
                self.abort();
                poisoned.into_inner()
            })
    }

    fn record_error(&self, error: MilpError) {
        let mut slot = self
            .error
            .lock()
            .unwrap_or_else(|poisoned| {
                self.abort();
                poisoned.into_inner()
            });
        if slot.is_none() {
            *slot = Some(error);
        }
    }

    fn take_error(&self) -> Option<MilpError> {
        let mut slot = self
            .error
            .lock()
            .unwrap_or_else(|poisoned| {
                self.abort();
                poisoned.into_inner()
            });
        slot.take()
    }

    fn read_incumbent_obj(&self) -> Option<f64> {
        let guard = self
            .incumbent
            .read()
            .unwrap_or_else(|poisoned| {
                self.abort();
                poisoned.into_inner()
            });
        guard.as_ref().map(|(obj, _)| *obj)
    }

    fn write_incumbent(&self, value: Option<(f64, Vec<f64>)>) {
        let mut guard = self
            .incumbent
            .write()
            .unwrap_or_else(|poisoned| {
                self.abort();
                poisoned.into_inner()
            });
        *guard = value;
    }

    /// 写锁双重检查更新 incumbent；返回是否成为新全局 incumbent。
    fn try_update_incumbent(&self, sense: Sense, inc: Incumbent, node_no: u64, elapsed_secs: f64) -> bool {
        if !inc.objective.is_finite() {
            return false; // 防御：非有限目标不入 incumbent
        }
        let mut guard = self
            .incumbent
            .write()
            .unwrap_or_else(|poisoned| {
                self.abort();
                poisoned.into_inner()
            });
        let better = match guard.as_ref() {
            None => true,
            Some((obj, _)) => match sense {
                Sense::Maximize => inc.objective > *obj,
                Sense::Minimize => inc.objective < *obj,
            },
        };
        if better {
            *guard = Some((inc.objective, inc.values));
            drop(guard);
            self.last_improve_node.store(node_no, Ordering::SeqCst);
            self.last_improve_millis
                .store((elapsed_secs * 1000.0) as u64, Ordering::SeqCst);
            true
        } else {
            false
        }
    }

    // —— stop（CAS 一次性置位 + 双 cv 唤醒）——

    fn set_stop(&self, code: u8) {
        if self
            .stop
            .compare_exchange(STOP_NONE, code, Ordering::SeqCst, Ordering::SeqCst)
            .is_ok()
        {
            self.notify_all();
        }
    }

    /// 异常中止（锁中毒等）：error 优先判定，stop 仅作通用终止信号。
    fn abort(&self) {
        self.set_stop(STOP_NODE);
    }

    fn stop_code(&self) -> u8 {
        self.stop.load(Ordering::SeqCst)
    }

    fn notify_all(&self) {
        self.cv_workers.notify_all();
        self.cv_progress.notify_all();
    }

    /// worker 等待节点（池空且 active > 0 且未 done/stop）。
    fn wait_workers<'a>(&'a self, guard: MutexGuard<'a, PoolState>) -> MutexGuard<'a, PoolState> {
        self.cv_workers
            .wait(guard)
            .unwrap_or_else(|poisoned| {
                self.abort();
                poisoned.into_inner()
            })
    }

    // —— 节点完成转移（(5)：push 持锁点 + 两段式 gap 检查）——

    /// 节点完成：压 children（如有）、释放持有槽（active−1、inflight=None）、
    /// 持锁聚合 db 快照、解锁后两段式 gap 检查。
    fn complete_node(
        &self,
        w: usize,
        mut children: Vec<PoolNode>,
        sense: Sense,
        mipgap: Option<f64>,
    ) {
        let db = {
            let mut pool = self.lock_pool();
            let grew = !children.is_empty();
            pool.nodes.append(&mut children);
            pool.active -= 1;
            pool.inflight[w] = None;
            if grew {
                // push ⇒ 唤醒一个等待 worker（完成者自身也会回头取节点）。
                self.cv_workers.notify_one();
            }
            pool_db_aggregate(&pool, sense)
        };
        self.maybe_gap(sense, mipgap, db);
        self.maybe_degrade(sense, db);
    }

    /// 两段式 gap 检查（调用方已释放全部锁）：读 incumbent → 无锁计算 →
    /// gap ≤ P 且 inc/db 有限 ⇒ stop CAS 0→Gap 胜者独占写证据 + 双 cv
    /// notify_all。
    fn maybe_gap(&self, sense: Sense, mipgap: Option<f64>, db: Option<f64>) {
        let Some(limit) = mipgap else {
            return;
        };
        if !limit.is_finite() || limit <= 0.0 {
            return; // None / 0 / 非有限 = 关闭
        }
        let (Some(db), Some(inc)) = (db, self.read_incumbent_obj()) else {
            return;
        };
        let Some(gap) = gap_percent(sense, inc, db) else {
            return;
        };
        if gap <= limit
            && self
                .stop
                .compare_exchange(STOP_NONE, STOP_GAP, Ordering::SeqCst, Ordering::SeqCst)
                .is_ok()
        {
            *self.lock_gap_evidence() = Some(GapEvidence {
                db,
                incumbent: inc,
                gap_percent: Some(gap),
            });
            self.notify_all();
        }
    }

    /// 策略级自动降级检查（与 maybe_gap 同在节点完成路径、同样两段式）：
    /// 先做无锁资格预检（degrade_budget_secs = None 即关闭），再读 incumbent
    /// （无嵌套锁）、持 degrade 锁串行判定（degrade 为叶子锁，不与池/incumbent
    /// 嵌套）；触发 ⇒ stop CAS 0→Degrade 胜者独占写 gap 证据（口径同
    /// STOP_GAP）+ 双 cv notify_all。终态映射见 [`ParallelStop::AutoDegrade`]。
    fn maybe_degrade(&self, sense: Sense, db: Option<f64>) {
        let Some(budget) = self.degrade_budget_secs else {
            return; // 未启用：无预算或非 auto 模式
        };
        let now = self.start.elapsed().as_secs_f64();
        if !AutoDegrade::eligible(self.degrade_cfg, now, Some(budget)) {
            return;
        }
        let snap = DegradeSnapshot {
            sense,
            elapsed_secs: now,
            budget_secs: Some(budget),
            has_incumbent: self.read_incumbent_obj().is_some(),
            last_improvement_secs: self.last_improve_millis.load(Ordering::SeqCst) as f64 / 1000.0,
            bound: db,
        };
        if !self.lock_degrade().decide(&snap) {
            return;
        }
        // decide = true 蕴含 db 有限且 incumbent 存在；incumbent 重读一次，
        // 期间再改进只会让 gap 更小（与 STOP_GAP 证据口径一致）。
        let Some(db) = db else {
            return; // 防御：不可达
        };
        if self
            .stop
            .compare_exchange(STOP_NONE, STOP_DEGRADE, Ordering::SeqCst, Ordering::SeqCst)
            .is_ok()
        {
            let inc = self.read_incumbent_obj();
            *self.lock_gap_evidence() = Some(GapEvidence {
                db,
                incumbent: inc.unwrap_or(f64::NAN),
                gap_percent: inc.and_then(|i| gap_percent(sense, i, db)),
            });
            self.notify_all();
        }
    }

    /// 退出路径（(6)）：仍持有则释放（active−1、inflight=None）；
    /// stop ≠ 0 或（空且 active == 0）⇒ done + 双 cv notify_all；
    /// 否则非空 ⇒ cv_workers.notify_one。
    fn worker_release(&self, w: usize) {
        let mut pool = self.lock_pool();
        if pool.inflight[w].is_some() {
            pool.active -= 1;
            pool.inflight[w] = None;
        }
        if self.stop_code() != STOP_NONE || (pool.nodes.is_empty() && pool.active == 0) {
            pool.done = true;
            self.notify_all();
        } else if !pool.nodes.is_empty() {
            self.cv_workers.notify_one();
        }
    }

    /// 决策快照（读锁取 incumbent、原子量取计数；无锁内嵌套）。
    fn snapshot(&self, time_budget_secs: Option<f64>, start: &Instant) -> SearchSnapshot {
        let nodes = self.performed.load(Ordering::SeqCst);
        SearchSnapshot {
            nodes,
            incumbent: self.read_incumbent_obj(),
            nodes_since_improvement: nodes.saturating_sub(self.last_improve_node.load(Ordering::SeqCst)),
            time_budget_secs,
            elapsed_secs: start.elapsed().as_secs_f64(),
        }
    }
}

/// 池对偶界聚合：池内 est ∪ inflight Some 槽（调用前滤除非有限值）。
fn pool_db_aggregate(pool: &PoolState, sense: Sense) -> Option<f64> {
    aggregate_bound(
        sense,
        pool.nodes
            .iter()
            .map(|n| n.est)
            .filter(|e| e.is_finite())
            .chain(
                pool.inflight
                    .iter()
                    .filter_map(|slot| *slot)
                    .filter(|e| e.is_finite()),
            ),
    )
}

// ===========================================================================
// 并行分支定界主入口
// ===========================================================================

/// 并行分支定界。前置：`config.threads >= 2`；`branch_rule`/`node_rule`
/// 非默认值 ⇒ `Err(Internal)`。
pub fn parallel_branch_and_bound(
    problem: &Problem,
    opts: &BbOptions,
    config: &ParallelConfig,
    tol: Tolerances,
) -> Result<ParallelOutcome> {
    let start = Instant::now();
    if config.threads < 2 {
        return Err(MilpError::Internal(format!(
            "parallel_branch_and_bound 前置失败：threads = {}（须 ≥ 2；threads = 1 由调用方走顺序路径）",
            config.threads
        )));
    }
    if opts.branch_rule != BranchRule::MostFractional || opts.node_rule != NodeRule::DepthFirst {
        return Err(MilpError::Internal(
            "parallel_branch_and_bound 前置失败：仅支持默认分支/节点规则（MostFractional/DepthFirst）".to_string(),
        ));
    }
    let int_vars = problem.integer_vars();
    let sense = problem.sense;
    let mut stats = BbStats::default();
    let mut simplex = Simplex::new(problem, tol)?;
    let bounds = Bounds::from_problem(problem);
    // incumbent = (目标, 取值)；目标一律按取整后的值评估（WP3 复核口径）。
    let mut incumbent: Option<(f64, Vec<f64>)> = None;
    // 根阶段施加的割清单（terms, rhs）——worker 按序 add_row 重放。
    let mut recorded_cuts: Vec<(Vec<(usize, f64)>, f64)> = Vec::new();
    let mut pure_lp = false;

    // —— 根 LP（计 1 节点）——
    let mut sol = simplex.solve(&bounds);
    stats.nodes = 1;

    let mut terminal: Option<SolveStatus> = None;
    if int_vars.is_empty() {
        // 纯 LP 捷径：无整数变量 ⇒ 单次 LP 直接按其状态返回。
        pure_lp = true;
        terminal = Some(match sol.status {
            LpStatus::Optimal => {
                incumbent = Some((sol.objective, sol.x.clone()));
                SolveStatus::Optimal
            }
            LpStatus::Infeasible => SolveStatus::Infeasible,
            LpStatus::Unbounded => SolveStatus::Unbounded,
            LpStatus::IterationLimit => SolveStatus::LimitReached,
        });
    } else {
        // 根 solve 非 Optimal ⇒ 不进割循环与启发式，直接终局。
        if let Some(status) = root_terminal_status(sol.status) {
            terminal = Some(status);
        }
        if terminal.is_none() && opts.root_cuts {
            let (applied, st, cuts_list) = cuts::cut_loop_recorded(
                &mut simplex,
                &bounds,
                problem,
                &int_vars,
                CutOptions::default(),
                &tol,
            )?;
            stats.cuts = applied;
            recorded_cuts = cuts_list
                .into_iter()
                .map(|cut| (cut.terms, cut.rhs))
                .collect();
            if let Some(status) = root_terminal_status(st) {
                terminal = Some(status);
            }
            if terminal.is_none() {
                // 取回割后的根 LP 解（记账性 resolve，不计节点）。
                sol = simplex.resolve(&bounds);
                if let Some(status) = root_terminal_status(sol.status) {
                    terminal = Some(status);
                }
            }
        }
        if terminal.is_none() && opts.heuristics {
            // 根启发式：round_and_check + dive 供应初始 incumbent。
            if let Some(inc) = heuristics::round_and_check(problem, &bounds, &sol.x, &tol) {
                maybe_update_local(sense, &mut incumbent, inc);
            }
            if let Some(inc) = heuristics::dive(
                problem,
                &mut simplex,
                &bounds,
                &sol.x,
                opts.dive_max_fixes,
                &tol,
            ) {
                maybe_update_local(sense, &mut incumbent, inc);
            }
            // dive 扰动了 Simplex 内部基：真实界 resolve 回取（记账性）。
            sol = simplex.resolve(&bounds);
            if let Some(status) = root_terminal_status(sol.status) {
                terminal = Some(status);
            }
        }
    }

    // —— 首行进度无条件发射点：根阶段结束的唯一汇合点（R6）——
    if let Some(callback) = &config.progress {
        let info = ProgressInfo {
            sense,
            elapsed_secs: start.elapsed().as_secs_f64(),
            nodes: 1,
            incumbent: incumbent
                .as_ref()
                .map(|(obj, _)| *obj)
                .filter(|v| v.is_finite()),
            bound: if sol.status == LpStatus::Optimal && sol.objective.is_finite() {
                Some(sol.objective)
            } else {
                None
            },
        };
        callback(&info);
    }

    // —— 根终局直返 ——
    if let Some(status) = terminal {
        stats.lp_iterations = simplex.iterations();
        stats.time_ms = wall_ms(&start);
        let (best_objective, best_values) = match incumbent {
            Some((obj, values)) => (Some(obj), Some(values)),
            None => (None, None),
        };
        let bound = if pure_lp && status == SolveStatus::Optimal {
            best_objective
        } else {
            None
        };
        let gap = gap_percent(
            sense,
            best_objective.unwrap_or(f64::NAN),
            bound.unwrap_or(f64::NAN),
        );
        return Ok(ParallelOutcome {
            root_terminal: Some(status),
            stop: None,
            best_objective,
            best_values,
            stats,
            bound,
            gap_percent: gap,
        });
    }

    // —— 根节点 (d)/(e)/(f) 判定（私有快照版，无 trail）——
    let mut children: Vec<PoolNode> = Vec::new();
    {
        // (d)：LP bound 不优于 incumbent ⇒ 根关闭。
        let mut closed = false;
        if let Some((obj, _)) = incumbent.as_ref() {
            if bound_no_better(
                sense,
                sol.objective,
                *obj,
                opts.prune_abs,
                opts.prune_rel,
            ) {
                closed = true;
            }
        }
        if !closed {
            // (e)：整数性满足 ⇒ round_and_check 复核；通过 ⇒ 记 incumbent
            // 后关闭；不通过（伪整数）⇒ 落入分支。
            let integral = is_integral(&int_vars, &sol.x, tol.integrality);
            let mut check_passed = false;
            if integral {
                if let Some(inc) = heuristics::round_and_check(problem, &bounds, &sol.x, &tol) {
                    maybe_update_local(sense, &mut incumbent, inc);
                    check_passed = true;
                }
            }
            if !(integral && check_passed) {
                // (f)：分支压 children（est = 根目标）。
                children = branch_children(&int_vars, &bounds, &sol);
            }
        }
    }

    // —— 共享状态初始化（根阶段末尾：reserved/performed 同步置 1）——
    let shared = Shared::new(
        config.threads,
        config.degrade,
        if config.adaptive_heuristics && opts.time_limit_ms > 0 {
            Some(opts.time_limit_ms as f64 / 1000.0)
        } else {
            None
        },
        start,
    );
    shared.write_incumbent(incumbent.clone());
    shared.reserved.store(1, Ordering::SeqCst);
    shared.performed.store(1, Ordering::SeqCst);
    shared
        .last_improve_node
        .store(u64::from(incumbent.is_some()), Ordering::SeqCst);
    if incumbent.is_some() {
        shared.last_improve_millis.store(
            u64::try_from(start.elapsed().as_millis()).unwrap_or(u64::MAX),
            Ordering::SeqCst,
        );
    }
    shared.lock_controller().notify_root_done(1);

    // —— 根子节点入池【持锁】+ 两段式 gap/降级检查（spawn 之前）——
    if !children.is_empty() {
        let db = {
            let mut pool = shared.lock_pool();
            pool.nodes.append(&mut children);
            pool_db_aggregate(&pool, sense)
        };
        shared.maybe_gap(sense, config.mipgap_percent, db);
        shared.maybe_degrade(sense, db);
    }

    // 根 gap/降级触发 ⇒ 直接收尾，不 spawn。
    if matches!(shared.stop_code(), STOP_GAP | STOP_DEGRADE) {
        stats.lp_iterations = simplex.iterations();
        stats.time_ms = wall_ms(&start);
        let (best_objective, best_values) = match incumbent {
            Some((obj, values)) => (Some(obj), Some(values)),
            None => (None, None),
        };
        let evidence = *shared.lock_gap_evidence();
        let (bound, gap) = match evidence {
            Some(ev) => (Some(ev.db), ev.gap_percent),
            None => (None, None), // 防御：CAS 胜者必写证据
        };
        return Ok(ParallelOutcome {
            root_terminal: None,
            stop: Some(if shared.stop_code() == STOP_GAP {
                ParallelStop::Gap
            } else {
                ParallelStop::AutoDegrade
            }),
            best_objective,
            best_values,
            stats,
            bound,
            gap_percent: gap,
        });
    }

    // 根即闭合（无子节点）⇒ 完备收尾，无需 spawn。
    if shared.lock_pool().nodes.is_empty() {
        stats.lp_iterations = simplex.iterations();
        stats.time_ms = wall_ms(&start);
        let (best_objective, best_values) = match incumbent {
            Some((obj, values)) => (Some(obj), Some(values)),
            None => (None, None),
        };
        return Ok(ParallelOutcome {
            root_terminal: None,
            stop: None,
            best_objective,
            best_values,
            stats,
            bound: best_objective,
            gap_percent: best_objective.map(|_| 0.0),
        });
    }

    // —— 并行树搜索（仅当 stop 仍为 None）——
    let mipgap = config.mipgap_percent;
    let adaptive = config.adaptive_heuristics;
    let time_budget_secs = if opts.time_limit_ms > 0 {
        Some(opts.time_limit_ms as f64 / 1000.0)
    } else {
        None
    };
    let worker_opts = *opts;
    let mut worker_iters: u64 = 0;
    let mut join_failed = false;
    let progress_cb = config.progress.clone();
    let log_period = config.log_period_secs;

    // 提升共享借用：move 闭包按值捕获这些 & 引用（引用本身 Copy），避免
    // 整个所有权被移入首个闭包。
    let problem_ref = problem;
    let int_vars_ref = &int_vars;
    let cuts_ref = &recorded_cuts;
    let shared_ref = &shared;

    thread::scope(|scope| {
        let mut handles = Vec::with_capacity(config.threads);
        for w in 0..config.threads {
            handles.push(scope.spawn(move || {
                worker_run(
                    w,
                    problem_ref,
                    int_vars_ref,
                    cuts_ref,
                    shared_ref,
                    worker_opts,
                    tol,
                    mipgap,
                    adaptive,
                    time_budget_secs,
                    start,
                )
            }));
        }
        // 协调线程（scope 内的调用线程）：进度节奏。
        if let (Some(callback), Some(period)) = (&progress_cb, log_period) {
            if period.is_finite() && period > 0.0 {
                progress_loop(shared_ref, callback, period, &start, sense);
            }
        }
        // 逐个 join：worker panic ⇒ Internal（不外泄 panic）。
        for handle in handles {
            match handle.join() {
                Ok(iters) => worker_iters = worker_iters.wrapping_add(iters),
                Err(_) => join_failed = true,
            }
        }
    });

    if join_failed {
        return Err(MilpError::Internal(
            "并行 worker 线程 panic（join 失败）".to_string(),
        ));
    }
    if let Some(error) = shared.take_error() {
        return Err(error);
    }

    // —— 组装（Gap 用粘滞证据，绝不在 join 后重新聚合）——
    stats.nodes = shared.performed.load(Ordering::SeqCst);
    stats.lp_iterations = simplex
        .iterations()
        .wrapping_add(worker_iters);
    stats.time_ms = wall_ms(&start);
    let (best_objective, best_values) = {
        let guard = shared
            .incumbent
            .read()
            .unwrap_or_else(|poisoned| {
                shared.abort();
                poisoned.into_inner()
            });
        match guard.as_ref() {
            Some((obj, values)) => (Some(*obj), Some(values.clone())),
            None => (None, None),
        }
    };
    let (stop_out, bound, gap) = match shared.stop_code() {
        code @ (STOP_GAP | STOP_DEGRADE) => {
            let evidence = *shared.lock_gap_evidence();
            (
                Some(if code == STOP_GAP {
                    ParallelStop::Gap
                } else {
                    ParallelStop::AutoDegrade
                }),
                evidence.map(|ev| ev.db),
                evidence.and_then(|ev| ev.gap_percent),
            )
        }
        code @ (STOP_NODE | STOP_TIME) => {
            let stop = Some(if code == STOP_NODE {
                ParallelStop::NodeLimit
            } else {
                ParallelStop::TimeLimit
            });
            let bound = pool_db_aggregate(&shared.lock_pool(), sense);
            (
                stop,
                bound,
                gap_percent(
                    sense,
                    best_objective.unwrap_or(f64::NAN),
                    bound.unwrap_or(f64::NAN),
                ),
            )
        }
        _ => {
            // 完备：有 incumbent ⇒ bound = best、gap = 0；Infeasible ⇒ 均无。
            (
                None,
                best_objective,
                best_objective.map(|_| 0.0),
            )
        }
    };
    Ok(ParallelOutcome {
        root_terminal: None,
        stop: stop_out,
        best_objective,
        best_values,
        stats,
        bound,
        gap_percent: gap,
    })
}

// ===========================================================================
// worker 主循环
// ===========================================================================

/// worker 线程主体：返回本 worker 的累计 LP 迭代数（join 汇总）。
#[allow(clippy::too_many_arguments)]
fn worker_run(
    w: usize,
    problem: &Problem,
    int_vars: &[usize],
    cuts: &[(Vec<(usize, f64)>, f64)],
    shared: &Shared,
    opts: BbOptions,
    tol: Tolerances,
    mipgap: Option<f64>,
    adaptive: bool,
    time_budget_secs: Option<f64>,
    start: Instant,
) -> u64 {
    let sense = problem.sense;
    let mut simplex = match Simplex::new(problem, tol) {
        Ok(s) => s,
        Err(error) => {
            shared.record_error(error);
            shared.set_stop(STOP_NODE);
            return 0;
        }
    };
    // 按序重放根阶段记录的割清单（与根线程 add_row 顺序一致）。
    for (terms, rhs) in cuts {
        if let Err(error) = simplex.add_row(terms, RelOp::Le, *rhs) {
            shared.record_error(error);
            shared.set_stop(STOP_NODE);
            return 0;
        }
    }
    loop {
        // —— 取节点（持池锁；空且 active > 0 且未 done ⇒ cv_workers.wait）——
        let node = {
            let mut pool = shared.lock_pool();
            loop {
                if shared.stop_code() != STOP_NONE {
                    drop(pool);
                    shared.worker_release(w);
                    return simplex.iterations();
                }
                if let Some(node) = pool.nodes.pop() {
                    pool.active += 1;
                    pool.inflight[w] = Some(node.est);
                    break node;
                }
                if pool.active == 0 {
                    // 池空且无在处理节点：搜索完备。
                    pool.done = true;
                    shared.notify_all();
                    drop(pool);
                    return simplex.iterations();
                }
                pool = shared.wait_workers(pool);
            }
        };
        // —— 预约制上限检查（成功弹出后；触发即丢节点：不 resolve、不计
        //    performed——镜像顺序『pop 后检查、触发即丢』语义）——
        let reserve = shared.reserved.fetch_add(1, Ordering::SeqCst) + 1;
        if opts.node_limit != 0 && reserve > opts.node_limit {
            shared.set_stop(STOP_NODE);
            shared.worker_release(w);
            return simplex.iterations();
        }
        if opts.time_limit_ms != 0 && wall_exceeded(&start, opts.time_limit_ms) {
            shared.set_stop(STOP_TIME);
            shared.worker_release(w);
            return simplex.iterations();
        }
        // —— resolve（计 1 performed）——
        let sol = simplex.resolve(&node.bounds);
        let node_no = shared.performed.fetch_add(1, Ordering::SeqCst) + 1;
        match sol.status {
            // (a) Infeasible ⇒ 剪枝。
            LpStatus::Infeasible => {
                shared.complete_node(w, Vec::new(), sense, mipgap);
                continue;
            }
            // (b) 非根 Unbounded（父有界时理论不发生的数值异常）⇒ 诚实 Err。
            LpStatus::Unbounded => {
                shared.record_error(nonroot_unbounded_error());
                shared.set_stop(STOP_NODE);
                shared.worker_release(w);
                return simplex.iterations();
            }
            // (c) IterationLimit ⇒ 置 stop=NodeLimit（粘滞）。
            LpStatus::IterationLimit => {
                shared.set_stop(STOP_NODE);
                shared.complete_node(w, Vec::new(), sense, mipgap);
                continue;
            }
            LpStatus::Optimal => {}
        }
        // —— (d)/(e)/adaptive/(f) ——
        let mut sol = sol;
        let mut children: Vec<PoolNode> = Vec::new();
        let mut closed = false;
        // (d)：LP bound 不优于 incumbent（读锁取当前值）⇒ 节点关闭。
        if let Some(obj) = shared.read_incumbent_obj() {
            if bound_no_better(sense, sol.objective, obj, opts.prune_abs, opts.prune_rel) {
                closed = true;
            }
        }
        // (e)：整数性满足 ⇒ round_and_check 复核；通过 ⇒ 写锁双重检查更新
        // incumbent + last_improve_node 复位；伪整数 ⇒ 继续分支。
        if !closed {
            let integral = is_integral(int_vars, &sol.x, tol.integrality);
            let mut check_passed = false;
            if integral {
                if let Some(inc) = heuristics::round_and_check(problem, &node.bounds, &sol.x, &tol)
                {
                    shared.try_update_incumbent(sense, inc, node_no, start.elapsed().as_secs_f64());
                    check_passed = true;
                }
            }
            closed = integral && check_passed;
        }
        // adaptive：控制器 decide = Run ⇒ dive + resolve 回取 + record(hit)。
        if !closed && adaptive {
            let decision = {
                let snap = shared.snapshot(time_budget_secs, &start);
                let controller = shared.lock_controller();
                controller.decide(&snap)
            };
            if decision == HeuristicDecision::Run {
                let hit = match heuristics::dive(
                    problem,
                    &mut simplex,
                    &node.bounds,
                    &sol.x,
                    opts.dive_max_fixes,
                    &tol,
                ) {
                    Some(inc) => {
                        shared.try_update_incumbent(sense, inc, node_no, start.elapsed().as_secs_f64())
                    }
                    None => false,
                };
                // dive 扰动了基：真实界 resolve 回取（记账性，不计节点）。
                let back = simplex.resolve(&node.bounds);
                match back.status {
                    LpStatus::Optimal => sol = back,
                    LpStatus::Infeasible => closed = true,
                    LpStatus::IterationLimit => shared.set_stop(STOP_NODE),
                    LpStatus::Unbounded => {
                        shared.record_error(nonroot_unbounded_error());
                        shared.set_stop(STOP_NODE);
                        shared.worker_release(w);
                        return simplex.iterations();
                    }
                }
                shared.lock_controller().record(hit);
            }
        }
        // (f)：分支压 children（est = 本节点目标）。
        if !closed {
            children = branch_children(int_vars, &node.bounds, &sol);
        }
        // —— 节点完成转移 + 两段式 gap 检查 ——
        shared.complete_node(w, children, sense, mipgap);
    }
}

// ===========================================================================
// 协调线程（进度节奏）
// ===========================================================================

/// 进度节奏：**接管即发首周期行**（`next_deadline` 初值 0，首轮流转立即
/// 到期），此后 cv_progress.wait_timeout(池锁, 至下一 deadline)；被唤醒后
/// 仅当未 done 且距上行 ≥ period 才快照（两段式）+ 回调；done/stop 后只
/// 退出不打印（快照时必有池非空或 active > 0，bound 不会空聚合）。
///
/// 入场即发的理由：若首周期行须等满一个 period，亚毫秒级搜索与极小
/// period（如 0.001s）将比拼定时精度——wait_timeout 的实际唤醒受 OS 定时
/// 粒度量化（Windows 默认 tick ≈ 15.6ms），短于一个 tick 的等待可能整体
/// 错过 deadline：搜索先完成 ⇒ notify_all 唤醒后只见 done 退出，『搜索期间
/// 每 F 秒一行』的契约在最需要观测的短搜索上反而一行皆无。接管时刻树搜索
/// 必然进行中（spawn 刚返回、根子节点已入池；done 未置 ⇒ 池非空或
/// active > 0），首行快照口径与后续周期行一致且不会空聚合出 bound=none
///（N5 不变）；首行契约（R6）的根阶段行在 spawn 之前发射，先后次序不变。
fn progress_loop(
    shared: &Shared,
    callback: &ProgressCallback,
    period: f64,
    start: &Instant,
    sense: Sense,
) {
    // 0 ⇒ 首轮流转立即到期：协调线程接管节奏即发首周期行。
    let mut next_deadline = 0.0f64;
    let mut pool = shared.lock_pool();
    loop {
        if pool.done || shared.stop_code() != STOP_NONE {
            break;
        }
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed + 1e-9 >= next_deadline {
            // 两段式快照：池锁内聚合 db → 解锁 → 读 incumbent → 回调
            //（回调期间不持有求解器任何锁）。
            let db = pool_db_aggregate(&pool, sense);
            drop(pool);
            let info = ProgressInfo {
                sense,
                elapsed_secs: elapsed,
                nodes: shared.performed.load(Ordering::SeqCst),
                incumbent: shared.read_incumbent_obj().filter(|v| v.is_finite()),
                bound: db,
            };
            // 回调 panic 不得外泄（scope 闭包内）：诚实记 error 并终止。
            if catch_unwind(AssertUnwindSafe(|| callback(&info))).is_err() {
                shared.record_error(MilpError::Internal(
                    "progress 回调 panic".to_string(),
                ));
                shared.abort();
                return;
            }
            next_deadline = elapsed + period;
            pool = shared.lock_pool();
            // 回调期间搜索可能已结束：重取锁后回到顶部 done 检查再排队等待。
            continue;
        }
        let wait_secs = (next_deadline - elapsed).max(0.0);
        let (guard, _) = shared
            .cv_progress
            .wait_timeout(pool, Duration::from_secs_f64(wait_secs))
            .unwrap_or_else(|poisoned| {
                shared.abort();
                poisoned.into_inner()
            });
        pool = guard;
    }
}

// ===========================================================================
// 根阶段辅助：终局分类与记录式割循环
// ===========================================================================

/// 根节点 LP 状态 → 整树终态（`None` = Optimal，继续树搜索）。
fn root_terminal_status(status: LpStatus) -> Option<SolveStatus> {
    match status {
        LpStatus::Optimal => None,
        LpStatus::Infeasible => Some(SolveStatus::Infeasible),
        LpStatus::Unbounded => Some(SolveStatus::Unbounded),
        LpStatus::IterationLimit => Some(SolveStatus::LimitReached),
    }
}

/// 非根节点 LP 无界的数值错误（镜像 branch_bound 同名语义）。
fn nonroot_unbounded_error() -> MilpError {
    MilpError::Numerical("node LP unbounded after bounded root".to_string())
}

// ===========================================================================
// 节点判定纯函数（自 branch_bound.rs 逐字移植的私有副本，测试向量钉死）
// ===========================================================================

/// (d) 剪枝判定（按原始 sense 方向）：
/// Max: bound ≤ inc + prune_abs + prune_rel·|inc|；Min 镜像。
fn bound_no_better(sense: Sense, bound: f64, inc: f64, prune_abs: f64, prune_rel: f64) -> bool {
    let slack = prune_abs + prune_rel * inc.abs();
    match sense {
        Sense::Maximize => bound <= inc + slack,
        Sense::Minimize => bound >= inc - slack,
    }
}

/// 整数性判定：全部 int_var 满足 |x−round(x)| ≤ integrality。
fn is_integral(int_vars: &[usize], x: &[f64], integrality: f64) -> bool {
    int_vars.iter().all(|&j| match x.get(j) {
        Some(v) => (v - v.round()).abs() <= integrality,
        None => false,
    })
}

/// MostFractional 分支变量选择：max min(f, 1−f)（f = x−⌊x⌋），并列取最小下标。
fn select_branch_var(int_vars: &[usize], x: &[f64]) -> Option<(usize, f64)> {
    let mut best: Option<(usize, f64, f64)> = None; // (var, 值, 得分)
    for &j in int_vars {
        let Some(&v) = x.get(j) else {
            continue;
        };
        if !v.is_finite() {
            continue;
        }
        let f = v - v.floor();
        let score = f.min(1.0 - f);
        let replace = match best {
            None => true,
            Some((_, _, s)) => score > s, // 严格大于 ⇒ 并列保留最小下标
        };
        if replace {
            best = Some((j, v, score));
        }
    }
    best.map(|(j, v, _)| (j, v))
}

/// 施加子节点分支界（快照版：对节点自有 Bounds 副本直接收紧，无 trail）：
/// 相对当前界严格收缩且新区间非空才应用，返回 true；空界（(f)）或无收缩
/// （无进展）返回 false，调用方不压栈（不 resolve、不计节点）。
fn tighten_child_snapshot(bounds: &mut Bounds, var: usize, tighten_lb: bool, bound: f64) -> bool {
    if tighten_lb {
        if bound <= bounds.lb[var] || bound > bounds.ub[var] {
            return false;
        }
        bounds.lb[var] = bound;
    } else {
        if bound >= bounds.ub[var] || bound < bounds.lb[var] {
            return false;
        }
        bounds.ub[var] = bound;
    }
    true
}

/// 分支压栈（(f)）：MostFractional 选变量，生成上下行两个子节点
/// （上行先压、下行后压——LIFO 弹出时下行先探索，与顺序栈序一致）；
/// est = 本节点 resolve 目标（对子树有效的乐观界）。
fn branch_children(int_vars: &[usize], bounds: &Bounds, sol: &LpSolution) -> Vec<PoolNode> {
    // 无可分支变量为内部不可达防御（不剪枝即丢弃节点）。
    let Some((var, val)) = select_branch_var(int_vars, &sol.x) else {
        return Vec::new();
    };
    let floor = val.floor();
    let mut children = Vec::new();
    for (tighten_lb, bound) in [(true, floor + 1.0), (false, floor)] {
        let mut child = bounds.clone();
        if tighten_child_snapshot(&mut child, var, tighten_lb, bound) {
            children.push(PoolNode {
                bounds: child,
                est: sol.objective,
            });
        }
    }
    children
}

/// 严格更优才替换 incumbent（按原始 sense 方向；目标与取值均为 WP3 取整后
/// 口径，防报不可行值）——根阶段本地版。
fn maybe_update_local(sense: Sense, incumbent: &mut Option<(f64, Vec<f64>)>, inc: Incumbent) {
    if !inc.objective.is_finite() {
        return; // 防御：非有限目标不入 incumbent
    }
    let better = match incumbent.as_ref() {
        None => true,
        Some((obj, _)) => match sense {
            Sense::Maximize => inc.objective > *obj,
            Sense::Minimize => inc.objective < *obj,
        },
    };
    if better {
        *incumbent = Some((inc.objective, inc.values));
    }
}

/// 时间上限检查（调用方保证 `limit_ms > 0`）。
fn wall_exceeded(start: &Instant, limit_ms: u64) -> bool {
    start.elapsed().as_millis() >= limit_ms as u128
}

/// 墙钟毫秒（溢出饱和，不 panic）。
fn wall_ms(start: &Instant) -> u64 {
    u64::try_from(start.elapsed().as_millis()).unwrap_or(u64::MAX)
}

// ===========================================================================
// 单元测试（除注明外全部确定性；定位算例文件用 env!("CARGO_MANIFEST_DIR")）
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::branch_bound::branch_and_bound;
    use crate::instance::load_instance;
    use crate::model::VarType;
    use std::path::Path;
    use std::sync::{Arc, Mutex as StdMutex};

    fn tol() -> Tolerances {
        Tolerances::default()
    }

    fn assert_close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 1e-6 * (1.0 + b.abs()),
            "数值不符: {a} vs {b}"
        );
    }

    fn instance_path(name: &str) -> String {
        format!("{}/instances/{name}.milp", env!("CARGO_MANIFEST_DIR"))
    }

    /// 并行配置（自适应关、无进度行）。
    fn pconfig(threads: usize, mipgap: Option<f64>) -> ParallelConfig {
        ParallelConfig {
            threads,
            mipgap_percent: mipgap,
            adaptive_heuristics: false,
            log_period_secs: None,
            degrade: DegradePolicy::default(),
            progress: None,
        }
    }

    /// dive 钉死问题（heuristics.rs:241-263）：max 10x1+x2 s.t. 20x1+3x2≤43，
    /// x∈[0,10] 整数；根 LP 21.5，根 dive 得 21。
    fn dive_problem() -> Problem {
        let mut p = Problem::new("dive", Sense::Maximize);
        p.add_var("x1", VarType::Integer).unwrap();
        p.add_var("x2", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_var_bounds(1, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 20.0), (1, 3.0)], RelOp::Le, 43.0)
            .unwrap();
        p
    }

    fn dive_opts() -> BbOptions {
        BbOptions {
            root_cuts: false,
            heuristics: true,
            ..BbOptions::default()
        }
    }

    /// 前置失败：threads < 2 ⇒ Err(Internal)。
    #[test]
    fn preconditions_rejected() {
        let p = dive_problem();
        for threads in [0usize, 1] {
            let r = parallel_branch_and_bound(&p, &BbOptions::default(), &pconfig(threads, None), tol());
            assert!(
                matches!(r, Err(MilpError::Internal(_))),
                "threads = {threads} 应拒绝"
            );
        }
    }

    /// (a) 根 gap 触发（spawn 之前完成 ⇒ 与线程数无关地确定）：
    /// mipgap=3.0 ⇒ gap≈2.381% ≤ 3 ⇒ stop=Gap、best=21、bound=21.5（证据）、
    /// gap_percent=≈2.38。
    #[test]
    fn gap_triggers_at_root_before_spawn() {
        let p = dive_problem();
        let expected_gap = 100.0 * 0.5 / 21.0;
        for threads in [2usize, 4] {
            let out = parallel_branch_and_bound(
                &p,
                &dive_opts(),
                &pconfig(threads, Some(3.0)),
                tol(),
            )
            .unwrap();
            assert_eq!(out.root_terminal, None, "threads = {threads}");
            assert_eq!(out.stop, Some(ParallelStop::Gap), "threads = {threads}");
            assert_close(out.best_objective.unwrap(), 21.0);
            assert_close(out.bound.unwrap(), 21.5);
            let g = out.gap_percent.unwrap();
            assert!(
                (g - expected_gap).abs() <= 1e-9,
                "gap = {g}，期望 ≈ {expected_gap}"
            );
            assert_eq!(out.stats.nodes, 1, "根阶段后即触发，无 worker 节点");
        }
    }

    /// (a) mipgap=2.0 ⇒ 根 gap≈2.381% > 2 不触发 ⇒ 树完备：
    /// stop=None、best=21、bound=Some(21)、gap=Some(0.0)。
    #[test]
    fn gap_not_triggered_completes() {
        let p = dive_problem();
        for threads in [2usize, 4] {
            let out = parallel_branch_and_bound(
                &p,
                &dive_opts(),
                &pconfig(threads, Some(2.0)),
                tol(),
            )
            .unwrap();
            assert_eq!(out.root_terminal, None, "threads = {threads}");
            assert_eq!(out.stop, None, "threads = {threads}：2.381% > 2% 不触发");
            assert_close(out.best_objective.unwrap(), 21.0);
            assert_close(out.bound.unwrap(), 21.0);
            assert_eq!(out.gap_percent, Some(0.0), "完备且有 incumbent ⇒ 0.0");
        }
    }

    /// (b) node_limit=2：child1 预约号 r=2（处理，performed=2）、child2 r=3
    /// （>2 丢弃 + stop=NodeLimit）——与哪颗子节点被预约无关地确定。
    #[test]
    fn node_limit_two_stops_deterministically() {
        let p = dive_problem();
        let opts = BbOptions {
            node_limit: 2,
            ..dive_opts()
        };
        for threads in [2usize, 4] {
            let out =
                parallel_branch_and_bound(&p, &opts, &pconfig(threads, None), tol()).unwrap();
            assert_eq!(out.stop, Some(ParallelStop::NodeLimit), "threads = {threads}");
            assert_close(
                out.best_objective.unwrap(),
                21.0,
            );
            assert_eq!(out.stats.nodes, 2, "根 + 恰一颗子节点");
        }
    }

    /// (c) 5 个原验收算例：threads=2/4 与顺序 branch_and_bound 最优值一致
    ///（≤1e-6），终局完备。
    #[test]
    fn acceptance_objectives_match_sequential() {
        const CASES: &[(&str, f64)] = &[
            ("lp-prod-3", 13.0),
            ("knapsack-10", 17.0),
            ("assign-4x4", 10.0),
            ("milp-wolsey-2", 40.0),
            ("milp-mixed-5", 29.0),
        ];
        for &(name, expected) in CASES {
            let problem =
                load_instance(Path::new(&instance_path(name))).expect("算例加载失败");
            let seq =
                branch_and_bound(&problem, &BbOptions::default(), tol()).expect("顺序求解失败");
            assert_eq!(seq.status, SolveStatus::Optimal, "算例 {name}");
            assert_close(seq.best_objective.unwrap(), expected);
            for threads in [2usize, 4] {
                let out = parallel_branch_and_bound(
                    &problem,
                    &BbOptions::default(),
                    &pconfig(threads, None),
                    tol(),
                )
                .unwrap_or_else(|e| panic!("算例 {name} threads={threads} 求解失败：{e}"));
                // lp-prod-3 无整数变量 ⇒ 纯 LP 捷径（root_terminal=Some(Optimal)）；
                // 其余 MILP 算例进入树搜索（root_terminal=None）。
                let expect_root_terminal = if problem.integer_vars().is_empty() {
                    Some(SolveStatus::Optimal)
                } else {
                    None
                };
                assert_eq!(
                    out.root_terminal, expect_root_terminal,
                    "算例 {name} threads={threads}"
                );
                assert_eq!(out.stop, None, "算例 {name} threads={threads} 应完备");
                let par_obj = out
                    .best_objective
                    .unwrap_or_else(|| panic!("算例 {name} threads={threads} 完备但无 incumbent"));
                assert!(
                    (par_obj - seq.best_objective.unwrap()).abs() <= 1e-6,
                    "算例 {name} threads={threads}：parallel = {par_obj}，sequential = {}",
                    seq.best_objective.unwrap()
                );
            }
        }
    }

    /// (c) infeasible/unbounded 根终局透传。
    #[test]
    fn status_instances_root_terminal_passthrough() {
        let problem = load_instance(Path::new(&instance_path("infeasible-tiny"))).unwrap();
        for threads in [2usize, 4] {
            let out = parallel_branch_and_bound(
                &problem,
                &BbOptions::default(),
                &pconfig(threads, None),
                tol(),
            )
            .unwrap();
            assert_eq!(out.root_terminal, Some(SolveStatus::Infeasible));
            assert!(out.best_objective.is_none());
            assert!(out.best_values.is_none());
            assert_eq!(out.stop, None);
        }
        let problem = load_instance(Path::new(&instance_path("unbounded-tiny"))).unwrap();
        for threads in [2usize, 4] {
            let out = parallel_branch_and_bound(
                &problem,
                &BbOptions::default(),
                &pconfig(threads, None),
                tol(),
            )
            .unwrap();
            assert_eq!(out.root_terminal, Some(SolveStatus::Unbounded));
            assert!(out.best_objective.is_none());
            assert!(out.best_values.is_none());
            assert_eq!(out.stop, None);
        }
    }

    /// 纯 LP Optimal 根终局：bound = Some(objective)、gap = Some(0.0)。
    #[test]
    fn pure_lp_root_terminal_bound_and_gap() {
        // max 3x+2y s.t. x+y≤4, x+3y≤6, [0,3]^2 → (3,1)，目标 11。
        let mut p = Problem::new("lp", Sense::Maximize);
        p.add_var("x", VarType::Continuous).unwrap();
        p.add_var("y", VarType::Continuous).unwrap();
        p.set_var_bounds(0, 0.0, 3.0).unwrap();
        p.set_var_bounds(1, 0.0, 3.0).unwrap();
        p.set_objective(vec![(0, 3.0), (1, 2.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 4.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 3.0)], RelOp::Le, 6.0)
            .unwrap();
        let out = parallel_branch_and_bound(&p, &BbOptions::default(), &pconfig(2, None), tol())
            .unwrap();
        assert_eq!(out.root_terminal, Some(SolveStatus::Optimal));
        assert_close(out.best_objective.unwrap(), 11.0);
        assert_close(out.bound.unwrap(), 11.0);
        assert_eq!(out.gap_percent, Some(0.0));
        assert_eq!(out.stats.nodes, 1);
    }

    /// (d) inflight None 槽不参与聚合（单元级：直接构造池状态）。
    #[test]
    fn pool_aggregation_skips_none_and_nonfinite_inflight() {
        let pool = PoolState {
            nodes: vec![
                PoolNode {
                    bounds: Bounds {
                        lb: vec![0.0],
                        ub: vec![1.0],
                    },
                    est: 21.5,
                },
                PoolNode {
                    bounds: Bounds {
                        lb: vec![0.0],
                        ub: vec![1.0],
                    },
                    est: 3.0,
                },
            ],
            active: 1,
            done: false,
            inflight: vec![None, Some(2.0), Some(f64::NAN)],
        };
        assert_eq!(pool_db_aggregate(&pool, Sense::Maximize), Some(21.5));
        assert_eq!(pool_db_aggregate(&pool, Sense::Minimize), Some(2.0));
        // 全空 ⇒ None。
        let empty = PoolState {
            nodes: Vec::new(),
            active: 0,
            done: false,
            inflight: vec![None, None],
        };
        assert_eq!(pool_db_aggregate(&empty, Sense::Maximize), None);
    }

    /// hard-knapsack-12（12 变量 0-1 背包，与 WP4 落盘文件同参）：
    /// 2^12 穷举钉死最优 523；threads=2/4 各重复 100 次，每次断言完备且
    /// |best−523| ≤ 1e-6（丢唤醒/死锁/池/incumbent/gap 竞态在重复下必现；
    /// 单次求解亚毫秒，循环耗时可忽略）。
    #[test]
    fn hard_knapsack_12_stress_100x() {
        let weights = [
            7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0, 31.0, 37.0, 41.0, 43.0, 47.0,
        ];
        let values = [
            11.0, 20.0, 32.0, 44.0, 55.0, 70.0, 88.0, 99.0, 111.0, 122.0, 133.0, 144.0,
        ];
        let mut p = Problem::new("hard-knapsack-12", Sense::Maximize);
        for i in 0..12 {
            p.add_var(&format!("x{i}"), VarType::Binary).unwrap();
        }
        p.set_objective((0..12).map(|i| (i, values[i])).collect(), 0.0)
            .unwrap();
        p.add_constraint(
            None,
            (0..12).map(|i| (i, weights[i])).collect(),
            RelOp::Le,
            172.0,
        )
        .unwrap();

        // 穷举钉死最优（4096 个子集）。
        let mut best = 0.0f64;
        for mask in 0u32..(1 << 12) {
            let (mut w, mut v) = (0.0f64, 0.0f64);
            for i in 0..12 {
                if mask & (1 << i) != 0 {
                    w += weights[i];
                    v += values[i];
                }
            }
            if w <= 172.0 {
                best = best.max(v);
            }
        }
        assert_eq!(best, 523.0, "穷举最优（与工作包说明一致）");

        let base = ParallelConfig {
            mipgap_percent: None,
            adaptive_heuristics: true, // auto 默认
            log_period_secs: None,
            progress: None,
            ..pconfig(2, None)
        };
        for threads in [2usize, 4] {
            let config = ParallelConfig {
                threads,
                mipgap_percent: base.mipgap_percent,
                adaptive_heuristics: base.adaptive_heuristics,
                log_period_secs: base.log_period_secs,
                degrade: DegradePolicy::default(),
                progress: base.progress.clone(),
            };
            for rep in 0..100 {
                let out = parallel_branch_and_bound(&p, &BbOptions::default(), &config, tol())
                    .unwrap_or_else(|e| panic!("threads={threads} rep={rep} 求解失败：{e}"));
                assert_eq!(
                    out.stop,
                    None,
                    "threads={threads} rep={rep} 应完备（丢唤醒/竞态会在此暴露）"
                );
                assert_eq!(out.root_terminal, None);
                let obj = out.best_objective.unwrap();
                assert!(
                    (obj - 523.0).abs() <= 1e-6,
                    "threads={threads} rep={rep}：objective = {obj}，期望 523"
                );
            }
        }
    }

    /// 首行进度（R6）：根阶段结束汇合点无条件发射——根 gap 触发（不
    /// spawn）时恰好一行，字段口径 nodes=1 / incumbent=21 / bound=21.5 /
    /// gap=2.38%。
    #[test]
    fn first_progress_line_emitted_at_root_junction() {
        let p = dive_problem();
        let lines: Arc<StdMutex<Vec<String>>> = Arc::new(StdMutex::new(Vec::new()));
        let sink = lines.clone();
        let config = ParallelConfig {
            threads: 2,
            mipgap_percent: Some(3.0),
            adaptive_heuristics: false,
            log_period_secs: Some(1.0),
            degrade: DegradePolicy::default(),
            progress: Some(Arc::new(move |info: &ProgressInfo| {
                sink.lock()
                    .unwrap_or_else(|poisoned| poisoned.into_inner())
                    .push(crate::search_ctrl::format_progress_line(info));
            })),
        };
        let out = parallel_branch_and_bound(&p, &dive_opts(), &config, tol()).unwrap();
        assert_eq!(out.stop, Some(ParallelStop::Gap));
        let lines = lines
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        assert_eq!(lines.len(), 1, "gap 触发于 spawn 前 ⇒ 恰好首行：{lines:?}");
        assert!(
            lines[0].starts_with("progress t="),
            "行前缀：{}",
            lines[0]
        );
        assert!(
            lines[0].contains("nodes=1 incumbent=21.000000 bound=21.500000 gap=2.38%"),
            "字段口径：{}",
            lines[0]
        );
    }

    /// 根割清单重放一致性（并行正确性前提）：cuts::cut_loop_recorded 在
    /// wolsey 上必产割；全新 Simplex 按序 add_row 记录的割后，根 LP 目标与
    /// 记录方一致（worker 表格复现根割表格）。
    #[test]
    fn cut_loop_recorded_replay_consistency() {
        let build = || {
            let mut p = Problem::new("wolsey", Sense::Maximize);
            p.add_var("x1", VarType::Integer).unwrap();
            p.add_var("x2", VarType::Integer).unwrap();
            p.set_objective(vec![(0, 5.0), (1, 8.0)], 0.0).unwrap();
            p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 6.0)
                .unwrap();
            p.add_constraint(None, vec![(0, 5.0), (1, 9.0)], RelOp::Le, 45.0)
                .unwrap();
            p
        };
        let problem = build();
        let bounds = Bounds::from_problem(&problem);
        let int_vars = problem.integer_vars();
        let t = tol();

        // 记录方：根 solve + 记录式割循环。
        let mut a = Simplex::new(&problem, t).unwrap();
        let root = a.solve(&bounds);
        assert_eq!(root.status, LpStatus::Optimal);
        assert!((root.objective - 41.25).abs() < 1e-9);
        let (applied, st, cuts_list) = cuts::cut_loop_recorded(
            &mut a,
            &bounds,
            &problem,
            &int_vars,
            CutOptions::default(),
            &t,
        )
        .unwrap();
        assert_eq!(st, LpStatus::Optimal);
        assert!(applied >= 1, "wolsey 根 LP 41.25 处必产割");
        assert_eq!(cuts_list.len() as u64, applied);
        let after_a = a.resolve(&bounds);
        assert_eq!(after_a.status, LpStatus::Optimal);

        // 重放方（worker 路径）：全新 Simplex + 按序 add_row + resolve。
        let mut b = Simplex::new(&problem, t).unwrap();
        b.solve(&bounds);
        for cut in &cuts_list {
            b.add_row(&cut.terms, RelOp::Le, cut.rhs).unwrap();
        }
        let after_b = b.resolve(&bounds);
        assert_eq!(after_b.status, LpStatus::Optimal);
        assert!(
            (after_a.objective - after_b.objective).abs() <= 1e-9,
            "重放目标 {} 与记录方 {} 不一致",
            after_b.objective,
            after_a.objective
        );
    }

    // —— 移植纯函数的测试向量（自 branch_bound.rs 既有向量钉死）——

    /// (d) 剪枝比较按原始 sense 方向（Max/Min 镜像）。
    #[test]
    fn bound_prune_respects_sense() {
        let (pa, pr) = (1e-9, 1e-9);
        assert!(bound_no_better(Sense::Maximize, 10.0, 10.0, pa, pr));
        assert!(!bound_no_better(Sense::Maximize, 10.5, 10.0, pa, pr));
        assert!(bound_no_better(Sense::Minimize, 2.0, 2.0, pa, pr));
        assert!(!bound_no_better(Sense::Minimize, 1.5, 2.0, pa, pr));
    }

    /// MostFractional 选择与并列最小下标破键。
    #[test]
    fn most_fractional_selection_and_tie_break() {
        let iv = [0usize, 1, 2];
        assert_eq!(select_branch_var(&iv, &[1.5, 2.2, 0.8]), Some((0, 1.5)));
        assert_eq!(select_branch_var(&iv, &[1.25, 3.75, 0.0]), Some((0, 1.25)));
        assert_eq!(select_branch_var(&iv, &[2.0, 1.3, 4.0]), Some((1, 1.3)));
        assert_eq!(select_branch_var(&iv, &[-1.5, 0.0, 0.0]), Some((0, -1.5)));
        assert_eq!(select_branch_var(&[], &[1.5]), None);
    }

    /// 整数性判定：容差内视为整数；下标越界视为不整数。
    #[test]
    fn is_integral_rules() {
        let iv = [0usize, 1];
        assert!(is_integral(&iv, &[2.0, 3.0], 1e-6));
        assert!(is_integral(&iv, &[2.0 + 5e-7, 0.0], 1e-6));
        assert!(!is_integral(&iv, &[2.5, 0.0], 1e-6));
        assert!(!is_integral(&iv, &[], 1e-6), "取值缺失 ⇒ 非整数");
        assert!(is_integral(&[], &[2.5], 1e-6), "无整数变量 ⇒ 整数");
    }

    /// tighten_child_snapshot：空界/无收缩拒绝，正常收紧接受。
    #[test]
    fn tighten_child_snapshot_rules() {
        let mut b = Bounds {
            lb: vec![0.0, 0.4],
            ub: vec![5.0, 0.6],
        };
        assert!(
            tighten_child_snapshot(&mut b, 0, false, 2.0),
            "下行收紧 ub 至 2"
        );
        assert_eq!(b.ub[0], 2.0);
        assert!(!tighten_child_snapshot(&mut b, 0, false, 2.0), "无收缩");
        // 新副本：上行收紧 lb 至 3。
        let mut c = Bounds {
            lb: vec![0.0],
            ub: vec![5.0],
        };
        assert!(
            tighten_child_snapshot(&mut c, 0, true, 3.0),
            "上行收紧 lb 至 3"
        );
        assert_eq!(c.lb[0], 3.0);
        assert!(!tighten_child_snapshot(&mut c, 0, true, 6.0), "越过 ub 区间空");
        assert!(
            tighten_child_snapshot(&mut c, 0, true, 5.0),
            "等于 ub 收敛到单点区间 [5,5]（与顺序路径同判：仅严格越过才空）"
        );
        assert_eq!(c.lb[0], 5.0);
        // 分数区间 [0.4, 0.6]：两侧分支界均空/无收缩（既有算例语义）。
        assert!(!tighten_child_snapshot(&mut b, 1, false, 0.0), "0 ≤ lb 无收缩");
        assert!(!tighten_child_snapshot(&mut b, 1, true, 1.0), "1 > ub 区间空");
    }

    /// 策略级自动降级（并行路径；激进参数 start=0/plateau=0 ⇒ 首个节点
    /// 完成即触发）：stop = AutoDegrade、携带 incumbent 21、gap 证据非负
    /// 有限（同 STOP_GAP 口径，由 CAS 胜者写入）。
    #[test]
    fn parallel_auto_degrade_stops_with_incumbent() {
        let p = dive_problem();
        let config = ParallelConfig {
            threads: 4,
            mipgap_percent: None,
            adaptive_heuristics: true,
            log_period_secs: None,
            degrade: DegradePolicy {
                start_frac: 0.0,
                plateau_frac: 0.0,
                ..DegradePolicy::default()
            },
            progress: None,
        };
        // 关根割：GMI 割会把该问题根节点割成整数解直接闭合（无子节点入
        // 池），降级检查点（根入池/节点完成）无从触发；关掉以复现「树搜索
        // 中途双停滞」形态（与顺序侧 solve_controlled_auto_degrade 同配置）。
        let mut opts = BbOptions {
            root_cuts: false,
            heuristics: true,
            ..BbOptions::default()
        };
        opts.time_limit_ms = 10_000;
        let out = parallel_branch_and_bound(&p, &opts, &config, tol())
            .unwrap_or_else(|e| panic!("求解失败：{e}"));
        assert_eq!(out.root_terminal, None);
        assert_eq!(
            out.stop,
            Some(ParallelStop::AutoDegrade),
            "激进参数 + 预算 ⇒ 必须触发降级（实际 {:?}）",
            out.stop
        );
        let obj = out.best_objective.expect("降级必携带 incumbent");
        assert!(
            (obj - 21.0).abs() <= 1e-6,
            "incumbent = {obj}，期望根 dive 供应的 21"
        );
        let gap = out.gap_percent.expect("降级 ⇒ gap 证据必为 Some");
        assert!(gap.is_finite() && gap >= 0.0, "gap = {gap}");
    }
}
