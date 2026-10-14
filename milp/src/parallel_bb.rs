//! 并行分支定界引擎（工作包：并行分支定界引擎与搜索控制层）。
//!
//! 架构：调用线程先完成与顺序 [`crate::branch_bound::branch_and_bound`]
//! 同构的根阶段（根 LP → 纯 LP 捷径 → 根割循环（记录割清单）→ 根启发式 →
//! resolve 回取），随后把根子节点（各自携带完整 [`Bounds`] 快照与 est =
//! 父节点目标）压入全局窃取池，再经 [`std::thread::scope`] 派生 N 个
//! worker。worker 各自持有独立 [`Simplex`]（Send、不共享），按序 `add_row`
//! 重放根阶段记录的割清单，首节点 resolve 自动冷启动，任意界跳转由
//! snap/Phase-1 重建处理。
//!
//! 节点归属（并行亲和，单一大锁保护）：每个存活节点恰归属一个容器——
//! 某 worker 的本地 `VecDeque` 队列、全局窃取池、或某 worker 的 inflight
//! 槽，全部转移在池锁临界区内原子完成（无丢无重）。worker 分支时把两个
//! 孩子按序压**自己队首**（owner 随后队首弹出即得近亲子代——同亲热启动
//! 使对偶修复迭代数贴近顺序 DFS 档），本地空时先窃取池尾、再从随机起点
//! 轮询各 victim **队尾**偷单个（最老 deferred 分支；owner 正在探索的近亲
//! 前沿不动，首版不拆半）。victim 本地 push 不发全局 notify——空闲 worker
//! 经 `wait_timeout(≤5ms)` 周期醒来重扫窃取；stop/done/窃取池 push 走
//! `notify_all` 快路径立即唤醒。
//!
//! 正确性骨架（与顺序路径逐字相同的判定）：
//! - 每节点恰处理一次（锁内原子弹出，无丢无重）；已处理节点的全部子节点
//!   必入 owner 本地队列（空界分支除外——与顺序 tighten_child 语义一致）
//!   ⇒ 树覆盖无遗漏；
//! - 剪枝判定 (d)/(e)/(f) 与 [`crate::branch_bound`] 同口径：bound_no_better
//!   / is_integral 为共享 pub(crate) 实现（含整数目标 cutoff 的
//!   effective_prune_bound 合成，三调用点同口径；原逐字私有副本已删除）；
//!   陈旧前沿 est 弹出剪枝在池弹出后、reserve 之前（不 reserve、不
//!   performed——与顺序侧 pop 后跳过同口径）；incumbent 单调 ⇒ 剪枝永不
//!   丢更优解；
//! - 根割循环的 IterationLimit 分流与顺序侧同款：割后 resolve 触发迭代
//!   上限 ⇒ 丢弃全部根割（worker 不重放）、重建无割 Simplex 恢复 Optimal
//!   根态后继续搜索（有效割不切任何整数可行解 ⇒ 回退可行域只会更大；
//!   重建前迭代经 base_iters 累计不丢，stats.cuts 如实清零）；
//! - Pseudocost 分支（默认规则）：worker 各持本地 [`Pseudocosts`]（无全局
//!   锁）——子节点 resolve 后按 PoolNode 记账三元组（var/tighten_lb/
//!   frac_dir）以顺序侧同口径被动记账；浅层节点经共享的
//!   [`crate::branch_bound::probe_reliability`] 冷启动（同口径预算/深度门，
//!   探测只读界并恢复，不参与节点判定）；子代入池序与顺序侧 P1-5
//!   best-estimate 压栈序同口径翻转；
//! - 自适应 dive 的记账性回取 resolve：状态不参与节点判定（dive 前的
//!   resolve 已证实节点 LP 可行/最优）——Infeasible/IterationLimit 弃用
//!   结果继续分支（绝不 closed 丢子树，见 [`dive_back_action`]），
//!   Unbounded 诚实 Err；
//! - stop 粘滞（CAS 一次性置位）与 LP 迭代上限剪枝标记 limit_hit（镜像
//!   顺序 limit_hit：剪枝 + 粘滞 + 树继续，终态诚实 LimitReached）：截断
//!   绝不冒充完备；gap 触发用触发时刻证据快照（粘滞、clamp 后非负有限）。
//!
//! gap 保守性：聚合 db 取『全部本地队列的全部节点 est ∪ 窃取池全部 est ∪
//! 在处理槽 est』（**全量口径**， [`pool_db_aggregate`] 唯一实现）的方向
//! 极值，est 为父谱系 LP 目标（对子树有效）；快照后 db 沿谱系只可能变严、
//! incumbent 只可能更优 ⇒ 计算得到的 gap ≥ 真实 gap（N1 的 clamp 只把负
//! gap 提到 0，而 db ≤ inc 蕴含最优 ≤ inc，0% 触发诚实）⇒ 触发只推迟、
//! 不虚假。全量扫描成本被 complete_node 门控控制（仅 mipgap 有限正值或
//! degrade 预算存在时聚合）；严禁退回『各队列尾+inflight』口径——漏窃取
//! 池且依赖未写明的 est 单调不变量，低估开放节点极值使 gap 偏小属不保守
//! 方向。
//!
//! 锁纪律：任一线程任一时刻至多持一把锁，绝不嵌套；单一大锁（池锁）覆盖
//! 全部本地队列 + 窃取池 + active/inflight/done（4 线程竞争可控，规避分锁
//! +终止检测复杂度）；gap 检查两段式（池锁内聚合 → 释放 → incumbent 读
//! → 释放 → 无锁计算 + CAS）；锁中毒 `into_inner` 恢复并置 stop；join Err
//! （worker panic）映射 `MilpError::Internal`，不外泄 panic。
//!
//! 与顺序路径的已知行为差异（不影响正确性与最优值，README 口径）：
//! nodes/迭代数/耗时/取值（多最优时）跨运行可变。根割循环经
//! `cuts::cut_loop_recorded`（带记账、与顺序 `cuts::cut_loop` 逐字节一致），
//! worker 按序 `add_row` 重放记录的割清单 ⇒ 根割行为与顺序路径一致。

use std::collections::VecDeque;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicU8, Ordering};
use std::sync::{Condvar, Mutex, MutexGuard, RwLock};
use std::thread;
use std::time::{Duration, Instant};

use crate::branch_bound::{bound_no_better, is_integral, probe_reliability, BbOptions, BbStats, BranchRule, NodeRule};
use crate::cuts::{self, CutOptions};
use crate::error::{MilpError, Result};
use crate::heuristics::{self, Incumbent};
use crate::model::{objective_is_integral, Problem, RelOp, Sense, SolveStatus, Tolerances};
use crate::search_ctrl::{
    aggregate_bound, gap_percent, AutoDegrade, DegradePolicy, DegradeSnapshot, HeuristicController,
    HeuristicDecision, ProgressCallback, ProgressInfo, Pseudocosts, SearchSnapshot,
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

/// worker 等待/窃取轮询周期上限：本地队列与窃取池均空且仍有在处理节点时，
/// `cv_workers.wait_timeout(5ms)` 周期醒来从随机起点重扫各 victim 队尾——
/// 解决『victim 本地队列 push 不发全局 notify 则偷不到』的永久睡眠风险；
/// stop/done/窃取池 push 走 `notify_all` 快路径立即唤醒，常态不依赖本超时。
const WORKER_WAIT_TIMEOUT: Duration = Duration::from_millis(5);

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
/// 镜像 branch_bound 的 limit_hit 粘滞语义；节点 LP 迭代上限经独立粘滞
/// 标记 limit_hit 剪枝续搜、不截断树，终局经 [`effective_stop_code`] 同样
/// 映射到 NodeLimit ⇒ LimitReached，见 M3）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParallelStop {
    /// 节点上限（含节点 LP 迭代上限：剪枝续搜，探索完备后仍诚实 LimitReached）。
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
    /// 停止原因；None = 树探索完备（池空 ∧ active == 0 且从未置 stop，也
    /// 从未置 LP 迭代上限剪枝标记 limit_hit——置位时按 NodeLimit 类截断
    /// 诚实报告 LimitReached，不冒充完备）。
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
/// 路径必有限）+ 分支记账三元组（var/tighten_lb/frac_dir，供子节点 resolve
/// 后按顺序侧同口径做伪代价被动记账）+ 深度（可靠性探测的深度门控）。
struct PoolNode {
    bounds: Bounds,
    est: f64,
    /// 父节点分支决策的变量。
    var: usize,
    /// true = 上行分支（父节点收紧 lb 至 ⌊x̄⌋+1）；false = 下行（ub 至 ⌊x̄⌋）。
    tighten_lb: bool,
    /// 该方向距被违反界的分数距离（与顺序侧 Pending.frac_dir 同口径）。
    frac_dir: f64,
    /// 节点深度（根子节点 = 1）。
    depth: u32,
}

/// 池状态（**单一大锁**保护：全部本地队列 + 窃取池 + active/inflight/done，
/// 4 线程竞争可控，规避分锁+终止检测复杂度）。不变量：`active` == inflight
/// 中 `Some` 槽位数；每个存活节点恰归属一个容器（某本地队列 / 窃取池 /
/// 某 worker 的 inflight 槽），全部转移在锁临界区内原子完成 ⇒ 无丢无重。
struct PoolState {
    /// **全局 DFS 栈模式**（R6-retry，轴3）：true ⇒ 子代不进本地队列而
    /// 直接压全局 steal 栈顶（Vec push），worker 恒从栈顶弹（Vec pop =
    /// LIFO）——全局探索序 ≈ 顺序 DFS 序（谁弹到谁延续 t1 的下一个节点），
    /// 对齐 t1 的 incumbent 时序与剪枝，抑制并行搜索膨胀（local 亲和模式
    /// 下 4 棵独立 DFS 树各自延迟发现 incumbent）。false = 本地亲和模式
    ///（默认，现状）。
    work_global: bool,
    /// per-worker 本地队列（双端：**队首** = owner 弹出端 = 最新子代——分支时
    /// 两个孩子按序 `push_front` 到自己队首，owner 随后 `pop_front` 即得近亲
    /// 子代，同亲热启动保亲和；**队尾** = 最老 deferred 分支，留给窃取者）。
    locals: Vec<VecDeque<PoolNode>>,
    /// 全局窃取池（Vec 尾弹出 = LIFO）：根子节点入池；worker 本地空时优先
    /// 从池尾窃取。push 走 `notify_all` 快路径。
    steal: Vec<PoolNode>,
    /// 在处理节点数。
    active: usize,
    /// 搜索结束标志（完备或 stop 后置位）。
    done: bool,
    /// Some(est) = worker w 在处理节点（est 为该节点父谱系目标）；
    /// None = 空闲（聚合时剔除，无哨兵毒化）。
    inflight: Vec<Option<f64>>,
}

impl PoolState {
    /// 全部待处理容器是否为空：**全部本地队列 ∧ 窃取池**（done 判定口径，
    /// 漏任一容器都会提前 done 丢树）。
    fn all_empty(&self) -> bool {
        self.steal.is_empty() && self.locals.iter().all(|q| q.is_empty())
    }

    /// 取节点（调用方持池锁；弹出即从容器原子移除，无丢无重）：
    /// ① 本地队首（owner 最新子代，保亲和 LIFO）→ ② 窃取池尾 →
    /// ③ 从随机起点轮询各 victim **队尾**偷**单个**（首版不拆半：victim 只
    ///    失去一个最老 deferred 分支，owner 正在探索的近亲前沿不动；被偷
    ///    节点成为窃取者自家子树的种子，亲和按 worker 各自积累）。
    /// 全部为空返回 None（调用方判定终止或 wait_timeout）。
    fn take_local_or_steal(&mut self, w: usize, rng: &mut u64) -> Option<PoolNode> {
        if let Some(node) = self.locals[w].pop_front() {
            return Some(node);
        }
        if let Some(node) = self.steal.pop() {
            return Some(node);
        }
        let threads = self.locals.len().max(1);
        let start = (splitmix_next(rng) % threads as u64) as usize;
        for k in 0..threads {
            let v = (start + k) % threads;
            if v == w {
                continue;
            }
            if let Some(node) = self.locals[v].pop_back() {
                return Some(node);
            }
        }
        None
    }
}

/// SplitMix64 单步（确定性伪随机，无墙钟依赖）：窃取扫描的随机起点。
fn splitmix_next(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// worker 取节点结果：Node = 成功取到（active/inflight 已在锁内记账）；
/// Finished = 搜索结束（stop 置位、或 done——含本线程在锁内判定的完备）。
enum AcquireOutcome {
    Node(PoolNode),
    Finished,
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
    /// worker 等待/窃取：本地与窃取池皆空且 active>0 ⇒ wait_timeout(≤5ms)
    /// 周期醒来重扫 victim 队尾；窃取池 push / done / stop → notify_all
    /// （快路径立即唤醒）；本地队列 push **静默**（亲和保护——唤醒的睡者
    /// 会与 owner 争夺刚压入的近亲子代，靠超时兜底窃取）。
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
    /// LP 迭代上限剪枝的粘滞标记（镜像顺序路径 limit_hit）：worker 节点
    /// resolve 迭代上限 ⇒ 剪枝该节点并置位，**不终止树搜索**；终局若 stop
    /// 仍为 NONE 且本标记置位 ⇒ 诚实映射 NodeLimit 类截断（LimitReached）。
    limit_hit: AtomicBool,
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
                work_global: false,
                locals: (0..threads).map(|_| VecDeque::new()).collect(),
                steal: Vec::new(),
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
            limit_hit: AtomicBool::new(false),
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

    /// worker 取节点（单一大锁内：本地队首 → 窃取池尾 → 随机起点扫各
    /// victim 队尾偷单个）。暂无节点且仍有在处理 ⇒ `cv_workers.
    /// wait_timeout(≤5ms)` 周期醒来重扫（victim 本地 push 不发全局 notify，
    /// 靠本超时兜底窃取，杜绝永久睡眠）；stop/done/窃取池 push 的
    /// `notify_all` 走快路径立即唤醒。
    ///
    /// 完备判定（全部本地队列空 ∧ 窃取池空 ∧ active == 0 ⇒ done）与空
    /// 聚合机会在同一池锁临界区内原子发生（M1 契约，见 [`complete_node`]）：
    /// 否则『末节点解锁』与『某 worker 醒来置 done』之间存在窗口，进度线程
    /// 恰在此快照会输出 incumbent=X bound=none gap=inf 的矛盾行。
    fn acquire_node(&self, w: usize, rng: &mut u64) -> AcquireOutcome {
        let mut pool = self.lock_pool();
        loop {
            if pool.done || self.stop_code() != STOP_NONE {
                return AcquireOutcome::Finished;
            }
            if let Some(node) = pool.take_local_or_steal(w, rng) {
                pool.active += 1;
                pool.inflight[w] = Some(node.est);
                return AcquireOutcome::Node(node);
            }
            if pool.active == 0 {
                // 全部容器空（take 失败已蕴含本地队列与窃取池皆空）且无在
                // 处理节点：搜索完备。done 与（门控开启时的）空聚合同临界区。
                pool.done = true;
                self.notify_all();
                return AcquireOutcome::Finished;
            }
            let (guard, _) = self
                .cv_workers
                .wait_timeout(pool, WORKER_WAIT_TIMEOUT)
                .unwrap_or_else(|poisoned| {
                    self.abort();
                    poisoned.into_inner()
                });
            pool = guard;
        }
    }

    // —— 节点完成转移（(5)：push 持锁点 + 两段式 gap 检查）——

    /// 节点完成：children 压入 owner 本地队列队首、释放持有槽（active−1、
    /// inflight=None）、（门控开启时）持锁聚合 db 快照、解锁后两段式 gap 检查。
    ///
    /// 完备判定（全部本地队列空 ∧ 窃取池空 ∧ active == 0 ⇒ done）与空聚合
    /// 在**同一池锁临界区**内原子发生（M1）：否则『末节点解锁』与『某
    /// worker 醒来置 done』之间存在窗口，进度线程恰在此快照会输出
    /// incumbent=X bound=none gap=inf 的矛盾行。
    ///
    /// 聚合门控：仅 (mipgap 有限正值 || degrade 预算存在) 才做 O(开放节点)
    /// 全量聚合；关闭时 db=None 传入 maybe_gap/maybe_degrade——两者对 None
    /// 各自早退（maybe_gap 的 limit 判空/非正防御、maybe_degrade 的预算判空
    /// 与 AutoDegrade::decide 步骤 3 界缺失防御），语义不变。
    fn complete_node(
        &self,
        w: usize,
        mut children: Vec<PoolNode>,
        sense: Sense,
        mipgap: Option<f64>,
        work_global_notify: bool,
    ) {
        let need_db = matches!(mipgap, Some(g) if g.is_finite() && g > 0.0)
            || self.degrade_budget_secs.is_some();
        let db = {
            let mut pool = self.lock_pool();
            if pool.work_global {
                // 全局 DFS 栈模式：子代压 steal 栈顶（append ⇒ 末元素 = 栈顶
                // = 下一弹出，与本地模式的 push_front 首探索序逐子一致）；
                // push 走 notify_all 快路径（全局栈无亲和可保护，等待者应
                // 立即被唤醒取新栈顶）。
                pool.steal.append(&mut children);
            } else {
                // 子代压自己队首（保同亲 LIFO）：按序 push_front ⇒ owner 队首弹出
                // 顺序与旧全局池『append 后尾弹出』逐子节点一致（children 末元素
                // 先探索）。本地队列 push **不发全局 notify**——等待 worker 经
                // wait_timeout 周期扫 victim 队尾窃取；若此处 notify_one，被唤醒
                // 的睡者会与 owner 争夺刚压入的近亲子代，恰好破坏亲和机制。
                for child in children.drain(..) {
                    pool.locals[w].push_front(child);
                }
            }
            pool.active -= 1;
            pool.inflight[w] = None;
            if pool.all_empty() && pool.active == 0 {
                // 全部容器空且无在处理节点：搜索完备。done 与上方（门控开启
                // 时的）空聚合同临界区，进度线程要么在窗口前快照（容器/在
                // 处理非空 ⇒ bound 必非空聚合），要么看到 done 直接退出——
                // 矛盾行窗口被关闭。
                pool.done = true;
                self.notify_all();
            }
            need_db.then(|| pool_db_aggregate(&pool, sense)).flatten()
        };
        // 全局 DFS 栈模式：栈顶有新子代 ⇒ 唤醒等待者（本地亲和模式不
        // notify——见上方亲和机制注释；此处无亲和可保护）。
        if work_global_notify {
            self.cv_workers.notify_all();
        }
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

    /// 退出路径：仍持有则释放（active−1、inflight=None）；
    /// stop ≠ 0 或（全部容器空且 active == 0）⇒ done + 双 cv notify_all；
    /// 否则窃取池非空 ⇒ cv_workers.notify_one（本地队列不为空不唤醒——
    /// 等待 worker 经 wait_timeout 自调度窃取）。
    fn worker_release(&self, w: usize) {
        let mut pool = self.lock_pool();
        if pool.inflight[w].is_some() {
            pool.active -= 1;
            pool.inflight[w] = None;
        }
        if self.stop_code() != STOP_NONE || (pool.all_empty() && pool.active == 0) {
            pool.done = true;
            self.notify_all();
        } else if !pool.steal.is_empty() {
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

/// 池对偶界聚合（**全量口径**，complete_node（门控后）与截断终界/进度行
/// 共用的唯一实现）：**所有本地队列的全部节点 est ∪ 窃取池全部 est ∪
/// inflight Some 槽**，先滤除非有限值，取方向极值。全量扫描保证不低估
/// 开放节点极值（gap 判定只可能被推迟、不偏提前终止——保守方向）；任何
/// 『只看各容器尾』的口径都漏掉容器内部节点与窃取池，属禁止项。
fn pool_db_aggregate(pool: &PoolState, sense: Sense) -> Option<f64> {
    aggregate_bound(
        sense,
        pool.locals
            .iter()
            .flat_map(|q| q.iter().map(|n| n.est))
            .filter(|e| e.is_finite())
            .chain(
                pool.steal
                    .iter()
                    .map(|n| n.est)
                    .filter(|e| e.is_finite()),
            )
            .chain(
                pool.inflight
                    .iter()
                    .filter_map(|slot| *slot)
                    .filter(|e| e.is_finite()),
            ),
    )
}

/// 终局停止码合成（M3，纯函数）：`limit_hit` 是 LP 迭代上限剪枝的粘滞标记
/// （镜像顺序 limit_hit），**不**截断树搜索；探索完备（stop 仍 NONE）而标记
/// 置位 ⇒ 映射 STOP_NODE（NodeLimit 类 ⇒ 终态诚实 LimitReached，bound 取
/// 剩余池聚合——完备后为空 ⇒ None，与顺序侧栈空聚合同口径）。粘滞 stop
/// （Gap/Degrade/NodeLimit/TimeLimit）优先于 limit_hit。
fn effective_stop_code(stop: u8, limit_hit: bool) -> u8 {
    if stop == STOP_NONE && limit_hit {
        STOP_NODE
    } else {
        stop
    }
}

// ===========================================================================
// 并行分支定界主入口
// ===========================================================================

/// 并行分支定界。前置：`config.threads >= 2`；`node_rule` 非默认值 ⇒
/// `Err(Internal)`。分支规则任意（MostFractional/Pseudocost：Pseudocost 下
/// worker 各持本地 [`Pseudocosts`]，无全局锁）。
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
    if opts.node_rule != NodeRule::DepthFirst {
        return Err(MilpError::Internal(
            "parallel_branch_and_bound 前置失败：仅支持默认节点规则（DepthFirst）".to_string(),
        ));
    }
    let int_vars = problem.integer_vars();
    let sense = problem.sense;
    // 整数目标门控（cutoff 剪枝 + est 弹出剪枝共用；纯问题属性，算一次）。
    let obj_integral = objective_is_integral(problem);
    let mut stats = BbStats::default();
    // 割丢弃回退（根割 IterationLimit 重建 Simplex）前累计的 LP 迭代数：
    // 终局统计以 base + 当前 Simplex 迭代数记账，重建不丢根 LP + 割循环迭代。
    let mut base_iters: u64 = 0;
    let mut simplex = Simplex::new(problem, tol)?;
    // W2 内部墙钟时限（与顺序侧同源合成；0 = None ⇒ 零墙钟读点）。重建
    // Simplex（预割轮 / IterationLimit 回滚）同款重注入。
    let time_deadline = crate::branch_bound::bb_time_deadline(&start, opts);
    simplex.set_time_deadline(time_deadline);
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
            // W2：纯 LP 单节点正是节点边界检查的盲区——根 LP 内部墙钟超限
            // ⇒ 诚实 LimitReached（不冒充 Optimal；纯 LP 无 incumbent 概念）。
            LpStatus::IterationLimit | LpStatus::TimeLimit => SolveStatus::LimitReached,
        });
    } else {
        // 根 solve 非 Optimal ⇒ 不进割循环与启发式，直接终局。
        if let Some(status) = root_terminal_status(sol.status) {
            terminal = Some(status);
        }
        // —— 根预割启发式（镜像顺序侧 branch_bound.rs 根预割轮；LP Optimal
        // 且 root_cuts && heuristics 时额外执行一轮）——
        // 割循环可能把根 LP 顶点从「单断点近整数」推向「多分数中面」（界更
        // 紧但取整/潜水质量下降——max 背包类：纯 LP 贪心顶点仅容量断点一个
        // 分数变量，dive 的断点回退链恰好在其上跨过；割后中面顶点上同一
        // dive 只得弱解）。割前顶点上先跑一轮 round_and_check + dive、与
        // 割后一轮取最优（maybe_update_local 只增不减）；incumbent 只增不减
        // 不影响割的生成。`root_cuts=false` 时不执行 ⇒ 行为逐字节不变。
        if terminal.is_none() && opts.root_cuts && opts.heuristics {
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
            if opts.mitm_windows > 0 && incumbent.is_none() {
                if let Some(inc) = heuristics::subset_sum_mitm(
                    problem,
                    &bounds,
                    &sol.x,
                    opts.mitm_windows,
                    &tol,
                ) {
                    maybe_update_local(sense, &mut incumbent, inc);
                }
            }
            // dive 扰动 Simplex 内部基：割循环入口的防御性 resolve 只恢复
            // Optimal 性，不恢复**同一个**最优基——重退化根 LP（如背包类）
            // 的备择最优基上 GMI 分离结果不同，recorded_cuts 会偏离无预割
            // 路径且与顺序侧冷基割集不一致（实测 11 变量退化背包：修复前
            // par(threads=4) heur on/off 割数 4/5，顺序侧（已重建冷基）5/5）。
            // 此处重建冷 Simplex + 重解根 LP（与顺序侧 branch_bound.rs 预割
            // 轮收尾、本函数 IterationLimit 回退同款模式），使割循环与
            // worker 的 recorded_cuts 重放恒在冷基上；incumbent 已由
            // maybe_update_local 收取，不受重建影响。迭代记账先进
            // base_iters，不丢数。
            base_iters = base_iters.saturating_add(simplex.iterations());
            simplex = Simplex::new(problem, tol)?;
            simplex.set_time_deadline(time_deadline);
            sol = simplex.solve(&bounds);
            if let Some(status) = root_terminal_status(sol.status) {
                terminal = Some(status);
            }
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
            // R8 回滚只对迭代上限触发（W2：TimeLimit 走 else 分支直接诚实
            // 终态——回滚重建+重放是又一次整段根 LP，会把超限点拖过
            // deadline）。
            if st == LpStatus::IterationLimit {
                // 割后 resolve 触发迭代上限：割循环停在非 Optimal 态，不可作
                // 为树搜索根态，也不作整树终态——丢弃全部根割（recorded_cuts
                // 保持空 ⇒ worker 不重放），重建无割 Simplex + solve 恢复
                // Optimal 根态后继续启发式与树搜索。正确性：无割原始松弛的
                // 可行域 ⊇ 含割版本（有效割不切任何整数可行解）⇒ 回退后的
                // 完备搜索覆盖全部整数可行解，最坏等同于关割行为。
                // 记账完整：重建前迭代累计进 base_iters；割已丢，stats.cuts
                // 如实清零。
                base_iters = base_iters.saturating_add(simplex.iterations());
                stats.cuts = 0;
                simplex = Simplex::new(problem, tol)?;
                simplex.set_time_deadline(time_deadline);
                sol = simplex.solve(&bounds);
                if let Some(status) = root_terminal_status(sol.status) {
                    terminal = Some(status);
                }
            } else {
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
            if opts.mitm_windows > 0 && incumbent.is_none() {
                if let Some(inc) = heuristics::subset_sum_mitm(
                    problem,
                    &bounds,
                    &sol.x,
                    opts.mitm_windows,
                    &tol,
                ) {
                    maybe_update_local(sense, &mut incumbent, inc);
                }
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
        stats.lp_iterations = base_iters.saturating_add(simplex.iterations());
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
        // (d)：LP bound 不优于 incumbent 的有效剪枝界（整数目标门控含
        // cutoff，与顺序 branch_or_close 同口径）⇒ 根关闭。
        let mut closed = false;
        if let Some((obj, _)) = incumbent.as_ref() {
            if bound_no_better(
                sense,
                sol.objective,
                *obj,
                obj_integral,
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
                // (f)：分支压 children（est = 根目标）。根阶段用空伪代价
                // 统计做规则化选择（Pseudocost 空表 = f⁻·f⁺ 乘积规则；worker
                // 各自冷启动本地统计——根的观测无法带给 worker，可靠性探测
                // 在根阶段无收益故不做，worker 从深度 1 节点起探测）。
                let root_pc = Pseudocosts::new(problem.vars.len());
                children = branch_children(
                    &int_vars,
                    &bounds,
                    &sol,
                    opts,
                    &root_pc,
                    tol.integrality,
                    1,
                );
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
    // R6-retry 全局 DFS 工作栈：模式写入口池（worker 的子代路由由
    // complete_node 读该标志，notify 参数由 opts 侧同源给出）。
    shared.lock_pool().work_global = opts.global_work_stack;
    shared.write_incumbent(incumbent.clone());
    shared.reserved.store(1, Ordering::SeqCst);
    shared.performed.store(1, Ordering::SeqCst);
    // 自适应 stall 锚点对齐顺序路径（AdaptiveState::new 的
    // last_improvement_node = 0）：根启发式供应的 incumbent 只刷新时刻锚点
    // （last_improve_millis，镜像顺序侧仅更新 last_improvement_secs），**不**
    // 记作节点 1 处的改进——否则 stall(400) 阈值在并行侧晚一个节点触发，
    // 与顺序侧的控制器决策序列分叉。
    if incumbent.is_some() {
        shared.last_improve_millis.store(
            u64::try_from(start.elapsed().as_millis()).unwrap_or(u64::MAX),
            Ordering::SeqCst,
        );
    }
    shared.lock_controller().notify_root_done(1);

    // —— 根子节点入窃取池【持锁】+ 两段式 gap/降级检查（spawn 之前）——
    if !children.is_empty() {
        let db = {
            let mut pool = shared.lock_pool();
            // 根子节点入窃取池（尾弹出 ⇒ 首个 worker 取到的子节点与旧全局
            // 池 append+尾弹出一致）；窃取池 push 走 notify_all 快路径（此刻
            // 尚无 waiter，纯协议统一——spawn 后各 worker 直接从池尾取得）。
            pool.steal.append(&mut children);
            shared.cv_workers.notify_all();
            pool_db_aggregate(&pool, sense)
        };
        shared.maybe_gap(sense, config.mipgap_percent, db);
        shared.maybe_degrade(sense, db);
    }

    // 根 gap/降级触发 ⇒ 直接收尾，不 spawn。
    if matches!(shared.stop_code(), STOP_GAP | STOP_DEGRADE) {
        stats.lp_iterations = base_iters.saturating_add(simplex.iterations());
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
    if shared.lock_pool().all_empty() {
        stats.lp_iterations = base_iters.saturating_add(simplex.iterations());
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
                    obj_integral,
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
        // 迭代汇总与 crate 其余统计同口径用 saturating（溢出饱和而非静默回绕）。
        for handle in handles {
            match handle.join() {
                Ok(iters) => worker_iters = worker_iters.saturating_add(iters),
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
    stats.lp_iterations = base_iters
        .saturating_add(simplex.iterations())
        .saturating_add(worker_iters);
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
    let (stop_out, bound, gap) = match effective_stop_code(
        shared.stop_code(),
        shared.limit_hit.load(Ordering::SeqCst),
    ) {
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
    obj_integral: bool,
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
    // W2：worker 的节点 resolve 同受内部墙钟时限（超限 ⇒ LpStatus::TimeLimit
    // ⇒ 剪枝 + limit_hit，随后节点边界的 wall_exceeded 检查置 STOP_TIME）。
    simplex.set_time_deadline(crate::branch_bound::bb_time_deadline(&start, &opts));
    // 按序重放根阶段记录的割清单（与根线程 add_row 顺序一致）。
    for (terms, rhs) in cuts {
        if let Err(error) = simplex.add_row(terms, RelOp::Le, *rhs) {
            shared.record_error(error);
            shared.set_stop(STOP_NODE);
            return 0;
        }
    }
    // 本地伪代价统计（仅 Pseudocost 规则读写；worker 各持一份，无全局锁——
    // 冷启动由 probe_reliability 在浅层节点供应可靠观测）。
    let mut pc = Pseudocosts::new(problem.vars.len());
    // 窃取扫描随机起点（确定性 SplitMix64，按 worker 错峰种子）。
    let mut steal_rng = (w as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    loop {
        // —— 取节点（单一大锁：本地队首 → 窃取池尾 → 随机起点扫各 victim
        //     队尾偷单个；全空且 active==0 ⇒ done；暂无且 active>0 ⇒
        //     cv_workers.wait_timeout(≤5ms) 周期醒来重扫）——
        let mut node = match shared.acquire_node(w, &mut steal_rng) {
            AcquireOutcome::Node(node) => node,
            AcquireOutcome::Finished => {
                shared.worker_release(w);
                return simplex.iterations();
            }
        };
        // —— est 弹出剪枝（仅整数目标启用，与顺序侧 pop 后跳过同口径）：
        //    父节点 LP 目标 est 已不优于当前 incumbent 的有效剪枝界 ⇒ 该
        //    子树不可能再出严格更优解，直接丢弃。位于 reserve fetch_add
        //    之前：不 reserve、不 performed、不 resolve——杜绝『同
        //    node_limit 下 t4 报 NodeLimit 而 t1 Optimal』的口径分裂。锁序
        //    与既有 Infeasible 剪枝路径一致（complete_node：池锁内归还槽位
        //    → 释放 → 两段式 gap/降级检查）。——
        if obj_integral {
            if let Some(obj) = shared.read_incumbent_obj() {
                if bound_no_better(sense, node.est, obj, true, opts.prune_abs, opts.prune_rel) {
                    shared.complete_node(w, Vec::new(), sense, mipgap, false);
                    continue;
                }
            }
        }
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
                shared.complete_node(w, Vec::new(), sense, mipgap, false);
                continue;
            }
            // (b) 非根 Unbounded（父有界时理论不发生的数值异常）⇒ 诚实 Err。
            LpStatus::Unbounded => {
                shared.record_error(nonroot_unbounded_error());
                shared.set_stop(STOP_NODE);
                shared.worker_release(w);
                return simplex.iterations();
            }
            // (c) IterationLimit ⇒ 剪枝本节点 + 置 limit_hit 粘滞标记（镜像
            //     顺序 PruneLimit ⇒ {limit_hit = true; continue}）：树继续
            //     探索其余节点，终态仍诚实映射 LimitReached——两路径终态与
            //     探索范围一致，不再以粘滞 stop 截断整树。
            // (c') TimeLimit（W2）⇒ 同处置（剪枝 + limit_hit）；超限后
            //     worker 下一节点边界的 wall_exceeded 检查立即置 STOP_TIME，
            //     终态经 TimeLimit 诚实映射 LimitReached。
            LpStatus::IterationLimit | LpStatus::TimeLimit => {
                shared.limit_hit.store(true, Ordering::SeqCst);
                shared.complete_node(w, Vec::new(), sense, mipgap, false);
                continue;
            }
            LpStatus::Optimal => {}
        }
        // 伪代价观测（与顺序侧主循环完全同口径）：子节点 LP 目标相对父 est
        // 的方向化增益 ÷ 该方向分数距离（node.var/tighten_lb/frac_dir 即父
        // 分支决策记账三元组）；非 Optimal 子节点已在上方剪枝路径 continue，
        // 此处必为 Optimal。
        if opts.branch_rule == BranchRule::Pseudocost {
            let gain = match sense {
                Sense::Maximize => node.est - sol.objective,
                Sense::Minimize => sol.objective - node.est,
            };
            if gain.is_finite() && gain > 0.0 && node.frac_dir > tol.integrality {
                pc.record(node.var, node.tighten_lb, gain / node.frac_dir);
            }
        }
        // —— (d)/(e)/adaptive/(f) ——
        let mut sol = sol;
        let mut children: Vec<PoolNode> = Vec::new();
        let mut closed = false;
        // (d)：LP bound 不优于 incumbent 的有效剪枝界（整数目标门控含
        // cutoff，与顺序 branch_or_close 同口径；读锁取当前值）⇒ 节点关闭。
        if let Some(obj) = shared.read_incumbent_obj() {
            if bound_no_better(
                sense,
                sol.objective,
                obj,
                obj_integral,
                opts.prune_abs,
                opts.prune_rel,
            ) {
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
                // P0-2：节点 LP 的可行性/最优性已由 dive 前的 resolve 证实，
                // 回取状态不参与节点判定——Infeasible 属 dive 扰基后的数值
                // 伪影，弃用结果继续分支（对齐顺序侧 branch_or_close 的
                // `let _ = resolve` 语义），**绝不据此 closed=true 丢子树**
                //（曾在此丢失含最优解的子树窗口）；Unbounded 保留诚实 Err。
                let back = simplex.resolve(&node.bounds);
                match dive_back_action(back.status) {
                    DiveBackAction::Refresh => sol = back,
                    DiveBackAction::Discard => {}
                    DiveBackAction::Fail => {
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
            // P1-4 可靠性冷启动探测（与顺序侧同一 pub(crate) 辅助：本地
            // Pseudocosts + 同口径预算/深度门控，无锁；探测只读界并恢复，
            // 不参与节点判定，子树不丢）。
            if opts.branch_rule == BranchRule::Pseudocost {
                probe_reliability(
                    sense,
                    &mut node.bounds,
                    int_vars,
                    &sol.x,
                    sol.objective,
                    tol,
                    &opts,
                    &mut pc,
                    &mut simplex,
                    node.depth as usize,
                );
            }
            children = branch_children(
                int_vars,
                &node.bounds,
                &sol,
                &opts,
                &pc,
                tol.integrality,
                node.depth + 1,
            );
        }
        // —— 节点完成转移 + 两段式 gap 检查 ——
        let notify_stack = opts.global_work_stack && !children.is_empty();
        shared.complete_node(w, children, sense, mipgap, notify_stack);
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

/// 根节点 LP 状态 → 整树终态（`None` = Optimal，继续树搜索）。W2：根 LP
/// 内部墙钟超限（TimeLimit）⇒ 诚实 LimitReached——独立于 IterationLimit，
/// 不触发根割的 R8 回滚（回滚重放是又一次整段根 LP，会把超限点拖过
/// deadline）。
fn root_terminal_status(status: LpStatus) -> Option<SolveStatus> {
    match status {
        LpStatus::Optimal => None,
        LpStatus::Infeasible => Some(SolveStatus::Infeasible),
        LpStatus::Unbounded => Some(SolveStatus::Unbounded),
        LpStatus::IterationLimit | LpStatus::TimeLimit => Some(SolveStatus::LimitReached),
    }
}

/// dive 后记账性回取 resolve 的状态分流（纯函数，语义供单元钉死）。
///
/// 节点 LP 的可行性/最优性由 dive **前**的 resolve 证实；回取只为把扰动的
/// 基拉回本节点界，其状态不参与节点判定：
/// - Optimal ⇒ 用回取解刷新 `sol`（同界重解，目标一致）；
/// - Infeasible ⇒ **弃用**：dive 扰基后对偶单纯形 Phase-1 的数值伪影，
///   不是本节点真实属性——保留 dive 前的 LP 解继续分支，绝不 closed 丢
///   子树（对齐顺序侧 branch_or_close 的 `let _ = resolve` 语义）；
/// - IterationLimit ⇒ 同上弃用（不截断树）；
/// - Unbounded ⇒ 数值异常（父有界时收紧/dive 不可能引出无界）⇒ 诚实 Err。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DiveBackAction {
    /// Optimal：用回取解刷新本节点 LP 解。
    Refresh,
    /// Infeasible / IterationLimit：弃用回取状态，继续按原 LP 解分支。
    Discard,
    /// Unbounded：数值异常，整树诚实终止。
    Fail,
}

fn dive_back_action(status: LpStatus) -> DiveBackAction {
    match status {
        LpStatus::Optimal => DiveBackAction::Refresh,
        LpStatus::Infeasible | LpStatus::IterationLimit | LpStatus::TimeLimit => {
            DiveBackAction::Discard
        }
        LpStatus::Unbounded => DiveBackAction::Fail,
    }
}

/// 非根节点 LP 无界的数值错误（镜像 branch_bound 同名语义）。
fn nonroot_unbounded_error() -> MilpError {
    MilpError::Numerical("node LP unbounded after bounded root".to_string())
}

// ===========================================================================
// 节点判定纯函数（bound_no_better / is_integral 原逐字私有副本已删除，
// 改导入 crate::branch_bound 的 pub(crate) 共享实现——整数目标 cutoff
// （effective_prune_bound 合成）三调用点同口径，杜绝第三份分叉；本模块
// 保留测试向量钉共享函数。select_branch_var 仍为本模块快照版路径专用副本。）
// ===========================================================================

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

/// 分支压栈（(f)）：按 `opts.branch_rule` 选变量（Pseudocost 用调用方本地
/// [`Pseudocosts`]——worker 各持一份，无全局锁；根阶段用空统计做规则化
/// 选择），生成上下行两个子节点（默认上行先产、下行后产——owner 本地
/// 队列按序 push_front 后队首弹出时下行先探索，与顺序栈序一致；P1-5
/// best-estimate 压栈序与顺序侧 branch_or_close 同口径翻转）；est = 本节点
/// resolve 目标（对子树有效的乐观界）；var/tighten_lb/frac_dir 供子节点
/// resolve 后按顺序侧同口径做伪代价被动记账。
fn branch_children(
    int_vars: &[usize],
    bounds: &Bounds,
    sol: &LpSolution,
    opts: &BbOptions,
    pc: &Pseudocosts,
    frac_tol: f64,
    child_depth: u32,
) -> Vec<PoolNode> {
    // 无可分支变量为内部不可达防御（不剪枝即丢弃节点）。
    let Some((var, val)) = (match opts.branch_rule {
        BranchRule::MostFractional => select_branch_var(int_vars, &sol.x),
        BranchRule::Pseudocost => pc.select(int_vars, &sol.x, frac_tol),
    }) else {
        return Vec::new();
    };
    let floor = val.floor();
    let frac = val - floor;
    // P1-5 best-estimate 子代压栈序（与顺序侧 branch_or_close 完全同口径）：
    // 两方向均有伪代价观测（n_obs ≥ 1）且上行估计退化严格更小时，上行后
    // 压 ⇒ 池尾弹出先探索；否则保持现序（下行先探索）。仅改变入池序，
    // 两分支域划分完备性与池互斥弹出『每节点恰处理一次』不变。
    let explore_up_first = pc.n_obs(var) > 0
        && pc.est_gain(var, true, 1.0 - frac) < pc.est_gain(var, false, frac);
    let up = (true, floor + 1.0, 1.0 - frac);
    let down = (false, floor, frac);
    let order = if explore_up_first {
        [down, up]
    } else {
        [up, down]
    };
    let mut children = Vec::new();
    for (tighten_lb, bound, frac_dir) in order {
        let mut child = bounds.clone();
        if tighten_child_snapshot(&mut child, var, tighten_lb, bound) {
            children.push(PoolNode {
                bounds: child,
                est: sol.objective,
                var,
                tighten_lb,
                frac_dir,
                depth: child_depth,
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

    /// 测试用池节点：est/深度可控（var 复用作压力测试的节点 id）。
    fn test_node(est: f64, depth: u32) -> PoolNode {
        PoolNode {
            bounds: Bounds {
                lb: vec![0.0],
                ub: vec![1.0],
            },
            est,
            var: est as usize,
            tighten_lb: false,
            frac_dir: 0.5,
            depth,
        }
    }

    /// 测试用池状态：locals 按 est 序构造（Vec 首元素 = 队首）。
    fn test_pool(
        locals: Vec<Vec<f64>>,
        steal: Vec<f64>,
        inflight: Vec<Option<f64>>,
        active: usize,
    ) -> PoolState {
        PoolState {
            work_global: false,
            locals: locals
                .into_iter()
                .map(|ests| ests.into_iter().map(|e| test_node(e, 1)).collect())
                .collect(),
            steal: steal.into_iter().map(|e| test_node(e, 1)).collect(),
            active,
            done: false,
            inflight,
        }
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

    /// dive 钉死问题（heuristics.rs:241-263）：max 10x1+x2+0.5 s.t.
    /// 20x1+3x2≤43，x∈[0,10] 整数；根 LP 22.0，根 dive 得 21.5。
    /// offset=0.5 有意非整数化（整数目标 cutoff/est 弹出剪枝门控判否，
    /// 保住 gap/降级/node-limit 路径覆盖；见 branch_bound.rs 同名 fixture）。
    fn dive_problem() -> Problem {
        let mut p = Problem::new("dive", Sense::Maximize);
        p.add_var("x1", VarType::Integer).unwrap();
        p.add_var("x2", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_var_bounds(1, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 1.0)], 0.5).unwrap();
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


    /// （L1 回归钉）根预割启发式轮对并行割循环**基透明**：预割 dive 扰动
    /// Simplex 基后，若不在割循环前重建冷基，recorded_cuts（及 worker 重放
    /// 的表格）会在备择最优基上分离——与顺序侧冷基割集不一致且偏离无预割
    /// 路径（实测本算例修复前 par(threads=4) heur on/off 割数 4/5，顺序侧
    /// 修复后 5/5）。修复后 `heuristics` 开关不得改变并行根割施加数；
    /// incumbent 只增不减 ⇒ 两配置均须完备且同最优。算例与顺序侧
    /// branch_bound.rs::precut_round_is_basis_transparent_to_cut_loop 同源
    /// （随机族搜索定位的 11 变量退化背包，dive 恰落在备择最优基上）。
    #[test]
    fn precut_round_is_basis_transparent_to_cut_loop_parallel() {
        let obj = [22.0, 33.0, 29.0, 13.0, 18.0, 33.0, 10.0, 12.0, 15.0, 27.0, 9.0];
        let weights = [21.0, 41.0, 35.0, 25.0, 39.0, 11.0, 42.0, 28.0, 20.0, 36.0, 12.0];
        let cap = 145.07999999999998_f64;
        let mut p = Problem::new("knapsack-degen-par", Sense::Maximize);
        for i in 0..obj.len() {
            p.add_var(&format!("x{i}"), VarType::Binary).unwrap();
        }
        p.set_objective(obj.iter().enumerate().map(|(j, &c)| (j, c)).collect(), 0.0)
            .unwrap();
        p.add_constraint(
            None,
            weights.iter().enumerate().map(|(j, &w)| (j, w)).collect(),
            RelOp::Le,
            cap,
        )
        .unwrap();
        let run = |heuristics: bool| {
            parallel_branch_and_bound(
                &p,
                &BbOptions {
                    root_cuts: true,
                    heuristics,
                    ..BbOptions::default()
                },
                &pconfig(4, None),
                tol(),
            )
            .unwrap_or_else(|e| panic!("heuristics={heuristics} 求解失败：{e}"))
        };
        let with_heur = run(true);
        let without_heur = run(false);
        assert_eq!(
            with_heur.stats.cuts, without_heur.stats.cuts,
            "预割轮不得改变并行根割施加数（recorded_cuts 基透明）"
        );
        assert_eq!(with_heur.stop, None, "heuristics=true 应完备");
        assert_eq!(without_heur.stop, None, "heuristics=false 应完备");
        let a = with_heur.best_objective.expect("应有最优解");
        let b = without_heur.best_objective.expect("应有最优解");
        assert!(
            (a - b).abs() <= 1e-7 * (1.0 + b.abs()),
            "两配置最优应一致: {a} vs {b}"
        );
    }

    /// （P0-2 回归）dive 回取状态分流：节点 LP 可行性已由 dive 前 resolve
    /// 证实——回取 Infeasible 是扰基数值伪影，必须 Discard（弃用结果继续
    /// 分支，**不得 closed 丢子树**；曾在此丢失含最优解的子树窗口）；
    /// IterationLimit 同 Discard；TimeLimit（W2）同 Discard（回取不参与
    /// 节点判定，超限由节点边界检查停树）；Optimal ⇒ Refresh；
    /// Unbounded ⇒ Fail。
    #[test]
    fn dive_back_action_discards_infeasible_artifact() {
        assert_eq!(
            dive_back_action(LpStatus::Infeasible),
            DiveBackAction::Discard,
            "回取 Infeasible = 数值伪影 ⇒ 弃用结果，不 closed、不丢子树"
        );
        assert_eq!(dive_back_action(LpStatus::IterationLimit), DiveBackAction::Discard);
        assert_eq!(dive_back_action(LpStatus::TimeLimit), DiveBackAction::Discard);
        assert_eq!(dive_back_action(LpStatus::Optimal), DiveBackAction::Refresh);
        assert_eq!(dive_back_action(LpStatus::Unbounded), DiveBackAction::Fail);
    }

    /// （P0-2 行为级）dive 后回取 Infeasible 场景子树不丢：hard-knapsack-12
    /// 自适应 dive 路径 + threads=4 重复求解恒收敛真最优 523（若回取伪影
    /// 据此 closed，最优解所在子树被丢即失最优；重复下必现）。
    #[test]
    fn dive_back_infeasible_never_loses_subtree_stress() {
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
        let config = pconfig(4, None);
        let config = ParallelConfig {
            adaptive_heuristics: true,
            ..config
        };
        for rep in 0..50 {
            let out = parallel_branch_and_bound(&p, &BbOptions::default(), &config, tol())
                .unwrap_or_else(|e| panic!("rep {rep} 求解失败：{e}"));
            assert_eq!(out.stop, None, "rep {rep} 应完备");
            assert_close(out.best_objective.unwrap(), 523.0);
        }
    }

    /// (a) 根 gap 触发（spawn 之前完成 ⇒ 与线程数无关地确定）：
    /// mipgap=3.0 ⇒ gap≈2.326% ≤ 3 ⇒ stop=Gap、best=21.5、bound=22.0（证据）、
    /// gap_percent=≈2.33。
    #[test]
    fn gap_triggers_at_root_before_spawn() {
        let p = dive_problem();
        let expected_gap = 100.0 * 0.5 / 21.5;
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
            assert_close(out.best_objective.unwrap(), 21.5);
            assert_close(out.bound.unwrap(), 22.0);
            let g = out.gap_percent.unwrap();
            assert!(
                (g - expected_gap).abs() <= 1e-9,
                "gap = {g}，期望 ≈ {expected_gap}"
            );
            assert_eq!(out.stats.nodes, 1, "根阶段后即触发，无 worker 节点");
        }
    }

    /// (a) mipgap=2.0 ⇒ 根 gap≈2.326% > 2 不触发 ⇒ 树完备：
    /// stop=None、best=21.5、bound=Some(21.5)、gap=Some(0.0)。
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
            assert_eq!(out.stop, None, "threads = {threads}：2.326% > 2% 不触发");
            assert_close(out.best_objective.unwrap(), 21.5);
            assert_close(out.bound.unwrap(), 21.5);
            assert_eq!(out.gap_percent, Some(0.0), "完备且有 incumbent ⇒ 0.0");
        }
    }

    /// (b) node_limit=2：child1 预约号 r=2（处理，performed=2）、child2 r=3
    /// （>2 丢弃 + stop=NodeLimit）——与哪颗子节点被预约无关地确定
    ///（offset=0.5 非整数目标 ⇒ est 弹出剪枝关闭，口径与机制引入前一致）。
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
                21.5,
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

    /// (d) 全量口径聚合（单元级：直接构造池状态）：遍历**所有本地队列的
    /// 全部节点 + 窃取池全部 + inflight Some 槽**取方向极值；inflight
    /// None 槽/非有限值剔除。极值节点位于本地队列**中部**与窃取池时必须
    /// 取到——『各队列尾+inflight』口径会漏掉队中/池中节点（低估开放节点
    /// 极值使 gap 偏小，属不保守方向），本向量钉死防口径回归。
    #[test]
    fn pool_aggregation_full_scan_skips_none_and_nonfinite() {
        let pool = test_pool(
            // 本地队列 0：[front=5.0, middle=21.5, back=3.0]；队列 1 空。
            vec![vec![5.0, 21.5, 3.0], Vec::new()],
            // 窃取池：[7.0]。
            vec![7.0],
            vec![None, Some(2.0), Some(f64::NAN)],
            1,
        );
        // Max 极值 21.5 在队列中部：队尾口径只能取 max(3.0, 7.0, 2.0) = 7.0。
        assert_eq!(pool_db_aggregate(&pool, Sense::Maximize), Some(21.5));
        assert_eq!(pool_db_aggregate(&pool, Sense::Minimize), Some(2.0));
        // 全空 ⇒ None。
        let empty = test_pool(
            vec![Vec::new(), Vec::new()],
            Vec::new(),
            vec![None, None],
            0,
        );
        assert_eq!(pool_db_aggregate(&empty, Sense::Maximize), None);
    }

    /// 取节点顺序契约：① 本地队首（LIFO，最新子代）→ ② 窃取池尾（LIFO）
    /// → ③ victim 队尾偷单个（最老 deferred 分支；owner 的近亲前沿不动）。
    #[test]
    fn take_local_then_steal_pool_then_victim_tail() {
        let mut pool = test_pool(
            vec![vec![11.0, 12.0], vec![31.0, 32.0]],
            vec![21.0, 22.0],
            vec![None, None],
            0,
        );
        let mut rng = 0xDEAD_BEEFu64;
        // w=0：本地队首优先（11 = 最后 push_front 的最新子代）。
        assert_eq!(pool.take_local_or_steal(0, &mut rng).unwrap().est, 11.0);
        assert_eq!(pool.take_local_or_steal(0, &mut rng).unwrap().est, 12.0);
        // 本地空 ⇒ 窃取池尾（22 先于 21，LIFO）。
        assert_eq!(pool.take_local_or_steal(0, &mut rng).unwrap().est, 22.0);
        assert_eq!(pool.take_local_or_steal(0, &mut rng).unwrap().est, 21.0);
        // 窃取池空 ⇒ victim w=1 队尾（32 = 最老 deferred 分支；任何随机起点
        // 扫描下唯一非空 victim 必被命中，结果确定）。
        assert_eq!(pool.take_local_or_steal(0, &mut rng).unwrap().est, 32.0);
        assert_eq!(pool.take_local_or_steal(0, &mut rng).unwrap().est, 31.0);
        // 全容器空 ⇒ None。
        assert!(pool.take_local_or_steal(0, &mut rng).is_none());
        assert!(pool.all_empty());
    }

    /// complete_node：子代压 owner 本地队首且探索序与旧全局池口径一致
    /// （children 末元素先探索）；done 判定要求**全部本地队列 ∧ 窃取池皆空**
    /// （漏任一容器会提前 done 丢树）。
    #[test]
    fn complete_node_pushes_children_to_owner_front_and_done_requires_all_empty() {
        let shared = Shared::new(2, DegradePolicy::default(), None, Instant::now());
        // 种子入窃取池，w=0 取得（active/inflight 锁内记账）。
        shared.lock_pool().steal.push(test_node(9.0, 0));
        let mut rng = 1u64;
        let node = match shared.acquire_node(0, &mut rng) {
            AcquireOutcome::Node(n) => n,
            AcquireOutcome::Finished => panic!("窃取池有节点，不应判定结束"),
        };
        assert_eq!(node.est, 9.0);
        assert_eq!(shared.lock_pool().active, 1);
        // 完成：子代 [up, down] 压 owner 队首 ⇒ owner 下一弹出 = down（与旧
        // 全局池 append[up,down]+尾弹出完全同探索序）。
        shared.complete_node(
            0,
            vec![test_node(1.0, 1), test_node(2.0, 1)],
            Sense::Maximize,
            None,
            false,
        );
        {
            let pool = shared.lock_pool();
            assert_eq!(pool.active, 0);
            assert!(!pool.done, "子代在 owner 本地队列 ⇒ 未完备");
            assert_eq!(pool.locals[0].len(), 2);
        }
        {
            let mut pool = shared.lock_pool();
            let first = pool.take_local_or_steal(0, &mut rng).unwrap();
            let second = pool.take_local_or_steal(0, &mut rng).unwrap();
            assert_eq!(first.est, 2.0, "children 末元素先探索（down-first 保持）");
            assert_eq!(second.est, 1.0);
        }
        // 模拟 w=0 在处理最后一颗节点（叶子）：无子代完成 ⇒ 同临界区 done。
        {
            let mut pool = shared.lock_pool();
            pool.active = 1;
            pool.inflight[0] = Some(2.0);
        }
        shared.complete_node(0, Vec::new(), Sense::Maximize, None, false);
        assert!(shared.lock_pool().done, "全容器空且 active==0 ⇒ done");
    }

    /// 窃取原子性压力（多 worker 窄树耗尽场景）：真实 Shared/PoolState +
    /// 4 线程驱动 acquire_node/complete_node（worker_run 同一代码路径），
    /// 确定性满二叉树（depth < 7 分裂）。守恒断言：每个创建的节点恰被处理
    /// 一次（无重复、无丢失 = 255 = 2^8−1），全部线程经 Finished 退出且
    /// done 置位——等待/窃取协议无永久睡眠、无提前 done、无节点损耗。
    #[test]
    fn steal_atomicity_stress_node_conservation() {
        const THREADS: usize = 4;
        const MAX_DEPTH: u32 = 7;
        let shared = Shared::new(THREADS, DegradePolicy::default(), None, Instant::now());
        // 种子节点（id=0, depth=0）入窃取池。
        shared.lock_pool().steal.push(test_node(0.0, 0));
        let processed: Arc<StdMutex<Vec<u64>>> = Arc::new(StdMutex::new(Vec::new()));
        std::thread::scope(|scope| {
            for w in 0..THREADS {
                let processed = processed.clone();
                let shared = &shared;
                scope.spawn(move || {
                    let mut rng = (w as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
                    loop {
                        match shared.acquire_node(w, &mut rng) {
                            AcquireOutcome::Finished => {
                                shared.worker_release(w);
                                break;
                            }
                            AcquireOutcome::Node(node) => {
                                processed
                                    .lock()
                                    .unwrap_or_else(|p| p.into_inner())
                                    .push(node.var as u64);
                                let id = node.var as u64;
                                let children = if node.depth < MAX_DEPTH {
                                    vec![
                                        test_node((id * 2 + 1) as f64, node.depth + 1),
                                        test_node((id * 2 + 2) as f64, node.depth + 1),
                                    ]
                                } else {
                                    Vec::new()
                                };
                                shared.complete_node(w, children, Sense::Maximize, None, false);
                            }
                        }
                    }
                });
            }
        });
        assert!(shared.lock_pool().done, "窄树耗尽后必须置 done（永久睡眠/提前退出回归）");
        let mut ids = processed
            .lock()
            .unwrap_or_else(|p| p.into_inner())
            .clone();
        ids.sort_unstable();
        let expected: Vec<u64> = (0..(1u64 << (MAX_DEPTH + 1)) - 1).collect();
        assert_eq!(ids.len(), expected.len(), "处理总数守恒（无丢失）");
        assert_eq!(ids, expected, "每节点恰处理一次（无重复、无丢失）");
    }

    /// 等待协议（防永久睡眠）：本地空 + 窃取池空且仍有在处理节点时，
    /// 等待 worker 经 wait_timeout(≤5ms) 周期醒来，偷到 victim 后续
    /// push_front 的本地节点【victim 不发任何 notify】——窃取目标是
    /// victim 队尾（最老），owner 的最新子代留在队首不动。
    #[test]
    fn wait_timeout_steals_victim_local_push_without_notify() {
        let shared = Shared::new(2, DegradePolicy::default(), None, Instant::now());
        // w=0 在处理（active=1）⇒ w=1 取不到节点时进入 wait_timeout。
        {
            let mut pool = shared.lock_pool();
            pool.active = 1;
            pool.inflight[0] = Some(1.0);
        }
        let (est, latency) = std::thread::scope(|scope| {
            let shared = &shared;
            let handle = scope.spawn(move || {
                let mut rng = 42u64;
                match shared.acquire_node(1, &mut rng) {
                    AcquireOutcome::Node(n) => n.est,
                    AcquireOutcome::Finished => panic!("active=1 时不应判定结束"),
                }
            });
            // 让 waiter 进入等待（覆盖 OS 调度余量）后，victim w=0 压两个
            // 子代到队首【不发全局 notify——协议规定本地 push 静默】。
            std::thread::sleep(Duration::from_millis(50));
            let push_t0 = Instant::now();
            {
                let mut pool = shared.lock_pool();
                pool.locals[0].push_front(test_node(5.0, 2)); // 较老（队尾）
                pool.locals[0].push_front(test_node(6.0, 2)); // 最新（队首）
            }
            let est = handle.join().unwrap();
            (est, push_t0.elapsed())
        });
        assert_eq!(est, 5.0, "窃取者取 victim 队尾（最老 deferred 分支）");
        assert!(
            latency < Duration::from_secs(2),
            "wait_timeout 必须周期醒来窃取（实测 {latency:?}）——永久睡眠回归"
        );
    }

    /// 等待协议快路径：窃取池 push 走 notify_all，等待 worker 立即被唤醒
    /// 取得节点（不依赖 wait_timeout 超时兜底）。
    #[test]
    fn steal_pool_push_notifies_waiters_immediately() {
        let shared = Shared::new(2, DegradePolicy::default(), None, Instant::now());
        {
            let mut pool = shared.lock_pool();
            pool.active = 1;
            pool.inflight[0] = Some(1.0);
        }
        let (est, latency) = std::thread::scope(|scope| {
            let shared = &shared;
            let handle = scope.spawn(move || {
                let mut rng = 7u64;
                match shared.acquire_node(1, &mut rng) {
                    AcquireOutcome::Node(n) => n.est,
                    AcquireOutcome::Finished => panic!("active=1 时不应判定结束"),
                }
            });
            std::thread::sleep(Duration::from_millis(50));
            let push_t0 = Instant::now();
            shared.lock_pool().steal.push(test_node(77.0, 1));
            shared.cv_workers.notify_all();
            let est = handle.join().unwrap();
            (est, push_t0.elapsed())
        });
        assert_eq!(est, 77.0);
        assert!(
            latency < Duration::from_secs(2),
            "notify_all 快路径必须唤醒等待 worker（实测 {latency:?}）"
        );
    }

    /// t4 与 t1 同界收敛（bench 真实算例）：knapsack-200/500（.opt 为 DP
    /// 独立真值）与 mknap-3x100，threads=4 并行完备搜索（重复 3 次覆盖
    /// 窃取调度的节点归属漂移）与 t1 顺序路径最优值逐次一致。mknap-5x300
    /// （90s 量级）归 bench 集成回归，不入单测。
    #[test]
    fn t4_matches_t1_objectives_on_bench_knapsacks() {
        const CASES: &[(&str, f64)] = &[
            ("knapsack-200", 3409.0),
            ("knapsack-500", 8975.0),
            ("mknap-3x100", 2331.0),
        ];
        for &(name, expected) in CASES {
            let path = format!(
                "{}/../bench/instances/{name}.milp",
                env!("CARGO_MANIFEST_DIR")
            );
            let problem = load_instance(Path::new(&path))
                .unwrap_or_else(|e| panic!("算例 {name} 加载失败：{e}"));
            let seq = branch_and_bound(&problem, &BbOptions::default(), tol())
                .unwrap_or_else(|e| panic!("算例 {name} 顺序求解失败：{e}"));
            assert_eq!(seq.status, SolveStatus::Optimal, "算例 {name}");
            assert_close(seq.best_objective.unwrap(), expected);
            for rep in 0..3 {
                let out = parallel_branch_and_bound(
                    &problem,
                    &BbOptions::default(),
                    &pconfig(4, None),
                    tol(),
                )
                .unwrap_or_else(|e| panic!("算例 {name} rep{rep} 求解失败：{e}"));
                assert_eq!(out.stop, None, "算例 {name} rep{rep} 应完备");
                assert_close(
                    out.best_objective.unwrap(),
                    seq.best_objective.unwrap(),
                );
            }
        }
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
    /// spawn）时恰好一行，字段口径 nodes=1 / incumbent=21.5 / bound=22.0 /
    /// gap=2.33%（offset=0.5 非整数目标 ⇒ cutoff 关闭，根不闭合、gap 路径
    /// 照常触发）。
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
            lines[0].contains("nodes=1 incumbent=21.500000 bound=22.000000 gap=2.33%"),
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

    // —— 共享纯函数的测试向量（原逐字副本已删除，向量改钉
    //    crate::branch_bound 的 pub(crate) 共享实现）——

    /// (d) 剪枝比较按原始 sense 方向（Max/Min 镜像；integral=false 退化
    /// 口径——与 branch_bound.rs 同名向量一致，钉共享函数不分叉）。
    #[test]
    fn bound_prune_respects_sense() {
        let (pa, pr) = (1e-9, 1e-9);
        assert!(bound_no_better(Sense::Maximize, 10.0, 10.0, false, pa, pr));
        assert!(!bound_no_better(Sense::Maximize, 10.5, 10.0, false, pa, pr));
        assert!(bound_no_better(Sense::Minimize, 2.0, 2.0, false, pa, pr));
        assert!(!bound_no_better(Sense::Minimize, 1.5, 2.0, false, pa, pr));
        // 整数目标 cutoff 口径（共享 effective_prune_bound 合成）：
        // est ∈ (inc+slack, inc+1−ε] ⇒ 剪枝；est = inc+1 恰在网格 ⇒ 保留。
        assert!(bound_no_better(Sense::Maximize, 21.5, 21.0, true, pa, pr));
        assert!(!bound_no_better(Sense::Maximize, 22.0, 21.0, true, pa, pr));
        assert!(bound_no_better(Sense::Minimize, 29.5, 30.0, true, pa, pr));
        assert!(!bound_no_better(Sense::Minimize, 29.0, 30.0, true, pa, pr));
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

    /// 并行 t=2 与 t=1 同 node_limit 同终态/objective（est 弹出剪枝在
    /// reserve 之前 ⇒ 不 reserve、不 performed：陈旧前沿不烧节点预算——
    /// 杜绝『同 node_limit 下 t4 报 NodeLimit 而 t1 Optimal』的口径分裂）。
    ///
    /// 载体与 branch_bound::chain_knapsack 同构：1 大件（w=7,v=10）+ 8 个
    /// 相同小件（w=2,v=1），容量 12；整数目标（cutoff 启用），穷举最优 12。
    /// 预算取『非整数孪生（= 旧口径）的完备需求』——该值对两引擎都是充裕
    /// 上界（无剪枝全树 ≤ 此值；cutoff 生效后 ≪ 此值）；孪生自身在同值 −1
    /// 的预算下必截断（LimitReached），证明预算确实卡在旧口径的需求上。
    #[test]
    fn parallel_t2_same_node_limit_same_terminal_as_t1() {
        let chain = |offset: f64| {
            let n = 8usize;
            let mut p = Problem::new("chain-knap", Sense::Maximize);
            p.add_var("b", VarType::Binary).unwrap();
            for i in 0..n {
                p.add_var(&format!("s{i}"), VarType::Binary).unwrap();
            }
            let mut terms = vec![(0usize, 10.0)];
            terms.extend((0..n).map(|i| (i + 1, 1.0)));
            p.set_objective(terms, offset).unwrap();
            let mut weights = vec![(0usize, 7.0)];
            weights.extend((0..n).map(|i| (i + 1, 2.0)));
            p.add_constraint(None, weights, RelOp::Le, 12.0).unwrap();
            p
        };
        let plain = BbOptions {
            root_cuts: false,
            heuristics: false,
            ..BbOptions::default()
        };
        // 旧口径（offset=0.5 孪生）的完备需求 = 预算上界（顺序侧确定性）。
        let need_old = branch_and_bound(&chain(0.5), &plain, tol())
            .unwrap()
            .stats
            .nodes;
        let limit = need_old;
        // 预算确实卡在旧口径需求上：−1 即截断。
        let tight = BbOptions {
            node_limit: limit - 1,
            ..plain
        };
        let twin = branch_and_bound(&chain(0.5), &tight, tol()).unwrap();
        assert_eq!(
            twin.status,
            SolveStatus::LimitReached,
            "孪生在同预算 −1 下必截断（实际 nodes={}）",
            twin.stats.nodes
        );
        // 同 node_limit：t=1（受控顺序）与 t=2（并行）同终态 Optimal、同目标。
        let opts = BbOptions {
            node_limit: limit,
            ..plain
        };
        let t1 = crate::branch_bound::branch_and_bound_controlled(
            &chain(0.0),
            &opts,
            &crate::branch_bound::SearchControl::default(),
            tol(),
            None,
        )
        .unwrap();
        assert_eq!(t1.status, SolveStatus::Optimal);
        assert_close(t1.best_objective.unwrap(), 12.0);
        let out = parallel_branch_and_bound(&chain(0.0), &opts, &pconfig(2, None), tol()).unwrap();
        assert_eq!(out.root_terminal, None);
        assert_eq!(
            out.stop,
            None,
            "同预算下 t2 不得 NodeLimit（est-skip 不 reserve；实际 {:?}，nodes={}）",
            out.stop,
            out.stats.nodes
        );
        assert_close(out.best_objective.unwrap(), t1.best_objective.unwrap());
        assert_close(out.best_objective.unwrap(), 12.0);
        assert!(
            out.stats.nodes <= limit,
            "performed 含根 =1；实际 {}",
            out.stats.nodes
        );
    }

    /// 策略级自动降级（并行路径；激进参数 start=0/plateau=0 ⇒ 首个节点
    /// 完成即触发）：stop = AutoDegrade、携带 incumbent 21.5（offset=0.5
    /// 非整数目标 ⇒ cutoff 关闭，树搜索路径与机制引入前一致）、gap 证据
    /// 非负有限（同 STOP_GAP 口径，由 CAS 胜者写入）。
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
            (obj - 21.5).abs() <= 1e-6,
            "incumbent = {obj}，期望根 dive 供应的 21.5"
        );
        let gap = out.gap_percent.expect("降级 ⇒ gap 证据必为 Some");
        assert!(gap.is_finite() && gap >= 0.0, "gap = {gap}");
    }

    /// （M2）并行自适应 stall 锚点初始化与顺序路径对齐：根启发式供应的
    /// incumbent 不把 last_improve_node 记作节点 1——顺序侧
    /// AdaptiveState::new 恒 0、根阶段只刷新时刻锚点；只有树内真实改进
    ///（try_update_incumbent）才推进节点锚点。否则 stall(400) 阈值在并行
    /// 侧晚一个节点触发，两路径的控制器决策序列分叉。
    #[test]
    fn parallel_stall_anchor_init_matches_sequential() {
        let shared = Shared::new(2, DegradePolicy::default(), None, Instant::now());
        assert_eq!(shared.last_improve_node.load(Ordering::SeqCst), 0);
        // 模拟根阶段末尾的共享状态写入（incumbent 存在、节点锚点不动）。
        shared.write_incumbent(Some((21.5, vec![2.0, 1.0])));
        assert_eq!(
            shared.last_improve_node.load(Ordering::SeqCst),
            0,
            "根 incumbent 不计为节点 1 的改进（与顺序侧 last_improvement_node = 0 对齐）"
        );
        // 树内真实改进按 node_no 推进锚点。
        let inc = Incumbent {
            objective: 22.0,
            values: vec![2.0, 2.0],
        };
        assert!(shared.try_update_incumbent(Sense::Maximize, inc, 7, 0.5));
        assert_eq!(shared.last_improve_node.load(Ordering::SeqCst), 7);
    }

    /// （M3）终局停止码合成：LP 迭代上限剪枝标记（limit_hit，镜像顺序
    /// PruneLimit 的 limit_hit）在树探索完备（stop 仍 NONE）时映射
    /// STOP_NODE ⇒ 终态诚实 LimitReached，不冒充完备；粘滞 stop 优先。
    #[test]
    fn effective_stop_code_limit_hit_composition() {
        assert_eq!(effective_stop_code(STOP_NONE, false), STOP_NONE);
        assert_eq!(
            effective_stop_code(STOP_NONE, true),
            STOP_NODE,
            "迭代上限剪枝后探索完备 ⇒ NodeLimit 类诚实截断"
        );
        assert_eq!(effective_stop_code(STOP_GAP, true), STOP_GAP, "粘滞 stop 优先");
        assert_eq!(effective_stop_code(STOP_DEGRADE, true), STOP_DEGRADE);
        assert_eq!(effective_stop_code(STOP_TIME, true), STOP_TIME);
        assert_eq!(effective_stop_code(STOP_NODE, false), STOP_NODE);
    }

    /// （M1 压测）并行进度行永不出 bound=none：threads=4 + 极小 log 周期
    ///（协调线程近似忙转快照）+ hard-knapsack-12 循环 100 次，断言 0 条
    /// bound=none 行。修复前『末节点完成解锁池锁』与『置 done』之间的窗口
    /// 会被空聚合快照出 incumbent=X bound=none gap=inf 矛盾行（原评审
    /// 800 次循环 3.1M 行中复现 1804 行；顺序路径同参数 0 例）；修复后
    /// done 置位与空聚合在同一池锁临界区内原子发生，窗口关闭。
    #[test]
    fn parallel_progress_lines_never_report_bound_none() {
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
        for rep in 0..100 {
            let lines: Arc<StdMutex<Vec<String>>> = Arc::new(StdMutex::new(Vec::new()));
            let sink = lines.clone();
            let config = ParallelConfig {
                threads: 4,
                mipgap_percent: None,
                adaptive_heuristics: true,
                log_period_secs: Some(1e-9),
                degrade: DegradePolicy::default(),
                progress: Some(Arc::new(move |info: &ProgressInfo| {
                    sink.lock()
                        .unwrap_or_else(|poisoned| poisoned.into_inner())
                        .push(crate::search_ctrl::format_progress_line(info));
                })),
            };
            let out = parallel_branch_and_bound(&p, &BbOptions::default(), &config, tol())
                .unwrap_or_else(|e| panic!("rep {rep} 求解失败：{e}"));
            assert_eq!(out.stop, None, "rep {rep} 应完备");
            assert_close(out.best_objective.unwrap(), 523.0);
            let lines = lines
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            assert!(!lines.is_empty(), "rep {rep}：极小周期下必有进度行");
            for l in lines.iter() {
                assert!(!l.contains("bound=none"), "rep {rep} 矛盾进度行：{l}");
            }
        }
    }
}
