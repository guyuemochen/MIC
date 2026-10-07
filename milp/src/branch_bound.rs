//! 深度优先分支定界（WP4 契约，本文件为 WP4 实现）。
//!
//! 架构：单一 [`crate::simplex::Simplex`] 实例 + Bounds 界栈（trail）+ 深度优先，
//! 节点零拷贝；下行分支先探索；回溯严格按 trail 逆序恢复旧界（恢复后由
//! resolve 的松弛失贴 snap 处理非基内点）。
//!
//! 节点处理（resolve 返回后逐条判定，顺序固定）：
//! (a) Infeasible ⇒ 剪枝回溯；
//! (b) Unbounded：根节点 → 整树终止返回 [`SolveStatus::Unbounded`]；非根节点
//!     （父有界时理论不发生，数值异常时可能出现）⇒ 立即整树终止返回
//!     `Err(MilpError::Numerical)`，不得静默继续；
//! (c) IterationLimit ⇒ 剪枝并置 limit_hit 标记；
//! (d) LP bound 不优于 incumbent（Max: bound ≤ inc + prune_abs + prune_rel·|inc|；
//!     Min 镜像）⇒ 剪枝；
//! (e) 整数性满足（全部 int_var |x−round| ≤ tol.integrality）⇒ round_and_check
//!     复核：通过且严格更优 ⇒ 更新 incumbent 并回溯；复核不通过（伪整数解）⇒
//!     不得剪枝、不得记 incumbent，必须按 MostFractional 继续分支；
//! (f) 空界（tighten 返回 false）⇒ 剪枝。
//!
//! 根节点顺序：int_vars 空 ⇒ 单次 LP 直接返回；根 solve 非 Optimal ⇒ 不进割循环
//! 与启发式，直接按 (a)/(b)/(c) 终局；根 Optimal ⇒ （可选）cut_loop ⇒ （可选）
//! round_and_check + dive ⇒ DFS 主循环。
//!
//! 终态判定（与上限诚实性整合）：主循环因栈空自然退出或因 node_limit/
//! time_limit 提前退出后统一判定——limit_hit 为真或上限触发 ⇒ LimitReached
//! （附 incumbent，若有）；否则有 incumbent ⇒ Optimal、无 ⇒ Infeasible。
//!
//! # 实现说明（WP4）
//! - 节点记账：根 LP 求解与每个树节点的 resolve 各计 1 节点；割循环后取回
//!   根 LP 解、dive 后恢复根热状态属记账性 resolve，不计节点（其迭代数仍
//!   全部计入 `stats.lp_iterations`）。
//! - 树结构：显式 `stack`（待应用的分支决策 [`Pending`]）+ `trail`（已应用
//!   界修改的回放记录）。每个 `Pending` 记录应用前 trail 应回退到的长度
//!   `unwind_to`，保证兄弟节点在父节点界状态上施加自己的分支界；回退严格
//!   按 trail 逆序恢复旧值。
//! - 进展保证：子节点分支界相对当前界严格收缩且新区间非空时才压栈
//!   （[`tighten_child`]）。非分数 x̄ 的两个分支界中至多一侧无收缩，伪整数
//!   回退分支因此在最坏情形下也单调收缩，节点上限兜底。
//! - 上限检查位于栈弹出之后、resolve 之前：栈空 ⇒ 搜索自然耗尽（先于上限
//!   判定，完整搜索不因恰好用满预算而谎报 LimitReached）；弹出待处理决策后
//!   触发上限 ⇒ 搜索确实未完成 ⇒ 诚实置 limit_hit。
//! - 可选优化（默认关闭）：纯整数系数目标可用 cutoff=⌊inc⌋−1 加强剪枝，
//!   因浮点边界错杀风险本实现不启用。
//!
//! # 受控扩展（默认全关 ⇒ 与原实现逐字节一致）
//! [`branch_and_bound_controlled`] 经 [`SearchControl`]（threads / mipgap /
//! adaptive_heuristics / log_period_secs）与 [`SearchHooks`]（进度回调）注入
//! 新能力，既有 [`branch_and_bound`] 签名与行为完全不变：
//! - 顺序路径（threads ≤ 1）全部改动门控：`Pending.est`（压栈时父节点 LP
//!   目标）为纯数据不参与决策；首行进度在根阶段结束的汇合点发射（含纯 LP
//!   捷径与根终局，先发射再返回）；主循环每轮顶部（pop 之前）仅在栈非空时
//!   ①发射周期行（bound = est 聚合，gap 经 [`crate::search_ctrl::gap_percent`]
//!   clamp ≥ 0）与 ②做 mipgap 判定（达标即粘滞跳出并记触发时刻证据）；
//!   自适应启发式在 branch_or_close 的分支路径上（(d) 未剪枝且非整数记录、
//!   选分支变量之前）由控制器决定是否补跑 dive。
//! - threads ≥ 2 组装 [`crate::parallel_bb::ParallelConfig`] 委托
//!   [`crate::parallel_bb::parallel_branch_and_bound`]，并按统一口径映射
//!   ParallelOutcome
//!   （根终局直返；NodeLimit/TimeLimit/Gap ⇒ LimitReached，Gap 附
//!   `gap_terminated = true` 与证据 bound；无 stop ⇒ incumbent?Optimal:Infeasible）。

use std::time::Instant;

use crate::cuts::{self, CutOptions};
use crate::error::{MilpError, Result};
use crate::heuristics::{self, Incumbent};
use crate::model::{Problem, Sense, SolveStatus, Tolerances};
use crate::search_ctrl::{
    aggregate_bound, gap_percent, AutoDegrade, DegradePolicy, DegradeSnapshot, HeuristicConfig,
    HeuristicController, HeuristicDecision, ProgressInfo, SearchSnapshot,
};
use crate::simplex::{Bounds, LpSolution, LpStatus, Simplex};

/// 分支变量选择规则（首版仅此一档）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BranchRule {
    /// 最分数变量（max min(f, 1−f)，并列取最小下标）。
    MostFractional,
}

/// 节点选择规则（首版仅 DFS；best-bound 需按节点存基或冷重启，代价大）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeRule {
    /// 深度优先 + 界栈热启动 + 回溯恢复旧界。
    DepthFirst,
}

/// 分支定界参数。
#[derive(Debug, Clone, Copy)]
pub struct BbOptions {
    /// 分支变量规则。默认 MostFractional。
    pub branch_rule: BranchRule,
    /// 节点选择规则。默认 DepthFirst。
    pub node_rule: NodeRule,
    /// 节点上限（每次 resolve 计 1 node）。默认 100_000。
    pub node_limit: u64,
    /// 时间上限毫秒；0 = 不限。默认 0。
    pub time_limit_ms: u64,
    /// 是否启用根节点割循环。默认 true。
    pub root_cuts: bool,
    /// 是否启用启发式。默认 true。
    pub heuristics: bool,
    /// dive 修复变量数上限。默认 8。
    pub dive_max_fixes: usize,
    /// 剪枝绝对容差。默认 1e-9。
    pub prune_abs: f64,
    /// 剪枝相对容差。默认 1e-9。
    pub prune_rel: f64,
}

impl Default for BbOptions {
    /// 默认值：MostFractional / DepthFirst / 100_000 / 0 / true / true / 8 / 1e-9 / 1e-9。
    fn default() -> Self {
        BbOptions {
            branch_rule: BranchRule::MostFractional,
            node_rule: NodeRule::DepthFirst,
            node_limit: 100_000,
            time_limit_ms: 0,
            root_cuts: true,
            heuristics: true,
            dive_max_fixes: 8,
            prune_abs: 1e-9,
            prune_rel: 1e-9,
        }
    }
}

/// 搜索控制旋钮（全部默认 = 现行为；纯数据，Debug+Clone+Copy 安全）。
///
/// 经 [`branch_and_bound_controlled`] 注入；[`SearchControl::default`] 与
/// `hooks = None` 的组合下顺序路径与既有 [`branch_and_bound`] 逐字节一致
/// （所有新能力均门控，默认全关）。
#[derive(Debug, Clone, Copy)]
pub struct SearchControl {
    /// 线程数。默认 1；1 = 纯顺序路径（与现版本一致），≥2 组装
    /// [`crate::parallel_bb::ParallelConfig`] 委托并行引擎。
    pub threads: usize,
    /// MIP gap 百分比阈值。默认 None = 关闭；Some(0.0) 等价关闭；Some(p) 且
    /// p > 0 时剩余对偶界与 incumbent 的 gap（一律 clamp ≥ 0，见
    /// [`crate::search_ctrl::gap_percent`]）≤ p% 即提前终止，终态诚实标注
    /// LimitReached 且 `gap_terminated = true`。
    pub mipgap_percent: Option<f64>,
    /// 是否启用自适应启发式（树节点分支路径上按控制器决策补跑 dive）。
    /// 默认 false（现行为：仅根节点启发式）。
    pub adaptive_heuristics: bool,
    /// 周期进度行的间隔秒数。默认 None = 不输出周期行（首行仍随钩子发出）。
    pub log_period_secs: Option<f64>,
    /// 策略级自动降级参数（「继续证最优」vs「提前返回最好解」）。仅在
    /// `adaptive_heuristics = true` 且给了时间预算（`BbOptions::time_limit_ms
    /// != 0`）时评估；触发时与 mipgap 同径提前终止（LimitReached +
    /// gap_terminated），另附 `auto_degraded = true` 标记。默认
    /// [`DegradePolicy::default`]。
    pub degrade: DegradePolicy,
}

impl Default for SearchControl {
    /// 默认值：threads=1 / mipgap_percent=None / adaptive_heuristics=false /
    /// log_period_secs=None / degrade=DegradePolicy::default()（= 现版本
    /// 行为，零变化）。
    fn default() -> Self {
        SearchControl {
            threads: 1,
            mipgap_percent: None,
            adaptive_heuristics: false,
            log_period_secs: None,
            degrade: DegradePolicy::default(),
        }
    }
}

/// 进度回调注入（仅 derive(Clone)：回调闭包不可 Copy，derive(Copy) 触发 E0277）。
#[derive(Clone)]
pub struct SearchHooks {
    /// 进度回调（`Arc<dyn Fn(&ProgressInfo) + Send + Sync>`）。
    pub progress: crate::search_ctrl::ProgressCallback,
}

/// 分支定界统计。
#[derive(Debug, Clone, Copy, Default)]
pub struct BbStats {
    /// 处理节点数（每次 resolve 计 1）。
    pub nodes: u64,
    /// 累计单纯形迭代数。
    pub lp_iterations: u64,
    /// 施加割数。
    pub cuts: u64,
    /// 累计耗时毫秒。
    pub time_ms: u64,
}

/// 分支定界结果。
#[derive(Debug, Clone)]
pub struct BbResult {
    /// 求解终态。gap 提前终止时 = LimitReached（诚实：未证明最优）。
    pub status: SolveStatus,
    /// 最优目标（原始 sense + obj_offset 口径；LimitReached 时为已有 incumbent）。
    pub best_objective: Option<f64>,
    /// 最优解（与入参 problem 同维同序，未经 postsolve；整数变量已取整）。
    pub best_values: Option<Vec<f64>>,
    /// 统计信息。
    pub stats: BbStats,
    /// 求解结束时的对偶界（原始 sense + obj_offset 口径）。gap 终止 =
    /// 触发时刻证据 db（粘滞快照，不依赖后续栈态）；截断（上限触发）= 剩余
    /// 待处理子树的 est 聚合（栈空时可 None）；完备且有 incumbent =
    /// best_objective；Infeasible/Unbounded = None。
    pub bound: Option<f64>,
    /// 是否因达到 mipgap 阈值提前终止（true 时 status 必为 LimitReached，
    /// bound 为触发时刻证据）。
    pub gap_terminated: bool,
    /// 是否因策略级自动降级提前终止（auto 模式 + 时间预算 + 双停滞触发）。
    /// true ⇒ `gap_terminated` 必为 true、status 必为 LimitReached——与
    /// mipgap 同径的诚实提前终止，此字段区分触发原因。
    pub auto_degraded: bool,
}

/// 深度优先分支定界主入口（签名与行为完全不变；默认受控路径的薄委托）。
pub fn branch_and_bound(problem: &Problem, opts: &BbOptions, tol: Tolerances) -> Result<BbResult> {
    branch_and_bound_controlled(problem, opts, &SearchControl::default(), tol, None)
}

/// 受控入口：threads ≤ 1 走纯顺序路径；threads ≥ 2 组装
/// [`crate::parallel_bb::ParallelConfig`] 委托
/// [`crate::parallel_bb::parallel_branch_and_bound`] 并
/// 按统一口径映射 ParallelOutcome：
/// - 根阶段终局 ⇒ 直返；
/// - stop ∈ {NodeLimit, TimeLimit}（或顺序侧 IterationLimit 类截断）⇒
///   LimitReached（附 incumbent，若有）；
/// - stop = Gap ⇒ LimitReached 且 `gap_terminated = true`（bound/gap_percent
///   用触发时刻证据口径）；
/// - stop = None ⇒ 有 incumbent 则 Optimal、无则 Infeasible（粘滞 stop 优先
///   于完备判定）。
///
/// 默认（[`SearchControl::default`] + `hooks = None`）与 [`branch_and_bound`]
/// 逐字节一致。
pub fn branch_and_bound_controlled(
    problem: &Problem,
    opts: &BbOptions,
    ctrl: &SearchControl,
    tol: Tolerances,
    hooks: Option<&SearchHooks>,
) -> Result<BbResult> {
    if ctrl.threads >= 2 {
        return parallel_entry(problem, opts, ctrl, tol, hooks);
    }
    sequential_bb(problem, opts, ctrl, tol, hooks)
}

/// 顺序分支定界（既有 [`branch_and_bound`] 的实现体 + 全部门控的受控扩展；
/// 门控默认全关 ⇒ 行为与原实现逐字节一致）。
fn sequential_bb(
    problem: &Problem,
    opts: &BbOptions,
    ctrl: &SearchControl,
    tol: Tolerances,
    hooks: Option<&SearchHooks>,
) -> Result<BbResult> {
    let start = Instant::now();
    let int_vars = problem.integer_vars();
    let sense = problem.sense;
    let mut stats = BbStats::default();
    let mut simplex = Simplex::new(problem, tol)?;
    let mut bounds = Bounds::from_problem(problem);
    // incumbent = (目标, 取值)；目标一律按取整后的值评估（WP3 复核口径）。
    let mut incumbent: Option<(f64, Vec<f64>)> = None;
    // —— 受控扩展状态（默认全关 ⇒ 零行为变化）——
    // 上行进度行的发射时刻（秒）；仅钩子启用时被写入。
    let mut last_log: Option<f64> = None;
    // gap 达标提前终止（粘滞）与触发时刻证据 db（与并行证据口径一致，不依赖
    // 后续栈态）。
    let mut gap_hit = false;
    let mut gap_evidence: Option<f64> = None;
    // 自适应启发式状态（adaptive_heuristics=false ⇒ 恒 None，零开销零行为）。
    let mut adaptive = ctrl
        .adaptive_heuristics
        .then(|| AdaptiveState::new(HeuristicConfig::default(), ctrl.degrade));

    // —— 根 LP（计 1 节点）——
    let mut sol = simplex.solve(&bounds);
    stats.nodes = 1;

    // 纯 LP 捷径：无整数变量 ⇒ 单次 LP 直接按其状态返回。
    if int_vars.is_empty() {
        let status = match sol.status {
            LpStatus::Optimal => SolveStatus::Optimal,
            LpStatus::Infeasible => SolveStatus::Infeasible,
            LpStatus::Unbounded => SolveStatus::Unbounded,
            LpStatus::IterationLimit => SolveStatus::LimitReached,
        };
        let lp_obj = sol.objective;
        let incumbent = if sol.status == LpStatus::Optimal {
            Some((sol.objective, sol.x))
        } else {
            None
        };
        // 完备（Optimal）⇒ bound = 最优值；其余终局无对偶界可报。
        let bound = if status == SolveStatus::Optimal {
            Some(lp_obj)
        } else {
            None
        };
        // 首行进度：根阶段终局汇合点之一（先发射再返回）。
        if let Some(h) = hooks {
            emit_progress(
                h,
                &start,
                stats.nodes,
                incumbent.as_ref().map(|(o, _)| *o),
                bound,
                sense,
            );
        }
        return Ok(finish(
            status,
            incumbent,
            &mut stats,
            &simplex,
            start,
            bound,
            false,
            false,
        ));
    }

    // 根 solve 非 Optimal ⇒ 不进割循环与启发式（gomory_cuts 的前提是 Optimal
    // 表格），直接按 (a)/(b)/(c) 终局。
    if let Some(status) = root_terminal(classify_node(sol.status, true))? {
        // 首行进度：根阶段终局汇合点（先发射再返回；此时尚无前沿栈 ⇒ bound=none）。
        if let Some(h) = hooks {
            emit_progress(h, &start, stats.nodes, incumbent.as_ref().map(|(o, _)| *o), None, sense);
        }
        return Ok(finish(status, None, &mut stats, &simplex, start, None, false, false));
    }

    // —— 根割循环（可选；割行全局保留至树结束）——
    if opts.root_cuts {
        let (applied, st) = cuts::cut_loop(
            &mut simplex,
            &bounds,
            problem,
            &int_vars,
            CutOptions::default(),
            &tol,
        )?;
        stats.cuts = applied;
        if let Some(status) = root_terminal(classify_node(st, true))? {
            // 首行进度：根阶段终局汇合点（先发射再返回；此时尚无前沿栈 ⇒ bound=none）。
            if let Some(h) = hooks {
                emit_progress(h, &start, stats.nodes, incumbent.as_ref().map(|(o, _)| *o), None, sense);
            }
            return Ok(finish(
                status,
                None,
                &mut stats,
                &simplex,
                start,
                None,
                false,
                false,
            ));
        }
        // 取回割后的根 LP 解（记账性 resolve，不计节点）。
        sol = simplex.resolve(&bounds);
        if let Some(status) = root_terminal(classify_node(sol.status, true))? {
            // 首行进度：根阶段终局汇合点（先发射再返回；此时尚无前沿栈 ⇒ bound=none）。
            if let Some(h) = hooks {
                emit_progress(h, &start, stats.nodes, incumbent.as_ref().map(|(o, _)| *o), None, sense);
            }
            return Ok(finish(
                status,
                None,
                &mut stats,
                &simplex,
                start,
                None,
                false,
                false,
            ));
        }
    }

    // —— 根启发式（可选）：round_and_check + dive 供应初始 incumbent ——
    if opts.heuristics {
        if let Some(inc) = heuristics::round_and_check(problem, &bounds, &sol.x, &tol) {
            maybe_update(sense, &mut incumbent, inc);
        }
        if let Some(inc) = heuristics::dive(
            problem,
            &mut simplex,
            &bounds,
            &sol.x,
            opts.dive_max_fixes,
            &tol,
        ) {
            maybe_update(sense, &mut incumbent, inc);
        }
        // dive 扰动了 Simplex 内部基：用真实界恢复根热状态（记账性 resolve）。
        sol = simplex.resolve(&bounds);
        if let Some(status) = root_terminal(classify_node(sol.status, true))? {
            // 首行进度：根阶段终局汇合点（先发射再返回；此时尚无前沿栈 ⇒ bound=none）。
            if let Some(h) = hooks {
                emit_progress(h, &start, stats.nodes, incumbent.as_ref().map(|(o, _)| *o), None, sense);
            }
            return Ok(finish(
                status,
                None,
                &mut stats,
                &simplex,
                start,
                None,
                false,
                false,
            ));
        }
    }

    // 根启发式供应的 incumbent 记为根阶段末尾时刻的改进（降级停滞锚点）。
    if let Some(state) = adaptive.as_mut() {
        if incumbent.is_some() {
            state.last_improvement_secs = start.elapsed().as_secs_f64();
        }
    }

    // —— DFS 主循环状态 ——
    let mut stack: Vec<Pending> = Vec::new();
    let mut trail: Vec<TrailEntry> = Vec::new();
    let mut limit_hit = false;
    // 策略级自动降级触发标记（与 gap_hit 同径终止，报告侧另附原因）。
    let mut degraded = false;

    // 根节点按 (d)/(e)/(f) 判定（无界等已终局，此处的 Unbounded 恒为根语义）。
    branch_or_close(
        problem,
        &int_vars,
        &bounds,
        &sol,
        tol,
        opts,
        &mut incumbent,
        &mut stack,
        &mut trail,
        &mut simplex,
        stats.nodes,
        start.elapsed().as_secs_f64(),
        adaptive.as_mut(),
    );

    // 根阶段完成（根启发式刚运行过）：控制器节流锚点 = 根节点数（=1）。
    if let Some(state) = adaptive.as_mut() {
        state.controller.notify_root_done(stats.nodes);
    }

    // —— 根阶段结束：首行进度发射点（钩子启用时无条件发射，口径与并行一致；
    // 先发射再进入主循环。栈空时若已有 incumbent，搜索事实上完备 ⇒ bound =
    // incumbent（gap=0），杜绝 incumbent=X bound=none gap=inf 矛盾行）——
    if let Some(h) = hooks {
        let inc_obj = incumbent.as_ref().map(|(o, _)| *o);
        let bound = frontier_bound(sense, &stack, inc_obj);
        emit_progress(h, &start, stats.nodes, inc_obj, bound, sense);
        last_log = Some(start.elapsed().as_secs_f64());
    }

    loop {
        // —— 受控扩展（默认全关 ⇒ 与原实现逐字节一致）。①/② 仅在栈非空时
        // 参与：空栈轮次直接落到下方 pop 的自然 break（Optimal/Infeasible），
        // 杜绝 incumbent=X bound=none gap=inf 矛盾周期行——
        if !stack.is_empty() && !gap_hit {
            // ① 周期进度行：log 启用且距上行 ≥ log_period_secs 才发射。
            if let (Some(h), Some(period)) = (hooks, ctrl.log_period_secs) {
                let now = start.elapsed().as_secs_f64();
                if period > 0.0 && last_log.is_some_and(|t| now - t >= period) {
                    let inc_obj = incumbent.as_ref().map(|(o, _)| *o);
                    let db = aggregate_bound(sense, stack.iter().map(|q| q.est));
                    emit_progress(h, &start, stats.nodes, inc_obj, db, sense);
                    last_log = Some(now);
                }
            }
            // ② mipgap：gap（经统一 clamp ≥ 0）≤ P ⇒ 粘滞置位并跳出，记下
            // 触发时刻证据 (db, inc)——与并行证据口径一致，不依赖后续栈态。
            if let Some(p) = ctrl.mipgap_percent {
                if p > 0.0 {
                    let inc_obj = incumbent.as_ref().map(|(o, _)| *o);
                    let db = aggregate_bound(sense, stack.iter().map(|q| q.est));
                    if let (Some(inc), Some(db)) = (inc_obj, db) {
                        if let Some(g) = gap_percent(sense, inc, db) {
                            if g <= p {
                                gap_hit = true;
                                gap_evidence = Some(db);
                                break;
                            }
                        }
                    }
                }
            }
            // ③ 策略级自动降级（仅 auto 模式且给了时间预算）：资格预检通过
            // 才聚合前沿界（避免逐节点 O(栈) 开销）；双停滞（incumbent 与界
            // 都不再前进）⇒ 与 ② 同径提前终止，证据口径一致。
            if let (Some(state), true) = (adaptive.as_mut(), opts.time_limit_ms != 0) {
                let now = start.elapsed().as_secs_f64();
                let budget = opts.time_limit_ms as f64 / 1000.0;
                if AutoDegrade::eligible(state.degrade.cfg(), now, Some(budget)) {
                    let db = aggregate_bound(sense, stack.iter().map(|q| q.est));
                    let snap = DegradeSnapshot {
                        sense,
                        elapsed_secs: now,
                        budget_secs: Some(budget),
                        has_incumbent: incumbent.is_some(),
                        last_improvement_secs: state.last_improvement_secs,
                        bound: db,
                    };
                    if state.degrade.decide(&snap) {
                        gap_hit = true;
                        degraded = true;
                        gap_evidence = db;
                        break;
                    }
                }
            }
        }
        // 栈空 ⇒ 搜索自然耗尽（先于上限判定，见模块实现说明）。
        let Some(p) = stack.pop() else {
            break;
        };
        // 节点/时间上限：每次 resolve 前检查；0 = 不限。
        if opts.node_limit != 0 && stats.nodes >= opts.node_limit {
            limit_hit = true;
            break;
        }
        if opts.time_limit_ms != 0 && time_exceeded(&start, opts.time_limit_ms) {
            limit_hit = true;
            break;
        }
        // 回退到父路径界状态后施加本分支界；空界/无收缩 ⇒ (f) 剪枝（不 resolve）。
        unwind(&mut trail, &mut bounds, p.unwind_to);
        if !tighten_child(&mut bounds, &mut trail, p.var, p.tighten_lb, p.bound) {
            continue;
        }
        let sol = simplex.resolve(&bounds);
        stats.nodes += 1;

        // (a)/(b)/(c)：resolve 状态判定（顺序固定）。
        match classify_node(sol.status, false) {
            NodeAction::Prune => continue,
            NodeAction::PruneLimit => {
                limit_hit = true;
                continue;
            }
            NodeAction::FailUnbounded => return Err(nonroot_unbounded_error()),
            NodeAction::TreeUnbounded => {
                return Err(MilpError::Internal(
                    "非根节点分类不可能返回根无界动作".to_string(),
                ));
            }
            NodeAction::Process => {}
        }

        // (d)/(e)/(f) 与分支压栈。
        branch_or_close(
            problem,
            &int_vars,
            &bounds,
            &sol,
            tol,
            opts,
            &mut incumbent,
            &mut stack,
            &mut trail,
            &mut simplex,
            stats.nodes,
            start.elapsed().as_secs_f64(),
            adaptive.as_mut(),
        );
    }

    // —— 终态判定（与上限诚实性整合，禁止字面分裂）：gap 达标与上限触发
    // 一律 LimitReached（诚实：未证明最优），粘滞标记优先于完备判定——
    let status = if limit_hit || gap_hit {
        SolveStatus::LimitReached
    } else if incumbent.is_some() {
        SolveStatus::Optimal
    } else {
        SolveStatus::Infeasible
    };
    // —— bound 口径：gap 终止 = 触发时刻证据 db；截断 = 剩余栈 est 聚合
    //（可 None）；完备有 incumbent = best_objective；Infeasible = None——
    let best_objective = incumbent.as_ref().map(|(o, _)| *o);
    let bound = if gap_hit {
        gap_evidence
    } else if limit_hit {
        aggregate_bound(sense, stack.iter().map(|q| q.est))
    } else if status == SolveStatus::Optimal {
        best_objective
    } else {
        None
    };
    Ok(finish(
        status,
        incumbent,
        &mut stats,
        &simplex,
        start,
        bound,
        gap_hit,
        degraded,
    ))
}

// ===========================================================================
// 节点分类与判定（纯函数，便于单元级直接测试）
// ===========================================================================

/// 节点 resolve 状态 → 动作（模块内纯函数）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum NodeAction {
    /// (a) LP 不可行 ⇒ 剪枝回溯。
    Prune,
    /// (c) 迭代上限 ⇒ 剪枝并置 limit_hit 标记。
    PruneLimit,
    /// (b) 根节点无界 ⇒ 整树终止返回 [`SolveStatus::Unbounded`]。
    TreeUnbounded,
    /// (b) 非根节点无界（父有界时理论不发生的数值异常）⇒ 整树终止
    /// `Err(MilpError::Numerical)`，不得静默继续。
    FailUnbounded,
    /// Optimal ⇒ 继续 (d)/(e)/(f) 判定。
    Process,
}

/// 节点 resolve 状态分类：(a) 不可行、(b) 无界（根/非根分流）、(c) 迭代上限。
fn classify_node(status: LpStatus, is_root: bool) -> NodeAction {
    match status {
        LpStatus::Infeasible => NodeAction::Prune,
        LpStatus::IterationLimit => NodeAction::PruneLimit,
        LpStatus::Unbounded => {
            if is_root {
                NodeAction::TreeUnbounded
            } else {
                NodeAction::FailUnbounded
            }
        }
        LpStatus::Optimal => NodeAction::Process,
    }
}

/// 非根节点 LP 无界的数值错误：子节点仅为父界的收紧，收紧不可能使有界 LP
/// 变无界；出现即数值异常，整树终止。
fn nonroot_unbounded_error() -> MilpError {
    MilpError::Numerical("node LP unbounded after bounded root".to_string())
}

/// 根节点动作 → 整树终态（`None` = 继续树搜索）。
fn root_terminal(action: NodeAction) -> Result<Option<SolveStatus>> {
    match action {
        NodeAction::Process => Ok(None),
        NodeAction::Prune => Ok(Some(SolveStatus::Infeasible)),
        NodeAction::PruneLimit => Ok(Some(SolveStatus::LimitReached)),
        NodeAction::TreeUnbounded => Ok(Some(SolveStatus::Unbounded)),
        NodeAction::FailUnbounded => Err(MilpError::Internal(
            "根节点分类不可能返回非根无界动作".to_string(),
        )),
    }
}

/// (e) 整数节点判定结果。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum NodeEval {
    /// 复核通过 ⇒ 调用方更新 incumbent（若严格更优）后回溯。
    Record,
    /// 非整数，或伪整数复核失败 ⇒ 继续分支。
    Branch,
}

/// (e) 纯判定：`integral`（整数性满足）且 `check_passed`（round_and_check
/// 复核通过）⇒ [`NodeEval::Record`]；否则 [`NodeEval::Branch`]。
///
/// 伪整数回退：复核不通过时不得剪枝、不得记 incumbent，必须继续分支——
/// 堵死“丢最优/可行问题误报 Infeasible”路径。
fn integer_node_eval(integral: bool, check_passed: bool) -> NodeEval {
    if integral && check_passed {
        NodeEval::Record
    } else {
        NodeEval::Branch
    }
}

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

// ===========================================================================
// 界栈（trail）与分支决策
// ===========================================================================

/// trail 记录：某次界收紧前的旧值（逆序回放恢复）。
struct TrailEntry {
    var: usize,
    is_lb: bool,
    old: f64,
}

/// 待应用的分支决策（stack 元素）。
struct Pending {
    var: usize,
    /// true = 上行分支（收紧 lb 至 ⌊x̄⌋+1）；false = 下行（收紧 ub 至 ⌊x̄⌋）。
    tighten_lb: bool,
    bound: f64,
    /// 应用本决策前 trail 应回退到的长度（即父路径的 trail 状态）。
    unwind_to: usize,
    /// 压栈时父节点 LP 目标（该子树最优值的上界估计）。纯数据：不参与任何
    /// 决策，仅供受控口径（对偶界聚合 [`crate::search_ctrl::aggregate_bound`] /
    /// gap 判定 / 进度行）读取。
    est: f64,
}

/// 回退 trail 至长度 `to`（严格逆序恢复旧界）。
fn unwind(trail: &mut Vec<TrailEntry>, bounds: &mut Bounds, to: usize) {
    while trail.len() > to {
        let Some(e) = trail.pop() else {
            break;
        };
        if e.is_lb {
            bounds.lb[e.var] = e.old;
        } else {
            bounds.ub[e.var] = e.old;
        }
    }
}

/// 施加子节点分支界：相对当前界严格收缩且新区间非空才应用并记账 trail，
/// 返回 true；空界（(f)）或无收缩（无进展）返回 false，调用方跳过该分支
/// （不 resolve、不计节点）。
fn tighten_child(
    bounds: &mut Bounds,
    trail: &mut Vec<TrailEntry>,
    var: usize,
    tighten_lb: bool,
    bound: f64,
) -> bool {
    if tighten_lb {
        if bound <= bounds.lb[var] || bound > bounds.ub[var] {
            return false;
        }
        trail.push(TrailEntry {
            var,
            is_lb: true,
            old: bounds.lb[var],
        });
        bounds.lb[var] = bound;
    } else {
        if bound >= bounds.ub[var] || bound < bounds.lb[var] {
            return false;
        }
        trail.push(TrailEntry {
            var,
            is_lb: false,
            old: bounds.ub[var],
        });
        bounds.ub[var] = bound;
    }
    true
}

// ===========================================================================
// Optimal 节点处理与 incumbent 维护
// ===========================================================================

/// Optimal 节点的 (d)/(e) 判定与分支压栈：
/// - (d) 界不优于 incumbent ⇒ 节点关闭；
/// - (e) 整数性满足 ⇒ round_and_check 复核，通过 ⇒ 更新 incumbent（若严格
///   更优）后节点关闭；不通过（伪整数）⇒ 落入分支；
/// - （受控扩展，默认关闭）自适应启发式：(d) 未剪枝且本节点未按整数记录
///   （即落入分支路径）、选分支变量之前，控制器 decide = Run ⇒ 在本节点
///   补跑 dive（是否产生新全局 incumbent 作为 hit record）并用真实界
///   resolve 回取（记账性，不计节点）；
/// - 分支：MostFractional 选变量，下行（x ≤ ⌊x̄⌋）后压 ⇒ 先探索。
#[allow(clippy::too_many_arguments)]
fn branch_or_close(
    problem: &Problem,
    int_vars: &[usize],
    bounds: &Bounds,
    sol: &LpSolution,
    tol: Tolerances,
    opts: &BbOptions,
    incumbent: &mut Option<(f64, Vec<f64>)>,
    stack: &mut Vec<Pending>,
    trail: &mut Vec<TrailEntry>,
    simplex: &mut Simplex,
    nodes: u64,
    elapsed_secs: f64,
    mut adaptive: Option<&mut AdaptiveState>,
) {
    // (d)：LP bound 不优于 incumbent ⇒ 剪枝（按原始 sense 方向比较）。
    if let Some((obj, _)) = incumbent.as_ref() {
        if bound_no_better(
            problem.sense,
            sol.objective,
            *obj,
            opts.prune_abs,
            opts.prune_rel,
        ) {
            return;
        }
    }
    // (e)：整数性满足 ⇒ round_and_check 复核（目标按取整后的值评估）。
    let integral = is_integral(int_vars, &sol.x, tol.integrality);
    let mut check_passed = false;
    let mut improved = false;
    if integral {
        if let Some(inc) = heuristics::round_and_check(problem, bounds, &sol.x, &tol) {
            improved = maybe_update(problem.sense, incumbent, inc);
            check_passed = true;
        }
    }
    if integer_node_eval(integral, check_passed) != NodeEval::Branch {
        // 复核通过 ⇒ 回溯（节点关闭）；incumbent 改进刷新自适应锚点。
        if let Some(state) = adaptive {
            if improved {
                state.last_improvement_node = nodes;
                state.last_improvement_secs = elapsed_secs;
            }
        }
        return;
    }
    // 自适应启发式（仅分支路径、选分支变量之前；adaptive=None ⇒ 零行为）。
    if let Some(state) = adaptive.as_deref_mut() {
        let snap = SearchSnapshot {
            nodes,
            incumbent: incumbent.as_ref().map(|(o, _)| *o),
            nodes_since_improvement: nodes.saturating_sub(state.last_improvement_node),
            time_budget_secs: if opts.time_limit_ms == 0 {
                None
            } else {
                Some(opts.time_limit_ms as f64 / 1000.0)
            },
            elapsed_secs,
        };
        if state.controller.decide(&snap) == HeuristicDecision::Run {
            let hit = match heuristics::dive(
                problem,
                simplex,
                bounds,
                &sol.x,
                opts.dive_max_fixes,
                &tol,
            ) {
                Some(inc) => maybe_update(problem.sense, incumbent, inc),
                None => false,
            };
            state.controller.record(hit);
            if hit {
                state.last_improvement_node = nodes;
                state.last_improvement_secs = elapsed_secs;
            }
            // dive 扰动了 Simplex 内部基：用本节点真实界 resolve 回取
            //（记账性，不计节点；分支决策仍基于传入的本节点 LP 解 sol）。
            let _ = simplex.resolve(bounds);
        }
    }
    // 分支：MostFractional；无可分支变量为内部不可达防御（不剪枝即丢弃节点）。
    let Some((var, val)) = select_branch_var(int_vars, &sol.x) else {
        return;
    };
    let floor = val.floor();
    let unwind_to = trail.len();
    // 上行先压（后弹出）、下行后压（先探索）；est = 压栈时父节点 LP 目标。
    stack.push(Pending {
        var,
        tighten_lb: true,
        bound: floor + 1.0,
        unwind_to,
        est: sol.objective,
    });
    stack.push(Pending {
        var,
        tighten_lb: false,
        bound: floor,
        unwind_to,
        est: sol.objective,
    });
}

/// 顺序路径的自适应启发式状态：WP1 控制器 + incumbent 改进锚点（供
/// [`SearchSnapshot::nodes_since_improvement`] 与降级停滞判定使用）。
struct AdaptiveState {
    controller: HeuristicController,
    last_improvement_node: u64,
    /// 最近一次 incumbent 严格改进的时刻（秒；根启发式供应的解记根阶段
    /// 末尾时刻）。策略级自动降级的双停滞信号之一。
    last_improvement_secs: f64,
    /// 策略级自动降级判定器（窗口界基准有状态）。
    degrade: AutoDegrade,
}

impl AdaptiveState {
    fn new(cfg: HeuristicConfig, degrade_policy: DegradePolicy) -> Self {
        AdaptiveState {
            controller: HeuristicController::new(cfg),
            last_improvement_node: 0,
            last_improvement_secs: 0.0,
            degrade: AutoDegrade::new(degrade_policy),
        }
    }
}

/// 严格更优才替换 incumbent（按原始 sense 方向；目标与取值均为 WP3 取整后
/// 口径，防报不可行值）。返回是否发生了替换（自适应控制器的 hit 语义）。
fn maybe_update(sense: Sense, incumbent: &mut Option<(f64, Vec<f64>)>, inc: Incumbent) -> bool {
    if !inc.objective.is_finite() {
        return false; // 防御：非有限目标不入 incumbent
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
    better
}

/// 时间上限检查（调用方保证 `limit_ms > 0`）。
fn time_exceeded(start: &Instant, limit_ms: u64) -> bool {
    start.elapsed().as_millis() >= limit_ms as u128
}

/// 组装结果并补全统计（耗时与累计迭代）；bound/gap_terminated 口径由调用方
/// 按终态类别给出（见 [`BbResult::bound`] 字段说明）。
fn finish(
    status: SolveStatus,
    incumbent: Option<(f64, Vec<f64>)>,
    stats: &mut BbStats,
    simplex: &Simplex,
    start: Instant,
    bound: Option<f64>,
    gap_terminated: bool,
    auto_degraded: bool,
) -> BbResult {
    stats.lp_iterations = simplex.iterations();
    stats.time_ms = u64::try_from(start.elapsed().as_millis()).unwrap_or(u64::MAX);
    let (best_objective, best_values) = match incumbent {
        Some((obj, values)) => (Some(obj), Some(values)),
        None => (None, None),
    };
    BbResult {
        status,
        best_objective,
        best_values,
        stats: *stats,
        bound,
        gap_terminated,
        auto_degraded,
    }
}

// ===========================================================================
// 受控口径的纯辅助（进度行与对偶界聚合；gap 统一走 search_ctrl::gap_percent）
// ===========================================================================

/// 首行进度的 bound 口径：栈非空 ⇒ 剩余聚合（每个开放子树的最优值不会优于
/// 其压栈时的父节点 LP 目标 est）；栈空且已有 incumbent ⇒ 搜索事实上完备，
/// bound = incumbent（gap=0）。
fn frontier_bound(sense: Sense, stack: &[Pending], incumbent: Option<f64>) -> Option<f64> {
    aggregate_bound(sense, stack.iter().map(|p| p.est)).or(incumbent)
}

/// 组装并发射一行进度（字段口径与并行引擎一致；gap 由展示层经
/// [`gap_percent`] 现算，缺任一侧或 inf 语义渲染为 inf）。
fn emit_progress(
    hooks: &SearchHooks,
    start: &Instant,
    nodes: u64,
    incumbent: Option<f64>,
    bound: Option<f64>,
    sense: Sense,
) {
    (hooks.progress)(&ProgressInfo {
        sense,
        elapsed_secs: start.elapsed().as_secs_f64(),
        nodes,
        incumbent,
        bound,
    });
}

// ===========================================================================
// 并行委托（threads ≥ 2；引擎本体见 crate::parallel_bb，本文件只做组装与映射）
// ===========================================================================

/// threads ≥ 2 的受控入口：由 [`SearchControl`] 组装
/// [`crate::parallel_bb::ParallelConfig`]（BbOptions 原样透传），委托
/// [`crate::parallel_bb::parallel_branch_and_bound`]，再按本文件统一口径把
/// [`crate::parallel_bb::ParallelOutcome`] 映射回 [`BbResult`]。
fn parallel_entry(
    problem: &Problem,
    opts: &BbOptions,
    ctrl: &SearchControl,
    tol: Tolerances,
    hooks: Option<&SearchHooks>,
) -> Result<BbResult> {
    let config = crate::parallel_bb::ParallelConfig {
        threads: ctrl.threads,
        mipgap_percent: ctrl.mipgap_percent,
        adaptive_heuristics: ctrl.adaptive_heuristics,
        log_period_secs: ctrl.log_period_secs,
        degrade: ctrl.degrade,
        progress: hooks.map(|h| h.progress.clone()),
    };
    let outcome = crate::parallel_bb::parallel_branch_and_bound(problem, opts, &config, tol)?;
    Ok(map_parallel_outcome(outcome))
}

/// 并行终态 → [`BbResult`]（与顺序路径口径逐一对应）：
/// - 根阶段终局（root_terminal）⇒ 直返（bound 引擎已按口径给出：纯 LP
///   Optimal = Some(objective)，其余 None）；
/// - stop ∈ {NodeLimit, TimeLimit}（及其他截断类）⇒ LimitReached，bound =
///   引擎上报的剩余聚合（截断不可知时诚实 None）；
/// - stop = Gap ⇒ LimitReached 且 gap_terminated = true，bound = 触发时刻
///   证据快照（引擎以 CAS 胜者独占写入，绝不join后重聚合）；
/// - stop = None ⇒ incumbent? Optimal : Infeasible（粘滞 stop 优先于完备），
///   完备有 incumbent 时 bound = best_objective。
fn map_parallel_outcome(outcome: crate::parallel_bb::ParallelOutcome) -> BbResult {
    use crate::parallel_bb::ParallelStop;
    // 根阶段终局 ⇒ 直返。
    if let Some(status) = outcome.root_terminal {
        return BbResult {
            status,
            best_objective: outcome.best_objective,
            best_values: outcome.best_values,
            stats: outcome.stats,
            bound: outcome.bound,
            gap_terminated: false,
            auto_degraded: false,
        };
    }
    // 粘滞 stop 优先于完备判定。NodeLimit/TimeLimit（含迭代上限类截断，
    // 引擎侧统一映射为 NodeLimit 停止码）⇒ 诚实 LimitReached。
    let (status, gap_terminated, auto_degraded, bound) = match outcome.stop {
        Some(ParallelStop::Gap) => (SolveStatus::LimitReached, true, false, outcome.bound),
        Some(ParallelStop::AutoDegrade) => (SolveStatus::LimitReached, true, true, outcome.bound),
        Some(ParallelStop::NodeLimit | ParallelStop::TimeLimit) => {
            (SolveStatus::LimitReached, false, false, outcome.bound)
        }
        None => {
            if outcome.best_objective.is_some() {
                (SolveStatus::Optimal, false, false, outcome.bound)
            } else {
                (SolveStatus::Infeasible, false, false, None)
            }
        }
    };
    BbResult {
        status,
        best_objective: outcome.best_objective,
        best_values: outcome.best_values,
        stats: outcome.stats,
        bound,
        gap_terminated,
        auto_degraded,
    }
}

// ===========================================================================
// 单元测试
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RelOp, VarType};

    /// 默认容差。
    fn tol() -> Tolerances {
        Tolerances::default()
    }

    /// 关闭根割与根启发式的选项：隔离分支定界本体行为
    /// （整数节点复核 round_and_check 属跨模块联调路径，见相关用例）。
    fn plain_opts() -> BbOptions {
        BbOptions {
            root_cuts: false,
            heuristics: false,
            ..BbOptions::default()
        }
    }

    fn assert_close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 1e-7 * (1.0 + b.abs()),
            "数值不符: {a} vs {b}"
        );
    }

    /// BbOptions 默认值逐字段符合契约。
    #[test]
    fn bb_options_default_contract() {
        let o = BbOptions::default();
        assert_eq!(o.branch_rule, BranchRule::MostFractional);
        assert_eq!(o.node_rule, NodeRule::DepthFirst);
        assert_eq!(o.node_limit, 100_000);
        assert_eq!(o.time_limit_ms, 0);
        assert!(o.root_cuts);
        assert!(o.heuristics);
        assert_eq!(o.dive_max_fixes, 8);
        assert_eq!(o.prune_abs, 1e-9);
        assert_eq!(o.prune_rel, 1e-9);
    }

    /// 节点分类纯函数：(Unbounded, is_root=false) → FailUnbounded（整树终止
    /// Err(Numerical) 分支）；(Unbounded, is_root=true) → TreeUnbounded
    /// （→ SolveStatus::Unbounded 分支）。非根 Unbounded 数学上不可用算例
    /// 构造（父有界时收紧界不可能无界），故只对纯函数断言。
    #[test]
    fn classify_node_unbounded_root_vs_nonroot() {
        assert_eq!(
            classify_node(LpStatus::Unbounded, true),
            NodeAction::TreeUnbounded
        );
        assert_eq!(
            classify_node(LpStatus::Unbounded, false),
            NodeAction::FailUnbounded
        );
        // FailUnbounded 动作到错误的映射。
        assert!(
            matches!(nonroot_unbounded_error(), MilpError::Numerical(msg) if msg.contains("unbounded"))
        );
        // 其余映射。
        assert_eq!(
            classify_node(LpStatus::Infeasible, false),
            NodeAction::Prune
        );
        assert_eq!(classify_node(LpStatus::Infeasible, true), NodeAction::Prune);
        assert_eq!(
            classify_node(LpStatus::IterationLimit, true),
            NodeAction::PruneLimit
        );
        assert_eq!(classify_node(LpStatus::Optimal, false), NodeAction::Process);
        // 根动作 → 终态。
        assert_eq!(
            root_terminal(NodeAction::TreeUnbounded).unwrap(),
            Some(SolveStatus::Unbounded)
        );
        assert_eq!(
            root_terminal(NodeAction::Prune).unwrap(),
            Some(SolveStatus::Infeasible)
        );
        assert_eq!(
            root_terminal(NodeAction::PruneLimit).unwrap(),
            Some(SolveStatus::LimitReached)
        );
        assert_eq!(root_terminal(NodeAction::Process).unwrap(), None);
    }

    /// 伪整数复核失败回退判定（纯函数）：integral 且复核不通过 ⇒ 必须分支
    /// （不得剪枝、不得记 incumbent）；复核通过 ⇒ Record；非整数 ⇒ Branch。
    #[test]
    fn pseudo_integer_fallback_branches() {
        assert_eq!(integer_node_eval(true, false), NodeEval::Branch);
        assert_eq!(integer_node_eval(true, true), NodeEval::Record);
        assert_eq!(integer_node_eval(false, false), NodeEval::Branch);
        assert_eq!(integer_node_eval(false, true), NodeEval::Branch);
    }

    /// (d) 剪枝比较按原始 sense 方向（Max/Min 镜像）。
    #[test]
    fn bound_prune_respects_sense() {
        let (pa, pr) = (1e-9, 1e-9);
        // Max：bound ≤ inc + slack ⇒ 剪枝。
        assert!(bound_no_better(Sense::Maximize, 10.0, 10.0, pa, pr));
        assert!(!bound_no_better(Sense::Maximize, 10.5, 10.0, pa, pr));
        // Min：bound ≥ inc − slack ⇒ 剪枝。
        assert!(bound_no_better(Sense::Minimize, 2.0, 2.0, pa, pr));
        assert!(!bound_no_better(Sense::Minimize, 1.5, 2.0, pa, pr));
    }

    /// MostFractional 选择与并列最小下标破键。
    #[test]
    fn most_fractional_selection_and_tie_break() {
        let iv = [0usize, 1, 2];
        // 得分 0.5 / 0.2 / 0.2 → 变量 0。
        assert_eq!(select_branch_var(&iv, &[1.5, 2.2, 0.8]), Some((0, 1.5)));
        // 1.25/3.75 二进制精确表示：得分恰为 0.25 / 0.25 / 0 的真并列
        // （1.2/3.8 等十进制值的分数性受 fl 噪声影响并非精确并列）。
        assert_eq!(select_branch_var(&iv, &[1.25, 3.75, 0.0]), Some((0, 1.25)));
        // 得分 0 / 0.3 / 0 → 变量 1。
        assert_eq!(select_branch_var(&iv, &[2.0, 1.3, 4.0]), Some((1, 1.3)));
        // 负值：floor(−1.5) = −2，f = 0.5。
        assert_eq!(select_branch_var(&iv, &[-1.5, 0.0, 0.0]), Some((0, -1.5)));
        assert_eq!(select_branch_var(&[], &[1.5]), None);
    }

    /// 空整数集（纯 LP）走捷径：单节点、无割，直接按 LP 状态返回。
    #[test]
    fn empty_int_vars_pure_lp_shortcut() {
        // max 3x+2y s.t. x+y<=4, x+3y<=6, [0,3]^2 → (3,1)，目标 11。
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
        let r = branch_and_bound(&p, &BbOptions::default(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_close(r.best_objective.unwrap(), 11.0);
        let v = r.best_values.unwrap();
        assert_eq!(v.len(), 2);
        assert_close(v[0], 3.0);
        assert_close(v[1], 1.0);
        assert_eq!(r.stats.nodes, 1);
        assert_eq!(r.stats.cuts, 0);
    }

    /// 根 Infeasible 透传：不触割循环与启发式（stats.cuts == 0）。
    #[test]
    fn root_infeasible_passthrough_no_cuts() {
        let mut p = Problem::new("inf", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Ge, 2.5)
            .unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 2.0)
            .unwrap();
        let r = branch_and_bound(&p, &BbOptions::default(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Infeasible);
        assert!(r.best_objective.is_none());
        assert!(r.best_values.is_none());
        assert_eq!(r.stats.nodes, 1);
        assert_eq!(r.stats.cuts, 0, "根非 Optimal 不得进入割循环");
    }

    /// 根 Unbounded 透传（整数变量存在而根 LP 松弛无界 ⇒ 按契约近似直接判
    /// Unbounded）：不触割循环与启发式。
    #[test]
    fn root_unbounded_passthrough_no_cuts() {
        let mut p = Problem::new("unb", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        let r = branch_and_bound(&p, &BbOptions::default(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Unbounded);
        assert!(r.best_objective.is_none());
        assert!(r.best_values.is_none());
        assert_eq!(r.stats.nodes, 1);
        assert_eq!(r.stats.cuts, 0, "根非 Optimal 不得进入割循环");
    }

    /// 极小 node_limit 驱动 limit_hit ⇒ LimitReached（附 incumbent，若有）：
    /// 绝不把截断搜索谎报为 Optimal/Infeasible。
    #[test]
    fn tiny_node_limit_reports_limit_reached() {
        let make = || {
            let mut p = Problem::new("lim", Sense::Maximize);
            p.add_var("x", VarType::Integer).unwrap();
            p.add_var("y", VarType::Integer).unwrap();
            p.set_objective(vec![(0, 1.0), (1, 1.0)], 0.0).unwrap();
            p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 3.5)
                .unwrap();
            p
        };
        // node_limit = 1：仅根节点（LP 最优 x=3.5 分数）。
        let opts = BbOptions {
            node_limit: 1,
            ..plain_opts()
        };
        let r = branch_and_bound(&make(), &opts, tol()).unwrap();
        assert_eq!(r.status, SolveStatus::LimitReached);
        assert_eq!(r.stats.nodes, 1);
        assert!(
            r.best_objective.is_none(),
            "未取得 incumbent 时不带目标返回"
        );
        // node_limit = 2：根 + 下行子节点（x=3, y=0.5 仍分数）后截断。
        let opts = BbOptions {
            node_limit: 2,
            ..plain_opts()
        };
        let r = branch_and_bound(&make(), &opts, tol()).unwrap();
        assert_eq!(r.status, SolveStatus::LimitReached);
        assert_eq!(r.stats.nodes, 2);
    }

    /// 整数变量界夹在分数区间内（无可行整数点）：两个分支界均空/无收缩 ⇒
    /// 不压栈、不 resolve，Infeasible 且仅 1 节点。
    #[test]
    fn infeasible_integer_band_no_children() {
        let mut p = Problem::new("band", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.4, 0.6).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Infeasible);
        assert_eq!(r.stats.nodes, 1, "空界分支不得 resolve");
    }

    /// 兄弟节点 LP 不可行剪枝 + trail 回溯正确性：根 2.6 分数，下行 [0,2] 与
    /// 上行 [3,5] 的 LP 均不可行 ⇒ 全树 3 节点后 Infeasible。
    #[test]
    fn sibling_lp_infeasible_prune_and_backtrack() {
        let mut p = Problem::new("sib", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 5.0).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Ge, 2.5)
            .unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 2.6)
            .unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Infeasible);
        assert_eq!(r.stats.nodes, 3, "根 + 两个 LP 不可行子节点");
        assert!(r.best_values.is_none());
    }

    /// （跨模块联调：整数节点复核经 WP3 heuristics::round_and_check）
    /// 单变量 MILP 全路径：max x s.t. 2x≤7 → 最优 3，3 节点。
    #[test]
    fn small_milp_optimal_integration() {
        let mut p = Problem::new("m1", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 2.0)], RelOp::Le, 7.0)
            .unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_close(r.best_objective.unwrap(), 3.0);
        let v = r.best_values.unwrap();
        assert_eq!(v.len(), 1);
        assert_close(v[0], 3.0);
        assert_eq!(r.stats.nodes, 3);
    }

    /// （跨模块联调）Min 方向镜像：min x s.t. x≥2.5 → 最优 3，3 节点
    /// （下行子节点 LP 不可行，上行取得整数 incumbent）。
    #[test]
    fn minimize_sense_integration() {
        let mut p = Problem::new("min1", Sense::Minimize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Ge, 2.5)
            .unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_close(r.best_objective.unwrap(), 3.0);
        assert_close(r.best_values.unwrap()[0], 3.0);
        assert_eq!(r.stats.nodes, 3);
    }

    /// （跨模块联调）两层分支 + incumbent 两次更新 + 深层回溯：
    /// max 10a+7b s.t. 5a+3b≤7, a,b∈{0,1} → LP 松弛 15（a=0.8,b=1），
    /// 整数最优 10（a=1,b=0），共 5 节点。
    #[test]
    fn binary_knapsack_two_level_integration() {
        let mut p = Problem::new("k2", Sense::Maximize);
        p.add_var("a", VarType::Binary).unwrap();
        p.add_var("b", VarType::Binary).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 7.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 5.0), (1, 3.0)], RelOp::Le, 7.0)
            .unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_close(r.best_objective.unwrap(), 10.0);
        let v = r.best_values.unwrap();
        assert_close(v[0], 1.0);
        assert_close(v[1], 0.0);
        assert_eq!(r.stats.nodes, 5);
    }

    /// （跨模块联调）行范数大（‖row‖₁=1e6）、整数偏差≈5e-7 的近整数根节点：
    /// 取整后行违反 0，复核应通过 ⇒ 根节点直接收敛（1 节点、目标 3）。
    #[test]
    fn near_integer_root_accept_integration() {
        let mut p = Problem::new("near1", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0e6)], RelOp::Le, 3_000_000.5)
            .unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_eq!(r.stats.nodes, 1, "近整数根节点经复核通过后不再分支");
        assert_close(r.best_objective.unwrap(), 3.0);
        assert_close(r.best_values.unwrap()[0], 3.0);
    }

    /// （跨模块联调）行范数大、取整造成行违反的伪整数边缘算例：
    /// LP 解 2.9999995（偏差≈5e-7 ≤ integrality），取整为 3 使行违反 0.5。
    /// - WP3 阶梯复核（violation ≤ 10feas + integrality×‖row‖₁ ≈ 1.0）通过
    ///   ⇒ 根节点记录 incumbent=3；
    /// - 复核若拒绝 ⇒ 伪整数回退分支，收敛到真整数最优 2。
    /// 两条路径都必须诚实给出 Optimal —— 关键回归点：复核失败时绝不能因
    /// 剪枝而误报 Infeasible。
    #[test]
    fn near_integer_rounding_edge_integration() {
        let mut p = Problem::new("near2", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0e6)], RelOp::Le, 2_999_999.5)
            .unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(
            r.status,
            SolveStatus::Optimal,
            "伪整数复核失败必须回退分支，不得误报 Infeasible"
        );
        let obj = r.best_objective.unwrap();
        assert!(
            obj == 2.0 || obj == 3.0,
            "阶梯语义两种结局：回退分支得 2；复核通过得 3。实际 {obj}"
        );
    }

    // =====================================================================
    // 受控入口（SearchControl / SearchHooks / gap / 进度行 / 并行委托）
    // =====================================================================

    /// SearchControl 默认值逐字段符合契约（= 现版本行为，零变化）。
    #[test]
    fn search_control_default_contract() {
        let c = SearchControl::default();
        assert_eq!(c.threads, 1, "默认 1 = 纯顺序路径");
        assert_eq!(c.mipgap_percent, None);
        assert!(!c.adaptive_heuristics, "默认关闭（现行为：仅根启发式）");
        assert_eq!(c.log_period_secs, None);
    }

    /// 受控默认（SearchControl::default + hooks=None）与既有
    /// branch_and_bound 结果全等（含 nodes 计数与 bound/gap_terminated 口径）。
    #[test]
    fn controlled_default_matches_branch_and_bound() {
        // 纯 LP。
        let mut lp = Problem::new("lp", Sense::Maximize);
        lp.add_var("x", VarType::Continuous).unwrap();
        lp.add_var("y", VarType::Continuous).unwrap();
        lp.set_var_bounds(0, 0.0, 3.0).unwrap();
        lp.set_var_bounds(1, 0.0, 3.0).unwrap();
        lp.set_objective(vec![(0, 3.0), (1, 2.0)], 0.0).unwrap();
        lp.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 4.0)
            .unwrap();
        lp.add_constraint(None, vec![(0, 1.0), (1, 3.0)], RelOp::Le, 6.0)
            .unwrap();
        // 整数两变量 + 极小 node_limit（LimitReached 路径）。
        let mut lim = Problem::new("lim", Sense::Maximize);
        lim.add_var("x", VarType::Integer).unwrap();
        lim.add_var("y", VarType::Integer).unwrap();
        lim.set_objective(vec![(0, 1.0), (1, 1.0)], 0.0).unwrap();
        lim.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 3.5)
            .unwrap();
        // 根不可行。
        let mut inf = Problem::new("inf", Sense::Maximize);
        inf.add_var("x", VarType::Integer).unwrap();
        inf.set_var_bounds(0, 0.0, 10.0).unwrap();
        inf.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        inf.add_constraint(None, vec![(0, 1.0)], RelOp::Ge, 2.5).unwrap();
        inf.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 2.0).unwrap();

        let cases: Vec<Problem> = vec![lp, lim, inf];
        for p in &cases {
            for opts in [plain_opts(), BbOptions::default()] {
                let a = branch_and_bound(p, &opts, tol()).unwrap();
                let b = branch_and_bound_controlled(
                    p,
                    &opts,
                    &SearchControl::default(),
                    tol(),
                    None,
                )
                .unwrap();
                assert_eq!(a.status, b.status, "{:?} {:?}", p.name, opts.root_cuts);
                assert_eq!(a.best_objective, b.best_objective);
                assert_eq!(a.best_values, b.best_values);
                assert_eq!(a.stats.nodes, b.stats.nodes, "nodes 计数必须全等");
                assert_eq!(a.stats.lp_iterations, b.stats.lp_iterations);
                assert_eq!(a.stats.cuts, b.stats.cuts);
                assert_eq!(a.bound, b.bound);
                assert_eq!(a.gap_terminated, b.gap_terminated);
                assert!(!b.gap_terminated, "默认路径绝不可 gap 终止");
            }
        }
        // LimitReached 子句：node_limit=2 时两路径同为 LimitReached 且 nodes=2。
        let opts = BbOptions {
            node_limit: 2,
            ..plain_opts()
        };
        let a = branch_and_bound(&cases[1], &opts, tol()).unwrap();
        let b = branch_and_bound_controlled(
            &cases[1],
            &opts,
            &SearchControl::default(),
            tol(),
            None,
        )
        .unwrap();
        assert_eq!(a.status, SolveStatus::LimitReached);
        assert_eq!(a.stats.nodes, b.stats.nodes);
        assert_eq!(a.stats.nodes, 2);
    }

    /// dive 问题：max 10x1+x2 s.t. 20x1+3x2≤43，整数 [0,10]。
    /// 根 LP 21.5（x1=2.15 分数），根启发式 dive 得 (2,1)=21，整数最优 21。
    fn dive_problem() -> Problem {
        let mut p = Problem::new("dive-gap", Sense::Maximize);
        p.add_var("x1", VarType::Integer).unwrap();
        p.add_var("x2", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_var_bounds(1, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 20.0), (1, 3.0)], RelOp::Le, 43.0)
            .unwrap();
        p
    }

    /// dive 问题的受控选项：关根割（割可能改变根 LP 界，破坏用例确定性），
    /// 保留根启发式（dive 供应 incumbent=21）。
    fn dive_opts() -> BbOptions {
        BbOptions {
            root_cuts: false,
            heuristics: true,
            ..BbOptions::default()
        }
    }

    /// 顺序确定性 gap 用例（dive 问题）：
    /// - mipgap=Some(3.0)：主循环首轮（栈=[根的两个子节点，est=21.5]）即
    ///   gap=(21.5−21)/21×100≈2.381 ≤ 3 ⇒ LimitReached + gap_terminated，
    ///   best=21，bound=触发时刻证据 21.5；
    /// - mipgap=Some(2.0)：2.381 > 2 不触发，全树 3 节点收敛 Optimal 21。
    #[test]
    fn sequential_gap_termination_dive_problem() {
        let p = dive_problem();
        // Some(3.0) ⇒ gap 终止。
        let ctrl = SearchControl {
            mipgap_percent: Some(3.0),
            ..SearchControl::default()
        };
        let r = branch_and_bound_controlled(&p, &dive_opts(), &ctrl, tol(), None).unwrap();
        assert_eq!(r.status, SolveStatus::LimitReached, "gap 终止诚实标注截断");
        assert!(r.gap_terminated);
        assert_close(r.best_objective.unwrap(), 21.0);
        assert_close(r.bound.expect("gap 终止必有证据 bound"), 21.5);
        assert_eq!(r.stats.nodes, 1, "首轮即触发：仅根节点");
        // gap_percent 以触发时刻证据重算：Some(≈2.38) 且非负有限。
        let g = crate::search_ctrl::gap_percent(
            Sense::Maximize,
            r.best_objective.unwrap(),
            r.bound.unwrap(),
        )
        .expect("gap_terminated ⇒ gap_percent 必为 Some");
        assert!(g.is_finite() && g >= 0.0);
        assert!((g - 2.3809523809523813).abs() < 1e-9, "实际 {g}");

        // Some(2.0) ⇒ 不触发，完备最优。
        let ctrl = SearchControl {
            mipgap_percent: Some(2.0),
            ..SearchControl::default()
        };
        let r = branch_and_bound_controlled(&p, &dive_opts(), &ctrl, tol(), None).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert!(!r.gap_terminated);
        assert_close(r.best_objective.unwrap(), 21.0);
        assert_close(r.bound.unwrap(), 21.0);
        assert_eq!(r.stats.nodes, 3);
    }

    /// gap 负值防护（N1）：gap_percent 一律 clamp ≥ 0——db ≤ inc（incumbent
    /// 跳跃越过陈旧 est，或 est=39.999999999999996 < inc=40.0 的浮点噪声）时
    /// gap = 0 而非负值（负 gap 会击穿 'gap=\d+\.\d+%' 正则与
    /// 『gap_reached 行永不矛盾』契约）。单元级直接调用覆盖为主。
    #[test]
    fn gap_percent_clamps_negative_to_zero() {
        use crate::search_ctrl::gap_percent;
        // Max：db ≤ inc ⇒ Some(0.0)。
        assert_eq!(
            gap_percent(Sense::Maximize, 40.0, 39.999999999999996),
            Some(0.0),
            "浮点噪声方向的 db<inc 必须 clamp 为 0"
        );
        assert_eq!(gap_percent(Sense::Maximize, 21.0, 21.0), Some(0.0));
        // 微小负 raw（约 -2.5e-12）clamp 后恰为 0（不会打印 -0.00）。
        assert_eq!(
            gap_percent(Sense::Maximize, 40.0, 40.0 - 1.0e-11),
            Some(0.0)
        );
        // 正常正值：Max=(bound−inc)/|inc|×100；Min 镜像=(inc−bound)/|inc|×100。
        let g = gap_percent(Sense::Maximize, 21.0, 21.5).unwrap();
        assert!((g - 50.0 / 21.0).abs() < 1e-9, "Max 方向 (21.5-21)/21×100，实际 {g}");
        let g = gap_percent(Sense::Minimize, 30.0, 27.0).unwrap();
        assert!((g - 10.0).abs() < 1e-9, "Min 镜像 (30-27)/30×100，实际 {g}");
        // Min 方向 db ≥ inc（无可改进）⇒ 0。
        assert_eq!(gap_percent(Sense::Minimize, 27.0, 30.0), Some(0.0));
        // 非有限输入（缺 incumbent/bound 的表达方式）⇒ None（inf 语义）。
        assert_eq!(gap_percent(Sense::Maximize, f64::NAN, 21.5), None);
        assert_eq!(gap_percent(Sense::Maximize, 21.0, f64::NAN), None);
    }

    /// API 级顺序侧 clamp 场景：大阈值 mipgap 触发时 gap_percent 必为非负
    /// 有限 Some（触发即有 incumbent + 证据 bound，统一 clamp 保证不 inf/负）。
    #[test]
    fn sequential_gap_trigger_yields_nonnegative_finite_gap() {
        // binary_knapsack 两层问题：root LP 15 分数；搜索中 incumbent=10、
        // 剩余栈 est=14.666… ⇒ gap=(14.667-10)/10×100≈46.67% ≤ 50 触发。
        let mut p = Problem::new("k2", Sense::Maximize);
        p.add_var("a", VarType::Binary).unwrap();
        p.add_var("b", VarType::Binary).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 7.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 5.0), (1, 3.0)], RelOp::Le, 7.0)
            .unwrap();
        let ctrl = SearchControl {
            mipgap_percent: Some(50.0),
            ..SearchControl::default()
        };
        let r = branch_and_bound_controlled(&p, &plain_opts(), &ctrl, tol(), None).unwrap();
        assert_eq!(r.status, SolveStatus::LimitReached);
        assert!(r.gap_terminated);
        let g = crate::search_ctrl::gap_percent(
            Sense::Maximize,
            r.best_objective.unwrap(),
            r.bound.unwrap(),
        )
        .expect("gap_terminated ⇒ Some");
        assert!(g.is_finite() && g >= 0.0, "gap 必须非负有限，实际 {g}");
        assert!((g - 14.0 / 3.0 * 10.0).abs() < 1e-6, "实际 {g}");
    }

    /// mipgap=Some(0.0) ≡ None（结果全等；Some(0.0)=关闭）。
    #[test]
    fn controlled_mipgap_zero_equals_none() {
        let p = dive_problem();
        let none_r =
            branch_and_bound_controlled(&p, &dive_opts(), &SearchControl::default(), tol(), None)
                .unwrap();
        let zero_r = branch_and_bound_controlled(
            &p,
            &dive_opts(),
            &SearchControl {
                mipgap_percent: Some(0.0),
                ..SearchControl::default()
            },
            tol(),
            None,
        )
        .unwrap();
        assert_eq!(none_r.status, zero_r.status);
        assert_eq!(none_r.best_objective, zero_r.best_objective);
        assert_eq!(none_r.best_values, zero_r.best_values);
        assert_eq!(none_r.stats.nodes, zero_r.stats.nodes);
        assert_eq!(none_r.stats.lp_iterations, zero_r.stats.lp_iterations);
        assert_eq!(none_r.bound, zero_r.bound);
        assert!(!zero_r.gap_terminated, "Some(0.0) 不得触发 gap 终止");
    }

    /// bound 字段口径：完备 Optimal ⇒ best_objective；Infeasible ⇒ None；
    /// 截断（node_limit）⇒ 剩余栈 est 聚合（此处 = 根 LP 21.5）。
    #[test]
    fn controlled_bound_accounting() {
        // 完备 Optimal。
        let p = dive_problem();
        let r =
            branch_and_bound_controlled(&p, &dive_opts(), &SearchControl::default(), tol(), None)
                .unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_eq!(r.bound, Some(21.0));
        // Infeasible ⇒ None。
        let mut q = Problem::new("band", Sense::Maximize);
        q.add_var("x", VarType::Integer).unwrap();
        q.set_var_bounds(0, 0.4, 0.6).unwrap();
        q.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        let r = branch_and_bound_controlled(&q, &plain_opts(), &SearchControl::default(), tol(), None)
            .unwrap();
        assert_eq!(r.status, SolveStatus::Infeasible);
        assert_eq!(r.bound, None);
        // 截断：node_limit=1 ⇒ 根后栈非空即断 ⇒ bound=剩余聚合=根 est 21.5。
        let opts = BbOptions {
            node_limit: 1,
            ..dive_opts()
        };
        let r = branch_and_bound_controlled(
            &p,
            &opts,
            &SearchControl::default(),
            tol(),
            None,
        )
        .unwrap();
        assert_eq!(r.status, SolveStatus::LimitReached);
        assert!(!r.gap_terminated);
        assert_close(r.bound.unwrap(), 21.5);
        assert_close(r.best_objective.unwrap(), 21.0);
    }

    // —— hard-knapsack-12（N3 合成算例：并发强度与进度行载体）——

    /// hard-knapsack-12：12 变量 0-1 背包（权重/价值/容量见常量）。
    fn hard_knapsack_12() -> Problem {
        const W: [f64; 12] =
            [7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0, 31.0, 37.0, 41.0, 43.0, 47.0];
        const V: [f64; 12] =
            [11.0, 20.0, 32.0, 44.0, 55.0, 70.0, 88.0, 99.0, 111.0, 122.0, 133.0, 144.0];
        let mut p = Problem::new("hard-knapsack-12", Sense::Maximize);
        for j in 0..12 {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        p.set_objective((0..12).map(|j| (j, V[j])).collect(), 0.0)
            .unwrap();
        p.add_constraint(None, (0..12).map(|j| (j, W[j])).collect(), RelOp::Le, 172.0)
            .unwrap();
        p
    }

    /// 2^12 穷举背包最优（真值基准）。
    fn brute_force_knapsack(w: &[f64], v: &[f64], cap: f64) -> f64 {
        let n = w.len();
        assert_eq!(v.len(), n);
        let mut best = 0.0;
        for mask in 0u32..(1u32 << n) {
            let (mut ww, mut vv) = (0.0, 0.0);
            for j in 0..n {
                if mask & (1 << j) != 0 {
                    ww += w[j];
                    vv += v[j];
                }
            }
            if ww <= cap + 1e-9 && vv > best {
                best = vv;
            }
        }
        best
    }

    /// 顺序默认路径与 2^12 穷举一致（最优 523）。
    #[test]
    fn hard_knapsack_12_sequential_matches_bruteforce() {
        const W: [f64; 12] =
            [7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0, 31.0, 37.0, 41.0, 43.0, 47.0];
        const V: [f64; 12] =
            [11.0, 20.0, 32.0, 44.0, 55.0, 70.0, 88.0, 99.0, 111.0, 122.0, 133.0, 144.0];
        assert_eq!(brute_force_knapsack(&W, &V, 172.0), 523.0, "穷举真值须为 523");
        let p = hard_knapsack_12();
        let r = branch_and_bound(&p, &BbOptions::default(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_close(r.best_objective.unwrap(), 523.0);
    }

    // —— 进度行（首行契约 + 周期行仅栈非空）——

    /// 进度行快照（从 ProgressInfo 提取的纯数据，避免依赖 ProgressInfo: Clone）。
    #[derive(Debug, Clone, Copy)]
    struct Line {
        nodes: u64,
        incumbent: Option<f64>,
        bound: Option<f64>,
        gap: Option<f64>,
    }

    /// 收集器钩子：回调把每行快照压入共享 Vec（回调须 Send+Sync）。gap 由
    /// ProgressInfo 的 incumbent/bound 经 search_ctrl::gap_percent 现算（与
    /// format_progress_line 同口径）。
    fn collector() -> (std::sync::Arc<std::sync::Mutex<Vec<Line>>>, SearchHooks) {
        use std::sync::{Arc, Mutex};
        let lines = Arc::new(Mutex::new(Vec::new()));
        let sink = Arc::clone(&lines);
        let hooks = SearchHooks {
            progress: Arc::new(move |i: &crate::search_ctrl::ProgressInfo| {
                if let Ok(mut v) = sink.lock() {
                    v.push(Line {
                        nodes: i.nodes,
                        incumbent: i.incumbent,
                        bound: i.bound,
                        gap: crate::search_ctrl::gap_percent(
                            i.sense,
                            i.incumbent.unwrap_or(f64::NAN),
                            i.bound.unwrap_or(f64::NAN),
                        ),
                    });
                }
            }),
        };
        (lines, hooks)
    }

    /// hard-knapsack-12 + 极小 log_period：首行在根阶段结束（nodes=1）；
    /// 此后周期行仅在栈非空时出现 ⇒ 全部行 bound 均 Some（无
    /// 'incumbent=X bound=none gap=inf' 矛盾行）；行数 ≥ 2（首行后仍有行）。
    #[test]
    fn progress_lines_first_after_root_and_never_bound_none() {
        let p = hard_knapsack_12();
        let (lines, hooks) = collector();
        let ctrl = SearchControl {
            log_period_secs: Some(1e-9),
            ..SearchControl::default()
        };
        let r = branch_and_bound_controlled(&p, &BbOptions::default(), &ctrl, tol(), Some(&hooks))
            .unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_close(r.best_objective.unwrap(), 523.0);
        let lines = lines.lock().unwrap();
        assert!(!lines.is_empty(), "必须收到首行进度");
        assert_eq!(
            lines[0].nodes, 1,
            "首行在根阶段结束的汇合点发射（nodes=1）"
        );
        assert!(
            lines.len() >= 2,
            "极小周期下首行之后至少一行周期行，实际 {} 行",
            lines.len()
        );
        for (i, l) in lines.iter().enumerate() {
            assert!(l.bound.is_some(), "第 {i} 行（nodes={}）bound 不得为 none", l.nodes);
        }
        // nodes 单调不减；行内容不钉数值。
        for w in lines.windows(2) {
            assert!(w[1].nodes >= w[0].nodes, "nodes 单调不减");
        }
    }

    /// 纯 LP 算例也必须收到首行：恰好一行、nodes=1、bound=incumbent=LP 最优、
    /// gap=Some(0)（完备搜索，clamp 后无 -0.00）。
    #[test]
    fn progress_first_line_pure_lp() {
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
        let (lines, hooks) = collector();
        let ctrl = SearchControl {
            log_period_secs: Some(1e-9),
            ..SearchControl::default()
        };
        let r = branch_and_bound_controlled(&p, &BbOptions::default(), &ctrl, tol(), Some(&hooks))
            .unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        let lines = lines.lock().unwrap();
        assert_eq!(lines.len(), 1, "纯 LP 只有根，恰一行");
        assert_eq!(lines[0].nodes, 1);
        assert_close(lines[0].incumbent.unwrap(), 11.0);
        assert_close(lines[0].bound.unwrap(), 11.0);
        assert_eq!(lines[0].gap, Some(0.0));
    }

    /// 根不可行终局：仍先发射首行再返回（nodes=1；incumbent/bound 均 none）。
    #[test]
    fn progress_first_line_root_infeasible() {
        let mut p = Problem::new("inf", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Ge, 2.5).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 2.0).unwrap();
        let (lines, hooks) = collector();
        let r = branch_and_bound_controlled(
            &p,
            &BbOptions::default(),
            &SearchControl::default(),
            tol(),
            Some(&hooks),
        )
        .unwrap();
        assert_eq!(r.status, SolveStatus::Infeasible);
        let lines = lines.lock().unwrap();
        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0].nodes, 1);
        assert_eq!(lines[0].incumbent, None);
        assert_eq!(lines[0].bound, None);
    }

    // —— 并行委托与自适应（对接 crate::parallel_bb / crate::search_ctrl）——

    /// 5 个原验收算例（已知最优）。
    fn acceptance_problems() -> Vec<(String, Problem)> {
        use std::path::Path;
        let cases = [
            ("lp-prod-3", 13.0),
            ("knapsack-10", 17.0),
            ("assign-4x4", 10.0),
            ("milp-wolsey-2", 40.0),
            ("milp-mixed-5", 29.0),
        ];
        cases
            .iter()
            .map(|&(name, _)| {
                let path = format!("{}/instances/{name}.milp", env!("CARGO_MANIFEST_DIR"));
                let p = crate::instance::load_instance(Path::new(&path))
                    .unwrap_or_else(|e| panic!("算例 {name} 加载失败：{e}"));
                (name.to_string(), p)
            })
            .collect()
    }

    /// threads=4 于 5 个原验收算例：最优值与默认（threads=1）完全一致
    ///（池互斥弹出每节点恰处理一次 + incumbent 单调 + 粘滞 stop ⇒ 不丢最优）。
    #[test]
    fn threads4_matches_default_on_acceptance_instances() {
        for (name, p) in acceptance_problems() {
            let seq = branch_and_bound(&p, &BbOptions::default(), tol())
                .unwrap_or_else(|e| panic!("{name} 顺序求解失败：{e}"));
            let par = branch_and_bound_controlled(
                &p,
                &BbOptions::default(),
                &SearchControl {
                    threads: 4,
                    ..SearchControl::default()
                },
                tol(),
                None,
            )
            .unwrap_or_else(|e| panic!("{name} 并行求解失败：{e}"));
            assert_eq!(par.status, SolveStatus::Optimal, "算例 {name}");
            assert!(
                seq.status.is_optimal() && par.best_objective == seq.best_objective,
                "算例 {name}：并行最优值必须与顺序一致（seq={:?} par={:?}）",
                seq.best_objective,
                par.best_objective
            );
        }
    }

    /// threads=4 于 hard-knapsack-12（最优 523）重复 ≥100 次：每次断言
    /// objective=523 且 Optimal（丢唤醒/死锁/竞态在重复下必现；不依赖计时）。
    #[test]
    fn threads4_hard_knapsack_12_repeated_100() {
        let p = hard_knapsack_12();
        let ctrl = SearchControl {
            threads: 4,
            ..SearchControl::default()
        };
        for i in 0..100 {
            let r = branch_and_bound_controlled(&p, &BbOptions::default(), &ctrl, tol(), None)
                .unwrap_or_else(|e| panic!("第 {i} 次并行运行失败：{e}"));
            assert_eq!(r.status, SolveStatus::Optimal, "第 {i} 次运行终态");
            assert_close(r.best_objective.unwrap(), 523.0);
        }
    }

    /// adaptive_heuristics=true 于 5 个原验收算例：仍 Optimal 且目标不变。
    #[test]
    fn adaptive_heuristics_keeps_acceptance_optima() {
        for (name, p) in acceptance_problems() {
            let base = branch_and_bound(&p, &BbOptions::default(), tol())
                .unwrap_or_else(|e| panic!("{name} 默认求解失败：{e}"));
            let ctrl = SearchControl {
                adaptive_heuristics: true,
                ..SearchControl::default()
            };
            let r = branch_and_bound_controlled(&p, &BbOptions::default(), &ctrl, tol(), None)
                .unwrap_or_else(|e| panic!("{name} 自适应求解失败：{e}"));
            assert_eq!(r.status, SolveStatus::Optimal, "算例 {name}");
            assert_close(r.best_objective.unwrap(), base.best_objective.unwrap());
            assert!(!r.gap_terminated);
        }
    }
}
