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

use std::time::Instant;

use crate::cuts::{self, CutOptions};
use crate::error::{MilpError, Result};
use crate::heuristics::{self, Incumbent};
use crate::model::{Problem, Sense, SolveStatus, Tolerances};
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
    /// 求解终态。
    pub status: SolveStatus,
    /// 最优目标（原始 sense + obj_offset 口径；LimitReached 时为已有 incumbent）。
    pub best_objective: Option<f64>,
    /// 最优解（与入参 problem 同维同序，未经 postsolve；整数变量已取整）。
    pub best_values: Option<Vec<f64>>,
    /// 统计信息。
    pub stats: BbStats,
}

/// 深度优先分支定界主入口。
pub fn branch_and_bound(problem: &Problem, opts: &BbOptions, tol: Tolerances) -> Result<BbResult> {
    let start = Instant::now();
    let int_vars = problem.integer_vars();
    let sense = problem.sense;
    let mut stats = BbStats::default();
    let mut simplex = Simplex::new(problem, tol)?;
    let mut bounds = Bounds::from_problem(problem);
    // incumbent = (目标, 取值)；目标一律按取整后的值评估（WP3 复核口径）。
    let mut incumbent: Option<(f64, Vec<f64>)> = None;

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
        let incumbent = if sol.status == LpStatus::Optimal {
            Some((sol.objective, sol.x))
        } else {
            None
        };
        return Ok(finish(status, incumbent, &mut stats, &simplex, start));
    }

    // 根 solve 非 Optimal ⇒ 不进割循环与启发式（gomory_cuts 的前提是 Optimal
    // 表格），直接按 (a)/(b)/(c) 终局。
    if let Some(status) = root_terminal(classify_node(sol.status, true))? {
        return Ok(finish(status, None, &mut stats, &simplex, start));
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
            return Ok(finish(status, None, &mut stats, &simplex, start));
        }
        // 取回割后的根 LP 解（记账性 resolve，不计节点）。
        sol = simplex.resolve(&bounds);
        if let Some(status) = root_terminal(classify_node(sol.status, true))? {
            return Ok(finish(status, None, &mut stats, &simplex, start));
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
            return Ok(finish(status, None, &mut stats, &simplex, start));
        }
    }

    // —— DFS 主循环状态 ——
    let mut stack: Vec<Pending> = Vec::new();
    let mut trail: Vec<TrailEntry> = Vec::new();
    let mut limit_hit = false;

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
    );

    loop {
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
        );
    }

    // —— 终态判定（与上限诚实性整合，禁止字面分裂）——
    let status = if limit_hit {
        SolveStatus::LimitReached
    } else if incumbent.is_some() {
        SolveStatus::Optimal
    } else {
        SolveStatus::Infeasible
    };
    Ok(finish(status, incumbent, &mut stats, &simplex, start))
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
/// - 分支：MostFractional 选变量，下行（x ≤ ⌊x̄⌋）后压 ⇒ 先探索。
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
    if integral {
        if let Some(inc) = heuristics::round_and_check(problem, bounds, &sol.x, &tol) {
            maybe_update(problem.sense, incumbent, inc);
            check_passed = true;
        }
    }
    if integer_node_eval(integral, check_passed) != NodeEval::Branch {
        return; // 复核通过 ⇒ 回溯（节点关闭）
    }
    // 分支：MostFractional；无可分支变量为内部不可达防御（不剪枝即丢弃节点）。
    let Some((var, val)) = select_branch_var(int_vars, &sol.x) else {
        return;
    };
    let floor = val.floor();
    let unwind_to = trail.len();
    // 上行先压（后弹出）、下行后压（先探索）。
    stack.push(Pending {
        var,
        tighten_lb: true,
        bound: floor + 1.0,
        unwind_to,
    });
    stack.push(Pending {
        var,
        tighten_lb: false,
        bound: floor,
        unwind_to,
    });
}

/// 严格更优才替换 incumbent（按原始 sense 方向；目标与取值均为 WP3 取整后
/// 口径，防报不可行值）。
fn maybe_update(sense: Sense, incumbent: &mut Option<(f64, Vec<f64>)>, inc: Incumbent) {
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
fn time_exceeded(start: &Instant, limit_ms: u64) -> bool {
    start.elapsed().as_millis() >= limit_ms as u128
}

/// 组装结果并补全统计（耗时与累计迭代）。
fn finish(
    status: SolveStatus,
    incumbent: Option<(f64, Vec<f64>)>,
    stats: &mut BbStats,
    simplex: &Simplex,
    start: Instant,
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
}
