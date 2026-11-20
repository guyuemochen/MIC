//! 基础启发式与二元等式的折半枚举。
//!
//! - [`round_and_check`]：整数变量取最近整数并夹回界，经与 integrality 联动
//!   的容差阶梯复核后评估目标。
//! - [`dive`]：fix-and-dive 潜水，逐步固定“最接近整数的非整整数变量”。
//! - [`sub_mip`]（方向 1）：根节点 sub-MIP 启发式——根 LP 解中已取整数值的
//!   整数变量按该值固定（RENS 型受限子问题），对剩余分数变量跑一次封顶的
//!   嵌套分支定界（确定性深限/节点上限）。
//! - [`subset_sum_mitm`]（R5-retry）：子集和 meet-in-the-middle 可行性
//!   启发式——全 Binary + 全 Eq + 整系数类的**精确**可行点发现（n ≤ 42
//!   全变量半枚举；n ≤ 84 窗口半枚举 + LCG 确定性窗口），脱离 LP 信息源
//!   （round/dive/泵均锚定 LP 顶点邻域，在 market-split 类上零产出——
//!   planted 解与 LP 取整点汉明距离 ~15，LP 驱动机制原理上跨不过）。
//!   顺序根搜索还可要求全列匹配检查所有目标，只有通过精确性检查才给出证明。

use std::hash::{DefaultHasher, Hasher};
use std::time::Instant;

use crate::branch_bound::{is_integral, select_branch_var};
use crate::model::{Problem, RelOp, Sense, Tolerances, VarType};
use crate::simplex::{Bounds, LpSolution, LpStatus};
use crate::simplex_sparse::LpEngine;

/// 启发式得到的可行整数解。
#[derive(Debug, Clone)]
pub struct Incumbent {
    /// 目标值（按取整后的值评估，原始 sense + obj_offset 口径）。
    pub objective: f64,
    /// 结构变量取值（整数变量已取整）。
    pub values: Vec<f64>,
}

/// 启发式共享的剩余 LP 工作预算；不用于数学剪枝，不修改引擎的时限或迭代提示。
///
/// 旧特性在resolve边界扣累计迭代，允许超支一次调用；新增
/// root_heuristic_work_limit才改为调用内部的工作步（每次换基/界翻转），
/// 不改变原迭代统计。两种路径都计失败、反向重试和嵌套节点；均不保证
/// 给主搜索留下固定墙钟时间，重构等操作不是可抢断的工作步。
#[derive(Debug, Clone, Copy)]
pub(crate) struct HeuristicWorkBudget {
    pub remaining_work: u64,
    pub lp_calls: u64,
}

/// 在实际调用边界统一记账，避免快照恢复或引擎重建影响跨启发式余额。
fn resolve_with_work_budget(
    simplex: &mut dyn LpEngine,
    bounds: &Bounds,
    budget: &mut Option<&mut HeuristicWorkBudget>,
) -> Option<LpSolution> {
    let Some(work) = budget.as_deref_mut() else {
        return Some(simplex.resolve(bounds));
    };
    if work.remaining_work == 0 || work.lp_calls == 0 {
        return None;
    }
    #[cfg(feature = "root_heuristic_work_limit")]
    let sol = simplex.resolve_with_work_limit(bounds, &mut work.remaining_work);
    #[cfg(not(feature = "root_heuristic_work_limit"))]
    let sol = {
        let before = simplex.iterations();
        let sol = simplex.resolve(bounds);
        work.remaining_work = work.remaining_work.saturating_sub(simplex.iterations().saturating_sub(before));
        sol
    };
    work.lp_calls -= 1;
    Some(sol)
}

/// 简单取整 + 可行性复核。
///
/// 整数变量取 nearest 并夹回 bounds（取整值违反界超过 tol.integrality 则返回
/// `None`），连续变量不动；复核容差与 integrality 联动（容差阶梯一致性）：
/// 约束 violation ≤ 10×tol.feasibility + tol.integrality×‖row‖₁（‖row‖₁ = 行
/// 系数绝对值和），界检查同理放宽（10×tol.feasibility + tol.integrality）；
/// 可行则按取整值评估目标（原始 sense + obj_offset 口径）返回。
///
/// 该公式保证“整数判定阈值 1e-6”的节点解在任意行范数下都能过复核；复核
/// 失败返回 `None`（伪整数解不得作为 incumbent）。
pub fn round_and_check(
    problem: &Problem,
    bounds: &Bounds,
    lp_x: &[f64],
    tol: &Tolerances,
) -> Option<Incumbent> {
    let n = problem.vars.len();
    if lp_x.len() != n || bounds.lb.len() != n || bounds.ub.len() != n {
        return None;
    }
    if lp_x.iter().any(|v| !v.is_finite()) {
        return None;
    }
    let mut x = lp_x.to_vec();
    // 整数变量取 nearest 并夹回界；夹回距离超过 integrality 即放弃。
    for j in 0..n {
        if !matches!(problem.vars[j].var_type, VarType::Binary | VarType::Integer) {
            continue;
        }
        let rounded = x[j].round();
        let clamped = rounded.max(bounds.lb[j]).min(bounds.ub[j]);
        if (rounded - clamped).abs() > tol.integrality {
            return None; // 取整值违反界超过 integrality
        }
        x[j] = clamped;
    }
    // 约束复核（容差阶梯）。
    for c in &problem.constraints {
        let lhs = c.expr.eval(&x);
        if !lhs.is_finite() {
            return None;
        }
        let norm1: f64 = c.expr.terms.iter().map(|&(_, a)| a.abs()).sum();
        let allow = 10.0 * tol.feasibility + tol.integrality * norm1;
        let violation = match c.op {
            RelOp::Le => lhs - c.rhs,
            RelOp::Ge => c.rhs - lhs,
            RelOp::Eq => (lhs - c.rhs).abs(),
        };
        if violation > allow {
            return None;
        }
    }
    // 界复核（同理放宽；夹回后整数列精确在界上，连续列在 LP 容差内）。
    let allow_bounds = 10.0 * tol.feasibility + tol.integrality;
    for j in 0..n {
        if x[j] < bounds.lb[j] - allow_bounds || x[j] > bounds.ub[j] + allow_bounds {
            return None;
        }
    }
    let objective = problem.obj.eval(&x) + problem.obj_offset;
    if !objective.is_finite() {
        return None;
    }
    Some(Incumbent {
        objective,
        values: x,
    })
}

/// 潜水（fix-and-dive）启发式。
///
/// 复制 bounds，迭代固定“最接近整数的非整整数变量”（并列取最小下标）到
/// round 值（夹回界）并 `simplex.resolve(副本)`；Infeasible/Unbounded/
/// IterationLimit 即对当前点做 [`round_and_check`] 复核（可行则作为较弱的
/// incumbent 返回——部分修复的解仍可能取整可行，丢弃是纯浪费）；
/// 修复预算（`max_fixes`）用尽同样退化为取整复核而非直接放弃。得整数解
/// 后经 [`round_and_check`] 复核返回（复核用调用方传入的原始 bounds，
/// 保证 incumbent 全局有效）。
///
/// 预算语义：对偶热修复使每次 resolve ≈ 数次迭代，深潜（预算 ≫ 变量数）
/// 成本可控；调用方（`BbOptions::dive_max_fixes`）按规模给预算。
///
/// 注意：dive 会改变 Simplex 内部基，调用方事后必须用真实 bounds 重新 resolve。
pub fn dive(
    problem: &Problem,
    simplex: &mut dyn LpEngine,
    bounds: &Bounds,
    lp_x: &[f64],
    max_fixes: usize,
    tol: &Tolerances,
) -> Option<Incumbent> {
    dive_with_work_budget(problem, simplex, bounds, lp_x, max_fixes, tol, None)
}

/// 带可选共享工作预算的潜水；None 严格保留旧调用顺序，耗尽时复核最后已知点。
pub(crate) fn dive_with_work_budget(
    problem: &Problem,
    simplex: &mut dyn LpEngine,
    bounds: &Bounds,
    lp_x: &[f64],
    max_fixes: usize,
    tol: &Tolerances,
    mut budget: Option<&mut HeuristicWorkBudget>,
) -> Option<Incumbent> {
    let n = problem.vars.len();
    if lp_x.len() != n || bounds.lb.len() != n || bounds.ub.len() != n {
        return None;
    }
    if lp_x.iter().any(|v| !v.is_finite()) {
        return None;
    }
    let int_vars = problem.integer_vars();
    let mut dive_bounds = bounds.clone();
    let mut cur_x = lp_x.to_vec();
    let mut fixes = 0usize;
    // 每变量至多一次反向回退（ceil 不可行改试 floor / 反之），防振荡。
    let mut flipped = vec![false; n];
    loop {
        // 选最接近整数的非整整数变量（并列取最小下标）。
        let mut pick: Option<usize> = None;
        let mut best_dist = f64::INFINITY;
        for &j in &int_vars {
            let dist = (cur_x[j] - cur_x[j].round()).abs();
            if dist > tol.integrality && dist < best_dist {
                best_dist = dist;
                pick = Some(j);
            }
        }
        let j = match pick {
            None => {
                // 全部整数变量已满足整数性：复核后返回。
                return round_and_check(problem, bounds, &cur_x, tol);
            }
            Some(j) => j,
        };
        if fixes >= max_fixes {
            // 预算用尽：部分修复点退化为取整复核（可行则仍是一个合法
            // incumbent——复核用原始 bounds，全局有效性有契约保证）。
            return round_and_check(problem, bounds, &cur_x, tol);
        }
        // 固定到 round 值（夹回界）。
        let target = cur_x[j]
            .round()
            .max(dive_bounds.lb[j])
            .min(dive_bounds.ub[j]);
        dive_bounds.lb[j] = target;
        dive_bounds.ub[j] = target;
        fixes += 1;
        let Some(sol) = resolve_with_work_budget(simplex, &dive_bounds, &mut budget) else {
            return round_and_check(problem, bounds, &cur_x, tol);
        };
        match sol.status {
            LpStatus::Optimal => cur_x = sol.x,
            LpStatus::Infeasible if !flipped[j] => {
                // 一步反向回退：按 cur_x 在 target 的哪一侧选**对侧**整数
                // 邻域（cur_x < target ⇒ target 是向上取整，回退取
                // target−1；cur_x > target ⇒ 向下取整，回退取 target+1），
                // 落在原始界内才改固定再试（系数潜水的廉价回溯——max 背包
                // 类断点向上取整超容量、向下取整即恢复可行的形态由此跨过）。
                let other = if cur_x[j] < target { target - 1.0 } else { target + 1.0 };
                if other >= bounds.lb[j] - tol.integrality
                    && other <= bounds.ub[j] + tol.integrality
                {
                    flipped[j] = true;
                    dive_bounds.lb[j] = other;
                    dive_bounds.ub[j] = other;
                    let Some(retry) = resolve_with_work_budget(simplex, &dive_bounds, &mut budget) else {
                        return round_and_check(problem, bounds, &cur_x, tol);
                    };
                    match retry.status {
                        LpStatus::Optimal => cur_x = retry.x,
                        _ => {
                            return round_and_check(problem, bounds, &cur_x, tol);
                        }
                    }
                } else {
                    return round_and_check(problem, bounds, &cur_x, tol);
                }
            }
            LpStatus::Unbounded
            | LpStatus::IterationLimit
            | LpStatus::TimeLimit
            | LpStatus::Infeasible => {
                // 修复失败（含 W2 墙钟超限的 TimeLimit：时限本就非确定，
                // 中途弃潜不损正确性——dive 只是启发式供源）：退化为对失败
                // 前点的取整复核。
                return round_and_check(problem, bounds, &cur_x, tol);
            }
        }
    }
}

// ===========================================================================
// 根节点 sub-MIP 启发式（方向 1：RENS 型受限子问题嵌套分支定界）
// ===========================================================================

/// sub-MIP 嵌套分支定界的深限（分支层数上限）。
///
/// 取值依据（setcover-300x600 割后根上、阈值 = dive 的 162，release 实测，
/// 经带参内核 [`sub_mip_with_caps`] 的临时标定探针取证、数字落档后探针
/// 删除）：深限 3（HiGHS 同款 sub-MIP 启发式深限——`milp/target/hsc_log.txt`
/// 的 Solving report 明示 `Max sub-MIP depth 3`）在 67 个分数列的受限子问题
/// 上**零产出**（≤15 节点探不到整数叶，260ms 白烧）；深限 8 ⇒ 节点预算内
/// 即取得 145（最优 144，HiGHS 同窗 0.9s 得 144），且比深限放开（d=32/∞）
/// **同质量更快**（3.3s vs 6.2-7.0s——受限 setcover 的深链 LP 界爬升慢，
/// 深潜是白烧）。8 同时与仓库既有浅层深限惯例对齐（Pseudocost 探测 / W6
/// 局部割同为 8）。深限 8 下二叉子树至多 2^9−1 = 511 节点，节点上限是
/// 纯防御性兜底。
const SUB_MIP_MAX_DEPTH: u32 = 8;

/// sub-MIP 嵌套分支定界的节点上限（每次 resolve 计 1，与
/// [`crate::branch_bound::BbOptions::node_limit`] 同记账口径）。
///
/// 深限 8 已隐含 ≤511 节点，本上限在当前深限下永不触发，是深限未来调整
/// 时的防御性兜底（防嵌套搜索成本失控）；两上限合起来保证启发式预算恒为
/// 确定性常数——无墙钟读点、无随机数（成本随算例规模线性缩放：实测
/// setcover-300x600 全预算 ≈ 3.3s、setcover-100x150 在 dive 已命中最优的
/// 阈值剪枝下 ≈ 6ms 自限）。
const SUB_MIP_NODE_LIMIT: u64 = 1_000;

/// 根节点 sub-MIP 启发式（RENS：Relaxation Enforced Neighborhood Search）。
///
/// # 时机（调用方契约）
/// 根 LP 与根割循环求解完成、廉价启发式（[`round_and_check`] / [`dive`]）
/// 之后、MITM 与 DFS 主循环之前——即“根 LP（及根割循环）求解完成后”的
/// 根阶段；仅 `--heuristic-mode` 非 off 档（`BbOptions::heuristics = true`，
/// auto/on 均为 true）运行，off 时调用方整块跳过。
///
/// # 机制
/// 以根 LP 解 x* 为锚构造受限子问题：整数性已满足
///（|x*_j − round(x*_j)| ≤ tol.integrality）的整数变量按该整数值**固定**
///（lb = ub = round 值夹回原界），其余分数整数变量保持自由；对受限子问题
/// 跑一次封顶的嵌套 DFS 分支定界——分支变量选择与主循环同规则
///（most-fractional、并列最小下标破键，[`select_branch_var`]），下行先
/// 探索；确定性上限 = 深限 [`SUB_MIP_MAX_DEPTH`] + 节点上限
/// [`SUB_MIP_NODE_LIMIT`]。
///
/// # incumbent 纪律（诚实性红线）
/// 每个嵌套节点的候选解一律经 [`round_and_check`] 复核——整数变量取整、
/// 逐行可行性核对、目标按取整后的 x 重算——复核通过且**严格优于**阈值
///（外部传入的当前 incumbent 目标，无则任一可行解）才作为候选返回；采纳
/// 与否由调用方按严格更优判定（`maybe_update`）。本启发式只供应
/// incumbent，**绝不影响状态判定**：Optimal/limit 的语义、退出码、目标
/// 输出契约全部不动。
///
/// # 引擎纪律（对齐 [`dive`]）
/// 会扰动 Simplex 内部基（嵌套 resolve 热启动），调用方事后必须用真实根界
/// resolve 恢复根热状态（记账性，不计节点；嵌套迭代数如实累计在引擎
/// 统计里）。根 LP 已全整时立即返回 `None`（不 resolve、零扰动）。传入
/// 引擎可以是稠密 `Simplex` 或稀疏 revised 引擎（`&mut dyn LpEngine`），
/// 两引擎共用同一实现。
pub fn sub_mip(
    problem: &Problem,
    simplex: &mut dyn LpEngine,
    bounds: &Bounds,
    root_lp_x: &[f64],
    incumbent_obj: Option<f64>,
    tol: &Tolerances,
) -> Option<Incumbent> {
    sub_mip_with_work_budget(problem, simplex, bounds, root_lp_x, incumbent_obj, tol, None)
}

/// 带可选共享工作预算的受限搜索；所有嵌套节点共用同一余额。
pub(crate) fn sub_mip_with_work_budget(
    problem: &Problem,
    simplex: &mut dyn LpEngine,
    bounds: &Bounds,
    root_lp_x: &[f64],
    incumbent_obj: Option<f64>,
    tol: &Tolerances,
    budget: Option<&mut HeuristicWorkBudget>,
) -> Option<Incumbent> {
    let n = problem.vars.len();
    if root_lp_x.len() != n || bounds.lb.len() != n || bounds.ub.len() != n {
        return None;
    }
    if root_lp_x.iter().any(|v| !v.is_finite()) {
        return None;
    }
    let int_vars = problem.integer_vars();
    // 受限子问题（RENS 邻域）：整数性已满足的整数变量固定到 round 值
    //（夹回原界），分数整数变量保持自由。
    let mut sub = bounds.clone();
    let mut fractional = 0usize;
    for &j in &int_vars {
        let v = root_lp_x[j];
        if (v - v.round()).abs() > tol.integrality {
            fractional += 1;
        } else {
            let fixed = v.round().max(bounds.lb[j]).min(bounds.ub[j]);
            sub.lb[j] = fixed;
            sub.ub[j] = fixed;
        }
    }
    if fractional == 0 {
        return None; // 根 LP 已全整：无受限邻域可搜（主流程按 (e) 记录）
    }
    if budget.is_none() {
        return sub_mip_with_caps(problem, simplex, bounds, &sub, incumbent_obj, tol,
            SUB_MIP_MAX_DEPTH, SUB_MIP_NODE_LIMIT);
    }
    sub_mip_with_caps_and_budget(
        problem,
        simplex,
        bounds,
        &sub,
        incumbent_obj,
        tol,
        SUB_MIP_MAX_DEPTH,
        SUB_MIP_NODE_LIMIT,
        budget,
        Some(root_lp_x),
    )
}

/// [`sub_mip`] 的带参内核（深限/节点上限显式传入）：生产入口 [`sub_mip`]
/// 以具名常量调用；单元测试与上限标定取证以候选常量调用，与生产共享
/// 同一实现。`base` = 原始根界（候选复核口径，保证 incumbent 全局有效）；
/// `sub` = 已施加整数变量固定的受限界（DFS 起点）。
#[allow(clippy::too_many_arguments)] // 参数即调用契约（同 branch_or_close 先例）
pub(crate) fn sub_mip_with_caps(
    problem: &Problem,
    simplex: &mut dyn LpEngine,
    base: &Bounds,
    sub: &Bounds,
    incumbent_obj: Option<f64>,
    tol: &Tolerances,
    max_depth: u32,
    node_limit: u64,
) -> Option<Incumbent> {
    sub_mip_with_caps_and_budget(problem, simplex, base, sub, incumbent_obj, tol,
        max_depth, node_limit, None, None)
}

#[allow(clippy::too_many_arguments)]
fn sub_mip_with_caps_and_budget(
    problem: &Problem,
    simplex: &mut dyn LpEngine,
    base: &Bounds,
    sub: &Bounds,
    incumbent_obj: Option<f64>,
    tol: &Tolerances,
    max_depth: u32,
    node_limit: u64,
    budget: Option<&mut HeuristicWorkBudget>,
    fallback_x: Option<&[f64]>,
) -> Option<Incumbent> {
    let int_vars = problem.integer_vars();
    let mut ctx = SubMipCtx {
        problem,
        base,
        int_vars: &int_vars,
        sense: problem.sense,
        tol,
        // 阈值含外部 incumbent 目标：仅接受严格更优候选（防御性滤掉非有限值）。
        threshold: incumbent_obj.filter(|o| o.is_finite()),
        nodes_left: node_limit,
        max_depth,
        best: None,
        budget,
    };
    // 新增预算路径先复核最后已知根点，首次调用耗尽或失败也能保留候选。
    // None 不增加候选或阈值更新，严格保持历史搜索顺序。
    if ctx.budget.is_some() {
        if let Some(x) = fallback_x { ctx.offer(x); }
    }
    sub_mip_dfs(&mut ctx, simplex, sub, 0);
    ctx.best
}

/// 嵌套 DFS 的共享上下文（问题与容差只读引用 + 上限/阈值/预算/候选状态）。
struct SubMipCtx<'a, 'b> {
    problem: &'a Problem,
    /// 原始根界（候选复核用——保证 incumbent 全局有效，与 dive 同契约）。
    base: &'a Bounds,
    int_vars: &'a [usize],
    sense: Sense,
    tol: &'a Tolerances,
    /// 候选采纳阈值（外部 incumbent 目标或已找到候选的目标）：仅严格更优入池。
    threshold: Option<f64>,
    /// 剩余节点预算（每次 resolve 计 1；耗尽即整体停止）。
    nodes_left: u64,
    /// 深限（生产路径 = [`SUB_MIP_MAX_DEPTH`]；测试/取证可覆盖标定）。
    max_depth: u32,
    /// 已找到的最优候选。
    best: Option<Incumbent>,
    budget: Option<&'b mut HeuristicWorkBudget>,
}

impl SubMipCtx<'_, '_> {
    /// 候选提取与采纳：对 LP 点做 [`round_and_check`] 复核（整数变量取整 +
    /// 逐行可行性核对 + 目标按 x 重算），通过且严格优于阈值才更新候选并
    /// 收紧阈值。
    fn offer(&mut self, x: &[f64]) {
        let Some(inc) = round_and_check(self.problem, self.base, x, self.tol) else {
            return;
        };
        let better = match self.threshold {
            None => true,
            Some(t) => match self.sense {
                Sense::Maximize => inc.objective > t,
                Sense::Minimize => inc.objective < t,
            },
        };
        if better {
            self.threshold = Some(inc.objective);
            self.best = Some(inc);
        }
    }

    /// 界剪枝判定：节点 LP 界不严格优于阈值 ⇒ 该子树不可能出更优候选
    ///（LP 界是子树最优值的乐观界；相等即无严格改进空间）。
    fn bound_pruned(&self, bound: f64) -> bool {
        match self.threshold {
            None => false,
            Some(t) => match self.sense {
                Sense::Maximize => bound <= t,
                Sense::Minimize => bound >= t,
            },
        }
    }
}

/// 嵌套 DFS 单节点：resolve → 界剪枝 → 取整复核（候选提取）→ 深限/整数性
/// 封闭 → most-fractional 分支（下行先探索，与主循环 DFS 序一致）。
///
/// 分支界施加用逐子节点克隆（深限封顶下节点数 ≤ 2^(深限+1)−1，克隆成本
/// 可忽略；免主循环的 trail 回溯机制）。非 Optimal resolve（不可行/迭代
/// 上限/内部墙钟超限——后者属既有 W2 时限机制，本启发式不新增墙钟读点）
/// 一律按子树放弃处理；预算耗尽即整体停止（封顶诚实：放弃剩余邻域，不影
/// 响正确性）。
fn sub_mip_dfs(ctx: &mut SubMipCtx, simplex: &mut dyn LpEngine, node_bounds: &Bounds, depth: u32) {
    if ctx.nodes_left == 0 {
        return; // 节点预算耗尽：整体停止（确定性封顶）
    }
    let Some(sol) = resolve_with_work_budget(simplex, node_bounds, &mut ctx.budget) else {
        return; // 已求解节点均已复核，保留现有候选，不将预算当成不可行证据。
    };
    ctx.nodes_left -= 1;
    if sol.status != LpStatus::Optimal {
        return; // 不可行 / 上限截断 / 数值异常：子树放弃
    }
    if ctx.bound_pruned(sol.objective) {
        return;
    }
    ctx.offer(&sol.x);
    if depth >= ctx.max_depth || is_integral(ctx.int_vars, &sol.x, ctx.tol.integrality) {
        return; // 深限叶：邻域边界；整数性满足：无可分支变量（防御同判）
    }
    let Some((var, val)) = select_branch_var(ctx.int_vars, &sol.x) else {
        return; // 无可分支变量（防御）：节点丢弃
    };
    let floor = val.floor();
    // 下行（ub → ⌊x̄⌋）先探索；可应用性判定与主循环 tighten_child 同口径
    //（严格收缩且新区间非空）。
    if floor >= node_bounds.lb[var] && floor < node_bounds.ub[var] {
        let mut child = node_bounds.clone();
        child.ub[var] = floor;
        sub_mip_dfs(ctx, simplex, &child, depth + 1);
    }
    // 上行（lb → ⌊x̄⌋+1）。
    if floor + 1.0 > node_bounds.lb[var] && floor + 1.0 <= node_bounds.ub[var] {
        let mut child = node_bounds.clone();
        child.lb[var] = floor + 1.0;
        sub_mip_dfs(ctx, simplex, &child, depth + 1);
    }
}

// ===========================================================================
// 子集和 MITM 可行性启发式（R5-retry）
// ===========================================================================

/// v 的「2 的幂分母指数」：最小 q ≥ 0 使 v·2^q 为整数（|v| < 2^52 内
/// 精确）；v 非有限/分母含非 2 素因子/超出精度返回 None。
fn pow2_denom(v: f64) -> Option<u32> {
    if !v.is_finite() {
        return None;
    }
    let a = v.abs();
    if a >= 9.0e15 {
        return None; // 超出精确整数范围
    }
    let mut x = a;
    let mut q = 0u32;
    while x.fract() != 0.0 {
        x *= 2.0;
        q += 1;
        if q > 52 {
            return None; // 分母含奇因子（乘 2 永不整）
        }
    }
    Some(q)
}

/// 欧几里得 gcd（非负输入；g=0 时返回另一侧）。
fn gcd_i64(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a.abs()
}

/// 每半至多 2²¹ 个掩码，开放寻址表最多占 16 MiB。
const MITM_WINDOW: usize = 42;
const MITM_MAX_VARS: usize = 2 * MITM_WINDOW;

/// 一个子集的等式行和与目标贡献。每条行占 16 位，按模 2¹²⁸ 加减。
/// 窗口内每行的活动范围宽度小于 65536，目标也必须在其范围内。
/// 因而最低一个不等的行差绝对值在 1..65535，后面的行无法抵消它。
#[derive(Clone, Copy, Default)]
struct SubsetSum {
    key: u128,
    cost: i64,
}

impl SubsetSum {
    fn plus(self, other: Self) -> Self {
        Self {
            key: self.key.wrapping_add(other.key),
            cost: self.cost + other.cost,
        }
    }
}

/// 将半边再分成两组，按掩码查两个小表即可恢复行和及目标。
/// 一个 42 列窗口每个小表至多 2048 项，大表只需保存左半掩码。
struct HalfSums {
    low: Vec<SubsetSum>,
    high: Vec<SubsetSum>,
}

impl HalfSums {
    fn new(columns: &[SubsetSum]) -> Self {
        fn sums(columns: &[SubsetSum]) -> Vec<SubsetSum> {
            let mut table = Vec::with_capacity(1 << columns.len());
            table.push(SubsetSum::default());
            for &column in columns {
                for i in 0..table.len() {
                    table.push(table[i].plus(column));
                }
            }
            table
        }
        let (low, high) = columns.split_at(columns.len() / 2);
        Self {
            low: sums(low),
            high: sums(high),
        }
    }

    fn get(&self, mask: u32) -> SubsetSum {
        let mask = mask as usize;
        self.low[mask % self.low.len()].plus(self.high[mask / self.low.len()])
    }
    // 表至多半满。散列只确定初始位置，命中时重建完整的行和键再比较。
    fn slot(&self, slots: &[u32], key: u128) -> usize {
        let mut hash = DefaultHasher::new();
        hash.write_u128(key);
        let mut i = (hash.finish() as usize) & (slots.len() - 1);
        while slots[i] != u32::MAX && self.get(slots[i]).key != key {
            i = (i + 1) & (slots.len() - 1);
        }
        i
    }
}

fn mitm_expired(deadline: Option<Instant>) -> bool {
    deadline.is_some_and(|end| Instant::now() >= end)
}

fn mitm_match(
    rows: &[Vec<i64>],
    vars: &[usize],
    target: &[i64],
    costs: Option<&[i64]>,
    deadline: Option<Instant>,
) -> Option<(Vec<u8>, bool)> {
    let mut columns: Vec<_> = vars
        .iter()
        .map(|&j| SubsetSum {
            key: 0,
            cost: costs.map_or(0, |c| c[j]),
        })
        .collect();
    let mut target_key = 0u128;
    for (i, row) in rows.iter().enumerate() {
        let (mut min, mut max) = (0, 0);
        for (column, &j) in columns.iter_mut().zip(vars) {
            let a = row[j];
            min += a.min(0);
            max += a.max(0);
            column.key = column.key.wrapping_add((a as u128) << (16 * i));
        }
        if max - min >= 65536 || !(min..=max).contains(&target[i]) {
            return None;
        }
        target_key = target_key.wrapping_add((target[i] as u128) << (16 * i));
    }
    let half = vars.len() / 2;
    let left = HalfSums::new(&columns[..half]);
    let right = HalfSums::new(&columns[half..]);
    let empty = u32::MAX;
    // 左半只存掩码，同键时保留费用最低的一个。
    let mut slots = vec![empty; 2 << half];
    for mask in 0..1u32 << half {
        if mask % 4096 == 0 && mitm_expired(deadline) {
            return None;
        }
        let sum = left.get(mask);
        let pos = left.slot(&slots, sum.key);
        let old = slots[pos];
        if old == empty || sum.cost < left.get(old).cost {
            // 升序枚举，目标相等时保留较小掩码。
            slots[pos] = mask;
        }
    }
    let mut best = None;
    let mut complete = costs.is_some();
    for mask in 0..1u32 << (vars.len() - half) {
        if mask % 4096 == 0 && mitm_expired(deadline) {
            complete = false;
            break;
        }
        let sum = right.get(mask);
        let other = slots[left.slot(&slots, target_key.wrapping_sub(sum.key))];
        if other != empty {
            let assignment = (u64::from(mask) << half) | u64::from(other);
            let candidate = (sum.cost + left.get(other).cost, assignment);
            if best.is_none_or(|old| candidate < old) {
                best = Some(candidate);
            }
            if costs.is_none() {
                break;
            }
        }
    }
    let (_, mask) = best?;
    Some((
        (0..vars.len()).map(|j| ((mask >> j) & 1) as u8).collect(),
        complete,
    ))
}

// 分成四组。每次只存左侧两个分组之和落在指定余数桶的组合，
// 右侧取互补余数。每个完整匹配恰在一个桶内，桶号本身不参与等式判定。
const RESIDUES: usize = 256;
const WIDE_MITM_MAX: usize = 50;

struct ResidueHalf {
    sums: HalfSums,
    high_masks: [Vec<u32>; RESIDUES],
}

impl ResidueHalf {
    fn new(columns: &[SubsetSum]) -> Self {
        let sums = HalfSums::new(columns);
        let mut high_masks: [Vec<u32>; RESIDUES] = std::array::from_fn(|_| Vec::new());
        for (mask, s) in sums.high.iter().enumerate() {
            high_masks[s.key as u8 as usize].push(mask as u32);
        }
        Self { sums, high_masks }
    }

    fn count(&self, residue: u8) -> usize {
        self.sums
            .low
            .iter()
            .map(|s| self.high_masks[residue.wrapping_sub(s.key as u8) as usize].len())
            .sum()
    }

    fn visit(
        &self,
        residue: u8,
        deadline: Option<Instant>,
        mut f: impl FnMut(SubsetSum, u32),
    ) -> bool {
        let width = self.sums.low.len() as u32;
        for (low_mask, &low) in self.sums.low.iter().enumerate() {
            if low_mask % 16 == 0 && mitm_expired(deadline) {
                return false;
            }
            for &hi in &self.high_masks[residue.wrapping_sub(low.key as u8) as usize] {
                f(
                    low.plus(self.sums.high[hi as usize]),
                    hi * width + low_mask as u32,
                );
            }
        }
        true
    }
}

// 若恰选 k 列，活动度必在最小 k 项与最大 k 项的和之间。每个原行和
// 所有原行之和都可排除不可能的 k；每个保留的 k 至少要付最便宜的 k 笔费用。
fn restrict_cardinalities(coeff: &mut [i64], target: i64, allowed: &mut [bool]) {
    coeff.sort_unstable();
    let (mut lo, mut hi) = (0, 0);
    for (k, possible) in allowed.iter_mut().enumerate() {
        if k > 0 {
            lo += coeff[k - 1];
            hi += coeff[coeff.len() - k];
        }
        *possible &= (lo..=hi).contains(&target);
    }
}

/// 在完整二元模型的原尺度上匹配，只有查完所有余数桶才报告证优或不可行。
pub(crate) fn bucketed_equalities(
    problem: &Problem,
    tol: &Tolerances,
    deadline: Option<Instant>,
    accept: &mut dyn FnMut(&[f64]),
) -> crate::error::Result<Option<(crate::model::SolveStatus, Option<f64>)>> {
    use crate::model::SolveStatus;
    let n = problem.vars.len();
    if !(MITM_WINDOW + 1..=WIDE_MITM_MAX).contains(&n)
        || problem.constraints.is_empty()
        || problem.constraints.len() > 8
        || problem.constraints.iter().any(|r| r.op != RelOp::Eq)
        || problem
            .vars
            .iter()
            .any(|v| v.var_type != VarType::Binary || v.lb != 0. || v.ub != 1.)
    {
        return Ok(None);
    }
    let Some(costs) = exact_mitm_costs(problem) else {
        return Ok(None);
    };
    let mut cols: Vec<_> = costs
        .iter()
        .map(|&cost| SubsetSum { key: 0, cost })
        .collect();
    let mut target = 0u128;
    let mut allowed = vec![true; n + 1];
    let mut combined = vec![0; n];
    let mut total_rhs = 0;
    for (i, row) in problem.constraints.iter().enumerate() {
        // 每行都用整数，范围限制保证打包键无歧义，也让原模型浮点求和精确。
        if row.rhs.fract() != 0. || row.rhs.abs() >= 65536. {
            return Ok(None);
        }
        let mut lo = 0;
        let mut hi = 0;
        let mut coeff = vec![0; n];
        for &(j, a) in &row.expr.terms {
            if a.fract() != 0. || a.abs() >= 65536. {
                return Ok(None);
            }
            let v = a as i64;
            lo += v.min(0);
            hi += v.max(0);
            cols[j].key = cols[j].key.wrapping_add((v as u128) << (16 * i));
            coeff[j] = v;
            combined[j] += v;
        }
        let norm = (hi - lo) as f64;
        if norm >= 65536.
            || !(lo..=hi).contains(&(row.rhs as i64))
            || !(10. * tol.feasibility + tol.integrality * norm < 1.)
            || !(crate::solver::FINAL_VIOLATION_REL * (1. + row.rhs.abs() + norm) < 1.)
        {
            return Ok(None);
        }
        let rhs = row.rhs as i64;
        restrict_cardinalities(&mut coeff, rhs, &mut allowed);
        total_rhs += rhs;
        target = target.wrapping_add((rhs as u128) << (16 * i));
    }
    restrict_cardinalities(&mut combined, total_rhs, &mut allowed);
    let mut sorted = costs;
    sorted.sort_unstable();
    let mut sum = 0;
    let floor = allowed
        .iter()
        .enumerate()
        .filter_map(|(k, &possible)| {
            if k > 0 {
                sum += sorted[k - 1];
            }
            possible.then_some(sum)
        })
        .min();
    let Some(floor) = floor else {
        return Ok(Some((SolveStatus::Infeasible, None)));
    };
    let bound = problem.obj_offset
        + if problem.sense == Sense::Minimize {
            floor as f64
        } else {
            -floor as f64
        };
    let mut objective = None;
    let mut invalid = false;
    let status = bucket_join(&cols, target, deadline, |mask| {
        let point: Vec<_> = (0..n).map(|j| ((mask >> j) & 1) as f64).collect();
        if problem
            .constraints
            .iter()
            .any(|r| r.expr.eval(&point) != r.rhs)
        {
            invalid = true;
            return;
        }
        objective = Some(problem.obj_offset + problem.obj.eval(&point));
        accept(&point);
    });
    if invalid {
        return Err(crate::error::MilpError::Numerical(
            "分桶匹配解未通过原等式核验".into(),
        ));
    }
    Ok(status.map(|status| {
        (
            status,
            match status {
                SolveStatus::Optimal => objective,
                SolveStatus::LimitReached => Some(bound),
                _ => None,
            },
        )
    }))
}

fn bucket_join(
    cols: &[SubsetSum],
    target: u128,
    deadline: Option<Instant>,
    mut accept: impl FnMut(u64),
) -> Option<crate::model::SolveStatus> {
    use crate::model::SolveStatus;
    if mitm_expired(deadline) {
        return Some(SolveStatus::LimitReached);
    }
    let half = cols.len() / 2;
    let left = ResidueHalf::new(&cols[..half]);
    let right = ResidueHalf::new(&cols[half..]);
    let largest = (0..RESIDUES).map(|b| left.count(b as u8)).max()?;
    // 极端集中的余数仍遵守内存上限，交回已有窗口搜索。表只存左半掩码。
    if largest > 1 << 20 {
        return None;
    }
    let mut slots = vec![u32::MAX; 2 * largest.next_power_of_two()];
    let mut best = None;
    for b in 0..RESIDUES {
        if mitm_expired(deadline) {
            return Some(SolveStatus::LimitReached);
        }
        slots.fill(u32::MAX);
        if !left.visit(b as u8, deadline, |s, mask| {
            let pos = left.sums.slot(&slots, s.key);
            let old = slots[pos];
            if old == u32::MAX || (s.cost, mask) < (left.sums.get(old).cost, old) {
                slots[pos] = mask;
            }
        }) {
            return Some(SolveStatus::LimitReached);
        }
        let residue = (target as u8).wrapping_sub(b as u8);
        if !right.visit(residue, deadline, |s, mask| {
            let left_mask = slots[left.sums.slot(&slots, target.wrapping_sub(s.key))];
            if left_mask == u32::MAX {
                return;
            }
            let cost = s.cost + left.sums.get(left_mask).cost;
            let assignment = (u64::from(mask) << half) | u64::from(left_mask);
            if best.is_none_or(|old| (cost, assignment) < old) {
                best = Some((cost, assignment));
                accept(assignment);
            }
        }) {
            return Some(SolveStatus::LimitReached);
        }
    }
    Some(if mitm_expired(deadline) {
        SolveStatus::LimitReached
    } else if best.is_some() {
        SolveStatus::Optimal
    } else {
        SolveStatus::Infeasible
    })
}

/// 只在整数目标的每个部分和以及带常数的最终值都能用 f64 精确表示时证优。
fn exact_mitm_costs(problem: &Problem) -> Option<Vec<i64>> {
    const LIMIT: i64 = 1 << 53;
    if !crate::model::objective_is_integral(problem) || problem.obj_offset.abs() >= LIMIT as f64 {
        return None;
    }
    let mut magnitude = problem.obj_offset.abs() as i64;
    let mut costs = vec![0; problem.vars.len()];
    for &(j, value) in &problem.obj.terms {
        let cost = value as i64;
        magnitude += cost.abs();
        if magnitude >= LIMIT {
            return None;
        }
        costs[j] = if problem.sense == Sense::Maximize {
            -cost
        } else {
            cost
        };
    }
    Some(costs)
}

/// 在至多八条二元等式中寻找一个可行解。各行须能精确换算成小整数。
/// 42 列以内检查完整模型，43 至 84 列从根 LP 的取整点出发，固定窗口外的列，
/// 按固定种子重试至多 `windows` 个窗口。返回前用当前界和约束复核候选。
pub fn subset_sum_mitm(
    problem: &Problem,
    bounds: &Bounds,
    lp_x: &[f64],
    windows: usize,
    tol: &Tolerances,
) -> Option<Incumbent> {
    subset_sum_search(problem, bounds, lp_x, windows, tol, false, None).map(|s| s.incumbent)
}

pub(crate) struct SubsetSolution {
    pub incumbent: Incumbent,
    pub proven: bool,
}

/// 全列枚举可在核准目标和等式格距后同时证优，局部窗口只供应可行解。
pub(crate) fn subset_sum_search(
    problem: &Problem,
    bounds: &Bounds,
    lp_x: &[f64],
    windows: usize,
    tol: &Tolerances,
    prove: bool,
    deadline: Option<Instant>,
) -> Option<SubsetSolution> {
    let n = problem.vars.len();
    if windows == 0 || n == 0 || n > MITM_MAX_VARS || mitm_expired(deadline) {
        return None;
    }
    if !problem.vars.iter().all(|v| v.var_type == VarType::Binary) {
        return None;
    }
    if problem.constraints.is_empty() || problem.constraints.len() > 8 {
        return None;
    }
    if !problem.constraints.iter().all(|c| c.op == RelOp::Eq) {
        return None;
    }
    let m = problem.constraints.len();
    // 有理 2 的幂分母容忍 + i16 打包范围：R4 缩放使系数/右端带 2 的幂
    // 分母（如 127/128，f64 精确）——逐行乘最小 2 的幂恢复整系数（正因子
    // 乘法不改等式解集，f64 精确无舍入），非 2 的幂分母或有理不可表示仍拒。
    let mut a_rows: Vec<Vec<i64>> = vec![vec![0i64; n]; m];
    let mut b: Vec<i64> = Vec::with_capacity(m);
    let mut certify =
        prove
            && n <= MITM_WINDOW
            && bounds.lb.len() == n
            && bounds.ub.len() == n
            && problem.vars.iter().enumerate().all(|(j, v)| {
                v.lb == 0.0 && v.ub == 1.0 && bounds.lb[j] == 0.0 && bounds.ub[j] == 1.0
            })
            && [tol.feasibility, tol.integrality, tol.zero]
                .iter()
                .all(|t| t.is_finite() && *t >= 0.0);
    for (i, c) in problem.constraints.iter().enumerate() {
        // 行内全部量（含 rhs）的最小整数化 2 的幂指数。
        let mut power = pow2_denom(c.rhs)?;
        for &(_, a) in &c.expr.terms {
            power = power.max(pow2_denom(a)?);
        }
        let scale = 2f64.powi(power as i32);
        let rb = c.rhs * scale;
        if !rb.is_finite() || rb.fract() != 0.0 || rb.abs() >= 32000.0 {
            return None;
        }
        b.push(rb as i64);
        for &(j, a) in &c.expr.terms {
            if j >= n {
                return None;
            }
            let aj = a * scale;
            if !aj.is_finite() || aj.fract() != 0.0 || aj.abs() >= 32000.0 {
                return None;
            }
            a_rows[i][j] = aj as i64;
        }
        // 同除公因数后每个非零残差至少为 g/scale。验收容差须严格小于该格距。
        let mut g = b[i].abs();
        for &a in &a_rows[i] {
            g = gcd_i64(g, a.abs());
            if g == 1 {
                break;
            }
        }
        if certify {
            let norm: f64 = c.expr.terms.iter().map(|&(_, a)| a.abs()).sum();
            let allowance = 10.0 * tol.feasibility + tol.integrality * norm;
            certify = allowance < g as f64 / scale
                && c.expr.terms.iter().any(|&(_, a)| a.abs() > tol.zero);
        }
        if g > 1 {
            for a in &mut a_rows[i] {
                *a /= g;
            }
            b[i] /= g;
        }
    }
    // 打包安全：半侧 |a| 和上界（半变量数 × max|a|）与 |b| 都须在 i16 内。
    let half = n.div_ceil(2) as i64;
    for i in 0..m {
        let mx = a_rows[i].iter().map(|v| v.abs()).max().unwrap_or(0);
        if mx * half > 32767 || b[i].abs() > 32767 {
            return None;
        }
    }
    let (x01, proven) = if n <= MITM_WINDOW {
        let vars: Vec<_> = (0..n).collect();
        let costs = if certify {
            exact_mitm_costs(problem)
        } else {
            None
        };
        mitm_match(&a_rows, &vars, &b, costs.as_deref(), deadline)?
    } else {
        // 窗口模式：基点 = round(lp_x)（需有限）。
        if lp_x.len() != n || lp_x.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let x0: Vec<u8> = (0..n)
            .map(|j| lp_x[j].round().max(0.0).min(1.0) as u8)
            .collect();
        let wsize = MITM_WINDOW;
        let mut lcg: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut result = None;
        for _ in 0..windows {
            if mitm_expired(deadline) {
                break;
            }
            lcg = lcg
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let mut idx: Vec<usize> = (0..n).collect();
            for k in (1..idx.len()).rev() {
                lcg = lcg
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let r = ((lcg >> 33) as usize) % (k + 1);
                idx.swap(k, r);
            }
            idx.truncate(wsize);
            idx.sort_unstable();
            // 目标 = b − A·x0（窗口外），即窗口内新赋值须补足的量。
            let mut target = vec![0i64; m];
            for i in 0..m {
                let mut outside = 0i64;
                for j in 0..n {
                    if !idx.contains(&j) {
                        outside += a_rows[i][j] * x0[j] as i64;
                    }
                }
                target[i] = b[i] - outside;
            }
            if let Some((xw, _)) = mitm_match(&a_rows, &idx, &target, None, deadline) {
                let mut x = x0.clone();
                for (t, &j) in idx.iter().enumerate() {
                    x[j] = xw[t];
                }
                result = Some(x);
                break;
            }
        }
        (result?, false)
    };
    let x: Vec<f64> = x01.into_iter().map(|v| v as f64).collect();
    let incumbent = round_and_check(problem, bounds, &x, tol)?;
    // 直接调用分支定界也能到这里。公开字段可写，证优前确认没有重复列号等非规范输入。
    let proven = proven && problem.validate().is_ok() && !mitm_expired(deadline);
    Some(SubsetSolution { incumbent, proven })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simplex::Simplex;
    use crate::model::Sense;

    fn tol() -> Tolerances {
        Tolerances::default()
    }

    /// （R5-retry）MITM 生产入口在 marketsplit-4x40（全变量精确）与
    /// 5x50（窗口半枚举）上的零→有转变 + 逐行精确验证 + 同输入确定性。
    ///（5x50 窗口模式 debug 下 ~60s、release ~5s——重负载测试。）
    #[test]
    fn subset_sum_mitm_unlocks_marketsplit() {
        use crate::instance::load_instance;
        use std::path::Path;
        let t = Tolerances::default();
        for (name, expect_some) in [("marketsplit-4x40", true), ("marketsplit-5x50", true)] {
            let p = load_instance(Path::new(&format!("../bench/instances/{name}.milp"))).unwrap();
            let mut s = Simplex::new(&p, t).unwrap();
            let b = Bounds::from_problem(&p);
            let sol = s.solve(&b);
            assert_eq!(sol.status, LpStatus::Optimal);
            let inc1 = subset_sum_mitm(&p, &b, &sol.x, 8, &t);
            if expect_some {
                let inc = inc1.expect("{name} 应找到可行点");
                for c in &p.constraints {
                    let lhs: f64 = c.expr.terms.iter().map(|&(j, a)| a * inc.values[j]).sum();
                    assert!(
                        (lhs - c.rhs).abs() <= 1e-6 * (1.0 + c.rhs.abs()),
                        "{name} 行不满足：{lhs} vs {}",
                        c.rhs
                    );
                }
                // 确定性：同输入两次调用逐字节一致。
                let inc2 = subset_sum_search(&p, &b, &sol.x, 8, &t, name == "marketsplit-5x50", None).unwrap();
                assert!(!inc2.proven, "窗口或仅找可行解的调用都不能证优");
                let inc2 = inc2.incumbent;
                assert_eq!(inc.values, inc2.values);
                assert_eq!(inc.objective, inc2.objective);
            }
        }
    }

    #[test]
    fn subset_sum_join_matches_exhaustive_signed_rows_costs_and_ties() {
        // 左半的掩码 3 和 4 在首行同为 2；同键要分别验证价格更优与价格持平。
        let n = 9;
        let vars: Vec<_> = (0..n).collect();
        let row = vec![1, 1, 2, 0, -1, -1, 0, 2, -2];
        let second = vec![0, 0, 0, -3, 0, -3, -2, 0, 5];
        // 检查负行、最高 16 位进位、非连续窗口以及无解目标。
        let rows = vec![
            row.clone(),
            second.clone(),
            row.clone(),
            second.clone(),
            row.clone(),
            second.clone(),
            row.clone(),
            second.clone(),
        ];
        let window = [0, 1, 2, 3, 4, 5, 7, 8];
        for vars in [&vars[..], &window[..]] {
            for costs in [
                vec![0; n],
                vec![1; n],
                vec![-1; n],
                vec![1, 1, -4, 3, -2, 0, 0, -1, 5],
            ] {
                for a in -3..=4 {
                    for b in [-5, 0, 2, 4] {
                        let target = [a, b, a, b, a, b, a, b];
                        let expected = (0u32..1 << vars.len())
                            .filter_map(|mask| {
                                let chosen = || {
                                    vars.iter()
                                        .enumerate()
                                        .filter(|(bit, _)| mask & (1 << bit) != 0)
                                        .map(|(_, &j)| j)
                                };
                                (chosen().map(|j| row[j]).sum::<i64>() == a
                                    && chosen().map(|j| second[j]).sum::<i64>() == b)
                                    .then(|| (chosen().map(|j| costs[j]).sum::<i64>(), mask))
                            })
                            .min()
                            .map(|(_, mask)| mask);
                        let actual = mitm_match(&rows, vars, &target, Some(&costs), None).map(
                            |(x, complete)| {
                                assert!(complete);
                                x.iter()
                                    .enumerate()
                                    .map(|(j, &v)| u32::from(v) << j)
                                    .sum::<u32>()
                            },
                        );
                        assert_eq!(actual, expected);
                    }
                }
            }
        }
    }

    #[test]
    fn residue_buckets_match_direct_enumeration() {
        use crate::model::SolveStatus;
        for seed in 0..20i64 {
            let cols: Vec<_> = (0..12i64)
                .map(|j| {
                    let a = (3 * j + seed) % 9 - 4;
                    let b = (j * j + seed) % 7 - 3;
                    SubsetSum {
                        key: (a as u128).wrapping_add((b as u128) << 112),
                        cost: (seed * j) % 11 - 5,
                    }
                })
                .collect();
            let target = if seed % 5 == 0 {
                100u128
            } else {
                (seed - 8) as u128
            };
            let mut expected = None;
            for mask in 0..1u64 << cols.len() {
                let sum = cols
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| mask & (1 << j) != 0)
                    .fold(SubsetSum::default(), |s, (_, &c)| s.plus(c));
                if sum.key == target && expected.is_none_or(|old| (sum.cost, mask) < old) {
                    expected = Some((sum.cost, mask));
                }
            }
            let mut last = None;
            let status = bucket_join(&cols, target, None, |mask| last = Some(mask));
            assert_eq!(last, expected.map(|(_, mask)| mask));
            assert_eq!(
                status,
                Some(if expected.is_some() {
                    SolveStatus::Optimal
                } else {
                    SolveStatus::Infeasible
                })
            );
            assert_eq!(
                bucket_join(&cols, target, Some(Instant::now()), |_| panic!(
                    "到期后不应发出候选"
                )),
                Some(SolveStatus::LimitReached)
            );
        }
    }

    #[test]
    fn cardinality_proof_needs_every_original_row_and_its_integer_grid() {
        use crate::model::SolveStatus;
        let mut p = Problem::new("矛盾的选取数", Sense::Maximize);
        for j in 0..43 {
            p.add_var(&format!("v{j}"), VarType::Binary).unwrap();
        }
        let ones: Vec<_> = (0..43).map(|j| (j, 1.0)).collect();
        p.set_objective(vec![(0, -2.0), (1, 7.0)], -2.0).unwrap();
        for rhs in [2.0, 3.0] {
            p.add_constraint(None, ones.clone(), RelOp::Eq, rhs)
                .unwrap();
        }
        let check = |p: &Problem, tol: Tolerances| {
            bucketed_equalities(p, &tol, None, &mut |_| panic!("矛盾模型不应有解")).unwrap()
        };
        assert_eq!(check(&p, tol()), Some((SolveStatus::Infeasible, None)));
        assert!(check(
            &p,
            Tolerances {
                integrality: 0.1,
                ..tol()
            }
        )
        .is_none());
        p.constraints[1].rhs = 2.5;
        assert!(check(&p, tol()).is_none());
    }

    #[test]
    fn subset_sum_proof_requires_free_bounds_distinct_grid_and_exact_objective() {
        let mut p = Problem::new("折半证明门槛", Sense::Maximize);
        for j in 0..3 {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        p.add_constraint(None, vec![(0, 1.0), (1, 1.0), (2, 2.0)], RelOp::Eq, 2.0)
            .unwrap();
        p.set_objective(vec![(0, -3.0), (1, 1.0), (2, 5.0)], -7.0)
            .unwrap();
        let mut bounds = Bounds::from_problem(&p);
        let search =
            |p: &Problem, b: &Bounds| subset_sum_search(p, b, &[], 1, &tol(), true, None).unwrap();
        let result = search(&p, &bounds);
        assert!(result.proven);
        assert_eq!(result.incumbent.objective, -2.0);
        p.sense = Sense::Minimize;
        let result = search(&p, &bounds);
        assert!(result.proven);
        assert_eq!(result.incumbent.objective, -9.0);

        assert!(
            subset_sum_search(&p, &bounds, &[], 1, &tol(), true, Some(Instant::now())).is_none()
        );
        bounds.lb[0] = 1.0;
        assert!(!search(&p, &bounds).proven);
        bounds.lb[0] = 0.0;
        p.vars[1].lb = 1.0;
        assert!(!search(&p, &bounds).proven);
        p.vars[1].lb = 0.0;
        p.obj.terms.push((0, -1.0));
        assert!(!search(&p, &bounds).proven);
        p.obj.terms.pop();
        p.obj_offset = 0.5;
        assert!(!search(&p, &bounds).proven);
        p.obj_offset = ((1_i64 << 53) - 2) as f64;
        assert!(!search(&p, &bounds).proven);
        p.obj_offset = 0.0;
        let scale = 2.0_f64.powi(-40);
        let row = &mut p.constraints[0];
        row.rhs *= scale;
        for (_, a) in &mut row.expr.terms {
            *a *= scale;
        }
        assert!(!search(&p, &bounds).proven);
    }

    /// knapsack LP 解的取整复核：max 9a+5b+3c s.t. 7a+4b+3c ≤ 8（0-1）。
    /// LP 解：a 先翻到上界 1，b 补 1/4 → (1, 0.25, 0)，目标 10.25；
    /// 取整 → (1,0,0) 可行，目标 9。
    #[test]
    fn round_and_check_knapsack_lp_solution() {
        let mut p = Problem::new("knapsack", Sense::Maximize);
        let a = p.add_var("a", VarType::Binary).unwrap();
        let b = p.add_var("b", VarType::Binary).unwrap();
        let c = p.add_var("c", VarType::Binary).unwrap();
        p.set_objective(vec![(a, 9.0), (b, 5.0), (c, 3.0)], 0.0)
            .unwrap();
        p.add_constraint(None, vec![(a, 7.0), (b, 4.0), (c, 3.0)], RelOp::Le, 8.0)
            .unwrap();
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let bnds = Bounds::from_problem(&p);
        let sol = s.solve(&bnds);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.objective - 10.25).abs() < 1e-9, "LP 应为 10.25");
        assert!((sol.x[a] - 1.0).abs() < 1e-9 && (sol.x[b] - 0.25).abs() < 1e-9);
        let inc = round_and_check(&p, &bnds, &sol.x, &t).expect("取整解应可行");
        assert!((inc.objective - 9.0).abs() < 1e-9);
        assert_eq!(inc.values, vec![1.0, 0.0, 0.0]);
    }

    /// 取整后违反约束 / 取整值违反界超过 integrality → None。
    #[test]
    fn round_and_check_rejects_infeasible_or_out_of_bounds() {
        let t = tol();
        // 取整后违反行：a+b ≤ 0.9，LP 解 (0.9, 0) → (1,0) 违反 0.1。
        let mut p = Problem::new("r1", Sense::Maximize);
        let a = p.add_var("a", VarType::Binary).unwrap();
        let b = p.add_var("b", VarType::Binary).unwrap();
        p.set_objective(vec![(a, 1.0), (b, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(a, 1.0), (b, 1.0)], RelOp::Le, 0.9)
            .unwrap();
        let bnds = Bounds::from_problem(&p);
        assert!(round_and_check(&p, &bnds, &[0.9, 0.0], &t).is_none());
        // 取整值违反界：整数变量界 [0.2, 0.4]，LP 解 0.2 → round 0 距界 0.2。
        let mut q = Problem::new("r2", Sense::Minimize);
        let i = q.add_var("i", VarType::Integer).unwrap();
        q.set_var_bounds(i, 0.2, 0.4).unwrap();
        q.set_objective(vec![(i, 1.0)], 0.0).unwrap();
        let bnds2 = Bounds::from_problem(&q);
        assert!(round_and_check(&q, &bnds2, &[0.2], &t).is_none());
    }

    /// 连续变量不动、目标按取整值评估（含 obj_offset）。
    #[test]
    fn round_and_check_keeps_continuous_and_evaluates_objective() {
        let mut p = Problem::new("mixed", Sense::Minimize);
        let i = p.add_var("i", VarType::Integer).unwrap();
        let c = p.add_var("c", VarType::Continuous).unwrap();
        p.set_var_bounds(i, 0.0, 5.0).unwrap();
        p.set_var_bounds(c, 0.0, 10.0).unwrap();
        p.set_objective(vec![(i, 2.0), (c, 3.0)], 0.5).unwrap();
        p.add_constraint(None, vec![(i, 1.0), (c, 1.0)], RelOp::Ge, 3.4)
            .unwrap();
        let bnds = Bounds::from_problem(&p);
        let inc = round_and_check(&p, &bnds, &[2.6, 0.8], &tol()).expect("应可行");
        assert_eq!(inc.values, vec![3.0, 0.8]);
        assert!((inc.objective - (2.0 * 3.0 + 3.0 * 0.8 + 0.5)).abs() < 1e-9);
        assert!((inc.objective - 8.9).abs() < 1e-9);
    }

    /// 维度不符 / 非有限输入 → None（防御，不 panic）。
    #[test]
    fn round_and_check_rejects_bad_dimensions() {
        let mut p = Problem::new("dim", Sense::Minimize);
        let i = p.add_var("i", VarType::Integer).unwrap();
        p.set_objective(vec![(i, 1.0)], 0.0).unwrap();
        let bnds = Bounds::from_problem(&p);
        assert!(round_and_check(&p, &bnds, &[], &tol()).is_none());
        assert!(round_and_check(&p, &bnds, &[f64::NAN], &tol()).is_none());
    }

    /// 潜水找到整数解 (2,1)：max 10x1+x2 s.t. 20x1+3x2 ≤ 43（根 LP (2.15,0)）；
    /// dive 后调用方用真实 bounds 重新 resolve 恢复根 LP 21.5。
    #[test]
    fn dive_finds_incumbent_and_caller_resolves_back() {
        let mut p = Problem::new("dive", Sense::Maximize);
        let x1 = p.add_var("x1", VarType::Integer).unwrap();
        let x2 = p.add_var("x2", VarType::Integer).unwrap();
        p.set_var_bounds(x1, 0.0, 10.0).unwrap();
        p.set_var_bounds(x2, 0.0, 10.0).unwrap();
        p.set_objective(vec![(x1, 10.0), (x2, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(x1, 20.0), (x2, 3.0)], RelOp::Le, 43.0)
            .unwrap();
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let bnds = Bounds::from_problem(&p);
        let sol = s.solve(&bnds);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.objective - 21.5).abs() < 1e-9);
        let inc = dive(&p, &mut s, &bnds, &sol.x, 2, &t).expect("dive 应找到 (2,1)");
        assert!((inc.objective - 21.0).abs() < 1e-9);
        assert_eq!(inc.values, vec![2.0, 1.0]);
        // 契约注意：dive 改变了内部基，调用方必须用真实 bounds 重新 resolve。
        let back = s.resolve(&bnds);
        assert_eq!(back.status, LpStatus::Optimal);
        assert!((back.objective - 21.5).abs() < 1e-9);
    }

    /// wolsey 上潜水固定 x1=2、x2=4 触发不可行 ⇒ 一步反向回退改固定
    /// x2=3 ⇒ 找到可行整数解 (2,3)=34（升级前直接 None）；根 LP 41.25
    /// 仍可恢复（热启动不受潜水副本污染）。
    #[test]
    fn dive_backtracks_on_infeasible_fix_and_root_recovers() {
        let mut p = Problem::new("wolsey", Sense::Maximize);
        p.add_var("x1", VarType::Integer).unwrap();
        p.add_var("x2", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 5.0), (1, 8.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 6.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 5.0), (1, 9.0)], RelOp::Le, 45.0)
            .unwrap();
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let bnds = Bounds::from_problem(&p);
        let sol = s.solve(&bnds);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.objective - 41.25).abs() < 1e-9);
        let inc = dive(&p, &mut s, &bnds, &sol.x, 4, &t).expect("回退后应找到可行解");
        assert!((inc.objective - 34.0).abs() < 1e-9);
        assert_eq!(inc.values, vec![2.0, 3.0]);
        let back = s.resolve(&bnds);
        assert_eq!(back.status, LpStatus::Optimal);
        assert!((back.objective - 41.25).abs() < 1e-9);
    }

    /// 一步反向回退的方向按 cur_x 在 target 的侧别取**对侧**整数邻域
    /// （M8）：round-down（LP 2.2 → 固定 2 不可行）回退到 3，round-up
    ///（LP 2.8 → 固定 3 不可行）回退到 2。旧实现只看 target 相对 lb 的
    /// 位置，round-down 情形会误取 1（更远离 LP 值的一侧）⇒ 回退失败。
    #[test]
    fn dive_backtrack_picks_opposite_side_neighborhood() {
        let t = tol();
        // round-down：min x s.t. x ≥ 2.2（LP 解 2.2，固定 2 不可行、对侧 3
        // 可行；固定 1 不可行——旧方向两次失败后 dive 返回 None）。
        let mut p = Problem::new("bt-down", Sense::Minimize);
        let x = p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(x, 0.0, 5.0).unwrap();
        p.set_objective(vec![(x, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(x, 1.0)], RelOp::Ge, 2.2)
            .unwrap();
        let mut s = Simplex::new(&p, t).unwrap();
        let bnds = Bounds::from_problem(&p);
        let sol = s.solve(&bnds);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.x[0] - 2.2).abs() < 1e-9, "根 LP 应为 2.2");
        let inc = dive(&p, &mut s, &bnds, &sol.x, 2, &t)
            .expect("round-down 回退必须取对侧 3（旧实现取 1 失败）");
        assert!((inc.objective - 3.0).abs() < 1e-9);
        assert_eq!(inc.values, vec![3.0]);

        // round-up：max x s.t. x ≤ 2.8（LP 解 2.8，固定 3 不可行、对侧 2
        // 可行）。
        let mut q = Problem::new("bt-up", Sense::Maximize);
        let y = q.add_var("y", VarType::Integer).unwrap();
        q.set_var_bounds(y, 0.0, 5.0).unwrap();
        q.set_objective(vec![(y, 1.0)], 0.0).unwrap();
        q.add_constraint(None, vec![(y, 1.0)], RelOp::Le, 2.8)
            .unwrap();
        let mut s2 = Simplex::new(&q, t).unwrap();
        let bnds2 = Bounds::from_problem(&q);
        let sol2 = s2.solve(&bnds2);
        assert_eq!(sol2.status, LpStatus::Optimal);
        assert!((sol2.x[0] - 2.8).abs() < 1e-9, "根 LP 应为 2.8");
        let inc = dive(&q, &mut s2, &bnds2, &sol2.x, 2, &t)
            .expect("round-up 回退必须取对侧 2");
        assert!((inc.objective - 2.0).abs() < 1e-9);
        assert_eq!(inc.values, vec![2.0]);
    }

    /// max_fixes = 0：分数起点直接 None；整数起点经复核返回。
    #[test]
    fn dive_zero_fixes_behavior() {
        let t = tol();
        // 分数起点、零预算 → None。
        let mut p = Problem::new("frac", Sense::Maximize);
        p.add_var("x1", VarType::Integer).unwrap();
        p.add_var("x2", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 5.0), (1, 8.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 6.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 5.0), (1, 9.0)], RelOp::Le, 45.0)
            .unwrap();
        let mut s = Simplex::new(&p, t).unwrap();
        let bnds = Bounds::from_problem(&p);
        let sol = s.solve(&bnds);
        assert!(dive(&p, &mut s, &bnds, &sol.x, 0, &t).is_none());
        // 整数起点、零预算 → round_and_check 直接返回。
        let mut q = Problem::new("int", Sense::Maximize);
        let x = q.add_var("x", VarType::Integer).unwrap();
        q.set_var_bounds(x, 0.0, 5.0).unwrap();
        q.set_objective(vec![(x, 1.0)], 0.0).unwrap();
        q.add_constraint(None, vec![(x, 1.0)], RelOp::Le, 3.0)
            .unwrap();
        let mut s2 = Simplex::new(&q, t).unwrap();
        let bnds2 = Bounds::from_problem(&q);
        let sol2 = s2.solve(&bnds2);
        assert_eq!(sol2.status, LpStatus::Optimal);
        let inc = dive(&q, &mut s2, &bnds2, &sol2.x, 0, &t).expect("整数解应直接返回");
        assert!((inc.objective - 3.0).abs() < 1e-9);
        assert_eq!(inc.values, vec![3.0]);
    }

    /// 使用真实 LP 引擎与剩余调用数核对重试和递归，避免以迭代为零误判未调用。
    #[test]
    fn work_budget_counts_dive_retry_and_checks_candidate() {
        let mut p = Problem::new("预算反向重试", Sense::Minimize);
        let x = p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(x, 0.0, 5.0).unwrap();
        p.set_objective(vec![(x, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(x, 1.0)], RelOp::Ge, 2.2).unwrap();
        let b = Bounds::from_problem(&p);
        let t = tol();
        for calls in [0, 1, 2] {
            let mut engine = Simplex::new(&p, t).unwrap();
            let root = engine.solve(&b);
            let before = engine.iterations();
            let mut budget = HeuristicWorkBudget { remaining_work: 100, lp_calls: calls };
            let inc = dive_with_work_budget(&p, &mut engine, &b, &root.x, 3, &t, Some(&mut budget));
            assert_eq!(budget.lp_calls, 0);
            assert_eq!(budget.remaining_work, 100u64.saturating_sub(engine.iterations() - before));
            if calls < 2 {
                assert!(inc.is_none(), "最后已知点取整为2违反原约束，不得返回");
            } else {
                assert_eq!(inc.unwrap().values, vec![3.0]);
            }
        }
    }

    #[test]
    fn work_budget_zero_iterations_never_resolves_and_validates_root() {
        let p = sub_mip_problem();
        let b = Bounds::from_problem(&p);
        let t = tol();
        let mut engine = Simplex::new(&p, t).unwrap();
        let root = engine.solve(&b);
        let before = engine.iterations();
        let mut budget = HeuristicWorkBudget { remaining_work: 0, lp_calls: 10 };
        assert!(dive_with_work_budget(&p, &mut engine, &b, &root.x, 10, &t, Some(&mut budget)).is_none());
        assert!(sub_mip_with_work_budget(&p, &mut engine, &b, &root.x, None, &t, Some(&mut budget)).is_none());
        assert_eq!(budget.lp_calls, 10);
        assert_eq!(engine.iterations(), before);
        // 零预算仍允许经过原模型校验的取整解，而非伪造“未找到”。
        let inc = sub_mip_with_work_budget(&p, &mut engine, &b, &[1.0, 1.2], None, &t, Some(&mut budget)).unwrap();
        assert_eq!(inc.values, vec![1.0, 1.0]);
        assert_eq!(budget.lp_calls, 10);
    }

    #[test]
    fn work_budget_is_shared_by_nested_nodes_and_separate_heuristics() {
        let p = sub_mip_problem();
        let b = Bounds::from_problem(&p);
        let t = tol();
        let mut engine = Simplex::new(&p, t).unwrap();
        let root = engine.solve(&b);
        let mut budget = HeuristicWorkBudget { remaining_work: 100, lp_calls: 2 };
        let inc = sub_mip_with_work_budget(&p, &mut engine, &b, &root.x, None, &t, Some(&mut budget)).unwrap();
        assert_eq!(inc.values, vec![1.0, 1.0]);
        assert_eq!(budget.lp_calls, 0, "根及下行节点分别记账，上行不再调用");
        let before = engine.iterations();
        assert!(dive_with_work_budget(&p, &mut engine, &b, &root.x, 10, &t, Some(&mut budget)).is_none());
        assert_eq!(engine.iterations(), before);
    }

    #[test]
    #[cfg(not(feature = "root_heuristic_work_limit"))]
    fn work_budget_positive_iterations_stop_next_call() {
        let p = sub_mip_problem();
        let b = Bounds::from_problem(&p);
        let mut engine = Simplex::new(&p, tol()).unwrap();
        let mut work = HeuristicWorkBudget { remaining_work: 1, lp_calls: 10 };
        let mut budget = Some(&mut work);
        let before = engine.iterations();
        let sol = resolve_with_work_budget(&mut engine, &b, &mut budget).unwrap();
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!(engine.iterations() > before, "冷引擎确实消费正迭代数");
        assert_eq!(budget.as_ref().unwrap().remaining_work, 0);
        assert_eq!(budget.as_ref().unwrap().lp_calls, 9);
        let after = engine.iterations();
        assert!(resolve_with_work_budget(&mut engine, &b, &mut budget).is_none());
        assert_eq!(engine.iterations(), after);
        assert_eq!(budget.as_ref().unwrap().lp_calls, 9);
    }

    #[test]
    fn work_budget_none_preserves_legacy_results_and_iterations() {
        let p = sub_mip_problem();
        let b = Bounds::from_problem(&p);
        let t = tol();
        for nested in [false, true] {
            let mut old = Simplex::new(&p, t).unwrap();
            let root = old.solve(&b);
            let mut optional = old.clone();
            let (a, z) = if nested {
                (sub_mip(&p, &mut old, &b, &root.x, None, &t),
                 sub_mip_with_work_budget(&p, &mut optional, &b, &root.x, None, &t, None))
            } else {
                (dive(&p, &mut old, &b, &root.x, 10, &t),
                 dive_with_work_budget(&p, &mut optional, &b, &root.x, 10, &t, None))
            };
            assert_eq!(a.as_ref().map(|i| (&i.values, i.objective)), z.as_ref().map(|i| (&i.values, i.objective)));
            assert_eq!(old.iterations(), optional.iterations());
            assert_eq!(old.basis_columns(), optional.basis_columns());
        }
    }

    #[test]
    #[cfg(feature = "root_heuristic_work_limit")]
    fn work_limit_box_scope_snapshot_and_none_path() {
        use crate::simplex_sparse::{new_engine, LpEngineKind};
        let p = sub_mip_problem();
        let b = Bounds::from_problem(&p);
        let t = tol();
        for kind in [LpEngineKind::Tableau, LpEngineKind::Revised] {
            let mut engine = new_engine(&p, t, kind).unwrap();
            let deadline = std::time::Instant::now() + std::time::Duration::from_secs(3600);
            engine.set_time_deadline(Some(deadline));
            engine.set_iter_limit_hint(Some(1234));
            let mut work = HeuristicWorkBudget { remaining_work: 1, lp_calls: 10 };
            let limited = resolve_with_work_budget(&mut engine, &b, &mut Some(&mut work)).unwrap();
            assert_eq!(limited.status, LpStatus::IterationLimit);
            assert!(limited.objective.is_nan());
            assert_eq!(work.remaining_work, 0);
            assert_eq!(work.lp_calls, 9);
            let mut snapshot = engine.snapshot();
            for candidate in [&mut engine, &mut snapshot] {
                assert_eq!(candidate.time_deadline(), Some(deadline));
                let recovered = candidate.resolve(&b);
                assert_eq!(recovered.status, LpStatus::Optimal);
                assert_eq!(recovered.objective, 6.0);
                assert_eq!(candidate.time_deadline(), Some(deadline));
            }
            let root = engine.resolve(&b);
            let mut ordinary = engine.snapshot();
            let mut optional = engine.snapshot();
            let a = dive(&p, &mut ordinary, &b, &root.x, 10, &t);
            let z = dive_with_work_budget(&p, &mut optional, &b, &root.x, 10, &t, None);
            assert_eq!(a.as_ref().map(|v| (v.objective, &v.values)), z.as_ref().map(|v| (v.objective, &v.values)));
            assert_eq!(ordinary.iterations(), optional.iterations());
            assert_eq!(ordinary.basis_columns(), optional.basis_columns());
            assert_eq!(ordinary.time_deadline(), optional.time_deadline());
        }
    }

    #[test]
    #[cfg(feature = "root_heuristic_work_limit")]
    fn work_limit_nested_nodes_share_balance_and_root_restore_is_outside() {
        use crate::simplex_sparse::{new_engine, LpEngineKind};
        let p = sub_mip_problem();
        let b = Bounds::from_problem(&p);
        let t = tol();
        for kind in [LpEngineKind::Tableau, LpEngineKind::Revised] {
            let mut engine = new_engine(&p, t, kind).unwrap();
            assert_eq!(engine.solve(&b).status, LpStatus::Optimal);
            engine.add_row(&[(1, 1.0)], RelOp::Le, 1.75).unwrap();
            let root = engine.resolve(&b);
            assert_eq!(root.objective, 6.0);
            let rows = engine.num_rows();
            let before = engine.iterations();
            let mut work = HeuristicWorkBudget { remaining_work: 1, lp_calls: 5 };
            let inc = sub_mip_with_work_budget(&p, &mut engine, &b, &root.x, None, &t, Some(&mut work));
            assert!(inc.is_none(), "中断子节点的未证最优点不能被当作已验候选");
            assert_eq!(work.remaining_work, 0);
            assert_eq!(work.lp_calls, 3, "根邻域零步调用及第一个子节点分别计费");
            assert_eq!(engine.iterations() - before, 1);
            let calls = work.lp_calls;
            assert!(dive_with_work_budget(&p, &mut engine, &b, &root.x, 10, &t, Some(&mut work)).is_none());
            assert_eq!(work.lp_calls, calls);
            let restored = engine.resolve(&b);
            assert_eq!(restored.status, LpStatus::Optimal);
            assert_eq!(restored.objective, root.objective);
            assert_eq!(engine.num_rows(), rows, "恢复原根界时保留已提交割");
            assert_eq!(work.remaining_work, 0, "根恢复不借用或充值启发式余额");
        }
    }

    #[test]
    #[cfg(feature = "root_heuristic_work_limit")]
    fn work_limit_primal_flips_and_zero_budget_do_not_forge_proofs() {
        use crate::simplex_sparse::{new_engine, LpEngineKind};
        let mut flips = Problem::new("纯界翻转", Sense::Maximize);
        for name in ["x", "y"] {
            let j = flips.add_var(name, VarType::Continuous).unwrap();
            flips.set_var_bounds(j, 0.0, 1.0).unwrap();
        }
        flips.set_objective(vec![(0, 1.0), (1, 1.0)], 0.0).unwrap();
        let mut unbounded = Problem::new("真实无界", Sense::Maximize);
        let x = unbounded.add_var("x", VarType::Continuous).unwrap();
        unbounded.set_objective(vec![(x, 1.0)], 0.0).unwrap();
        let mut infeasible = Problem::new("需Phase1证明不可行", Sense::Minimize);
        let x = infeasible.add_var("x", VarType::Continuous).unwrap();
        infeasible.set_var_bounds(x, 0.0, 1.0).unwrap();
        infeasible.add_constraint(None, vec![(x, 1.0)], RelOp::Ge, 2.0).unwrap();
        for kind in [LpEngineKind::Tableau, LpEngineKind::Revised] {
            let bounds = Bounds::from_problem(&flips);
            let mut engine = new_engine(&flips, tol(), kind).unwrap();
            let mut remaining = 1;
            let stopped = engine.resolve_with_work_limit(&bounds, &mut remaining);
            assert_eq!(stopped.status, LpStatus::IterationLimit);
            assert_eq!(remaining, 0);
            assert_eq!(engine.iterations(), 1);
            assert_eq!(stopped.x, vec![1.0, 0.0]);
            assert_eq!(engine.resolve(&bounds).objective, 2.0);
            for (p, expected) in [(&unbounded, LpStatus::Unbounded), (&infeasible, LpStatus::Infeasible)] {
                let bounds = Bounds::from_problem(p);
                let mut engine = new_engine(p, tol(), kind).unwrap();
                let mut remaining = 0;
                let stopped = engine.resolve_with_work_limit(&bounds, &mut remaining);
                assert_eq!(stopped.status, LpStatus::IterationLimit);
                assert!(stopped.objective.is_nan());
                assert_eq!(engine.iterations(), 0);
                let mut remaining = 10;
                assert_eq!(engine.resolve_with_work_limit(&bounds, &mut remaining).status, expected);
            }
        }
    }

    // —— 根节点 sub-MIP 启发式（方向 1）——

    /// 载体：max 3a + 2b s.t. a + b ≤ 2.5（a Binary、b Integer [0,10]）。
    /// 根 LP (a=1, b=1.5) 目标 6——a 整数性满足（受限子问题固定 1）、b 分数
    /// 自由；受限子问题内深限可达的整数最优点 (1,1)=5（全空间穷举真值亦为
    /// 5：(1,1)=5 > (0,2)=4 > (1,0)=3）。
    fn sub_mip_problem() -> Problem {
        let mut p = Problem::new("sub-mip", Sense::Maximize);
        p.add_var("a", VarType::Binary).unwrap();
        p.add_var("b", VarType::Integer).unwrap();
        p.set_var_bounds(1, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 3.0), (1, 2.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 2.5)
            .unwrap();
        p
    }

    /// （方向 1）命中：受限子问题嵌套 DFS 在深限内找到 (1,1)=5；候选经
    /// **独立**复核（界/整数性/逐行可行性/目标按 x 重算——不复用生产校验
    /// 路径）全部通过；同构造双跑逐字节一致（确定性）。
    #[test]
    fn sub_mip_hits_restricted_optimum_validated_and_deterministic() {
        let p = sub_mip_problem();
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let bnds = Bounds::from_problem(&p);
        let sol = s.solve(&bnds);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.objective - 6.0).abs() < 1e-9, "根 LP 应为 6");
        assert!(
            (sol.x[0] - 1.0).abs() < 1e-9 && (sol.x[1] - 1.5).abs() < 1e-9,
            "根 LP 顶点 (1, 1.5)：a 整数性满足、b 分数"
        );
        let inc = sub_mip(&p, &mut s, &bnds, &sol.x, None, &t)
            .expect("受限子问题深限内应命中 (1,1)=5");
        assert!((inc.objective - 5.0).abs() < 1e-9);
        // 独立有效性复核：界、整数性、逐行可行、目标按 x 重算。
        assert_eq!(inc.values.len(), 2);
        for (j, v) in inc.values.iter().enumerate() {
            assert!(
                *v >= bnds.lb[j] - 1e-9 && *v <= bnds.ub[j] + 1e-9,
                "变量 j={j} 取值 {v} 越界"
            );
            if matches!(p.vars[j].var_type, VarType::Binary | VarType::Integer) {
                assert!((v - v.round()).abs() < 1e-9, "整数变量 j={j} 须取整");
            }
        }
        for c in &p.constraints {
            let lhs = c.expr.eval(&inc.values);
            let ok = match c.op {
                RelOp::Le => lhs <= c.rhs + 1e-9,
                RelOp::Ge => lhs >= c.rhs - 1e-9,
                RelOp::Eq => (lhs - c.rhs).abs() <= 1e-9,
            };
            assert!(ok, "行可行性：lhs={lhs} rhs={}", c.rhs);
        }
        let obj = p.obj.eval(&inc.values) + p.obj_offset;
        assert!((obj - inc.objective).abs() < 1e-9, "目标按 x 重算一致");
        // 确定性：同构造新引擎双跑逐字节一致。
        let mut s2 = Simplex::new(&p, t).unwrap();
        let sol2 = s2.solve(&bnds);
        let inc2 = sub_mip(&p, &mut s2, &bnds, &sol2.x, None, &t)
            .expect("双跑第二次也应命中");
        assert_eq!(inc.values, inc2.values);
        assert_eq!(inc.objective, inc2.objective);
    }

    /// （方向 1）不命中：整数不可行带（x ∈ [2.2, 2.8] 内无整数点）——根 LP
    /// 2.8 分数，下行（⌊2.8⌋=2 < lb）与上行（3 > ub）分支界均空 ⇒ None
    ///（诚实 miss，不得伪造候选）。
    #[test]
    fn sub_mip_misses_on_integer_infeasible_band() {
        let mut p = Problem::new("band", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Ge, 2.2).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 2.8).unwrap();
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let bnds = Bounds::from_problem(&p);
        let sol = s.solve(&bnds);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.x[0] - 2.8).abs() < 1e-9, "根 LP 应取上界 2.8");
        assert!(sub_mip(&p, &mut s, &bnds, &sol.x, None, &t).is_none());
    }

    /// （方向 1）阈值剪枝：外部 incumbent 目标恰为受限子问题最优 5 ⇒ 唯一
    /// 整数叶 (1,1) 的 LP 界 5 不严格更优 ⇒ 剪枝不采纳，返回 None；弱阈值
    ///（4.0）⇒ 5 严格更优 ⇒ 命中（只在优于当前 incumbent 时采纳）。
    #[test]
    fn sub_mip_threshold_prunes_non_improving_candidates() {
        let p = sub_mip_problem();
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let bnds = Bounds::from_problem(&p);
        let sol = s.solve(&bnds);
        assert!(
            sub_mip(&p, &mut s, &bnds, &sol.x, Some(5.0), &t).is_none(),
            "阈值 = 受限最优 ⇒ 无严格更优候选"
        );
        let mut s2 = Simplex::new(&p, t).unwrap();
        let sol2 = s2.solve(&bnds);
        let inc = sub_mip(&p, &mut s2, &bnds, &sol2.x, Some(4.0), &t)
            .expect("弱阈值应命中 5");
        assert!((inc.objective - 5.0).abs() < 1e-9);
    }

    /// （方向 1）根 LP 已全整：无受限邻域可搜 ⇒ 立即 None 且零引擎扰动
    ///（迭代数不变——早退路径不 resolve，不扰动基）。
    #[test]
    fn sub_mip_skips_and_does_not_touch_engine_when_root_integral() {
        let mut p = Problem::new("int-root", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 5.0).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 3.0).unwrap();
        let t = tol();
        let mut s = Simplex::new(&p, t).unwrap();
        let bnds = Bounds::from_problem(&p);
        let sol = s.solve(&bnds);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.x[0] - 3.0).abs() < 1e-9, "根 LP 应全整 (3)");
        let iters = s.iterations();
        assert!(sub_mip(&p, &mut s, &bnds, &sol.x, None, &t).is_none());
        assert_eq!(s.iterations(), iters, "全整早退不 resolve（零扰动）");
    }

    /// （方向 1）确定性上限语义（带参内核 [`sub_mip_with_caps`]）：深限 0
    ///（只解析受限根、不分支）与节点上限 1（恰一次 resolve）在受限根取整
    /// miss 的载体上都 ⇒ None 且引擎只多花一次 resolve 的迭代；深限放开到
    /// 1（下行叶 (1,1) 可达）⇒ 命中 5——两上限各自独立封顶搜索。
    #[test]
    fn sub_mip_caps_bound_the_nested_search() {
        let p = sub_mip_problem();
        let t = tol();
        // 受限界的构造与生产入口同口径（a=1 固定、b 自由）。
        let mk = || {
            let mut s = Simplex::new(&p, t).unwrap();
            let bnds = Bounds::from_problem(&p);
            let sol = s.solve(&bnds);
            assert_eq!(sol.status, LpStatus::Optimal);
            let mut sub = bnds.clone();
            sub.lb[0] = 1.0;
            sub.ub[0] = 1.0;
            (s, bnds, sub)
        };
        // 深限 0：只解析受限根（round miss）⇒ None（受限根固定 a=1 后同顶点
        //（1,1.5）仍最优，resolve 可零主元——迭代数不作断言）。
        {
            let (mut s, bnds, sub) = mk();
            assert!(sub_mip_with_caps(&p, &mut s, &bnds, &sub, None, &t, 0, 1_000).is_none());
        }
        // 节点上限 1：同深限 0 的语义（预算只够根）⇒ None。
        {
            let (mut s, bnds, sub) = mk();
            assert!(sub_mip_with_caps(&p, &mut s, &bnds, &sub, None, &t, 8, 1).is_none());
        }
        // 深限 1：下行叶 (1,1)=5 可达 ⇒ 命中。
        {
            let (mut s, bnds, sub) = mk();
            let inc = sub_mip_with_caps(&p, &mut s, &bnds, &sub, None, &t, 1, 1_000)
                .expect("深限 1 内下行叶应命中 (1,1)=5");
            assert!((inc.objective - 5.0).abs() < 1e-9);
        }
    }
}
