//! 统一求解编排与终检（WP4 契约，本文件为 WP4 实现）。
//!
//! 流程：`validate` →（可选）时间窗路线求解或 presolve → branch_and_bound
//! → postsolve → 安全网终检 → 报告。
//!
//! 受控入口 [`solve_controlled`] 在同一编排上注入 [`crate::branch_bound::SearchControl`]
//!（threads/mipgap/自适应启发式/日志周期）与 [`crate::branch_bound::SearchHooks`]
//!（进度回调），输出 [`ControlledSolveReport`]（常规报告 + gap/对偶界证据）；
//! [`solve`] 签名不变，内部为默认受控路径的薄委托。
//!
//! 安全网终检：对原始 Problem 逐约束复核 violation ≤ 1e-6×(1+|rhs|+‖row‖₁)、
//! 整数性 |x−round| ≤ 1e-5（超差的整数值已替换为 round 值）、目标以最终 values
//! 重算（值与目标自洽）；违反即 `Err(MilpError::Numerical)` 而非输出错解。
//! objective 口径 = 原 sense、含 obj_offset。
//!
//! # 实现说明（WP4）
//! - 开启预处理和启发式的单线程入口会先尝试可识别的时间窗路径覆盖。
//!   每条原行必须能由恰选一、时间窗或活动行程解释，候选要通过原模型复核。
//!   专用求解只有在精确对偶与覆盖搜索共同完成证明时才报最优。
//!   结构不适用时走通用流程，接手前已核验的解与已消耗的预算仍有效。
//! - 单线程预处理也可将独立的二元短行组编成选项。每个选项保存原变量赋值，
//!   改写后仍走通用分支定界，回调与终解都还原到原模型并检查。
//! - 若固定共享决策后各条需求可独立求解，可在根节点用整数成本表直接证优。
//! - 少量整数容量的二元装载可用单行背包上界精确搜索。非负特征分类可从冲突费用
//!   推出固定开关；周期等式可先修出一个原模型可行解，交给后续分支定界。
//! - 可由二层逻辑门完整解释的二元网络只搜索主列。每个可行赋值回填原列后终检。
//! - presolve 关闭时不做约束化简（`values` 仍在原列顺序）；
//!   开启时 `presolve_reductions` = 被固定删除的列数 + 被删除的行数 +
//!   被聚合合并的重复列数（W5）。
//! - 时间窗路线可在缩放前直接返回已核验的原空间解。通用路径在 B&B 带回取值
//!   （Optimal 或有 incumbent 的 LimitReached）时还原并终检，未找到解时
//!   报告 `objective` 与 `values` 均为 `None`。
//! - 终检容差为契约固定安全网值（相对 1e-6、整数 1e-5），宽于 WP3 复核
//!   阶梯，正常路径不会被误伤；超差说明数值已失真，宁可报错也不输出错解。

use std::time::{Duration, Instant};

use crate::branch_bound::{
    self, BbOptions, BbStats, NodeRule, SearchControl, SearchHooks, SearchPrelude,
};
use crate::conditional_rows;
use crate::error::{MilpError, Result};
use crate::heuristics::{self, Incumbent};
use crate::model::{Problem, RelOp, Sense, SolveStatus, Tolerances, VarType};
use crate::presolve;
use crate::search_ctrl::{gap_percent, Determinism};
use crate::simplex::Bounds;
use crate::time_window_routes::{sum_bounds, RouteModel};

/// 终检约束违反相对容差（契约固定：violation ≤ 1e-6×(1+|rhs|+‖row‖₁)）。
pub(crate) const FINAL_VIOLATION_REL: f64 = 1e-6;
/// 终检整数性容差（契约固定：|x−round(x)| ≤ 1e-5，超差即 Err）。
const FINAL_INTEGER_TOL: f64 = 1e-5;

/// 求解选项。
#[derive(Debug, Clone, Copy)]
pub struct SolveOptions {
    /// 数值容差。
    pub tol: Tolerances,
    /// 分支定界参数。
    pub bb: BbOptions,
    /// 是否启用预处理。默认 true。
    pub use_presolve: bool,
    /// 是否启用问题级数值缩放（R4，轴4：行列迭代几何均值，2 的幂系数，
    /// scale→presolve→solve→postsolve→unscale）。默认 true；CLI
    /// `--scaling=off` 关闭供 A/B。目标/解还原原尺度（终检与目标重算在
    /// 原问题上做），输出契约不变。
    pub use_scaling: bool,
}

impl Default for SolveOptions {
    /// 默认值：Tolerances::default() + BbOptions::default() + use_presolve=true
    /// + use_scaling=true。
    fn default() -> Self {
        SolveOptions {
            tol: Tolerances::default(),
            bb: BbOptions::default(),
            use_presolve: true,
            use_scaling: true,
        }
    }
}

/// 求解报告。
#[derive(Debug, Clone)]
pub struct SolveReport {
    /// 求解终态。
    pub status: SolveStatus,
    /// 目标值（原始 sense + obj_offset 口径；非 Optimal/LimitReached-无incumbent
    /// 时为 None）。
    pub objective: Option<f64>,
    /// 变量取值（原始 Problem 变量空间；无可行解时为 None）。
    pub values: Option<Vec<f64>>,
    /// 求解统计。
    pub stats: BbStats,
    /// 预处理化简量（删除的列/行计数）。
    pub presolve_reductions: usize,
}

/// 受控求解报告：在既有 [`SolveReport`] 之上附加 gap/对偶界证据字段。
#[derive(Debug, Clone)]
pub struct ControlledSolveReport {
    /// 常规求解报告（形状与 [`solve`] 的输出完全一致；gap 终止时
    /// status = LimitReached，携带 incumbent 解）。
    pub report: SolveReport,
    /// 是否因达到 mipgap 阈值提前终止（true ⇒ report.status = LimitReached）。
    pub gap_terminated: bool,
    /// 是否因策略级自动降级提前终止（auto 模式 + 时间预算 + 双停滞）。
    /// true ⇒ `gap_terminated` 必为 true——与 mipgap 同径的诚实提前终止，
    /// 此字段区分触发原因（「继续证最优」被自动放弃，返回当前最好解）。
    pub auto_degraded: bool,
    /// 触发时刻的 gap 百分比（经统一 clamp ≥ 0）。gap_terminated = true 时
    /// 必为非负有限的 Some；否则为 None。
    pub gap_percent: Option<f64>,
    /// 对偶界（口径见 [`crate::branch_bound::BbResult::bound`]）。
    pub bound: Option<f64>,
}

/// 受控统一求解入口：在 [`solve`] 的编排（validate → 可选 presolve →
/// 分支定界 → postsolve → 安全网终检）之上，经 [`SearchControl`] 与
/// [`SearchHooks`] 注入线程数 / mipgap / 自适应启发式 / 进度回调。
pub fn solve_controlled(
    problem: &Problem,
    opts: &SolveOptions,
    ctrl: &SearchControl,
    hooks: Option<&SearchHooks>,
) -> Result<ControlledSolveReport> {
    solve_controlled_impl(problem, opts, ctrl, hooks, None)
}

/// 单线程受控求解，并在搜索接受新解时报告接受时刻、原尺度目标和变量值。
/// 回调只读。工作模型的解按原列顺序还原，最终报告仍按原模型终检。
pub fn solve_controlled_with_incumbent(
    problem: &Problem,
    opts: &SolveOptions,
    ctrl: &SearchControl,
    hooks: Option<&SearchHooks>,
    observe: &dyn Fn(Instant, f64, &[f64]),
) -> Result<ControlledSolveReport> {
    if ctrl.threads != 1 {
        return Err(MilpError::InvalidModel("新解观测仅支持单线程".into()));
    }
    solve_controlled_impl(problem, opts, ctrl, hooks, Some(observe))
}

type IncumbentObserver<'a> = Option<&'a dyn Fn(Instant, f64, &[f64])>;

fn solve_controlled_impl(
    problem: &Problem,
    opts: &SolveOptions,
    ctrl: &SearchControl,
    hooks: Option<&SearchHooks>,
    incumbent_observer: IncumbentObserver<'_>,
) -> Result<ControlledSolveReport> {
    let start = Instant::now();
    problem.validate()?;
    let deadline = if ctrl.determinism == Determinism::Budget {
        None
    } else {
        branch_bound::bb_time_deadline(&start, &opts.bb)
    };
    let eligible = opts.use_presolve
        && ctrl.threads == 1
        && !(opts.bb.node_rule == NodeRule::HybridBestBound && opts.bb.local_cuts);
    let form = if eligible {
        if opts.bb.heuristics && !ctrl.mipgap_percent.is_some_and(|gap| gap > 0.0) {
            if let Some(report) = solve_integer_structure(
                problem,
                opts,
                ctrl,
                start,
                deadline,
                hooks,
                incumbent_observer,
            )? {
                return Ok(report);
            }
            if let Some(form) =
                crate::classification_conflicts::reformulate(problem, opts.tol, deadline)
            {
                return solve_classification(
                    problem,
                    form,
                    opts,
                    ctrl,
                    hooks,
                    incumbent_observer,
                    start,
                    deadline,
                );
            }
        }
        let form = conditional_rows::reformulate(problem, opts.tol, deadline);
        if deadline.is_some_and(|end| Instant::now() >= end) {
            let stopped = PresolveOutcome {
                start,
                deadline,
                stats: BbStats::default(),
                best: None,
                bound: None,
                status: Some(SolveStatus::LimitReached),
            };
            return Ok(stopped.report(problem, hooks, 0));
        }
        form
    } else {
        None
    };
    let Some(form) = form else {
        return solve_prepared(problem, opts, ctrl, hooks, incumbent_observer, start, None);
    };
    form.problem.validate()?;
    // 每个选项只代表一个整数组合。原空间和选项空间的目标（含常数）须相同。
    let restore = |point: &[f64], objective: f64| -> Result<Vec<f64>> {
        let x = form
            .restore(point)
            .ok_or_else(|| MilpError::Numerical("二元选项解无法还原为原变量".into()))?;
        if x.len() != problem.vars.len()
            || x.iter()
                .zip(&problem.vars)
                .any(|(&v, var)| !(v == 0.0 || v == 1.0) || v < var.lb || v > var.ub)
            || problem
                .constraints
                .iter()
                .any(|row| !(row.expr.eval(&x) <= row.rhs))
            || problem.obj.eval(&x) + problem.obj_offset != objective
        {
            return Err(MilpError::Numerical("二元选项解未通过原模型复核".into()));
        }
        Ok(x)
    };
    let callback_error = std::cell::Cell::new(None);
    let lift_event = |when, objective, point: &[f64]| {
        if let Some(observe) = incumbent_observer {
            match restore(point, objective) {
                Ok(x) => observe(when, objective, &x),
                Err(error) => callback_error.set(Some(error)),
            }
        }
    };
    let observer = incumbent_observer.map(|_| &lift_event as &dyn Fn(Instant, f64, &[f64]));
    let mut result = solve_prepared(&form.problem, opts, ctrl, hooks, observer, start, None)?;
    if let Some(error) = callback_error.take() {
        return Err(error);
    }
    result.report.values = match (result.report.values.take(), result.report.objective) {
        (Some(point), Some(objective)) => Some(restore(&point, objective)?),
        (None, None) => None,
        _ => return Err(MilpError::Internal("二元选项目标与解向量不对应".into())),
    };
    if deadline.is_some_and(|end| Instant::now() >= end) {
        result.report.status = SolveStatus::LimitReached;
    }
    result.report.stats.time_ms = start.elapsed().as_millis().min(u64::MAX as u128) as u64;
    Ok(result)
}

// 分类模型只保留可能改善已知解的分支。首解、后续回调、终解都在原变量下计价。
fn solve_classification(
    original: &Problem,
    form: crate::classification_conflicts::Reformulation,
    opts: &SolveOptions,
    ctrl: &SearchControl,
    hooks: Option<&SearchHooks>,
    observe: IncumbentObserver<'_>,
    start: Instant,
    deadline: Option<Instant>,
) -> Result<ControlledSolveReport> {
    form.problem.validate()?;
    let checked = |mut values: Vec<f64>| -> Result<Incumbent> {
        if values.len() != original.vars.len()
            || values.iter().zip(&original.vars).any(|(&v, b)| {
                !v.is_finite() || v < b.lb - opts.tol.feasibility || v > b.ub + opts.tol.feasibility
            })
        {
            return Err(MilpError::Numerical("分类解未通过原变量界检查".into()));
        }
        final_check(original, &mut values)?;
        let objective = original.obj.eval(&values) + original.obj_offset;
        Ok(Incumbent { objective, values })
    };
    let restore = |point: &[f64], objective: f64| -> Result<Incumbent> {
        let values = form
            .restore(original, point)
            .ok_or_else(|| MilpError::Numerical("分类模型的误差变量无法还原".into()))?;
        let candidate = checked(values)?;
        if (candidate.objective - objective).abs() > 1e-8 * (1.0 + objective.abs()) {
            return Err(MilpError::Numerical("分类解还原后目标不一致".into()));
        }
        Ok(candidate)
    };
    let seed = checked(form.seed.clone())?;
    if let Some(observe) = observe {
        observe(Instant::now(), seed.objective, &seed.values);
    }
    let best = std::cell::RefCell::new(seed);
    let error = std::cell::Cell::new(None);
    let lift = |when: Instant, objective: f64, point: &[f64]| match restore(point, objective) {
        Ok(candidate) if better(original.sense, candidate.objective, best.borrow().objective) => {
            if let Some(observe) = observe {
                observe(when, candidate.objective, &candidate.values);
            }
            best.replace(candidate);
        }
        Ok(_) => {}
        Err(e) => error.set(Some(e)),
    };
    let prior = PresolveOutcome {
        start,
        deadline,
        stats: BbStats {
            lp_iterations: form.lp_iterations,
            ..BbStats::default()
        },
        best: form.project_seed().and_then(|point| {
            heuristics::round_and_check(
                &form.problem,
                &Bounds::from_problem(&form.problem),
                &point,
                &opts.tol,
            )
        }),
        bound: None,
        status: None,
    };
    // 分类冲突的两个分支常有很不相同的代价，自动模式利用已有分支的目标变化选列。
    let mut work_opts = *opts;
    if work_opts.bb.branch_rule_auto {
        work_opts.bb.branch_rule_auto = false;
        work_opts.bb.branch_rule = branch_bound::BranchRule::Pseudocost;
    }
    let mut result = solve_prepared(
        &form.problem,
        &work_opts,
        ctrl,
        hooks,
        Some(&lift),
        start,
        Some(prior),
    )?;
    if let Some(e) = error.take() {
        return Err(e);
    }
    let mut best = best.into_inner();
    match (result.report.values.take(), result.report.objective) {
        (Some(point), Some(objective)) => {
            let candidate = restore(&point, objective)?;
            if result.report.status == SolveStatus::Optimal
                && best.objective + 1e-8 * (1.0 + objective.abs()) < candidate.objective
            {
                return Err(MilpError::Numerical("分类最优目标与已知可行解冲突".into()));
            }
            if better(original.sense, candidate.objective, best.objective) {
                best = candidate;
            }
        }
        (None, None) if result.report.status == SolveStatus::LimitReached => {}
        _ => return Err(MilpError::Numerical("分类模型终态与已知可行解冲突".into())),
    }
    if result.report.status == SolveStatus::Optimal {
        result.bound = Some(best.objective);
    } else if result.bound.is_some_and(|b| b > best.objective) {
        result.bound = None;
    }
    result.report.objective = Some(best.objective);
    result.report.values = Some(best.values);
    if deadline.is_some_and(|d| Instant::now() >= d) {
        result.report.status = SolveStatus::LimitReached;
    }
    result.report.stats.time_ms = start.elapsed().as_millis().min(u64::MAX as u128) as u64;
    Ok(result)
}

// 输入已验证，结构改写只在外层执行。已有的解、工作量和起点沿用到分支定界。
fn solve_prepared(
    problem: &Problem,
    opts: &SolveOptions,
    ctrl: &SearchControl,
    hooks: Option<&SearchHooks>,
    incumbent_observer: IncumbentObserver<'_>,
    start: Instant,
    prior: Option<PresolveOutcome>,
) -> Result<ControlledSolveReport> {
    let num_orig_vars = problem.vars.len();
    let mut route =
        prior.or_else(|| structured_prelude(problem, opts, ctrl, start, incumbent_observer));
    if let Some(r) = route.take_if(|r| r.status.is_some()) {
        return Ok(r.report(problem, hooks, 0));
    }

    // （可选）presolve → 缩减问题 + 还原映射。R4 起最前端先做问题级数值
    // 缩放（scale→presolve→solve→postsolve→unscale）：缩放使容差判定
    // 作用于均衡后量级（标准做法）；目标值不变（ĉ·ŷ = c·x）、正缩放保序
    // ⇒ incumbent 比较合法、割在缩放空间生成对原问题有效（线性双射保
    // 可行集）；终检与目标重算在**原问题**上做（unscale 后的 x），输出
    // 契约不变。2 的幂因子保证 unscale 除法逐位精确。
    let (pre_work, smap) = if opts.use_scaling {
        let (scaled, map) = crate::scaling::scale_problem(problem);
        (scaled, Some(map))
    } else {
        (problem.clone(), None)
    };
    let (work, map, reductions, preserved) = if opts.use_presolve {
        // 按最外层请求门控；并行小树前缀不能因临时 threads=1 捕获旧行。
        let capture = cfg!(feature = "certified_row_fallback") && ctrl.threads == 1
            && opts.bb.binary_row_propagation.resolve(ctrl.threads)
                == branch_bound::BinaryPropagationMode::On;
        let (reduced, map, preserved) = presolve::presolve_for_solve(&pre_work, &opts.tol, capture)?;
        // W5：聚合掉的重复列（组成员数 − 1）一并计入消减量。
        let reductions = map.fixed.len()
            + map.dropped_rows.len()
            + map
                .aggregations
                .iter()
                .map(|(members, _, _)| members.len().saturating_sub(1))
                .sum::<usize>();
        (reduced, Some(map), reductions, preserved)
    } else {
        (pre_work, None, 0, crate::binary_propagation::PreservedBinaryRows::default())
    };

    let restore_values = |reduced_x: Vec<f64>| -> Vec<f64> {
        let x_pre = match &map {
            Some(m) => presolve::postsolve(m, &reduced_x, num_orig_vars, &opts.tol),
            None => reduced_x,
        };
        match &smap {
            Some(s) => s.unscale_x(&x_pre),
            None => x_pre,
        }
    };
    let mut prelude = route
        .as_ref()
        .map_or_else(SearchPrelude::default, |r| SearchPrelude {
            start: Some(r.start),
            nodes: r.stats.nodes,
            lp_iterations: r.stats.lp_iterations,
            incumbent: r.best.as_ref().and_then(|incumbent| {
                work_seed(
                    &work, map.as_ref(), smap.as_ref(), incumbent, opts.tol, &restore_values,
                )
            }),
        });
    // 结构识别、编制选项及常规预处理已经使用的墙钟时间也在单线程总预算内。
    if ctrl.threads == 1 && opts.use_presolve {
        prelude.start = Some(start);
    }
    if let Some(mut r) = route.take_if(|r| r.time_exceeded()) {
        r.status = Some(SolveStatus::LimitReached);
        return Ok(r.report(problem, hooks, reductions));
    }
    let observe_work = |reduced_x: &[f64]| {
        let Some(observe) = incumbent_observer else { return };
        let accepted_at = Instant::now();
        let x = restore_values(reduced_x.to_vec());
        let objective = problem.obj.eval(&x) + problem.obj_offset;
        let known = route.as_ref().and_then(|r| r.best.as_ref());
        if known.is_none_or(|inc| !better(problem.sense, inc.objective, objective)) {
            observe(accepted_at, objective, &x);
        }
    };
    let work_observer = incumbent_observer.map(|_| &observe_work as &dyn Fn(&[f64]));
    let mut bb = branch_bound::branch_and_bound_with_preserved(
        &work, &opts.bb, ctrl, opts.tol, hooks, preserved, work_observer, prelude,
    )?;

    // postsolve → 缩放前变量空间 →（可选）unscale → 原问题空间；B&B 未
    // 带回取值则不带解返回。
    let (mut objective, mut values) = match bb.best_values.take() {
        Some(reduced_x) => {
            let mut x = restore_values(reduced_x);
            if x.len() != num_orig_vars {
                return Err(MilpError::Numerical(format!(
                    "解向量维数 {} 与原始变量数 {num_orig_vars} 不符",
                    x.len()
                )));
            }
            // 安全网终检（整数取整替换 + 逐约束复核），失败即 Err。
            final_check(problem, &mut x)?;
            // 目标以最终 values 重算（原 sense、含 obj_offset，值与目标自洽）。
            (Some(problem.obj.eval(&x) + problem.obj_offset), Some(x))
        }
        None => (None, None),
    };

    if let Some(r) = route {
        if let Some(inc) = r.best {
            if objective.is_none_or(|old| better(problem.sense, inc.objective, old)) {
                objective = Some(inc.objective);
                values = Some(inc.values);
                if bb.status != SolveStatus::LimitReached {
                    // 可信的原模型候选与通用终态冲突，保留候选和独立路线界。
                    bb.status = SolveStatus::LimitReached;
                    bb.bound = r.bound;
                } else if bb.bound.is_some_and(|bound| better(problem.sense, inc.objective, bound)) {
                    bb.bound = r.bound;
                    bb.gap_terminated = false;
                    bb.auto_degraded = false;
                }
                bb.best_objective = Some(inc.objective);
            }
        }
        if bb.status == SolveStatus::LimitReached && !bb.gap_terminated {
            if let Some(bound) = r.bound {
                if bb.bound.is_none_or(|old| better(problem.sense, old, bound)) {
                    bb.bound = Some(bound);
                }
            }
        }
    }

    // gap 证据：以 (最终 incumbent, 证据 bound) 经统一 clamp 重算——
    // BbResult.bound 在 gap 终止时即触发时刻证据 db（顺序路径触发即断，最终
    // incumbent == 触发时 incumbent，重算值 == 触发值；并行路径若停机收尾期
    // 间 incumbent 又有改进，重算值只会更小、仍非负有限）。故
    // gap_terminated = true 时必得非负有限 Some（触发判定即同一 db 上
    // gap_percent = Some(g ≤ P)）。
    let gap_percent = if bb.gap_terminated {
        gap_percent(
            problem.sense,
            bb.best_objective.unwrap_or(f64::NAN),
            bb.bound.unwrap_or(f64::NAN),
        )
    } else {
        None
    };

    Ok(ControlledSolveReport {
        report: SolveReport {
            status: bb.status,
            objective,
            values,
            stats: bb.stats,
            presolve_reductions: reductions,
        },
        gap_terminated: bb.gap_terminated,
        auto_degraded: bb.auto_degraded,
        gap_percent,
        bound: bb.bound,
    })
}

fn better(sense: Sense, new: f64, old: f64) -> bool {
    match sense {
        Sense::Minimize => new < old,
        Sense::Maximize => new > old,
    }
}

struct PresolveOutcome {
    start: Instant,
    deadline: Option<Instant>,
    stats: BbStats,
    best: Option<Incumbent>,
    bound: Option<f64>,
    status: Option<SolveStatus>,
}

impl PresolveOutcome {
    fn time_exceeded(&self) -> bool {
        self.deadline
            .is_some_and(|deadline| Instant::now() >= deadline)
    }

    fn report(
        mut self,
        problem: &Problem,
        hooks: Option<&SearchHooks>,
        reductions: usize,
    ) -> ControlledSolveReport {
        let (objective, values) = self
            .best
            .map_or((None, None), |inc| (Some(inc.objective), Some(inc.values)));
        self.stats.time_ms = self.start.elapsed().as_millis().min(u64::MAX as u128) as u64;
        if let Some(hooks) = hooks {
            branch_bound::emit_progress(
                hooks,
                &self.start,
                self.stats.nodes,
                objective,
                self.bound,
                problem.sense,
            );
        }
        ControlledSolveReport {
            report: SolveReport {
                status: self.status.unwrap_or(SolveStatus::LimitReached),
                objective,
                values,
                stats: self.stats,
                presolve_reductions: reductions,
            },
            gap_terminated: false,
            auto_degraded: false,
            gap_percent: None,
            bound: self.bound,
        }
    }
}

fn solve_integer_structure(
    problem: &Problem,
    opts: &SolveOptions,
    ctrl: &SearchControl,
    start: Instant,
    deadline: Option<Instant>,
    hooks: Option<&SearchHooks>,
    observe: IncumbentObserver<'_>,
) -> Result<Option<ControlledSolveReport>> {
    let mut best: Option<Incumbent> = None;
    let mut accept = |values: &[f64]| {
        // 分解器已按原行核验赋值；所有对外目标都在原列顺序下重算。
        let objective = problem.obj.eval(values) + problem.obj_offset;
        if best
            .as_ref()
            .is_none_or(|inc| better(problem.sense, objective, inc.objective))
        {
            if let Some(observe) = observe {
                observe(Instant::now(), objective, values);
            }
            best = Some(Incumbent {
                objective,
                values: values.to_vec(),
            });
        }
    };
    let mut root_result = crate::shared_demand::solve(problem, opts.tol, deadline, &mut accept)?;
    if root_result.is_none()
        && opts.bb.mitm_windows > 0
        && opts.bb.branch_rule_auto
        && opts.bb.node_rule == NodeRule::DepthFirst
    {
        root_result = heuristics::bucketed_equalities(problem, &opts.tol, deadline, &mut accept)?;
    }
    let outcome = if let Some((status, bound)) = root_result {
        Some((
            status,
            bound,
            BbStats {
                nodes: 1,
                ..BbStats::default()
            },
        ))
    } else if opts.bb.branch_rule_auto
        && opts.bb.node_rule == NodeRule::DepthFirst
        && opts.bb.lp_engine == crate::simplex_sparse::LpEngineKind::Auto
    {
        let limit = structural_node_limit(opts, ctrl);
        let mut result = None;
        for solve in [
            crate::boolean_network::solve,
            crate::packing::solve,
            crate::grouped_dp::solve,
            crate::weighted_cover::solve,
        ] {
            result = solve(problem, opts.tol, deadline, limit, &mut accept)?;
            if result.is_some() {
                break;
            }
        }
        result
    } else {
        None
    };
    let Some((mut status, bound, stats)) = outcome else {
        return Ok(None);
    };
    if status == SolveStatus::Optimal
        && (best.is_none() || best.as_ref().map(|inc| inc.objective) != bound)
    {
        return Err(MilpError::Numerical(
            "原模型目标与整数搜索证明不一致".into(),
        ));
    }
    if deadline.is_some_and(|end| Instant::now() >= end) {
        status = SolveStatus::LimitReached;
    }
    let result = PresolveOutcome {
        start,
        deadline,
        best,
        bound,
        status: Some(status),
        stats,
    };
    Ok(Some(result.report(problem, hooks, 0)))
}

// 节点上限为零表示不限。纯预算模式同时遵守墙钟参数折算的节点额度。
fn structural_node_limit(opts: &SolveOptions, ctrl: &SearchControl) -> u64 {
    let limit = if opts.bb.node_limit == 0 {
        u64::MAX
    } else {
        opts.bb.node_limit
    };
    if ctrl.determinism == Determinism::Budget && opts.bb.time_limit_ms > 0 {
        limit.min(branch_bound::budget_node_cap_of(&opts.bb))
    } else {
        limit
    }
}

fn structured_prelude(
    problem: &Problem,
    opts: &SolveOptions,
    ctrl: &SearchControl,
    start: Instant,
    observe: IncumbentObserver<'_>,
) -> Option<PresolveOutcome> {
    if !opts.use_presolve
        || !opts.bb.heuristics
        || ctrl.threads != 1
        || (opts.bb.node_rule == NodeRule::HybridBestBound && opts.bb.local_cuts)
    {
        return None;
    }
    let budget_mode = ctrl.determinism == Determinism::Budget;
    let deadline = if !budget_mode && opts.bb.time_limit_ms > 0 {
        start.checked_add(Duration::from_millis(opts.bb.time_limit_ms))
    } else {
        None
    };
    let node_limit = structural_node_limit(opts, ctrl);
    let mut r = PresolveOutcome {
        start,
        deadline,
        stats: BbStats::default(),
        best: None,
        bound: None,
        status: None,
    };
    let Some(model) = RouteModel::detect(problem) else {
        let mut report_cover = |values: &[f64]| {
            if let Some(observe) = observe {
                observe(
                    Instant::now(),
                    problem.obj.eval(values) + problem.obj_offset,
                    values,
                );
            }
        };
        if let Some((values, bound)) =
            crate::cover_repair::search(problem, deadline, &mut report_cover)
        {
            r.bound = bound;
            r.best = Some(Incumbent {
                objective: problem.obj.eval(&values) + problem.obj_offset,
                values,
            });
            if r.time_exceeded() {
                r.status = Some(SolveStatus::LimitReached);
            }
            return Some(r);
        }
        if let Some(bound) = crate::binary_tables::seed(
            problem, opts.tol, deadline, node_limit, &mut r.stats.nodes, |values| {
                let objective = problem.obj.eval(&values) + problem.obj_offset;
                if let Some(observe) = observe {
                    observe(Instant::now(), objective, &values);
                }
                r.best = Some(Incumbent { objective, values });
            },
        ) {
            r.bound = Some(bound);
            if r.time_exceeded() || r.stats.nodes >= node_limit {
                r.status = Some(SolveStatus::LimitReached);
            } else if r.best.as_ref().is_some_and(|inc| inc.objective == bound) {
                r.status = Some(SolveStatus::Optimal);
            }
            return Some(r);
        }
        let values = crate::periodic_equalities::seed(problem, deadline)?;
        let objective = problem.obj.eval(&values) + problem.obj_offset;
        if let Some(observe) = observe {
            observe(Instant::now(), objective, &values);
        }
        let improved = crate::periodic_equalities::polish(
            problem,
            &values,
            opts,
            deadline,
            &mut r.stats.lp_iterations,
        );
        let mut best = Incumbent { objective, values };
        if let Some(values) = improved {
            let new_cost = problem.obj.eval(&values) + problem.obj_offset;
            if better(problem.sense, new_cost, objective) {
                if let Some(observe) = observe {
                    observe(Instant::now(), new_cost, &values);
                }
                best = Incumbent {
                    objective: new_cost,
                    values,
                };
            }
        }
        let current = best.values.clone();
        let mut accept = |values: &[f64]| {
            let objective = problem.obj.eval(values) + problem.obj_offset;
            if better(problem.sense, objective, best.objective) {
                if let Some(observe) = observe {
                    observe(Instant::now(), objective, values);
                }
                best = Incumbent {
                    objective,
                    values: values.to_vec(),
                };
            }
        };
        r.bound = crate::periodic_equalities::improve(
            problem,
            &current,
            opts,
            deadline,
            node_limit,
            &mut r.stats,
            &mut accept,
        );
        r.best = Some(best);
        if r.time_exceeded() || r.stats.nodes >= node_limit {
            r.status = Some(SolveStatus::LimitReached);
        }
        return Some(r);
    };
    if !r.time_exceeded() {
        if let Some(values) = model.greedy_candidate(problem, deadline).or_else(|| {
            (!r.time_exceeded())
                .then(|| model.idle_candidate(problem))
                .flatten()
        }) {
            let objective = problem.obj.eval(&values) + problem.obj_offset;
            if let Some(observe) = observe {
                observe(Instant::now(), objective, &values);
            }
            r.best = Some(Incumbent { objective, values });
        }
    }
    // 有显式 gap 时保留通用搜索的证据与停机语义；贪心结果可供其使用。
    if !r.time_exceeded() && !ctrl.mipgap_percent.is_some_and(|gap| gap > 0.0) {
        let result = crate::route_cover::optimize(
            problem,
            &model,
            opts.tol,
            deadline,
            node_limit,
            r.best.as_ref().map(|inc| inc.values.as_slice()),
        );
        r.stats.nodes = result.nodes;
        r.stats.lp_iterations = result.lp_iterations;
        r.bound = result.bound.and_then(|bound| match problem.sense {
            Sense::Minimize => sum_bounds(problem.obj_offset, bound).map(|(lower, _)| lower),
            Sense::Maximize => sum_bounds(problem.obj_offset, -bound).map(|(_, upper)| upper),
        });
        if let Some(values) = result.best {
            let objective = problem.obj.eval(&values) + problem.obj_offset;
            // 证书对应的是精确比较选出的向量。即使浮点重算目标相同，
            // 也必须使用证书对应的向量。
            if result.proven_optimal
                || r.best
                    .as_ref()
                    .is_none_or(|inc| better(problem.sense, objective, inc.objective))
            {
                if let Some(observe) = observe {
                    if r.best.as_ref().is_none_or(|inc| inc.values != values) {
                        observe(Instant::now(), objective, &values);
                    }
                }
                r.best = Some(Incumbent { objective, values });
            }
            if result.proven_optimal {
                r.status = Some(SolveStatus::Optimal);
                r.bound = Some(objective);
            }
        }
        if r.status.is_none() && (result.interrupted || r.stats.nodes >= node_limit) {
            r.status = Some(SolveStatus::LimitReached);
        }
    }
    if r.time_exceeded() {
        r.status = Some(SolveStatus::LimitReached);
    }
    Some(r)
}

/// 缩放和预处理后先往返核验。原空间的候选由调用方单独保留，以防
/// 工作映射中的舍入或交换重复列改变了已核验的向量。
fn work_seed(
    work: &Problem,
    map: Option<&presolve::PostsolveMap>,
    scale: Option<&crate::scaling::ScaleMap>,
    incumbent: &Incumbent,
    tol: Tolerances,
    restore: &impl Fn(Vec<f64>) -> Vec<f64>,
) -> Option<Incumbent> {
    let scaled = match scale {
        Some(scale) => incumbent
            .values
            .iter()
            .zip(&scale.col)
            .map(|(x, factor)| x * factor)
            .collect(),
        None => incumbent.values.clone(),
    };
    let reduced = if let Some(map) = map {
        let mut reduced: Vec<_> = map
            .var_map
            .iter()
            .map(|&j| scaled.get(j).copied().unwrap_or(0.0))
            .collect();
        for (members, j, fixed) in &map.aggregations {
            if fixed.is_none() {
                *reduced.get_mut(*j)? = members.iter().map(|&k| scaled[k]).sum();
            }
        }
        reduced
    } else {
        scaled
    };
    let candidate = heuristics::round_and_check(work, &Bounds::from_problem(work), &reduced, &tol)?;
    (restore(candidate.values.clone()) == incumbent.values).then_some(candidate)
}

/// 统一求解入口（签名不变；默认受控路径的薄委托——
/// SearchControl::default + hooks=None，与原实现逐字节一致）。
pub fn solve(problem: &Problem, opts: &SolveOptions) -> Result<SolveReport> {
    Ok(solve_controlled(problem, opts, &SearchControl::default(), None)?.report)
}

/// 便捷入口：默认选项求解。
pub fn solve_default(problem: &Problem) -> Result<SolveReport> {
    solve(problem, &SolveOptions::default())
}

/// 安全网终检（就地修改 `x`）：
/// - 整数取整替换：Integer/Binary 变量 |x−round| ≤ [`FINAL_INTEGER_TOL`] 时
///   替换为 round 值，超差即 `Err(Numerical)`；
/// - 逐约束复核（对替换后的最终值）：violation > [`FINAL_VIOLATION_REL`]×
///   (1+|rhs|+‖row‖₁) 即 `Err(Numerical)`（‖row‖₁ = 行系数绝对值和）。
pub(crate) fn final_check(problem: &Problem, x: &mut [f64]) -> Result<()> {
    for (i, v) in problem.vars.iter().enumerate() {
        if matches!(v.var_type, VarType::Binary | VarType::Integer) {
            let val = x[i];
            let rounded = val.round();
            if (val - rounded).abs() > FINAL_INTEGER_TOL {
                return Err(MilpError::Numerical(format!(
                    "变量 `{}` 终值 {val} 距最近整数超过 {FINAL_INTEGER_TOL}",
                    v.name
                )));
            }
            x[i] = rounded;
        }
    }
    for (row, c) in problem.constraints.iter().enumerate() {
        let lhs = c.expr.eval(x);
        let violation = match c.op {
            RelOp::Le => (lhs - c.rhs).max(0.0),
            RelOp::Ge => (c.rhs - lhs).max(0.0),
            RelOp::Eq => (lhs - c.rhs).abs(),
        };
        let row_norm1: f64 = c.expr.terms.iter().map(|&(_, a)| a.abs()).sum();
        let allowed = FINAL_VIOLATION_REL * (1.0 + c.rhs.abs() + row_norm1);
        if !violation.is_finite() || violation > allowed {
            return Err(MilpError::Numerical(format!(
                "约束 {row}（{}）终检违反量 {violation} 超过容许值 {allowed}",
                c.name.as_deref().unwrap_or("<未命名>")
            )));
        }
    }
    Ok(())
}

// ===========================================================================
// 单元测试
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::branch_bound::{BranchRule, NodeRule};
    use crate::model::Sense;

    /// 关闭 presolve 的选项：隔离编排与终检本体（presolve 属跨模块路径）。
    fn opts_no_presolve() -> SolveOptions {
        SolveOptions {
            use_presolve: false,
            ..SolveOptions::default()
        }
    }

    /// 关闭 presolve 与根割/启发式：最小依赖路径（整数复核仍经 WP3）。
    fn opts_plain_bb() -> SolveOptions {
        SolveOptions {
            use_presolve: false,
            bb: crate::branch_bound::BbOptions {
                root_cuts: false,
                heuristics: false,
                ..crate::branch_bound::BbOptions::default()
            },
            ..SolveOptions::default()
        }
    }

    fn assert_close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 1e-7 * (1.0 + b.abs()),
            "数值不符: {a} vs {b}"
        );
    }

    #[test]
    fn incumbent_observer_recovers_original_columns_without_changing_result() {
        let mut p = Problem::new("观测原尺度", Sense::Maximize);
        p.add_var("x", VarType::Binary).unwrap();
        p.add_var("y", VarType::Continuous).unwrap();
        p.vars[1].lb = 4.0;
        p.vars[1].ub = 4.0;
        p.set_objective(vec![(0, 1.0), (1, 128.0)], 2.25).unwrap();
        let opts = SolveOptions::default();
        let ctrl = SearchControl::default();
        let normal = solve_controlled(&p, &opts, &ctrl, None).unwrap();
        let start = Instant::now();
        let events = std::cell::RefCell::new(Vec::new());
        let observed = solve_controlled_with_incumbent(&p, &opts, &ctrl, None, &|time, objective, x| {
            events.borrow_mut().push((time, objective, x.to_vec()));
        }).unwrap();
        let events = events.into_inner();
        assert!(!events.is_empty());
        assert!(events.iter().all(|(time, _, _)| *time >= start));
        assert_eq!(observed.report.status, normal.report.status);
        assert_eq!(observed.report.values, normal.report.values);
        let (_, goal, values) = events.last().unwrap();
        assert_eq!(Some(*goal), normal.report.objective);
        assert_eq!(values, &[1.0, 4.0]);
    }

    #[test]
    fn binary_public_bounds_cannot_produce_out_of_domain_solution() {
        let mut p = Problem::new("二元公开边界入口", Sense::Maximize);
        p.add_var("x", VarType::Binary).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        for use_presolve in [true, false] {
            let opts = SolveOptions { use_presolve, ..SolveOptions::default() };
            p.vars[0].lb = 2.0;
            p.vars[0].ub = 2.0;
            let result = solve(&p, &opts);
            assert!(result.is_err(), "非法二元界应在求解前拒绝，实际{result:?}");
            p.vars[0].lb = 0.2;
            p.vars[0].ub = 0.8;
            assert_eq!(solve(&p, &opts).unwrap().status, SolveStatus::Infeasible);
            p.vars[0].lb = 1.0;
            p.vars[0].ub = 1.0;
            let report = solve(&p, &opts).unwrap();
            assert_eq!(report.status, SolveStatus::Optimal);
            assert_eq!(report.objective, Some(1.0));
        }
    }

    /// SolveOptions 默认值逐字段符合契约。
    #[test]
    fn solve_options_default_contract() {
        let o = SolveOptions::default();
        let t = Tolerances::default();
        assert_eq!(o.tol.feasibility, t.feasibility);
        assert_eq!(o.tol.optimality, t.optimality);
        assert_eq!(o.tol.pivot, t.pivot);
        assert_eq!(o.tol.integrality, t.integrality);
        assert_eq!(o.tol.zero, t.zero);
        assert_eq!(o.bb.branch_rule, BranchRule::MostFractional);
        assert_eq!(o.bb.node_rule, NodeRule::DepthFirst);
        assert_eq!(o.bb.node_limit, 2_000_000);
        assert_eq!(o.bb.time_limit_ms, 0);
        assert!(o.bb.root_cuts);
        assert!(o.bb.heuristics);
        assert_eq!(o.bb.dive_max_fixes, 4096);
        assert_eq!(o.bb.pc_reliability_probes, 8);
        assert!(o.use_presolve);
    }

    /// 纯 LP（无整数变量、关 presolve）：Optimal、目标含 obj_offset、values
    /// 与原始问题同维同序、目标按最终值重算（自洽）。
    #[test]
    fn solve_lp_no_presolve() {
        // max 3x+2y s.t. x+y<=4, x+3y<=6, [0,3]^2 → (3,1)，目标 11 + 偏移 100。
        let mut p = Problem::new("lp", Sense::Maximize);
        p.add_var("x", VarType::Continuous).unwrap();
        p.add_var("y", VarType::Continuous).unwrap();
        p.set_var_bounds(0, 0.0, 3.0).unwrap();
        p.set_var_bounds(1, 0.0, 3.0).unwrap();
        p.set_objective(vec![(0, 3.0), (1, 2.0)], 100.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 4.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 3.0)], RelOp::Le, 6.0)
            .unwrap();
        let rep = solve(&p, &opts_no_presolve()).unwrap();
        assert_eq!(rep.status, SolveStatus::Optimal);
        assert_close(rep.objective.unwrap(), 111.0);
        let v = rep.values.as_ref().unwrap();
        assert_eq!(v.len(), 2, "values 与原始问题同维");
        assert_close(v[0], 3.0);
        assert_close(v[1], 1.0);
        // 目标自洽：以报告值重算。
        assert_close(p.obj.eval(v) + p.obj_offset, rep.objective.unwrap());
        assert_eq!(rep.presolve_reductions, 0);
        assert_eq!(rep.stats.nodes, 1);
    }

    /// 纯 LP 状态透传：Infeasible/Unbounded 时 objective 与 values 均为 None。
    #[test]
    fn solve_lp_status_passthrough() {
        let mut inf = Problem::new("inf", Sense::Minimize);
        inf.add_var("x", VarType::Continuous).unwrap();
        inf.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        inf.add_constraint(None, vec![(0, 1.0)], RelOp::Ge, 2.0)
            .unwrap();
        inf.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 1.0)
            .unwrap();
        let rep = solve(&inf, &opts_no_presolve()).unwrap();
        assert_eq!(rep.status, SolveStatus::Infeasible);
        assert!(rep.objective.is_none());
        assert!(rep.values.is_none());

        let mut unb = Problem::new("unb", Sense::Maximize);
        unb.add_var("x", VarType::Continuous).unwrap();
        unb.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        let rep = solve(&unb, &opts_no_presolve()).unwrap();
        assert_eq!(rep.status, SolveStatus::Unbounded);
        assert!(rep.objective.is_none());
        assert!(rep.values.is_none());
    }

    /// 非法模型（lb>ub）→ Err(InvalidModel)，不进求解。
    #[test]
    fn solve_invalid_model_returns_err() {
        let mut p = Problem::new("bad", Sense::Minimize);
        p.add_var("x", VarType::Continuous).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.vars[0].lb = 5.0;
        p.vars[0].ub = 1.0;
        assert!(matches!(
            solve(&p, &opts_no_presolve()),
            Err(MilpError::InvalidModel(_))
        ));
    }

    /// 终检的整数替换与约束复核直接可用纯数据验证：
    /// 偏差在 1e-5 内的整数值被替换为 round，violation 超相对容差即 Err。
    #[test]
    fn final_check_rounds_and_rejects() {
        let mut p = Problem::new("chk", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 10.0)
            .unwrap();
        // 偏差 1e-6 ≤ 1e-5：替换为 round。
        let mut x = vec![3.0 + 1e-6];
        assert!(final_check(&p, &mut x).is_ok());
        assert_eq!(x[0], 3.0);
        // 偏差超过 1e-5：Err(Numerical)。
        let mut x = vec![3.0 + 1e-4];
        assert!(matches!(
            final_check(&p, &mut x),
            Err(MilpError::Numerical(_))
        ));
        // 约束违反超过 1e-6×(1+|rhs|+‖row‖₁)：Err(Numerical)。
        let mut x = vec![12.0];
        assert!(matches!(
            final_check(&p, &mut x),
            Err(MilpError::Numerical(_))
        ));
    }

    /// （跨模块联调：整数节点复核经 WP3 heuristics::round_and_check）
    /// MILP（关 presolve、关割/启发式）：max x s.t. 2x≤7 → 3.0，整数已取整、
    /// 目标自洽。
    #[test]
    fn solve_milp_no_presolve_integration() {
        let mut p = Problem::new("m", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 2.0)], RelOp::Le, 7.0)
            .unwrap();
        let rep = solve(&p, &opts_plain_bb()).unwrap();
        assert_eq!(rep.status, SolveStatus::Optimal);
        assert_close(rep.objective.unwrap(), 3.0);
        let v = rep.values.as_ref().unwrap();
        assert_eq!(v.len(), 1);
        assert_close(v[0], 3.0);
        assert_close(p.obj.eval(v) + p.obj_offset, rep.objective.unwrap());
        assert_eq!(rep.presolve_reductions, 0);
    }

    /// （跨模块联调：presolve/cut_loop/heuristics 全链路）
    /// 默认选项（presolve + 根割 + 启发式）。
    #[test]
    fn solve_milp_default_integration() {
        let mut p = Problem::new("m", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 2.0)], RelOp::Le, 7.0)
            .unwrap();
        let rep = solve_default(&p).unwrap();
        assert_eq!(rep.status, SolveStatus::Optimal);
        assert_close(rep.objective.unwrap(), 3.0);
        assert_close(rep.values.as_ref().unwrap()[0], 3.0);
    }

    /// （跨模块联调：默认选项下的多变量背包）
    /// max 10a+7b+6c+4d+2e s.t. 5a+4b+3c+2d+e ≤ 9，0-1 变量 → 最优 18
    /// （a=c=e=1，重量恰 9）。
    #[test]
    fn solve_default_knapsack5_integration() {
        let mut p = Problem::new("k5", Sense::Maximize);
        for name in ["a", "b", "c", "d", "e"] {
            p.add_var(name, VarType::Binary).unwrap();
        }
        p.set_objective(vec![(0, 10.0), (1, 7.0), (2, 6.0), (3, 4.0), (4, 2.0)], 0.0)
            .unwrap();
        p.add_constraint(
            None,
            vec![(0, 5.0), (1, 4.0), (2, 3.0), (3, 2.0), (4, 1.0)],
            RelOp::Le,
            9.0,
        )
        .unwrap();
        let rep = solve_default(&p).unwrap();
        assert_eq!(rep.status, SolveStatus::Optimal);
        assert_close(rep.objective.unwrap(), 18.0);
        let v = rep.values.as_ref().unwrap();
        for (j, expect) in [(0usize, 1.0), (1, 0.0), (2, 1.0), (3, 0.0), (4, 1.0)] {
            assert_close(v[j], expect);
        }
    }

    /// （跨模块联调：经 WP3 presolve/postsolve/round_and_check）
    /// postsolve 还原维数：固定列 x≡1 被删除后，报告必须回到原始维数并回填
    /// 固定值（x=1，y=6，目标 6）。
    #[test]
    fn postsolve_restores_dimension_integration() {
        let mut p = Problem::new("fix", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.add_var("y", VarType::Integer).unwrap();
        p.set_var_bounds(0, 1.0, 1.0).unwrap();
        p.set_objective(vec![(1, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 7.7)
            .unwrap();
        let rep = solve(&p, &SolveOptions::default()).unwrap();
        assert_eq!(rep.status, SolveStatus::Optimal);
        let v = rep.values.as_ref().unwrap();
        assert_eq!(v.len(), 2, "postsolve 必须还原原始维数");
        assert!((v[0] - 1.0).abs() < 1e-9, "固定列回填固定值");
        assert!((v[1] - 6.0).abs() < 1e-9);
        assert_close(rep.objective.unwrap(), 6.0);
        assert!(rep.presolve_reductions >= 1, "固定列删除应计入化简量");
    }

    // =====================================================================
    // 受控入口（solve_controlled / ControlledSolveReport）
    // =====================================================================

    /// solve() = solve_controlled(默认) 的薄委托：报告全等。
    #[test]
    fn solve_controlled_default_matches_solve() {
        let mut p = Problem::new("m", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 2.0)], RelOp::Le, 7.0).unwrap();
        for opts in [opts_no_presolve(), SolveOptions::default()] {
            let a = solve(&p, &opts).unwrap();
            let b = solve_controlled(&p, &opts, &SearchControl::default(), None).unwrap();
            assert_eq!(a.status, b.report.status);
            assert_eq!(a.objective, b.report.objective);
            assert_eq!(a.values, b.report.values);
            assert_eq!(a.stats.nodes, b.report.stats.nodes);
            assert_eq!(a.presolve_reductions, b.report.presolve_reductions);
            assert!(!b.gap_terminated);
            assert_eq!(b.gap_percent, None);
        }
    }

    /// 受控顺序 gap 用例（dive 问题，关 presolve 与根割保证确定性；
    /// offset=0.5 有意非整数化——整数目标 cutoff 门控判否，保住 gap 路径
    /// 覆盖；offset=0.0 时根 bound 21.5 ≤ 21+1−ε 会被 cutoff 在 (d) 直接
    /// 闭合，终态翻转为完备 Optimal）：
    /// - mipgap=Some(3.0) ⇒ report.status=LimitReached（携带 incumbent 解）、
    ///   gap_terminated=true、gap_percent=Some(≈2.33)（非负有限）、bound≈22.0；
    /// - mipgap=Some(2.0) ⇒ Optimal 21.5、gap_percent=None。
    #[test]
    fn solve_controlled_sequential_gap_dive_problem() {
        let mut p = Problem::new("dive-gap", Sense::Maximize);
        p.add_var("x1", VarType::Integer).unwrap();
        p.add_var("x2", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_var_bounds(1, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 1.0)], 0.5).unwrap();
        p.add_constraint(None, vec![(0, 20.0), (1, 3.0)], RelOp::Le, 43.0)
            .unwrap();
        let opts = SolveOptions {
            use_presolve: false,
            bb: crate::branch_bound::BbOptions {
                root_cuts: false,
                heuristics: true,
                ..crate::branch_bound::BbOptions::default()
            },
            ..SolveOptions::default()
        };
        // Some(3.0)：gap=(22.0−21.5)/21.5×100≈2.326 ≤ 3 ⇒ LimitReached + 证据。
        let ctrl = SearchControl {
            mipgap_percent: Some(3.0),
            ..SearchControl::default()
        };
        let r = solve_controlled(&p, &opts, &ctrl, None).unwrap();
        assert_eq!(r.report.status, SolveStatus::LimitReached);
        assert!(r.gap_terminated);
        assert_close(r.report.objective.unwrap(), 21.5);
        assert_close(r.report.values.as_ref().unwrap()[0], 2.0);
        assert_close(r.bound.unwrap(), 22.0);
        let g = r.gap_percent.expect("gap_terminated ⇒ gap_percent 必为 Some");
        assert!(g.is_finite() && g >= 0.0, "gap_percent 非负有限，实际 {g}");
        assert!((g - 2.3255813953488373).abs() < 1e-9, "实际 {g}");
        // Some(2.0)：不触发 ⇒ Optimal。
        let ctrl = SearchControl {
            mipgap_percent: Some(2.0),
            ..SearchControl::default()
        };
        let r = solve_controlled(&p, &opts, &ctrl, None).unwrap();
        assert_eq!(r.report.status, SolveStatus::Optimal);
        assert!(!r.gap_terminated);
        assert_eq!(r.gap_percent, None);
        assert_close(r.report.objective.unwrap(), 21.5);
        assert_close(r.bound.unwrap(), 21.5);
    }

    /// 策略级自动降级（顺序路径；dive 问题关 presolve/根割保证确定性；
    /// offset=0.5 非整数化保住降级路径覆盖，见上用例注释）：
    /// - 激进参数（start_frac=0 / plateau_frac=0）+ auto + 时间预算 ⇒ 双停滞
    ///   判定即时触发：LimitReached + gap_terminated + auto_degraded，携带
    ///   incumbent 21.5 与 gap 证据（≈2.33%，同 mipgap 口径）；
    /// - 同参数但**无时间预算** ⇒ 契约红线：auto 仍证明最优（Optimal 21.5）；
    /// - 时间预算但**非 auto 模式** ⇒ 降级不评估：Optimal 21.5。
    #[test]
    fn solve_controlled_auto_degrade_sequential() {
        let mut p = Problem::new("dive-gap", Sense::Maximize);
        p.add_var("x1", VarType::Integer).unwrap();
        p.add_var("x2", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_var_bounds(1, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 1.0)], 0.5).unwrap();
        p.add_constraint(None, vec![(0, 20.0), (1, 3.0)], RelOp::Le, 43.0)
            .unwrap();
        let aggressive = crate::search_ctrl::DegradePolicy {
            start_frac: 0.0,
            plateau_frac: 0.0,
            ..crate::search_ctrl::DegradePolicy::default()
        };
        let base = SolveOptions {
            use_presolve: false,
            bb: crate::branch_bound::BbOptions {
                root_cuts: false,
                heuristics: true,
                ..crate::branch_bound::BbOptions::default()
            },
            ..SolveOptions::default()
        };
        // ① auto + 预算 + 激进参数 ⇒ 触发降级（诚实 LimitReached + 标记）。
        let mut o = base;
        o.bb.time_limit_ms = 10_000;
        let ctrl = SearchControl {
            adaptive_heuristics: true,
            degrade: aggressive,
            ..SearchControl::default()
        };
        let r = solve_controlled(&p, &o, &ctrl, None).unwrap();
        assert_eq!(r.report.status, SolveStatus::LimitReached);
        assert!(r.gap_terminated, "降级 ⇒ gap_terminated");
        assert!(r.auto_degraded, "降级 ⇒ auto_degraded");
        assert_close(r.report.objective.unwrap(), 21.5);
        assert_close(r.report.values.as_ref().unwrap()[0], 2.0);
        let g = r.gap_percent.expect("降级 ⇒ gap_percent 必为 Some");
        assert!(g.is_finite() && g >= 0.0, "gap 非负有限，实际 {g}");
        // ② 同参数但无预算 ⇒ 证明最优（契约红线：无需求表达不降级）。
        let ctrl = SearchControl {
            adaptive_heuristics: true,
            degrade: aggressive,
            ..SearchControl::default()
        };
        let r = solve_controlled(&p, &base, &ctrl, None).unwrap();
        assert_eq!(r.report.status, SolveStatus::Optimal);
        assert!(!r.gap_terminated);
        assert!(!r.auto_degraded);
        assert_close(r.report.objective.unwrap(), 21.5);
        // ③ 预算但非 auto 模式 ⇒ 降级不评估。
        let ctrl = SearchControl {
            adaptive_heuristics: false,
            degrade: aggressive,
            ..SearchControl::default()
        };
        let r = solve_controlled(&p, &o, &ctrl, None).unwrap();
        assert_eq!(r.report.status, SolveStatus::Optimal);
        assert!(!r.gap_terminated);
        assert!(!r.auto_degraded);
    }

    // =====================================================================
    // S3：revised 引擎全路径（--lp-engine revised 的库层口径）
    // =====================================================================

    /// revised 选项（threads=1 默认；关 presolve 隔离 B&B 本体）。
    fn opts_revised() -> SolveOptions {
        SolveOptions {
            use_presolve: false,
            bb: crate::branch_bound::BbOptions {
                lp_engine: crate::simplex_sparse::LpEngineKind::Revised,
                ..crate::branch_bound::BbOptions::default()
            },
            ..SolveOptions::default()
        }
    }

    fn bench_instance(name: &str) -> std::path::PathBuf {
        std::path::PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/../bench/instances"))
            .join(format!("{name}.milp"))
    }

    fn milp_instance(name: &str) -> std::path::PathBuf {
        std::path::PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/instances"))
            .join(format!("{name}.milp"))
    }

    /// revised 全路径小 MILP：背包 5 项 → 18（同 tableau 钉值）。
    #[test]
    fn revised_full_path_knapsack5() {
        let mut p = Problem::new("k5", Sense::Maximize);
        for name in ["a", "b", "c", "d", "e"] {
            p.add_var(name, VarType::Binary).unwrap();
        }
        p.set_objective(vec![(0, 10.0), (1, 7.0), (2, 6.0), (3, 4.0), (4, 2.0)], 0.0)
            .unwrap();
        p.add_constraint(
            None,
            vec![(0, 5.0), (1, 4.0), (2, 3.0), (3, 2.0), (4, 1.0)],
            RelOp::Le,
            9.0,
        )
        .unwrap();
        let rep = solve(&p, &opts_revised()).unwrap();
        assert_eq!(rep.status, SolveStatus::Optimal);
        assert_close(rep.objective.unwrap(), 18.0);
    }

    /// revised × milp/instances 全路径对账（presolve/根割/启发式/局部割全开，
    /// 默认参数）：目标值与 tableau 引擎恒同（状态恒同；节点/lpi 允许不同）。
    #[test]
    fn revised_matches_tableau_on_milp_instances() {
        for name in [
            "knapsack-10",
            "lp-prod-3",
            "milp-mixed-5",
            "milp-wolsey-2",
            "assign-4x4",
            "hard-knapsack-12",
        ] {
            let p = crate::instance::load_instance(&milp_instance(name)).unwrap();
            let tab = solve_default(&p).unwrap();
            let rev = solve(&p, &opts_revised_with_presolve()).unwrap();
            assert_eq!(
                format!("{:?}", tab.status),
                format!("{:?}", rev.status),
                "{name}: 终态分叉"
            );
            if tab.status == SolveStatus::Optimal {
                let a = tab.objective.unwrap();
                let b = rev.objective.unwrap();
                assert!(
                    (a - b).abs() <= 1e-6 * (1.0 + a.abs()),
                    "{name}: 目标分叉 tableau={a} revised={b}"
                );
            }
        }
    }

    /// revised 选项（presolve 开启——与默认 tableau 完全同参数，仅引擎不同）。
    fn opts_revised_with_presolve() -> SolveOptions {
        SolveOptions {
            bb: crate::branch_bound::BbOptions {
                lp_engine: crate::simplex_sparse::LpEngineKind::Revised,
                ..crate::branch_bound::BbOptions::default()
            },
            ..SolveOptions::default()
        }
    }

    /// bench 轻负载对账：setcover-100x150 已知最优 276（revised 全路径）。
    #[test]
    fn revised_setcover_100x150_known_optimum_276() {
        let p = crate::instance::load_instance(&bench_instance("setcover-100x150")).unwrap();
        let rep = solve(&p, &opts_revised_with_presolve()).unwrap();
        assert_eq!(rep.status, SolveStatus::Optimal);
        assert!(
            (rep.objective.unwrap() - 276.0).abs() <= 1e-6,
            "obj = {}",
            rep.objective.unwrap()
        );
    }

    /// F5 守卫：revised × threads≥2 在库层显式报错（不静默回 tableau）。
    #[test]
    fn revised_with_threads_rejected_loudly() {
        let mut p = Problem::new("m", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 2.0)], RelOp::Le, 7.0).unwrap();
        let ctrl = SearchControl {
            threads: 2,
            ..SearchControl::default()
        };
        let r = solve_controlled(&p, &opts_revised_with_presolve(), &ctrl, None);
        assert!(matches!(r, Err(MilpError::InvalidModel(_))), "实际 {r:?}");
        // threads=1 同参数正常。
        let r = solve_controlled(&p, &opts_revised_with_presolve(), &SearchControl::default(), None);
        assert!(r.is_ok());
    }

    /// t1 确定性（revised）：同引擎同参数双跑全部行为字段逐字节相同
    ///（目标/nodes/lpi/values 按位；mknap-3x100 带 B&B 深树路径）。
    #[test]
    fn revised_deterministic_double_run() {
        let p = crate::instance::load_instance(&bench_instance("mknap-3x100")).unwrap();
        let a = solve(&p, &opts_revised_with_presolve()).unwrap();
        let b = solve(&p, &opts_revised_with_presolve()).unwrap();
        assert_eq!(a.status, b.status);
        assert_eq!(
            a.objective.unwrap().to_bits(),
            b.objective.unwrap().to_bits()
        );
        assert_eq!(a.stats.nodes, b.stats.nodes);
        assert_eq!(a.stats.lp_iterations, b.stats.lp_iterations);
        assert_eq!(a.stats.cuts, b.stats.cuts);
        assert_eq!(a.values, b.values);
        // 目标为已知最优（对账锚：mknap-3x100 = 2331，tableau 口径与
        // python_bench micpython/highs/scip 三方记录一致）。
        assert!(
            (a.objective.unwrap() - 2331.0).abs() <= 1e-6,
            "obj = {}",
            a.objective.unwrap()
        );
    }

    /// W2 时限注入兼容（revised 全路径）：极短时限 ⇒ LimitReached 诚实终态
    ///（不冒充 Optimal）；无时限 ⇒ Optimal。
    #[test]
    fn revised_time_limit_truncates_honestly() {
        let p = crate::instance::load_instance(&bench_instance("setcover-100x150")).unwrap();
        let opts = SolveOptions {
            bb: crate::branch_bound::BbOptions {
                lp_engine: crate::simplex_sparse::LpEngineKind::Revised,
                time_limit_ms: 1,
                ..crate::branch_bound::BbOptions::default()
            },
            ..SolveOptions::default()
        };
        let rep = solve(&p, &opts).unwrap();
        assert_eq!(rep.status, SolveStatus::LimitReached);
        let full = solve(&p, &opts_revised_with_presolve()).unwrap();
        assert_eq!(full.status, SolveStatus::Optimal);
    }

    /// S4 类自适应规则解析：m×(n+m) ≥ 10^6 ⇒ revised，否则 tableau；
    /// 显式指定覆盖优先；SolveOptions/BbOptions 默认 Auto。
    #[test]
    fn auto_engine_resolution_rule() {
        use crate::simplex_sparse::{
            effective_engine, LpEngineKind, ENGINE_AUTO_CELL_THRESHOLD,
        };
        // 默认 Auto（三层同源）。
        assert_eq!(SolveOptions::default().bb.lp_engine, LpEngineKind::Auto);
        assert_eq!(
            crate::branch_bound::BbOptions::default().lp_engine,
            LpEngineKind::Auto
        );
        // 合成规模：阈下（m=300, n=600 ⇒ 270k）与阈上（m=800, n=4000 ⇒
        // 3.84M；lp-sparse 同款规模）。
        let mut small = Problem::new("s", Sense::Maximize);
        for j in 0..600 {
            small.add_var(&format!("x{j}"), VarType::Continuous).unwrap();
        }
        for _ in 0..300 {
            small
                .add_constraint(None, vec![(0, 1.0)], RelOp::Le, 1.0)
                .unwrap();
        }
        let mut big = Problem::new("b", Sense::Maximize);
        for j in 0..4000 {
            big.add_var(&format!("x{j}"), VarType::Continuous).unwrap();
        }
        for _ in 0..800 {
            big.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 1.0)
                .unwrap();
        }
        assert_eq!(
            effective_engine(&small, LpEngineKind::Auto),
            LpEngineKind::Tableau,
            "270k 格 < 10^6 ⇒ tableau"
        );
        assert_eq!(
            effective_engine(&big, LpEngineKind::Auto),
            LpEngineKind::Revised,
            "3.84M 格 ≥ 10^6 ⇒ revised（lp-sparse 同款规模）"
        );
        // 显式覆盖优先（大算例显式 tableau / 小算例显式 revised 均原样）。
        assert_eq!(
            effective_engine(&big, LpEngineKind::Tableau),
            LpEngineKind::Tableau
        );
        assert_eq!(
            effective_engine(&small, LpEngineKind::Revised),
            LpEngineKind::Revised
        );
        assert_eq!(ENGINE_AUTO_CELL_THRESHOLD, 1_000_000);
    }

    /// Auto 默认在小算例上走 tableau（与显式 tableau 同目标；默认组合
    /// 阈下行为不变的应用层证据）。
    #[test]
    fn auto_default_small_instance_matches_tableau() {
        let mut p = Problem::new("k5", Sense::Maximize);
        for name in ["a", "b", "c", "d", "e"] {
            p.add_var(name, VarType::Binary).unwrap();
        }
        p.set_objective(vec![(0, 10.0), (1, 7.0), (2, 6.0), (3, 4.0), (4, 2.0)], 0.0)
            .unwrap();
        p.add_constraint(
            None,
            vec![(0, 5.0), (1, 4.0), (2, 3.0), (3, 2.0), (4, 1.0)],
            RelOp::Le,
            9.0,
        )
        .unwrap();
        let auto = solve_default(&p).unwrap();
        let tab = solve(
            &p,
            &SolveOptions {
                bb: crate::branch_bound::BbOptions {
                    lp_engine: crate::simplex_sparse::LpEngineKind::Tableau,
                    ..crate::branch_bound::BbOptions::default()
                },
                ..SolveOptions::default()
            },
        )
        .unwrap();
        assert_eq!(auto.status, tab.status);
        assert_eq!(auto.objective, tab.objective);
        assert_eq!(auto.stats.nodes, tab.stats.nodes);
        assert_close(auto.objective.unwrap(), 18.0);
    }

    /// Auto × threads≥2 × 大算例（有效引擎 = revised）在库层显式报错
    ///（解析层不知规模不拦 auto；小算例 auto→tableau × threads≥2 正常）。
    #[test]
    fn auto_engine_threads_guard_by_effective_kind() {
        let mut big = Problem::new("b", Sense::Maximize);
        for j in 0..4000 {
            big.add_var(&format!("x{j}"), VarType::Continuous).unwrap();
        }
        for _ in 0..800 {
            big.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 1.0)
                .unwrap();
        }
        let ctrl2 = SearchControl {
            threads: 2,
            ..SearchControl::default()
        };
        // 关 presolve：合成大问题会被预处理清空（简约后 m=0 ⇒ 有效引擎
        // 回 tableau，守卫不触发）；关掉后问题原样到达派发层。
        let opts_big = SolveOptions {
            use_presolve: false,
            ..SolveOptions::default()
        };
        let r = solve_controlled(&big, &opts_big, &ctrl2, None);
        assert!(matches!(r, Err(MilpError::InvalidModel(_))), "实际 {r:?}");
        // 小算例（auto→tableau）× threads≥2 走并行路径正常求解。
        let mut small = Problem::new("m", Sense::Maximize);
        small.add_var("x", VarType::Integer).unwrap();
        small.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        small.add_constraint(None, vec![(0, 2.0)], RelOp::Le, 7.0).unwrap();
        let r = solve_controlled(&small, &SolveOptions::default(), &ctrl2, None);
        assert!(r.is_ok());
    }

    /// t1 确定性（Auto 默认在大算例上解析为 revised）：同参数双跑全部
    /// 行为字段逐字节相同（lp-wide-1000x2000，6M 格；重负载 #[ignore]，
    /// 与 acceptance 显式跑同款）。
    #[test]
    #[ignore]
    fn auto_engine_deterministic_double_run_large() {
        let p = crate::instance::load_instance(&bench_instance("lp-wide-1000x2000")).unwrap();
        assert_eq!(
            crate::simplex_sparse::effective_engine(
                &p,
                crate::simplex_sparse::LpEngineKind::Auto
            ),
            crate::simplex_sparse::LpEngineKind::Revised,
            "前提：该算例按规则解析为 revised"
        );
        let a = solve_default(&p).unwrap();
        let b = solve_default(&p).unwrap();
        assert_eq!(a.status, b.status);
        assert_eq!(a.objective.unwrap().to_bits(), b.objective.unwrap().to_bits());
        assert_eq!(a.stats.nodes, b.stats.nodes);
        assert_eq!(a.stats.lp_iterations, b.stats.lp_iterations);
        assert_eq!(a.values, b.values);
        assert!(
            (a.objective.unwrap() - 112.600268).abs() <= 1e-5,
            "obj = {}（S4 A/B 双引擎一致值 112.600268）",
            a.objective.unwrap()
        );
    }

    /// 验收目标（重负载，#[ignore]；`cargo test --release --lib solver::tests
    /// -- --ignored --test-threads=1` 显式跑）：revised 下 mknap-5x300 =
    /// 4762、setcover-300x600 = 144（已知最优，与 tableau 口径一致）。
    #[test]
    #[ignore]
    fn revised_acceptance_mknap_5x300_and_setcover_300x600() {
        let t0 = std::time::Instant::now();
        let p = crate::instance::load_instance(&bench_instance("mknap-5x300")).unwrap();
        let rep = solve(&p, &opts_revised_with_presolve()).unwrap();
        let mknap_secs = t0.elapsed().as_secs_f64();
        assert_eq!(rep.status, SolveStatus::Optimal, "mknap-5x300 终态");
        assert!(
            (rep.objective.unwrap() - 4762.0).abs() <= 1e-6,
            "mknap-5x300 obj = {}（钉 4762）",
            rep.objective.unwrap()
        );
        let t0 = std::time::Instant::now();
        let p = crate::instance::load_instance(&bench_instance("setcover-300x600")).unwrap();
        let rep = solve(&p, &opts_revised_with_presolve()).unwrap();
        let setcover_secs = t0.elapsed().as_secs_f64();
        assert_eq!(rep.status, SolveStatus::Optimal, "setcover-300x600 终态");
        assert!(
            (rep.objective.unwrap() - 144.0).abs() <= 1e-6,
            "setcover-300x600 obj = {}（钉 144）",
            rep.objective.unwrap()
        );
        eprintln!(
            "revised acceptance: mknap-5x300=4762 @{mknap_secs:.1}s, setcover-300x600=144 @{setcover_secs:.1}s"
        );
    }
}

#[cfg(all(test, feature = "certified_row_fallback"))]
mod certified_fallback_tests {
    use super::*;
    use crate::binary_propagation::{CaptureAudit, TestDeadlineGuard};
    use crate::branch_bound::{BinaryPropagationMode, BranchRule, NodeRule};
    use crate::presolve::certified_fallback_tests::branch_model;
    use crate::search_ctrl::Determinism;
    use crate::simplex_sparse::LpEngineKind;

    fn options(engine: LpEngineKind) -> SolveOptions {
        SolveOptions { use_scaling: false, bb: BbOptions { lp_engine: engine,
            binary_row_propagation: BinaryPropagationMode::Auto, branch_rule_auto: false,
            branch_rule: BranchRule::MostFractional, heuristics: false, root_cuts: false,
            local_cuts: false, root_rc_fixing: false, ..BbOptions::default() }, ..SolveOptions::default() }
    }

    #[test]
    fn certified_fallback_solver_matches_known_optimum_across_engines_and_domains() {
        let p = branch_model();
        let mut total_fixes = 0;
        for engine in [LpEngineKind::Tableau, LpEngineKind::Revised] {
            for scaling in [false, true] {
                for rule in [NodeRule::DepthFirst, NodeRule::HybridBestBound] {
                    for rc in [false, true] {
                        let mut opts = options(engine);
                        opts.use_scaling = scaling;
                        opts.bb.node_rule = rule;
                        opts.bb.root_rc_fixing = rc;
                        let auto = solve(&p, &opts).unwrap();
                        opts.bb.binary_row_propagation = BinaryPropagationMode::On;
                        let on = solve(&p, &opts).unwrap();
                        opts.bb.binary_row_propagation = BinaryPropagationMode::Off;
                        let off = solve(&p, &opts).unwrap();
                        for report in [&auto, &on, &off] {
                            assert_eq!(report.status, SolveStatus::Optimal);
                            assert_eq!(report.objective, Some(2.125));
                            let mut values = report.values.clone().unwrap();
                            final_check(&p, &mut values).unwrap();
                        }
                        assert_eq!(auto.stats.nodes, on.stats.nodes);
                        assert_eq!(auto.stats.lp_iterations, on.stats.lp_iterations);
                        assert_eq!(auto.values, on.values);
                        assert_eq!(auto.stats.propagation.fixes, on.stats.propagation.fixes);
                        assert_eq!(off.stats.propagation.accepted_rows, 0);
                        assert!(auto.stats.propagation.calls > 0);
                        assert_eq!(auto.stats.propagation.accepted_rows, 2);
                        assert_eq!(auto.stats.propagation.skipped_rows, 0);
                        total_fixes += auto.stats.propagation.fixes;
                    }
                }
            }
        }
        assert!(total_fixes > 0, "必须实际使用被救回的传播能力");
    }

    #[test]
    fn certified_fallback_capture_gate_and_root_lazy_statistics() {
        let p = branch_model();
        for presolve in [false, true] {
            for mode in [BinaryPropagationMode::Auto, BinaryPropagationMode::On, BinaryPropagationMode::Off] {
                let mut opts = options(LpEngineKind::Tableau);
                opts.use_presolve = presolve;
                opts.bb.binary_row_propagation = mode;
                opts.bb.node_limit = 1;
                let audit = CaptureAudit::new();
                let result = solve(&p, &opts).unwrap();
                assert_eq!(result.status, SolveStatus::LimitReached);
                let expected = if presolve && mode != BinaryPropagationMode::Off { 1 } else { 0 };
                assert_eq!(audit.count(), expected);
                assert_eq!(result.stats.propagation.calls, 0);
                assert_eq!(result.stats.propagation.accepted_rows, 0,
                    "未展开孩子时不能提前发布传播数据库统计");
            }
        }
        let mut root_closed = p;
        root_closed.set_objective(vec![(0, -1.0), (1, -2.0), (2, -3.0)], 0.25).unwrap();
        let audit = CaptureAudit::new();
        let closed = solve(&root_closed, &options(LpEngineKind::Tableau)).unwrap();
        assert_eq!(closed.status, SolveStatus::Optimal);
        assert_eq!(closed.objective, Some(0.25));
        assert_eq!(audit.count(), 1, "根部闭合也已支付旧行认证成本");
        assert_eq!(closed.stats.propagation.accepted_rows, 0);
    }

    #[test]
    fn certified_fallback_parallel_prefix_stays_off_and_explicit_on_errors() {
        let p = branch_model();
        let opts = options(LpEngineKind::Tableau);
        let ctrl = SearchControl { threads: 2, ..SearchControl::default() };
        let audit = CaptureAudit::new();
        let parallel = {
            // 若最外层Auto被顺序前缀重开，传播会立即中断并使本断言失败。
            let _deadline = TestDeadlineGuard::on_check(1);
            solve_controlled(&p, &opts, &ctrl, None).unwrap()
        };
        assert_eq!(parallel.report.status, SolveStatus::Optimal);
        assert_eq!(parallel.report.objective, Some(2.125));
        assert_eq!(parallel.report.stats.propagation.calls, 0);
        assert_eq!(audit.count(), 0);
        let mut explicit = opts;
        explicit.bb.binary_row_propagation = BinaryPropagationMode::On;
        let error = solve_controlled(&p, &explicit, &ctrl, None).unwrap_err();
        assert!(error.to_string().contains("二元静态行传播"));
        assert_eq!(audit.count(), 0);
    }

    #[test]
    fn certified_fallback_deadline_keeps_frontier_and_budget_run_is_deterministic() {
        let p = branch_model();
        for engine in [LpEngineKind::Tableau, LpEngineKind::Revised] {
            let opts = options(engine);
            let interrupted = {
                let _deadline = TestDeadlineGuard::on_check(1);
                solve_controlled(&p, &opts, &SearchControl::default(), None).unwrap()
            };
            assert_eq!(interrupted.report.status, SolveStatus::LimitReached);
            assert_eq!(interrupted.report.stats.nodes, 1);
            assert_eq!(interrupted.report.stats.propagation.accepted_rows, 2);
            assert_eq!(interrupted.report.stats.propagation.calls, 1);
            assert_eq!(interrupted.report.objective, None);
            assert!(interrupted.bound.unwrap() >= 2.125 - 1e-8);
            let ctrl = SearchControl { determinism: Determinism::Budget, ..SearchControl::default() };
            let a = solve_controlled(&p, &opts, &ctrl, None).unwrap();
            let b = solve_controlled(&p, &opts, &ctrl, None).unwrap();
            assert_eq!(a.report.status, SolveStatus::Optimal);
            assert_eq!(a.report.objective, Some(2.125));
            assert_eq!(a.bound, b.bound);
            assert_eq!(a.report.values, b.report.values);
            assert_eq!(a.report.stats.nodes, b.report.stats.nodes);
            assert_eq!(a.report.stats.lp_iterations, b.report.stats.lp_iterations);
            assert_eq!(format!("{:?}", a.report.stats.propagation), format!("{:?}", b.report.stats.propagation));
        }
    }
}
