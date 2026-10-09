//! 统一求解编排与终检（WP4 契约，本文件为 WP4 实现）。
//!
//! 流程：`validate` →（可选）presolve → branch_and_bound(缩减问题) → postsolve
//! → 安全网终检 → 报告。
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
//! - presolve 关闭时直接在原始问题上求解（`values` 天然同维）；
//!   开启时 `presolve_reductions` = 被固定删除的列数 + 被删除的行数。
//! - 只有 B&B 带回取值（Optimal 或有 incumbent 的 LimitReached）才走
//!   postsolve 与终检；Infeasible/Unbounded/无 incumbent 的 LimitReached
//!   报告 `objective`/`values` 均为 `None`。
//! - 终检容差为契约固定安全网值（相对 1e-6、整数 1e-5），宽于 WP3 复核
//!   阶梯，正常路径不会被误伤；超差说明数值已失真，宁可报错也不输出错解。

use crate::branch_bound::{self, BbOptions, BbStats, SearchControl, SearchHooks};
use crate::error::{MilpError, Result};
use crate::model::{Problem, RelOp, SolveStatus, Tolerances, VarType};
use crate::presolve;
use crate::search_ctrl::gap_percent;

/// 终检约束违反相对容差（契约固定：violation ≤ 1e-6×(1+|rhs|+‖row‖₁)）。
const FINAL_VIOLATION_REL: f64 = 1e-6;
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
}

impl Default for SolveOptions {
    /// 默认值：Tolerances::default() + BbOptions::default() + use_presolve=true。
    fn default() -> Self {
        SolveOptions {
            tol: Tolerances::default(),
            bb: BbOptions::default(),
            use_presolve: true,
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
    problem.validate()?;
    let num_orig_vars = problem.vars.len();

    // （可选）presolve → 缩减问题 + 还原映射。
    let (work, map, reductions) = if opts.use_presolve {
        let (reduced, map) = presolve::presolve(problem, &opts.tol)?;
        let reductions = map.fixed.len() + map.dropped_rows.len();
        (reduced, Some(map), reductions)
    } else {
        (problem.clone(), None, 0)
    };

    let bb = branch_bound::branch_and_bound_controlled(&work, &opts.bb, ctrl, opts.tol, hooks)?;

    // postsolve → 原始变量空间；B&B 未带回取值则不带解返回。
    let (objective, values) = match bb.best_values {
        Some(reduced_x) => {
            let mut x = match &map {
                Some(m) => presolve::postsolve(m, &reduced_x, num_orig_vars, &opts.tol),
                None => reduced_x,
            };
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
fn final_check(problem: &Problem, x: &mut [f64]) -> Result<()> {
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
}
