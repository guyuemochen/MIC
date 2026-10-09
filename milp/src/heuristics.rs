//! 基础启发式（只供应 incumbent，不影响正确性）（WP3 实现）。
//!
//! - [`round_and_check`]：整数变量取最近整数并夹回界，经与 integrality 联动
//!   的容差阶梯复核后评估目标。
//! - [`dive`]：fix-and-dive 潜水，逐步固定“最接近整数的非整整数变量”。

use crate::model::{Problem, RelOp, Tolerances, VarType};
use crate::simplex::{Bounds, LpStatus, Simplex};

/// 启发式得到的可行整数解。
#[derive(Debug, Clone)]
pub struct Incumbent {
    /// 目标值（按取整后的值评估，原始 sense + obj_offset 口径）。
    pub objective: f64,
    /// 结构变量取值（整数变量已取整）。
    pub values: Vec<f64>,
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
    simplex: &mut Simplex,
    bounds: &Bounds,
    lp_x: &[f64],
    max_fixes: usize,
    tol: &Tolerances,
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
        let sol = simplex.resolve(&dive_bounds);
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
                    let retry = simplex.resolve(&dive_bounds);
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
            LpStatus::Unbounded | LpStatus::IterationLimit | LpStatus::Infeasible => {
                // 修复失败：退化为对失败前点的取整复核。
                return round_and_check(problem, bounds, &cur_x, tol);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Sense;

    fn tol() -> Tolerances {
        Tolerances::default()
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
}
