//! 保守预处理与还原（WP3 实现）。
//!
//! 规则仅四类且可证明安全，每轮严格按 R1→R2→R3→R4 顺序施加，禁止重排
//! （反例：整数变量的非整数单点区间 [2.3,2.3] 若先被 R2 固定列规则处理，
//! 会把整数变量固定在非整数值 2.3，终检报 `Err(Numerical)` 而非正确的
//! Infeasible；先 R1 取整则 lb←3 > ub←2 得空界，Simplex 判 Infeasible）：
//! - R1：整数变量（Integer/Binary）界取整 lb←ceil(lb−tol.feasibility)、
//!   ub←floor(ub+tol.feasibility)；取整后 lb>ub 时保留空界（Simplex 将判
//!   Infeasible），不报错也不放宽。
//! - R2：固定列 |ub−lb| ≤ tol.feasibility（空界列不固定）：代入各行右端与
//!   目标（常数移入 obj_offset）后删列。经 R1 后整数变量的单点区间必为整数
//!   值或已空界，因此整数固定值精确。
//! - R3：空行（所有 |a_j| ≤ tol.zero）：满足（判定式与 Simplex 构建期近零行
//!   完全一致：Le rhs ≥ −feas、Ge rhs ≤ feas、Eq |rhs| ≤ feas）则删行；
//!   违反则保留该行——Simplex 构建期近零行判定会置 infeasible_flag 判
//!   Infeasible。
//! - R4：单系数行 a·x_j (op) rhs：换算为 x_j 的界（a<0 时翻转 op，Eq 双侧），
//!   与现界取交集（只收紧）后删行；整数变量交集后再按 R1 取整。
//!
//! 最多迭代 5 轮，无 reductions 即停。数值防御：代入/换算产生非有限右端、
//! 目标常数溢出或无穷收紧界时返回 `Err(MilpError::Numerical)`。

use crate::error::{MilpError, Result};
use crate::model::{Constraint, LinearExpr, Problem, RelOp, Tolerances, VarType, Variable};

/// 预处理轮数上限（轮内无 reductions 即提前停止）。
const MAX_PRESOLVE_ROUNDS: usize = 5;

/// postsolve 还原映射。
#[derive(Debug, Clone)]
pub struct PostsolveMap {
    /// 被固定删除的 (原始变量下标, 固定值)。
    pub fixed: Vec<(usize, f64)>,
    /// 被删除的原始约束下标。
    pub dropped_rows: Vec<usize>,
    /// 缩减问题变量下标 → 原始变量下标。
    pub var_map: Vec<usize>,
}

/// 预处理工作行（原始变量下标空间；alive 标记是否仍存续）。
struct WorkRow {
    /// 该行对应的原始约束下标。
    origin: usize,
    /// 稀疏项（原始变量下标，升序、无零系数——继承 LinearExpr 不变量）。
    terms: Vec<(usize, f64)>,
    /// 关系运算符。
    op: RelOp,
    /// 右端（随 R2 固定列代入更新）。
    rhs: f64,
    /// 是否仍在缩减问题中。
    alive: bool,
}

/// 变量类型是否要求取整（Integer/Binary）。
fn is_integer_type(var_type: VarType) -> bool {
    matches!(var_type, VarType::Integer | VarType::Binary)
}

/// 把 −0.0 规范化为 +0.0（ceil/floor 在零附近可能产生负零，仅影响显示）。
fn normalize_zero(v: f64) -> f64 {
    if v == 0.0 {
        0.0
    } else {
        v
    }
}

/// 保守预处理：绝不切掉整数可行最优；返回 (缩减问题, 还原映射)。
///
/// 输入须已通过 [`Problem::validate`]（有限系数/右端/有限 lb、ub 有限或 +∞）。
/// 空界（lb>ub）会原样保留在缩减问题里，由 Simplex 判 Infeasible——预处理
/// 自身不直接下不可行结论（与 R1/R3 的“交给 Simplex”语义一致）。
pub fn presolve(p: &Problem, tol: &Tolerances) -> Result<(Problem, PostsolveMap)> {
    let n = p.vars.len();
    let mut lb: Vec<f64> = p.vars.iter().map(|v| v.lb).collect();
    let mut ub: Vec<f64> = p.vars.iter().map(|v| v.ub).collect();
    let mut var_alive = vec![true; n];
    let mut obj_terms: Vec<(usize, f64)> = p.obj.terms.clone();
    let mut obj_offset = p.obj_offset;
    let mut rows: Vec<WorkRow> = p
        .constraints
        .iter()
        .enumerate()
        .map(|(origin, c)| WorkRow {
            origin,
            terms: c.expr.terms.clone(),
            op: c.op,
            rhs: c.rhs,
            alive: true,
        })
        .collect();
    let mut fixed: Vec<(usize, f64)> = Vec::new();
    let mut dropped_rows: Vec<usize> = Vec::new();

    for _round in 0..MAX_PRESOLVE_ROUNDS {
        let mut reductions = 0usize;

        // —— R1：整数变量界取整（每轮首先施加，禁止与 R2 换序）——
        for j in 0..n {
            if !var_alive[j] || !is_integer_type(p.vars[j].var_type) {
                continue;
            }
            let new_lb = normalize_zero((lb[j] - tol.feasibility).ceil());
            let new_ub = normalize_zero((ub[j] + tol.feasibility).floor()); // ub=+∞ 时保持 +∞
            if new_lb != lb[j] || new_ub != ub[j] {
                reductions += 1;
            }
            // 取整后 lb>ub 保留空界（Simplex 判 Infeasible）。
            lb[j] = new_lb;
            ub[j] = new_ub;
        }

        // —— R2：固定列（|ub−lb| ≤ feas；空界列不固定）——
        for j in 0..n {
            if !var_alive[j] || lb[j] > ub[j] {
                continue;
            }
            if (ub[j] - lb[j]).abs() > tol.feasibility {
                continue;
            }
            let value = normalize_zero(0.5 * (lb[j] + ub[j]));
            if !value.is_finite() {
                return Err(MilpError::Numerical(format!(
                    "变量 {j} 固定值非有限（lb = {}，ub = {}）",
                    lb[j], ub[j]
                )));
            }
            // 代入目标：常数移入 offset。
            if let Some(pos) = obj_terms.iter().position(|&(t, _)| t == j) {
                obj_offset += obj_terms[pos].1 * value;
                if !obj_offset.is_finite() {
                    return Err(MilpError::Numerical(format!(
                        "固定变量 {j} 代入目标后 obj_offset 溢出"
                    )));
                }
                obj_terms.remove(pos);
            }
            // 代入各行右端并删项。
            for row in rows.iter_mut() {
                if !row.alive {
                    continue;
                }
                if let Some(pos) = row.terms.iter().position(|&(t, _)| t == j) {
                    row.rhs -= row.terms[pos].1 * value;
                    if !row.rhs.is_finite() {
                        return Err(MilpError::Numerical(format!(
                            "固定变量 {j} 代入行 {} 后右端溢出",
                            row.origin
                        )));
                    }
                    row.terms.remove(pos);
                }
            }
            var_alive[j] = false;
            fixed.push((j, value));
            reductions += 1;
        }

        // —— R3：空行（所有 |a_j| ≤ tol.zero）——
        for row in rows.iter_mut() {
            if !row.alive {
                continue;
            }
            if !row.terms.iter().all(|&(_, a)| a.abs() <= tol.zero) {
                continue;
            }
            // 判定式与 Simplex 构建期近零行逐字一致，保证两边结论不漂移。
            let ok = match row.op {
                RelOp::Le => row.rhs >= -tol.feasibility,
                RelOp::Ge => row.rhs <= tol.feasibility,
                RelOp::Eq => row.rhs.abs() <= tol.feasibility,
            };
            if ok {
                row.alive = false;
                dropped_rows.push(row.origin);
                reductions += 1;
            }
            // 违反的空行保留：Simplex 构建期会置 infeasible_flag 判 Infeasible。
        }

        // —— R4：单系数行 → 变量界（交集只收紧）后删行 ——
        for row in rows.iter_mut() {
            if !row.alive || row.terms.len() != 1 {
                continue;
            }
            let (j, a) = row.terms[0];
            if !var_alive[j] {
                continue; // 防御：固定列的项已在 R2 移除，不应出现
            }
            let r = row.rhs;
            let (cand_lb, cand_ub): (Option<f64>, Option<f64>) = match row.op {
                RelOp::Le if a > 0.0 => (None, Some(r / a)),
                RelOp::Le => (Some(r / a), None),
                RelOp::Ge if a > 0.0 => (Some(r / a), None),
                RelOp::Ge => (None, Some(r / a)),
                RelOp::Eq => (Some(r / a), Some(r / a)),
            };
            if let Some(v) = cand_lb {
                if v == f64::INFINITY {
                    return Err(MilpError::Numerical(format!(
                        "行 {} 换算出无穷下界（右端 {r} / 系数 {a} 溢出）",
                        row.origin
                    )));
                }
                if v.is_finite() && v > lb[j] {
                    lb[j] = v;
                    reductions += 1;
                }
            }
            if let Some(v) = cand_ub {
                if v == f64::NEG_INFINITY {
                    return Err(MilpError::Numerical(format!(
                        "行 {} 换算出负无穷上界（右端 {r} / 系数 {a} 溢出）",
                        row.origin
                    )));
                }
                if v.is_finite() && v < ub[j] {
                    ub[j] = v;
                    reductions += 1;
                }
            }
            // 整数变量交集后再按 R1 取整（只收紧；空界保持空界）。
            if is_integer_type(p.vars[j].var_type) {
                lb[j] = normalize_zero((lb[j] - tol.feasibility).ceil());
                ub[j] = normalize_zero((ub[j] + tol.feasibility).floor());
            }
            row.alive = false;
            dropped_rows.push(row.origin);
            reductions += 1;
        }

        if reductions == 0 {
            break;
        }
    }

    // —— 构建缩减问题（存活变量重编号；行保持原始相对顺序）——
    let mut new_index = vec![usize::MAX; n];
    let mut var_map: Vec<usize> = Vec::new();
    let mut vars: Vec<Variable> = Vec::new();
    for j in 0..n {
        if !var_alive[j] {
            continue;
        }
        new_index[j] = var_map.len();
        var_map.push(j);
        vars.push(Variable {
            name: p.vars[j].name.clone(),
            var_type: p.vars[j].var_type,
            lb: normalize_zero(lb[j]),
            ub: normalize_zero(ub[j]),
        });
    }
    let obj = LinearExpr {
        terms: obj_terms
            .iter()
            .filter_map(|&(j, c)| {
                if var_alive[j] {
                    Some((new_index[j], c))
                } else {
                    None
                }
            })
            .collect(),
    };
    let constraints: Vec<Constraint> = rows
        .iter()
        .filter(|row| row.alive)
        .map(|row| Constraint {
            name: p.constraints[row.origin].name.clone(),
            expr: LinearExpr {
                terms: row.terms.iter().map(|&(j, c)| (new_index[j], c)).collect(),
            },
            op: row.op,
            rhs: row.rhs,
        })
        .collect();
    let reduced = Problem {
        name: p.name.clone(),
        sense: p.sense,
        vars,
        obj,
        obj_offset,
        constraints,
    };
    Ok((
        reduced,
        PostsolveMap {
            fixed,
            dropped_rows,
            var_map,
        },
    ))
}

/// 还原到原始变量空间；整数变量值 |x−round| ≤ tol.integrality 时替换为
/// round 值。
///
/// 说明：本函数签名不含问题模型，无法区分整数列，故对全部“距整数不超过
/// tol.integrality”的分量统一取整（含恰在整数附近的连续分量）。该替换偏差
/// ≤ tol.integrality，处于容差阶梯（复核预算 ~1e-6×‖row‖₁）之内，且取整后
/// 目标由上层按最终值重算，值与目标仍自洽。维度不符的缺失分量以 NaN 示警
/// （上层终检会拒绝），不 panic。
pub fn postsolve(
    map: &PostsolveMap,
    reduced_x: &[f64],
    num_orig_vars: usize,
    tol: &Tolerances,
) -> Vec<f64> {
    let mut x = vec![f64::NAN; num_orig_vars];
    for (reduced_idx, &orig) in map.var_map.iter().enumerate() {
        if orig < num_orig_vars {
            x[orig] = reduced_x.get(reduced_idx).copied().unwrap_or(f64::NAN);
        }
    }
    for &(orig, value) in &map.fixed {
        if orig < num_orig_vars {
            x[orig] = value;
        }
    }
    for v in x.iter_mut() {
        if v.is_finite() {
            let r = v.round();
            if (r - *v).abs() <= tol.integrality {
                *v = r;
            }
        }
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Sense, Tolerances, VarType};
    use crate::simplex::{Bounds, LpStatus, Simplex};

    fn tol() -> Tolerances {
        Tolerances::default()
    }

    fn assert_close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 1e-7 * (1.0 + b.abs()),
            "数值不符: {a} vs {b}"
        );
    }

    /// 对原始问题逐约束复核（violation ≤ 1e-6×(1+|rhs|+‖row‖₁)）。
    fn original_feasible(p: &Problem, x: &[f64]) -> bool {
        for c in &p.constraints {
            let lhs = c.expr.eval(x);
            let norm1: f64 = c.expr.terms.iter().map(|&(_, a)| a.abs()).sum();
            let allow = 1e-6 * (1.0 + c.rhs.abs() + norm1);
            let ok = match c.op {
                RelOp::Le => lhs <= c.rhs + allow,
                RelOp::Ge => lhs >= c.rhs - allow,
                RelOp::Eq => (lhs - c.rhs).abs() <= allow,
            };
            if !ok {
                return false;
            }
        }
        true
    }

    /// 固定列 / 空行 / 单系数行三类化简同时命中，且 postsolve 还原后对原始
    /// 约束可行、目标一致。
    #[test]
    fn presolve_hits_fixed_empty_single_and_restores() {
        let mut p = Problem::new("hit", Sense::Maximize);
        let x0 = p.add_var("x0", VarType::Continuous).unwrap();
        let x1 = p.add_var("x1", VarType::Integer).unwrap();
        let x2 = p.add_var("x2", VarType::Continuous).unwrap();
        let x3 = p.add_var("x3", VarType::Integer).unwrap();
        p.set_var_bounds(x0, 0.0, 10.0).unwrap();
        p.set_var_bounds(x1, 0.0, 10.0).unwrap();
        p.set_var_bounds(x2, 1.5, 1.5).unwrap(); // 固定列（连续）
        p.set_objective(vec![(x0, 3.0), (x1, 2.0), (x2, 1.0), (x3, 4.0)], 0.0)
            .unwrap();
        p.add_constraint(
            Some("r0"),
            vec![(x0, 1.0), (x1, 1.0), (x2, 1.0), (x3, 1.0)],
            RelOp::Le,
            10.0,
        )
        .unwrap();
        p.add_constraint(Some("r1"), vec![(x0, 2.0), (x1, -1.0)], RelOp::Ge, 1.0)
            .unwrap();
        // 近零行（|a| = 1e-13 ≤ zero）：满足 → 删除。
        p.add_constraint(Some("r2"), vec![(x0, 1e-13)], RelOp::Le, 5.0)
            .unwrap();
        // 单系数行一对：x3 的界经交集成为 [2,2] → 次轮 R2 固定删除。
        p.add_constraint(Some("r3"), vec![(x3, 1.0)], RelOp::Ge, 2.0)
            .unwrap();
        p.add_constraint(Some("r4"), vec![(x3, 1.0)], RelOp::Le, 2.0)
            .unwrap();

        let (q, map) = presolve(&p, &tol()).unwrap();
        // 结构断言。
        assert_eq!(q.vars.len(), 2);
        assert_eq!(map.var_map, vec![x0, x1]);
        assert_eq!(map.fixed, vec![(x2, 1.5), (x3, 2.0)]);
        assert_eq!(map.dropped_rows, vec![2, 3, 4]);
        assert_eq!(q.constraints.len(), 2);
        assert_eq!(
            q.constraints[0].expr.terms,
            vec![(0, 1.0), (1, 1.0)],
            "r0 应为 x0 + x1 ≤ 6.5"
        );
        assert_close(q.constraints[0].rhs, 6.5);
        assert_eq!(q.constraints[1].expr.terms, vec![(0, 2.0), (1, -1.0)]);
        assert_eq!(q.constraints[1].op, RelOp::Ge);
        assert_close(q.constraints[1].rhs, 1.0);
        assert_eq!(q.obj.terms, vec![(0, 3.0), (1, 2.0)]);
        assert_close(q.obj_offset, 1.5 * 1.0 + 2.0 * 4.0);

        // 求解缩减问题并还原。
        let mut s = Simplex::new(&q, tol()).unwrap();
        let sol = s.solve(&Bounds::from_problem(&q));
        assert_eq!(sol.status, LpStatus::Optimal);
        // 缩减 LP 最优 (6.5, 0)：目标 19.5 + offset 9.5 = 29。
        assert_close(sol.objective, 29.0);
        let full = postsolve(&map, &sol.x, p.vars.len(), &tol());
        assert_eq!(full.len(), 4);
        assert!(original_feasible(&p, &full), "还原解须对原始约束可行");
        let obj_full = p.obj.eval(&full) + p.obj_offset;
        assert_close(obj_full, sol.objective);
        assert_close(full[x3], 2.0);
    }

    /// 整数变量非整数单点区间（z>=2.3 且 z<=2.3）：经 R1→R4 得空界而非固定在
    /// 2.3，Simplex 判 Infeasible（轮内次序 R1 在前的契约反例）。
    #[test]
    fn presolve_integer_singleton_yields_infeasible_not_fixed() {
        let mut p = Problem::new("singleton", Sense::Minimize);
        let z = p.add_var("z", VarType::Integer).unwrap();
        p.set_objective(vec![(z, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(z, 1.0)], RelOp::Ge, 2.3)
            .unwrap();
        p.add_constraint(None, vec![(z, 1.0)], RelOp::Le, 2.3)
            .unwrap();
        let (q, map) = presolve(&p, &tol()).unwrap();
        assert_eq!(q.vars.len(), 1, "空界列不得被 R2 固定删除");
        assert!(
            q.vars[0].lb > q.vars[0].ub,
            "整数单点区间 [2.3,2.3] 经 R1 应为空界"
        );
        assert_close(q.vars[0].lb, 3.0);
        assert_close(q.vars[0].ub, 2.0);
        assert!(map.fixed.is_empty(), "不得把整数变量固定在非整数值 2.3");
        assert_eq!(map.dropped_rows, vec![0, 1]);
        let mut s = Simplex::new(&q, tol()).unwrap();
        let sol = s.solve(&Bounds::from_problem(&q));
        assert_eq!(sol.status, LpStatus::Infeasible);
    }

    /// 违反的空行必须保留（交给 Simplex 判 Infeasible），不得删行。
    #[test]
    fn presolve_violated_empty_row_kept_for_simplex() {
        let mut p = Problem::new("empty", Sense::Minimize);
        let x = p.add_var("x", VarType::Continuous).unwrap();
        p.set_objective(vec![(x, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![], RelOp::Ge, 0.5).unwrap(); // 0 ≥ 0.5 违反
        p.add_constraint(None, vec![(x, 1.0)], RelOp::Le, 3.0)
            .unwrap();
        let (q, map) = presolve(&p, &tol()).unwrap();
        // 行 1（x ≤ 3）是单系数行，经 R4 化为界后删除是正确行为；
        // 行 0（空 Ge 行，违反）必须保留。
        assert_eq!(map.dropped_rows, vec![1], "只允许删单系数行");
        assert_eq!(q.constraints.len(), 1);
        assert!(q.constraints[0].expr.terms.is_empty());
        assert_eq!(q.constraints[0].op, RelOp::Ge);
        let mut s = Simplex::new(&q, tol()).unwrap();
        assert_eq!(
            s.solve(&Bounds::from_problem(&q)).status,
            LpStatus::Infeasible
        );
    }

    /// R1 只取整整数变量界，连续变量界不动；Binary 默认界不受影响。
    #[test]
    fn presolve_rounds_integer_bounds_only() {
        let mut p = Problem::new("round", Sense::Minimize);
        let i = p.add_var("i", VarType::Integer).unwrap();
        let c = p.add_var("c", VarType::Continuous).unwrap();
        let b = p.add_var("b", VarType::Binary).unwrap();
        p.set_var_bounds(i, 0.4, 3.7).unwrap();
        p.set_var_bounds(c, 0.4, 3.7).unwrap();
        p.set_objective(vec![(i, 1.0), (c, 1.0), (b, 1.0)], 0.0)
            .unwrap();
        p.add_constraint(None, vec![(i, 1.0), (c, 1.0)], RelOp::Ge, 0.1)
            .unwrap();
        let (q, _) = presolve(&p, &tol()).unwrap();
        assert_eq!((q.vars[0].lb, q.vars[0].ub), (1.0, 3.0));
        assert_eq!((q.vars[1].lb, q.vars[1].ub), (0.4, 3.7), "连续界不动");
        assert_eq!((q.vars[2].lb, q.vars[2].ub), (0.0, 1.0));
    }

    /// 无可化简的问题：原样返回（恒等映射、无固定/删行）。
    #[test]
    fn presolve_noop_on_clean_problem() {
        let mut p = Problem::new("clean", Sense::Minimize);
        let x = p.add_var("x", VarType::Continuous).unwrap();
        let y = p.add_var("y", VarType::Integer).unwrap();
        p.set_var_bounds(x, -1.25, 4.75).unwrap();
        p.set_var_bounds(y, -2.0, 7.0).unwrap();
        p.set_objective(vec![(x, 1.0), (y, 3.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(x, 1.0), (y, 1.0)], RelOp::Le, 5.0)
            .unwrap();
        let (q, map) = presolve(&p, &tol()).unwrap();
        assert_eq!(map.var_map, vec![0, 1]);
        assert!(map.fixed.is_empty());
        assert!(map.dropped_rows.is_empty());
        assert_eq!(q.constraints.len(), 1);
        assert_eq!(q.vars.len(), 2);
        assert_eq!((q.vars[0].lb, q.vars[0].ub), (-1.25, 4.75));
        assert_eq!((q.vars[1].lb, q.vars[1].ub), (-2.0, 7.0));
        assert_eq!(q.obj.terms, vec![(0, 1.0), (1, 3.0)]);
        assert_close(q.obj_offset, 0.0);
    }

    /// postsolve：近整数分量取整、分数分量与固定值不动、缺失分量 NaN 示警。
    #[test]
    fn postsolve_snaps_near_integers_only() {
        let map = PostsolveMap {
            fixed: vec![(2, 1.5)],
            dropped_rows: vec![],
            var_map: vec![0, 1],
        };
        let x = postsolve(&map, &[0.9999999995, 2.5], 3, &tol());
        assert_eq!(x, vec![1.0, 2.5, 1.5]);
        // 距整数超过 integrality 的分量不动。
        let x = postsolve(&map, &[1.0000049, 2.5], 3, &tol());
        assert!((x[0] - 1.0000049).abs() < 1e-12);
        // 维度不符：缺失分量 NaN（示警而非 panic）。
        let bad = postsolve(&map, &[1.0], 3, &tol());
        assert!(bad[1].is_nan());
    }
}
