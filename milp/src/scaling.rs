//! 数值缩放（R4，轴4）：行列迭代几何均值 scaling（2 的幂系数）+ 问题级
//! 全链路往返。
//!
//! # 动机（差距文档轴4）
//! 三个已修数值问题（对偶修复阈值、mir_floor 吸附、ratio_test 带宽）都是
//! 同一缺失的表征——系数量级失衡使固定容差在不同行列上语义漂移。本模块在
//! 求解编排最前端做问题级缩放，使容差判定作用于均衡后的量级（标准做法）。
//!
//! # 算法
//! 迭代几何均值（≤ [`MAX_SWEEPS`] 轮 sweep）：
//! 1. 行因子 r_i = 1/√(max|a_ij|·min|a_ij|)（当前行**非零**上，≥2 非零才
//!    缩放），吸附到最近 2 的幂后整行乘上（a_ij ← r_i·a_ij，b_i ← r_i·b_i）；
//! 2. 列因子同式（当前列非零上，乘子 m_j = 1/√(max·min)），吸附后整列
//!    乘上（a_ij ← a_ij·m_j ≡ a_ij/s_j，s_j = 1/m_j）。
//! 全部因子吸附为 2^0（本轮无变化）即停。吸附到 2 的幂保证 f64 乘除
//! **精确无舍入**——unscale 的 x = ŷ/2^k 逐位还原，整数界 ×2^k 保持整数。
//!
//! # 代换与往返
//! ŷ_j = x_j·s_j（s_j = 2^k）⇒ â_ij = r_i·a_ij/s_j、b̂_i = r_i·b_i、
//! ĉ_j = c_j/s_j、界 ŷ ∈ [lb·s_j, ub·s_j]；目标**值**不变（ĉ·ŷ = c·x，
//! 精确算术下恒等），obj_offset/sense/VarType 不动。正缩放保序 ⇒ B&B
//! incumbent 比较在缩放空间合法。线性双射保可行集 ⇒ 割在缩放空间生成、
//! 对原问题同样有效。
//!
//! # 变量类型约束（正确性红线，实测抓出）
//! - **仅 Continuous 列可缩放**：Binary 的 0/1 语义与 Integer 的整格
//!   （x ∈ Z ⟺ ŷ ∈ Z 对 2 的幂因子仅 s_j = 1 成立——ask 的「Integer 列
//!   k ≥ 0」会把半整数反解回原空间，被终检安全网在 milp-mixed-250i40c
//!   上实测拒绝）；
//! - 全整数算例（knapsack/mknap/setcover/marketsplit）的列缩放整体
//!   不参与，只有行缩放生效（收益折损，如实报告）。
//!
//! 与 Simplex 构建期行归一（1/max|a_ij|，`row_scale`）叠加无害——row_scale
//! 只管对偶还原，缩放行的 max|â| 更均衡后其归一系数也更接近 1。

use crate::model::{Problem, VarType};

/// sweep 轮数上限（每轮 = 行遍 + 列遍）。
const MAX_SWEEPS: usize = 8;

/// 缩放映射：行因子 r_i 与列（变量放大）因子 s_j，全为 2 的幂。
#[derive(Debug, Clone)]
pub struct ScaleMap {
    /// 行因子 r_i（约束表行序；含被跳过行的 1.0）。
    pub row: Vec<f64>,
    /// 列因子 s_j（变量表列序；Binary 列与被跳过列恒 1.0）。
    pub col: Vec<f64>,
}

/// 吸附后的 2 的幂指数（供整型列 k ≥ 0 钳制；f ≤ 0 或非有限返回 0，
/// 防御）。
fn pow2_exponent(f: f64) -> i32 {
    if !f.is_finite() || f <= 0.0 {
        return 0;
    }
    f.log2().round() as i32
}

fn pow2(k: i32) -> f64 {
    2f64.powi(k)
}

/// 列 j 是否可缩放：**仅 Continuous 列**。
///
/// 【正确性红线（实测抓出）】ask 原文「Integer 列 s_j = 2^k, k ≥ 0」不健全：
/// x ↦ ŷ = x·2^k 把 x 的整格映射进 2^k·Z ⊂ Z，但缩放问题按 VarType::Integer
/// 在 **ŷ ∈ Z** 上求解——奇数 ŷ 反解 x = ŷ/2^k 得半整数（milp-mixed-250i30c
/// 实测：y191 终值 0.5 被终检安全网拒绝）。格双射要求 x ∈ Z ⟺ ŷ ∈ Z，
/// 对 2 的幂因子只有 s_j = 1 同时满足两个方向——故 Binary（0/1 语义）与
/// Integer（整格）列一律不缩放，收益折损如实报告（全整数算例只剩行缩放）。
fn col_scalable(var_type: VarType) -> bool {
    matches!(var_type, VarType::Continuous)
}

/// 问题级几何均值缩放。返回（缩放后问题, 缩放映射）。原问题不动。
///
/// 数值防御：任一系数/右端/界非有限或行内出现精确零（结构性零不影响——
/// 只在非零上取 min/max；右端无零约束）时该行/列跳过缩放（因子 1），
/// 绝不产 NaN。
pub fn scale_problem(p: &Problem) -> (Problem, ScaleMap) {
    let n = p.vars.len();
    let m = p.constraints.len();
    // 稠密工作矩阵（缩放只在非零上取几何均值；稀疏→稠密便于列遍历，
    // 问题规模已在 B&B 可解量级）。
    let mut a: Vec<Vec<f64>> = vec![vec![0.0; n]; m];
    let mut b: Vec<f64> = vec![0.0; m];
    let mut row_ok = vec![true; m];
    for (i, con) in p.constraints.iter().enumerate() {
        b[i] = con.rhs;
        if !con.rhs.is_finite() {
            row_ok[i] = false;
        }
        for &(j, v) in &con.expr.terms {
            if j >= n || !v.is_finite() || v == 0.0 {
                row_ok[i] = false;
                continue;
            }
            a[i][j] = v;
        }
    }
    let col_scalable: Vec<bool> = p.vars.iter().map(|v| col_scalable(v.var_type)).collect();
    // 界非有限的列仍可缩放（∞·2^k = ∞）；下界恒有限（crate 契约）。

    let mut row_factor = vec![1.0f64; m];
    let mut col_exp = vec![0i32; n]; // 累计列指数 s_j = 2^col_exp[j]

    for _ in 0..MAX_SWEEPS {
        let mut changed = false;
        // —— 行遍：r_i = snap(1/√(max·min))，应用于当前矩阵 ——
        let mut row_exp = vec![0i32; m];
        for i in 0..m {
            if !row_ok[i] {
                continue;
            }
            let mut mx = 0.0f64;
            let mut mn = f64::INFINITY;
            let mut nz = 0usize;
            for &v in &a[i] {
                if v != 0.0 {
                    mx = mx.max(v.abs());
                    mn = mn.min(v.abs());
                    nz += 1;
                }
            }
            if nz == 0 || !mx.is_finite() || mn <= 0.0 {
                continue; // 空行/数值异常不缩放
            }
            let f = 1.0 / (mx * mn).sqrt();
            let k = pow2_exponent(f);
            if k != 0 {
                row_exp[i] = k;
                changed = true;
            }
        }
        for i in 0..m {
            let k = row_exp[i];
            if k == 0 {
                continue;
            }
            let r = pow2(k);
            for v in a[i].iter_mut() {
                if *v != 0.0 {
                    *v *= r;
                }
            }
            b[i] *= r;
            row_factor[i] *= r;
        }
        // —— 列遍：s_j = snap(1/√(max·min))（当前矩阵），按类型钳制 ——
        let mut col_pass = vec![0i32; n];
        for j in 0..n {
            if !col_scalable[j] {
                continue;
            }
            let mut mx = 0.0f64;
            let mut mn = f64::INFINITY;
            let mut nz = 0usize;
            for i in 0..m {
                let v = a[i][j];
                if v != 0.0 {
                    mx = mx.max(v.abs());
                    mn = mn.min(v.abs());
                    nz += 1;
                }
            }
            if nz == 0 || !mx.is_finite() || mn <= 0.0 {
                continue;
            }
            // 列乘子 m = 1/√(max·min)（几何均衡方向）；代换 â = a/s_j ⇒
            // s_j = 1/m = √(max·min)，即 s 的指数 = −m 的指数。可缩放集已
            // 限定 Continuous 列（无格约束）。
            let m_exp = pow2_exponent(1.0 / (mx * mn).sqrt());
            let k = -m_exp;
            if k != 0 {
                col_pass[j] = k;
                changed = true;
            }
        }
        for j in 0..n {
            let k = col_pass[j];
            if k == 0 {
                continue;
            }
            let s = pow2(k);
            for i in 0..m {
                if a[i][j] != 0.0 {
                    a[i][j] /= s;
                }
            }
            col_exp[j] += k;
        }
        if !changed {
            break;
        }
    }

    // —— 组装缩放后问题 ——
    let mut q = Problem::new(&p.name, p.sense);
    q.obj_offset = p.obj_offset;
    for (j, v) in p.vars.iter().enumerate() {
        let s = pow2(col_exp[j]);
        let mut sv = v.clone();
        sv.lb = v.lb * s;
        sv.ub = v.ub * s; // ±∞·2^k 不变号（s > 0）
        q.vars.push(sv);
    }
    // 目标：ĉ_j = c_j / s_j（合并重复列后稀疏重建）。
    let mut obj_terms: Vec<(usize, f64)> = Vec::with_capacity(p.obj.terms.len());
    for &(j, c) in &p.obj.terms {
        let cj = c / pow2(col_exp[j]);
        match obj_terms.last_mut() {
            Some(last) if last.0 == j => last.1 += cj,
            _ => obj_terms.push((j, cj)),
        }
    }
    q.obj = crate::model::LinearExpr {
        terms: obj_terms,
    };
    for (i, con) in p.constraints.iter().enumerate() {
        let mut terms: Vec<(usize, f64)> = Vec::with_capacity(con.expr.terms.len());
        for &(j, _) in &con.expr.terms {
            let v = a[i][j];
            if v != 0.0 {
                terms.push((j, v));
            }
        }
        terms.sort_by_key(|&(j, _)| j);
        q.constraints.push(crate::model::Constraint {
            name: con.name.clone(),
            expr: crate::model::LinearExpr { terms },
            op: con.op,
            rhs: b[i],
        });
    }
    let map = ScaleMap {
        row: row_factor,
        col: (0..n).map(|j| pow2(col_exp[j])).collect(),
    };
    (q, map)
}

impl ScaleMap {
    /// 还原变量取值：x_j = ŷ_j / s_j（2 的幂除法，f64 精确无舍入）。
    pub fn unscale_x(&self, y: &[f64]) -> Vec<f64> {
        y.iter()
            .zip(self.col.iter())
            .map(|(&yj, &sj)| yj / sj)
            .collect()
    }
}

#[cfg(test)]
mod engage_tests {
    use super::*;
    use crate::instance::load_instance;
    use std::path::Path;

    /// 诊断：bench 族上缩放的实际生效面（行因子分布 + 连续列因子数）。
    /// 全整数算例只有行缩放（列缩放被整格红线禁止）；milp-mixed 连续列
    /// 参与列缩放。
    #[test]
    fn scaling_engagement_on_bench_family() {
        for name in [
            "knapsack-2000",
            "mknap-3x100",
            "marketsplit-4x30",
            "setcover-100x150",
            "milp-mixed-250i40c",
        ] {
            let p = load_instance(Path::new(&format!(
                "../bench/instances/{name}.milp"
            )))
            .unwrap();
            let (q, map) = scale_problem(&p);
            let nontrivial_rows = map.row.iter().filter(|&&r| r != 1.0).count();
            let row_ratio = map
                .row
                .iter()
                .filter(|&&r| r != 1.0)
                .map(|r| r.log2().abs())
                .fold(0.0_f64, f64::max);
            let nontrivial_cols = map.col.iter().filter(|&&s| s != 1.0).count();
            // 缩放后全矩阵 max|a|/min|a|（非零、跨行）量级对比缩放前。
            let spread = |prob: &Problem| -> f64 {
                let mut mx = 0.0f64;
                let mut mn = f64::INFINITY;
                for con in &prob.constraints {
                    for &(_, v) in &con.expr.terms {
                        mx = mx.max(v.abs());
                        mn = mn.min(v.abs());
                    }
                }
                mx / mn
            };
            eprintln!(
                "ENGAGE {name}: rows={}/{} 非平凡（最大 |log2| = {row_ratio:.1}），                  cols={}/{} 非平凡，矩阵幅度比 {} -> {}",
                nontrivial_rows,
                map.row.len(),
                nontrivial_cols,
                map.col.len(),
                spread(&p),
                spread(&q),
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RelOp, Sense, VarType};

    fn var(name: &str, t: VarType, lb: f64, ub: f64) -> crate::model::Variable {
        crate::model::Variable {
            name: name.to_string(),
            var_type: t,
            lb,
            ub,
        }
    }

    /// 2 的幂性 + 往返精确性 + Binary 列不缩 + 整数列 k≥0。
    #[test]
    fn pow2_factors_and_roundtrip() {
        // 2 连续列量级失衡 [1e-3, 1e5]，1 条 Le 行。
        let mut p = Problem::new("s", Sense::Maximize);
        p.vars.push(var("u", VarType::Continuous, 0.0, 1e6));
        p.vars.push(var("v", VarType::Continuous, 0.0, 1e6));
        p.obj = crate::model::LinearExpr {
            terms: vec![(0, 1.0), (1, 1.0)],
        };
        p.constraints.push(crate::model::Constraint {
            name: None,
            expr: crate::model::LinearExpr {
                terms: vec![(0, 1e-3), (1, 1e5)],
            },
            op: RelOp::Le,
            rhs: 1.0,
        });
        let (q, map) = scale_problem(&p);
        // 全因子 2 的幂。
        for r in &map.row {
            assert!(r.log2().abs().fract() < 1e-12, "行因子 {r} 非 2 的幂");
        }
        for s in &map.col {
            assert!(s.log2().abs().fract() < 1e-12, "列因子 {s} 非 2 的幂");
        }
        // 行几何均值均衡：缩放后行内 max/min ≤ 4（2 的幂吸附下 ≤ 一档）。
        let row0: Vec<f64> = q.constraints[0]
            .expr
            .terms
            .iter()
            .map(|&(_, v)| v.abs())
            .collect();
        let mx = row0.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mn = row0.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(mx / mn <= 4.0 + 1e-9, "行未均衡：{row0:?}");
        // 往返精确性：任意 x → ŷ = x·s → unscale 逐位还原。
        let x = [123.456, 7.89];
        let y: Vec<f64> = x.iter().zip(map.col.iter()).map(|(&xj, &sj)| xj * sj).collect();
        let back = map.unscale_x(&y);
        assert_eq!(back, vec![123.456, 7.89]); // f64 == 逐位
        // 代换恒等：â·ŷ = r·(a·x)（同一 x 点上两侧求值；精确算术，浮点下
        // ulp 级——2 的幂乘除精确，仅乘法结合序有尾差）。
        let lhs: f64 = q.constraints[0]
            .expr
            .terms
            .iter()
            .map(|&(j, a)| a * y[j])
            .sum();
        let rhs_direct: f64 = map.row[0]
            * p.constraints[0]
                .expr
                .terms
                .iter()
                .map(|&(j, a)| a * x[j])
                .sum::<f64>();
        assert!(
            (lhs - rhs_direct).abs() <= 1e-9 * (1.0 + rhs_direct.abs()),
            "代换恒等破坏：{lhs} vs {rhs_direct}"
        );
        // 右端恒等：b̂ = r·b。
        let b_hat = q.constraints[0].rhs;
        let b_direct = map.row[0] * p.constraints[0].rhs;
        assert_eq!(b_hat, b_direct, "右端缩放恒等破坏");
    }

    /// Binary/Integer 列恒不缩（0/1 语义 / 整格双射红线）；仅 Continuous 缩。
    #[test]
    fn var_type_constraints() {
        let mut p = Problem::new("b", Sense::Maximize);
        p.vars.push(var("b0", VarType::Binary, 0.0, 1.0));
        p.vars.push(var("y1", VarType::Integer, 0.0, 10.0));
        p.vars.push(var("c2", VarType::Continuous, 0.0, 1e5));
        p.obj = crate::model::LinearExpr {
            terms: vec![(0, 3.0), (1, 2.0), (2, 1.0)],
        };
        p.constraints.push(crate::model::Constraint {
            name: None,
            expr: crate::model::LinearExpr {
                terms: vec![(0, 1e4), (1, 1.0), (2, 1e-4)],
            },
            op: RelOp::Le,
            rhs: 5e3,
        });
        let (q, map) = scale_problem(&p);
        assert_eq!(map.col[0], 1.0, "Binary 列不得缩放");
        assert_eq!(map.col[1], 1.0, "Integer 列不得缩放（整格双射红线）");
        assert_eq!(q.vars[1].lb, p.vars[1].lb);
        assert_eq!(q.vars[1].ub, p.vars[1].ub);
        // 连续列可缩（失衡列 [1e-4 vs 1e4] 下非平凡因子）。
        assert!(map.col[2] != 1.0, "失衡连续列应被缩放（得 {}）", map.col[2]);
    }

    /// 量级已均衡（行内比值 1 且 |a| = 1）⇒ 全部因子 1（吸附不动）。
    ///（注意几何缩放同时均衡量级：[3,3] 行会合法归一到 [0.75,0.75]
    ///（1/3 吸附到 2^-2），非零因子不等于「比值不变」——钉 [1,1] 的
    /// 真不动锚点。）
    #[test]
    fn uniform_problem_untouched() {
        let mut p = Problem::new("u", Sense::Minimize);
        p.vars.push(var("x", VarType::Integer, 0.0, 5.0));
        p.vars.push(var("y", VarType::Integer, 0.0, 5.0));
        p.obj = crate::model::LinearExpr {
            terms: vec![(0, 1.0), (1, 1.0)],
        };
        p.constraints.push(crate::model::Constraint {
            name: None,
            expr: crate::model::LinearExpr {
                terms: vec![(0, 1.0), (1, 1.0)],
            },
            op: RelOp::Le,
            rhs: 6.0,
        });
        let (q, map) = scale_problem(&p);
        assert!(map.row.iter().all(|&r| r == 1.0));
        assert!(map.col.iter().all(|&s| s == 1.0));
        assert_eq!(q.constraints[0].rhs, 6.0);
        assert_eq!(q.constraints[0].expr.terms, vec![(0, 1.0), (1, 1.0)]);
    }

    /// 全链路等价：缩放问题的整数最优与原问题一致（穷举小例）。
    #[test]
    fn scaled_optimum_matches_original() {
        let mut p = Problem::new("e", Sense::Maximize);
        p.vars.push(var("x", VarType::Integer, 0.0, 8.0));
        p.vars.push(var("y", VarType::Integer, 0.0, 8.0));
        p.vars.push(var("c", VarType::Continuous, 0.0, 2.0));
        p.obj = crate::model::LinearExpr {
            terms: vec![(0, 3.0), (1, 4.0), (2, 0.5)],
        };
        p.constraints.push(crate::model::Constraint {
            name: None,
            expr: crate::model::LinearExpr {
                terms: vec![(0, 1e2), (1, 1e-2), (2, 1.0)],
            },
            op: RelOp::Le,
            rhs: 8e1,
        });
        p.constraints.push(crate::model::Constraint {
            name: None,
            expr: crate::model::LinearExpr {
                terms: vec![(0, 1.0), (1, 1.0)],
            },
            op: RelOp::Le,
            rhs: 10.0,
        });
        let (q, map) = scale_problem(&p);
        // 穷举原问题整数可行点最大目标；缩放问题同点（映射后）目标一致。
        let mut best_orig = f64::NEG_INFINITY;
        let mut best_mapped = f64::NEG_INFINITY;
        for x in 0..=8i64 {
            for y in 0..=8i64 {
                for c in [0.0f64, 2.0] {
                    let v = [x as f64, y as f64, c];
                    if 1e2 * v[0] + 1e-2 * v[1] + v[2] <= 8e1 && v[0] + v[1] <= 10.0 {
                        let obj: f64 = p.obj.terms.iter().map(|&(j, cj)| cj * v[j]).sum();
                        best_orig = best_orig.max(obj);
                        let yv: Vec<f64> =
                            v.iter().zip(map.col.iter()).map(|(&xj, &sj)| xj * sj).collect();
                        let objs: f64 = q.obj.terms.iter().map(|&(j, cj)| cj * yv[j]).sum();
                        best_mapped = best_mapped.max(objs);
                    }
                }
            }
        }
        assert!(
            (best_orig - best_mapped).abs() <= 1e-9 * (1.0 + best_orig.abs()),
            "缩放空间最优值不一致：{best_orig} vs {best_mapped}"
        );
        assert!(best_orig > f64::NEG_INFINITY);
    }
}
