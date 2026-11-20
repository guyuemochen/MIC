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
    // 工作副本只保存输入非零项；列统计按行升序累积，不再展开 m×n 矩阵。
    // 合法模型的项已严格升序，因此每个系数仍按原轮序先乘行因子、再除列因子。
    let mut a: Vec<Vec<(usize, f64)>> = Vec::with_capacity(m);
    let mut b: Vec<f64> = vec![0.0; m];
    let mut row_ok = vec![true; m];
    let mut canonical_rows = vec![true; m];
    for (i, con) in p.constraints.iter().enumerate() {
        b[i] = con.rhs;
        if !con.rhs.is_finite() {
            row_ok[i] = false;
        }
        let mut row = Vec::with_capacity(con.expr.terms.len());
        let mut previous = None;
        for &(j, v) in &con.expr.terms {
            if previous.is_some_and(|prev| prev >= j) {
                canonical_rows[i] = false;
            }
            previous = Some(j);
            if j >= n || !v.is_finite() || v == 0.0 {
                row_ok[i] = false;
                canonical_rows[i] = false;
                continue;
            }
            row.push((j, v));
        }
        if !canonical_rows[i] {
            // 求解入口会先 validate；此分支仅保留直接调用时旧实现的有限项
            // 覆盖规则：乱序排序，同列最后一个有效值覆盖前值，不扩张模型契约。
            row.sort_by_key(|&(j, _)| j);
            row.dedup_by(|later, earlier| {
                if later.0 == earlier.0 {
                    earlier.1 = later.1;
                    true
                } else {
                    false
                }
            });
        }
        a.push(row);
    }
    let col_scalable: Vec<bool> = p.vars.iter().map(|v| col_scalable(v.var_type)).collect();
    // 界非有限的列仍可缩放（∞·2^k = ∞）；下界恒有限（crate 契约）。

    let mut row_factor = vec![1.0f64; m];
    let mut col_exp = vec![0i32; n]; // 累计列指数 s_j = 2^col_exp[j]
    let mut row_exp = vec![0i32; m];
    // 全整数模型原本就不做列缩放，不为这类模型增加逐项列统计。
    let has_scalable_columns = col_scalable.iter().any(|&scalable| scalable);
    let column_workspace_len = if has_scalable_columns { n } else { 0 };
    let mut col_max = vec![0.0f64; column_workspace_len];
    let mut col_min = vec![f64::INFINITY; column_workspace_len];
    let mut col_nz = vec![0usize; column_workspace_len];
    let mut col_divisor = vec![1.0f64; column_workspace_len];

    for _ in 0..MAX_SWEEPS {
        let mut changed = false;
        // —— 行遍：r_i = snap(1/√(max·min))，应用于当前矩阵 ——
        row_exp.fill(0);
        for i in 0..m {
            if !row_ok[i] {
                continue;
            }
            let mut mx = 0.0f64;
            let mut mn = f64::INFINITY;
            let mut nz = 0usize;
            for &(_, v) in &a[i] {
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
            for (_, v) in a[i].iter_mut() {
                if *v != 0.0 {
                    *v *= r;
                }
            }
            b[i] *= r;
            row_factor[i] *= r;
        }
        if has_scalable_columns {
            // —— 列遍：各列仍按原始行号升序取 min/max；零项不参与 ——
            col_max.fill(0.0);
            col_min.fill(f64::INFINITY);
            col_nz.fill(0);
            col_divisor.fill(1.0);
            for row in &a {
                for &(j, v) in row {
                    if col_scalable[j] && v != 0.0 {
                        col_max[j] = col_max[j].max(v.abs());
                        col_min[j] = col_min[j].min(v.abs());
                        col_nz[j] += 1;
                    }
                }
            }
            for j in 0..n {
                if !col_scalable[j] {
                    continue;
                }
                let mx = col_max[j];
                let mn = col_min[j];
                if col_nz[j] == 0 || !mx.is_finite() || mn <= 0.0 {
                    continue;
                }
                // 列乘子 m = 1/√(max·min)（几何均衡方向）；代换 â = a/s_j ⇒
                // s_j = 1/m = √(max·min)，即 s 的指数 = −m 的指数。可缩放集已
                // 限定 Continuous 列（无格约束）。
                let m_exp = pow2_exponent(1.0 / (mx * mn).sqrt());
                let k = -m_exp;
                if k != 0 {
                    col_divisor[j] = pow2(k);
                    col_exp[j] += k;
                    changed = true;
                }
            }
            // 同一轮各列互不依赖，按稀疏行访问不会改变任一系数的运算顺序。
            for row in &mut a {
                for (j, v) in row {
                    let s = col_divisor[*j];
                    if s != 1.0 && *v != 0.0 {
                        *v /= s;
                    }
                }
            }
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
    q.obj = crate::model::LinearExpr { terms: obj_terms };
    for (i, (con, mut row)) in p.constraints.iter().zip(a).enumerate() {
        let terms = if canonical_rows[i] {
            // 合法输入的列序未变，直接转移非零工作项，避免复制整行。
            row.retain(|&(_, v)| v != 0.0);
            row
        } else {
            // 仅兼容未经 validate 的直接调用：保留旧输出的重复项和排序，
            // 越界下标也维持原来的拒绝方式，不将非法模型静默修成可用模型。
            let mut terms = Vec::with_capacity(con.expr.terms.len());
            for &(j, _) in &con.expr.terms {
                assert!(j < n, "缩放约束项变量下标越界");
                let v = row
                    .binary_search_by_key(&j, |&(col, _)| col)
                    .map(|pos| row[pos].1)
                    .unwrap_or(0.0);
                if v != 0.0 {
                    terms.push((j, v));
                }
            }
            terms.sort_by_key(|&(j, _)| j);
            terms
        };
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
            let p = load_instance(Path::new(&format!("../bench/instances/{name}.milp"))).unwrap();
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

    // 改动前的稠密算法只供小模型逐位对照；生产路径不得使用。
    fn dense_reference_scale_problem(p: &Problem) -> (Problem, ScaleMap) {
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
        q.obj = crate::model::LinearExpr { terms: obj_terms };
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
        let y: Vec<f64> = x
            .iter()
            .zip(map.col.iter())
            .map(|(&xj, &sj)| xj * sj)
            .collect();
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
                        let yv: Vec<f64> = v
                            .iter()
                            .zip(map.col.iter())
                            .map(|(&xj, &sj)| xj * sj)
                            .collect();
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

    fn assert_float_bits(actual: f64, expected: f64, context: &str) {
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "{context} 的位模式不一致：{actual:?} / {expected:?}"
        );
    }

    fn assert_terms_bits(actual: &[(usize, f64)], expected: &[(usize, f64)], context: &str) {
        assert_eq!(actual.len(), expected.len(), "{context} 的项数不一致");
        for (&(aj, av), &(ej, ev)) in actual.iter().zip(expected) {
            assert_eq!(aj, ej, "{context} 的列号不一致");
            assert_float_bits(av, ev, context);
        }
    }

    fn assert_scaled_bits(p: &Problem) {
        let (actual, actual_map) = scale_problem(p);
        let (expected, expected_map) = dense_reference_scale_problem(p);
        let context = p.name.as_str();
        assert_eq!(actual.name, expected.name);
        assert_eq!(actual.sense, expected.sense);
        assert_float_bits(actual.obj_offset, expected.obj_offset, context);
        assert_terms_bits(&actual.obj.terms, &expected.obj.terms, context);
        assert_eq!(actual.vars.len(), expected.vars.len());
        for (av, ev) in actual.vars.iter().zip(&expected.vars) {
            assert_eq!(av.name, ev.name);
            assert_eq!(av.var_type, ev.var_type);
            assert_float_bits(av.lb, ev.lb, context);
            assert_float_bits(av.ub, ev.ub, context);
        }
        assert_eq!(actual.constraints.len(), expected.constraints.len());
        for (ac, ec) in actual.constraints.iter().zip(&expected.constraints) {
            assert_eq!(ac.name, ec.name);
            assert_eq!(ac.op, ec.op);
            assert_float_bits(ac.rhs, ec.rhs, context);
            assert_terms_bits(&ac.expr.terms, &ec.expr.terms, context);
        }
        assert_eq!(actual_map.row.len(), expected_map.row.len());
        assert_eq!(actual_map.col.len(), expected_map.col.len());
        for (&a, &e) in actual_map.row.iter().zip(&expected_map.row) {
            assert_float_bits(a, e, context);
        }
        for (&a, &e) in actual_map.col.iter().zip(&expected_map.col) {
            assert_float_bits(a, e, context);
        }
    }

    struct ScalingRng(u64);

    impl ScalingRng {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1);
            self.0
        }

        fn coefficient(&mut self) -> f64 {
            let sign = if self.next() & (1 << 20) == 0 {
                -1.0
            } else {
                1.0
            };
            let exponent = (self.next() % 161) as i32 - 80;
            let mantissa = ((self.next() % 511) + 1) as f64 / 32.0;
            sign * mantissa * pow2(exponent)
        }
    }

    /// 逐位对照旧算法，覆盖宽/窄、空行/列、全整数和混合列、正负系数，
    /// 也检查名称、上下界、目标常数及关系符号未因稀疏化改变。
    #[test]
    fn sparse_scaling_matches_dense_reference_bitwise() {
        for seed in 0..12 {
            let mut rng = ScalingRng(seed + 1);
            for (m, n) in [
                (0, 0),
                (0, 7),
                (7, 0),
                (1, 1),
                (1, 31),
                (31, 1),
                (5, 17),
                (17, 5),
                (23, 31),
                (31, 23),
                (9, 9),
            ] {
                for density in [0, 1, 6, 11] {
                    let mut p = Problem::new(
                        &format!("逐位_{seed}_{m}_{n}_{density}"),
                        if seed % 2 == 0 {
                            Sense::Minimize
                        } else {
                            Sense::Maximize
                        },
                    );
                    p.obj_offset = if seed % 2 == 0 { -0.0 } else { -3.25 };
                    for j in 0..n {
                        let kind = match seed % 3 {
                            0 => {
                                if j % 2 == 0 {
                                    VarType::Binary
                                } else {
                                    VarType::Integer
                                }
                            }
                            1 => VarType::Continuous,
                            _ => match j % 3 {
                                0 => VarType::Binary,
                                1 => VarType::Integer,
                                _ => VarType::Continuous,
                            },
                        };
                        let (lb, ub) = if kind == VarType::Binary {
                            (0.0, 1.0)
                        } else {
                            (
                                if j % 7 == 0 { -0.0 } else { (j % 3) as f64 },
                                if j % 5 == 0 { f64::INFINITY } else { 30.25 },
                            )
                        };
                        p.vars.push(var(&format!("x{j}"), kind, lb, ub));
                        if j % 4 != 0 {
                            p.obj.terms.push((j, rng.coefficient()));
                        }
                    }
                    for i in 0..m {
                        let mut terms = Vec::new();
                        for j in 0..n {
                            if rng.next() % 11 < density {
                                terms.push((j, rng.coefficient()));
                            }
                        }
                        p.constraints.push(crate::model::Constraint {
                            name: if i % 2 == 0 {
                                Some(format!("r{i}"))
                            } else {
                                None
                            },
                            expr: crate::model::LinearExpr { terms },
                            op: match i % 3 {
                                0 => RelOp::Le,
                                1 => RelOp::Ge,
                                _ => RelOp::Eq,
                            },
                            rhs: if i % 7 == 0 { -0.0 } else { rng.coefficient() },
                        });
                    }
                    p.validate().unwrap();
                    assert_scaled_bits(&p);
                }
            }
        }
    }

    /// 乘积上溢/下溢时仍按旧指数防御跳过，不顺手改变数值策略。
    #[test]
    fn sparse_scaling_preserves_extreme_finite_coefficients() {
        let mut p = Problem::new("极端有限系数", Sense::Minimize);
        for j in 0..5 {
            p.vars.push(var(
                &format!("x{j}"),
                VarType::Continuous,
                0.0,
                f64::INFINITY,
            ));
        }
        for (i, values) in [[1e200, -2e200], [1e-200, -2e-200]].iter().enumerate() {
            p.constraints.push(crate::model::Constraint {
                name: None,
                expr: crate::model::LinearExpr {
                    terms: vec![(i * 2, values[0]), (i * 2 + 1, values[1])],
                },
                op: RelOp::Le,
                rhs: -0.0,
            });
        }
        p.validate().unwrap();
        assert_scaled_bits(&p);
    }

    /// 只维持公开函数已有防御，不把未校验的重复/零/非有限项升级成合法模型。
    #[test]
    fn direct_scaling_keeps_existing_invalid_row_behavior() {
        let mut p = Problem::new("未校验行兼容", Sense::Minimize);
        for j in 0..3 {
            p.vars
                .push(var(&format!("x{j}"), VarType::Continuous, 0.0, 10.0));
        }
        p.constraints.push(crate::model::Constraint {
            name: None,
            expr: crate::model::LinearExpr {
                terms: vec![(2, 3.0), (0, 32.0), (2, -8.0), (1, 0.0), (2, f64::NAN)],
            },
            op: RelOp::Ge,
            rhs: 0.5,
        });
        assert!(p.validate().is_err());
        assert_scaled_bits(&p);
        p.constraints[0].expr.terms = vec![(2, 3.0), (0, 32.0), (2, -8.0)];
        assert!(p.validate().is_err());
        assert_scaled_bits(&p);
    }

    /// 十万行列但每行至多一个非零项：旧实现会试图分配约 80 GB，
    /// 本例只跑稀疏实现并按单项缩放的解析结果逐项验收，禁止调用稠密参考。
    #[test]
    fn sparse_scaling_handles_one_hundred_thousand_rows_and_columns() {
        let n = 100_007usize;
        let m = 100_003usize;
        let mut p = Problem::new("十万级稀疏容量", Sense::Minimize);
        p.vars = (0..n)
            .map(|j| var(&format!("x{j}"), VarType::Continuous, 0.0, 16.0))
            .collect();
        p.constraints = (0..m)
            .map(|i| {
                let exponent = (i % 41) as i32 - 20;
                let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                crate::model::Constraint {
                    name: None,
                    expr: crate::model::LinearExpr {
                        terms: if i % 17 == 0 {
                            Vec::new()
                        } else {
                            vec![(i, sign * pow2(exponent))]
                        },
                    },
                    op: RelOp::Le,
                    rhs: 2.0,
                }
            })
            .collect();
        p.validate().unwrap();
        let (q, map) = scale_problem(&p);
        assert_eq!(q.vars.len(), n);
        assert_eq!(q.constraints.len(), m);
        assert_eq!(map.col.len(), n);
        assert_eq!(map.row.len(), m);
        assert!(map.col.iter().all(|&s| s.to_bits() == 1.0f64.to_bits()));
        for (i, con) in q.constraints.iter().enumerate() {
            if i % 17 == 0 {
                assert!(con.expr.terms.is_empty());
                assert_float_bits(map.row[i], 1.0, "空行因子");
                assert_float_bits(con.rhs, 2.0, "空行右端");
            } else {
                let r = pow2(20 - (i % 41) as i32);
                let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                assert_terms_bits(&con.expr.terms, &[(i, sign)], "十万级单项缩放");
                assert_float_bits(map.row[i], r, "十万级行因子");
                assert_float_bits(con.rhs, 2.0 * r, "十万级右端");
            }
        }
        q.validate().unwrap();
    }
}
