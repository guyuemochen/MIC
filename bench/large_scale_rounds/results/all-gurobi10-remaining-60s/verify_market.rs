// 以独立的两个排序数组给单位成本二元等式计算下界；每个真解必在标量投影上匹配。
// 每半的高 59 位存投影余数，低 5 位存取一数。投影碰撞只能让所得下界变松。
use milp::model::{RelOp, Sense, VarType};
use std::{env, path::Path};

fn enumerate(weights: &[u64], target: Option<u64>) -> Vec<u64> {
    let mut values = Vec::with_capacity(1 << weights.len());
    values.push(target.unwrap_or(0) << 5);
    for &w in weights {
        let step = if target.is_some() {
            w.wrapping_neg()
        } else {
            w
        };
        for i in 0..values.len() {
            values.push(values[i].wrapping_add((step << 5) | 1));
        }
    }
    values.sort_unstable();
    values
}

fn main() {
    let p = milp::instance::load_instance(Path::new(&env::args().nth(1).unwrap())).unwrap();
    let n = p.vars.len();
    assert!(n <= 60 && p.sense == Sense::Minimize && p.obj_offset == 0.);
    assert!(p
        .vars
        .iter()
        .all(|v| v.var_type == VarType::Binary && v.lb == 0. && v.ub == 1.));
    assert!(p.obj.terms.len() == n && p.obj.terms.iter().all(|&(_, c)| c == 1.));
    let mut weights = vec![0_u64; n];
    let mut target = 0_u64;
    let mut salt = 231_u64;
    for row in &p.constraints {
        assert!(row.op == RelOp::Eq && row.rhs == row.rhs.round() && row.rhs.abs() < 1e6);
        salt = salt
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        target = target.wrapping_add(salt.wrapping_mul(row.rhs as i64 as u64));
        for &(j, a) in &row.expr.terms {
            assert!(a == a.round() && a.abs() < 1e6);
            weights[j] = weights[j].wrapping_add(salt.wrapping_mul(a as i64 as u64));
        }
    }
    let left = enumerate(&weights[..n / 2], None);
    let right = enumerate(&weights[n / 2..], Some(target));
    let (mut i, mut j, mut lower) = (0, 0, None);
    while i < left.len() && j < right.len() {
        let a = left[i] >> 5;
        let b = right[j] >> 5;
        if a < b {
            i += 1;
        } else if a > b {
            j += 1;
        } else {
            // 排序让同一投影下取一数最少的组合排在最前。
            let count = (left[i] & 31) + (right[j] & 31);
            lower = Some(lower.map_or(count, |v: u64| v.min(count)));
            while i < left.len() && left[i] >> 5 == a {
                i += 1;
            }
            while j < right.len() && right[j] >> 5 == b {
                j += 1;
            }
        }
    }
    println!("{{\"left_subsets\":{},\"right_subsets\":{},\"projection_bits\":59,\"certified_lower_bound\":{}}}", left.len(),right.len(), lower.unwrap());
}
