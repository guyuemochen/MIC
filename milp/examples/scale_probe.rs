//! 缩放容量探针：模型生成、缩放和完整求解分别计时；不改变求解策略。
//! 用法：scale_probe scale <维数> <每行项数> <continuous|integer|mixed>
//! 或 scale_probe solve <维数>，或 scale_probe file <模型路径>。

use milp::instance::load_instance;
use milp::model::{Constraint, LinearExpr, Problem, RelOp, Sense, VarType, Variable};
use milp::scaling::{scale_problem, ScaleMap};
use milp::solver::{solve, SolveOptions};
use std::error::Error;
use std::path::Path;
use std::time::Instant;

fn generated(n: usize, width: usize, kind: &str, diagonal: bool) -> Problem {
    let mut p = Problem::new("稀疏容量验证", Sense::Maximize);
    for j in 0..n {
        let var_type = match kind {
            "integer" => VarType::Integer,
            "mixed" if j % 3 == 0 => VarType::Integer,
            _ => VarType::Continuous,
        };
        p.vars.push(Variable {
            name: format!("x{j}"),
            var_type,
            lb: 0.0,
            ub: 8.0,
        });
        p.obj.terms.push((j, -1.0));
    }
    for i in 0..n {
        let mut terms = Vec::with_capacity(width);
        for k in 0..width {
            let exponent = ((i * 13 + k * 7) % 25) as i32 - 12;
            let sign = if !diagonal && (i + k) % 3 == 0 {
                -1.0
            } else {
                1.0
            };
            terms.push(((i + k) % n, sign * 2.0_f64.powi(exponent)));
        }
        terms.sort_by_key(|&(j, _)| j);
        p.constraints.push(Constraint {
            name: None,
            expr: LinearExpr { terms },
            op: RelOp::Le,
            rhs: 4.0,
        });
    }
    p
}

fn mix(hash: &mut u64, value: u64) {
    for byte in value.to_le_bytes() {
        *hash ^= u64::from(byte);
        *hash = hash.wrapping_mul(0x100000001b3);
    }
}

/// 按位摘要覆盖系数、界、右端、目标和映射，用于同一工具链新旧对照。
fn fingerprint(p: &Problem, map: &ScaleMap) -> u64 {
    let mut h = 0xcbf29ce484222325;
    mix(&mut h, p.vars.len() as u64);
    mix(&mut h, p.constraints.len() as u64);
    mix(&mut h, p.obj_offset.to_bits());
    for v in &p.vars {
        mix(&mut h, v.var_type as u64);
        mix(&mut h, v.lb.to_bits());
        mix(&mut h, v.ub.to_bits());
    }
    for &(j, value) in &p.obj.terms {
        mix(&mut h, j as u64);
        mix(&mut h, value.to_bits());
    }
    for c in &p.constraints {
        mix(&mut h, c.op as u64);
        mix(&mut h, c.rhs.to_bits());
        mix(&mut h, c.expr.terms.len() as u64);
        for &(j, value) in &c.expr.terms {
            mix(&mut h, j as u64);
            mix(&mut h, value.to_bits());
        }
    }
    for &value in map.row.iter().chain(&map.col) {
        mix(&mut h, value.to_bits());
    }
    h
}

fn run() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let start = Instant::now();
    let (p, complete) = match args.first().map(String::as_str) {
        Some("scale") if args.len() == 4 => {
            let n: usize = args[1].parse()?;
            let width: usize = args[2].parse()?;
            if n == 0
                || width == 0
                || width > n
                || !["continuous", "integer", "mixed"].contains(&args[3].as_str())
            {
                return Err("维数、项数或变量类型非法".into());
            }
            (generated(n, width, &args[3], false), false)
        }
        Some("solve") if args.len() == 2 => {
            let n: usize = args[1].parse()?;
            if n == 0 {
                return Err("维数必须大于零".into());
            }
            // 独立单系数行由预处理消去，已知最优解全零；只验证完整入口容量。
            (generated(n, 1, "continuous", true), true)
        }
        Some("file") if args.len() == 2 => (load_instance(Path::new(&args[1]))?, false),
        _ => {
            return Err(
                "用法：scale_probe scale N W continuous|integer|mixed；solve N；file PATH".into(),
            )
        }
    };
    p.validate()?;
    let setup_secs = start.elapsed().as_secs_f64();
    let nnz: usize = p.constraints.iter().map(|c| c.expr.terms.len()).sum();
    let operation = Instant::now();
    if complete {
        let report = solve(&p, &SolveOptions::default())?;
        let operation_secs = operation.elapsed().as_secs_f64();
        let objective = report
            .objective
            .filter(|v| v.is_finite())
            .map_or_else(|| "null".to_string(), |v| format!("{v:?}"));
        let mut digest = 0xcbf29ce484222325;
        if let Some(values) = report.values.as_ref() {
            for &value in values {
                mix(&mut digest, value.to_bits());
            }
        }
        println!("{{\"mode\":\"solve\",\"rows\":{},\"cols\":{},\"nnz\":{nnz},\"setup_secs\":{setup_secs},\"operation_secs\":{operation_secs},\"total_secs\":{},\"status\":\"{:?}\",\"objective\":{objective},\"nodes\":{},\"lp_iterations\":{},\"digest\":\"{digest:016x}\"}}",
            p.constraints.len(), p.vars.len(), start.elapsed().as_secs_f64(),
            report.status, report.stats.nodes, report.stats.lp_iterations);
    } else {
        let (q, map) = scale_problem(&p);
        let operation_secs = operation.elapsed().as_secs_f64();
        q.validate()?;
        let digest = fingerprint(&q, &map);
        println!("{{\"mode\":\"scale\",\"rows\":{},\"cols\":{},\"nnz\":{nnz},\"setup_secs\":{setup_secs},\"operation_secs\":{operation_secs},\"total_secs\":{},\"digest\":\"{digest:016x}\"}}",
            p.constraints.len(), p.vars.len(), start.elapsed().as_secs_f64());
    }
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}
