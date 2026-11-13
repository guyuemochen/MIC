//! 只执行原生产求解一次，捕获首个真实第一阶段；不额外定价、重构或求解。
use milp::branch_bound::SearchControl;
use milp::instance::load_instance;
use milp::simplex::phase1_progress_diagnostics;
use milp::simplex_sparse::LpEngineKind;
use milp::solver::{solve_controlled, SolveOptions};
use std::error::Error;
use std::io::{BufWriter, Write};
use std::path::Path;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() < 2 {
        return Err(
            "用法：phase1_progress_probe 输入.milp 输出.json [秒数] [dense|revised]".into(),
        );
    }
    let problem = load_instance(Path::new(&args[0]))?;
    let mut options = SolveOptions::default();
    options.bb.lp_engine = match args.get(3).map(String::as_str).unwrap_or("revised") {
        "dense" => LpEngineKind::Tableau,
        "revised" => LpEngineKind::Revised,
        _ => return Err("未知LP引擎".into()),
    };
    options.bb.time_limit_ms = args
        .get(2)
        .map(|s| s.parse::<u64>())
        .transpose()?
        .unwrap_or(30)
        .checked_mul(1000)
        .ok_or("时限溢出")?;
    phase1_progress_diagnostics::reset();
    let result = solve_controlled(
        &problem,
        &options,
        &SearchControl {
            threads: 1,
            adaptive_heuristics: true,
            ..SearchControl::default()
        },
        None,
    );
    let capture = phase1_progress_diagnostics::take();
    let mut out = BufWriter::new(std::fs::File::create(&args[1])?);
    // 外部原文件与engine_input_lp可能经历预处理/缩放，绝不能混用变量坐标。
    write!(out,"{{\"float_encoding\":\"ieee754_binary64_hex_bits\",\"external_input_not_engine_coordinates\":")?;
    phase1_progress_diagnostics::write_problem_json(&problem, &mut out)?;
    match result {
        Ok(report) => write!(
            out,
            ",\"status\":\"{:?}\",\"nodes\":{},\"lp_iterations\":{}",
            report.report.status, report.report.stats.nodes, report.report.stats.lp_iterations
        )?,
        Err(_) => write!(out, ",\"solver_error\":true")?,
    }
    write!(out, ",\"capture\":")?;
    if let Some(capture) = capture {
        capture.write_json(&mut out)?;
    } else {
        write!(out, "null")?;
    }
    writeln!(out, "}}")?;
    out.flush()?;
    Ok(())
}
