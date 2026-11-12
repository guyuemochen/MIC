//! 首次LU失败证据捕获。默认只执行原生产求解；--recheck额外做离线式基分解，
//! 不把复核结果送回求解器，不据模型已知标签判不可行。
use milp::branch_bound::SearchControl;
use milp::instance::load_instance;
use milp::simplex_sparse::{lu_failure_diagnostics, LpEngineKind};
use milp::solver::{solve_controlled, SolveOptions};
use std::error::Error;
use std::io::{BufWriter, Write};
use std::path::Path;

fn string(out: &mut impl Write, value: &str) -> std::io::Result<()> {
    write!(out, "\"")?;
    for c in value.chars() {
        match c {
            '"' => write!(out, "\\\"")?,
            '\\' => write!(out, "\\\\")?,
            '\u{0}'..='\u{1f}' => write!(out, "\\u{:04x}", u32::from(c))?,
            _ => write!(out, "{c}")?,
        }
    }
    write!(out, "\"")
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() < 2 {
        return Err("用法：lu_failure_probe 输入.milp 输出.json [--recheck] [--replay-history] [--time-limit-secs 秒]".into());
    }
    let mut opts = SolveOptions::default();
    opts.bb.lp_engine = LpEngineKind::Revised;
    opts.bb.time_limit_ms = 30_000;
    let mut recheck = false;
    let mut replay_history = false;
    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--recheck" => recheck = true,
            "--replay-history" => replay_history = true,
            "--time-limit-secs" => {
                i += 1;
                opts.bb.time_limit_ms = args
                    .get(i)
                    .ok_or("时限缺值")?
                    .parse::<u64>()?
                    .checked_mul(1000)
                    .ok_or("时限溢出")?;
            }
            value => return Err(format!("未知选项{value}").into()),
        }
        i += 1;
    }
    let problem = load_instance(Path::new(&args[0]))?;
    lu_failure_diagnostics::reset();
    let result = solve_controlled(
        &problem,
        &opts,
        &SearchControl {
            threads: 1,
            adaptive_heuristics: true,
            ..SearchControl::default()
        },
        None,
    );
    let failure = lu_failure_diagnostics::take();
    let mut out = BufWriter::new(std::fs::File::create(&args[1])?);
    write!(out, "{{\"schema\":1,\"problem\":")?;
    string(&mut out, &problem.name)?;
    match result {
        Ok(report) => write!(
            out,
            ",\"status\":\"{:?}\",\"nodes\":{},\"lp_iterations\":{}",
            report.report.status, report.report.stats.nodes, report.report.stats.lp_iterations
        )?,
        Err(error) => {
            write!(out, ",\"solver_error\":")?;
            string(&mut out, &error.to_string())?;
        }
    }
    write!(out, ",\"failure\":")?;
    if let Some(failure) = failure {
        failure.write_json(&mut out)?;
        if recheck {
            // take已关闭捕获；复核仅用保存的矩阵，不访问或修复生产引擎。
            let checks = failure.recheck();
            write!(out, ",\"fresh_checks\":[")?;
            for (i, check) in checks.iter().enumerate() {
                if i > 0 {
                    write!(out, ",")?;
                }
                check.write_json(&mut out)?;
            }
            write!(out, "]")?;
        }
        if replay_history {
            write!(out, ",\"history_checks\":[")?;
            for (i, (after, check)) in failure.recheck_history().iter().enumerate() {
                if i > 0 {
                    write!(out, ",")?;
                }
                write!(out, "{{\"after_pivots\":{after},\"check\":")?;
                check.write_json(&mut out)?;
                write!(out, "}}")?;
            }
            write!(out, "]")?;
        }
    } else {
        write!(out, "null")?;
    }
    writeln!(out, "}}")?;
    out.flush()?;
    Ok(())
}
