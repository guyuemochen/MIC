//! 单次求解探针与原模型导出。用于基准，不执行额外的根 LP。
use milp::branch_bound::{BranchRule, NodeRule, SearchControl};
use milp::instance::{load_instance, write_instance};
use milp::model::{Constraint, LinearExpr, Problem, RelOp, Sense, SolveStatus, VarType, Variable};
use milp::search_ctrl::{gap_percent, DegradePolicy};
use milp::simplex_sparse::LpEngineKind;
use milp::solver::{solve_controlled, SolveOptions};
use std::error::Error;
use std::io::{self, Write};
use std::path::Path;
use std::time::Instant;

fn number(value: Option<f64>) -> String {
    value
        .filter(|v| v.is_finite())
        .map_or_else(|| "null".into(), |v| format!("{v:?}"))
}

fn string(out: &mut dyn Write, value: &str) -> io::Result<()> {
    write!(out, "\"")?;
    for ch in value.chars() {
        match ch {
            '"' => write!(out, "\\\"")?,
            '\\' => write!(out, "\\\\")?,
            '\u{0}'..='\u{1f}' => write!(out, "\\u{:04x}", u32::from(ch))?,
            _ => write!(out, "{ch}")?,
        }
    }
    write!(out, "\"")
}

fn array<T>(
    out: &mut dyn Write,
    values: &[T],
    mut element: impl FnMut(&mut dyn Write, &T) -> io::Result<()>,
) -> io::Result<()> {
    write!(out, "[")?;
    for (index, value) in values.iter().enumerate() {
        if index != 0 {
            write!(out, ",")?;
        }
        element(out, value)?;
    }
    write!(out, "]")
}

fn export(out: &mut impl Write, problem: &Problem) -> io::Result<()> {
    write!(out, "{{\"name\":")?;
    string(out, &problem.name)?;
    let sense = match problem.sense {
        Sense::Minimize => "min",
        Sense::Maximize => "max",
    };
    write!(
        out,
        ",\"sense\":\"{sense}\",\"offset\":{:?},\"objective\":",
        problem.obj_offset
    )?;
    array(out, &problem.obj.terms, |w, (j, a)| {
        write!(w, "[{j},{a:?}]")
    })?;
    write!(out, ",\"variables\":")?;
    array(out, &problem.vars, |w, var| {
        let kind = match var.var_type {
            VarType::Continuous => "continuous",
            VarType::Integer => "integer",
            VarType::Binary => "binary",
        };
        write!(w, "{{\"name\":")?;
        string(w, &var.name)?;
        write!(
            w,
            ",\"kind\":\"{kind}\",\"lb\":{:?},\"ub\":{}}}",
            var.lb,
            number(Some(var.ub))
        )
    })?;
    write!(out, ",\"rows\":")?;
    array(out, &problem.constraints, |w, row| {
        let op = match row.op {
            RelOp::Le => "<=",
            RelOp::Ge => ">=",
            RelOp::Eq => "=",
        };
        write!(w, "{{\"op\":\"{op}\",\"rhs\":{:?},\"terms\":", row.rhs)?;
        array(w, &row.expr.terms, |w, (j, a)| write!(w, "[{j},{a:?}]"))?;
        write!(w, "}}")
    })?;
    writeln!(out, "}}")
}

fn toggle(arg: &str) -> Result<bool, Box<dyn Error>> {
    match arg {
        "on" => Ok(true),
        "off" => Ok(false),
        _ => Err(format!("期望 on 或 off，收到 {arg}").into()),
    }
}

fn solve(out: &mut impl Write, problem: &Problem, args: &[String]) -> Result<(), Box<dyn Error>> {
    let mut opts = SolveOptions::default();
    let mut ctrl = SearchControl {
        adaptive_heuristics: true,
        ..SearchControl::default()
    };
    let mut pairs = args.chunks_exact(2);
    if !pairs.remainder().is_empty() {
        return Err("求解选项必须成对给出".into());
    }
    for pair in &mut pairs {
        let flag = pair[0].as_str();
        let value = pair[1].as_str();
        match flag {
            "--node-limit" => opts.bb.node_limit = value.parse()?,
            "--time-limit-secs" => {
                opts.bb.time_limit_ms = value.parse::<u64>()?.checked_mul(1000).ok_or("秒数过大")?
            }
            "--threads" => ctrl.threads = value.parse()?,
            "--lp-engine" => {
                opts.bb.lp_engine = match value {
                    "auto" => LpEngineKind::Auto,
                    "tableau" => LpEngineKind::Tableau,
                    "revised" => LpEngineKind::Revised,
                    _ => return Err(format!("未知 LP 引擎 {value}").into()),
                }
            }
            "--node-rule" => {
                opts.bb.node_rule = match value {
                    "dfs" => NodeRule::DepthFirst,
                    "hybrid" => NodeRule::HybridBestBound,
                    _ => return Err(format!("未知节点规则 {value}").into()),
                }
            }
            "--branch-rule" => {
                opts.bb.branch_rule_auto = value == "auto";
                opts.bb.branch_rule = match value {
                    "auto" | "mostfrac" => BranchRule::MostFractional,
                    "pseudo" => BranchRule::Pseudocost,
                    _ => return Err(format!("未知分支规则 {value}").into()),
                };
            }
            "--pc-probes" => opts.bb.pc_reliability_probes = value.parse()?,
            "--probe-reuse" => opts.bb.reuse_probe_evidence = toggle(value)?,
            "--early-stop" => {
                ctrl.degrade = if toggle(value)? {
                    DegradePolicy::default()
                } else {
                    DegradePolicy::disabled()
                }
            }
            "--local-cuts" => opts.bb.local_cuts = toggle(value)?,
            "--root-cuts" => opts.bb.root_cuts = toggle(value)?,
            "--rc-fixing" => opts.bb.root_rc_fixing = toggle(value)?,
            _ => return Err(format!("未知求解选项 {flag}").into()),
        }
    }
    if ctrl.threads == 0 {
        return Err("线程数至少为 1".into());
    }
    let start = Instant::now();
    let result = solve_controlled(problem, &opts, &ctrl, None)?;
    let elapsed = start.elapsed().as_secs_f64();
    let report = &result.report;
    let stats = report.stats;
    let status = match report.status {
        SolveStatus::Optimal => "optimal",
        SolveStatus::Infeasible => "infeasible",
        SolveStatus::Unbounded => "unbounded",
        SolveStatus::LimitReached => "limit_reached",
    };
    let gap = report
        .objective
        .zip(result.bound)
        .and_then(|(objective, bound)| gap_percent(problem.sense, objective, bound));
    let absolute_gap = report
        .objective
        .zip(result.bound)
        .map(|(objective, bound)| (objective - bound).abs());
    write!(out,
        "{{\"status\":\"{status}\",\"raw_status\":\"{status}\",\"elapsed_secs\":{elapsed:?},\"objective\":{},\"bound\":{},\"gap_percent\":{},\"gap_absolute\":{},\"nodes\":{},\"lp_iterations\":{},\"cuts\":{},\"rc_fixed\":{},\"rc_reapplies\":{},\"rc_refixed\":{},\"presolve_reductions\":{},\"auto_degraded\":{},\"gap_terminated\":{},\"values\":",
        number(report.objective), number(result.bound), number(gap), number(absolute_gap), stats.nodes,
        stats.lp_iterations, stats.cuts, stats.rc_fixed, stats.rc_reapplies, stats.rc_refixed,
        report.presolve_reductions, result.auto_degraded, result.gap_terminated,
    )?;
    if let Some(values) = &report.values {
        array(out, values, |w, value| {
            write!(w, "{}", number(Some(*value)))
        })?;
    } else {
        write!(out, "null")?;
    }
    writeln!(out, "}}")?;
    Ok(())
}

fn generate(size: usize, width: usize, kind: &str) -> Result<Problem, Box<dyn Error>> {
    if size == 0
        || width == 0
        || width > size
        || !matches!(kind, "continuous" | "integer" | "mixed")
    {
        return Err("期望 1 ≤ 每行项数 ≤ 维数，类型为 continuous、integer 或 mixed".into());
    }
    let mut problem = Problem::new(&format!("generated-{size}-{width}-{kind}"), Sense::Maximize);
    for j in 0..size {
        let var_type = if kind == "integer" || (kind == "mixed" && j % 3 == 0) {
            VarType::Integer
        } else {
            VarType::Continuous
        };
        problem.vars.push(Variable {
            name: format!("x{j}"),
            var_type,
            lb: 0.0,
            ub: 8.0,
        });
        problem.obj.terms.push((j, (1 + j % 7) as f64));
    }
    for i in 0..size {
        let mut terms: Vec<_> = (0..width)
            .map(|k| ((i + k) % size, (1 + (i + 7 * k) % 11) as f64))
            .collect();
        terms.sort_unstable_by_key(|&(j, _)| j);
        problem.constraints.push(Constraint {
            name: None,
            expr: LinearExpr { terms },
            op: RelOp::Le,
            rhs: (3 * width) as f64,
        });
    }
    Ok(problem)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let mut out = io::BufWriter::new(io::stdout().lock());
    match args.as_slice() {
        [command, path] if command == "export" => export(&mut out, &load_instance(Path::new(path))?)?,
        [command, path, opts @ ..] if command == "solve" => solve(&mut out, &load_instance(Path::new(path))?, opts)?,
        [command, size, width, kind] if command == "generate" => {
            let model = generate(size.parse()?, width.parse()?, kind)?;
            out.write_all(write_instance(&model)?.as_bytes())?;
        }
        _ => return Err("用法：round_probe solve <模型> [选项] | export <模型> | generate <维数> <每行项数> <类型>".into()),
    }
    out.flush()?;
    Ok(())
}
