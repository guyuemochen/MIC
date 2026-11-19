//! 只运行候选缩减模型的首次根 LP，按原列下标输出报告解。
use milp::branch_bound::BbOptions;
use milp::instance::load_instance;
use milp::model::Tolerances;
use milp::presolve::{postsolve, presolve};
use milp::scaling::scale_problem;
use milp::simplex::{Bounds, LpStatus};
use milp::simplex_sparse::{effective_engine, new_engine};
use std::io::{BufWriter, Error, ErrorKind, Write};
use std::path::PathBuf;
use std::time::{Duration, Instant};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let paths: Vec<PathBuf> = std::env::args_os().skip(1).map(PathBuf::from).collect();
    let [model] = paths.as_slice() else {
        return Err(Error::new(ErrorKind::InvalidInput, "需要一个原始模型路径").into());
    };
    let start = Instant::now();
    let original = load_instance(model)?;
    original.validate()?;
    let tolerance = Tolerances::default();
    let (scaled, scale) = scale_problem(&original);
    let (reduced, mapping) = presolve(&scaled, &tolerance)?;
    if !mapping.aggregations.is_empty() {
        return Err(Error::new(ErrorKind::InvalidData, "二元聚合后的 LP 取值无法逐列还原").into());
    }
    let kind = effective_engine(&reduced, BbOptions::default().lp_engine);
    let mut engine = new_engine(&reduced, tolerance, kind)?;
    engine.set_time_deadline(Some(start + Duration::from_secs(30)));
    let solution = engine.solve(&Bounds::from_problem(&reduced));
    let diagnostics = engine.proof_diagnostics();
    let first_artificial = reduced.num_vars() + engine.column_layout().slack.len();
    let basic_artificials = engine
        .basis_columns()
        .iter()
        .filter(|&&j| j >= first_artificial)
        .count();
    eprintln!(
        "model={} engine={:?} status={:?} objective={:.17e} iterations={} phase1_pending={:?} basic_artificials={} numerical_error={:?} elapsed={:.6}",
        original.name,
        kind,
        solution.status,
        solution.objective,
        engine.iterations(),
        diagnostics.phase1_pending,
        basic_artificials,
        diagnostics.numerical_error,
        start.elapsed().as_secs_f64(),
    );
    if solution.status == LpStatus::Optimal {
        // 禁止对报告点自行取整。两例均无聚合，只恢复固定列和列尺度。
        let exact = Tolerances {
            integrality: 0.0,
            ..tolerance
        };
        let scaled = postsolve(&mapping, &solution.x, original.num_vars(), &exact);
        let values = scale.unscale_x(&scaled);
        let mut out = BufWriter::new(std::io::stdout().lock());
        for (j, (variable, value)) in original.vars.iter().zip(values).enumerate() {
            writeln!(out, "{j}\t{}\t{:.17e}", variable.name, value)?;
        }
    }
    Ok(())
}
