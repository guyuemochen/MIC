//! 只跑生产缩放和预处理，写出缩减模型及原始列固定记录。
use milp::instance::{load_instance, write_instance};
use milp::model::Tolerances;
use milp::presolve::presolve;
use milp::scaling::scale_problem;
use std::io::{Error, ErrorKind};
use std::path::PathBuf;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let paths: Vec<PathBuf> = std::env::args_os().skip(1).map(PathBuf::from).collect();
    let (input, output) = match paths.as_slice() {
        [input] => (input, None),
        [input, output] => (input, Some(output)),
        _ => {
            return Err(
                Error::new(ErrorKind::InvalidInput, "需要输入模型，可另给缩减输出路径").into(),
            )
        }
    };
    let problem = load_instance(input)?;
    problem.validate()?;
    let scaling_start = Instant::now();
    let (scaled, _) = scale_problem(&problem);
    let elapsed_scaling = scaling_start.elapsed();
    let presolve_start = Instant::now();
    let (reduced, map) = presolve(&scaled, &Tolerances::default())?;
    let elapsed_presolve = presolve_start.elapsed();
    if let Some(output) = output {
        std::fs::write(output, write_instance(&reduced)?)?;
    }
    eprintln!(
        "{} {}列/{}行 → {}列/{}行 固定{} 聚合{} 缩放{:.3}ms 预处理{:.3}ms",
        input.display(),
        scaled.num_vars(),
        scaled.num_constraints(),
        reduced.num_vars(),
        reduced.num_constraints(),
        map.fixed.len(),
        map.aggregations.len(),
        elapsed_scaling.as_secs_f64() * 1000.0,
        elapsed_presolve.as_secs_f64() * 1000.0,
    );
    for (j, value) in map.fixed {
        println!("{}\t{}\t{value}", problem.vars[j].name, j);
    }
    Ok(())
}
