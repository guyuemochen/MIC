//! 默认与诊断特性构建的公开LP定迭代对账；不设墙钟，输出供逐字节比较。
use milp::instance::load_instance;
use milp::model::Tolerances;
use milp::simplex::{Bounds, Simplex};
use milp::simplex_sparse::SparseSimplex;
use std::error::Error;
use std::path::Path;
fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 2 {
        return Err("用法：phase1_trajectory_probe 输入.milp 迭代上限".into());
    }
    let p = load_instance(Path::new(&args[0]))?;
    let cap = args[1].parse::<u64>()?;
    let bounds = Bounds::from_problem(&p);
    #[cfg(feature = "phase1_progress_diagnostics")]
    milp::simplex::phase1_progress_diagnostics::reset();
    let mut dense = Simplex::new(&p, Tolerances::default()).unwrap();
    dense.set_iter_limit_hint(Some(cap));
    let d = dense.solve(&bounds);
    println!(
        "Dense {:?} objective={:016x} iterations={} basis={:?} x={:?}",
        d.status,
        d.objective.to_bits(),
        dense.iterations(),
        dense.basis_columns(),
        d.x.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
    );
    #[cfg(feature = "phase1_progress_diagnostics")]
    {
        milp::simplex::phase1_progress_diagnostics::take();
        milp::simplex::phase1_progress_diagnostics::reset();
    }
    let mut sparse = SparseSimplex::new(&p, Tolerances::default()).unwrap();
    sparse.set_iter_limit_hint(Some(cap));
    let r = sparse.solve(&bounds);
    println!(
        "Revised {:?} objective={:016x} iterations={} basis={:?} x={:?}",
        r.status,
        r.objective.to_bits(),
        sparse.iterations(),
        sparse.basis_columns(),
        r.x.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
    );
    #[cfg(feature = "phase1_progress_diagnostics")]
    milp::simplex::phase1_progress_diagnostics::take();
    Ok(())
}
