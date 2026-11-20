//! 缩放与完整预处理后，调用生产二元传播认证，输出可认证行数；不运行 LP。
#[allow(dead_code)]
#[path = "../../milp/src/binary_propagation.rs"]
mod binary_propagation;

use milp::{model, simplex};
use std::{error::Error, path::Path};

fn main() -> Result<(), Box<dyn Error>> {
    let tol = model::Tolerances::default();
    for path in std::env::args_os().skip(1) {
        let path = Path::new(&path);
        let model = milp::instance::load_instance(path)?;
        model.validate()?;
        let (scaled, _) = milp::scaling::scale_problem(&model);
        let (reduced, _) = milp::presolve::presolve(&scaled, &tol)?;
        let stats = binary_propagation::BinaryPropagation::new(&reduced, tol).stats;
        println!(
            "{}\t{}\t{}",
            path.display(),
            stats.accepted_rows,
            stats.skipped_rows
        );
    }
    Ok(())
}
