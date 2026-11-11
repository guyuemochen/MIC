//! 稀疏 LP 成本探针：只解 LP 松弛，不做预处理或分支定界。
//! 用法：cargo run --release --features sparse_profile --example sparse_cost -- 算例 [迭代上限]
//! 正式端到端计时关闭 sparse_profile；扫描候选在独立实验分支验证。
use milp::instance::load_instance;
use milp::model::Tolerances;
use milp::simplex::Bounds;
use milp::simplex_sparse::SparseSimplex;
use std::path::Path;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let path = args.next().ok_or("缺少算例路径")?;
    let limit: Option<u64> = args.next().map(|v| v.parse()).transpose()?;
    if args.next().is_some() {
        return Err("参数过多".into());
    }
    let problem = load_instance(Path::new(&path))?;
    #[cfg(feature = "sparse_profile")]
    milp::simplex_sparse::cost_profile::reset();
    let start = Instant::now();
    let mut lp = SparseSimplex::new(&problem, Tolerances::default())?;
    lp.set_iter_limit_hint(limit);
    let result = lp.solve(&Bounds::from_problem(&problem));
    println!(
        "status={:?} objective={:?} iterations={} wall_secs={:.6}",
        result.status,
        result.objective,
        lp.iterations(),
        start.elapsed().as_secs_f64()
    );
    #[cfg(feature = "sparse_profile")]
    {
        let costs = milp::simplex_sparse::cost_profile::snapshot();
        for (k, label) in ["refactor", "pricing", "basis_solve", "basis_update"]
            .iter()
            .enumerate()
        {
            println!(
                "{label}_calls={} {label}_secs={:.6}",
                costs.calls[k],
                costs.nanos[k] as f64 / 1e9
            );
        }
    }
    Ok(())
}
