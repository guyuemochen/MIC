//! 父域保存方式的固定工作量对照；不改变生产CLI与求解配置。
//! 用法：domain_stats <实例> <dfs|hybrid|hybrid_rc> <节点预算>
use milp::branch_bound::{BbOptions, NodeRule, SearchControl};
use milp::instance::load_instance;
use milp::solver::{solve_controlled, SolveOptions};
use std::{path::Path, time::Instant};
fn number(x: Option<f64>) -> String {
    x.filter(|v| v.is_finite()).map_or("null".into(), |v| format!("{v:?}"))
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 { return Err("用法：domain_stats <实例> <dfs|hybrid|hybrid_rc> <节点预算>".into()); }
    let mut bb = BbOptions { node_limit: args[2].parse()?, ..BbOptions::default() };
    match args[1].as_str() {
        "dfs" => {},
        "hybrid" | "hybrid_rc" => {
            bb.node_rule = NodeRule::HybridBestBound;
            bb.local_cuts = false;
            bb.root_rc_fixing = args[1] == "hybrid_rc";
        },
        _ => return Err("未知模式".into()),
    }
    let p = load_instance(Path::new(&args[0]))?;
    let start = Instant::now();
    // 不设墙钟时限、不主动降级，便于逐字段验证完全相同的搜索工作量。
    let r = solve_controlled(&p, &SolveOptions { bb, ..SolveOptions::default() }, &SearchControl::default(), None)?;
    let elapsed = start.elapsed().as_secs_f64();
    let values = match r.report.values.as_ref() {
        Some(xs) if xs.iter().all(|x| x.is_finite()) => format!("{xs:?}"),
        Some(_) => return Err("求解返回非有限解".into()),
        None => "null".into(),
    };
    let s = &r.report.stats;
    println!("{{\"status\":\"{:?}\",\"objective\":{},\"bound\":{},\"gap_percent\":{},\"nodes\":{},\"lp_iterations\":{},\"cuts\":{},\"rc_fixed\":{},\"rc_reapplies\":{},\"rc_refixed\":{},\"auto_degraded\":{},\"values\":{},\"wall_secs\":{}}}", r.report.status, number(r.report.objective), number(r.bound), number(r.gap_percent), s.nodes, s.lp_iterations, s.cuts, s.rc_fixed, s.rc_reapplies, s.rc_refixed, r.auto_degraded, values, elapsed);
    Ok(())
}
