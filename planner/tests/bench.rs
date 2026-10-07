//! 性能基准（工作包 wp3）：全部 `#[ignore]` + `#[test]`，供手动运行，不断言数值。
//!
//! ```text
//! cargo test --release --manifest-path planner/Cargo.toml --test bench -- --ignored --nocapture
//! ```
//!
//! 对 6 个算例分别计时 `construct(FF/BF)` 与 `local_search(TabuSearch{tabu_size:8},
//! Steps(300_000))`，打印表格（算例、构造分、搜索分、耗时）。数值仅性能/质量参考，
//! 不构成验收断言；务必 `--release`。

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use planner::construct::{construct, ConstructionHeuristic};
use planner::format::parse_file;
use planner::localsearch::{local_search, LocalSearchConfig, Metaheuristic, Termination};
use planner::moves::MoveKind;
use planner::score::ScoreCalculator;

const INSTANCES: &[&str] = &[
    "roster-small",
    "cloud-small",
    "timetable-small",
    "assign-10",
    "roster-medium",
    "cloud-medium",
];

/// 算例绝对路径（经 CARGO_MANIFEST_DIR；测试进程 cwd 不定）。
fn instance_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("instances")
        .join(format!("{name}.txt"))
}

/// 基准表：每算例计时 construct(FF/BF) 与 TabuSearch{8} 300k 步搜索（seed 42、
/// Mixed{0.5}、samples 32、FF 起点——与验收测试同配置，仅步数加大）。
#[test]
#[ignore]
fn bench_construct_and_tabu_300k() {
    println!(
        "=== planner bench（release；FF/BF 构造 + TabuSearch{{tabu_size:8}} Steps(300_000) seed 42 Mixed{{0.5}} samples 32）==="
    );
    println!(
        "{:<16} {:>12} {:>12} {:>12} {:>10} {:>10} {:>10}",
        "instance", "FF (h,s)", "BF (h,s)", "best (h,s)", "FF ms", "BF ms", "search ms"
    );
    for &name in INSTANCES {
        let path = instance_path(name);
        let model = Arc::new(
            parse_file(&path).unwrap_or_else(|e| panic!("{name}: parse failed: {e}")),
        );
        let calc = ScoreCalculator::new(model.clone());

        let t0 = Instant::now();
        let ff = construct(&model, ConstructionHeuristic::FirstFit)
            .unwrap_or_else(|e| panic!("{name}: FF failed: {e}"));
        let ff_ms = t0.elapsed().as_millis();
        let ff_score = calc.total_score(&ff);

        let t1 = Instant::now();
        let bf = construct(&model, ConstructionHeuristic::BestFit)
            .unwrap_or_else(|e| panic!("{name}: BF failed: {e}"));
        let bf_ms = t1.elapsed().as_millis();
        let bf_score = calc.total_score(&bf);

        let cfg = LocalSearchConfig {
            metaheuristic: Metaheuristic::TabuSearch { tabu_size: 8 },
            seed: 42,
            termination: Termination::Steps(300_000),
            move_kind: MoveKind::Mixed { change_ratio: 0.5 },
            samples_per_step: 32,
        };
        let t2 = Instant::now();
        let r = local_search(ff, &cfg).unwrap_or_else(|e| panic!("{name}: search failed: {e}"));
        let search_ms = t2.elapsed().as_millis();

        println!(
            "{:<16} {:>12} {:>12} {:>12} {:>10} {:>10} {:>10}",
            name,
            format!("({},{})", ff_score.hard, ff_score.soft),
            format!("({},{})", bf_score.hard, bf_score.soft),
            format!("({},{})", r.best_score.hard, r.best_score.soft),
            ff_ms,
            bf_ms,
            search_ms
        );
    }
}
