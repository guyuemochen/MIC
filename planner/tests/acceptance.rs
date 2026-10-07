//! 统一验收测试（工作包 wp3）。
//!
//! ## 运行方式
//! **务必 `--release`**（`must_match_known_optima` 的 200k 步与 6 算例 × 多配置搜索在
//! debug 下过慢）：
//! ```text
//! cargo test --manifest-path planner/Cargo.toml --release
//! ```
//!
//! ## 约定
//! - 算例路径一律经 `env!("CARGO_MANIFEST_DIR")` 定位（测试进程 cwd 不定，禁止相对路径）。
//! - 本文件调用 moves / localsearch / solver 的**真实实现**，须在 wp2 落地后运行
//!   （此前运行会在 `todo!()` 处失败——这是任务书声明的先后依赖，不是本包缺陷）。
//! - 精确配置（TabuSearch{tabu_size:8}、LateAcceptance{history_len:30}）取自
//!   `localsearch` 契约的唯一定义，不另造数值。

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Arc;

use planner::construct::{construct, ConstructionHeuristic};
use planner::format::parse_file;
use planner::localsearch::{local_search, LocalSearchConfig, Metaheuristic, Termination};
use planner::model::HardSoftScore;
use planner::moves::{random_move, MoveKind};
use planner::parallel::{solve_parallel, solve_parallel_detailed};
use planner::rng::Rng;
use planner::score::{IncrementalScorer, ScoreCalculator};
use planner::solver::{solve, SolverConfig};

/// 单个算例的验收数据。
struct Case {
    name: &'static str,
    /// FF 构造分基线（hard 恒断言；medium 的 soft 不锁定 → None）。
    ff_hard: i64,
    ff_soft: Option<i64>,
    /// BF 构造分基线（medium 不锁定 → None）。
    bf: Option<(i64, i64)>,
    /// 必达最优（None = 仅参考值，无断言依赖）。
    optimum: Option<i64>,
}

/// 基线表（出处：rustc 1.95 原型剪枝 DFS + 评审方 python verify.py 两轮独立复核相等）。
const CASES: &[Case] = &[
    Case {
        name: "roster-small",
        ff_hard: 0,
        ff_soft: Some(36),
        bf: Some((0, 27)),
        optimum: Some(27),
    },
    Case {
        name: "cloud-small",
        ff_hard: 6, // FF 起点带硬违反，由局部搜索修复至 0
        ff_soft: Some(16),
        bf: Some((3, 13)),
        optimum: Some(13),
    },
    Case {
        name: "timetable-small",
        ff_hard: 0,
        ff_soft: Some(16),
        bf: Some((0, 17)),
        optimum: Some(14),
    },
    Case {
        name: "assign-10",
        ff_hard: 0,
        ff_soft: Some(113),
        bf: Some((0, 41)),
        optimum: Some(36),
    },
    // 中等算例：仅锁定 FF hard == 0；BF 只测确定性与完整性；122 / 77 为原型参考值
    // （无断言依赖，见 README「算例清单与已知值」）。
    Case {
        name: "roster-medium",
        ff_hard: 0,
        ff_soft: None,
        bf: None,
        optimum: None,
    },
    Case {
        name: "cloud-medium",
        ff_hard: 0,
        ff_soft: None,
        bf: None,
        optimum: None,
    },
];

/// 算例绝对路径（经 CARGO_MANIFEST_DIR；测试进程 cwd 不定）。
fn instance_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("instances")
        .join(format!("{name}.txt"))
}

fn load(name: &str) -> Arc<planner::model::ProblemModel> {
    let path = instance_path(name);
    Arc::new(parse_file(&path).unwrap_or_else(|e| panic!("{name}: parse failed: {e}")))
}

/// 1) 增量一致性：每算例随机游走 5000 步（seed 42、Mixed{0.5}），每步断言
///    preview == 应用后 inc.score() == ScoreCalculator::total_score（全量真值）。
#[test]
fn incremental_consistency() {
    for case in CASES {
        let model = load(case.name);
        let calc = ScoreCalculator::new(model.clone());
        let mut sol = construct(&model, ConstructionHeuristic::FirstFit)
            .unwrap_or_else(|e| panic!("{}: construct failed: {e}", case.name));
        let mut inc = IncrementalScorer::new(model.clone());
        inc.reset(&sol);
        // 入口不变式：影子分 == 全量重算
        assert_eq!(
            inc.score(),
            calc.total_score(&sol),
            "{}: reset score != full recompute",
            case.name
        );
        let mut rng = Rng::new(42);
        let mut applied = 0usize;
        for step in 0..5000 {
            let Some(mv) = random_move(&sol, &inc, &mut rng, MoveKind::Mixed { change_ratio: 0.5 })
            else {
                continue; // 全域退化才返回 None（6 算例不可达）
            };
            let preview = mv.preview(&inc);
            mv.apply(&mut sol, &mut inc);
            assert_eq!(
                preview,
                inc.score(),
                "{}: step {step}: preview != shadow score after apply",
                case.name
            );
            assert_eq!(
                preview,
                calc.total_score(&sol),
                "{}: step {step}: preview != full recompute",
                case.name
            );
            applied += 1;
        }
        assert!(applied > 0, "{}: no move applied in 5000 steps", case.name);
    }
}

/// 2) 构造确定性：FF/BF 各跑两次赋值逐位相等、is_complete、构造分等于基线表
///    （此处从算例文件读取，与 wp1 在 format.rs 内嵌字符串的测试互为独立复核）。
#[test]
fn construction_deterministic() {
    for case in CASES {
        let model = load(case.name);
        let calc = ScoreCalculator::new(model.clone());

        // ---- First Fit ----
        let a = construct(&model, ConstructionHeuristic::FirstFit)
            .unwrap_or_else(|e| panic!("{}: FF failed: {e}", case.name));
        let b = construct(&model, ConstructionHeuristic::FirstFit)
            .unwrap_or_else(|e| panic!("{}: FF #2 failed: {e}", case.name));
        assert!(a.is_complete(), "{}: FF solution incomplete", case.name);
        assert!(
            a.assignment() == b.assignment(),
            "{}: FF not deterministic",
            case.name
        );
        let s = calc.total_score(&a);
        assert_eq!(
            s.hard, case.ff_hard,
            "{}: FF hard {} != baseline {}",
            case.name, s.hard, case.ff_hard
        );
        if let Some(want) = case.ff_soft {
            assert_eq!(
                s.soft, want,
                "{}: FF soft {} != baseline {}",
                case.name, s.soft, want
            );
        }

        // ---- Best Fit ----
        let a = construct(&model, ConstructionHeuristic::BestFit)
            .unwrap_or_else(|e| panic!("{}: BF failed: {e}", case.name));
        let b = construct(&model, ConstructionHeuristic::BestFit)
            .unwrap_or_else(|e| panic!("{}: BF #2 failed: {e}", case.name));
        assert!(a.is_complete(), "{}: BF solution incomplete", case.name);
        assert!(
            a.assignment() == b.assignment(),
            "{}: BF not deterministic",
            case.name
        );
        if let Some((want_hard, want_soft)) = case.bf {
            let s = calc.total_score(&a);
            assert_eq!(
                (s.hard, s.soft),
                (want_hard, want_soft),
                "{}: BF score {:?} != baseline ({}, {})",
                case.name,
                (s.hard, s.soft),
                want_hard,
                want_soft
            );
        }
    }
}

/// 搜索配置（契约默认：seed 42、Steps(50_000)、Mixed{0.5}、samples 32、FF 起点）。
fn ls_config(metaheuristic: Metaheuristic, steps: usize, seed: u64) -> LocalSearchConfig {
    LocalSearchConfig {
        metaheuristic,
        seed,
        termination: Termination::Steps(steps),
        move_kind: MoveKind::Mixed { change_ratio: 0.5 },
        samples_per_step: 32,
    }
}

/// 3) 搜索不劣于构造：全部 6 算例 × 两种元启发式（TabuSearch{8} / LateAcceptance{30}，
///    seed=42、Steps(50_000)、Mixed{0.5}、FF 起点）断言 best_score ≤ initial_score（字典序）。
#[test]
fn search_no_worse_than_construction() {
    let metas = [
        Metaheuristic::TabuSearch { tabu_size: 8 },
        Metaheuristic::LateAcceptance { history_len: 30 },
    ];
    for case in CASES {
        let model = load(case.name);
        let init = construct(&model, ConstructionHeuristic::FirstFit)
            .unwrap_or_else(|e| panic!("{}: FF failed: {e}", case.name));
        for meta in &metas {
            let name = match meta {
                Metaheuristic::TabuSearch { .. } => "tabu",
                Metaheuristic::LateAcceptance { .. } => "la",
                Metaheuristic::SimulatedAnnealing { .. } => "sa",
            };
            let cfg = ls_config(meta.clone(), 50_000, 42);
            let r = local_search(init.clone(), &cfg)
                .unwrap_or_else(|e| panic!("{}: {name} search failed: {e}", case.name));
            assert!(
                r.best_score <= r.initial_score,
                "{}: {name} best {:?} > initial {:?}",
                case.name,
                r.best_score,
                r.initial_score
            );
            assert!(
                r.best.is_complete(),
                "{}: {name} best solution incomplete",
                case.name
            );
        }
    }
}

/// 4) 可行且可复现：全部 6 算例以默认配置（SolverConfig::default()）求解两次，
///    断言 best_score.hard == 0 且两次 best.assignment() 逐位相等、best_score 相等。
#[test]
fn feasible_and_reproducible() {
    for case in CASES {
        let model = load(case.name);
        let r1 = solve(model.clone(), &SolverConfig::default())
            .unwrap_or_else(|e| panic!("{}: solve #1 failed: {e}", case.name));
        let r2 = solve(model, &SolverConfig::default())
            .unwrap_or_else(|e| panic!("{}: solve #2 failed: {e}", case.name));
        assert_eq!(
            r1.best_score.hard, 0,
            "{}: default config must reach feasible (got {:?})",
            case.name, r1.best_score
        );
        assert_eq!(
            r1.best_score, r2.best_score,
            "{}: best_score differs between two identical runs",
            case.name
        );
        assert!(
            r1.best.assignment() == r2.best.assignment(),
            "{}: assignment differs between two identical runs",
            case.name
        );
        assert!(r1.best.is_complete(), "{}: best incomplete", case.name);
    }
}

/// 5) 必达最优：4 个必达算例按种子表 [42, 7, 2026] 依次以 Steps(200_000)、
///    FF + TabuSearch{tabu_size:8} 求解，任一种子命中 (hard == 0 && soft == 期望值) 即通过
///    （两套独立实现均实测 seed 42/7 在 50k 步已命中，200k 留裕量；跨实现不承诺轨迹一致，
///    故不限定单一种子）。
#[test]
fn must_match_known_optima() {
    let seeds = [42u64, 7, 2026];
    for case in CASES.iter().filter(|c| c.optimum.is_some()) {
        let want = case.optimum.unwrap();
        let model = load(case.name);
        let mut hit: Option<(u64, i64)> = None;
        for &seed in &seeds {
            let cfg = SolverConfig {
                construction: ConstructionHeuristic::FirstFit,
                local_search: ls_config(Metaheuristic::TabuSearch { tabu_size: 8 }, 200_000, seed),
                run_local_search: true,
            };
            let r = solve(model.clone(), &cfg)
                .unwrap_or_else(|e| panic!("{}: solve(seed {seed}) failed: {e}", case.name));
            if r.best_score.hard == 0 && r.best_score.soft == want {
                hit = Some((seed, r.best_score.soft));
                break;
            }
        }
        assert!(
            hit.is_some(),
            "{}: no seed in [42, 7, 2026] hit (hard=0, soft={want}) within 200k steps",
            case.name
        );
    }
}

/// 6) CLI 冒烟：默认参数运行 6 算例，断言 exit 0、stdout 恰两行、首行 "hard = 0"、
///    末行形如 "objective = <非负整数>"；必达算例另断言 objective == 期望最优；
///    非法输入（--time-secs NaN / -1）断言 exit 1 且无 panic 输出。
#[test]
fn cli_smoke() {
    let exe = env!("CARGO_BIN_EXE_planner");
    for case in CASES {
        let out = Command::new(exe)
            .arg(instance_path(case.name))
            .output()
            .unwrap_or_else(|e| panic!("{}: spawn planner failed: {e}", case.name));
        assert!(
            out.status.success(),
            "{}: planner exited with {:?}",
            case.name,
            out.status.code()
        );
        let stdout = String::from_utf8_lossy(&out.stdout);
        let lines: Vec<&str> = stdout.lines().collect();
        assert_eq!(
            lines.len(),
            2,
            "{}: stdout must be exactly 2 lines, got {stdout:?}",
            case.name
        );
        assert_eq!(lines[0], "hard = 0", "{}: first line {:?}", case.name, lines[0]);
        let obj_text = lines[1]
            .strip_prefix("objective = ")
            .unwrap_or_else(|| panic!("{}: last line not 'objective = <int>': {:?}", case.name, lines[1]));
        let objective: i64 = obj_text
            .parse()
            .unwrap_or_else(|_| panic!("{}: objective not an integer: {obj_text:?}", case.name));
        assert!(objective >= 0, "{}: negative objective {objective}", case.name);
        if let Some(want) = case.optimum {
            assert_eq!(
                objective, want,
                "{}: CLI objective {objective} != known optimum {want}",
                case.name
            );
        }
    }

    // 非法 --time-secs：exit 1，stderr "error: ..."，无 panic
    for bad in ["NaN", "-1"] {
        let out = Command::new(exe)
            .arg(instance_path("roster-small"))
            .arg("--time-secs")
            .arg(bad)
            .output()
            .unwrap_or_else(|e| panic!("spawn planner failed: {e}"));
        assert_eq!(
            out.status.code(),
            Some(1),
            "--time-secs {bad}: exit code {:?}",
            out.status.code()
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        let stdout = String::from_utf8_lossy(&out.stdout);
        assert!(
            stderr.starts_with("error:"),
            "--time-secs {bad}: stderr {stderr:?}"
        );
        assert!(
            !stderr.contains("panicked") && !stdout.contains("panicked"),
            "--time-secs {bad}: panic output (stderr {stderr:?}, stdout {stdout:?})"
        );
    }
}

/// 7) 并行必达最优：4 个必达算例 × `solve_parallel(…, 4)`（FF + TabuSearch{8} +
///    Steps(200_000) + seed 42 + Mixed{0.5} + samples 32）。worker 0 = 主 seed 42 =
///    `must_match_known_optima` 的首种子路径（确定性命中），聚合取 min 而 `(0, 最优)`
///    已是下界 ⇒ 断言 `best_score == (0, 期望最优)`。
#[test]
fn parallel_matches_known_optima() {
    for case in CASES.iter().filter(|c| c.optimum.is_some()) {
        let want = case.optimum.unwrap();
        let model = load(case.name);
        let cfg = SolverConfig {
            construction: ConstructionHeuristic::FirstFit,
            local_search: ls_config(Metaheuristic::TabuSearch { tabu_size: 8 }, 200_000, 42),
            run_local_search: true,
        };
        let r = solve_parallel(model, &cfg, 4)
            .unwrap_or_else(|e| panic!("{}: solve_parallel failed: {e}", case.name));
        assert_eq!(
            r.best_score,
            HardSoftScore::new(0, want),
            "{}: parallel best {:?} != known optimum (0, {want})",
            case.name,
            r.best_score
        );
    }
}

/// 8) 并行不劣于单线程（结构性质）：roster-medium / cloud-medium 以同 cfg
///    （Steps(50_000)、seed 42、tabu）对比 `solve` 与 `solve_parallel_detailed(…, 4)`：
///    (a) 聚合 ≤ 单线程（字典序）；(b) 恰 4 个 worker；(c) worker 0 分数 == `solve` 的
///    best_score（端到端钉死 worker-0 = 主 seed 规则）；(d) 两者 hard == 0
///    （与 `feasible_and_reproducible` 同基线）。
#[test]
fn parallel_no_worse_than_single_thread() {
    for name in ["roster-medium", "cloud-medium"] {
        let model = load(name);
        let cfg = SolverConfig {
            construction: ConstructionHeuristic::FirstFit,
            local_search: ls_config(Metaheuristic::TabuSearch { tabu_size: 8 }, 50_000, 42),
            run_local_search: true,
        };
        let t1 = solve(model.clone(), &cfg)
            .unwrap_or_else(|e| panic!("{name}: solve failed: {e}"));
        let t4 = solve_parallel_detailed(model, &cfg, 4)
            .unwrap_or_else(|e| panic!("{name}: solve_parallel_detailed failed: {e}"));
        assert!(
            t4.best.best_score <= t1.best_score,
            "{name}: parallel {:?} > single {:?}",
            t4.best.best_score,
            t1.best_score
        );
        assert_eq!(t4.worker_scores.len(), 4, "{name}: worker count");
        assert_eq!(
            t4.worker_scores[0], t1.best_score,
            "{name}: worker 0 must run with the main seed (worker-0 rule)"
        );
        assert_eq!(t1.best_score.hard, 0, "{name}: single-thread must be feasible");
        assert_eq!(t4.best.best_score.hard, 0, "{name}: parallel must be feasible");
    }
}

/// 9) 固定 seed + 固定 threads + Steps 终止 ⇒ 逐位可复现：roster-small / cloud-medium
///    （threads=3、seed 42、Steps(50_000)、tabu）两次运行 best_score 相等、
///    best.assignment() 逐位相等、winner 相等。
#[test]
fn parallel_reproducible_fixed_seed_and_threads() {
    for name in ["roster-small", "cloud-medium"] {
        let model = load(name);
        let cfg = SolverConfig {
            construction: ConstructionHeuristic::FirstFit,
            local_search: ls_config(Metaheuristic::TabuSearch { tabu_size: 8 }, 50_000, 42),
            run_local_search: true,
        };
        let a = solve_parallel_detailed(model.clone(), &cfg, 3)
            .unwrap_or_else(|e| panic!("{name}: run #1 failed: {e}"));
        let b = solve_parallel_detailed(model, &cfg, 3)
            .unwrap_or_else(|e| panic!("{name}: run #2 failed: {e}"));
        assert_eq!(
            a.best.best_score, b.best.best_score,
            "{name}: best_score differs between two identical runs"
        );
        assert!(
            a.best.best.assignment() == b.best.best.assignment(),
            "{name}: assignment differs between two identical runs"
        );
        assert_eq!(a.winner, b.winner, "{name}: winner differs");
    }
}

/// 10) 重复运行冒烟：roster-small（tabu）与 cloud-medium（la）以 threads=4 各连跑
///     5 轮——每轮 Ok、无 panic（join 全部由库保证，此处端到端冒烟）；Steps 终止下
///     各轮 best_score 逐位一致。
#[test]
fn parallel_repeated_runs_no_crash() {
    let setups = [
        ("roster-small", Metaheuristic::TabuSearch { tabu_size: 8 }),
        ("cloud-medium", Metaheuristic::LateAcceptance { history_len: 30 }),
    ];
    for (name, meta) in setups {
        let model = load(name);
        let cfg = SolverConfig {
            construction: ConstructionHeuristic::FirstFit,
            local_search: ls_config(meta, 50_000, 42),
            run_local_search: true,
        };
        let first = solve_parallel(model.clone(), &cfg, 4)
            .unwrap_or_else(|e| panic!("{name}: round 1 failed: {e}"));
        for round in 2..=5 {
            let r = solve_parallel(model.clone(), &cfg, 4)
                .unwrap_or_else(|e| panic!("{name}: round {round} failed: {e}"));
            assert_eq!(
                r.best_score, first.best_score,
                "{name}: round {round} best_score differs from round 1"
            );
        }
    }
}

/// 11) CLI `--threads` 冒烟：`--threads 4` 跑 roster-small → exit 0、stdout 恰两行、
///     首行 `hard = 0`、objective == 27，且与 `--threads 1` 的 objective 相等（必达算例
///     双双命中 27）；`--threads 0` / `--threads abc` / `--threads 4 --threads 2` →
///     exit 1 且 stderr 以 `error:` 开头、无 panicked 字样。
#[test]
fn cli_threads_smoke() {
    let exe = env!("CARGO_BIN_EXE_planner");

    // ---- 合法：--threads 4 与 --threads 1（后者即旧版 solve() 路径）----
    let mut objectives: Vec<i64> = Vec::new();
    for t in ["4", "1"] {
        let out = Command::new(exe)
            .arg(instance_path("roster-small"))
            .arg("--threads")
            .arg(t)
            .output()
            .unwrap_or_else(|e| panic!("--threads {t}: spawn planner failed: {e}"));
        assert!(
            out.status.success(),
            "--threads {t}: exited with {:?}",
            out.status.code()
        );
        let stdout = String::from_utf8_lossy(&out.stdout);
        let lines: Vec<&str> = stdout.lines().collect();
        assert_eq!(
            lines.len(),
            2,
            "--threads {t}: stdout must be exactly 2 lines, got {stdout:?}"
        );
        assert_eq!(lines[0], "hard = 0", "--threads {t}: first line {:?}", lines[0]);
        let obj_text = lines[1].strip_prefix("objective = ").unwrap_or_else(|| {
            panic!("--threads {t}: last line not 'objective = <int>': {:?}", lines[1])
        });
        let objective: i64 = obj_text
            .parse()
            .unwrap_or_else(|_| panic!("--threads {t}: objective not an integer: {obj_text:?}"));
        assert_eq!(objective, 27, "--threads {t}: roster-small objective != 27");
        objectives.push(objective);
    }
    assert_eq!(
        objectives[0], objectives[1],
        "--threads 4 and --threads 1 objectives differ"
    );

    // ---- 非法：越界 / 非整数 / 重复（进程断言模式同 cli_smoke）----
    let bad_cases: [&[&str]; 3] = [
        &["--threads", "0"],
        &["--threads", "abc"],
        &["--threads", "4", "--threads", "2"],
    ];
    for bad in bad_cases {
        let mut cmd = Command::new(exe);
        cmd.arg(instance_path("roster-small"));
        for a in bad {
            cmd.arg(a);
        }
        let out = cmd.output().unwrap_or_else(|e| panic!("{bad:?}: spawn planner failed: {e}"));
        assert_eq!(
            out.status.code(),
            Some(1),
            "{bad:?}: exit code {:?}",
            out.status.code()
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        let stdout = String::from_utf8_lossy(&out.stdout);
        assert!(stderr.starts_with("error:"), "{bad:?}: stderr {stderr:?}");
        assert!(
            !stderr.contains("panicked") && !stdout.contains("panicked"),
            "{bad:?}: panic output (stderr {stderr:?}, stdout {stdout:?})"
        );
    }
}
