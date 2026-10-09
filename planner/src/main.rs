//! planner CLI（工作包 wp3「CLI、算例、基准与验收」）。
//!
//! ## 用法
//! ```text
//! planner <算例路径> [--steps N] [--time-secs F] [--seed S]
//!          [--algo tabu|sa|la] [--construct ff|bf] [--samples K]
//!          [--threads N] [--verbose]
//! ```
//! 默认 steps=50000、seed=42、algo=tabu、construct=ff、samples=32、threads=1。
//!
//! ## 参数语义
//! - `--steps` 与 `--time-secs` 互斥（同给 exit 1）。
//! - `--threads N`：多起点并行 worker 数（`1 ≤ N ≤ 1024`，默认 1）。threads=1 走与
//!   旧版**完全相同**的 `solve()` 调用路径（行为与输出逐字节一致）；N > 1 时
//!   `solve_parallel` 由父线程 validate+construct **恰一次**得同一初始解，N 个 worker
//!   各自对该初始解独立 local_search（worker 0 用主 seed、i ≥ 1 用 `rng::derive_seed`
//!   派生），结束取最优。与 `--steps`/`--time-secs`
//!   无互斥，stdout 恰两行契约与 exit code 契约不变。
//! - `--algo tabu` → `TabuSearch { tabu_size: 8 }`；
//!   `--algo sa` 固定 `SimulatedAnnealing { start_temp: 1000.0, end_temp: 1.0 }`；
//!   `--algo la` 固定 `LateAcceptance { history_len: 30 }`（参数取自 `localsearch`
//!   契约的唯一默认定义，不新造数值；无额外调参旗标）。
//! - `--samples` 仅 tabu 读取（每步采样数）；sa/la 下接受但无效果（不耗随机数）。
//! - 移动配比固定 `Mixed { change_ratio: 0.5 }`（无 CLI 旗标）。
//!
//! ## 数值安全（normative，全程不 panic）
//! - `--time-secs`：解析为 f64 后先 `v.is_finite() && 0.0 <= v <= 36000.0`，
//!   再 `Duration::try_from_secs_f64(v)`，任一不过 → stderr `error: invalid --time-secs`
//!   + exit 1。**禁止 `Duration::from_secs_f64`**（NaN/负/非有限/溢出直接 panic，std 文档明示）。
//! - `--steps`/`--seed`/`--samples`：解析为非负整数（u64），非法或超界
//!   （steps > 10^9、samples ∉ [1, 10^6]）exit 1。
//!
//! ## 输出契约
//! 成功时 stdout **恰好两行**且在结束时输出（均取 `SearchResult.best_score`，`hard` 后无
//! 多余空白、`objective` 后无内容）：
//! ```text
//! hard = <i64>
//! objective = <i64>
//! ```
//! 错误：stderr `error: <消息>` + exit 1；正常求解（无论 hard 是否为 0）exit 0。
//!
//! ## --verbose
//! 进度只走 stderr：构造分（求解前单独构造一次打印，与 `solve` 内部构造确定性一致）、
//! threads 数、搜索结束后的初始分 / 实际步数 / 历史最优；threads > 1 时另追加胜出
//! worker 与各 worker 的 `best = (h, s)` 一行。**步级进度（每 5000 步当前最优）**：lib 公共
//! API（`solver::solve` / `localsearch::local_search`）不提供逐步回调，故 CLI 无法在搜索
//! 过程中报告中间最优——这是当前 API 的固有限制，见 README「可复现性与 verbose」一节。

#![forbid(unsafe_code)]

use std::path::Path;
use std::sync::Arc;
use std::time::Duration;

use planner::construct::{construct, ConstructionHeuristic};
use planner::format::parse_file;
use planner::localsearch::{LocalSearchConfig, Metaheuristic, Termination};
use planner::model::HardSoftScore;
use planner::moves::MoveKind;
use planner::parallel::solve_parallel_detailed;
use planner::score::ScoreCalculator;
use planner::solver::{solve, SolverConfig};

/// 元启发式选择（CLI 旗标 → 固定参数，映射全文唯一）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Algo {
    /// `--algo tabu` → `TabuSearch { tabu_size: 8 }`。
    Tabu,
    /// `--algo sa` → `SimulatedAnnealing { start_temp: 1000.0, end_temp: 1.0 }`。
    Sa,
    /// `--algo la` → `LateAcceptance { history_len: 30 }`。
    La,
}

/// 解析后的 CLI 参数（None = 未给，用默认值）。
#[derive(Debug, Clone, PartialEq)]
struct Args {
    path: String,
    steps: Option<u64>,
    time_limit: Option<Duration>,
    seed: u64,
    algo: Algo,
    construction: ConstructionHeuristic,
    samples: u64,
    threads: u64,
    verbose: bool,
}

/// 硬界：steps ≤ 10^9、samples ∈ [1, 10^6]（越界按非法参数 exit 1）。
const MAX_STEPS: u64 = 1_000_000_000;
const MAX_SAMPLES: u64 = 1_000_000;
/// threads 合法上界（与 lib 侧 `parallel::MAX_THREADS` 数值一致，单测自查；
/// 合法域 [1, 1024]）。
const MAX_THREADS: u64 = 1024;
/// 默认值（与 `LocalSearchConfig::default` 的契约默认同源）。
const DEFAULT_STEPS: u64 = 50_000;
const DEFAULT_SEED: u64 = 42;
const DEFAULT_SAMPLES: u64 = 32;
/// 默认单线程（threads=1 走与旧版完全相同的 `solve()` 路径）。
const DEFAULT_THREADS: u64 = 1;

/// 读取旗标的值参数（消耗 `argv[i+1]`）；缺失返回 `Err`（exit 1 路径，不 panic）。
fn flag_value<'a>(argv: &'a [String], i: &mut usize, flag: &str) -> Result<&'a str, String> {
    *i += 1;
    argv.get(*i)
        .map(|s| s.as_str())
        .ok_or_else(|| format!("{flag} expects a value"))
}

/// 手写参数解析（std only）。任何用户输入只产生 `Err(String)`，不 panic。
fn parse_args(argv: &[String]) -> Result<Args, String> {
    const USAGE: &str = "usage: planner <算例路径> [--steps N] [--time-secs F] [--seed S] \
                         [--algo tabu|sa|la] [--construct ff|bf] [--samples K] \
                         [--threads N] [--verbose]";

    let mut path: Option<String> = None;
    let mut steps: Option<u64> = None;
    let mut time_limit: Option<Duration> = None;
    let mut seed: Option<u64> = None;
    let mut algo: Option<Algo> = None;
    let mut construction: Option<ConstructionHeuristic> = None;
    let mut samples: Option<u64> = None;
    let mut threads: Option<u64> = None;
    let mut verbose = false;

    let mut i = 0usize;
    while i < argv.len() {
        let arg = argv[i].as_str();
        match arg {
            "--verbose" => verbose = true,
            "--steps" => {
                if steps.is_some() {
                    return Err("duplicate --steps".to_string());
                }
                if time_limit.is_some() {
                    return Err("--steps and --time-secs are mutually exclusive".to_string());
                }
                let raw = flag_value(argv, &mut i, "--steps")?;
                let v: u64 = raw
                    .parse()
                    .map_err(|_| format!("invalid --steps '{raw}' (expected integer in [0, 10^9])"))?;
                if v > MAX_STEPS {
                    return Err(format!("invalid --steps '{raw}': exceeds 10^9"));
                }
                steps = Some(v);
            }
            "--time-secs" => {
                if time_limit.is_some() {
                    return Err("duplicate --time-secs".to_string());
                }
                if steps.is_some() {
                    return Err("--steps and --time-secs are mutually exclusive".to_string());
                }
                let raw = flag_value(argv, &mut i, "--time-secs")?;
                // 数值安全：先有限性与范围预检，再 try_from_secs_f64（禁 from_secs_f64，
                // 其对 NaN/负/溢出直接 panic）。两道检查任一不过均为同一错误消息。
                let v: f64 = raw
                    .parse()
                    .map_err(|_| "invalid --time-secs".to_string())?;
                if !(v.is_finite() && (0.0..=36000.0).contains(&v)) {
                    return Err("invalid --time-secs".to_string());
                }
                let d = Duration::try_from_secs_f64(v)
                    .map_err(|_| "invalid --time-secs".to_string())?;
                time_limit = Some(d);
            }
            "--seed" => {
                if seed.is_some() {
                    return Err("duplicate --seed".to_string());
                }
                let raw = flag_value(argv, &mut i, "--seed")?;
                let v: u64 = raw
                    .parse()
                    .map_err(|_| format!("invalid --seed '{raw}' (expected u64)"))?;
                seed = Some(v);
            }
            "--algo" => {
                if algo.is_some() {
                    return Err("duplicate --algo".to_string());
                }
                let raw = flag_value(argv, &mut i, "--algo")?;
                algo = match raw {
                    "tabu" => Some(Algo::Tabu),
                    "sa" => Some(Algo::Sa),
                    "la" => Some(Algo::La),
                    other => {
                        return Err(format!("invalid --algo '{other}' (expected tabu|sa|la)"))
                    }
                };
            }
            "--construct" => {
                if construction.is_some() {
                    return Err("duplicate --construct".to_string());
                }
                let raw = flag_value(argv, &mut i, "--construct")?;
                construction = match raw {
                    "ff" => Some(ConstructionHeuristic::FirstFit),
                    "bf" => Some(ConstructionHeuristic::BestFit),
                    other => {
                        return Err(format!("invalid --construct '{other}' (expected ff|bf)"))
                    }
                };
            }
            "--samples" => {
                if samples.is_some() {
                    return Err("duplicate --samples".to_string());
                }
                let raw = flag_value(argv, &mut i, "--samples")?;
                let v: u64 = raw
                    .parse()
                    .map_err(|_| format!("invalid --samples '{raw}' (expected integer in [1, 10^6])"))?;
                if v < 1 {
                    // 与 lib 侧运行期 `samples_per_step ≥ 1` 校验同域：解析期即拒绝 0，
                    // 声明域与强制域一致（此前 0 可通过解析、运行期才被拒）。
                    return Err(format!("invalid --samples '{raw}': must be >= 1"));
                }
                if v > MAX_SAMPLES {
                    return Err(format!("invalid --samples '{raw}': exceeds 10^6"));
                }
                samples = Some(v);
            }
            "--threads" => {
                if threads.is_some() {
                    return Err("duplicate --threads".to_string());
                }
                let raw = flag_value(argv, &mut i, "--threads")?;
                let v: u64 = raw.parse().map_err(|_| {
                    format!("invalid --threads '{raw}' (expected integer in [1, 1024])")
                })?;
                if v < 1 || v > MAX_THREADS {
                    return Err(format!("invalid --threads '{raw}': must be within [1, 1024]"));
                }
                threads = Some(v);
            }
            other if other.starts_with("--") => {
                return Err(format!("unknown option '{other}'"));
            }
            other => {
                if path.is_some() {
                    return Err(format!("unexpected extra argument '{other}'"));
                }
                path = Some(other.to_string());
            }
        }
        i += 1;
    }

    let path = path.ok_or_else(|| format!("missing instance path ({USAGE})"))?;

    Ok(Args {
        path,
        steps,
        time_limit,
        seed: seed.unwrap_or(DEFAULT_SEED),
        algo: algo.unwrap_or(Algo::Tabu),
        construction: construction.unwrap_or(ConstructionHeuristic::FirstFit),
        samples: samples.unwrap_or(DEFAULT_SAMPLES),
        threads: threads.unwrap_or(DEFAULT_THREADS),
        verbose,
    })
}

/// 主流程：parse_file → solve → 输出。所有错误统一 `Err(String)` → `error: <消息>` + exit 1。
fn run(args: &Args) -> Result<(), String> {
    let model = Arc::new(parse_file(Path::new(&args.path)).map_err(|e| e.to_string())?);

    if args.verbose {
        // 构造是确定性的（不耗随机数）：此处单独构造一次仅为打印构造分，
        // 与 solve() 内部的构造结果逐位一致。
        let sol = construct(&model, args.construction).map_err(|e| e.to_string())?;
        let s = ScoreCalculator::new(model.clone()).total_score(&sol);
        eprintln!("[verbose] instance   : {}", args.path);
        eprintln!("[verbose] construct  : hard = {}, objective = {}", s.hard, s.soft);
        eprintln!("[verbose] threads    : {}", args.threads);
    }

    let termination = match args.time_limit {
        Some(d) => Termination::TimeLimit(d),
        None => Termination::Steps(
            usize::try_from(args.steps.unwrap_or(DEFAULT_STEPS))
                .map_err(|_| "invalid --steps".to_string())?,
        ),
    };
    let metaheuristic = match args.algo {
        Algo::Tabu => Metaheuristic::TabuSearch { tabu_size: 8 },
        Algo::Sa => Metaheuristic::SimulatedAnnealing {
            start_temp: 1000.0,
            end_temp: 1.0,
        },
        Algo::La => Metaheuristic::LateAcceptance { history_len: 30 },
    };
    let samples = usize::try_from(args.samples).map_err(|_| "invalid --samples".to_string())?;
    let cfg = SolverConfig {
        construction: args.construction,
        local_search: LocalSearchConfig {
            metaheuristic,
            seed: args.seed,
            termination,
            move_kind: MoveKind::Mixed { change_ratio: 0.5 },
            samples_per_step: samples,
        },
        run_local_search: true,
    };
    // ---- 求解分派 ----
    // threads=1：与旧版完全相同的 solve() 调用路径（行为与输出逐字节一致，结构性保证）。
    // threads>1：多起点并行；solve_parallel 即 solve_parallel_detailed(...)?.best
    // （lib 中前者逐字委托后者），此处直接调 detailed 以便 --verbose 打印 winner /
    // worker_scores——非 verbose 下结果与 solve_parallel 逐位一致。
    let mut parallel_detail: Option<(usize, Vec<HardSoftScore>)> = None;
    let result = if args.threads == 1 {
        solve(model, &cfg).map_err(|e| e.to_string())?
    } else {
        // parse_args 已校验 threads ∈ [1, 1024]，usize::try_from 在该域下恒成功；
        // 仍按"全程不 panic"映射 Err。
        let threads_usize = usize::try_from(args.threads).map_err(|_| {
            format!("invalid --threads '{}': must be within [1, {MAX_THREADS}]", args.threads)
        })?;
        let report =
            solve_parallel_detailed(model, &cfg, threads_usize).map_err(|e| e.to_string())?;
        parallel_detail = Some((report.winner, report.worker_scores.clone()));
        report.best
    };

    if args.verbose {
        eprintln!(
            "[verbose] initial    : hard = {}, objective = {}",
            result.initial_score.hard, result.initial_score.soft
        );
        eprintln!("[verbose] steps done : {}", result.steps);
        if let Some((winner, worker_scores)) = parallel_detail {
            eprintln!("[verbose] winner     : worker {winner} of {}", args.threads);
            // HardSoftScore 无 Display：逐字段格式化为 (h, s)。
            for (i, s) in worker_scores.iter().enumerate() {
                eprintln!("[verbose] worker {i:<3} : best = ({}, {})", s.hard, s.soft);
            }
        }
    }
    // stdout 恰两行，结束时输出；progress 只走 stderr。
    println!("hard = {}", result.best_score.hard);
    println!("objective = {}", result.best_score.soft);
    Ok(())
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    if let Err(msg) = parse_args(&argv).and_then(|args| run(&args)) {
        eprintln!("error: {msg}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn defaults_when_only_path_given() {
        let a = parse_args(&args(&["inst.txt"])).expect("bare path must parse");
        assert_eq!(a.path, "inst.txt");
        assert_eq!(a.steps, None);
        assert_eq!(a.time_limit, None);
        assert_eq!(a.seed, 42);
        assert_eq!(a.algo, Algo::Tabu);
        assert_eq!(a.construction, ConstructionHeuristic::FirstFit);
        assert_eq!(a.samples, 32);
        assert!(!a.verbose);
    }

    #[test]
    fn all_flags_in_any_order() {
        let a = parse_args(&args(&[
            "i.txt", "--steps", "1000", "--seed", "7", "--algo", "sa", "--construct", "bf",
            "--samples", "64", "--verbose",
        ]))
        .expect("full option set must parse");
        assert_eq!(a.steps, Some(1000));
        assert_eq!(a.time_limit, None);
        assert_eq!(a.seed, 7);
        assert_eq!(a.algo, Algo::Sa);
        assert_eq!(a.construction, ConstructionHeuristic::BestFit);
        assert_eq!(a.samples, 64);
        assert!(a.verbose);

        // 旗标与路径顺序无关
        let b = parse_args(&args(&["--verbose", "--algo", "la", "i.txt"])).expect("flags before path");
        assert_eq!(b.algo, Algo::La);
        assert_eq!(b.path, "i.txt");
        assert!(b.verbose);
    }

    #[test]
    fn steps_and_time_secs_are_mutually_exclusive() {
        let e = parse_args(&args(&["i", "--steps", "10", "--time-secs", "1.5"]))
            .expect_err("steps + time-secs must fail");
        assert_eq!(e, "--steps and --time-secs are mutually exclusive");
        let e = parse_args(&args(&["i", "--time-secs", "1.5", "--steps", "10"]))
            .expect_err("time-secs + steps must fail (either order)");
        assert_eq!(e, "--steps and --time-secs are mutually exclusive");
    }

    #[test]
    fn time_secs_numeric_safety() {
        // 非法：NaN / inf / 负 / 超 36000 / 非数字字面量 / 解析溢出为 inf
        for bad in ["NaN", "nan", "-1", "-0.001", "inf", "-inf", "Infinity", "1e309", "36000.1", "abc", ""] {
            let e = parse_args(&args(&["i", "--time-secs", bad]))
                .err()
                .unwrap_or_else(|| panic!("case {bad:?} must be rejected"));
            assert_eq!(e, "invalid --time-secs", "case {bad:?} -> {e}");
        }
        // 合法：0 / 0.0 / 1.5 / 36000 / 科学计数
        for good in ["0", "0.0", "1.5", "36000", "3.5e2"] {
            let a = parse_args(&args(&["i", "--time-secs", good]))
                .unwrap_or_else(|e| panic!("case {good:?} must parse: {e}"));
            assert!(a.time_limit.is_some(), "case {good:?}");
            assert_eq!(a.steps, None);
        }
        // 1.5 秒精确转换为 Duration
        let a = parse_args(&args(&["i", "--time-secs", "1.5"])).expect("1.5");
        assert_eq!(
            a.time_limit,
            Some(Duration::try_from_secs_f64(1.5).expect("1.5 must convert"))
        );
        // 0 秒合法（TimeLimit(0)：立即返回初始解）
        let a = parse_args(&args(&["i", "--time-secs", "0"])).expect("0");
        assert_eq!(a.time_limit, Some(Duration::ZERO));
    }

    #[test]
    fn integer_flag_bounds_and_literals() {
        // --steps：非负有限字面量，≤ 10^9
        for bad in ["x", "-1", "1.5", "1e3", "", "1000000001", "99999999999999999999999"] {
            assert!(
                parse_args(&args(&["i", "--steps", bad])).is_err(),
                "--steps {bad:?} must be rejected"
            );
        }
        assert_eq!(parse_args(&args(&["i", "--steps", "0"])).unwrap().steps, Some(0));
        assert_eq!(
            parse_args(&args(&["i", "--steps", "1000000000"])).unwrap().steps,
            Some(1_000_000_000)
        );

        // --samples：有限整数字面量，1 ≤ K ≤ 10^6（0 在解析期即拒绝）
        for bad in ["x", "-5", "2.5", "0", "1000001"] {
            assert!(
                parse_args(&args(&["i", "--samples", bad])).is_err(),
                "--samples {bad:?} must be rejected"
            );
        }
        assert_eq!(
            parse_args(&args(&["i", "--samples", "1"])).unwrap().samples,
            1
        );
        assert_eq!(
            parse_args(&args(&["i", "--samples", "1000000"])).unwrap().samples,
            1_000_000
        );

        // --seed：u64 全域
        for bad in ["x", "-1", "1.0", ""] {
            assert!(
                parse_args(&args(&["i", "--seed", bad])).is_err(),
                "--seed {bad:?} must be rejected"
            );
        }
        assert_eq!(parse_args(&args(&["i", "--seed", "0"])).unwrap().seed, 0);
        assert_eq!(
            parse_args(&args(&["i", "--seed", "18446744073709551615"])).unwrap().seed,
            u64::MAX
        );
    }

    /// P6 回归：--samples 0 在**解析期**即拒绝（消息注明须 ≥1），不再拖到运行期由
    /// lib 侧 `samples_per_step must be >= 1` 校验兜底；声明域与强制域一致。
    #[test]
    fn samples_zero_rejected_at_parse_with_ge1_message() {
        let e = parse_args(&args(&["i", "--samples", "0"]))
            .err()
            .unwrap_or_else(|| panic!("--samples 0 must be rejected at parse time"));
        assert_eq!(e, "invalid --samples '0': must be >= 1");
        // 解析失败消息的取值域说明同步为 [1, 10^6]
        let e = parse_args(&args(&["i", "--samples", "x"]))
            .err()
            .unwrap_or_else(|| panic!("--samples x must be rejected"));
        assert_eq!(e, "invalid --samples 'x' (expected integer in [1, 10^6])");
    }

    #[test]
    fn algo_construct_and_usage_errors() {        assert!(parse_args(&args(&["i", "--algo", "SA"])).is_err(), "大小写敏感");
        assert!(parse_args(&args(&["i", "--algo", "greedy"])).is_err());
        assert!(parse_args(&args(&["i", "--construct", "zz"])).is_err());
        assert!(parse_args(&args(&["i", "--unknown"])).is_err(), "未知旗标");
        assert!(parse_args(&args(&["i", "--steps"])).is_err(), "旗标缺值");
        assert!(parse_args(&args(&[])).is_err(), "缺算例路径");
        assert!(parse_args(&args(&["a", "b"])).is_err(), "多余位置参数");
        assert!(parse_args(&args(&["i", "--steps", "1", "--steps", "2"])).is_err(), "重复旗标");
    }

    #[test]
    fn run_reports_missing_file_as_error_without_panic() {
        // IO 错误路径：run 返回 Err（main 打印 error: + exit 1），全程不 panic。
        let a = parse_args(&args(&["no-such-instance-file.txt"])).expect("parses");
        match run(&a) {
            Err(msg) => assert!(
                msg.contains("io error") || msg.contains("os error"),
                "unexpected message: {msg}"
            ),
            Ok(()) => panic!("missing file must be an error"),
        }
    }

    // ---- --threads（wp4 新增；上方 7 个既有用例逐字未动）----

    #[test]
    fn threads_defaults_to_one() {
        let a = parse_args(&args(&["inst.txt"])).expect("bare path must parse");
        assert_eq!(a.threads, 1, "default threads must be 1 (single-thread path)");
        // 显式 1 同样合法（仍是 solve() 路径）
        let a = parse_args(&args(&["inst.txt", "--threads", "1"])).expect("--threads 1");
        assert_eq!(a.threads, 1);
        // 合法上界 1024
        let a = parse_args(&args(&["inst.txt", "--threads", "1024"])).expect("1024 in bounds");
        assert_eq!(a.threads, 1024);
    }

    #[test]
    fn threads_combines_with_other_flags_in_any_order() {
        // 与 --steps/--seed/--algo 任意顺序组合（无互斥）
        let a = parse_args(&args(&[
            "i.txt", "--threads", "4", "--steps", "1000", "--seed", "7", "--algo", "sa",
        ]))
        .expect("--threads first must parse");
        assert_eq!(a.threads, 4);
        assert_eq!(a.steps, Some(1000));
        assert_eq!(a.seed, 7);
        assert_eq!(a.algo, Algo::Sa);

        let b = parse_args(&args(&[
            "--algo", "la", "--seed", "9", "--steps", "50", "--threads", "8", "i.txt",
        ]))
        .expect("--threads last must parse");
        assert_eq!(b.algo, Algo::La);
        assert_eq!(b.seed, 9);
        assert_eq!(b.steps, Some(50));
        assert_eq!(b.threads, 8);
        assert_eq!(b.path, "i.txt");

        // 与 --time-secs 也可组合（--steps 与 --time-secs 互斥规则不变）
        let c = parse_args(&args(&["i.txt", "--threads", "4", "--time-secs", "1.5"]))
            .expect("--threads + --time-secs");
        assert_eq!(c.threads, 4);
        assert!(c.time_limit.is_some());
    }

    #[test]
    fn threads_rejects_invalid_literals() {
        // 0（低于下界）/ 1.5 / x（非整数）/ 1025（超上界）/ u64 溢出字面量
        for bad in ["0", "1.5", "x", "1025", "99999999999999999999999"] {
            assert!(
                parse_args(&args(&["i", "--threads", bad])).is_err(),
                "--threads {bad:?} must be rejected"
            );
        }
        // 错误消息风格与 --steps 一致（parse 失败 / 越界两类）
        let e = parse_args(&args(&["i", "--threads", "x"]))
            .err()
            .unwrap_or_else(|| panic!("'x' must be rejected"));
        assert_eq!(e, "invalid --threads 'x' (expected integer in [1, 1024])");
        let e = parse_args(&args(&["i", "--threads", "0"]))
            .err()
            .unwrap_or_else(|| panic!("'0' must be rejected"));
        assert_eq!(e, "invalid --threads '0': must be within [1, 1024]");
        let e = parse_args(&args(&["i", "--threads", "1025"]))
            .err()
            .unwrap_or_else(|| panic!("'1025' must be rejected"));
        assert_eq!(e, "invalid --threads '1025': must be within [1, 1024]");
    }

    #[test]
    fn threads_duplicate_and_missing_value_rejected() {
        let e = parse_args(&args(&["i", "--threads", "4", "--threads", "2"]))
            .err()
            .unwrap_or_else(|| panic!("duplicate --threads must be rejected"));
        assert_eq!(e, "duplicate --threads");
        // 缺值：flag_value 返回 Err
        assert!(
            parse_args(&args(&["i", "--threads"])).is_err(),
            "--threads without a value must be rejected"
        );
    }

    /// CLI 侧自查：本文件 MAX_THREADS 与 lib 侧 parallel::MAX_THREADS 数值一致。
    #[test]
    fn threads_constants_match_lib() {
        assert_eq!(MAX_THREADS as usize, planner::parallel::MAX_THREADS);
        assert_eq!(DEFAULT_THREADS, 1);
    }
}
