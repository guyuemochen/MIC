//! milp CLI 入口（WP5 契约实现）。
//!
//! 用法：`milp <算例路径> [--time-limit-secs N] [--verbose]`。
//!
//! 输出契约：Optimal 时 stdout 最后一行必须精确是
//! `objective = <format_objective(v)>`；此前输出 `status = optimal`、
//! `nodes = N`，`--verbose` 时逐行 `name = value`。Infeasible/Unbounded/
//! LimitReached 输出 `status = ...`（无 objective 行），LimitReached 且有
//! incumbent 时附 `best = <v>`。错误到 stderr `error: <MilpError>`。
//! 退出码：0=最优 2=不可行 3=无界 4=上限触发 1=错误/用法。
//!
//! 实现纪律：用户输入路径禁止 unwrap/expect/panic；`--time-limit-secs N`
//! 映射 `BbOptions::time_limit_ms` 使用 `N.saturating_mul(1000)`（裸乘法在
//! debug 极值下溢出 panic，违反铁律）；写 stdout 失败（如管道断裂）静默
//! 忽略而非 panic。

use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

use milp::branch_bound::{BbOptions, BranchRule, NodeRule};
use milp::instance::load_instance;
use milp::model::{Problem, SolveStatus, Tolerances};
use milp::solver::{solve, SolveOptions, SolveReport};

/// 用法说明（打印到 stderr）。
const USAGE: &str = "用法：milp <算例路径> [--time-limit-secs N] [--verbose]";

/// 命令行配置。
#[derive(Debug)]
struct Config {
    /// 算例文件路径。
    instance_path: String,
    /// 时间上限（秒）；映射 BbOptions.time_limit_ms 必须用
    /// `N.saturating_mul(1000)`（禁止裸乘法：debug 极值溢出 panic 违反
    /// “输入路径禁 panic”铁律）。
    time_limit_secs: Option<u64>,
    /// 是否输出详细变量取值。
    verbose: bool,
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let config = match parse_args(&args) {
        Ok(config) => config,
        Err(message) => {
            eprintln!("error: {message}");
            return ExitCode::from(1);
        }
    };
    let problem = match load_instance(Path::new(&config.instance_path)) {
        Ok(problem) => problem,
        Err(error) => {
            eprintln!("error: {error}");
            return ExitCode::from(1);
        }
    };
    let options = solve_options_from(&config);
    let report = match solve(&problem, &options) {
        Ok(report) => report,
        Err(error) => {
            eprintln!("error: {error}");
            return ExitCode::from(1);
        }
    };
    // 安全网：Optimal 报告必须同时携带目标值与变量取值（诚实性要求，
    // 求解器侧违约时不冒充成功，按错误退出）。
    if report.status.is_optimal() && (report.objective.is_none() || report.values.is_none()) {
        eprintln!("error: 求解器报告 Optimal 但缺少目标值/变量取值（内部不一致）");
        return ExitCode::from(1);
    }
    print_report(&report, &problem, config.verbose);
    exit_code_for(report.status)
}

/// 解析命令行参数。
///
/// 用法：`milp <算例路径> [--time-limit-secs N] [--verbose]`；
/// 无参打印用法到 stderr 并返回 `Err`；`N` 非数字或为负 → `Err`。
fn parse_args(args: &[String]) -> Result<Config, String> {
    if args.is_empty() {
        eprintln!("{USAGE}");
        return Err("缺少算例路径".to_string());
    }
    let mut instance_path: Option<String> = None;
    let mut time_limit_secs: Option<u64> = None;
    let mut verbose = false;
    let mut index = 0usize;
    while index < args.len() {
        let arg = args[index].as_str();
        if arg == "--verbose" {
            verbose = true;
        } else if arg == "--time-limit-secs" {
            if time_limit_secs.is_some() {
                return Err(format!("`--time-limit-secs` 重复出现（{USAGE}）"));
            }
            let Some(value) = args.get(index + 1) else {
                return Err(format!("`--time-limit-secs` 缺少数值参数 N（{USAGE}）"));
            };
            let secs = match value.parse::<u64>() {
                Ok(secs) => secs,
                Err(_) => {
                    let reason = if value.starts_with('-') {
                        "为负数"
                    } else {
                        "无法解析为非负整数"
                    };
                    return Err(format!(
                        "`--time-limit-secs` 的参数 `{value}` {reason}（期望非负整数秒数）"
                    ));
                }
            };
            time_limit_secs = Some(secs);
            index += 1; // 额外消费数值 token
        } else if arg.starts_with("--") {
            return Err(format!("未知参数 `{arg}`（{USAGE}）"));
        } else if instance_path.is_some() {
            return Err(format!("多余的位置参数 `{arg}`（{USAGE}）"));
        } else {
            instance_path = Some(arg.to_string());
        }
        index += 1;
    }
    let Some(instance_path) = instance_path else {
        return Err(format!("缺少算例路径（{USAGE}）"));
    };
    Ok(Config {
        instance_path,
        time_limit_secs,
        verbose,
    })
}

/// 由 [`Config`] 构造求解选项。
///
/// 数值取契约默认值（`Tolerances::default()` 与 `BbOptions::default()` 的
/// 文档值），仅 `time_limit_ms` 按 CLI 覆盖；映射必须用
/// `N.saturating_mul(1000)`，`u64::MAX` 秒饱和为 `u64::MAX` 毫秒而非溢出。
fn solve_options_from(config: &Config) -> SolveOptions {
    SolveOptions {
        tol: Tolerances::default(),
        bb: BbOptions {
            branch_rule: BranchRule::MostFractional,
            node_rule: NodeRule::DepthFirst,
            node_limit: 100_000,
            time_limit_ms: config
                .time_limit_secs
                .map_or(0, |secs| secs.saturating_mul(1000)),
            root_cuts: true,
            heuristics: true,
            dive_max_fixes: 8,
            prune_abs: 1e-9,
            prune_rel: 1e-9,
        },
        use_presolve: true,
    }
}

/// 目标值格式化：0<|v|<1e-6 或 |v|≥1e15 时 `format!("{:e}", v)`，
/// 否则 `format!("{:.6}", v)`。
fn format_objective(v: f64) -> String {
    let magnitude = v.abs();
    if (0.0 < magnitude && magnitude < 1e-6) || magnitude >= 1e15 {
        format!("{v:e}")
    } else {
        format!("{v:.6}")
    }
}

/// 状态 → 输出 token。
fn status_token(status: SolveStatus) -> &'static str {
    match status {
        SolveStatus::Optimal => "optimal",
        SolveStatus::Infeasible => "infeasible",
        SolveStatus::Unbounded => "unbounded",
        SolveStatus::LimitReached => "limit_reached",
    }
}

/// 按输出契约打印求解报告（Optimal 时最后一行精确为
/// `objective = <format_objective(v)>`）。
fn print_report(report: &SolveReport, problem: &Problem, verbose: bool) {
    let mut stdout = std::io::stdout().lock();
    write_report(&mut stdout, report, problem, verbose);
}

/// [`print_report`] 的可测试内核：写入任意 `Write`（写失败静默忽略，不 panic）。
fn write_report<W: Write>(out: &mut W, report: &SolveReport, problem: &Problem, verbose: bool) {
    let _ = writeln!(out, "status = {}", status_token(report.status));
    if report.status.is_optimal() {
        let _ = writeln!(out, "nodes = {}", report.stats.nodes);
        if verbose {
            if let Some(values) = &report.values {
                for (variable, value) in problem.vars.iter().zip(values.iter()) {
                    let _ = writeln!(out, "{} = {}", variable.name, format_objective(*value));
                }
            }
        }
        // 契约：最后一行必须精确是 objective 行（main 已安全网保证 Some）。
        if let Some(objective) = report.objective {
            let _ = writeln!(out, "objective = {}", format_objective(objective));
        }
    } else if report.status == SolveStatus::LimitReached {
        if let Some(best) = report.objective {
            let _ = writeln!(out, "best = {}", format_objective(best));
        }
    }
}

/// 终态 → 退出码：0=最优 2=不可行 3=无界 4=上限触发。
fn exit_code_for(status: SolveStatus) -> ExitCode {
    match status {
        SolveStatus::Optimal => ExitCode::SUCCESS,
        SolveStatus::Infeasible => ExitCode::from(2),
        SolveStatus::Unbounded => ExitCode::from(3),
        SolveStatus::LimitReached => ExitCode::from(4),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use milp::branch_bound::BbStats;
    use milp::model::{Sense, VarType};

    /// 构造测试用报告。
    fn report(
        status: SolveStatus,
        objective: Option<f64>,
        values: Option<Vec<f64>>,
        nodes: u64,
    ) -> SolveReport {
        SolveReport {
            status,
            objective,
            values,
            stats: BbStats {
                nodes,
                ..BbStats::default()
            },
            presolve_reductions: 0,
        }
    }

    /// 构造两个变量的测试问题（x1 integer、x2 continuous）。
    fn two_var_problem() -> Problem {
        let mut problem = Problem::new("t", Sense::Maximize);
        problem.add_var("x1", VarType::Integer).unwrap();
        problem.add_var("x2", VarType::Continuous).unwrap();
        problem
            .set_objective(vec![(0, 3.0), (1, 2.0)], 0.0)
            .unwrap();
        problem
    }

    /// args 便捷构造。
    fn args(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    // —— parse_args ——

    /// 仅算例路径：默认无时间上限、非 verbose。
    #[test]
    fn parse_args_path_only() {
        let config = parse_args(&args(&["a.milp"])).unwrap();
        assert_eq!(config.instance_path, "a.milp");
        assert_eq!(config.time_limit_secs, None);
        assert!(!config.verbose);
    }

    /// 全部合法组合（含 flags 先于路径的顺序）。
    #[test]
    fn parse_args_all_flags() {
        let config = parse_args(&args(&["a.milp", "--time-limit-secs", "5", "--verbose"]))
            .expect("合法组合必须解析成功");
        assert_eq!(config.instance_path, "a.milp");
        assert_eq!(config.time_limit_secs, Some(5));
        assert!(config.verbose);

        let config = parse_args(&args(&["--verbose", "a.milp"])).unwrap();
        assert_eq!(config.instance_path, "a.milp");
        assert!(config.verbose);
        assert_eq!(config.time_limit_secs, None);

        // 0 秒合法（time_limit_ms=0 即“不限”）。
        assert_eq!(
            parse_args(&args(&["a.milp", "--time-limit-secs", "0"]))
                .unwrap()
                .time_limit_secs,
            Some(0)
        );
        // u64 上限合法（映射毫秒时由 saturating_mul 饱和）。
        assert_eq!(
            parse_args(&args(&[
                "a.milp",
                "--time-limit-secs",
                &u64::MAX.to_string()
            ]))
            .unwrap()
            .time_limit_secs,
            Some(u64::MAX)
        );
    }

    /// 无参 → Err（用法已打印到 stderr，此处只断言 Err 与信息）。
    #[test]
    fn parse_args_empty_rejected() {
        match parse_args(&[]) {
            Err(message) => assert!(message.contains("缺少算例路径"), "{message}"),
            Ok(config) => panic!("无参应报错，得到 {config:?}"),
        }
    }

    /// 各类非法参数一律 Err：N 缺失/非数字/为负、未知 flag、多余位置参数、
    /// 重复 --time-limit-secs、只有 flag 无路径。
    #[test]
    fn parse_args_invalid_rejected() {
        let bad: &[&[&str]] = &[
            &["--time-limit-secs"],                                          // N 缺失
            &["a.milp", "--time-limit-secs"],                                // N 缺失
            &["a.milp", "--time-limit-secs", "abc"],                         // 非数字
            &["a.milp", "--time-limit-secs", "3.5"],                         // 非整数
            &["a.milp", "--time-limit-secs", "-5"],                          // 负数
            &["a.milp", "--time-limit-secs", "5", "--time-limit-secs", "7"], // 重复
            &["a.milp", "--bogus"],                                          // 未知 flag
            &["a.milp", "b.milp"],                                           // 多余位置参数
            &["--verbose"],                                                  // 无路径
            &["--time-limit-secs", "5"],                                     // 无路径
        ];
        for case in bad {
            match parse_args(&args(case)) {
                Err(message) => assert!(!message.is_empty()),
                Ok(config) => panic!("参数 {:?} 应报错，得到 {config:?}", case),
            }
        }
    }

    /// N 为负的错误信息应说明“负数”（与“无法解析”区分）。
    #[test]
    fn parse_args_negative_message() {
        let message =
            parse_args(&args(&["a.milp", "--time-limit-secs", "-5"])).expect_err("负数应报错");
        assert!(message.contains("负数"), "{message}");
    }

    // —— time_limit_secs → time_limit_ms 映射 ——

    /// 映射必须 saturating：u64::MAX 秒 → u64::MAX 毫秒（不 panic、不回绕）。
    #[test]
    fn time_limit_mapping_saturates() {
        let config = Config {
            instance_path: "a.milp".to_string(),
            time_limit_secs: None,
            verbose: false,
        };
        assert_eq!(solve_options_from(&config).bb.time_limit_ms, 0);
        let config = Config {
            instance_path: "a.milp".to_string(),
            time_limit_secs: Some(2),
            verbose: false,
        };
        assert_eq!(solve_options_from(&config).bb.time_limit_ms, 2000);
        let config = Config {
            instance_path: "a.milp".to_string(),
            time_limit_secs: Some(u64::MAX),
            verbose: false,
        };
        assert_eq!(
            solve_options_from(&config).bb.time_limit_ms,
            u64::MAX,
            "u64::MAX 秒 × 1000 必须饱和为 u64::MAX 而非回绕"
        );
    }

    // —— format_objective ——

    /// 常规量级用 {:.6}；科学记数触发边界（0<|v|<1e-6 或 |v|≥1e15）。
    #[test]
    fn format_objective_branches() {
        assert_eq!(format_objective(0.0), "0.000000");
        assert_eq!(format_objective(13.0), "13.000000");
        assert_eq!(format_objective(-2.5), "-2.500000");
        assert_eq!(format_objective(1e-6), "0.000001", "1e-6 本身不走科学记数");
        assert_eq!(format_objective(9.999999e-7), "9.999999e-7");
        assert_eq!(format_objective(-1e-7), "-1e-7");
        assert_eq!(format_objective(9.99e14), "999000000000000.000000");
        assert_eq!(format_objective(1e15), "1e15");
        assert_eq!(format_objective(-2.5e18), "-2.5e18");
    }

    // —— write_report / print_report 输出契约 ——

    /// Optimal 非 verbose：status、nodes、objective 三行，objective 精确为末行。
    #[test]
    fn report_optimal_last_line_is_objective() {
        let problem = two_var_problem();
        let report = report(SolveStatus::Optimal, Some(17.0), Some(vec![0.0, 2.0]), 3);
        let mut out = Vec::new();
        write_report(&mut out, &report, &problem, false);
        let text = String::from_utf8(out).unwrap();
        assert_eq!(
            text, "status = optimal\nnodes = 3\nobjective = 17.000000\n",
            "输出必须逐行精确"
        );
        assert_eq!(text.lines().last(), Some("objective = 17.000000"));
    }

    /// Optimal + verbose：逐行 name = value 出现在 objective 行之前。
    #[test]
    fn report_optimal_verbose_values() {
        let problem = two_var_problem();
        let report = report(SolveStatus::Optimal, Some(5.5), Some(vec![1.0, 2.25]), 7);
        let mut out = Vec::new();
        write_report(&mut out, &report, &problem, true);
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "status = optimal");
        assert_eq!(lines[1], "nodes = 7");
        assert_eq!(lines[2], "x1 = 1.000000");
        assert_eq!(lines[3], "x2 = 2.250000");
        assert_eq!(lines.last(), Some(&"objective = 5.500000"));
    }

    /// 非最优终态：只有 status 行（Infeasible/Unbounded/无 incumbent 的
    /// LimitReached），绝无 objective 行。
    #[test]
    fn report_non_optimal_no_objective_line() {
        let problem = two_var_problem();
        for (status, token) in [
            (SolveStatus::Infeasible, "infeasible"),
            (SolveStatus::Unbounded, "unbounded"),
        ] {
            let report = report(status, None, None, 1);
            let mut out = Vec::new();
            write_report(&mut out, &report, &problem, false);
            let text = String::from_utf8(out).unwrap();
            assert_eq!(text, format!("status = {token}\n"), "{status:?}");
            assert!(!text.contains("objective"));
        }
        // LimitReached 无 incumbent：仅 status 行。
        let report = report(SolveStatus::LimitReached, None, None, 100_000);
        let mut out = Vec::new();
        write_report(&mut out, &report, &problem, false);
        assert_eq!(String::from_utf8(out).unwrap(), "status = limit_reached\n");
    }

    /// LimitReached 有 incumbent：status 行后附 best 行，无 objective 行。
    #[test]
    fn report_limit_reached_with_incumbent_prints_best() {
        let problem = two_var_problem();
        let report = report(
            SolveStatus::LimitReached,
            Some(12.5),
            Some(vec![1.0, 4.0]),
            100_000,
        );
        let mut out = Vec::new();
        write_report(&mut out, &report, &problem, false);
        let text = String::from_utf8(out).unwrap();
        assert_eq!(text, "status = limit_reached\nbest = 12.500000\n");
        assert!(!text.contains("objective"));
    }

    /// 退出码映射：0/2/3/4。
    #[test]
    fn exit_codes_map() {
        assert_eq!(
            exit_code_for(SolveStatus::Optimal),
            ExitCode::SUCCESS,
            "Optimal → 0"
        );
        assert_eq!(exit_code_for(SolveStatus::Infeasible), ExitCode::from(2));
        assert_eq!(exit_code_for(SolveStatus::Unbounded), ExitCode::from(3));
        assert_eq!(exit_code_for(SolveStatus::LimitReached), ExitCode::from(4));
    }

    // —— 本包交付的 7 个算例文件可解析（WP1 已实现，纯解析层面即可验证）——

    /// 7 个算例全部可被 load_instance 解析且模型校验通过、算例名正确。
    #[test]
    fn all_seven_instances_parse() {
        for name in [
            "lp-prod-3",
            "knapsack-10",
            "assign-4x4",
            "milp-wolsey-2",
            "milp-mixed-5",
            "infeasible-tiny",
            "unbounded-tiny",
        ] {
            let path = format!("{}/instances/{name}.milp", env!("CARGO_MANIFEST_DIR"));
            match load_instance(Path::new(&path)) {
                Ok(problem) => {
                    assert_eq!(problem.name, name, "算例名应与文件名一致");
                    assert!(problem.validate().is_ok(), "{name} 校验失败");
                }
                Err(error) => panic!("算例 {name} 解析失败：{error}"),
            }
        }
    }
}
