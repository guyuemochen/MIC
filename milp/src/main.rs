//! milp CLI 入口（WP5 契约实现 + 受控搜索参数包）。
//!
//! 用法：`milp <算例路径> [--time-limit-secs N] [--threads N] [--log-secs F]
//! [--mipgap P] [--heuristic-mode auto|on|off] `
//! [--branch-rule most_fractional|pseudocost] [--scaling=on|off] [--mitm=on|off] [--work-stack local|global] [--verbose]`。
//!
//! 输出契约：
//! - **Optimal**：stdout 依次 `status = optimal`、`nodes = N`（`--verbose` 时
//!   逐行 `name = value`），**最后一行精确为** `objective = <format_objective(v)>`；
//! - **Infeasible / Unbounded**：仅 `status = ...`（无 objective 行）；
//! - **LimitReached**（时间/节点上限触发）：`status = limit_reached`，有
//!   incumbent 时附 `best = <v>`（无 objective 行）；
//! - **gap 终止**（`--mipgap P` 且 P>0 触发）：依次 `status = gap_reached`、
//!   `nodes = N`、`gap = <%.2f>`（库层触发时刻证据 + 统一 clamp 双保证为
//!   非负有限，永不 inf/负）、`best = <数值>`（**无 objective 行**——未证明
//!   最优，绝不冒充）；
//! - **自动降级终止**（`--heuristic-mode auto` 且给了 `--time-limit-secs N`，
//!   预算过半后解与对偶界双停滞时触发）：依次 `status = auto_degraded`、
//!   `nodes = N`、`gap = <%.2f>`、`best = <数值>`（行形状同 gap 终止、
//!   status 区分原因；同样**无 objective 行**——未证明最优）。
//!
//! 错误到 stderr `error: <MilpError>`。退出码：0=最优 2=不可行 3=无界
//! 4=上限触发与 gap 终止（均未证明最优） 1=错误/用法。
//!
//! 实现纪律：用户输入路径禁止 unwrap/expect/panic；`--time-limit-secs N`
//! 映射 `BbOptions::time_limit_ms` 使用 `N.saturating_mul(1000)`（裸乘法在
//! debug 极值下溢出 panic，违反铁律）；写 stdout 失败（如管道断裂）静默
//! 忽略而非 panic。

use std::io::Write;
use std::path::Path;
use std::process::ExitCode;
use std::sync::Arc;

use milp::branch_bound::{BbOptions, BranchRule, NodeRule, SearchControl, SearchHooks};
use milp::instance::load_instance;
use milp::model::{Problem, SolveStatus, Tolerances};
use milp::search_ctrl::{format_progress_line, ProgressInfo};
use milp::solver::{solve_controlled, ControlledSolveReport, SolveOptions, SolveReport};

/// 用法说明（打印到 stderr）。
const USAGE: &str = "用法：milp <算例路径> [--time-limit-secs N] [--threads N] \
[--log-secs F] [--mipgap P] [--heuristic-mode auto|on|off] [--branch-rule most_fractional|pseudocost] [--scaling=on|off] [--mitm=on|off] [--verbose]";

/// `--threads` 允许的最大线程数（1..=1024；1 = 纯顺序路径）。
const MAX_THREADS: usize = 1024;

/// 启发式模式（`--heuristic-mode`；默认 Auto）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HeuristicModeCli {
    /// auto：heuristics=true 且自适应（控制器按需在树节点补跑 dive）。
    Auto,
    /// on：heuristics=true 且不自适应（升级前行为：仅根启发式）。
    On,
    /// off：heuristics=false。
    Off,
}

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
    /// 线程数（`--threads N`，1..=1024）；默认 1 = 纯顺序路径（输出与
    /// 升级前逐字节一致），≥2 走并行引擎。
    threads: usize,
    /// 进度行周期秒数（`--log-secs F`，有限且 F>0）；默认 None = 不输出
    /// 任何进度行（保持既有 stdout 契约逐字节稳定）。
    log_secs: Option<f64>,
    /// MIP gap 百分比阈值（`--mipgap P`，有限且 P≥0）；默认 None = 关闭；
    /// Some(0.0) 等价关闭。
    mipgap_percent: Option<f64>,
    /// 启发式模式；默认 Auto。
    heuristic_mode: HeuristicModeCli,
    /// 分支规则（`--branch-rule`）；默认 None = 类自适应（R3：求解入口按
    /// 结构类检测解析——全 Binary + 全正系数 Le 行 + 行数 ≤ 64 ⇒ Pseudocost
    /// + 可靠性冷启动 8，否则 MostFractional）；显式 most_fractional/
    /// pseudocost 覆盖优先（pseudocost 含可靠性冷启动，`--pc-probes` 仅
    /// bench 探针暴露）。
    branch_rule: Option<BranchRule>,
    /// 数值缩放（`--scaling=on|off`，R4）；默认 on（行列迭代几何均值，
    /// 2 的幂系数；目标与解还原原尺度，输出契约不变）。
    scaling: bool,
    /// 子集和 MITM 可行性启发式（`--mitm=on|off`，R5-retry）；默认 on
    ///（全 Eq + 全 Binary 类的最后手段，廉价启发式 miss 时根一次求解）。
    mitm: bool,
    /// 并行工作栈模式（`--work-stack local|global`，R6-retry）；默认 local
    ///（本地亲和，现状）；global = 全局 DFS 栈（并行探索序对齐顺序 DFS，
    /// 抑制搜索膨胀；仅 threads ≥ 2 生效）。
    work_stack_global: bool,
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
    let ctrl = search_control_from(&config);
    // 进度回调仅在显式给出 --log-secs 时注入（默认零进度行，stdout 契约
    // 逐字节稳定）。回调内部自行短暂持有 stdout 锁，单行写入、失败静默。
    let hooks = config.log_secs.map(|_| progress_hooks());
    let controlled = match solve_controlled(&problem, &options, &ctrl, hooks.as_ref()) {
        Ok(controlled) => controlled,
        Err(error) => {
            eprintln!("error: {error}");
            return ExitCode::from(1);
        }
    };
    let report = &controlled.report;
    // 安全网：Optimal 报告必须同时携带目标值与变量取值（诚实性要求，
    // 求解器侧违约时不冒充成功，按错误退出）。
    if report.status.is_optimal() && (report.objective.is_none() || report.values.is_none()) {
        eprintln!("error: 求解器报告 Optimal 但缺少目标值/变量取值（内部不一致）");
        return ExitCode::from(1);
    }
    // gap 终止：独立输出契约（status = gap_reached / nodes / gap / best，
    // 无 objective 行），退出码 4。
    if controlled.gap_terminated {
        return gap_report_and_exit(&controlled);
    }
    print_report(report, &problem, config.verbose);
    exit_code_for(report.status)
}

/// 解析命令行参数。
///
/// 用法见 [`USAGE`]；无参打印用法到 stderr 并返回 `Err`；各数值参数缺失/
/// 非法（非整数、越界、非有限、非正、未知 token）或重复出现一律 `Err`
/// （退出码 1），风格与错误信息与既有 `--time-limit-secs` 分支一致。
fn parse_args(args: &[String]) -> Result<Config, String> {
    if args.is_empty() {
        eprintln!("{USAGE}");
        return Err("缺少算例路径".to_string());
    }
    let mut instance_path: Option<String> = None;
    let mut time_limit_secs: Option<u64> = None;
    let mut verbose = false;
    let mut threads: Option<usize> = None;
    let mut log_secs: Option<f64> = None;
    let mut mipgap_percent: Option<f64> = None;
    let mut heuristic_mode: Option<HeuristicModeCli> = None;
    let mut branch_rule: Option<BranchRule> = None;
    let mut scaling: Option<bool> = None;
    let mut mitm: Option<bool> = None;
    let mut work_stack_global: Option<bool> = None;
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
        } else if arg == "--threads" {
            if threads.is_some() {
                return Err(format!("`--threads` 重复出现（{USAGE}）"));
            }
            let Some(value) = args.get(index + 1) else {
                return Err(format!("`--threads` 缺少数值参数 N（{USAGE}）"));
            };
            let n = match value.parse::<usize>() {
                Ok(n) => n,
                Err(_) => {
                    let reason = if value.starts_with('-') {
                        "为负数"
                    } else {
                        "无法解析为整数"
                    };
                    return Err(format!(
                        "`--threads` 的参数 `{value}` {reason}（期望 1..={MAX_THREADS} 的整数）"
                    ));
                }
            };
            if n == 0 || n > MAX_THREADS {
                return Err(format!(
                    "`--threads` 的参数 `{value}` 超出范围（期望 1..={MAX_THREADS} 的整数）"
                ));
            }
            threads = Some(n);
            index += 1; // 额外消费数值 token
        } else if arg == "--log-secs" {
            if log_secs.is_some() {
                return Err(format!("`--log-secs` 重复出现（{USAGE}）"));
            }
            let Some(value) = args.get(index + 1) else {
                return Err(format!("`--log-secs` 缺少数值参数 F（{USAGE}）"));
            };
            let secs = match value.parse::<f64>() {
                Ok(secs) => secs,
                Err(_) => {
                    return Err(format!(
                        "`--log-secs` 的参数 `{value}` 无法解析为数值（期望有限且大于 0 的秒数）"
                    ));
                }
            };
            if !secs.is_finite() || secs <= 0.0 {
                return Err(format!(
                    "`--log-secs` 的参数 `{value}` 必须为有限且大于 0 的秒数（{USAGE}）"
                ));
            }
            log_secs = Some(secs);
            index += 1; // 额外消费数值 token
        } else if arg == "--mipgap" {
            if mipgap_percent.is_some() {
                return Err(format!("`--mipgap` 重复出现（{USAGE}）"));
            }
            let Some(value) = args.get(index + 1) else {
                return Err(format!("`--mipgap` 缺少数值参数 P（{USAGE}）"));
            };
            let percent = match value.parse::<f64>() {
                Ok(percent) => percent,
                Err(_) => {
                    return Err(format!(
                        "`--mipgap` 的参数 `{value}` 无法解析为数值（期望 ≥0 的百分比）"
                    ));
                }
            };
            if !percent.is_finite() || percent < 0.0 {
                return Err(format!(
                    "`--mipgap` 的参数 `{value}` 必须为有限且 ≥0 的百分比（{USAGE}）"
                ));
            }
            mipgap_percent = Some(percent);
            index += 1; // 额外消费数值 token
        } else if arg == "--heuristic-mode" {
            if heuristic_mode.is_some() {
                return Err(format!("`--heuristic-mode` 重复出现（{USAGE}）"));
            }
            let Some(value) = args.get(index + 1) else {
                return Err(format!(
                    "`--heuristic-mode` 缺少参数 auto|on|off（{USAGE}）"
                ));
            };
            let mode = match value.as_str() {
                "auto" => HeuristicModeCli::Auto,
                "on" => HeuristicModeCli::On,
                "off" => HeuristicModeCli::Off,
                unknown => {
                    return Err(format!(
                        "`--heuristic-mode` 的参数 `{unknown}` 未知（期望 auto|on|off）（{USAGE}）"
                    ));
                }
            };
            heuristic_mode = Some(mode);
            index += 1; // 额外消费模式 token
        } else if arg == "--branch-rule" {
            if branch_rule.is_some() {
                return Err(format!("`--branch-rule` 重复出现（{USAGE}）"));
            }
            let Some(value) = args.get(index + 1) else {
                return Err(format!(
                    "`--branch-rule` 缺少参数 most_fractional|pseudocost（{USAGE}）"
                ));
            };
            let rule = match value.as_str() {
                "most_fractional" => BranchRule::MostFractional,
                "pseudocost" => BranchRule::Pseudocost,
                other => {
                    return Err(format!(
                        "未知 --branch-rule `{other}`（期望 most_fractional|pseudocost，{USAGE}）"
                    ))
                }
            };
            branch_rule = Some(rule);
            index += 1; // 额外消费规则 token
        } else if arg == "--scaling" || arg.starts_with("--scaling=") {
            if scaling.is_some() {
                return Err(format!("`--scaling` 重复出现（{USAGE}）"));
            }
            // 兼容 `--scaling=on|off`（等号连写，ask 原文）与
            // `--scaling on|off`（空格）两种形式。
            let value: String = if let Some(v) = arg.strip_prefix("--scaling=") {
                v.to_string()
            } else {
                let Some(v) = args.get(index + 1) else {
                    return Err(format!("`--scaling` 缺少参数 on|off（{USAGE}）"));
                };
                index += 1; // 额外消费 on/off token
                v.clone()
            };
            let v = match value.as_str() {
                "on" => true,
                "off" => false,
                other => {
                    return Err(format!(
                        "未知 --scaling 值 `{other}`（期望 on|off，{USAGE}）"
                    ))
                }
            };
            scaling = Some(v);
        } else if arg == "--mitm" || arg.starts_with("--mitm=") {
            if mitm.is_some() {
                return Err(format!("`--mitm` 重复出现（{USAGE}）"));
            }
            let value: String = if let Some(v) = arg.strip_prefix("--mitm=") {
                v.to_string()
            } else {
                let Some(v) = args.get(index + 1) else {
                    return Err(format!("`--mitm` 缺少参数 on|off（{USAGE}）"));
                };
                index += 1;
                v.clone()
            };
            let v = match value.as_str() {
                "on" => true,
                "off" => false,
                other => {
                    return Err(format!(
                        "未知 --mitm 值 `{other}`（期望 on|off，{USAGE}）"
                    ))
                }
            };
            mitm = Some(v);
        } else if arg == "--work-stack" || arg.starts_with("--work-stack=") {
            if work_stack_global.is_some() {
                return Err(format!("`--work-stack` 重复出现（{USAGE}）"));
            }
            let value: String = if let Some(v) = arg.strip_prefix("--work-stack=") {
                v.to_string()
            } else {
                let Some(v) = args.get(index + 1) else {
                    return Err(format!("`--work-stack` 缺少参数 local|global（{USAGE}）"));
                };
                index += 1;
                v.clone()
            };
            let v = match value.as_str() {
                "local" => false,
                "global" => true,
                other => {
                    return Err(format!(
                        "未知 --work-stack 值 `{other}`（期望 local|global，{USAGE}）"
                    ))
                }
            };
            work_stack_global = Some(v);
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
        // 默认 1 = 纯顺序路径（输出与升级前逐字节一致）。
        threads: threads.unwrap_or(1),
        log_secs,
        mipgap_percent,
        // 默认 auto：heuristics=true + 自适应（控制器阈值高于原 7 算例的
        // 节点数 ⇒ 原算例 auto ≡ on，见 README）。
        heuristic_mode: heuristic_mode.unwrap_or(HeuristicModeCli::Auto),
        branch_rule, // None = 类自适应默认（R3）
        scaling: scaling.unwrap_or(true), // 默认 on（R4）
        mitm: mitm.unwrap_or(true),
        work_stack_global: work_stack_global.unwrap_or(false),      // 默认 on（R5-retry）
    })
}

/// 由 [`Config`] 构造求解选项。
///
/// 数值取契约默认值（`Tolerances::default()` 与 `BbOptions::default()` 的
/// 文档值），仅 `time_limit_ms`（经 `N.saturating_mul(1000)` 饱和映射）与
/// `heuristics`（按 `--heuristic-mode`：off ⇒ false，on/auto ⇒ true）按 CLI
/// 覆盖；其余字段仍用契约默认值字面量（本包内维护）。分支规则默认
/// Pseudocost + 可靠性冷启动探测预算 8（2026-10 起；CLI 不设旗标，默认
/// 路径即启用——纯被动复现走 bench 探针 `--branch-rule pseudo
/// --pc-probes 0`）。
fn solve_options_from(config: &Config) -> SolveOptions {
    SolveOptions {
        tol: Tolerances::default(),
        bb: BbOptions {
            branch_rule: config.branch_rule.unwrap_or(BranchRule::MostFractional),
            branch_rule_auto: config.branch_rule.is_none(), // None ⇒ 类自适应（R3）
            node_rule: NodeRule::DepthFirst,
            node_limit: 2_000_000,
            time_limit_ms: config
                .time_limit_secs
                .map_or(0, |secs| secs.saturating_mul(1000)),
            root_cuts: true,
            heuristics: !matches!(config.heuristic_mode, HeuristicModeCli::Off),
            dive_max_fixes: 4096,
            prune_abs: 1e-9,
            prune_rel: 1e-9,
            pc_reliability_probes: 8,
            mitm_windows: if config.mitm { 8 } else { 0 },
            global_work_stack: config.work_stack_global,
        },
        use_presolve: true,
        use_scaling: config.scaling,
    }
}

/// 由 [`Config`] 构造搜索控制（`--threads` / `--mipgap` / `--log-secs` /
/// `--heuristic-mode auto` 直通 [`SearchControl`]）。
///
/// adaptive_heuristics 仅在 auto 模式开启（on/off 均为 false）；默认
/// （threads=1、无 mipgap、无 log、auto）与升级前行为逐字节一致。
fn search_control_from(config: &Config) -> SearchControl {
    SearchControl {
        degrade: milp::search_ctrl::DegradePolicy::default(),
        threads: config.threads,
        mipgap_percent: config.mipgap_percent,
        adaptive_heuristics: config.heuristic_mode == HeuristicModeCli::Auto,
        log_period_secs: config.log_secs,
    }
}

/// 进度回调钩子（仅 `--log-secs` 显式给出时注入）。
///
/// 每次回调经 [`format_progress_line`] 格式化为单行写入 stdout：
/// `stdout().lock()` 短暂持锁（回调内不获取求解器任何锁、立即返回），
/// 写失败（如管道断裂）静默忽略而非 panic。
fn progress_hooks() -> SearchHooks {
    SearchHooks {
        progress: Arc::new(|info: &ProgressInfo| {
            let mut stdout = std::io::stdout().lock();
            let _ = writeln!(stdout, "{}", format_progress_line(info));
        }),
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

/// gap 终止的报告输出与退出码。
///
/// 输出契约：依次 `status = gap_reached`、`nodes = N`、`gap = <%.2f>`、
/// `best = <数值>`（无 objective 行——gap 终止未证明最优，绝不冒充）。
/// gap 数值由库层触发时刻证据快照 + 统一 clamp（≥0）双保证为非负有限；
/// 本函数仍做防御复核（证据缺失/非有限/负 ⇒ 内部不一致，按错误退出 1），
/// 绝不编造或原样输出矛盾数值（如 inf/负百分比）。退出码 4。
fn gap_report_and_exit(controlled: &ControlledSolveReport) -> ExitCode {
    let Some(gap) = controlled.gap_percent.filter(|g| g.is_finite() && *g >= 0.0) else {
        eprintln!("error: gap 终止但缺少非负有限的 gap 证据（内部不一致）");
        return ExitCode::from(1);
    };
    let mut stdout = std::io::stdout().lock();
    // 自动降级与 mipgap 达标共用证据与行形状，仅 status token 区分原因。
    if controlled.auto_degraded {
        write_degrade_report(
            &mut stdout,
            controlled.report.stats.nodes,
            gap,
            controlled.report.objective,
        );
    } else {
        write_gap_report(
            &mut stdout,
            controlled.report.stats.nodes,
            gap,
            controlled.report.objective,
        );
    }
    ExitCode::from(4)
}

/// [`gap_report_and_exit`] 的可测试内核：写入任意 `Write`（写失败静默忽略）。
fn write_gap_report<W: Write>(out: &mut W, nodes: u64, gap_percent: f64, best: Option<f64>) {
    let _ = writeln!(out, "status = gap_reached");
    let _ = writeln!(out, "nodes = {nodes}");
    let _ = writeln!(out, "gap = {gap_percent:.2}%");
    if let Some(best) = best {
        let _ = writeln!(out, "best = {}", format_objective(best));
    }
}

/// 自动降级终止的报告输出（行形状同 [`write_gap_report`]，status token
/// 区分触发原因：`auto_degraded` = auto 模式下「继续证最优」被策略性放弃）。
fn write_degrade_report<W: Write>(out: &mut W, nodes: u64, gap_percent: f64, best: Option<f64>) {
    let _ = writeln!(out, "status = auto_degraded");
    let _ = writeln!(out, "nodes = {nodes}");
    let _ = writeln!(out, "gap = {gap_percent:.2}%");
    if let Some(best) = best {
        let _ = writeln!(out, "best = {}", format_objective(best));
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

    /// 构造 gap 终止的受控报告（write/gap 路径测试用）。
    fn controlled_gap(nodes: u64, objective: Option<f64>, gap: Option<f64>) -> ControlledSolveReport {
        ControlledSolveReport {
            report: report(SolveStatus::LimitReached, objective, None, nodes),
            gap_terminated: true,
            auto_degraded: false,
            gap_percent: gap,
            bound: None,
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

    /// Config 便捷构造（新字段取默认值；单测内集中补齐，避免每个用例重复）。
    fn config_with(time_limit_secs: Option<u64>) -> Config {
        Config {
            instance_path: "a.milp".to_string(),
            time_limit_secs,
            verbose: false,
            threads: 1,
            log_secs: None,
            mipgap_percent: None,
            heuristic_mode: HeuristicModeCli::Auto,
            branch_rule: None,
            scaling: true,
            mitm: true,
            work_stack_global: false,
        }
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

    // —— 新 flag 解析：默认值 ——

    /// 仅路径时新字段全为默认：threads=1、log_secs=None、mipgap=None、
    /// heuristic_mode=Auto（默认 = 升级前行为 + 自适应；原 7 算例 auto≡on）。
    #[test]
    fn parse_args_new_flags_defaults() {
        let config = parse_args(&args(&["a.milp"])).unwrap();
        assert_eq!(config.threads, 1, "默认 1 = 纯顺序路径");
        assert_eq!(config.log_secs, None, "默认不输出任何进度行");
        assert_eq!(config.mipgap_percent, None, "默认关闭 mipgap");
        assert_eq!(config.heuristic_mode, HeuristicModeCli::Auto);
    }

    // —— 新 flag 解析：合法全表 ——

    /// --threads / --log-secs / --mipgap / --heuristic-mode 的合法组合
    ///（含边界值、flag 先于路径、全组合一次给出）。
    #[test]
    fn parse_args_new_flags_legal() {
        // --threads：1 与 1024 边界、普通值、先于路径。
        assert_eq!(
            parse_args(&args(&["a.milp", "--threads", "1"])).unwrap().threads,
            1
        );
        assert_eq!(
            parse_args(&args(&["a.milp", "--threads", "1024"]))
                .unwrap()
                .threads,
            1024
        );
        assert_eq!(
            parse_args(&args(&["--threads", "4", "a.milp"]))
                .unwrap()
                .threads,
            4
        );
        // --log-secs：小数/科学记数（有限且 >0）。
        assert_eq!(
            parse_args(&args(&["a.milp", "--log-secs", "0.001"]))
                .unwrap()
                .log_secs,
            Some(0.001)
        );
        assert_eq!(
            parse_args(&args(&["a.milp", "--log-secs", "1e-3"]))
                .unwrap()
                .log_secs,
            Some(0.001)
        );
        assert_eq!(
            parse_args(&args(&["a.milp", "--log-secs", "5"]))
                .unwrap()
                .log_secs,
            Some(5.0)
        );
        // --mipgap：0（≡关闭）与正值。
        assert_eq!(
            parse_args(&args(&["a.milp", "--mipgap", "0"]))
                .unwrap()
                .mipgap_percent,
            Some(0.0)
        );
        assert_eq!(
            parse_args(&args(&["a.milp", "--mipgap", "3.5"]))
                .unwrap()
                .mipgap_percent,
            Some(3.5)
        );
        // --heuristic-mode：三档 token。
        for (token, mode) in [
            ("auto", HeuristicModeCli::Auto),
            ("on", HeuristicModeCli::On),
            ("off", HeuristicModeCli::Off),
        ] {
            assert_eq!(
                parse_args(&args(&["a.milp", "--heuristic-mode", token]))
                    .unwrap()
                    .heuristic_mode,
                mode
            );
        }
        // 全部新 flag + 既有 flag 一次组合（顺序任意、路径居中）。
        let config = parse_args(&args(&[
            "--threads",
            "4",
            "a.milp",
            "--time-limit-secs",
            "9",
            "--log-secs",
            "0.5",
            "--mipgap",
            "2.5",
            "--heuristic-mode",
            "on",
            "--verbose",
        ]))
        .expect("全组合必须解析成功");
        assert_eq!(config.instance_path, "a.milp");
        assert_eq!(config.threads, 4);
        assert_eq!(config.time_limit_secs, Some(9));
        assert_eq!(config.log_secs, Some(0.5));
        assert_eq!(config.mipgap_percent, Some(2.5));
        assert_eq!(config.heuristic_mode, HeuristicModeCli::On);
        assert!(config.verbose);
    }

    // —— 新 flag 解析：非法全表 ——

    /// 新 flag 各类非法一律 Err：threads 0/越界/非整数/负/缺失/重复；
    /// log-secs 0/负/非数/NaN/inf/缺失/重复；mipgap 负/NaN/inf/非数/缺失/
    /// 重复；heuristic-mode 未知 token（大小写敏感）/缺失/重复。
    #[test]
    fn parse_args_new_flags_invalid() {
        let bad: &[&[&str]] = &[
            // --threads。
            &["a.milp", "--threads"],                                // N 缺失
            &["a.milp", "--threads", "0"],                          // 下越界
            &["a.milp", "--threads", "1025"],                       // 上越界
            &["a.milp", "--threads", "99999999999999999999999"],    // 上越界（溢出量级）
            &["a.milp", "--threads", "abc"],                        // 非整数
            &["a.milp", "--threads", "3.5"],                        // 非整数
            &["a.milp", "--threads", "-4"],                         // 负数
            &["a.milp", "--threads", "2", "--threads", "3"],        // 重复
            // --log-secs。
            &["a.milp", "--log-secs"],                              // F 缺失
            &["a.milp", "--log-secs", "0"],                         // 非正
            &["a.milp", "--log-secs", "-1.5"],                      // 负
            &["a.milp", "--log-secs", "abc"],                       // 非数
            &["a.milp", "--log-secs", "NaN"],                       // NaN
            &["a.milp", "--log-secs", "inf"],                       // inf
            &["a.milp", "--log-secs", "1e999"],                     // 溢出为 inf
            &["a.milp", "--log-secs", "0.5", "--log-secs", "1"],    // 重复
            // --mipgap。
            &["a.milp", "--mipgap"],                                // P 缺失
            &["a.milp", "--mipgap", "-1"],                          // 负
            &["a.milp", "--mipgap", "abc"],                         // 非数
            &["a.milp", "--mipgap", "NaN"],                         // NaN
            &["a.milp", "--mipgap", "inf"],                         // inf
            &["a.milp", "--mipgap", "1e999"],                       // 溢出为 inf
            &["a.milp", "--mipgap", "1", "--mipgap", "2"],          // 重复
            // --heuristic-mode。
            &["a.milp", "--heuristic-mode"],                        // token 缺失
            &["a.milp", "--heuristic-mode", "yes"],                 // 未知 token
            &["a.milp", "--heuristic-mode", "AUTO"],                // 大小写敏感
            &["a.milp", "--heuristic-mode", "on", "--heuristic-mode", "off"], // 重复
        ];
        for case in bad {
            match parse_args(&args(case)) {
                Err(message) => assert!(!message.is_empty(), "参数 {case:?}"),
                Ok(config) => panic!("参数 {:?} 应报错，得到 {config:?}", case),
            }
        }
        // 数值 flag 缺参时不应把后续 flag 误当数值（--threads 后跟 --verbose）。
        assert!(parse_args(&args(&["a.milp", "--threads", "--verbose"])).is_err());
    }

    // —— heuristic-mode → BbOptions / SearchControl 映射 ——

    /// off ⇒ heuristics=false；on/auto ⇒ true；adaptive 仅 auto 开启；
    /// threads/mipgap/log 直通 SearchControl。
    #[test]
    fn heuristic_mode_maps_to_options_and_control() {
        // off：启发式关闭。
        let config = Config {
            heuristic_mode: HeuristicModeCli::Off,
            ..config_with(None)
        };
        assert!(!solve_options_from(&config).bb.heuristics);
        assert!(!search_control_from(&config).adaptive_heuristics);
        // on：启发式开、不自适应（升级前行为）。
        let config = Config {
            heuristic_mode: HeuristicModeCli::On,
            ..config_with(None)
        };
        assert!(solve_options_from(&config).bb.heuristics);
        assert!(!search_control_from(&config).adaptive_heuristics);
        // auto：启发式开 + 自适应（默认）。
        let config = Config {
            heuristic_mode: HeuristicModeCli::Auto,
            ..config_with(None)
        };
        assert!(solve_options_from(&config).bb.heuristics);
        assert!(search_control_from(&config).adaptive_heuristics);
        // 直通字段。
        let config = Config {
            threads: 4,
            log_secs: Some(0.25),
            mipgap_percent: Some(3.0),
            ..config_with(None)
        };
        let ctrl = search_control_from(&config);
        assert_eq!(ctrl.threads, 4);
        assert_eq!(ctrl.log_period_secs, Some(0.25));
        assert_eq!(ctrl.mipgap_percent, Some(3.0));
        // 其余 BbOptions 字段仍为契约默认值字面量。
        let opts = solve_options_from(&config);
        assert_eq!(opts.bb.node_limit, 2_000_000);
        assert!(opts.bb.root_cuts);
        assert_eq!(opts.bb.dive_max_fixes, 4096);
        assert_eq!(opts.bb.prune_abs, 1e-9);
        assert_eq!(opts.bb.prune_rel, 1e-9);
        // 默认分支规则：类自适应（R3）——字段基线 MostFractional + auto=true
        //（求解入口按结构类检测翻转；显式 --branch-rule 关 auto）。
        assert_eq!(opts.bb.branch_rule, BranchRule::MostFractional);
        assert!(opts.bb.branch_rule_auto);
        assert_eq!(opts.bb.pc_reliability_probes, 8);
        assert!(opts.use_presolve);
    }

    // —— time_limit_secs → time_limit_ms 映射 ——

    /// 映射必须 saturating：u64::MAX 秒 → u64::MAX 毫秒（不 panic、不回绕）。
    #[test]
    fn time_limit_mapping_saturates() {
        let config = config_with(None);
        assert_eq!(solve_options_from(&config).bb.time_limit_ms, 0);
        let config = config_with(Some(2));
        assert_eq!(solve_options_from(&config).bb.time_limit_ms, 2000);
        let config = config_with(Some(u64::MAX));
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

    // —— gap_terminated 报告输出契约 ——

    /// gap 终止输出逐行精确：status = gap_reached、nodes、gap（非负有限
    /// 百分比，{:.2}% 且绝无 -0.00/inf）、best（无 objective 行）；退出码 4；
    /// 证据缺失/非有限/负 ⇒ 防御性按错误退出 1（不编造数值）。
    #[test]
    fn gap_report_exact_lines_and_exit_code() {
        // 可测试内核：常规触发（gap≈2.38）。
        let mut out = Vec::new();
        write_gap_report(&mut out, 17, 2.3809523809523813, Some(21.0));
        let text = String::from_utf8(out).unwrap();
        assert_eq!(
            text,
            "status = gap_reached\nnodes = 17\ngap = 2.38%\nbest = 21.000000\n",
            "gap 终止输出必须逐行精确"
        );
        assert!(
            !text.contains("objective"),
            "gap 终止绝不输出 objective 行（未证明最优）"
        );
        // clamp 后的 0 gap：输出 0.00% 而非 -0.00%。
        let mut out = Vec::new();
        write_gap_report(&mut out, 5, 0.0, Some(40.0));
        let text = String::from_utf8(out).unwrap();
        assert_eq!(text, "status = gap_reached\nnodes = 5\ngap = 0.00%\nbest = 40.000000\n");
        // 大 gap 百分比仍为有限非负。
        let mut out = Vec::new();
        write_gap_report(&mut out, 2, 46.666_666_666_666_664, Some(10.0));
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "status = gap_reached\nnodes = 2\ngap = 46.67%\nbest = 10.000000\n"
        );
        // 自动降级：行形状同 gap 终止，仅 status token 区分触发原因；
        // 同样绝不输出 objective 行。
        let mut out = Vec::new();
        write_degrade_report(&mut out, 7, 2.3809523809523813, Some(21.0));
        let degrade_text = String::from_utf8(out).unwrap();
        assert_eq!(
            degrade_text,
            "status = auto_degraded\nnodes = 7\ngap = 2.38%\nbest = 21.000000\n"
        );
        assert!(
            !degrade_text.contains("objective ="),
            "降级终止绝不输出 objective 行（未证明最优）"
        );
        // 端到端（stdout 由测试框架捕获）：合法证据 ⇒ 退出码 4。
        assert_eq!(
            gap_report_and_exit(&controlled_gap(17, Some(21.0), Some(2.3809523809523813))),
            ExitCode::from(4),
            "gap_terminated ⇒ 退出码 4"
        );
        // 防御复核：证据缺失/NaN/inf/负 ⇒ 退出码 1，不输出矛盾数值。
        for bad_gap in [None, Some(f64::NAN), Some(f64::INFINITY), Some(-1e-11)] {
            assert_eq!(
                gap_report_and_exit(&controlled_gap(9, Some(1.0), bad_gap)),
                ExitCode::from(1),
                "非法 gap 证据 {bad_gap:?} 必须按内部不一致退出"
            );
        }
    }

    // —— 进度回调钩子 ——

    /// 回调可调用且不 panic（单行写 stdout；写失败静默；stdout 由测试框架
    /// 捕获）。行格式的逐字段口径由库层 search_ctrl::format_progress_line
    /// 的单元测试钉死。
    #[test]
    fn progress_hook_callable_without_panic() {
        let hooks = progress_hooks();
        let info = ProgressInfo {
            sense: Sense::Maximize,
            elapsed_secs: 0.0,
            nodes: 1,
            incumbent: Some(13.0),
            bound: Some(13.0),
        };
        (hooks.progress)(&info);
        let info = ProgressInfo {
            sense: Sense::Maximize,
            elapsed_secs: 1.25,
            nodes: 7,
            incumbent: None,
            bound: None,
        };
        (hooks.progress)(&info);
    }

    // —— 本包交付的算例文件可解析（WP1 已实现，纯解析层面即可验证）——

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

    /// 新算例 hard-knapsack-12 可被 load_instance 解析且校验通过、变量数 12。
    #[test]
    fn hard_knapsack_12_instance_parses() {
        let path = format!(
            "{}/instances/hard-knapsack-12.milp",
            env!("CARGO_MANIFEST_DIR")
        );
        let problem = load_instance(Path::new(&path))
            .unwrap_or_else(|error| panic!("算例 hard-knapsack-12 解析失败：{error}"));
        assert_eq!(problem.name, "hard-knapsack-12");
        assert!(problem.validate().is_ok());
        assert_eq!(problem.num_vars(), 12);
        assert_eq!(problem.constraints.len(), 1);
    }
}
