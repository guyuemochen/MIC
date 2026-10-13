//! profile_stats —— 库级基准测量探针（bench 计时矩阵的测量工具，M 包重建版）。
//!
//! 用法：`profile_stats <算例路径> [--threads N] [--time-limit-secs S]
//! [--branch-rule mostfrac|pseudo] [--pc-probes N] [--root-cuts on|off]`。
//!
//! # 定位与调用契约
//! - 纯测量工具：只经库 pub API（[`milp::solver::solve_controlled`]）求解，
//!   **不修改求解语义、不触碰 `milp/src` 任何文件**；探针存在与否不影响求解器。
//! - 调用方 `bench/run_timings.py`：`profile_stats.exe <instance.milp>
//!   --threads N --time-limit-secs S`，stdout **末行**为单行 JSON
//!   （本探针 stdout 仅此一行）。
//! - `wall_secs` 为 solve_controlled 周期内计时（不含算例加载与 root_bound
//!   旁路，与 bench 既有口径一致）；`stats_time_ms` = `BbStats::time_ms`。
//! - 退出码：0 = 完成（**含 limit_reached/auto_degraded 等诚实终态**——
//!   bench 侧 `subprocess.run(check=True)` 要求非错误终态也退出 0）；
//!   1 = 用法错误 / 算例加载失败 / 求解返回 Err。
//!
//! # JSON 字段集（对照 bench/analyze_baseline.py 硬消费点逐字段给出）
//! - `instance`（file_stem）、`threads`、`nodes`、`objective`、
//!   `lp_iterations`、`wall_secs`、`status`、`cuts`、`stats_time_ms`、
//!   `bound`、`gap_percent`（注意：字段名是 `gap_percent`；输出侧键
//!   `gap_percent_at_stop` 是 analyze_baseline.py 的重命名）；
//! - `presolve_reductions`（诊断用，旧探针遗留字段，非 analyze 消费点）；
//! - `root_bound`（新增，见下）。
//! - 浮点一律舍入到 6 位小数后以 Rust `{:?}`（Debug）输出（值级兼容现存
//!   baseline.json 的浮点口径，如 `165.0`/`112.668615`）；非有限值
//!   （inf/NaN，防御路径）落 `null`。
//!
//! # 求解选项：复刻 CLI 默认 + 构造纪律
//! - 选项 = CLI 默认（MostFractional/DepthFirst/node_limit=2M/root_cuts/
//!   heuristics=auto 即 `adaptive_heuristics=true`/dive 4096/presolve on）
//!   叠加旗标差异：`--time-limit-secs`（saturating_mul(1000) 映射毫秒）、
//!   `--branch-rule`、`--root-cuts off`。
//! - **构造纪律**：分支/割/预算选项一律 `..BbOptions::default()` /
//!   `..SolveOptions::default()` / `..SearchControl::default()` 展开，只覆盖
//!   差异字段——包 B 将给 `BbOptions` 新增 `pc_reliability_probes` 字段，
//!   全字段字面量构造会在 B 合入时把本文件编译断裂；展开式则只需 B 在
//!   旗标映射处补一行。
//! - 旗标间**无隐含联动**：纯被动 Pseudocost = `--branch-rule pseudo
//!   --pc-probes 0`（即现状 `BranchRule::Pseudocost`）；默认可靠性 =
//!   `--branch-rule pseudo --pc-probes 8`。B 合入前 `--pc-probes` 非 0
//!   时本探针**直接报错退出**（`pc_reliability_probes` 字段尚不存在，
//!   接受该值即静默无效），B 合入后由 B 接线到新字段。
//!
//! # root_bound 旁路（正确性论证）
//! - `root_bound` = **根割后的根 LP 界**（原始 sense + obj_offset 口径）：
//!   与求解器同款的 pub API 复现根割循环——`presolve::presolve`（与
//!   solve_controlled 默认编排一致）→ `Simplex::new` + `solve` →
//!   `cuts::cut_loop_recorded(CutOptions::default())` → `resolve` 读
//!   objective。该循环与 `branch_bound` 根阶段逐字节同路径（`cut_loop`
//!   即 `cut_loop_recorded` 丢弃割清单后的委托），无私有状态访问。
//! - 旁路在主求解（计时窗口）**之后**、以全新 Simplex 实例运行，不影响
//!   `wall_secs` 与主求解统计；确定性：同输入同结果。
//! - 无合法界时为 `null`：根 LP 非 Optimal、割循环终态非 Optimal（如
//!   setcover-100x150 的根割 IterationLimit 传播）、或 resolve 非 Optimal。
//! - 用途：供包 A（根割引擎）测量根界推进（改造前后 root_bound 对照）。

use std::io::Write;
use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use milp::branch_bound::{BbOptions, BranchRule, SearchControl};
use milp::cuts::{cut_loop_recorded, CutOptions};
use milp::instance::load_instance;
use milp::model::{Problem, SolveStatus, Tolerances};
use milp::presolve;
use milp::simplex::{Bounds, LpStatus, Simplex};
use milp::solver::{solve_controlled, ControlledSolveReport, SolveOptions};

/// 用法说明（打印到 stderr）。
const USAGE: &str = "用法：profile_stats <算例路径> [--threads N] [--time-limit-secs S] \
[--branch-rule mostfrac|pseudo] [--pc-probes N] [--root-cuts on|off] \
[--work-stack auto|local|global]";

/// `--threads` 允许的最大线程数（与 CLI main.rs 同口径）。
const MAX_THREADS: usize = 1024;

/// 命令行配置（默认值 = CLI 默认的测量复刻：threads=1、无时限、mostfrac、
/// pc-probes=8（复刻库默认）、root-cuts on、work-stack auto（W1 按类选栈，
/// 复刻 CLI 多线程默认））。
#[derive(Debug)]
struct Config {
    /// 算例文件路径。
    instance_path: String,
    /// 线程数（1..=1024）；默认 1。
    threads: usize,
    /// 时间上限秒（None = 不限）；经 saturating_mul(1000) 映射毫秒。
    time_limit_secs: Option<u64>,
    /// 分支规则旗标；默认 None = 类自适应（R3，复刻 CLI 默认；求解入口按
    /// 结构类检测解析）。显式 mostfrac/pseudo 覆盖优先。
    branch_rule: Option<BranchRule>,
    /// Pseudocost 可靠性探针数；默认 8（复刻 BbOptions::pc_reliability_probes
    /// 库默认；0 = 纯被动）。R3 接线到 BbOptions（此前误留包 B 未合入守卫）。
    pc_probes: usize,
    /// 根割开关；默认 on（CLI 默认 root_cuts=true）。
    root_cuts: bool,
    /// 并行工作栈模式（`--work-stack auto|local|global`，R6 接线 + W1 按类
    /// 选栈；仅 threads ≥ 2 生效）；默认 auto（复刻 CLI 多线程默认：类检测
    /// 命中多维背包类 ⇒ global，否则 local）。global = 全局 DFS 栈，压并行
    /// 节点膨胀（P1 取证 §5 M1 提案的接线点）。
    work_stack: WorkStack,
}

/// 并行工作栈模式（与 CLI main.rs 同口径的三态）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WorkStack {
    /// 类自适应（求解入口按结构类检测解析）。
    Auto,
    /// 本地亲和 + 偷取（显式覆盖）。
    Local,
    /// 全局 DFS 栈（显式覆盖）。
    Global,
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
    // wall_secs = solve_controlled 周期内计时（不含算例加载与 root_bound 旁路）。
    let start = Instant::now();
    let controlled = match solve_controlled(&problem, &options, &ctrl, None) {
        Ok(controlled) => controlled,
        Err(error) => {
            eprintln!("error: {error}");
            return ExitCode::from(1);
        }
    };
    let wall_secs = start.elapsed().as_secs_f64();
    // root_bound 旁路：计时窗口之外，全新 Simplex，不影响主求解统计。
    let root_bound = root_bound_after_root_cuts(&problem);
    let report = &controlled.report;
    let instance = Path::new(&config.instance_path)
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| config.instance_path.clone());
    let mut stdout = std::io::stdout().lock();
    let _ = writeln!(
        stdout,
        "{{\"instance\": {}, \"threads\": {}, \"wall_secs\": {}, \
\"stats_time_ms\": {}, \"status\": \"{}\", \"objective\": {}, \
\"nodes\": {}, \"lp_iterations\": {}, \"cuts\": {}, \
\"presolve_reductions\": {}, \"bound\": {}, \"gap_percent\": {}, \
\"root_bound\": {}}}",
        json_string(&instance),
        config.threads,
        json_f64(Some(wall_secs)),
        report.stats.time_ms,
        status_token(&controlled),
        json_f64(report.objective),
        report.stats.nodes,
        report.stats.lp_iterations,
        report.stats.cuts,
        report.presolve_reductions,
        json_f64(controlled.bound),
        json_f64(controlled.gap_percent),
        json_f64(root_bound),
    );
    // 诚实终态（limit_reached/auto_degraded/…）也是完成的测量 ⇒ 退出 0。
    ExitCode::SUCCESS
}

/// 解析命令行参数（风格与 CLI main.rs 一致：缺失/非法/重复/未知一律 Err；
/// 用户输入路径不 unwrap/expect/panic）。
fn parse_args(args: &[String]) -> Result<Config, String> {
    if args.is_empty() {
        eprintln!("{USAGE}");
        return Err("缺少算例路径".to_string());
    }
    let mut instance_path: Option<String> = None;
    let mut threads: Option<usize> = None;
    let mut time_limit_secs: Option<u64> = None;
    let mut branch_rule: Option<BranchRule> = None;
    let mut pc_probes: Option<usize> = None;
    let mut root_cuts: Option<bool> = None;
    let mut work_stack: Option<WorkStack> = None;
    let mut index = 0usize;
    while index < args.len() {
        let arg = args[index].as_str();
        if arg == "--threads" {
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
            index += 1;
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
            index += 1;
        } else if arg == "--branch-rule" {
            if branch_rule.is_some() {
                return Err(format!("`--branch-rule` 重复出现（{USAGE}）"));
            }
            let Some(value) = args.get(index + 1) else {
                return Err(format!("`--branch-rule` 缺少参数 mostfrac|pseudo（{USAGE}）"));
            };
            let rule = match value.as_str() {
                "mostfrac" => BranchRule::MostFractional,
                "pseudo" => BranchRule::Pseudocost,
                unknown => {
                    return Err(format!(
                        "`--branch-rule` 的参数 `{unknown}` 未知（期望 mostfrac|pseudo）（{USAGE}）"
                    ));
                }
            };
            branch_rule = Some(rule);
            index += 1;
        } else if arg == "--pc-probes" {
            if pc_probes.is_some() {
                return Err(format!("`--pc-probes` 重复出现（{USAGE}）"));
            }
            let Some(value) = args.get(index + 1) else {
                return Err(format!("`--pc-probes` 缺少数值参数 N（{USAGE}）"));
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
                        "`--pc-probes` 的参数 `{value}` {reason}（期望 ≥0 的整数）"
                    ));
                }
            };
            pc_probes = Some(n);
            index += 1;
        } else if arg == "--root-cuts" {
            if root_cuts.is_some() {
                return Err(format!("`--root-cuts` 重复出现（{USAGE}）"));
            }
            let Some(value) = args.get(index + 1) else {
                return Err(format!("`--root-cuts` 缺少参数 on|off（{USAGE}）"));
            };
            let on = match value.as_str() {
                "on" => true,
                "off" => false,
                unknown => {
                    return Err(format!(
                        "`--root-cuts` 的参数 `{unknown}` 未知（期望 on|off）（{USAGE}）"
                    ));
                }
            };
            root_cuts = Some(on);
            index += 1;
        } else if arg == "--work-stack" {
            if work_stack.is_some() {
                return Err(format!("`--work-stack` 重复出现（{USAGE}）"));
            }
            let Some(value) = args.get(index + 1) else {
                return Err(format!(
                    "`--work-stack` 缺少参数 auto|local|global（{USAGE}）"
                ));
            };
            let mode = match value.as_str() {
                "auto" => WorkStack::Auto,
                "local" => WorkStack::Local,
                "global" => WorkStack::Global,
                unknown => {
                    return Err(format!(
                        "`--work-stack` 的参数 `{unknown}` 未知（期望 auto|local|global）（{USAGE}）"
                    ));
                }
            };
            work_stack = Some(mode);
            index += 1;
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
        threads: threads.unwrap_or(1),
        time_limit_secs,
        branch_rule, // None = 类自适应（R3，复刻 CLI 默认）
        pc_probes: pc_probes.unwrap_or(8),
        root_cuts: root_cuts.unwrap_or(true),
        // 默认 auto（W1，复刻 CLI 多线程默认：类检测命中 ⇒ global，否则
        // local）；显式 local/global 供 A/B 消融。
        work_stack: work_stack.unwrap_or(WorkStack::Auto),
    })
}

/// 由 [`Config`] 构造求解选项：CLI 默认 + 旗标差异。
///
/// 构造纪律：`..BbOptions::default()` / `..SolveOptions::default()` 展开，
/// 只覆盖差异字段（branch_rule / root_cuts / time_limit_ms）——B 给
/// `BbOptions` 新增 `pc_reliability_probes` 时本函数不编译断裂。
/// CLI 默认复刻：node_limit=2M / DepthFirst / heuristics=true（auto 模式，
/// BbOptions 默认即 true）/ dive 4096 / prune 1e-9 / presolve on。
fn solve_options_from(config: &Config) -> SolveOptions {
    SolveOptions {
        bb: BbOptions {
            branch_rule: config.branch_rule.unwrap_or(BranchRule::MostFractional),
            branch_rule_auto: config.branch_rule.is_none(), // None ⇒ 类自适应（R3）
            pc_reliability_probes: config.pc_probes,
            root_cuts: config.root_cuts,
            // W1 按类选栈：auto ⇒ 求解入口按结构类检测解析（复刻 CLI 默认）；
            // 显式 local/global ⇒ 覆盖优先。
            global_work_stack: config.work_stack == WorkStack::Global,
            work_stack_auto: config.work_stack == WorkStack::Auto,
            time_limit_ms: config
                .time_limit_secs
                .map_or(0, |secs| secs.saturating_mul(1000)),
            ..BbOptions::default()
        },
        ..SolveOptions::default()
    }
}

/// 由 [`Config`] 构造搜索控制：threads 直通；`adaptive_heuristics = true`
/// 复刻 CLI 默认 `--heuristic-mode auto`（bench 基线口径）。
///
/// 构造纪律：`..SearchControl::default()` 展开，只覆盖差异字段。
fn search_control_from(config: &Config) -> SearchControl {
    SearchControl {
        threads: config.threads,
        adaptive_heuristics: true,
        ..SearchControl::default()
    }
}

/// 根割后的根 LP 界（旁路；正确性论证见模块文档）。
///
/// 与 solve_controlled 默认编排同款 pub API 复现根阶段：presolve →
/// Simplex::new + solve → `cut_loop_recorded(CutOptions::default())` →
/// resolve 读 objective（原始 sense + obj_offset 口径，Simplex 内部已含
/// offset；presolve 把固定列目标常数折叠进 offset，口径保持）。
/// 任一阶段非 Optimal（含割循环 IterationLimit 传播）⇒ None ⇒ JSON null。
fn root_bound_after_root_cuts(problem: &Problem) -> Option<f64> {
    let tol = Tolerances::default();
    let (work, _map) = presolve::presolve(problem, &tol).ok()?;
    let int_vars = work.integer_vars();
    let mut simplex = Simplex::new(&work, tol).ok()?;
    let bounds = Bounds::from_problem(&work);
    let root = simplex.solve(&bounds);
    if root.status != LpStatus::Optimal {
        return None;
    }
    let (_applied, cut_status, _recorded) =
        cut_loop_recorded(&mut simplex, &bounds, &work, &int_vars, CutOptions::default(), &tol)
            .ok()?;
    if cut_status != LpStatus::Optimal {
        return None;
    }
    let after = simplex.resolve(&bounds);
    (after.status == LpStatus::Optimal).then_some(after.objective)
}

/// 终态 token：与 CLI 输出契约同款区分（auto_degraded / gap_reached），
/// 其余按 [`SolveStatus`]（optimal / infeasible / unbounded / limit_reached）。
fn status_token(controlled: &ControlledSolveReport) -> &'static str {
    if controlled.gap_terminated && controlled.auto_degraded {
        "auto_degraded"
    } else if controlled.gap_terminated {
        "gap_reached"
    } else {
        match controlled.report.status {
            SolveStatus::Optimal => "optimal",
            SolveStatus::Infeasible => "infeasible",
            SolveStatus::Unbounded => "unbounded",
            SolveStatus::LimitReached => "limit_reached",
        }
    }
}

/// JSON 数值（f64）：与旧探针同款口径——先舍入到 6 位小数（值级兼容现存
/// baseline.json 的浮点字段，如 112.668615/165.0），再 Rust `{:?}` 输出
/// （恒含小数点或指数，合法 JSON 数值）；None、非有限（inf/NaN，防御路径）
/// 或超出舍入安全量级（|v| ≥ 1e9，乘 1e6 丢有效位）⇒ 原值/`null` 如实处理。
fn json_f64(value: Option<f64>) -> String {
    match value {
        Some(v) if !v.is_finite() => "null".to_string(),
        Some(v) => {
            let rounded = if v.abs() < 1e9 {
                (v * 1e6).round() / 1e6
            } else {
                v
            };
            format!("{rounded:?}")
        }
        None => "null".to_string(),
    }
}

/// JSON 字符串：转义 `"`、`\` 与控制字符（算例名/status 常规为 ASCII，
/// 转义仅防御外部路径名）。
fn json_string(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for c in value.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}
