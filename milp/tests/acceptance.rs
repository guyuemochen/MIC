//! 验收测试（WP5 契约）。
//!
//! - 5 个数值型算例走 Optimal 主路径：终态、目标值、维数、整数性、
//!   约束满足性（容差 1e-6×(1+|rhs|+‖row‖₁)）、目标按 values 重算一致；
//! - 2 个状态型算例：Infeasible / Unbounded（均无 objective/values）；
//! - CLI 层：objective 末行数值解析（容差 1e-6，不依赖具体小数位数）、
//!   不可行/无界退出码 2/3 且无 objective 行。
//!
//! 算例路径拼接用 `concat!(env!("CARGO_MANIFEST_DIR"), "/instances/", name, ".milp")`
//! 的等价形式（name 非字面量处用 format!）。
//!
//! # 受控搜索参数包扩展（追加测试）
//! - CLI 默认输出逐字节基线（原 7 算例；threads=1 与默认 auto 红线端到端
//!   钉死）；
//! - 库级/CLI 级 threads=4 并发一致性（含新算例 hard-knapsack-12：200+ 节点
//!   真实压满节点池；CLI 级重复 ≥10 次不依赖计时）；
//! - `--log-secs` 进度行格式（首行契约钉死 + 全行正则 + 空栈矛盾行防护）；
//! - `--mipgap`（0≡关闭、大阈值 gap_reached 择一断言、库级确定性用例）；
//! - 自适应启发式 auto 于原算例与新算例仍精确最优；`--heuristic-mode on`
//!   在 200+ 节点算例上钉死 on 路径基线。
//!
//! （时间/节点上限触发的 LimitReached 仍不做 CLI 级验收——依赖计时不稳定；
//! 由 WP4 单元测试以极小 node_limit 覆盖。gap 终止（`--mipgap`）不依赖计时，
//! 在 acceptance_cli_mipgap 做 CLI 级验收。）
//!
//! # 重负载验收（L4 门禁拆分）
//! setcover-100x150 依赖的两项验收（`acceptance_setcover_100x150_unlocked_optimal`
//! / `acceptance_default_path_deterministic_two_runs`）debug 构建下合计
//! ~130s（release 亚秒~秒级），已标 `#[ignore]` 移出默认快路径：默认
//! `cargo test` 只跑轻验收；门禁以 `cargo test --test acceptance --
//! --ignored` 显式补跑（脚本门禁另有 release 全量验收兜底），覆盖不丢。

use std::path::Path;
use std::process::Command;

use milp::branch_bound::{BbOptions, SearchControl};
use milp::instance::load_instance;
use milp::model::{Problem, RelOp, Sense, SolveStatus, VarType};
use milp::solver::{solve, solve_controlled, solve_default, SolveOptions};

/// 数值型验收算例与已知最优目标值（mustMatch，须精确达到，容差 1e-6）。
const CASES: &[(&str, f64)] = &[
    ("lp-prod-3", 13.0),
    ("knapsack-10", 17.0),
    ("assign-4x4", 10.0),
    ("milp-wolsey-2", 40.0),
    ("milp-mixed-5", 29.0),
];

/// 算例文件路径：`instances/<name>.milp`（相对 crate 根）。
fn instance_path(name: &str) -> String {
    format!("{}/instances/{name}.milp", env!("CARGO_MANIFEST_DIR"))
}

/// 每个 CASE：load_instance → solve_default → assert status==Optimal
/// → |objective - expected| <= 1e-6
/// → values 维数=num_vars、整数变量 |x-round|<=1e-5、全部约束满足
///   （容差 1e-6×(1+|rhs|+行L1范数)）、目标按 values 重算一致（<=1e-6）。
#[test]
fn acceptance_all_instances_optimal() {
    for &(name, expected) in CASES {
        let path = instance_path(name);
        let problem = load_instance(Path::new(&path))
            .unwrap_or_else(|error| panic!("算例 {name} 加载失败：{error}"));
        let report =
            solve_default(&problem).unwrap_or_else(|error| panic!("算例 {name} 求解失败：{error}"));

        // 终态必须为 Optimal（诚实性：不得把截断结果冒充最优）。
        assert_eq!(report.status, SolveStatus::Optimal, "算例 {name} 终态");

        // 目标值精确达到（容差 1e-6）。
        let objective = report
            .objective
            .unwrap_or_else(|| panic!("算例 {name} 报告 Optimal 但无 objective"));
        assert!(
            (objective - expected).abs() <= 1e-6,
            "算例 {name}：objective = {objective}，期望 {expected}"
        );

        // values 维数 = num_vars。
        let values = report
            .values
            .unwrap_or_else(|| panic!("算例 {name} 报告 Optimal 但无 values"));
        assert_eq!(values.len(), problem.num_vars(), "算例 {name} values 维数");

        // 整数变量 |x - round(x)| <= 1e-5。
        for &idx in &problem.integer_vars() {
            assert!(
                (values[idx] - values[idx].round()).abs() <= 1e-5,
                "算例 {name} 变量 {}（{}）整数性超差",
                problem.vars[idx].name,
                values[idx]
            );
        }

        // 全部约束满足：violation <= 1e-6 × (1 + |rhs| + ‖row‖₁)。
        for (row, constraint) in problem.constraints.iter().enumerate() {
            let lhs = constraint.expr.eval(&values);
            let row_l1: f64 = constraint.expr.terms.iter().map(|&(_, c)| c.abs()).sum();
            let tol = 1e-6 * (1.0 + constraint.rhs.abs() + row_l1);
            let violation = match constraint.op {
                RelOp::Le => lhs - constraint.rhs,
                RelOp::Ge => constraint.rhs - lhs,
                RelOp::Eq => (lhs - constraint.rhs).abs(),
            };
            assert!(
                violation <= tol,
                "算例 {name} 约束 {row}（{}）违反量 {violation} 超容差 {tol}",
                constraint.name.as_deref().unwrap_or("<未命名>")
            );
        }

        // 目标按 values 重算一致（含 obj_offset）。
        let recomputed = problem.obj.eval(&values) + problem.obj_offset;
        assert!(
            (recomputed - objective).abs() <= 1e-6,
            "算例 {name} 目标重算 {recomputed} 与报告值 {objective} 不一致"
        );
    }
}

/// infeasible-tiny ⇒ Infeasible；unbounded-tiny ⇒ Unbounded
/// （均无 objective/values）。
#[test]
fn acceptance_status_instances() {
    let path = instance_path("infeasible-tiny");
    let problem = load_instance(Path::new(&path))
        .unwrap_or_else(|error| panic!("算例 infeasible-tiny 加载失败：{error}"));
    let report = solve_default(&problem)
        .unwrap_or_else(|error| panic!("算例 infeasible-tiny 求解失败：{error}"));
    assert_eq!(report.status, SolveStatus::Infeasible);
    assert!(report.objective.is_none(), "不可行不得携带 objective");
    assert!(report.values.is_none(), "不可行不得携带 values");

    let path = instance_path("unbounded-tiny");
    let problem = load_instance(Path::new(&path))
        .unwrap_or_else(|error| panic!("算例 unbounded-tiny 加载失败：{error}"));
    let report = solve_default(&problem)
        .unwrap_or_else(|error| panic!("算例 unbounded-tiny 求解失败：{error}"));
    assert_eq!(report.status, SolveStatus::Unbounded);
    assert!(report.objective.is_none(), "无界不得携带 objective");
    assert!(report.values.is_none(), "无界不得携带 values");
}

/// CLI 求解 knapsack-10：退出码 0，stdout 最后一行 strip_prefix
/// "objective = " 后 parse::<f64>() 与 17.0 相差 <= 1e-6
/// （数值解析比较，不依赖具体小数位数）。
#[test]
fn acceptance_cli_objective_line() {
    let output = Command::new(env!("CARGO_BIN_EXE_milp"))
        .arg(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/instances/knapsack-10.milp"
        ))
        .output()
        .expect("启动 milp CLI 失败");
    assert_eq!(
        output.status.code(),
        Some(0),
        "knapsack-10 应以退出码 0 报最优，stderr：{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    let last_line = stdout
        .trim_end_matches(['\r', '\n'])
        .lines()
        .next_back()
        .unwrap_or_default();
    let value_text = last_line.strip_prefix("objective = ").unwrap_or_else(|| {
        panic!(
            "stdout 最后一行应为 `objective = <数值>`，得到 {last_line:?}（完整 stdout：{stdout}）"
        )
    });
    let value: f64 = value_text
        .parse()
        .unwrap_or_else(|_| panic!("objective 行无法数值解析：{value_text:?}"));
    assert!(
        (value - 17.0).abs() <= 1e-6,
        "CLI objective = {value}，期望 17.0"
    );
}

/// CLI 状态型算例退出码：
/// - infeasible-tiny ⇒ 退出码 2、stdout 含 "status = infeasible"、无 "objective = " 行；
/// - unbounded-tiny ⇒ 退出码 3、stdout 含 "status = unbounded"、无 "objective = " 行。
///
/// （LimitReached 退出码 4 不做 CLI 级验收——依赖计时不稳定；由 WP4 单元
/// 测试以极小 node_limit 覆盖。）
#[test]
fn acceptance_cli_exit_codes() {
    let output = Command::new(env!("CARGO_BIN_EXE_milp"))
        .arg(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/instances/infeasible-tiny.milp"
        ))
        .output()
        .expect("启动 milp CLI 失败");
    assert_eq!(
        output.status.code(),
        Some(2),
        "infeasible-tiny 应以退出码 2 报不可行，stdout：{}",
        String::from_utf8_lossy(&output.stdout)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("status = infeasible"),
        "stdout 应含状态行：{stdout}"
    );
    assert!(
        !stdout.contains("objective = "),
        "不可行不得输出 objective 行：{stdout}"
    );

    let output = Command::new(env!("CARGO_BIN_EXE_milp"))
        .arg(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/instances/unbounded-tiny.milp"
        ))
        .output()
        .expect("启动 milp CLI 失败");
    assert_eq!(
        output.status.code(),
        Some(3),
        "unbounded-tiny 应以退出码 3 报无界，stdout：{}",
        String::from_utf8_lossy(&output.stdout)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("status = unbounded"),
        "stdout 应含状态行：{stdout}"
    );
    assert!(
        !stdout.contains("objective = "),
        "无界不得输出 objective 行：{stdout}"
    );
}

// ===========================================================================
// 受控搜索参数包扩展（既有 4 测试不动，以下全部追加；全部不依赖计时——
// 首行进度由『根阶段结束必发』契约保证，并发断言只看 objective/status）
// ===========================================================================

/// 运行 milp CLI（args 首元素为算例路径）。
fn cli(args: &[&str]) -> (Option<i32>, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_milp"))
        .args(args)
        .output()
        .expect("启动 milp CLI 失败");
    (output.status.code(), String::from_utf8_lossy(&output.stdout).into_owned())
}

/// stdout 最后一行 strip_prefix "objective = " 后数值解析。
fn last_objective(stdout: &str) -> f64 {
    let last_line = stdout
        .trim_end_matches(['\r', '\n'])
        .lines()
        .next_back()
        .unwrap_or_default();
    let value_text = last_line.strip_prefix("objective = ").unwrap_or_else(|| {
        panic!(
            "stdout 最后一行应为 `objective = <数值>`，得到 {last_line:?}（完整 stdout：{stdout}）"
        )
    });
    value_text
        .parse()
        .unwrap_or_else(|_| panic!("objective 行无法数值解析：{value_text:?}"))
}

/// 数值字段校验（无正则 crate，手写等价判定）：
/// `none` 或 `-?\d+(\.\d+)?([eE][+-]?\d+)?`。
fn number_token_ok(token: &str) -> bool {
    if token == "none" {
        return true;
    }
    let body = token.strip_prefix('-').unwrap_or(token);
    let mut split_exp = body.splitn(2, ['e', 'E']);
    let mantissa = split_exp.next().unwrap_or("");
    let exponent = split_exp.next();
    let mut split_dot = mantissa.splitn(2, '.');
    let int_part = split_dot.next().unwrap_or("");
    let frac_part = split_dot.next();
    let digits = |s: &str| !s.is_empty() && s.chars().all(|c| c.is_ascii_digit());
    if !digits(int_part) {
        return false;
    }
    if let Some(frac) = frac_part {
        if !digits(frac) {
            return false;
        }
    }
    if let Some(exp) = exponent {
        let exp_digits = exp.strip_prefix(['+', '-']).unwrap_or(exp);
        if !digits(exp_digits) {
            return false;
        }
    }
    true
}

/// 进度行整行校验，等价于
/// `^progress t=\d+\.\d+ nodes=\d+ incumbent=(none|-?\d+(\.\d+)?([eE][+-]?\d+)?) bound=(none|-?\d+(\.\d+)?([eE][+-]?\d+)?) gap=(inf|\d+\.\d+%)$`。
fn progress_line_ok(line: &str) -> bool {
    let digits = |s: &str| !s.is_empty() && s.chars().all(|c| c.is_ascii_digit());
    // progress t=<\d+\.\d+>
    let Some(rest) = line.strip_prefix("progress t=") else {
        return false;
    };
    let Some((t, rest)) = rest.split_once(' ') else {
        return false;
    };
    let mut t_parts = t.splitn(2, '.');
    if !digits(t_parts.next().unwrap_or("")) || !digits(t_parts.next().unwrap_or("")) {
        return false;
    }
    // nodes=<\d+>
    let Some(rest) = rest.strip_prefix("nodes=") else {
        return false;
    };
    let Some((nodes, rest)) = rest.split_once(' ') else {
        return false;
    };
    if !digits(nodes) {
        return false;
    }
    // incumbent=<number>
    let Some(rest) = rest.strip_prefix("incumbent=") else {
        return false;
    };
    let Some((incumbent, rest)) = rest.split_once(' ') else {
        return false;
    };
    if !number_token_ok(incumbent) {
        return false;
    }
    // bound=<number>
    let Some(rest) = rest.strip_prefix("bound=") else {
        return false;
    };
    let Some((bound, rest)) = rest.split_once(' ') else {
        return false;
    };
    if !number_token_ok(bound) {
        return false;
    }
    // gap=(inf|\d+\.\d+%)
    let Some(gap) = rest.strip_prefix("gap=") else {
        return false;
    };
    if gap == "inf" {
        return true;
    }
    let Some(percent) = gap.strip_suffix('%') else {
        return false;
    };
    let mut gap_parts = percent.splitn(2, '.');
    digits(gap_parts.next().unwrap_or("")) && digits(gap_parts.next().unwrap_or(""))
}

/// 对一份 CLI stdout 的全部 progress 行做统一校验：
/// ≥1 行、每行匹配 [`progress_line_ok`]、首 progress 行先于 status 行。
fn assert_progress_lines_wellformed(stdout: &str) -> Vec<&str> {
    let lines: Vec<&str> = stdout.lines().collect();
    let progress: Vec<&str> = lines
        .iter()
        .copied()
        .filter(|l| l.starts_with("progress "))
        .collect();
    assert!(!progress.is_empty(), "必须至少一行进度：{stdout}");
    for line in &progress {
        assert!(progress_line_ok(line), "进度行格式非法：{line:?}（stdout：{stdout}）");
    }
    let first_progress = lines
        .iter()
        .position(|l| l.starts_with("progress "))
        .expect("上方已保证存在 progress 行");
    let status_line = lines
        .iter()
        .position(|l| l.starts_with("status "))
        .unwrap_or(usize::MAX);
    assert!(
        first_progress < status_line,
        "首 progress 行必须先于 status 行：{stdout}"
    );
    progress
}

/// 首 progress 行去掉 `progress t=<数值>` 前缀后的剩余字段（逐字段精确断言用；
/// t 为墙钟数值不钉）。
fn first_progress_tail(stdout: &str) -> &str {
    let first = stdout
        .lines()
        .find(|l| l.starts_with("progress "))
        .unwrap_or_else(|| panic!("无 progress 行：{stdout}"));
    let idx = first
        .find(" nodes=")
        .unwrap_or_else(|| panic!("progress 行缺 nodes 字段：{first:?}"));
    assert!(
        first.starts_with("progress t=") && {
            let t = &first["progress t=".len()..idx];
            !t.is_empty()
                && t.chars().all(|c| c.is_ascii_digit() || c == '.')
                && t.contains('.')
        },
        "t 字段应为 \\d+\\.\\d+：{first:?}"
    );
    &first[idx..]
}

/// dive 问题（与 solver/branch_bound 单元测试同构；关 presolve 与根割保证
/// gap 触发确定性）：max 10·x1 + x2 + 0.5 s.t. 20·x1 + 3·x2 ≤ 43，整数
/// [0,10]²。根 LP 22.0 分数，根启发式 dive 得 (2,1)=21.5，整数最优 21.5。
/// offset=0.5 有意非整数化（整数目标 cutoff 门控判否）：offset=0.0 时根
/// bound 21.5 ≤ 21+1−ε 会被 cutoff 在 (d) 直接闭合，gap 触发路径消失。
fn dive_gap_problem() -> Problem {
    let mut p = Problem::new("dive-gap", Sense::Maximize);
    p.add_var("x1", VarType::Integer).unwrap();
    p.add_var("x2", VarType::Integer).unwrap();
    p.set_var_bounds(0, 0.0, 10.0).unwrap();
    p.set_var_bounds(1, 0.0, 10.0).unwrap();
    p.set_objective(vec![(0, 10.0), (1, 1.0)], 0.5).unwrap();
    p.add_constraint(None, vec![(0, 20.0), (1, 3.0)], RelOp::Le, 43.0)
        .unwrap();
    p
}

/// (1) CLI 默认参数输出逐字节基线：**仅原 7 算例**（threads=1 与默认 auto
/// 红线端到端钉死）。hard-knapsack-12 超自适应阈值、并行节点数可变，不入
/// 逐字节基线（见 (6)/(6b) 与 README）。
#[test]
fn acceptance_cli_default_output_identical() {
    let baselines: &[(&str, &str, i32)] = &[
        ("lp-prod-3", "status = optimal\nnodes = 1\nobjective = 13.000000\n", 0),
        ("knapsack-10", "status = optimal\nnodes = 1\nobjective = 17.000000\n", 0),
        ("assign-4x4", "status = optimal\nnodes = 1\nobjective = 10.000000\n", 0),
        // wolsey-2：整数目标（系数 5/8 × Integer）⇒ cutoff 生效；默认分支
        // 规则换 Pseudocost（含可靠性冷启动探测）后节点数随分支/压栈序
        // 变化（A+B 叠加实测钉值），目标 40 精确不变。
        // （WP6 包 A 复核：本包根割改动（4×40 + c-MIR/隐含界/提升）后实测
        // 仍为 7；禁用新割复算（同轮数纯 GMI+覆盖）亦为 7 ⇒ 该值当时对
        // 割族不敏感。R2-retry 聚合 MIR（「表格行⊕原始行」GMI 化）入批后
        // 根割界直接收紧到整数最优 40 ⇒ 根节点即证最优（nodes 7→1），
        // 目标 40 恒为已知整数最优不变——钉值随本轮割族增强重钉。）
        ("milp-wolsey-2", "status = optimal\nnodes = 1\nobjective = 40.000000\n", 0),
        ("milp-mixed-5", "status = optimal\nnodes = 1\nobjective = 29.000000\n", 0),
        ("infeasible-tiny", "status = infeasible\n", 2),
        ("unbounded-tiny", "status = unbounded\n", 3),
    ];
    for &(name, expected_stdout, expected_code) in baselines {
        let (code, stdout) = cli(&[&instance_path(name)]);
        assert_eq!(
            code,
            Some(expected_code),
            "算例 {name} 退出码，stdout：{stdout}"
        );
        assert_eq!(
            stdout, expected_stdout,
            "算例 {name} 默认参数输出必须与升级前逐字节一致"
        );
    }
}

/// (2) 库级 threads=4 与 threads=1（solve_default）一致：5 个原数值算例
/// status 相同且 |objective 差| ≤ 1e-6；infeasible/unbounded status 相同；
/// hard-knapsack-12（新算例文件）：Optimal 且 objective=523±1e-6
///（195 节点 ≫ 原算例的 1/7 节点，真实压满节点池路径）。
#[test]
fn acceptance_threads4_matches_threads1() {
    let ctrl4 = SearchControl {
        threads: 4,
        ..SearchControl::default()
    };
    for &(name, expected) in CASES {
        let problem = load_instance(Path::new(&instance_path(name)))
            .unwrap_or_else(|error| panic!("算例 {name} 加载失败：{error}"));
        let seq = solve_default(&problem)
            .unwrap_or_else(|error| panic!("算例 {name} 顺序求解失败：{error}"));
        let par = solve_controlled(&problem, &SolveOptions::default(), &ctrl4, None)
            .unwrap_or_else(|error| panic!("算例 {name} 并行求解失败：{error}"));
        assert_eq!(
            par.report.status, seq.status,
            "算例 {name}：threads=4 与 threads=1 终态必须相同"
        );
        assert_eq!(seq.status, SolveStatus::Optimal, "算例 {name} 顺序基线");
        let diff = (par
            .report
            .objective
            .unwrap_or_else(|| panic!("算例 {name} 并行无 objective"))
            - seq
                .objective
                .unwrap_or_else(|| panic!("算例 {name} 顺序无 objective")))
        .abs();
        assert!(diff <= 1e-6, "算例 {name}：objective 差 {diff} 超 1e-6");
        assert!(
            (par.report.objective.unwrap() - expected).abs() <= 1e-6,
            "算例 {name} 并行 objective 须达期望 {expected}"
        );
    }
    // 状态型算例：status 相同（Infeasible / Unbounded）。
    for name in ["infeasible-tiny", "unbounded-tiny"] {
        let problem = load_instance(Path::new(&instance_path(name)))
            .unwrap_or_else(|error| panic!("算例 {name} 加载失败：{error}"));
        let seq = solve_default(&problem)
            .unwrap_or_else(|error| panic!("算例 {name} 顺序求解失败：{error}"));
        let par = solve_controlled(&problem, &SolveOptions::default(), &ctrl4, None)
            .unwrap_or_else(|error| panic!("算例 {name} 并行求解失败：{error}"));
        assert_eq!(par.report.status, seq.status, "算例 {name}：终态必须相同");
        assert!(par.report.objective.is_none(), "算例 {name} 不得携带 objective");
    }
    // hard-knapsack-12：Optimal 且 523（2^12 穷举已独立验证，见 README）。
    let problem = load_instance(Path::new(&instance_path("hard-knapsack-12")))
        .unwrap_or_else(|error| panic!("算例 hard-knapsack-12 加载失败：{error}"));
    let par = solve_controlled(&problem, &SolveOptions::default(), &ctrl4, None)
        .unwrap_or_else(|error| panic!("算例 hard-knapsack-12 并行求解失败：{error}"));
    assert_eq!(par.report.status, SolveStatus::Optimal);
    let objective = par
        .report
        .objective
        .unwrap_or_else(|| panic!("hard-knapsack-12 报告 Optimal 但无 objective"));
    assert!(
        (objective - 523.0).abs() <= 1e-6,
        "hard-knapsack-12 objective = {objective}，期望 523"
    );
    assert!(!par.gap_terminated);
}

/// (3) CLI 级并发：knapsack-10 --threads 4 ⇒ 退出码 0、末行 objective=17±1e-6；
/// hard-knapsack-12 --threads 4 ⇒ 退出码 0、末行 objective=523±1e-6，重复
/// ≥10 次同断言（195 节点压满池；丢唤醒/死锁/竞态在重复下必现，单次亚毫秒
/// 耗时可忽略；断言不依赖计时）。
#[test]
fn acceptance_cli_threads4() {
    let (code, stdout) = cli(&[&instance_path("knapsack-10"), "--threads", "4"]);
    assert_eq!(code, Some(0), "knapsack-10 --threads 4 应报最优：{stdout}");
    assert!(
        (last_objective(&stdout) - 17.0).abs() <= 1e-6,
        "knapsack-10 --threads 4 objective：{stdout}"
    );

    for run in 1..=10 {
        let (code, stdout) = cli(&[&instance_path("hard-knapsack-12"), "--threads", "4"]);
        assert_eq!(
            code, Some(0),
            "hard-knapsack-12 --threads 4 第 {run} 次运行应报最优：{stdout}"
        );
        assert!(
            (last_objective(&stdout) - 523.0).abs() <= 1e-6,
            "hard-knapsack-12 --threads 4 第 {run} 次运行 objective：{stdout}"
        );
    }
}

/// (4) `--log-secs` 进度日志格式：
/// - lp-prod-3 --log-secs 5：首行逐字段精确
///   ` nodes=1 incumbent=13.000000 bound=13.000000 gap=0.00%`（t 数值不钉），
///   末行仍 objective = 13.000000；
/// - knapsack-10：首行 nodes=1、incumbent=17.000000；
/// - milp-wolsey-2：首行 incumbent=none、gap=inf；
/// - infeasible-tiny / unbounded-tiny：首行 ` nodes=1 incumbent=none bound=none gap=inf`；
/// - 全部算例：每个 progress 行匹配统一正则、≥1 行、首 progress 行先于
///   status 行；
/// - hard-knapsack-12 --log-secs 0.001 --threads 4 变体：周期行 ≥1（总行
///   数 ≥2 = 首行 + ≥1 周期行）、全部匹配正则、且无 `bound=none`
///   （空栈/空池矛盾行防护，N5）。
#[test]
fn acceptance_cli_progress_log_format() {
    // lp-prod-3：首行逐字段精确 + 末行 objective 不受进度行影响。
    let (code, stdout) = cli(&[&instance_path("lp-prod-3"), "--log-secs", "5"]);
    assert_eq!(code, Some(0), "{stdout}");
    assert_eq!(
        first_progress_tail(&stdout),
        " nodes=1 incumbent=13.000000 bound=13.000000 gap=0.00%",
        "lp-prod-3 首行字段：{stdout}"
    );
    assert_eq!(
        stdout.trim_end_matches(['\r', '\n']).lines().next_back(),
        Some("objective = 13.000000"),
        "末行必须仍是 objective 行：{stdout}"
    );
    assert_progress_lines_wellformed(&stdout);

    // knapsack-10：首行 nodes=1、incumbent=17.000000。
    let (code, stdout) = cli(&[&instance_path("knapsack-10"), "--log-secs", "5"]);
    assert_eq!(code, Some(0), "{stdout}");
    let tail = first_progress_tail(&stdout);
    assert!(tail.starts_with(" nodes=1 "), "knapsack-10 首行 nodes：{stdout}");
    assert!(
        tail.contains(" incumbent=17.000000 "),
        "knapsack-10 首行 incumbent：{stdout}"
    );
    assert_progress_lines_wellformed(&stdout);

    // milp-wolsey-2：根 dive 现在成功（深预算 + 一步回退）⇒ 首行即有
    // incumbent（2026-10 启发式增强前为 none/gap=inf——当时 8 次修复
    // 预算在断点链上必耗尽；数值 34→37 随根割族增强漂移，A+B 叠加实测
    // 钉值；R2-retry 聚合 MIR 入批后根割界收紧到整数最优 40 ⇒ 首行
    // incumbent 37→40，仍为已知整数最优、随本轮割族增强重钉）。
    let (code, stdout) = cli(&[&instance_path("milp-wolsey-2"), "--log-secs", "5"]);
    assert_eq!(code, Some(0), "{stdout}");
    let tail = first_progress_tail(&stdout);
    assert!(
        tail.contains(" incumbent=40.000000 "),
        "milp-wolsey-2 首行 incumbent：{stdout}"
    );
    assert!(
        !tail.ends_with(" gap=inf"),
        "milp-wolsey-2 首行 gap 应为数值（有 incumbent）：{stdout}"
    );
    assert_progress_lines_wellformed(&stdout);

    // infeasible-tiny / unbounded-tiny：首行四字段全钉。
    for name in ["infeasible-tiny", "unbounded-tiny"] {
        let (code, stdout) = cli(&[&instance_path(name), "--log-secs", "5"]);
        assert!(code == Some(2) || code == Some(3), "算例 {name} 退出码：{stdout}");
        assert_eq!(
            first_progress_tail(&stdout),
            " nodes=1 incumbent=none bound=none gap=inf",
            "算例 {name} 首行字段：{stdout}"
        );
        assert_progress_lines_wellformed(&stdout);
    }

    // hard-knapsack-12 变体：--log-secs 0.001 --threads 4（并行协调线程）。
    let (code, stdout) = cli(&[
        &instance_path("hard-knapsack-12"),
        "--log-secs",
        "0.001",
        "--threads",
        "4",
    ]);
    assert_eq!(code, Some(0), "{stdout}");
    let progress = assert_progress_lines_wellformed(&stdout);
    assert!(
        progress.len() >= 2,
        "极小周期下首行之后至少一行周期行，实际 {} 行：{stdout}",
        progress.len()
    );
    for line in &progress {
        assert!(
            !line.contains("bound=none"),
            "周期行不得出现 bound=none（空栈/空池矛盾行，N5）：{line:?}"
        );
    }
    assert!(
        (last_objective(&stdout) - 523.0).abs() <= 1e-6,
        "hard-knapsack-12 变体末行 objective：{stdout}"
    );
}

/// (5) `--mipgap`：
/// - `--mipgap 0` 于 5 个原数值算例 ⇒ 退出码 0 + objective==期望（≡关闭）；
/// - `--mipgap 100` 于 milp-wolsey-2 ⇒ 退出码 ∈ {0,4} 择一断言：0 ⇒ 末行
///   objective=40±1e-6；4 ⇒ 含 status = gap_reached、gap 行匹配
///   `gap=\d+\.\d+%`（非 inf、非负——clamp+证据双保证，捕获 -0.00%/inf 矛盾）、
///   best 行可解析且 ≤ 40+1e-6、无 objective 行；
/// - 库级确定性补充（dive 问题，关 presolve 与根割）：mipgap=Some(3.0) ⇒
///   LimitReached + gap_terminated + objective=21 + gap_percent≈2.38；
///   Some(2.0) ⇒ Optimal 21。
#[test]
fn acceptance_cli_mipgap() {
    // --mipgap 0 ≡ 关闭。
    for &(name, expected) in CASES {
        let (code, stdout) = cli(&[&instance_path(name), "--mipgap", "0"]);
        assert_eq!(code, Some(0), "算例 {name} --mipgap 0 应报最优：{stdout}");
        let objective = last_objective(&stdout);
        assert!(
            (objective - expected).abs() <= 1e-6,
            "算例 {name} --mipgap 0 objective = {objective}，期望 {expected}"
        );
    }

    // --mipgap 100 于 milp-wolsey-2：择一断言。
    let (code, stdout) = cli(&[&instance_path("milp-wolsey-2"), "--mipgap", "100"]);
    match code {
        Some(0) => {
            assert!(
                (last_objective(&stdout) - 40.0).abs() <= 1e-6,
                "退出码 0 分支：objective 须为 40：{stdout}"
            );
        }
        Some(4) => {
            assert!(
                stdout.contains("status = gap_reached"),
                "退出码 4 分支须含 gap_reached 状态行：{stdout}"
            );
            let gap_line = stdout
                .lines()
                .find(|l| l.starts_with("gap = "))
                .unwrap_or_else(|| panic!("退出码 4 分支须含 gap 行：{stdout}"));
            let gap_text = gap_line
                .strip_prefix("gap = ")
                .and_then(|v| v.strip_suffix('%'))
                .unwrap_or_else(|| panic!("gap 行格式非法：{gap_line:?}"));
            // 匹配 \d+\.\d+%：无 '-'、非 inf、恰一小数点、其余全数字。
            assert!(
                !gap_text.starts_with('-')
                    && gap_text.contains('.')
                    && gap_text
                        .chars()
                        .all(|c| c.is_ascii_digit() || c == '.')
                    && gap_text.matches('.').count() == 1,
                "gap 值须为非负有限百分比（如 9.12%），得到 {gap_line:?}"
            );
            let gap: f64 = gap_text
                .parse()
                .unwrap_or_else(|_| panic!("gap 无法数值解析：{gap_line:?}"));
            assert!(
                gap.is_finite() && gap >= 0.0,
                "gap 须非负有限（clamp+证据双保证）：{gap_line:?}"
            );
            let best_line = stdout
                .lines()
                .find(|l| l.starts_with("best = "))
                .unwrap_or_else(|| panic!("退出码 4 分支须含 best 行：{stdout}"));
            let best: f64 = best_line
                .strip_prefix("best = ")
                .unwrap_or_else(|| panic!("best 行格式非法：{best_line:?}"))
                .parse()
                .unwrap_or_else(|_| panic!("best 无法数值解析：{best_line:?}"));
            assert!(
                best <= 40.0 + 1e-6,
                "gap 终止的 incumbent 不得超过真最优 40：{best_line:?}"
            );
            assert!(
                !stdout.contains("objective = "),
                "gap 终止不得输出 objective 行：{stdout}"
            );
        }
        other => panic!("--mipgap 100 退出码应为 0 或 4，得到 {other:?}：{stdout}"),
    }

    // 库级确定性：dive 问题（gap=(22.0−21.5)/21.5×100≈2.326）。
    let problem = dive_gap_problem();
    let opts = SolveOptions {
        use_presolve: false,
        bb: BbOptions {
            root_cuts: false,
            heuristics: true,
            ..BbOptions::default()
        },
        ..SolveOptions::default()
    };
    let ctrl = SearchControl {
        mipgap_percent: Some(3.0),
        ..SearchControl::default()
    };
    let r = solve_controlled(&problem, &opts, &ctrl, None)
        .unwrap_or_else(|error| panic!("mipgap=3.0 求解失败：{error}"));
    assert_eq!(r.report.status, SolveStatus::LimitReached, "gap 终止诚实标注截断");
    assert!(r.gap_terminated);
    let objective = r
        .report
        .objective
        .expect("gap 终止须携带 incumbent 解的目标值");
    assert!(
        (objective - 21.5).abs() <= 1e-9,
        "mipgap=3.0 触发时 incumbent = {objective}，期望 21.5"
    );
    let gap = r.gap_percent.expect("gap_terminated ⇒ gap_percent 必为 Some");
    assert!(gap.is_finite() && gap >= 0.0, "gap 须非负有限：{gap}");
    assert!(
        (gap - 100.0 * 0.5 / 21.5).abs() < 1e-9,
        "gap_percent = {gap}，期望 ≈2.33"
    );
    // Some(2.0)：2.326 > 2 不触发 ⇒ Optimal 21.5。
    let ctrl = SearchControl {
        mipgap_percent: Some(2.0),
        ..SearchControl::default()
    };
    let r = solve_controlled(&problem, &opts, &ctrl, None)
        .unwrap_or_else(|error| panic!("mipgap=2.0 求解失败：{error}"));
    assert_eq!(r.report.status, SolveStatus::Optimal);
    assert!(!r.gap_terminated);
    assert_eq!(r.gap_percent, None);
    assert!((r.report.objective.unwrap() - 21.5).abs() <= 1e-9);
}

/// (6) 库层 adaptive_heuristics=true（auto 模式）：5 个原数值算例
/// **及 hard-knapsack-12** ⇒ Optimal + 期望目标（auto 于超阈值算例仍精确
/// 最优——dive 只供应 incumbent，不改剪枝/分支语义；nodes 不钉）。
#[test]
fn acceptance_heuristic_auto_optimal() {
    let ctrl = SearchControl {
        adaptive_heuristics: true,
        ..SearchControl::default()
    };
    for &(name, expected) in CASES {
        let problem = load_instance(Path::new(&instance_path(name)))
            .unwrap_or_else(|error| panic!("算例 {name} 加载失败：{error}"));
        let r = solve_controlled(&problem, &SolveOptions::default(), &ctrl, None)
            .unwrap_or_else(|error| panic!("算例 {name} auto 求解失败：{error}"));
        assert_eq!(r.report.status, SolveStatus::Optimal, "算例 {name}");
        assert!(!r.gap_terminated, "算例 {name} 不得 gap 终止");
        let objective = r
            .report
            .objective
            .unwrap_or_else(|| panic!("算例 {name} auto 无 objective"));
        assert!(
            (objective - expected).abs() <= 1e-6,
            "算例 {name} auto objective = {objective}，期望 {expected}"
        );
    }
    // hard-knapsack-12：195 节点超自适应阈值（32/50），auto 仍精确最优。
    let problem = load_instance(Path::new(&instance_path("hard-knapsack-12")))
        .unwrap_or_else(|error| panic!("算例 hard-knapsack-12 加载失败：{error}"));
    let r = solve_controlled(&problem, &SolveOptions::default(), &ctrl, None)
        .unwrap_or_else(|error| panic!("算例 hard-knapsack-12 auto 求解失败：{error}"));
    assert_eq!(r.report.status, SolveStatus::Optimal);
    assert!(!r.gap_terminated);
    let objective = r
        .report
        .objective
        .expect("hard-knapsack-12 auto 无 objective");
    assert!(
        (objective - 523.0).abs() <= 1e-6,
        "hard-knapsack-12 auto objective = {objective}，期望 523"
    );
}

/// (6b) CLI `--heuristic-mode on` 在 hard-knapsack-12 上的确定性基线
///（on ≡ heuristics=true 且不自适应）。整数目标（整数系数 × 全 Binary）⇒
/// cutoff=inc+1−ε 收紧 (d) 界 + 陈旧前沿 est 弹出剪枝；根覆盖割（WP4 增补，
/// cuts.rs）对全 Binary 背包行再收紧；默认分支规则换 Pseudocost（含可靠性
/// 冷启动探测 + best-estimate 压栈序）改变树形；objective 523.000000 钉死
/// 不变（2^12 穷举真值）。默认 auto 的 nodes 不钉（见 README：阈值超限后
/// 控制器可补跑 dive；目标不变）。nodes 为 A+B 叠加实测钉值。
#[test]
fn acceptance_cli_heuristic_on_identical() {
    let (code, stdout) = cli(&[&instance_path("hard-knapsack-12"), "--heuristic-mode", "on"]);
    assert_eq!(code, Some(0), "{stdout}");
    // 演进：根 dive 供应 incumbent 227→235；整数目标 cutoff 235→211；根
    // 覆盖割 211→195；默认 Pseudocost + 可靠性探测 + best-estimate 压栈序
    // （A+B 叠加实测）195→201（目标不变——树形受割与分支序交互影响）。
    // （WP6 包 A 复核归因：本包割轮 5→4 + c-MIR/提升族后实测仍为 201；
    // 禁用新割复算（同 4 轮纯 GMI+覆盖）得 93，割轮保持 5 时新割族下为
    // 73——树形对割集混沌敏感，201 为 A+B 叠加终态实测钉值。）
    // （2026-10 根预割启发式轮：割后顶点中面化使根 incumbent 453 质量塌
    // 缩，割前顶点补跑一轮 round_and_check + dive 后根 incumbent 511 ⇒
    // 剪枝提前，201→184（目标不变）。）
    // （2026-10 复核后默认规则回退 MostFractional：184→159。）
    // （R3 类自适应默认：hard-knapsack-12 命中背包类检测（全 Binary +
    // 正系数 Le 单行）⇒ 默认翻 Pseudocost+probes8（定档矩阵实测该口径
    // 在背包类上净优：mknap-5x300 −25% 证最优墙钟；本例树形 159→196，
    // 目标 523 恒为已知最优——钉值随默认规则演进重钉；显式复现旧行：
    // --branch-rule most_fractional。）
    assert_eq!(
        stdout,
        "status = optimal
nodes = 196
objective = 523.000000
",
        "hard-knapsack-12 --heuristic-mode on 输出（R3 类自适应默认确定性基线）"
    );
}

/// 契约：auto 模式 + 充足时间预算（默认降级参数 0.5/0.25/0.01）在全部
/// 验收算例上仍证明最优——默认策略不会在秒级完成的算例上误触发降级
/// （资格预检要求预算过半 + 双停滞窗口）；降级只应在真实长停滞时放弃
/// 证明。同时顺序（threads=1）与并行（threads=4）口径一致。
#[test]
fn acceptance_auto_with_budget_still_optimal() {
    for threads in [1usize, 4] {
        let ctrl = SearchControl {
            threads,
            adaptive_heuristics: true,
            ..SearchControl::default()
        };
        for &(name, expected) in CASES {
            let path = instance_path(name);
            let problem = load_instance(std::path::Path::new(&path))
                .unwrap_or_else(|error| panic!("算例 {name} 加载失败：{error}"));
            let mut opts = SolveOptions::default();
            opts.bb.time_limit_ms = 60_000;
            let r = solve_controlled(&problem, &opts, &ctrl, None)
                .unwrap_or_else(|error| panic!("算例 {name}（threads={threads}）求解失败：{error}"));
            assert_eq!(
                r.report.status,
                SolveStatus::Optimal,
                "算例 {name}（threads={threads}）：默认降级参数不得误触发"
            );
            assert!(!r.auto_degraded, "算例 {name}（threads={threads}）不应降级");
            let objective = r
                .report
                .objective
                .unwrap_or_else(|| panic!("算例 {name} 报告 Optimal 但无 objective"));
            assert!(
                (objective - expected).abs() <= 1e-6,
                "算例 {name}（threads={threads}）：objective = {objective}，期望 {expected}"
            );
        }
    }
}

/// bench 算例路径（`bench/instances/<name>.milp`，milp crate 根的上一级；
/// bench/ 归画像与复测所有，本文件只读引用，不修改）。
fn bench_instance_path(name: &str) -> String {
    format!("{}/../bench/instances/{name}.milp", env!("CARGO_MANIFEST_DIR"))
}

/// （P0-1 验收 + 割轮 IterationLimit 分流）setcover-100x150：根割循环
/// 的割后 resolve 触发迭代上限时不再作整树终态——丢弃全部根割、重建
/// 无割 Simplex 恢复 Optimal 根态后继续搜索。修复前 6/6 次运行 1 节点
/// limit_reached 无解（bench/raw_runs.jsonl 第 3 轮基线）；修复后：
/// - t1/t4 均 Optimal 且目标 276（节点数 > 1：搜索越过根节点）；
/// - threads≥2 + adaptive 与 t1 目标一致（子树不丢、不误终态）；
/// - 记账完整：默认路径 lp_iterations **严格大于**关割等价路径（回退后
///   树与关割一致，多出的恰是根 LP + 割循环的重建前迭代）；
/// - 割丢弃如实：stats.cuts == 0（割已随 Simplex 一起丢弃）。
///
/// 重负载门禁（L4）：debug 构建下 8 次 setcover-100x150 求解（t1 + 关割
/// + threads∈{2,4}×3 重复，含并行根预割启发式轮）实测 ~100s、release
/// ~1s/次——默认 `cargo test` 跳过（#[ignore]），独立门禁显式运行：
/// `cargo test --test acceptance -- --ignored`；脚本门禁另跑 release 全量
/// 验收兜底。验收覆盖不丢，只是移出默认快路径。
#[ignore = "setcover-100x150 重负载：debug ~100s；门禁显式 `cargo test -- --ignored` 运行（release 验收兜底）"]
#[test]
fn acceptance_setcover_100x150_unlocked_optimal() {
    let path = bench_instance_path("setcover-100x150");
    let problem = load_instance(Path::new(&path))
        .unwrap_or_else(|error| panic!("算例 setcover-100x150 加载失败：{error}"));

    // t1 默认路径：Optimal 276。
    let seq = solve_default(&problem)
        .unwrap_or_else(|error| panic!("setcover-100x150 t1 求解失败：{error}"));
    assert_eq!(seq.status, SolveStatus::Optimal, "t1 终态");
    assert!(
        (seq.objective.unwrap() - 276.0).abs() <= 1e-6,
        "t1 objective = {:?}，期望 276",
        seq.objective
    );
    assert!(seq.stats.nodes > 1, "搜索越过根节点（实际 {} 节点）", seq.stats.nodes);
    // 割保留如实（R8 重钉：割循环收紧预算 + 回滚语义——失败轮弃增量割、
    // 已提交割保留，不再触发 branch_bound 的整树回退。旧钉值 0 对应
    // 『回退丢弃全部根割』，R8 后实测 cuts=23（已提交割数）、树 45→13
    // 节点；目标 276 恒为已知最优不变）。
    assert_eq!(seq.stats.cuts, 23, "回滚语义：已提交割保留（R8 重钉）");

    // 记账完整：lp_iterations 含重建前迭代（根 LP + 割循环）。
    let opts_nocuts = SolveOptions {
        bb: BbOptions {
            root_cuts: false,
            ..BbOptions::default()
        },
        ..SolveOptions::default()
    };
    let nocuts = solve(&problem, &opts_nocuts)
        .unwrap_or_else(|error| panic!("setcover-100x150 关割求解失败：{error}"));
    assert_eq!(nocuts.status, SolveStatus::Optimal);
    assert!((nocuts.objective.unwrap() - 276.0).abs() <= 1e-6);
    assert!(
        seq.stats.lp_iterations > nocuts.stats.lp_iterations,
        "默认（回退）lp_iterations={} 须 > 关割={}（含重建前根割迭代）",
        seq.stats.lp_iterations,
        nocuts.stats.lp_iterations
    );

    // threads≥2 + adaptive 与 t1 目标一致（重复 3 次：并行口径下每次都
    // 不得丢子树/误终态）。
    for threads in [2usize, 4] {
        let ctrl = SearchControl {
            threads,
            adaptive_heuristics: true,
            ..SearchControl::default()
        };
        for rep in 0..3 {
            let r = solve_controlled(&problem, &SolveOptions::default(), &ctrl, None)
                .unwrap_or_else(|error| {
                    panic!("setcover-100x150 threads={threads} rep={rep} 求解失败：{error}")
                });
            assert_eq!(
                r.report.status,
                SolveStatus::Optimal,
                "threads={threads} rep={rep} 终态"
            );
            assert!(
                (r.report.objective.unwrap() - 276.0).abs() <= 1e-6,
                "threads={threads} rep={rep}：objective 须与 t1 一致（276）"
            );
            assert!(!r.gap_terminated && !r.auto_degraded);
        }
    }
}

/// （确定性验收）固定 node_limit（默认 2M）、time_limit=0 下默认路径
/// （Pseudocost + 可靠性探测预算 8）两次运行逐项一致：探测预算按
/// 节点/迭代计、绝不读墙钟——无时间上限时自适应/降级的时间面也不
/// 参与决策，搜索全确定。
///
/// 重负载门禁（L4，与 setcover 主验收同因）：2 次 setcover-100x150
/// 默认求解在 debug 下实测 ~33s（release 亚秒）——#[ignore] 移出默认
/// 快路径，门禁 `cargo test --test acceptance -- --ignored` 显式运行。
#[ignore = "setcover-100x150 重负载（2 次默认求解 debug ~33s）；门禁显式 `cargo test -- --ignored` 运行"]
#[test]
fn acceptance_default_path_deterministic_two_runs() {
    let path = bench_instance_path("setcover-100x150");
    let problem = load_instance(Path::new(&path)).unwrap();
    let r1 = solve_default(&problem).unwrap();
    let r2 = solve_default(&problem).unwrap();
    assert_eq!(r1.status, r2.status);
    assert_eq!(r1.objective, r2.objective);
    assert_eq!(r1.stats.nodes, r2.stats.nodes, "节点数跨次一致");
    assert_eq!(r1.stats.lp_iterations, r2.stats.lp_iterations);

    // 本地算例同断言（不同规模/结构）。
    let problem = load_instance(Path::new(&instance_path("hard-knapsack-12"))).unwrap();
    let r1 = solve_default(&problem).unwrap();
    let r2 = solve_default(&problem).unwrap();
    assert_eq!(r1.stats.nodes, r2.stats.nodes);
    assert_eq!(r1.stats.lp_iterations, r2.stats.lp_iterations);
    assert_eq!(r1.objective, r2.objective);
}

/// （P1-4 纯被动等价验收）pc_reliability_probes=0 ≡ 纯被动 Pseudocost：
/// 探测零动作 ⇒ 与无探测实现逐节点一致（两次运行 + 两个入口逐项相等，
/// 机制性保证的行为级复验）；默认预算 8 仅改变统计冷启动与树形——
/// 目标值精确不变（搜索完整性：探测只读界并恢复，不剪任何可行解）。
#[test]
fn acceptance_pc_probes_zero_passive_parity() {
    let problem = load_instance(Path::new(&instance_path("hard-knapsack-12"))).unwrap();
    let passive = SolveOptions {
        bb: BbOptions {
            branch_rule: milp::branch_bound::BranchRule::Pseudocost,
            pc_reliability_probes: 0,
            ..BbOptions::default()
        },
        ..SolveOptions::default()
    };
    // 纯被动两次运行逐项一致（确定性）。
    let p1 = solve(&problem, &passive).unwrap();
    let p2 = solve(&problem, &passive).unwrap();
    assert_eq!(p1.status, SolveStatus::Optimal);
    assert_eq!(p1.stats.nodes, p2.stats.nodes, "纯被动逐节点一致");
    assert_eq!(p1.stats.lp_iterations, p2.stats.lp_iterations);
    assert_eq!(p1.objective, p2.objective);
    // 入口等价：solve（默认受控委托）与 solve_controlled(默认) 同逐项。
    let p3 = solve_controlled(&problem, &passive, &SearchControl::default(), None).unwrap();
    assert_eq!(p3.report.stats.nodes, p1.stats.nodes);
    assert_eq!(p3.report.stats.lp_iterations, p1.stats.lp_iterations);
    // 默认预算 8：目标值与纯被动精确一致（树形允许不同）。
    let probed = solve_default(&problem).unwrap();
    assert_eq!(probed.status, SolveStatus::Optimal);
    assert!(
        (probed.objective.unwrap() - p1.objective.unwrap()).abs() <= 1e-6,
        "探测预算只改统计不剪枝 ⇒ 目标必须一致"
    );
}
