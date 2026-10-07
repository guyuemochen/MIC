//! 验收测试（WP5 契约）。
//!
//! - 5 个数值型算例走 Optimal 主路径：终态、目标值、维数、整数性、
//!   约束满足性（容差 1e-6×(1+|rhs|+‖row‖₁)）、目标按 values 重算一致；
//! - 2 个状态型算例：Infeasible / Unbounded（均无 objective/values）；
//! - CLI 层：objective 末行数值解析（容差 1e-6，不依赖具体小数位数）、
//!   不可行/无界退出码 2/3 且无 objective 行。
//!
//! 算例路径拼接用 `concat!(env!("CARGO_MANIFEST_DIR"), "/instances/", name, ".milp")`
//! 的等价形式（name 非字面量处用 format!）。LimitReached 退出码 4 不做 CLI
//! 级验收（依赖计时不稳定；由 WP4 单元测试以极小 node_limit 覆盖）。

use std::path::Path;
use std::process::Command;

use milp::instance::load_instance;
use milp::model::{RelOp, SolveStatus};
use milp::solver::solve_default;

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
