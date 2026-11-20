//! 仅用公开LP API复现比值测试边界，记录原始约束残差；不运行MILP搜索。
use milp::model::{Problem, RelOp, Sense, Tolerances, VarType};
use milp::simplex::{Bounds, LpSolution, Simplex};
use milp::simplex_sparse::SparseSimplex;
use std::panic::{catch_unwind, AssertUnwindSafe};

fn model(sense: Sense, reverse: bool, overflow: bool) -> Problem {
    let mut p = Problem::new("比值安全自然反例", sense);
    p.add_var("x", VarType::Continuous).unwrap();
    let sign = if sense == Sense::Maximize { 1.0 } else { -1.0 };
    if overflow {
        p.vars[0].lb = -1e308;
        p.vars[0].ub = 1e308;
        p.set_objective(vec![(0, sign)], 0.0).unwrap();
    } else {
        p.add_var("y", VarType::Continuous).unwrap();
        p.set_objective(vec![(0, sign), (1, sign * 0.1)], 0.0).unwrap();
        let mut rows = vec![(vec![(0, 1.0)], 1000.0 + 2.0_f64.powi(-21)),
                            (vec![(0, 1.0), (1, 0.5)], 1000.0)];
        if reverse { rows.reverse(); }
        for (terms, rhs) in rows {
            p.add_constraint(None, terms, RelOp::Le, rhs).unwrap();
        }
    }
    p.validate().unwrap();
    p
}

fn run(p: &Problem, sparse: bool) -> LpSolution {
    let bounds = Bounds::from_problem(p);
    if sparse {
        SparseSimplex::new(p, Tolerances::default()).unwrap().solve(&bounds)
    } else {
        Simplex::new(p, Tolerances::default()).unwrap().solve(&bounds)
    }
}

fn report(p: &Problem, sparse: bool, reverse: bool, overflow: bool) {
    let engine = if sparse { "revised" } else { "tableau" };
    let case = if overflow { "finite_width_overflow" } else { "phase2_negative_slack" };
    let expected_x = if overflow { 1e308 } else { 1000.0 };
    let expected_objective = if p.sense == Sense::Maximize { expected_x } else { -expected_x };
    println!("CASE case={case} engine={engine} sense={:?} reversed={reverse} debug_assertions={} expected_objective={expected_objective:.17e}", p.sense, cfg!(debug_assertions));
    match catch_unwind(AssertUnwindSafe(|| run(p, sparse))) {
        Ok(sol) => {
            let recomputed = p.obj.terms.iter().map(|&(j, a)| a * sol.x[j]).sum::<f64>() + p.obj_offset;
            println!("RESULT caught_panic=false status={:?} objective={:.17e} recomputed={recomputed:.17e} x={:?} x_bits={:?}",
                sol.status, sol.objective, sol.x, sol.x.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
            for (i, row) in p.constraints.iter().enumerate() {
                let activity: f64 = row.expr.terms.iter().map(|&(j, a)| a * sol.x[j]).sum();
                println!("ROW index={i} activity={activity:.17e} rhs={:.17e} le_violation={:.17e}", row.rhs, activity - row.rhs);
            }
            for (j, var) in p.vars.iter().enumerate() {
                println!("BOUND index={j} lower_violation={:.17e} upper_violation={:.17e}", var.lb - sol.x[j], sol.x[j] - var.ub);
            }
        }
        Err(payload) => {
            let message = payload.downcast_ref::<String>().cloned()
                .or_else(|| payload.downcast_ref::<&str>().map(|value| (*value).to_string()))
                .unwrap_or_else(|| "非字符串panic载荷".to_string());
            println!("RESULT caught_panic=true message={message:?}");
        }
    }
    println!("END");
}

fn main() {
    for overflow in [false, true] {
        for sense in [Sense::Maximize, Sense::Minimize] {
            for sparse in [false, true] {
                for reverse in [false, true] {
                    if overflow && reverse { continue; }
                    report(&model(sense, reverse, overflow), sparse, reverse, overflow);
                }
            }
        }
    }
}
