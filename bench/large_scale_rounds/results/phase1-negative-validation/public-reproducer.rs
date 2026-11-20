// 有限合法输入的第一阶段负人工变量公开接口复现；不写内部状态。
use milp::model::{Problem, RelOp, Sense, Tolerances, VarType};
use milp::simplex::{Bounds, LpSolution, Simplex};
use milp::simplex_sparse::SparseSimplex;

fn report(engine: &str, sense: Sense, larger_first: bool, problem: &Problem,
          solution: &LpSolution, iterations: u64, basis: &[usize], basic_values: &[f64]) {
    // 直接从原问题系数重算等式残差；不复用求解器表格或可行性函数。
    let residuals: Vec<f64> = problem.constraints.iter().map(|row| {
        row.expr.terms.iter().map(|&(j, a)| a * solution.x[j]).sum::<f64>() - row.rhs
    }).collect();
    let max_residual = residuals.iter().fold(0.0_f64, |m, r| m.max(r.abs()));
    let bound_violation = problem.vars.iter().enumerate().fold(0.0_f64, |m, (j, var)| {
        m.max((var.lb - solution.x[j]).max(0.0)).max((solution.x[j] - var.ub).max(0.0))
    });
    println!("engine={engine} sense={sense:?} larger_first={larger_first} status={:?} objective={:?} x={:?} iterations={iterations} basis={basis:?} basic_values={basic_values:?} original_residuals={residuals:?} max_original_residual={max_residual:.18e} bound_violation={bound_violation:.18e}",
        solution.status, solution.objective, solution.x);
}

fn main() {
    let difference = 2.0_f64.powi(-21);
    let tolerance = Tolerances::default();
    println!("difference={difference:.18e} feasibility={:.18e} difference_over_feasibility={:.9} ratio_band={:.18e}",
        tolerance.feasibility, difference / tolerance.feasibility,
        (tolerance.feasibility * 1001.0).min(tolerance.feasibility * 1000.0));
    for sense in [Sense::Minimize, Sense::Maximize] {
        for larger_first in [true, false] {
            let mut problem = Problem::new("phase1_negative_artificial", sense);
            let x = problem.add_var("x", VarType::Continuous).unwrap();
            problem.set_objective(vec![], 0.0).unwrap();
            let rhs = if larger_first { [1000.0 + difference, 1000.0] }
                else { [1000.0, 1000.0 + difference] };
            for b in rhs {
                problem.add_constraint(None, vec![(x, 1.0)], RelOp::Eq, b).unwrap();
            }
            problem.validate().unwrap();
            let bounds = Bounds::from_problem(&problem);
            let mut dense = Simplex::new(&problem, tolerance).unwrap();
            let dense_solution = dense.solve(&bounds);
            let dense_values: Vec<f64> = (0..dense.num_rows()).map(|r| dense.tableau_row_terms(r).1).collect();
            report("Dense", sense, larger_first, &problem, &dense_solution, dense.iterations(), dense.basis_columns(), &dense_values);
            let mut sparse = SparseSimplex::new(&problem, tolerance).unwrap();
            let sparse_solution = sparse.solve(&bounds);
            let sparse_values: Vec<f64> = (0..sparse.num_rows()).map(|r| sparse.tableau_row_terms(r).1).collect();
            report("Revised", sense, larger_first, &problem, &sparse_solution, sparse.iterations(), sparse.basis_columns(), &sparse_values);
        }
    }
}
