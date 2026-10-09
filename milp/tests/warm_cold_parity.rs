//! 热启动（对偶修复）与冷启动全解的一致性差分（界收紧序列）。
//!
//! 背景：resolve 的对偶修复快路径在大规模分支下是主路径；本测试在随机
//! 背包/混合整数模型的 DFS 式界收紧序列上，逐点对比热解与冷解的状态与
//! 目标值——任何对偶修复的语义偏差（目标分叉/状态分叉）都会在此暴露。
use milp::model::{Problem, RelOp, Sense, Tolerances, VarType};
use milp::simplex::{Bounds, Simplex};

fn knapsack(n: usize, seed: u64) -> Problem {
    let mut lcg = seed;
    let mut rnd = || {
        lcg = lcg.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((lcg >> 33) % 41 + 7) as f64
    };
    let mut p = Problem::new("k", Sense::Maximize);
    for i in 0..n {
        p.add_var(&format!("x{i}"), if i % 3 == 0 { VarType::Integer } else { VarType::Binary }).unwrap();
        p.set_var_bounds(i, 0.0, if i % 3 == 0 { 5.0 } else { 1.0 }).unwrap();
    }
    let terms: Vec<(usize, f64)> = (0..n).map(|i| (i, rnd())).collect();
    let cap: f64 = terms.iter().map(|&(_, w)| w).sum::<f64>() * 0.4;
    p.set_objective((0..n).map(|i| (i, terms[i].1 + 3.0)).collect(), 0.0).unwrap();
    p.add_constraint(None, terms, RelOp::Le, cap).unwrap();
    p
}

#[test]
fn dual_repair_matches_cold_objective() {
    let tol = Tolerances::default();
    for seed in 1..=40u64 {
        let p = knapsack(60, seed);
        let mut warm = Simplex::new(&p, tol).unwrap();
        let mut bounds = Bounds::from_problem(&p);
        let root = warm.solve(&bounds);
        assert!(matches!(root.status, milp::simplex::LpStatus::Optimal), "seed {seed} 根非 Optimal");
        // DFS 式收紧 60 步：每步取最分数变量收紧一侧界。
        for step in 0..60 {
            let sol = warm.resolve(&bounds);
            let frac = (0..p.vars.len())
                .filter(|&j| p.vars[j].var_type != VarType::Continuous)
                .map(|j| (j, (sol.x[j] - sol.x[j].floor()).abs()))
                .filter(|&(_, f)| f > 1e-6 && f < 1.0 - 1e-6)
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
            let Some((j, _)) = frac else { break };
            let v = sol.x[j];
            if step % 2 == 0 {
                bounds.ub[j] = bounds.ub[j].min(v.floor());
            } else {
                bounds.lb[j] = bounds.lb[j].max(v.ceil());
            }
            if !bounds.lb[j].is_finite() { break; }
            // 热解 vs 冷解 对比。
            let warm_sol = warm.resolve(&bounds);
            let mut cold = Simplex::new(&p, tol).unwrap();
            let cold_sol = cold.solve(&bounds);
            match (&warm_sol.status, &cold_sol.status) {
                (a, b) if std::mem::discriminant(a) == std::mem::discriminant(b) => {}
                (a, b) => panic!("seed {seed} step {step}: 状态分叉 warm={a:?} cold={b:?}"),
            }
            if matches!(warm_sol.status, milp::simplex::LpStatus::Optimal) {
                let dw = warm_sol.objective;
                let dc = cold_sol.objective;
                assert!(
                    (dw - dc).abs() <= 1e-7 * (1.0 + dc.abs()),
                    "seed {seed} step {step}: 目标分叉 warm={dw} cold={dc}（差 {:.3e}）",
                    dw - dc
                );
            }
        }
    }
}
