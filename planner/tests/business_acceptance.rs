//! 从实际入口求解，用独立的业务计算验证可行性、费用、完工时间和终止状态。

#[path = "../examples/business_models/cases.rs"]
mod cases;

use planner::business::{optimize, Config, Expr, Model, Outcome, Rule, StopReason};
use planner::Metaheuristic;
use std::time::Duration;

fn delivery_cost(routes: &[Vec<usize>]) -> i64 {
    let positions = [1i64, 3, 8];
    let mut total = 0;
    for route in routes {
        let mut last = 0;
        for &i in route {
            total += (positions[i] - last).abs();
            last = positions[i];
        }
        total += last.abs();
    }
    total
}

fn delivery_optimum() -> i64 {
    let mut best = i64::MAX;
    for a in 0..3 {
        for b in 0..3 {
            if a == b {
                continue;
            }
            let c = 3 - a - b;
            let customers = [a, b, c];
            // 每车最多接两个订单，三单只有 1+2、2+1 两类分配。
            for k in 1..=2 {
                best = best.min(delivery_cost(&[
                    customers[..k].to_vec(),
                    customers[k..].to_vec(),
                ]));
            }
        }
    }
    best
}

fn factory_makespan(machines: &[Vec<usize>]) -> Option<i64> {
    // 用逐项放行来模拟生产：工艺前驱和同机前一项都已完成，才安排此工序。
    let durations = [2i64, 1, 2, 1];
    let job_prev = [None, Some(0), None, Some(2)];
    let mut ends = [None; 4];
    for _ in 0..4 {
        for machine in machines {
            for (p, &i) in machine.iter().enumerate() {
                let machine_prev = p.checked_sub(1).map(|q| machine[q]);
                let predecessors: Vec<_> = [machine_prev, job_prev[i]]
                    .iter()
                    .flatten()
                    .copied()
                    .collect();
                if predecessors.iter().all(|&j| ends[j].is_some()) {
                    let start = predecessors
                        .iter()
                        .filter_map(|&j| ends[j])
                        .max()
                        .unwrap_or(0);
                    ends[i] = Some(start + durations[i]);
                }
            }
        }
    }
    ends.iter()
        .copied()
        .collect::<Option<Vec<_>>>()?
        .into_iter()
        .max()
}

fn factory_optimum() -> i64 {
    let mut best = i64::MAX;
    for m0 in [vec![0, 3], vec![3, 0]] {
        for m1 in [vec![1, 2], vec![2, 1]] {
            if let Some(t) = factory_makespan(&[m0.clone(), m1]) {
                best = best.min(t);
            }
        }
    }
    best
}

fn check_result(name: &str, r: &Outcome) {
    assert!(r.complete, "{name} 未生成完整解");
    assert_eq!(r.score.hard, 0, "{name} 存在硬违反");
    match name {
        "roster" => {
            let shifts: Vec<_> = r.state.values.iter().map(|v| v.unwrap()).collect();
            assert_eq!(shifts.len(), 3);
            assert_eq!(shifts[2], 1);
            assert!(shifts.iter().all(|x| [0, 1].contains(x)));
            for id in [0, 1] {
                assert!(shifts.iter().filter(|&&v| v == id).count() <= 2);
            }
            let preferences = shifts
                .iter()
                .map(|&v| if v == 1 { 0 } else { 2 })
                .collect::<Vec<_>>();
            assert_eq!(
                r.score.penalties,
                [vec![0, 0], preferences.clone()].concat()
            );
            assert_eq!(preferences.iter().sum::<i64>(), r.score.soft);
            assert_eq!(r.score.soft, 2);
        }
        "timetable" => {
            let vars: Vec<_> = r.state.values.iter().map(|v| v.unwrap()).collect();
            assert_eq!(vars.len(), 6);
            let lessons = vars.chunks_exact(2).collect::<Vec<_>>();
            assert_ne!(lessons[0][0], lessons[1][0]);
            for i in 0..3 {
                for j in i + 1..3 {
                    assert_ne!(lessons[i], lessons[j]);
                }
            }
            assert_eq!(r.score.soft, lessons.iter().map(|l| l[0]).sum::<i64>());
            assert_eq!(r.score.soft, 1);
        }
        "delivery" => {
            let routes = &r.state.lists[0];
            assert_eq!(routes.len(), 2);
            assert!(routes.iter().all(|r| r.len() <= 2));
            let mut ids = routes.iter().flatten().copied().collect::<Vec<_>>();
            ids.sort_unstable();
            assert_eq!(ids, vec![0, 1, 2]);
            assert_eq!(delivery_cost(routes), r.score.soft);
            assert_eq!(r.score.soft, delivery_optimum());
        }
        "factory" => {
            let machines = &r.state.lists[0];
            assert_eq!(machines.len(), 2);
            let mut jobs = machines.clone();
            for m in &mut jobs {
                m.sort_unstable();
            }
            assert_eq!(jobs, vec![vec![0, 3], vec![1, 2]]);
            assert_eq!(factory_makespan(machines), Some(r.score.soft));
            assert_eq!(r.score.soft, factory_optimum());
        }
        other => panic!("未知算例 {other}"),
    }
}

#[test]
fn all_algorithms_solve_four_business_models() {
    assert_eq!(delivery_optimum(), 18);
    assert_eq!(factory_optimum(), 3);
    for algo in [
        Metaheuristic::TabuSearch { tabu_size: 8 },
        Metaheuristic::LateAcceptance { history_len: 10 },
        Metaheuristic::SimulatedAnnealing {
            start_temp: 5.0,
            end_temp: 0.1,
        },
    ] {
        for (name, m) in cases::models() {
            let cfg = Config {
                algorithm: algo.clone(),
                steps: 250,
                seed: 42,
                ..Config::default()
            };
            let r = optimize(&m, &cfg).unwrap();
            assert_eq!(r.steps, 250);
            assert_eq!(r.reason, StopReason::Steps);
            check_result(name, &r);
        }
    }
}

#[test]
fn step_budgets_are_reproducible_with_one_and_three_workers() {
    for (name, m) in cases::models() {
        let mut best_single = None;
        for threads in [1, 3] {
            let cfg = Config {
                threads,
                steps: 100,
                seed: 67,
                ..Config::default()
            };
            let a = optimize(&m, &cfg).unwrap();
            let b = optimize(&m, &cfg).unwrap();
            assert_eq!(a.state, b.state, "{name}、{threads} 个线程");
            assert_eq!(a.score, b.score);
            assert_eq!((a.steps, a.reason), (b.steps, b.reason));
            let rank = (a.score.hard, a.score.soft);
            if let Some(old) = best_single {
                assert!(rank <= old);
            } else {
                best_single = Some(rank);
            }
        }
    }
}

#[test]
fn stop_reason_and_solution_quality_are_independent() {
    for (_name, m) in cases::models() {
        let r = optimize(
            &m,
            &Config {
                time_limit: Some(Duration::ZERO),
                ..Config::default()
            },
        )
        .unwrap();
        assert_eq!(r.reason, StopReason::Time);
        assert!(!r.complete);
        assert_eq!(r.steps, 0);
        let r = optimize(
            &m,
            &Config {
                steps: 0,
                ..Config::default()
            },
        )
        .unwrap();
        assert_eq!(r.reason, StopReason::Steps);
        assert!(r.complete);
    }
    let mut m = Model::default();
    m.domains.push(vec![1]);
    let violation = m.expr(Expr::Const(4));
    m.rules.push(Rule {
        name: "无法消除的罚项".into(),
        hard: true,
        weight: 2,
        expr: violation,
    });
    let r = optimize(&m, &Config::default()).unwrap();
    assert_eq!((r.steps, r.reason), (0, StopReason::NoMoves));
    assert!(r.complete);
    assert_eq!(r.score.hard, 8);
    assert_eq!(r.score.penalties, vec![8]);
    assert!(optimize(
        &m,
        &Config {
            threads: 0,
            ..Config::default()
        }
    )
    .is_err());
    m.domains[0].clear();
    assert!(optimize(&m, &Config::default()).is_err());
}

#[test]
fn tied_business_values_follow_registration_order() {
    let m = Model {
        domains: vec![vec![20, 10]],
        ..Model::default()
    };
    for construction in [
        planner::ConstructionHeuristic::FirstFit,
        planner::ConstructionHeuristic::BestFit,
    ] {
        let r = optimize(
            &m,
            &Config {
                construction,
                steps: 0,
                ..Config::default()
            },
        )
        .unwrap();
        assert_eq!(r.state.values, [Some(20)]);
        assert!(r.feasible());
    }
}
