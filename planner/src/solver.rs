//! 编排：validate → 构造 →（可选）局部搜索 → [`SearchResult`]。

use std::sync::Arc;

use crate::construct::{construct, ConstructionHeuristic};
use crate::localsearch::{local_search, LocalSearchConfig, SearchResult};
use crate::model::{PlannerError, ProblemModel};
use crate::score::ScoreCalculator;

/// 求解器配置（默认值全文唯一来源）。
#[derive(Debug, Clone)]
pub struct SolverConfig {
    /// 构造启发式；默认 [`ConstructionHeuristic::FirstFit`]。
    pub construction: ConstructionHeuristic,
    /// 局部搜索配置；默认 [`LocalSearchConfig::default`]。
    pub local_search: LocalSearchConfig,
    /// 是否执行局部搜索；默认 true（false = 仅构造）。
    pub run_local_search: bool,
}

impl Default for SolverConfig {
    /// 默认：FirstFit / LocalSearchConfig::default() / true。
    fn default() -> Self {
        SolverConfig {
            construction: ConstructionHeuristic::FirstFit,
            local_search: LocalSearchConfig::default(),
            run_local_search: true,
        }
    }
}

/// 求解编排：`validate → construct`（Err 透传）`→（可选）local_search → SearchResult`。
///
/// 恒有 `best_score ≤ initial_score`（字典序），`best` 为完整解。
/// `run_local_search = false` 时返回构造产物本身（steps = 0，best == 初始解）。
pub fn solve(model: Arc<ProblemModel>, cfg: &SolverConfig) -> Result<SearchResult, PlannerError> {
    model.validate()?;
    let initial = construct(&model, cfg.construction)?;
    if cfg.run_local_search {
        return local_search(initial, &cfg.local_search);
    }
    // 仅构造：不读取 local_search 配置（含其参数合法性——未被使用即不校验）
    let score = ScoreCalculator::new(model).total_score(&initial);
    Ok(SearchResult {
        best: initial,
        best_score: score,
        initial_score: score,
        steps: 0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::localsearch::{Metaheuristic, Termination};
    use crate::model::{HardSoftScore, PlanningEntity};

    /// construct.rs 测试同文的内嵌算例（经 format::parse 全链路）。
    /// 已知构造分：FF → [0,1,2]（hard 0, soft 8）；BF → [0,2,2]（hard 0, soft 3）。
    const HAND_INSTANCE: &str = "\
NAME hand
VALUES 3
ENTITIES 3
COSTROW 0 2 9 4
COSTROW 1 0 0 0
COST 1 0 7
COST 1 1 6
COST 1 2 1
FORBID 2 1
CONF 0 2
CAP 0 1
DEMAND 0 1
DEMAND 1 1
DEMAND 2 1
";

    fn hand_model() -> Arc<ProblemModel> {
        Arc::new(crate::format::parse(HAND_INSTANCE).expect("parse hand instance"))
    }

    #[test]
    fn solver_config_defaults_are_canonical() {
        let cfg = SolverConfig::default();
        assert_eq!(cfg.construction, ConstructionHeuristic::FirstFit);
        // LocalSearchConfig 无 PartialEq（契约 derive 仅 Debug/Clone）→ 逐字段比对
        let ls = LocalSearchConfig::default();
        assert_eq!(cfg.local_search.metaheuristic, ls.metaheuristic);
        assert_eq!(cfg.local_search.seed, ls.seed);
        assert_eq!(cfg.local_search.termination, ls.termination);
        assert_eq!(cfg.local_search.move_kind, ls.move_kind);
        assert_eq!(cfg.local_search.samples_per_step, ls.samples_per_step);
        assert_eq!(cfg.local_search.metaheuristic, Metaheuristic::TabuSearch { tabu_size: 8 });
        assert!(cfg.run_local_search);
    }

    /// 编排默认路径：validate → FF 构造（soft 8）→ tabu 局部搜索（不劣于构造分）。
    #[test]
    fn solve_runs_construction_then_local_search() {
        let model = hand_model();
        let cfg = SolverConfig {
            local_search: LocalSearchConfig {
                metaheuristic: Metaheuristic::TabuSearch { tabu_size: 8 },
                termination: Termination::Steps(3000),
                ..LocalSearchConfig::default()
            },
            ..SolverConfig::default()
        };
        let r = solve(model.clone(), &cfg).expect("solve");
        assert_eq!(r.initial_score, HardSoftScore::new(0, 8), "FF construction score");
        assert!(r.best_score <= r.initial_score);
        assert!(r.best.is_complete());
        assert!(r.steps > 0);
        let calc = ScoreCalculator::new(model.clone());
        assert_eq!(r.best_score, calc.total_score(&r.best));
        for e in 0..model.entity_count() {
            assert!(model.in_domain(e, r.best.value(e)));
        }
    }

    /// run_local_search = false：仅构造，best == 构造产物、steps = 0。
    #[test]
    fn solve_without_local_search_returns_construction() {
        let model = hand_model();
        let cfg = SolverConfig { run_local_search: false, ..SolverConfig::default() };
        let r = solve(model, &cfg).expect("solve ff-only");
        assert_eq!(r.steps, 0);
        assert_eq!(r.initial_score, HardSoftScore::new(0, 8));
        assert_eq!(r.best_score, r.initial_score);
        assert_eq!(r.best.assignment(), &[0, 1, 2]);

        // BestFit 构造分同样固定
        let model = hand_model();
        let cfg = SolverConfig {
            construction: ConstructionHeuristic::BestFit,
            run_local_search: false,
            ..SolverConfig::default()
        };
        let r = solve(model, &cfg).expect("solve bf-only");
        assert_eq!(r.best.assignment(), &[0, 2, 2]);
        assert_eq!(r.best_score, HardSoftScore::new(0, 3));
    }

    /// 模型非法（空域）：validate 的 Err 经 solve 透传（Model，而非 panic）。
    #[test]
    fn solve_invalid_model_returns_model_error() {
        let model = Arc::new(ProblemModel {
            name: "bad".to_string(),
            value_count: 2,
            resources: vec![],
            entities: vec![PlanningEntity {
                domain: vec![],
                demands: vec![],
                conflicts: vec![],
            }],
            cost: vec![vec![0, 0]],
        });
        match solve(model, &SolverConfig::default()) {
            Err(PlannerError::Model(_)) => {}
            other => panic!("expected Model error, got {other:?}"),
        }
    }

    /// 同模型同配置整链路（validate → 构造 → 局搜）逐位可复现。
    #[test]
    fn solve_is_deterministic_end_to_end() {
        let model = hand_model();
        let cfg = SolverConfig {
            local_search: LocalSearchConfig {
                metaheuristic: Metaheuristic::LateAcceptance { history_len: 30 },
                termination: Termination::Steps(2000),
                ..LocalSearchConfig::default()
            },
            ..SolverConfig::default()
        };
        let a = solve(model.clone(), &cfg).expect("solve a");
        let b = solve(model, &cfg).expect("solve b");
        assert_eq!(a.best.assignment(), b.best.assignment());
        assert_eq!(a.best_score, b.best_score);
        assert_eq!(a.steps, b.steps);
    }

    /// sa 配置经编排同样可用且不劣于构造分。
    #[test]
    fn solve_with_sa_config() {
        let model = hand_model();
        let cfg = SolverConfig {
            local_search: LocalSearchConfig {
                metaheuristic: Metaheuristic::SimulatedAnnealing {
                    start_temp: 1000.0,
                    end_temp: 1.0,
                },
                termination: Termination::Steps(2000),
                ..LocalSearchConfig::default()
            },
            ..SolverConfig::default()
        };
        let r = solve(model, &cfg).expect("solve sa");
        assert!(r.best_score <= r.initial_score);
        assert!(r.best.is_complete());
    }
}
