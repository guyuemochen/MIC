//! 局部搜索：禁忌搜索 / 模拟退火 / 迟到接受 + 终止条件与搜索结果。
//!
//! ## 算法规则（normative）
//! - 通用：入口 `IncrementalScorer::reset(initial)`（debug_assert 与全量重算相等）；
//!   维护同步对 (sol, inc) 与历史最优快照（初始 = 初始解分数与赋值）；比较一律
//!   [`HardSoftScore`] 字典序；每步先耗随机数生成移动、后判定。
//! - **TabuSearch**：每步采样 `samples_per_step` 个 [`crate::moves::random_move`]
//!   （None 样本跳过但随机数照常消耗）；可采纳 = 涉及实体全部满足
//!   `step ≥ tabu_until[e]`，或预演分**严格** < 历史最优（特赦）；可采纳样本中取预演分
//!   最小者（平局取采样序最前，严格小于才替换）；无可采纳样本则空步。应用后对涉及实体
//!   置 `tabu_until[e] = step + tabu_size`；预演分 < best 时更新快照。
//! - **SimulatedAnnealing**：每步一个 random_move；标量 `W = 1 + Σ_{e,v} cost`
//!   （validate 保证 ≤ 2^40+1）；`s = hard as i128 * W as i128 + soft as i128`，
//!   `Δ = s(preview) − s(当前)`（**全程 i128，禁 i64 乘法**——实测上界积
//!   2^46·(2^40+1) = 77371252455406635925372928 已超 i64）；`Δ ≤ 0` 接受；否则
//!   `p = exp(−(Δ as f64) / T_k)`，`gen_ratio(p)` 为 true 接受；
//!   `T_k = max(end_temp, start_temp · 0.999^k)`；best 仍按真实 HardSoftScore 维护。
//! - **LateAcceptance**：每步一个 random_move；预演分 ≤ `ring[step mod history_len]`
//!   则应用；无论接受与否，随后 `ring[step mod history_len] =` 当前（决策后）分数；
//!   ring 初始化全为初始解分数。

use std::sync::Arc;
use std::time::{Duration, Instant};

use crate::model::{HardSoftScore, PlannerError, ProblemModel, Solution};
use crate::moves::{random_move, Move, MoveKind};
use crate::rng::Rng;
use crate::score::IncrementalScorer;

/// 终止条件。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Termination {
    /// 恰好执行 N 步（含被拒绝步 / 空步）——比特级可复现；
    /// `Steps(0)` 合法 = 直接返回初始解。
    Steps(usize),
    /// 每步开始前检查已耗时；结果是"已完成步"的确定函数
    /// （墙钟截断点随机器速度不同，跨机器不保证相同）。
    TimeLimit(Duration),
}

impl Default for Termination {
    /// 默认 `Steps(50_000)`。
    fn default() -> Self {
        Termination::Steps(50_000)
    }
}

/// 元启发式（参数默认值全文唯一定义，CLI/README/测试均引用此处，不得另造数值）。
#[derive(Debug, Clone, PartialEq)]
pub enum Metaheuristic {
    /// 实体禁忌：被移动实体在 `tabu_size` 步内不可再动（特赦除外）。
    /// CLI/测试默认 `tabu_size = 8`；须 ≥ 1。
    TabuSearch { tabu_size: usize },
    /// 模拟退火。CLI 默认 `{ start_temp: 1000.0, end_temp: 1.0 }`；
    /// 须 `start_temp > 0.0`、`end_temp > 0.0`、`end_temp ≤ start_temp`。
    SimulatedAnnealing { start_temp: f64, end_temp: f64 },
    /// 迟到接受。CLI/测试默认 `history_len = 30`；须 ≥ 1
    /// （0 会使 `ring[step % 0]` 除零 panic，入口校验拦截）。
    LateAcceptance { history_len: usize },
}

/// 局部搜索配置（默认值全文唯一来源）。
#[derive(Debug, Clone)]
pub struct LocalSearchConfig {
    /// 默认 `TabuSearch { tabu_size: 8 }`。
    pub metaheuristic: Metaheuristic,
    /// 随机种子；默认 42。
    pub seed: u64,
    /// 终止条件；默认 `Steps(50_000)`。
    pub termination: Termination,
    /// 移动配比；默认 `Mixed { change_ratio: 0.5 }`。
    pub move_kind: MoveKind,
    /// 每步采样数；默认 32；**仅 TabuSearch 读取**（sa/la 完全忽略，
    /// 不因它多耗一个随机数）。
    pub samples_per_step: usize,
}

impl Default for LocalSearchConfig {
    /// 默认：TabuSearch{tabu_size:8} / seed 42 / Steps(50_000) / Mixed{0.5} / samples 32。
    fn default() -> Self {
        LocalSearchConfig {
            metaheuristic: Metaheuristic::TabuSearch { tabu_size: 8 },
            seed: 42,
            termination: Termination::Steps(50_000),
            move_kind: MoveKind::default(),
            samples_per_step: 32,
        }
    }
}

/// 搜索结果。
#[derive(Debug, Clone)]
pub struct SearchResult {
    /// 历史最优快照（完整解）。
    pub best: Solution,
    /// 历史最优分数；不变式 `best_score ≤ initial_score`（HardSoftScore 字典序）。
    pub best_score: HardSoftScore,
    /// 初始解（构造产物）分数。
    pub initial_score: HardSoftScore,
    /// 实际执行步数。
    pub steps: usize,
}

/// SA 标量 `s = hard·W + soft`：全程 i128（禁 i64 乘法）。
/// validate 界保证 hard ≤ 2^46、W ≤ 2^40+1、soft ≤ 2^40，故
/// s ≤ 2^46·(2^40+1) + 2^40 = 77371252455407735437000704 > i64::MAX（实测
/// i64::try_from 均失败）；仅 Δ/T 转 f64 进 exp。
#[inline]
fn scalar(score: HardSoftScore, w: i128) -> i128 {
    i128::from(score.hard) * w + i128::from(score.soft)
}

/// 局部搜索主循环。
///
/// 入口校验（全部 `Err(PlannerError::Invalid(..))`，杜绝 panic 类配置）：
/// `initial.is_complete()`；`tabu_size ≥ 1`；`samples_per_step ≥ 1`；
/// `history_len ≥ 1`（LA）；`start_temp > 0.0 && end_temp > 0.0 && end_temp ≤ start_temp`（SA）；
/// `Mixed { r }` 须 `0.0 ≤ r ≤ 1.0`。
/// 不变式：返回 `best_score ≤ initial_score`（字典序），`best` 为完整解。
pub fn local_search(
    initial: Solution,
    cfg: &LocalSearchConfig,
) -> Result<SearchResult, PlannerError> {
    // ---- 入口校验（先于任何随机数消耗 / 状态构建） ----
    if !initial.is_complete() {
        return Err(PlannerError::Invalid(
            "initial solution is not complete".to_string(),
        ));
    }
    match &cfg.metaheuristic {
        Metaheuristic::TabuSearch { tabu_size } if *tabu_size < 1 => {
            return Err(PlannerError::Invalid(format!(
                "tabu_size must be >= 1 (got {tabu_size})"
            )));
        }
        Metaheuristic::SimulatedAnnealing { start_temp, end_temp }
            if !(*start_temp > 0.0 && *end_temp > 0.0 && *end_temp <= *start_temp) =>
        {
            return Err(PlannerError::Invalid(format!(
                "SA needs start_temp > 0, end_temp > 0 and end_temp <= start_temp \
                 (got start_temp={start_temp}, end_temp={end_temp})"
            )));
        }
        Metaheuristic::LateAcceptance { history_len } if *history_len < 1 => {
            return Err(PlannerError::Invalid(format!(
                "history_len must be >= 1 (got {history_len})"
            )));
        }
        _ => {}
    }
    if cfg.samples_per_step < 1 {
        return Err(PlannerError::Invalid(format!(
            "samples_per_step must be >= 1 (got {})",
            cfg.samples_per_step
        )));
    }
    if let MoveKind::Mixed { change_ratio } = cfg.move_kind {
        // NaN 在此一并拒绝（NaN 与任何比较均为 false）
        if !(change_ratio >= 0.0 && change_ratio <= 1.0) {
            return Err(PlannerError::Invalid(format!(
                "change_ratio must be within [0.0, 1.0] (got {change_ratio})"
            )));
        }
    }

    // ---- 同步对与历史最优初始化 ----
    // Solution 经 model_arc() 直接共享其 Arc<ProblemModel>（零拷贝；模型不可变，
    // 内容与深拷贝逐位相等 ⇒ 搜索结果与随机流不变，且多起点并行下各 worker
    // 可安全复用同一模型，消除每线程一次 O(n·m) 深拷贝）。
    let model: Arc<ProblemModel> = initial.model_arc().clone();
    let n = model.entity_count();
    let mut rng = Rng::new(cfg.seed);
    let mut sol = initial;
    let mut inc = IncrementalScorer::new(model.clone());
    let initial_score = inc.reset(&sol);
    #[cfg(debug_assertions)]
    {
        let calc = crate::score::ScoreCalculator::new(model.clone());
        debug_assert_eq!(
            initial_score,
            calc.total_score(&sol),
            "reset must equal full recompute"
        );
    }
    let mut best = sol.clone();
    let mut best_score = initial_score;
    let mut steps = 0usize;

    // 每步开始前的终止判定：Steps 恰 N 步；TimeLimit 已耗时 ≥ limit 即停。
    // TimeLimit 的墙钟截断点因机器而异，但结果仍是"已完成步"的确定函数。
    let t0 = Instant::now();
    let running = |steps: usize| match cfg.termination {
        Termination::Steps(max) => steps < max,
        Termination::TimeLimit(limit) => t0.elapsed() < limit,
    };

    match &cfg.metaheuristic {
        Metaheuristic::TabuSearch { tabu_size } => {
            let tabu_size = *tabu_size;
            let mut tabu_until = vec![0usize; n];
            while running(steps) {
                // 采样：None 样本跳过（其随机数已在 random_move 内照常消耗）
                let mut chosen: Option<(HardSoftScore, Move)> = None;
                for _ in 0..cfg.samples_per_step {
                    let Some(mv) = random_move(&sol, &inc, &mut rng, cfg.move_kind) else {
                        continue;
                    };
                    let preview = mv.preview(&inc);
                    let (ents, k) = mv.entities();
                    // 判定顺序固定：先禁忌后特赦（特赦为严格 <）
                    let tabu_free = (0..k).all(|i| steps >= tabu_until[ents[i]]);
                    let asylum = preview < best_score;
                    if tabu_free || asylum {
                        // 平局取采样序最前：严格小于才替换
                        if chosen.map_or(true, |(bp, _)| preview < bp) {
                            chosen = Some((preview, mv));
                        }
                    }
                }
                if let Some((preview, mv)) = chosen {
                    mv.apply(&mut sol, &mut inc);
                    let (ents, k) = mv.entities();
                    for i in 0..k {
                        tabu_until[ents[i]] = steps + tabu_size;
                    }
                    if preview < best_score {
                        best_score = preview;
                        best = sol.clone();
                    }
                }
                // 无可采纳样本：空步（sol/best 均不变，步数照计）
                steps += 1;
            }
        }
        Metaheuristic::SimulatedAnnealing { start_temp, end_temp } => {
            let (start_temp, end_temp) = (*start_temp, *end_temp);
            // W = 1 + Σ_{e,v} cost（validate 保证 ≤ 2^40+1；i128 求和无窄化）
            let w: i128 = 1
                + model
                    .cost
                    .iter()
                    .flat_map(|row| row.iter())
                    .map(|&c| i128::from(c))
                    .sum::<i128>();
            while running(steps) {
                let t_k = end_temp.max(start_temp * 0.999f64.powf(steps as f64));
                if let Some(mv) = random_move(&sol, &inc, &mut rng, cfg.move_kind) {
                    let preview = mv.preview(&inc);
                    // Δ 全程 i128；Δ ≤ 0 判定不过 f64
                    let delta = scalar(preview, w) - scalar(inc.score(), w);
                    let accept = delta <= 0 || {
                        let p = (-(delta as f64) / t_k).exp();
                        rng.gen_ratio(p)
                    };
                    if accept {
                        mv.apply(&mut sol, &mut inc);
                        // best 按真实 HardSoftScore 维护（已访问解的最优快照）
                        if preview < best_score {
                            best_score = preview;
                            best = sol.clone();
                        }
                    }
                }
                steps += 1;
            }
        }
        Metaheuristic::LateAcceptance { history_len } => {
            let history_len = *history_len;
            let mut ring = vec![initial_score; history_len];
            while running(steps) {
                let idx = steps % history_len;
                if let Some(mv) = random_move(&sol, &inc, &mut rng, cfg.move_kind) {
                    let preview = mv.preview(&inc);
                    if preview <= ring[idx] {
                        mv.apply(&mut sol, &mut inc);
                        if preview < best_score {
                            best_score = preview;
                            best = sol.clone();
                        }
                    }
                }
                // 无论接受与否（含 None 步）：写回决策后当前分数
                ring[idx] = inc.score();
                steps += 1;
            }
        }
    }

    debug_assert!(
        best_score <= initial_score,
        "best_score must never be worse than initial_score"
    );
    debug_assert!(best.is_complete(), "best snapshot must stay complete");
    Ok(SearchResult {
        best,
        best_score,
        initial_score,
        steps,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::construct::{construct, ConstructionHeuristic};
    use crate::model::{PlanningEntity, ProblemModel, Resource};
    use crate::score::ScoreCalculator;

    /// 与 score.rs / moves.rs 测试同构的玩具模型（m=4、2 资源、6 实体、冲突链 + 三角）。
    fn toy_model() -> Arc<ProblemModel> {
        let cap_sentinel = i64::MAX / 2;
        let model = ProblemModel {
            name: "ls-toy".to_string(),
            value_count: 4,
            resources: vec![
                Resource { name: "cpu".to_string(), capacity: vec![2, 2, 3, cap_sentinel] },
                Resource { name: "mem".to_string(), capacity: vec![3, 1, 2, 2] },
            ],
            entities: vec![
                PlanningEntity { domain: vec![0, 1, 2, 3], demands: vec![2, 1], conflicts: vec![1, 2] },
                PlanningEntity { domain: vec![0, 1, 3], demands: vec![1, 2], conflicts: vec![0, 3] },
                PlanningEntity { domain: vec![1, 2, 3], demands: vec![3, 1], conflicts: vec![0, 3, 4] },
                PlanningEntity { domain: vec![0, 2, 3], demands: vec![1, 1], conflicts: vec![1, 2] },
                PlanningEntity { domain: vec![0, 1, 2, 3], demands: vec![2, 3], conflicts: vec![2, 5] },
                PlanningEntity { domain: vec![0, 2, 3], demands: vec![1, 1], conflicts: vec![4] },
            ],
            cost: vec![
                vec![4, 1, 7, 0],
                vec![2, 9, 0, 3],
                vec![0, 5, 6, 2],
                vec![8, 0, 1, 4],
                vec![3, 2, 8, 5],
                vec![6, 0, 4, 7],
            ],
        };
        model.validate().expect("toy model must be valid");
        Arc::new(model)
    }

    /// 单实体二值模型：cost = [c0, c1]，无资源无冲突（n=1 ⇒ Swap 结构性不可能）。
    fn single_entity_model(c0: i64, c1: i64) -> Arc<ProblemModel> {
        let model = ProblemModel {
            name: "single-entity".to_string(),
            value_count: 2,
            resources: vec![],
            entities: vec![PlanningEntity {
                domain: vec![0, 1],
                demands: vec![],
                conflicts: vec![],
            }],
            cost: vec![vec![c0, c1]],
        };
        model.validate().expect("single entity model must be valid");
        Arc::new(model)
    }

    /// SA 标量 i128 宽度模型：hard = 2^32、W = 2^40+1 ⇒ s ≈ 2^72 远超 i64::MAX。
    /// 误用 i64 乘法会溢出（debug 下 panic）。
    fn wide_scalar_model() -> Arc<ProblemModel> {
        let model = ProblemModel {
            name: "wide".to_string(),
            value_count: 2,
            resources: vec![Resource { name: "r".to_string(), capacity: vec![0, 1 << 40] }],
            entities: vec![
                PlanningEntity {
                    domain: vec![0, 1],
                    demands: vec![1 << 31],
                    conflicts: vec![],
                },
                PlanningEntity {
                    domain: vec![0, 1],
                    demands: vec![1 << 31],
                    conflicts: vec![],
                },
            ],
            cost: vec![vec![1 << 38, 1 << 38], vec![1 << 38, 1 << 38]],
        };
        model.validate().expect("wide model within hard bounds");
        Arc::new(model)
    }

    /// 同 seed 两次运行逐位一致 + 返回结果自洽（完整、域内、分数与全量重算相等）。
    fn assert_deterministic_and_sound(
        model: &Arc<ProblemModel>,
        initial: &Solution,
        cfg: &LocalSearchConfig,
    ) {
        let a = local_search(initial.clone(), cfg).expect("run a");
        let b = local_search(initial.clone(), cfg).expect("run b");
        assert_eq!(
            a.best.assignment(),
            b.best.assignment(),
            "same seed must reproduce best assignment bit for bit"
        );
        assert_eq!(a.best_score, b.best_score);
        assert_eq!(a.initial_score, b.initial_score);
        assert_eq!(a.steps, b.steps);
        assert!(a.best.is_complete(), "best must be complete");
        for e in 0..model.entity_count() {
            assert!(
                model.in_domain(e, a.best.value(e)),
                "best value of entity {e} must be in domain"
            );
        }
        let calc = ScoreCalculator::new(model.clone());
        assert_eq!(a.best_score, calc.total_score(&a.best));
        assert_eq!(a.initial_score, calc.total_score(initial));
        assert!(a.best_score <= a.initial_score, "best must be <= initial");
    }

    #[test]
    fn tabu_same_seed_reproduces_bit_for_bit() {
        let model = toy_model();
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let cfg = LocalSearchConfig {
            metaheuristic: Metaheuristic::TabuSearch { tabu_size: 8 },
            termination: Termination::Steps(4000),
            ..LocalSearchConfig::default()
        };
        assert_deterministic_and_sound(&model, &initial, &cfg);
    }

    #[test]
    fn la_same_seed_reproduces_bit_for_bit() {
        let model = toy_model();
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let cfg = LocalSearchConfig {
            metaheuristic: Metaheuristic::LateAcceptance { history_len: 30 },
            termination: Termination::Steps(4000),
            ..LocalSearchConfig::default()
        };
        assert_deterministic_and_sound(&model, &initial, &cfg);
    }

    #[test]
    fn sa_same_seed_reproduces_bit_for_bit() {
        let model = toy_model();
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let cfg = LocalSearchConfig {
            metaheuristic: Metaheuristic::SimulatedAnnealing {
                start_temp: 1000.0,
                end_temp: 1.0,
            },
            termination: Termination::Steps(4000),
            ..LocalSearchConfig::default()
        };
        assert_deterministic_and_sound(&model, &initial, &cfg);
    }

    /// 入口校验：非法配置一律 `Err(PlannerError::Invalid(..))`，绝不 panic。
    #[test]
    fn invalid_configurations_return_invalid() {
        let model = toy_model();
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let cases: Vec<LocalSearchConfig> = vec![
            LocalSearchConfig {
                metaheuristic: Metaheuristic::TabuSearch { tabu_size: 0 },
                ..LocalSearchConfig::default()
            },
            LocalSearchConfig {
                metaheuristic: Metaheuristic::LateAcceptance { history_len: 0 },
                ..LocalSearchConfig::default()
            },
            // end_temp > start_temp / start_temp ≤ 0 / end_temp ≤ 0
            LocalSearchConfig {
                metaheuristic: Metaheuristic::SimulatedAnnealing {
                    start_temp: 1000.0,
                    end_temp: 1001.0,
                },
                ..LocalSearchConfig::default()
            },
            LocalSearchConfig {
                metaheuristic: Metaheuristic::SimulatedAnnealing {
                    start_temp: 0.0,
                    end_temp: 1.0,
                },
                ..LocalSearchConfig::default()
            },
            LocalSearchConfig {
                metaheuristic: Metaheuristic::SimulatedAnnealing {
                    start_temp: 5.0,
                    end_temp: -1.0,
                },
                ..LocalSearchConfig::default()
            },
            // samples 校验对三种元启发式一律生效（sa/la 不读取它也要 ≥ 1）
            LocalSearchConfig {
                samples_per_step: 0,
                ..LocalSearchConfig::default()
            },
            LocalSearchConfig {
                metaheuristic: Metaheuristic::LateAcceptance { history_len: 30 },
                samples_per_step: 0,
                ..LocalSearchConfig::default()
            },
            LocalSearchConfig {
                metaheuristic: Metaheuristic::SimulatedAnnealing {
                    start_temp: 1000.0,
                    end_temp: 1.0,
                },
                samples_per_step: 0,
                ..LocalSearchConfig::default()
            },
            // ratio 越界 / NaN
            LocalSearchConfig {
                move_kind: MoveKind::Mixed { change_ratio: 1.5 },
                ..LocalSearchConfig::default()
            },
            LocalSearchConfig {
                move_kind: MoveKind::Mixed { change_ratio: -0.1 },
                ..LocalSearchConfig::default()
            },
            LocalSearchConfig {
                move_kind: MoveKind::Mixed { change_ratio: f64::NAN },
                ..LocalSearchConfig::default()
            },
        ];
        for cfg in &cases {
            match local_search(initial.clone(), cfg) {
                Err(PlannerError::Invalid(_)) => {}
                other => panic!("expected Invalid for cfg {cfg:?}, got {other:?}"),
            }
        }
        // 不完整初始解（含全未分配空解）
        let incomplete = Solution::empty(model.clone());
        match local_search(incomplete, &LocalSearchConfig::default()) {
            Err(PlannerError::Invalid(_)) => {}
            other => panic!("expected Invalid for incomplete initial, got {other:?}"),
        }
    }

    /// Steps(0)：合法，直接返回初始解，0 步。
    #[test]
    fn steps_zero_returns_initial_unchanged() {
        let model = toy_model();
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let cfg = LocalSearchConfig {
            termination: Termination::Steps(0),
            ..LocalSearchConfig::default()
        };
        let r = local_search(initial.clone(), &cfg).expect("zero steps");
        assert_eq!(r.steps, 0);
        assert_eq!(r.best.assignment(), initial.assignment());
        assert_eq!(r.best_score, r.initial_score);
    }

    /// TimeLimit(0)：首步开始前已超时 ⇒ 0 步返回初始解（确定行为）。
    #[test]
    fn time_limit_zero_completes_zero_steps() {
        let model = toy_model();
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let cfg = LocalSearchConfig {
            termination: Termination::TimeLimit(Duration::ZERO),
            ..LocalSearchConfig::default()
        };
        let r = local_search(initial.clone(), &cfg).expect("zero time");
        assert_eq!(r.steps, 0);
        assert_eq!(r.best.assignment(), initial.assignment());
        assert_eq!(r.best_score, r.initial_score);
    }

    /// TimeLimit(非零) 冒烟：步数不确定，但结果自洽（best 完整、域内、分数一致、不劣于初始）。
    #[test]
    fn time_limit_positive_returns_sound_result() {
        let model = toy_model();
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let cfg = LocalSearchConfig {
            termination: Termination::TimeLimit(Duration::from_millis(5)),
            ..LocalSearchConfig::default()
        };
        let r = local_search(initial, &cfg).expect("bounded run");
        assert!(r.best.is_complete());
        assert!(r.best_score <= r.initial_score);
        let calc = ScoreCalculator::new(model.clone());
        assert_eq!(r.best_score, calc.total_score(&r.best));
    }

    /// 全部移动被拒的人造场景：单实体、初始已最优（soft 0）、唯一非平凡移动 soft 5。
    /// LA：preview (0,5) ≤ ring (0,0) 恒假 ⇒ 全部拒绝；tabu：唯一移动可被应用
    /// （禁忌允许恶化、特赦不触发），但 best 快照守住初始最优。
    #[test]
    fn best_le_initial_when_every_move_is_rejected() {
        let model = single_entity_model(0, 5); // 值 0 代价 0（最优）、值 1 代价 5
        let initial = Solution::from_assignment(model.clone(), &[0]).expect("in domain");
        assert_eq!(
            ScoreCalculator::new(model.clone()).total_score(&initial),
            HardSoftScore::new(0, 0)
        );

        let la = LocalSearchConfig {
            metaheuristic: Metaheuristic::LateAcceptance { history_len: 30 },
            termination: Termination::Steps(600),
            ..LocalSearchConfig::default()
        };
        let r = local_search(initial.clone(), &la).expect("la");
        assert_eq!(r.steps, 600, "rejected steps still count as steps");
        assert_eq!(r.best_score, r.initial_score);
        assert_eq!(r.best_score, HardSoftScore::new(0, 0));
        assert_eq!(r.best.assignment(), &[0]);

        let tabu = LocalSearchConfig {
            metaheuristic: Metaheuristic::TabuSearch { tabu_size: 8 },
            termination: Termination::Steps(600),
            ..LocalSearchConfig::default()
        };
        let r = local_search(initial, &tabu).expect("tabu");
        assert_eq!(r.steps, 600);
        assert_eq!(
            r.best_score,
            HardSoftScore::new(0, 0),
            "tabu may wander but best snapshot holds the initial optimum"
        );
        assert_eq!(r.best.assignment(), &[0]);
    }

    /// FF 起点非最优时，tabu 与 la 都能找到更便宜的值（确定性小算例）。
    #[test]
    fn tabu_and_la_improve_first_fit_construction() {
        let model = single_entity_model(5, 1); // FF 起点值 0（soft 5），最优值 1（soft 1）
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        assert_eq!(initial.assignment(), &[0]);
        for meta in [
            Metaheuristic::TabuSearch { tabu_size: 8 },
            Metaheuristic::LateAcceptance { history_len: 30 },
        ] {
            let cfg = LocalSearchConfig {
                metaheuristic: meta.clone(),
                termination: Termination::Steps(200),
                ..LocalSearchConfig::default()
            };
            let r = local_search(initial.clone(), &cfg).expect("search");
            assert_eq!(r.initial_score, HardSoftScore::new(0, 5), "FF start is soft 5");
            assert_eq!(
                r.best_score,
                HardSoftScore::new(0, 1),
                "{meta:?}: must find the cheaper value"
            );
            assert_eq!(r.best.assignment(), &[1]);
        }
    }

    /// 三个元启发式 × FF/BF 构造：best ≤ initial、best 完整且域内、分数与全量重算一致。
    #[test]
    fn best_never_worse_across_algorithms_and_heuristics() {
        let model = toy_model();
        let metas = [
            Metaheuristic::TabuSearch { tabu_size: 8 },
            Metaheuristic::SimulatedAnnealing {
                start_temp: 1000.0,
                end_temp: 1.0,
            },
            Metaheuristic::LateAcceptance { history_len: 30 },
        ];
        for h in [ConstructionHeuristic::FirstFit, ConstructionHeuristic::BestFit] {
            let initial = construct(&model, h).expect("construct");
            for meta in &metas {
                let cfg = LocalSearchConfig {
                    metaheuristic: meta.clone(),
                    termination: Termination::Steps(1500),
                    ..LocalSearchConfig::default()
                };
                let r = local_search(initial.clone(), &cfg).expect("search");
                assert!(
                    r.best_score <= r.initial_score,
                    "{h:?}/{meta:?}: best must be <= initial"
                );
                assert_eq!(r.steps, 1500);
                assert!(r.best.is_complete());
                for e in 0..model.entity_count() {
                    assert!(model.in_domain(e, r.best.value(e)));
                }
                let calc = ScoreCalculator::new(model.clone());
                assert_eq!(r.best_score, calc.total_score(&r.best));
            }
        }
    }

    /// SA 标量必须走 i128：hard=2^32、W=2^40+1 ⇒ s≈2^72 > i64::MAX。
    /// 若实现误用 i64 乘法，debug 构建直接 overflow panic。
    #[test]
    fn sa_survives_scalar_magnitudes_beyond_i64() {
        let model = wide_scalar_model();
        // 初始 [0,0]：load[0] = 2·2^31 = 2^32 > cap[0] = 0 ⇒ hard = 2^32；soft = 2^39
        let initial = Solution::from_assignment(model.clone(), &[0, 0]).expect("in domain");
        let cfg = LocalSearchConfig {
            metaheuristic: Metaheuristic::SimulatedAnnealing {
                start_temp: 1000.0,
                end_temp: 1.0,
            },
            termination: Termination::Steps(2000),
            ..LocalSearchConfig::default()
        };
        let r = local_search(initial, &cfg).expect("sa on wide scalar model");
        assert_eq!(r.initial_score, HardSoftScore::new(1 << 32, 2 * (1 << 38)));
        assert!(r.best_score <= r.initial_score);
        let calc = ScoreCalculator::new(model.clone());
        assert_eq!(r.best_score, calc.total_score(&r.best));
    }

    #[test]
    fn defaults_are_canonical() {
        assert_eq!(Termination::default(), Termination::Steps(50_000));
        let c = LocalSearchConfig::default();
        assert_eq!(c.metaheuristic, Metaheuristic::TabuSearch { tabu_size: 8 });
        assert_eq!(c.seed, 42);
        assert_eq!(c.termination, Termination::Steps(50_000));
        assert_eq!(c.move_kind, MoveKind::Mixed { change_ratio: 0.5 });
        assert_eq!(c.samples_per_step, 32);
    }

    /// Steps 终止的可复现承诺：不同种子允许不同结果，但重跑同种子必须一致
    /// （覆盖非默认 seed）。
    #[test]
    fn non_default_seed_is_reproducible_too() {
        let model = toy_model();
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let cfg = LocalSearchConfig {
            seed: 2026,
            termination: Termination::Steps(1000),
            ..LocalSearchConfig::default()
        };
        assert_deterministic_and_sound(&model, &initial, &cfg);
    }
}
