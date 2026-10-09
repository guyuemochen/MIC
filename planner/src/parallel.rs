//! 多起点并行编排（OptaPlanner multistart 模式，std-only scoped threads）。
//!
//! ## 并发模型（零共享可变状态）
//! - [`std::thread::scope`] + [`std::thread::Builder::spawn_scoped`]（rustc 1.63 稳定）：
//!   已 spawn 线程在 scope 退出时自动 join（含闭包提前 return 后被 drop 的 handle），
//!   无泄漏、无死锁——worker 只做纯计算，不等待父线程或彼此。
//! - 线程间**零共享可变状态**（无 Mutex/RwLock/Atomic/channel）：唯一共享对象是
//!   只读 `Arc<ProblemModel>`（父线程构造一次初始解，各 worker 经
//!   [`crate::model::Solution::model_arc`] 零拷贝复用同一不可变模型）；配置在派生前
//!   逐 worker 克隆、仅替换 seed 后按所有权移入。每 worker 私有持有 Solution 克隆、
//!   [`crate::score::IncrementalScorer`] 与 [`crate::rng::Rng`]。
//! - worker panic 不吞不锁死：handle 按下标序 join、记住最小下标 payload、仍 join 完
//!   其余线程，然后在 scope **外** [`std::panic::resume_unwind`]；`spawn_scoped` 的
//!   io 失败（如资源不足）映射 [`PlannerError::Io`](crate::model::PlannerError::Io)
//!   并停止继续 spawn——返回 Err 前先按下标序 join 已派生的全部线程，其中的 panic
//!   payload 照常暂存并于 scope 外重抛（「panic 不吞」不因 spawn 失败破例）。
//!
//! ## 确定性（Steps 终止下逐位可复现）
//! 固定 seed + 固定 threads + [`crate::localsearch::Termination::Steps`] 终止 ⇒
//! 聚合结果逐位可复现：
//! 1. worker 种子集合确定：worker 0 恒用 `cfg.local_search.seed` 原值；worker i ≥ 1 用
//!    [`crate::rng::derive_seed`](seed, i)，固定 seed 时派生种子两两互异
//!    （i ↦ seed_i 在 i ≥ 1 上单射；与 worker 0 原值的碰撞概率 ~threads/2^64 ≤ 2^-54，
//!    即便命中也仅两个 worker 重复同一搜索，min 聚合与确定性不受影响——见
//!    [`crate::rng::derive_seed`]）；
//! 2. [`crate::localsearch::local_search`] 是 (initial, cfg) 的纯函数（initial 为父线程
//!    一次性构造的同一解）；
//! 3. 聚合 = min(best_score)，**严格小于才替换** ⇒ 并列取最小 worker 下标，且聚合在
//!    全部 join 完成后按下标序进行 ⇒ 与线程完成顺序无关。
//! 推论：worker 0 = 主 seed ⇒ 任何 threads=N 的运行完整包含 threads=1 的那次搜索，
//! 聚合取 min ⇒ Steps 终止下 threads=N 结果恒 ≤ threads=1（HardSoftScore 字典序）
//! ——"不劣于"是结构性质而非概率承诺。
//!
//! ## TimeLimit
//! 每 worker 收到同一 [`std::time::Duration`]（零共享状态方案，不设共享 deadline），
//! 各自从进入搜索时计时；墙钟截断本就不承诺跨机一致
//! （[`crate::localsearch::Termination::TimeLimit`] 既定语义），此时"不劣于"为尽力
//! 而非结构保证。

use std::any::Any;
use std::sync::Arc;

use crate::construct::construct;
use crate::localsearch::{local_search, LocalSearchConfig, SearchResult};
use crate::model::{HardSoftScore, PlannerError, ProblemModel, Solution};
use crate::rng;
use crate::solver::{solve, SolverConfig};

/// threads 合法上界（合法域 [1, 1024]）。
pub const MAX_THREADS: usize = 1024;

/// 多起点聚合报告。winner/worker_scores 仅描述性，聚合语义见
/// [`solve_parallel_detailed`]。
#[derive(Debug, Clone)]
pub struct ParallelReport {
    /// 胜出 worker 的 [`SearchResult`]（原样返回，不变式随胜者保持：
    /// `best_score ≤ initial_score`、`best` 完整）。
    pub best: SearchResult,
    /// 胜出 worker 下标（best_score 并列取最小下标，与完成顺序无关）。
    pub winner: usize,
    /// 各 worker 的 best_score，按 worker 下标序。长度分两种情形：
    /// - 并行路径（threads > 1 且 `run_local_search`）：长度 == threads；
    /// - 委托路径（threads == 1 或 `!run_local_search`，见 [`solve_parallel_detailed`]
    ///   第 2 条语义）：恒为 `vec![best.best_score]`（长度 1、winner == 0——该路径
    ///   与 [`crate::solver::solve`] 完全等价，无并行增益，报告仅包装单次求解）。
    pub worker_scores: Vec<HardSoftScore>,
}

/// 多起点并行求解，返回聚合最优（= `solve_parallel_detailed(...)?.best`）。
///
/// 语义与错误处理见 [`solve_parallel_detailed`]。
pub fn solve_parallel(
    model: Arc<ProblemModel>,
    cfg: &SolverConfig,
    threads: usize,
) -> Result<SearchResult, PlannerError> {
    solve_parallel_detailed(model, cfg, threads).map(|report| report.best)
}

/// [`spawn_and_join`] 的 worker 体类型：`(worker 下标, 共享初始解引用, 该 worker 配置)`。
type WorkerBody<'a> =
    &'a (dyn Fn(usize, &Solution, &LocalSearchConfig) -> Result<SearchResult, PlannerError> + Sync);

/// 生产 worker 体：克隆共享初始解并以其专属配置跑局部搜索。worker 下标参数生产
/// 路径不使用——供 [`spawn_and_join`] 的注入式测试区分 worker。
fn run_local_search_worker(
    _worker: usize,
    initial: &Solution,
    cfg: &LocalSearchConfig,
) -> Result<SearchResult, PlannerError> {
    local_search(initial.clone(), cfg)
}

/// 按句柄序（= spawn 序 = worker 下标序）join 全部线程：
/// - `Ok(Ok(r))`：`results` 为 `Some` 时按下标序收集（spawn 失败清理路径传 `None` 丢弃）；
/// - `Ok(Err(pe))`：记录**首个**（即最小下标）worker Err，作为返回值透传；
/// - `Err(payload)`：暂存**首个** panic payload 到 `panic_slot`（仍 join 完其余线程）。
fn join_handles(
    handles: Vec<std::thread::ScopedJoinHandle<'_, Result<SearchResult, PlannerError>>>,
    mut results: Option<&mut Vec<SearchResult>>,
    panic_slot: &mut Option<Box<dyn Any + Send>>,
) -> Option<PlannerError> {
    let mut first_err: Option<PlannerError> = None;
    for handle in handles {
        match handle.join() {
            Ok(Ok(result)) => {
                if let Some(rs) = results.as_mut() {
                    rs.push(result);
                }
            }
            Ok(Err(pe)) => {
                // 同 cfg 下各 worker 的校验结果恒同；取最小下标者透传
                if first_err.is_none() {
                    first_err = Some(pe);
                }
            }
            Err(payload) => {
                if panic_slot.is_none() {
                    *panic_slot = Some(payload);
                }
            }
        }
    }
    first_err
}

/// [`solve_parallel_detailed`] 的派生 + join 内核（抽出为独立函数以提供注入式测试点：
/// worker 体与模拟 spawn 失败；生产经 [`run_local_search_worker`] 与 `fail_spawn_at =
/// None` 进入，两条注入路径生产二进制均不可达）。
///
/// 语义（normative）：
/// 1. 依次为 worker i ∈ [0, worker_cfgs.len()) 派线程执行
///    `body(i, initial, &worker_cfgs[i])`（initial 以只读引用借入、cfg 按所有权移入，
///    每 worker 在自己线程内取用）；
/// 2. 某次 `spawn_scoped` io 失败（资源不足等）：停止继续 spawn，**先按下标序 join
///    已派生的全部线程**——其中的 panic payload 照常暂存（见第 4 条，「worker panic
///    不吞」不因 spawn 失败破例），正常/Err 结果丢弃——然后返回
///    `Err(PlannerError::Io)`（该双故障角落下 worker 的 Err 被 Io 覆盖，panic 除外）；
/// 3. 全部派生成功：按下标序 join；有 panic → 返回占位 Err（不可观察，见第 4 条）；
///    否则首个（最小下标）worker Err 透传；再否则全部结果按下标序返回；
/// 4. 任一 worker panic：在 scope **外**（全部线程已 join）[`std::panic::resume_unwind`]
///    最小下标 payload——优先于任何 Err 返回值，不吞不锁死。
fn spawn_and_join(
    initial: &Solution,
    worker_cfgs: Vec<LocalSearchConfig>,
    body: WorkerBody<'_>,
    fail_spawn_at: Option<usize>,
) -> Result<Vec<SearchResult>, PlannerError> {
    // worker panic payload 的暂存格：仅父线程（scope 闭包运行于父线程）读写，
    // 非 worker 间共享状态；scope 返回后若已置位则先 resume_unwind 再检查结果。
    let mut panic_payload: Option<Box<dyn Any + Send>> = None;

    let joined: Result<Vec<SearchResult>, PlannerError> = std::thread::scope(|scope| {
        // initial 以引用借入 'scope（move 闭包拷贝的是 &Solution 本身）；
        // cfg_i 按所有权移入。
        let mut handles = Vec::with_capacity(worker_cfgs.len());
        for (i, cfg_i) in worker_cfgs.into_iter().enumerate() {
            // fail_spawn_at：注入式测试专用——Some(i) 时第 i 次派生模拟 io 失败
            // （不真正建线程），走与真实失败完全相同的清理/返回路径；生产恒 None。
            let spawned = if fail_spawn_at == Some(i) {
                Err(std::io::Error::other("injected spawn failure (test-only)"))
            } else {
                std::thread::Builder::new().spawn_scoped(scope, move || body(i, initial, &cfg_i))
            };
            match spawned {
                Ok(handle) => handles.push(handle),
                Err(io) => {
                    // 资源不足等 io 失败：映射 Io 并停止继续 spawn。返回前先按下标序
                    // join 已派生线程——无泄漏无死锁，且其中的 panic payload 照常暂存
                    // （scope 外统一 resume_unwind，不吞）；结果丢弃（本路径以 Io 为
                    // 返回值）。
                    let _ = join_handles(handles, None, &mut panic_payload);
                    return Err(PlannerError::Io(io));
                }
            }
        }

        // 按下标序 join（handles 即 spawn 序 = worker 下标序）：
        // 首个 panic / 首个 Err 即最小下标；仍 join 完其余线程后才返回。
        let mut results = Vec::with_capacity(handles.len());
        let first_err = join_handles(handles, Some(&mut results), &mut panic_payload);
        if panic_payload.is_some() {
            // 占位 Err：永不逃逸——panic_payload 已置位，scope 返回后先 resume_unwind
            // 再检查 joined（见下），本值不会被调用方观察到。
            return Err(PlannerError::Invalid(
                "worker thread panicked (payload stashed for resume_unwind)".to_string(),
            ));
        }
        if let Some(pe) = first_err {
            return Err(pe);
        }
        Ok(results)
    });

    // worker panic：在 scope 外（全部线程已 join）重抛最小下标 payload，不吞不锁死。
    if let Some(payload) = panic_payload {
        std::panic::resume_unwind(payload);
    }
    joined
}

/// 多起点并行求解（含各 worker 明细）。
///
/// 语义（normative）：
/// 1. `threads == 0 || threads > `[`MAX_THREADS`] → `Err(PlannerError::Invalid(
///    format!("threads must be within [1, {MAX_THREADS}] (got {threads})")))`，不 panic；
/// 2. `threads == 1 || !cfg.run_local_search` → 与 [`crate::solver::solve`] 完全等价
///    （直接委托）；本函数包装 `ParallelReport { best, winner: 0,
///    worker_scores: vec![best.best_score] }`（构造确定性，并行无增益）；
/// 3. 父线程先行 `model.validate()?` 与 `construct(&model, cfg.construction)?`
///    （与 [`crate::solver::solve`] 同序），Err 在派生任何线程之前返回；
/// 4. worker i ∈ [0, threads)：`cfg_i = cfg.local_search.clone()` 仅替换 seed ——
///    i == 0 保留 `cfg.local_search.seed` 原值；i ≥ 1 用
///    [`rng::derive_seed(cfg.local_search.seed, i as u64)`](crate::rng::derive_seed)；
///    各自调用 `local_search(initial.clone(), &cfg_i)`（initial 为父线程构造解）；
/// 5. 聚合：任一 worker `Err(PlannerError)` → 返回最小下标的 Err（同 cfg 下校验与
///    seed 无关 ⇒ 各 worker Err 恒同）；否则 winner = argmin_i(best_score_i)，
///    严格小于才替换 ⇒ 并列取最小下标，与线程完成顺序无关；
/// 6. 派生线程 io 失败（如资源不足）：返回 `Err(PlannerError::Io)`；返回前先 join
///    已派生线程——其中的 worker panic 仍优先重抛（优先级与清理细节见
///    [`spawn_and_join`] 第 2/4 条）。
pub fn solve_parallel_detailed(
    model: Arc<ProblemModel>,
    cfg: &SolverConfig,
    threads: usize,
) -> Result<ParallelReport, PlannerError> {
    // ---- threads 合法域（先于一切；不 panic）----
    if threads == 0 || threads > MAX_THREADS {
        return Err(PlannerError::Invalid(format!(
            "threads must be within [1, {MAX_THREADS}] (got {threads})"
        )));
    }

    // ---- 无并行增益路径：直接委托 solve（threads=1 行为与单线程版逐位一致）----
    if threads == 1 || !cfg.run_local_search {
        let best = solve(model, cfg)?;
        let best_score = best.best_score;
        return Ok(ParallelReport { best, winner: 0, worker_scores: vec![best_score] });
    }

    // ---- 与 solver::solve 同序：validate → construct（Err 在派生线程前返回）----
    model.validate()?;
    let initial = construct(&model, cfg.construction)?;

    // ---- worker 配置：克隆并仅替换 seed（worker 0 保留主 seed 原值）----
    let worker_cfgs: Vec<LocalSearchConfig> = (0..threads)
        .map(|i| {
            let mut cfg_i = cfg.local_search.clone();
            if i > 0 {
                cfg_i.seed = rng::derive_seed(cfg.local_search.seed, i as u64);
            }
            cfg_i
        })
        .collect();

    // ---- 派生 + join + panic 语义内聚于 spawn_and_join（并作为注入式测试点）：
    // 此处仅消费其按下标序返回的各 worker 结果。----
    let results = spawn_and_join(&initial, worker_cfgs, &run_local_search_worker, None)?;

    // ---- 聚合：min(best_score)，严格小于才替换 ⇒ 并列取最小下标 ----
    debug_assert_eq!(results.len(), threads, "every worker must report a result");
    let worker_scores: Vec<HardSoftScore> = results.iter().map(|r| r.best_score).collect();
    let mut winner = 0usize;
    for i in 1..results.len() {
        if results[i].best_score < results[winner].best_score {
            winner = i;
        }
    }
    let best = results
        .into_iter()
        .nth(winner)
        .expect("winner index is within results");
    Ok(ParallelReport { best, winner, worker_scores })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::construct::ConstructionHeuristic;
    use crate::localsearch::{Metaheuristic, Termination};
    use crate::model::{PlanningEntity, Resource};
    use crate::score::ScoreCalculator;

    /// 与 localsearch.rs 测试同构的玩具模型（m=4、2 资源、6 实体、冲突链 + 三角）。
    fn toy_model() -> Arc<ProblemModel> {
        let cap_sentinel = i64::MAX / 2;
        let model = ProblemModel {
            name: "par-toy".to_string(),
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

    /// SA 标量 i128 宽度模型（与 localsearch.rs 同构）：hard = 2^32、W = 2^40+1。
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

    /// 标准测试配置：tabu / Steps / FF。
    fn tabu_cfg(steps: usize) -> SolverConfig {
        SolverConfig {
            local_search: LocalSearchConfig {
                termination: Termination::Steps(steps),
                ..LocalSearchConfig::default()
            },
            ..SolverConfig::default()
        }
    }

    /// threads 越界：0 与 1025 → Err(Invalid)（不 panic）；1024 为合法值——经
    /// run_local_search=false 的委托路径只验证边界、不实际 spawn 1024 线程。
    #[test]
    fn threads_out_of_bounds_returns_invalid() {
        let model = toy_model();
        let cfg = tabu_cfg(10);
        for bad in [0usize, MAX_THREADS + 1] {
            match solve_parallel_detailed(model.clone(), &cfg, bad) {
                Err(PlannerError::Invalid(msg)) => {
                    assert_eq!(msg, format!("threads must be within [1, {MAX_THREADS}] (got {bad})"));
                }
                other => panic!("threads={bad}: expected Invalid, got {other:?}"),
            }
            match solve_parallel(model.clone(), &cfg, bad) {
                Err(PlannerError::Invalid(_)) => {}
                other => panic!("solve_parallel threads={bad}: expected Invalid, got {other:?}"),
            }
        }
        // 1024 合法（仅验证边界；run_local_search=false → 委托路径，零 spawn）
        let cfg_no_ls = SolverConfig { run_local_search: false, ..tabu_cfg(10) };
        let report =
            solve_parallel_detailed(model, &cfg_no_ls, MAX_THREADS).expect("1024 is in bounds");
        assert_eq!(report.winner, 0);
        assert_eq!(report.worker_scores.len(), 1);
    }

    /// worker-0 规则：detailed(…, 4).worker_scores[0] 与 solve 的 best_score 逐位相等
    /// （worker 0 = 主 seed + 同一构造起点 ⇒ 搜索完全相同）。
    #[test]
    fn worker_zero_score_equals_solve() {
        let model = toy_model();
        let cfg = tabu_cfg(1500);
        let serial = solve(model.clone(), &cfg).expect("solve");
        let report = solve_parallel_detailed(model, &cfg, 4).expect("parallel");
        assert_eq!(report.worker_scores.len(), 4);
        assert_eq!(
            report.worker_scores[0],
            serial.best_score,
            "worker 0 must run with the main seed verbatim"
        );
        // 胜者一致性：best 即胜出 worker 的原样结果
        assert_eq!(report.best.best_score, report.worker_scores[report.winner]);
        // Steps 终止下聚合不劣于单线程（结构性质）
        assert!(report.best.best_score <= serial.best_score);
    }

    /// 聚合 == worker_scores 的 min 且并列取最小下标：Steps(0) 下全体 worker 同分
    /// （都返回构造解），winner 必为 0。
    #[test]
    fn aggregate_is_min_with_min_index_on_ties() {
        let model = toy_model();
        let cfg = tabu_cfg(0);
        let report =
            solve_parallel_detailed(model.clone(), &cfg, 5).expect("steps(0) parallel");
        assert_eq!(report.worker_scores.len(), 5);
        for s in &report.worker_scores {
            assert_eq!(*s, report.worker_scores[0], "Steps(0): all workers tie");
        }
        assert_eq!(report.winner, 0, "ties must resolve to the minimal index");
        assert_eq!(report.best.steps, 0);
        let expect = ScoreCalculator::new(model.clone()).total_score(
            &construct(&model, ConstructionHeuristic::FirstFit).expect("construct"),
        );
        assert_eq!(report.best.best_score, expect);
        assert_eq!(report.best.best_score, *report.worker_scores.iter().min().unwrap());
    }

    /// 固定 seed + 固定 threads + Steps 终止 ⇒ 两次运行逐位一致
    /// （best.assignment / best_score / steps / winner / worker_scores）。
    #[test]
    fn parallel_is_bit_reproducible() {
        let model = toy_model();
        let cfg = tabu_cfg(1500);
        let a = solve_parallel_detailed(model.clone(), &cfg, 4).expect("run a");
        let b = solve_parallel_detailed(model, &cfg, 4).expect("run b");
        assert_eq!(a.best.best.assignment(), b.best.best.assignment());
        assert_eq!(a.best.best_score, b.best.best_score);
        assert_eq!(a.best.steps, b.best.steps);
        assert_eq!(a.winner, b.winner);
        assert_eq!(a.worker_scores, b.worker_scores);
    }

    /// 5 轮 threads=4 重复：无 panic 且逐位一致。
    #[test]
    fn five_rounds_threads4_are_repeatable() {
        let model = toy_model();
        let cfg = tabu_cfg(800);
        let first = solve_parallel_detailed(model.clone(), &cfg, 4).expect("round 1");
        for round in 2..=5 {
            let r = solve_parallel_detailed(model.clone(), &cfg, 4)
                .unwrap_or_else(|e| panic!("round {round} failed: {e}"));
            assert_eq!(r.best.best.assignment(), first.best.best.assignment(), "round {round}");
            assert_eq!(r.best.best_score, first.best.best_score, "round {round}");
            assert_eq!(r.best.steps, first.best.steps, "round {round}");
            assert_eq!(r.worker_scores, first.worker_scores, "round {round}");
        }
    }

    /// threads=1 与 run_local_search=false 两条路径均与 solve 完全等价。
    #[test]
    fn delegation_paths_equal_solve() {
        let model = toy_model();
        // threads = 1
        let cfg = tabu_cfg(1200);
        let serial = solve(model.clone(), &cfg).expect("solve");
        let one = solve_parallel(model.clone(), &cfg, 1).expect("threads=1");
        assert_eq!(one.best.assignment(), serial.best.assignment());
        assert_eq!(one.best_score, serial.best_score);
        assert_eq!(one.initial_score, serial.initial_score);
        assert_eq!(one.steps, serial.steps);
        // run_local_search = false（threads > 1 也走委托）
        let cfg_no_ls = SolverConfig { run_local_search: false, ..tabu_cfg(1200) };
        let serial = solve(model.clone(), &cfg_no_ls).expect("solve no-ls");
        let report =
            solve_parallel_detailed(model, &cfg_no_ls, 4).expect("parallel no-ls");
        assert_eq!(report.best.best.assignment(), serial.best.assignment());
        assert_eq!(report.best.best_score, serial.best_score);
        assert_eq!(report.best.steps, serial.steps);
        assert_eq!(report.winner, 0);
        assert_eq!(report.worker_scores, vec![serial.best_score]);
    }

    /// 三种元启发式小步数冒烟（threads=3）：无 panic、结果自洽
    /// （完整、域内、分数与全量重算相等、不劣于初始）。
    #[test]
    fn three_metaheuristics_smoke() {
        let model = toy_model();
        let metas = [
            Metaheuristic::TabuSearch { tabu_size: 8 },
            Metaheuristic::SimulatedAnnealing { start_temp: 1000.0, end_temp: 1.0 },
            Metaheuristic::LateAcceptance { history_len: 30 },
        ];
        for meta in &metas {
            let cfg = SolverConfig {
                local_search: LocalSearchConfig {
                    metaheuristic: meta.clone(),
                    termination: Termination::Steps(300),
                    ..LocalSearchConfig::default()
                },
                ..SolverConfig::default()
            };
            let report =
                solve_parallel_detailed(model.clone(), &cfg, 3)
                    .unwrap_or_else(|e| panic!("{meta:?}: parallel failed: {e}"));
            assert!(report.best.best.is_complete(), "{meta:?}: best incomplete");
            assert!(report.best.best_score <= report.best.initial_score, "{meta:?}: worse than initial");
            assert_eq!(report.best.steps, 300, "{meta:?}: Steps termination");
            for e in 0..model.entity_count() {
                assert!(
                    model.in_domain(e, report.best.best.value(e)),
                    "{meta:?}: value out of domain"
                );
            }
            let calc = ScoreCalculator::new(model.clone());
            assert_eq!(
                report.best.best_score,
                calc.total_score(&report.best.best),
                "{meta:?}: best_score != full recompute"
            );
        }
    }

    /// SA 在 worker 线程内同样必须走 i128 标量（复刻 localsearch 的
    /// sa_survives_scalar_magnitudes_beyond_i64 场景；误用 i64 乘法 debug 下直接 panic）。
    #[test]
    fn sa_parallel_survives_scalar_magnitudes_beyond_i64() {
        let model = wide_scalar_model();
        let cfg = SolverConfig {
            local_search: LocalSearchConfig {
                metaheuristic: Metaheuristic::SimulatedAnnealing {
                    start_temp: 1000.0,
                    end_temp: 1.0,
                },
                termination: Termination::Steps(500),
                ..LocalSearchConfig::default()
            },
            ..SolverConfig::default()
        };
        let report = solve_parallel_detailed(model.clone(), &cfg, 2).expect("sa parallel");
        // 多起点入口经 FF 构造起点（[1,1]，可行、soft 2^39 —— 与直接以 [0,0] 起跑的
        // localsearch 既有测试不同）。SA 每步预演 Change→值0 的移动时
        // s = 2^31·(2^40+1) ≈ 2^71 > i64::MAX：误用 i64 乘法会在 worker 线程内
        // overflow panic（debug），经 resume_unwind 令本测试失败。
        let initial =
            construct(&model, ConstructionHeuristic::FirstFit).expect("construct wide");
        let expect_initial = ScoreCalculator::new(model.clone()).total_score(&initial);
        assert_eq!(expect_initial, HardSoftScore::new(0, 2 * (1 << 38)), "FF start");
        assert_eq!(report.best.initial_score, expect_initial);
        assert!(report.best.best_score <= report.best.initial_score);
        let calc = ScoreCalculator::new(model);
        assert_eq!(report.best.best_score, calc.total_score(&report.best.best));
    }

    /// 非法局部搜索配置：worker 的 Err 以最小下标透传（与 solve 的 Err 一致）。
    #[test]
    fn invalid_worker_config_propagates_error() {
        let model = toy_model();
        let bad = SolverConfig {
            local_search: LocalSearchConfig {
                metaheuristic: Metaheuristic::TabuSearch { tabu_size: 0 },
                termination: Termination::Steps(100),
                ..LocalSearchConfig::default()
            },
            ..SolverConfig::default()
        };
        let par = solve_parallel_detailed(model.clone(), &bad, 3);
        let ser = solve(model, &bad);
        match (par, ser) {
            (Err(PlannerError::Invalid(_)), Err(PlannerError::Invalid(_))) => {}
            (p, s) => panic!("expected Invalid from both, got parallel={p:?} serial={s:?}"),
        }
    }

    /// 非法模型：validate 的 Err 在派生任何线程之前返回（Err(Model)，不 panic）。
    #[test]
    fn invalid_model_errors_before_spawning() {
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
        match solve_parallel_detailed(model, &tabu_cfg(10), 4) {
            Err(PlannerError::Model(_)) => {}
            other => panic!("expected Model error, got {other:?}"),
        }
    }

    // ---- P4 回归：spawn 失败 / worker panic 的单故障与双故障角落 ----
    // （经 spawn_and_join 的注入式测试点：fail_spawn_at 模拟派生 io 失败、
    //   自定义 body 模拟 worker panic / Err；生产体 run_local_search_worker 与
    //   fail_spawn_at=None 的胶水由上方全部 solve_parallel_detailed 测试覆盖。）

    /// 测试用 SearchResult：best = 共享初始解克隆、steps 编码 worker 下标
    /// （仅用于断言结果按下标序）。
    fn canned(initial: &Solution, worker: usize) -> SearchResult {
        SearchResult {
            best: initial.clone(),
            best_score: HardSoftScore::new(0, worker as i64),
            initial_score: HardSoftScore::new(0, 0),
            steps: worker,
        }
    }

    /// 三个默认 worker 配置（种子不参与注入式 body，内容任意合法即可）。
    fn three_cfgs() -> Vec<LocalSearchConfig> {
        (0..3).map(|_| LocalSearchConfig::default()).collect()
    }

    /// 双故障角落（P4 核心）：worker 0 panic 且第 2 次派生注入 io 失败 ⇒
    /// 已派生线程仍被按下标序 join、panic payload **不被吞**——spawn_and_join 在
    /// scope 外 resume_unwind（catch_unwind 捕获，payload 内容保持最小下标者）。
    #[test]
    fn spawn_failure_still_unwinds_panicking_worker() {
        let model = toy_model();
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let body = |worker: usize, init: &Solution, _cfg: &LocalSearchConfig| {
            if worker == 0 {
                panic!("injected worker panic (worker 0)");
            }
            Ok(canned(init, worker))
        };
        // 静音 worker 线程 panic 的默认打印（测试输出整洁；不影响 payload 传递）
        let prev_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            // 派生 worker 0、1 成功（0 会 panic），第 2 次派生注入失败
            spawn_and_join(&initial, three_cfgs(), &body, Some(2))
        }));
        std::panic::set_hook(prev_hook);
        let payload = match outcome {
            Err(payload) => payload,
            Ok(_) => panic!("double fault must unwind: worker panic was swallowed"),
        };
        let msg = payload
            .downcast_ref::<&str>()
            .map(|s| s.to_string())
            .or_else(|| payload.downcast_ref::<String>().cloned())
            .expect("panic payload must be &str or String");
        assert_eq!(msg, "injected worker panic (worker 0)");
    }

    /// 单故障（仅 worker panic，无 spawn 失败）：最小下标 payload 于 scope 外重抛
    /// （既有路径的回归钉死，确保修复未改变该优先级）。
    #[test]
    fn worker_panic_without_spawn_failure_unwinds() {
        let model = toy_model();
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let body = |worker: usize, init: &Solution, _cfg: &LocalSearchConfig| {
            if worker == 1 {
                panic!("injected worker panic (worker 1)");
            }
            Ok(canned(init, worker))
        };
        let prev_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            spawn_and_join(&initial, three_cfgs(), &body, None)
        }));
        std::panic::set_hook(prev_hook);
        let payload = match outcome {
            Err(payload) => payload,
            Ok(_) => panic!("worker panic without spawn failure must unwind"),
        };
        let msg = payload
            .downcast_ref::<&str>()
            .map(|s| s.to_string())
            .or_else(|| payload.downcast_ref::<String>().cloned())
            .expect("panic payload must be &str or String");
        assert_eq!(msg, "injected worker panic (worker 1)");
    }

    /// 单故障（仅 spawn 失败，无 panic）：注入第 1 次派生失败 ⇒ 已派生的 worker 0
    /// 被按下标序 join（结果丢弃、不泄漏不挂起），返回 Err(Io)、不 panic。
    #[test]
    fn injected_spawn_failure_returns_io() {
        let model = toy_model();
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let body =
            |worker: usize, init: &Solution, _cfg: &LocalSearchConfig| Ok(canned(init, worker));
        match spawn_and_join(&initial, three_cfgs(), &body, Some(1)) {
            Err(PlannerError::Io(_)) => {}
            Err(other) => panic!("expected Io, got {other:?}"),
            Ok(_) => panic!("injected spawn failure must return Err(Io)"),
        }
    }

    /// 正常路径：全部派生成功时结果按 worker 下标序返回（steps 编码下标）。
    #[test]
    fn spawn_and_join_returns_results_in_index_order() {
        let model = toy_model();
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let body =
            |worker: usize, init: &Solution, _cfg: &LocalSearchConfig| Ok(canned(init, worker));
        let results =
            spawn_and_join(&initial, three_cfgs(), &body, None).expect("all workers succeed");
        assert_eq!(
            results.iter().map(|r| r.steps).collect::<Vec<_>>(),
            vec![0, 1, 2],
            "results must be in worker index order"
        );
    }

    /// worker Err 以最小下标透传（经重构后的 join_handles 收集路径）。
    #[test]
    fn worker_err_propagates_min_index() {
        let model = toy_model();
        let initial = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let body = |worker: usize, init: &Solution, _cfg: &LocalSearchConfig| {
            if worker >= 1 {
                Err(PlannerError::Invalid(format!("injected worker {worker} error")))
            } else {
                Ok(canned(init, worker))
            }
        };
        match spawn_and_join(&initial, three_cfgs(), &body, None) {
            Err(PlannerError::Invalid(msg)) => {
                assert_eq!(msg, "injected worker 1 error", "min worker index wins")
            }
            other => panic!("expected Invalid, got {other:?}"),
        }
    }
}
