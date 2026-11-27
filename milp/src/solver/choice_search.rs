//! 大量班次共享互斥行时，试探候选搜索和全模型证明，未达标时交回完整树搜索。
use super::{paired_primal, restricted_primal, roster_proof};
use crate::branch_bound::{BbOptions, BbStats, SearchControl, WorkIncumbentObserver};
use crate::heuristics::Incumbent;
use crate::model::{Problem, Sense, Tolerances, VarType};
use crate::search_ctrl::gap_percent;
use crate::simplex::LpSolution;
use crate::simplex_sparse::{effective_engine, LpEngine, LpEngineKind};
use std::cell::{Cell, RefCell};
use std::time::{Duration, Instant};

pub(crate) struct ChoiceSearch {
    pub(crate) ordinary_deadline: Instant,
    solve_deadline: Option<Instant>,
    target: f64,
    costs: Vec<f64>,
    pub(crate) paired_rows: bool,
}

impl ChoiceSearch {
    pub(crate) fn prepare(
        p: &Problem,
        opts: &BbOptions,
        ctrl: &SearchControl,
        start: Instant,
        solve_deadline: Option<Instant>,
        lp: &dyn LpEngine,
    ) -> Option<Self> {
        let target = ctrl.mipgap_percent?;
        // 调用者已检查单线程、正缺口、墙钟模式及大组互斥列。
        // 候选搜索目前只支持有非负整数余量的稀疏最小化问题。
        let now = Instant::now();
        if p.sense != Sense::Minimize
            || effective_engine(p, opts.lp_engine) != LpEngineKind::Revised
            || !p
                .vars
                .iter()
                .any(|v| v.var_type == VarType::Integer && v.lb == 0. && v.ub == f64::INFINITY)
            || solve_deadline.is_some_and(|end| now >= end)
        {
            return None;
        }
        let ordinary_budget = Duration::from_secs(8);
        let ordinary_deadline = solve_deadline.map_or(now + ordinary_budget, |end| {
            start + (end.duration_since(start) / 6).min(ordinary_budget)
        });
        Some(Self {
            ordinary_deadline,
            solve_deadline,
            target,
            costs: lp.structural_reduced_costs(),
            paired_rows: false,
        })
    }

    /// 互斥潜水可能用完整个普通预算。先试共享余量的候选搜索，有解即可直接做证明。
    pub(crate) fn start_paired(
        &mut self,
        p: &Problem,
        root: &LpSolution,
        opts: &BbOptions,
        tol: Tolerances,
        upper: Option<f64>,
        used_nodes: u64,
        accept: &dyn Fn(Incumbent),
    ) -> BbStats {
        if !opts.heuristics
            || upper
                .and_then(|u| gap_percent(p.sense, u, root.objective))
                .is_some_and(|g| g <= 2. * self.target)
        {
            return BbStats::default();
        }
        let nodes = if opts.node_limit == 0 {
            u64::MAX
        } else {
            opts.node_limit.saturating_sub(used_nodes)
        };
        let end = self.candidate_deadline(2, Duration::from_secs(24));
        let found = Cell::new(false);
        let stats = paired_primal::search(
            p, root, &self.costs, opts, tol, end, nodes, self.target, &|inc| {
                found.set(true);
                accept(inc);
            },
        );
        self.paired_rows = found.get();
        stats
    }

    fn candidate_deadline(&self, fraction: u32, cap: Duration) -> Instant {
        let now = Instant::now();
        let budget = self.solve_deadline.map_or(cap, |end| {
            (end.saturating_duration_since(now) / fraction).min(cap)
        });
        now + budget
    }

    /// 缺口至多为目标的两倍时，先试根与一对完整子域；任一侧中断仍沿用父界。
    pub(crate) fn probe(
        &self,
        p: &Problem,
        root: &LpSolution,
        upper: f64,
        tol: &Tolerances,
        nodes: u64,
    ) -> Option<(Option<Incumbent>, f64, BbStats)> {
        let gap = gap_percent(p.sense, upper, root.objective)?;
        if gap <= self.target || gap > 2. * self.target {
            return None;
        }
        self.certify(
            p,
            root,
            upper,
            tol,
            self.ordinary_deadline,
            nodes.min(3),
            false,
        )
    }

    fn certify(
        &self,
        p: &Problem,
        root: &LpSolution,
        upper: f64,
        tol: &Tolerances,
        deadline: Instant,
        nodes: u64,
        paired_rows: bool,
    ) -> Option<(Option<Incumbent>, f64, BbStats)> {
        if nodes == 0 || Instant::now() >= deadline {
            return None;
        }
        let mut cuts = roster_proof::separate(p, root, &self.costs, tol, deadline);
        if paired_rows {
            // 初始根已经满足的组覆盖割先不入表。任何子集仍给出完整整数域的有效松弛。
            cuts.extend(
                roster_proof::group_cover::separate(p, tol, deadline)
                    .into_iter()
                    .filter(|cut| {
                        let activity: f64 = cut.terms.iter().map(|&(j, a)| a * root.x[j]).sum();
                        activity > cut.rhs + tol.feasibility
                    }),
            );
        }
        if cuts.is_empty() {
            return None;
        }
        let proof = roster_proof::reduction::prove_compact(
            p,
            &cuts,
            root,
            upper,
            self.target,
            tol,
            deadline,
            nodes,
            paired_rows,
        );
        let stats = BbStats {
            nodes: proof.nodes,
            lp_iterations: proof.lp_iterations,
            cuts: if proof.nodes > 0 {
                cuts.len() as u64
            } else {
                0
            },
            ..Default::default()
        };
        Some((proof.candidate, proof.bound, stats))
    }

    pub(crate) fn improve(
        self,
        p: &Problem,
        root: &LpSolution,
        opts: &BbOptions,
        tol: Tolerances,
        current: Option<f64>,
        mut bound: f64,
        used_nodes: u64,
        observer: WorkIncumbentObserver<'_>,
    ) -> (Option<Incumbent>, f64, BbStats) {
        let nodes = if opts.node_limit == 0 {
            u64::MAX
        } else {
            opts.node_limit.saturating_sub(used_nodes)
        };
        let best_value = Cell::new(current.unwrap_or(f64::INFINITY));
        let best = RefCell::new(None);
        let improvements = Cell::new(0u8);
        let accept = |candidate: Incumbent| {
            if candidate.objective < best_value.get() {
                if let Some(report) = observer {
                    report(&candidate.values);
                }
                improvements.set(improvements.get().saturating_add(1));
                best_value.set(candidate.objective);
                best.replace(Some(candidate));
            }
        };
        // 有总时限时，至少保留剩余预算的三分之二做证明。
        let primal_end = self.candidate_deadline(3, Duration::from_secs(16));
        let paired_rows = self.paired_rows;
        let mut stats = BbStats::default();
        // 近目标缺口已在 probe 试过。差距较大时，先用小额度的全域证明试探：
        // 有效割可能同时收紧界、找出新解，及早达到目标即可省去受限搜索。
        if !paired_rows {
            if let Some(upper) = current.filter(|&upper| {
                gap_percent(p.sense, upper, bound).is_some_and(|g| g > 2. * self.target)
            }) {
                let deadline = self.candidate_deadline(6, Duration::from_secs(2));
                if let Some((candidate, lower, spent)) =
                    self.certify(p, root, upper, &tol, deadline, nodes.min(32), false)
                {
                    if let Some(inc) = candidate {
                        accept(inc);
                    }
                    if lower <= best_value.get() {
                        bound = bound.max(lower);
                    }
                    stats = spent;
                }
            }
        }
        let meets_target = |bound| {
            gap_percent(p.sense, best_value.get(), bound).is_some_and(|g| g <= self.target)
        };
        if !paired_rows && !meets_target(bound) {
            let prior_improvements = improvements.get();
            let spent = restricted_primal::search(
                p,
                root,
                &self.costs,
                &BbOptions {
                    node_limit: nodes.saturating_sub(stats.nodes),
                    ..*opts
                },
                tol,
                primal_end,
                &accept,
                // 只用完整模型的界及已复核的原模型可行解判断是否达到目标。
                Some(&|_| meets_target(bound) || (current.is_some()
                    && improvements.get().saturating_sub(prior_improvements) >= 3)),
            );
            stats.nodes += spent.nodes;
            stats.cuts += spent.cuts;
            stats.lp_iterations += spent.lp_iterations;
        }
        let remaining_nodes = nodes.saturating_sub(stats.nodes);
        if best_value.get().is_finite() && !meets_target(bound) {
            // 无总时限时只限制这次证明试探，未证部分留给完整搜索。
            let proof_deadline = self.solve_deadline
                .unwrap_or_else(|| Instant::now() + Duration::from_secs(40));
            if let Some((candidate, lower, spent)) = self.certify(
                p,
                root,
                best_value.get(),
                &tol,
                proof_deadline,
                remaining_nodes,
                paired_rows,
            ) {
                if let Some(inc) = candidate {
                    accept(inc);
                }
                if lower <= best_value.get() {
                    bound = bound.max(lower);
                }
                stats.nodes += spent.nodes;
                stats.cuts += spent.cuts;
                stats.lp_iterations += spent.lp_iterations;
            }
        }
        (best.into_inner(), bound, stats)
    }
}

#[cfg(test)]
mod tests {
    use crate::branch_bound::{BbOptions, SearchControl};
    use crate::search_ctrl::Determinism;
    use crate::model::{Problem, RelOp, Sense, SolveStatus, VarType, Variable};
    use crate::simplex_sparse::LpEngineKind;
    use crate::solver::{solve_controlled, SolveOptions};

    #[test]
    fn full_model_clique_closes_gap_only_when_ordinary_search_needs_it() {
        let mut p = Problem::new("互斥三角形", Sense::Minimize);
        p.add_vars(
            (0..4096)
                .map(|i| Variable {
                    name: format!("x{i}"),
                    var_type: VarType::Binary,
                    lb: 0.,
                    ub: 1.,
                })
                .collect(),
        )
        .unwrap();
        p.add_var("slack", VarType::Integer).unwrap();
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            p.add_constraint(None, vec![(i, 1.), (j, 1.)], RelOp::Le, 1.)
                .unwrap();
        }
        p.add_constraint(None, (3..4096).map(|i| (i, 1.)).collect(), RelOp::Le, 1.)
            .unwrap();
        // 三个选项两两互斥，整数最优为 -1，连续松弛可取三个 0.5 得到 -1.5。
        p.set_objective(vec![(0, -1.), (1, -1.), (2, -1.)], 0.)
            .unwrap();
        let mut opts = SolveOptions {
            use_presolve: false,
            use_scaling: false,
            bb: BbOptions {
                lp_engine: LpEngineKind::Revised,
                ..Default::default()
            },
            ..Default::default()
        };
        for (time_limit_ms, determinism, target, expected_cuts) in [
            (10_000, Determinism::Wall, 1., 1),
            (0, Determinism::Wall, 1., 1),
            (0, Determinism::Wall, 60., 0),
            // 精确搜索和节点预算模式沿用常规割循环。
            (0, Determinism::Wall, 0., 4),
            (0, Determinism::Budget, 1., 4),
        ] {
            opts.bb.time_limit_ms = time_limit_ms;
            let ctrl = SearchControl {
                determinism,
                mipgap_percent: Some(target),
                ..Default::default()
            };
            let result = solve_controlled(&p, &opts, &ctrl, None).unwrap();
            assert_eq!(result.report.objective, Some(-1.));
            assert_eq!(result.report.status, SolveStatus::Optimal);
            assert_eq!(result.bound, Some(-1.));
            assert_eq!(result.report.stats.cuts, expected_cuts,
                "time_limit_ms={time_limit_ms}, determinism={determinism:?}, target={target}");
        }
    }
}
