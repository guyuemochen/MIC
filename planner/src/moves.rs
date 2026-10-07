//! 移动（Change/Swap）与随机邻域生成。随机数消耗顺序与跨类回退规则逐字钉死，
//! 保证同一二进制同种子产生相同移动序列。

use crate::model::{EntityId, HardSoftScore, Solution, ValueId, UNASSIGNED};
use crate::rng::Rng;
use crate::score::IncrementalScorer;

/// 单个分支的采样重试上限（normative：Change / Swap 分支各至多 8 次，
/// 首选分支与回退分支独立计数）。
const BRANCH_RETRIES: usize = 8;

/// 原子移动。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Move {
    /// 把单个实体的取值改为 `new_value`（须在域内且 ≠ 当前值）。
    Change { entity: EntityId, new_value: ValueId },
    /// 互换 `a`、`b` 的取值（双方已分配、值不同、val(b) ∈ dom(a) 且 val(a) ∈ dom(b)）。
    Swap { a: EntityId, b: EntityId },
}

/// 随机移动的种类配比。
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MoveKind {
    /// 只生成 Change。
    Change,
    /// 只生成 Swap。
    Swap,
    /// `change_ratio ∈ [0, 1]`：每次生成前先 `gen_ratio(change_ratio)`，
    /// true ⇒ Change 分支、false ⇒ Swap 分支。
    Mixed { change_ratio: f64 },
}

impl Default for MoveKind {
    /// 默认 `Mixed { change_ratio: 0.5 }`。
    fn default() -> Self {
        MoveKind::Mixed { change_ratio: 0.5 }
    }
}

/// 均匀随机生成一个合法非平凡移动；返回 `None` 仅当全域退化
/// （如所有实体共用单值且全已分配）。
///
/// 随机数消耗顺序与回退规则（逐字钉死，单一无条件规则，与首选分支来源无关）：
/// 1. 首选分支选择：`kind = Change` → Change 分支；`kind = Swap` → Swap 分支；
///    `Mixed { r }`：r ≤ 0 → Swap 分支，r ≥ 1 → Change 分支，否则先 `gen_ratio(r)`：
///    true → Change 分支、false → Swap 分支。
/// 2. Change 分支：至多 8 次重试，每次 `gen_range(0, n)` 得 e、
///    `gen_range(0, dom[e].len())` 得 v，若 `sol.value(e) ≠ v` 即返回 `Change { e, v }`。
/// 3. Swap 分支：至多 8 次重试，每次 `gen_range(0, n)`、`gen_range(0, n)` 得 (a, b)，
///    满足 a ≠ b 且双方已赋值、值不同、val(b) ∈ dom(a) 且 val(a) ∈ dom(b) 即返回
///    `Swap { a, b }`。
/// 4. 首选分支 8 次重试失败后，**无条件**尝试另一分支（按 2)/3) 规则，至多 8 次）——
///    无论首选分支来自 kind 直选还是 Mixed 的随机选择；另一分支也 8 败则返回 `None`。
///    （跨类回退是必须的：m=3 小算例上无回退实测约 0.3% 步生成失败。）
pub fn random_move(
    sol: &Solution,
    inc: &IncrementalScorer,
    rng: &mut Rng,
    kind: MoveKind,
) -> Option<Move> {
    // 契约签名保留 inc（生成判定仅依赖 sol；sol 与 inc 的同步由 Move::apply 结构性保证）
    let _ = inc;
    let n = sol.model().entity_count();
    if n == 0 {
        // 全域退化兜底（validate 层拦截 n=0 模型；此处短路避免触碰 gen_range 的
        // low < high 前置条件），不消耗任何随机数
        return None;
    }
    // 1) 首选分支选择：Mixed 两端点（r ≤ 0 / r ≥ 1）不耗随机数，中间恰耗一个 gen_ratio
    let prefer_change = match kind {
        MoveKind::Change => true,
        MoveKind::Swap => false,
        MoveKind::Mixed { change_ratio: r } => {
            if r <= 0.0 {
                false
            } else if r >= 1.0 {
                true
            } else {
                rng.gen_ratio(r)
            }
        }
    };
    // 2)/3) 首选分支至多 8 次；4) 8 败后**无条件**试另一分支（同样至多 8 次），
    //     与首选分支来自 kind 直选还是 Mixed 随机选择无关
    if prefer_change {
        try_change(sol, rng).or_else(|| try_swap(sol, rng))
    } else {
        try_swap(sol, rng).or_else(|| try_change(sol, rng))
    }
}

/// Change 分支：每次重试消耗 `gen_range(0, n)`（实体）+ `gen_range(0, dom.len())`（值）。
/// 空域模型属非法输入（parse/validate 层拦截），`gen_range(0, 0)` 的 assert 属逻辑错误路径。
fn try_change(sol: &Solution, rng: &mut Rng) -> Option<Move> {
    let n = sol.model().entity_count();
    for _ in 0..BRANCH_RETRIES {
        let e = rng.gen_range(0, n);
        let dom = sol.model().domain(e);
        let v = dom[rng.gen_range(0, dom.len())];
        if sol.value(e) != v {
            return Some(Move::Change { entity: e, new_value: v });
        }
    }
    None
}

/// Swap 分支：每次重试消耗 `gen_range(0, n)` 两次得 (a, b)；
/// 合法 ⇔ a ≠ b ∧ 双方已赋值 ∧ 值不同 ∧ val(b) ∈ dom(a) ∧ val(a) ∈ dom(b)。
fn try_swap(sol: &Solution, rng: &mut Rng) -> Option<Move> {
    let n = sol.model().entity_count();
    for _ in 0..BRANCH_RETRIES {
        let a = rng.gen_range(0, n);
        let b = rng.gen_range(0, n);
        if a == b {
            continue;
        }
        let (va, vb) = (sol.value(a), sol.value(b));
        if va == UNASSIGNED || vb == UNASSIGNED || va == vb {
            continue;
        }
        if sol.model().in_domain(a, vb) && sol.model().in_domain(b, va) {
            return Some(Move::Swap { a, b });
        }
    }
    None
}

impl Move {
    /// 预演：返回应用本移动后的总分（内部走 [`IncrementalScorer::preview_change`] /
    /// [`IncrementalScorer::preview_swap`]），不修改任何状态。
    pub fn preview(&self, inc: &IncrementalScorer) -> HardSoftScore {
        match *self {
            Move::Change { entity, new_value } => inc.preview_change(entity, new_value),
            Move::Swap { a, b } => inc.preview_swap(a, b),
        }
    }

    /// 提交：同一函数内先对涉及实体 `debug_assert!(inc.value(e) == sol.value(e))`，
    /// 再同步 `Solution::set` + `IncrementalScorer::apply_*`（crate 内唯一解变更路径）。
    ///
    /// Swap 前置条件：双方已分配；违例 assert + 立即 panic（与
    /// [`IncrementalScorer::apply_swap`] 同语义：逻辑错误路径，正常数据流不可达）。
    pub fn apply(&self, sol: &mut Solution, inc: &mut IncrementalScorer) {
        match *self {
            Move::Change { entity, new_value } => {
                debug_assert_eq!(
                    inc.value(entity),
                    sol.value(entity),
                    "shadow value out of sync with solution at entity {entity}"
                );
                sol.set(entity, new_value);
                inc.apply_change(entity, new_value);
            }
            Move::Swap { a, b } => {
                assert!(
                    sol.value(a) != UNASSIGNED && sol.value(b) != UNASSIGNED,
                    "Move::apply Swap on unassigned entity"
                );
                debug_assert_eq!(
                    inc.value(a),
                    sol.value(a),
                    "shadow value out of sync with solution at entity {a}"
                );
                debug_assert_eq!(
                    inc.value(b),
                    sol.value(b),
                    "shadow value out of sync with solution at entity {b}"
                );
                let (va, vb) = (sol.value(a), sol.value(b));
                sol.set(a, vb);
                sol.set(b, va);
                inc.apply_swap(a, b);
            }
        }
    }

    /// 涉及实体列表与个数（Change 1 个，Swap 2 个），供禁忌管理使用。
    /// Change 的第 2 槽位填充为同一实体（消费方只读前 `count` 个）。
    pub fn entities(&self) -> ([EntityId; 2], usize) {
        match *self {
            Move::Change { entity, .. } => ([entity, entity], 1),
            Move::Swap { a, b } => ([a, b], 2),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{PlanningEntity, ProblemModel, Resource};
    use crate::score::ScoreCalculator;
    use std::sync::Arc;

    /// 与 score.rs 测试同构的玩具模型（m=4、2 资源、6 实体、冲突链 + 三角、
    /// 非平凡代价与需求），保证 Swap/Change 都有充足合法候选。
    fn toy_model() -> Arc<ProblemModel> {
        let cap_sentinel = i64::MAX / 2;
        let model = ProblemModel {
            name: "moves-toy".to_string(),
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

    /// 跨类回退用例模型：n=2、域均 {0,1}、无约束、代价 0。
    fn all_same_model() -> Arc<ProblemModel> {
        let model = ProblemModel {
            name: "all-same".to_string(),
            value_count: 2,
            resources: vec![],
            entities: vec![
                PlanningEntity { domain: vec![0, 1], demands: vec![], conflicts: vec![] },
                PlanningEntity { domain: vec![0, 1], demands: vec![], conflicts: vec![] },
            ],
            cost: vec![vec![0, 0], vec![0, 0]],
        };
        model.validate().expect("all-same model must be valid");
        Arc::new(model)
    }

    /// 构造同步对 (sol, inc)：赋值向量 → 完整解 + 已 reset 的影子评分器。
    fn setup(model: &Arc<ProblemModel>, assignment: &[ValueId]) -> (Solution, IncrementalScorer) {
        let sol =
            Solution::from_assignment(model.clone(), assignment).expect("assignment in domain");
        let mut inc = IncrementalScorer::new(model.clone());
        inc.reset(&sol);
        (sol, inc)
    }

    /// 契约合法性逐条断言（含"非平凡"：Change 不得为 no-op、Swap 值必须不同）。
    fn assert_legal(sol: &Solution, mv: &Move) {
        match *mv {
            Move::Change { entity, new_value } => {
                assert!(
                    sol.model().in_domain(entity, new_value),
                    "Change value {new_value} must be in domain of entity {entity}"
                );
                assert_ne!(
                    sol.value(entity),
                    new_value,
                    "Change must be non-trivial (new value equals current)"
                );
            }
            Move::Swap { a, b } => {
                assert_ne!(a, b, "Swap entities must differ");
                let (va, vb) = (sol.value(a), sol.value(b));
                assert_ne!(va, UNASSIGNED, "Swap: entity {a} must be assigned");
                assert_ne!(vb, UNASSIGNED, "Swap: entity {b} must be assigned");
                assert_ne!(va, vb, "Swap must swap different values (no no-op)");
                assert!(
                    sol.model().in_domain(a, vb),
                    "Swap: val(b)={vb} must be in dom(a={a})"
                );
                assert!(
                    sol.model().in_domain(b, va),
                    "Swap: val(a)={va} must be in dom(b={b})"
                );
            }
        }
    }

    /// 随机游走 1000 步：每次生成的移动必须合法（Swap 域相容、无 no-op）；
    /// preview == apply 后影子分 == 全量重算；影子赋值与解逐位一致。
    #[test]
    fn random_walk_moves_legal_and_preview_matches_apply_1000() {
        let model = toy_model();
        let (mut sol, mut inc) = setup(&model, &[3, 0, 3, 2, 0, 2]);
        let calc = ScoreCalculator::new(model.clone());
        let mut rng = Rng::new(2026_1007);
        for step in 0..1000 {
            let Some(mv) = random_move(
                &sol,
                &inc,
                &mut rng,
                MoveKind::Mixed { change_ratio: 0.5 },
            ) else {
                // 两分支均 8 败：无移动，随机数已照常消耗
                continue;
            };
            assert_legal(&sol, &mv);
            let preview = mv.preview(&inc);
            mv.apply(&mut sol, &mut inc);
            assert_eq!(
                inc.score(),
                preview,
                "step {step}: score after apply must equal preview"
            );
            assert_eq!(
                inc.score(),
                calc.total_score(&sol),
                "step {step}: shadow score must equal full recompute"
            );
            for e in 0..sol.model().entity_count() {
                assert_eq!(inc.value(e), sol.value(e), "step {step}: shadow drift at {e}");
            }
        }
    }

    /// 跨类回退用例 (a)：kind=Swap、全部实体当前同值且各域 ≥ 2 值
    /// ⇒ Swap 分支 8 次全败（任一候选 va == vb）⇒ 无条件回退 Change 分支。
    /// 逐 seed 断言（每 seed 结果确定）。
    #[test]
    fn fallback_swap_to_change_when_all_values_equal() {
        let model = all_same_model();
        let (sol, inc) = setup(&model, &[0, 0]);
        for seed in [42u64, 7, 2026, 1234567, 987654321] {
            let mut rng = Rng::new(seed);
            let mv = random_move(&sol, &inc, &mut rng, MoveKind::Swap);
            assert!(
                matches!(mv, Some(Move::Change { .. })),
                "seed {seed}: expected fallback to Change, got {mv:?}"
            );
        }
    }

    /// 跨类回退用例 (b)：所有实体域单值且全已分配该值
    /// ⇒ Change 分支（v 恒等于当前值）与 Swap 分支（值恒等）皆 8 败 ⇒ None。
    #[test]
    fn degenerate_single_value_domains_return_none() {
        let model = Arc::new(ProblemModel {
            name: "single".to_string(),
            value_count: 1,
            resources: vec![],
            entities: vec![
                PlanningEntity { domain: vec![0], demands: vec![], conflicts: vec![] },
                PlanningEntity { domain: vec![0], demands: vec![], conflicts: vec![] },
            ],
            cost: vec![vec![0], vec![0]],
        });
        model.validate().expect("single-value model must be valid");
        let (sol, inc) = setup(&model, &[0, 0]);
        for kind in [
            MoveKind::Change,
            MoveKind::Swap,
            MoveKind::Mixed { change_ratio: 0.5 },
            MoveKind::default(),
        ] {
            let mut rng = Rng::new(42);
            assert_eq!(
                random_move(&sol, &inc, &mut rng, kind),
                None,
                "kind {kind:?}: both branches must fail on single-value degenerate world"
            );
        }
    }

    /// 钉死回退路径的随机数消耗顺序：用同种子 Rng 独立模拟契约序列
    /// （Swap 8 × (gen_range,gen_range) 全败 → Change 8 × (gen_range,gen_range)），
    /// 与 random_move 实际产出逐位一致。
    #[test]
    fn fallback_consumes_rng_in_contract_order() {
        let model = all_same_model();
        let (sol, inc) = setup(&model, &[0, 0]);
        let n = 2usize;
        let seed = 4242u64;

        // —— 契约独立模拟（不经过 random_move）——
        let mut sim = Rng::new(seed);
        let mut expected = None;
        for _ in 0..8 {
            let a = sim.gen_range(0, n);
            let b = sim.gen_range(0, n);
            let ok = a != b
                && sol.value(a) != UNASSIGNED
                && sol.value(b) != UNASSIGNED
                && sol.value(a) != sol.value(b)
                && sol.model().in_domain(a, sol.value(b))
                && sol.model().in_domain(b, sol.value(a));
            if ok {
                expected = Some(Move::Swap { a, b });
                break;
            }
        }
        if expected.is_none() {
            for _ in 0..8 {
                let e = sim.gen_range(0, n);
                let dom = sol.model().domain(e);
                let v = dom[sim.gen_range(0, dom.len())];
                if sol.value(e) != v {
                    expected = Some(Move::Change { entity: e, new_value: v });
                    break;
                }
            }
        }

        let mut rng = Rng::new(seed);
        assert_eq!(
            random_move(&sol, &inc, &mut rng, MoveKind::Swap),
            expected,
            "random_move must match the contract-order simulation bit for bit"
        );
        assert!(matches!(expected, Some(Move::Change { .. })), "Swap 必 8 败 ⇒ 回退 Change");
    }

    /// Mixed 端点规则：r ≤ 0 ⇒ Swap 分支、r ≥ 1 ⇒ Change 分支，且不消耗额外随机数
    /// （同种子下与直选 kind 的结果逐位一致）。
    #[test]
    fn mixed_endpoints_match_direct_kinds_without_extra_draws() {
        let model = all_same_model();
        let (sol, inc) = setup(&model, &[0, 1]); // 此时 Swap 也可行，两端点各自走首选分支
        let seed = 99u64;
        let mut r1 = Rng::new(seed);
        let a = random_move(&sol, &inc, &mut r1, MoveKind::Swap);
        let mut r2 = Rng::new(seed);
        let b = random_move(&sol, &inc, &mut r2, MoveKind::Mixed { change_ratio: 0.0 });
        assert_eq!(a, b, "Mixed ratio 0.0 must equal Swap kind stream-for-stream");

        let mut r3 = Rng::new(seed);
        let c = random_move(&sol, &inc, &mut r3, MoveKind::Change);
        let mut r4 = Rng::new(seed);
        let d = random_move(&sol, &inc, &mut r4, MoveKind::Mixed { change_ratio: 1.0 });
        assert_eq!(c, d, "Mixed ratio 1.0 must equal Change kind stream-for-stream");
    }

    /// 首选分支统计（确定种子）：r=0 绝大多数产出 Swap、r=1 绝大多数产出 Change；
    /// 偶发回退样本仍必须合法。
    #[test]
    fn ratio_extremes_prefer_their_own_branch() {
        let model = toy_model();
        let (sol, inc) = setup(&model, &[3, 0, 3, 2, 0, 2]);
        let mut rng = Rng::new(7);
        let mut swaps = 0usize;
        let mut total = 0usize;
        for _ in 0..500 {
            if let Some(mv) =
                random_move(&sol, &inc, &mut rng, MoveKind::Mixed { change_ratio: 0.0 })
            {
                assert_legal(&sol, &mv);
                total += 1;
                if matches!(mv, Move::Swap { .. }) {
                    swaps += 1;
                }
            }
        }
        assert!(
            swaps * 100 > total * 90,
            "ratio 0 should overwhelmingly produce Swaps (got {swaps}/{total})"
        );

        let mut changes = 0usize;
        let mut total = 0usize;
        for _ in 0..500 {
            if let Some(mv) =
                random_move(&sol, &inc, &mut rng, MoveKind::Mixed { change_ratio: 1.0 })
            {
                assert_legal(&sol, &mv);
                total += 1;
                if matches!(mv, Move::Change { .. }) {
                    changes += 1;
                }
            }
        }
        assert!(
            changes * 100 > total * 90,
            "ratio 1 should overwhelmingly produce Changes (got {changes}/{total})"
        );
    }

    /// Change 与 Swap 的 apply 均成对同步 sol + inc，preview 恒等于应用后分数。
    #[test]
    fn apply_keeps_solution_and_shadow_in_sync() {
        let model = toy_model();
        let (mut sol, mut inc) = setup(&model, &[3, 0, 3, 2, 0, 2]);
        let calc = ScoreCalculator::new(model);

        let mv = Move::Change { entity: 0, new_value: 1 };
        let preview = mv.preview(&inc);
        mv.apply(&mut sol, &mut inc);
        assert_eq!(sol.value(0), 1);
        assert_eq!(inc.value(0), 1);
        assert_eq!(inc.score(), preview);
        assert_eq!(inc.score(), calc.total_score(&sol));

        // v0=1, v1=0：交换后 e0=0（域内）、e1=1（域内）
        let mv = Move::Swap { a: 0, b: 1 };
        let preview = mv.preview(&inc);
        mv.apply(&mut sol, &mut inc);
        assert_eq!(sol.value(0), 0);
        assert_eq!(sol.value(1), 1);
        assert_eq!(inc.value(0), sol.value(0));
        assert_eq!(inc.value(1), sol.value(1));
        assert_eq!(inc.score(), preview);
        assert_eq!(inc.score(), calc.total_score(&sol));
    }

    /// Swap 前置条件违例（双方已赋值）：assert + panic，与 apply_swap 同语义。
    #[test]
    fn apply_swap_on_unassigned_panics() {
        let model = toy_model();
        let mut sol = Solution::empty(model.clone());
        let mut inc = IncrementalScorer::new(model);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            Move::Swap { a: 0, b: 1 }.apply(&mut sol, &mut inc);
        }));
        assert!(result.is_err(), "Move::apply Swap on unassigned must panic");
    }

    #[test]
    fn entities_pairs_and_counts() {
        let c = Move::Change { entity: 3, new_value: 1 };
        assert_eq!(c.entities(), ([3, 3], 1));
        let s = Move::Swap { a: 0, b: 5 };
        assert_eq!(s.entities(), ([0, 5], 2));
    }

    #[test]
    fn move_kind_default_is_mixed_half() {
        assert_eq!(MoveKind::default(), MoveKind::Mixed { change_ratio: 0.5 });
    }
}
