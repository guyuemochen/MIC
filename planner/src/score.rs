//! 评分器：全量重算（[`ScoreCalculator`]，真值）与增量评分（[`IncrementalScorer`]，影子状态）。
//!
//! ## 增量一致性保证
//! 1. 影子状态（赋值 / 负载表 / 当前分）是 `(model, assignment)` 的纯函数；
//!    变更只能经 `apply_*` 同步（[`crate::model::Solution::set`] 为 `pub(crate)` 且仅被
//!    [`crate::construct::construct`] 与 [`crate::moves::Move::apply`] 成对调用）。
//! 2. `preview_*` 与 `apply_*` 共用同一私有增量内核（`change_terms` / `swap_terms`），
//!    数学上 preview − score == 全量重算差（所有评分项只经受影响格 / 受影响冲突边依赖该实体）。
//! 3. 兜底 = 验收测试的随机游走属性测试（每步断言 preview == 应用后影子分 == 全量重算）。

use std::sync::Arc;

use crate::model::{EntityId, HardSoftScore, ProblemModel, Solution, ValueId, UNASSIGNED};

/// 全量评分器（真值）。对任意解（含未分配变量的不完整解）给出正确总分，O(实体·资源 + 冲突边)。
///
/// hard = Σ_{r,v} max(0, load[r][v] − cap[r][v]) + |{(e,c) ∈ 冲突边 : val(e) = val(c) ≠ UNASSIGNED}|
/// soft = Σ_e cost[e][val(e)]；未分配实体不计任何项。
pub struct ScoreCalculator {
    model: Arc<ProblemModel>,
}

impl ScoreCalculator {
    /// 构造全量评分器。
    pub fn new(model: Arc<ProblemModel>) -> ScoreCalculator {
        ScoreCalculator { model }
    }

    /// 对任意解给出正确总分（含未分配变量的不完整解）。
    pub fn total_score(&self, sol: &Solution) -> HardSoftScore {
        let m = self.model.value_count;
        let n = self.model.entity_count();
        let mut hard: i64 = 0;
        let mut soft: i64 = 0;
        // 负载类硬违反：逐资源重建整张负载表
        for (r, res) in self.model.resources.iter().enumerate() {
            let mut load = vec![0i64; m];
            for e in 0..n {
                let v = sol.value(e);
                if v != UNASSIGNED {
                    load[v as usize] += self.model.entities[e].demands[r];
                }
            }
            for v in 0..m {
                hard += (load[v] - res.capacity[v]).max(0);
            }
        }
        // 冲突硬违反：每条无向边恰计一次（conf 邻接对称，取 c > e 一侧）
        for e in 0..n {
            let ve = sol.value(e);
            if ve == UNASSIGNED {
                continue;
            }
            for &c in &self.model.entities[e].conflicts {
                if c > e && sol.value(c) == ve {
                    hard += 1;
                }
            }
        }
        // 软代价：未分配实体不计
        for e in 0..n {
            let v = sol.value(e);
            if v != UNASSIGNED {
                soft += self.model.cost[e][v as usize];
            }
        }
        HardSoftScore::new(hard, soft)
    }

    /// 所属模型。
    pub fn model(&self) -> &ProblemModel {
        &self.model
    }
}

/// 增量评分器：与解同步演化的影子状态（赋值、每资源每值负载、当前 hard/soft）。
pub struct IncrementalScorer {
    model: Arc<ProblemModel>,
    value: Vec<ValueId>,
    load: Vec<Vec<i64>>,
    hard: i64,
    soft: i64,
}

impl IncrementalScorer {
    /// 全 [`crate::model::UNASSIGNED`] 状态，`score() == HardSoftScore::ZERO`。
    pub fn new(model: Arc<ProblemModel>) -> IncrementalScorer {
        IncrementalScorer {
            load: vec![vec![0; model.value_count]; model.resources.len()],
            value: vec![UNASSIGNED; model.entity_count()],
            hard: 0,
            soft: 0,
            model,
        }
    }

    /// 重建影子状态并返回总分；返回值必须等于
    /// [`ScoreCalculator::total_score`] 对同一解的结果。
    pub fn reset(&mut self, sol: &Solution) -> HardSoftScore {
        self.value = sol.assignment().to_vec();
        let m = self.model.value_count;
        for (r, _res) in self.model.resources.iter().enumerate() {
            let mut load = vec![0i64; m];
            for e in 0..self.model.entity_count() {
                let v = self.value[e];
                if v != UNASSIGNED {
                    load[v as usize] += self.model.entities[e].demands[r];
                }
            }
            self.load[r] = load;
        }
        self.hard = 0;
        self.soft = 0;
        // 负载违反直接读影子负载表
        for (r, res) in self.model.resources.iter().enumerate() {
            for v in 0..m {
                self.hard += (self.load[r][v] - res.capacity[v]).max(0);
            }
        }
        // 冲突违反（无向边一次）
        for e in 0..self.model.entity_count() {
            let ve = self.value[e];
            if ve == UNASSIGNED {
                continue;
            }
            for &c in &self.model.entities[e].conflicts {
                if c > e && self.value[c] == ve {
                    self.hard += 1;
                }
            }
        }
        // 软代价
        for e in 0..self.model.entity_count() {
            let v = self.value[e];
            if v != UNASSIGNED {
                self.soft += self.model.cost[e][v as usize];
            }
        }
        self.score()
    }

    /// 当前总分。
    pub fn score(&self) -> HardSoftScore {
        HardSoftScore::new(self.hard, self.soft)
    }

    /// 影子赋值读取（与同步的 [`Solution::value`] 一致）。
    pub fn value(&self, e: EntityId) -> ValueId {
        self.value[e]
    }

    /// 预演：把实体 `e` 改为 `new` 后的**新总分**（不是差值，杜绝调用方加减漂移）。
    /// `old`/`new` 任一可为 [`crate::model::UNASSIGNED`]（未分配→分配、分配→未分配均合法）。
    pub fn preview_change(&self, e: EntityId, new: ValueId) -> HardSoftScore {
        let (dh, ds) = self.change_terms(e, new);
        HardSoftScore::new(self.hard + dh, self.soft + ds)
    }

    /// 预演：`a`、`b` 互换取值后的新总分（**联合计算**，不许拆成两次 change——
    /// a、b 互为冲突邻居时拆开会算错）。
    ///
    /// 前置条件（显式契约）：`a`、`b` 均已分配（value ≠ UNASSIGNED）。违例属调用方逻辑错误，
    /// 实现必须 `assert!(..., "preview_swap on unassigned entity")` 立即 panic（禁止静默返回错值）。
    /// 这不是"用户输入"路径：经 parse/construct/local_search 的正常数据流不可达。
    pub fn preview_swap(&self, a: EntityId, b: EntityId) -> HardSoftScore {
        assert!(
            self.value[a] != UNASSIGNED && self.value[b] != UNASSIGNED,
            "preview_swap on unassigned entity"
        );
        let (dh, ds) = self.swap_terms(a, b);
        HardSoftScore::new(self.hard + dh, self.soft + ds)
    }

    /// 提交 change：同步更新影子赋值 / 负载 / 总分（与 [`Self::preview_change`] 同一内核）。
    pub fn apply_change(&mut self, e: EntityId, new: ValueId) {
        let (dh, ds) = self.change_terms(e, new);
        let old = self.value[e];
        if old == new {
            return;
        }
        for r in 0..self.model.resources.len() {
            let d = self.model.entities[e].demands[r];
            if old != UNASSIGNED {
                self.load[r][old as usize] -= d;
            }
            if new != UNASSIGNED {
                self.load[r][new as usize] += d;
            }
        }
        self.value[e] = new;
        self.hard += dh;
        self.soft += ds;
    }

    /// 提交 swap：前置条件与 panic 语义同 [`Self::preview_swap`]（已分配断言）。
    pub fn apply_swap(&mut self, a: EntityId, b: EntityId) {
        assert!(
            self.value[a] != UNASSIGNED && self.value[b] != UNASSIGNED,
            "apply_swap on unassigned entity"
        );
        let (dh, ds) = self.swap_terms(a, b);
        let (va, vb) = (self.value[a], self.value[b]);
        if va == vb {
            return;
        }
        for r in 0..self.model.resources.len() {
            let da = self.model.entities[a].demands[r];
            let db = self.model.entities[b].demands[r];
            self.load[r][va as usize] += db - da;
            self.load[r][vb as usize] += da - db;
        }
        self.value[a] = vb;
        self.value[b] = va;
        self.hard += dh;
        self.soft += ds;
    }

    // ---- 私有增量内核（preview_* 与 apply_* 共用，数值同源；UNASSIGNED 不占负载不计代价） ----

    /// (Δhard, Δsoft)：把 `e` 从当前值改为 `new`。
    ///
    /// Δhard_load = Σ_r [max(0, load[new]+d−cap[new]) − max(0, load[new]−cap[new])]
    /// （old ≠ UNA 时对称减 old 格的一项）；Δhard_conf =
    /// |{c ∈ conf(e) : val(c) == new}| − (old ≠ UNA ? |{c : val(c) == old}| : 0)；
    /// Δsoft = cost[e][new] − (old ≠ UNA ? cost[e][old] : 0)。
    ///
    /// `new == old`（含双双 UNASSIGNED）为 no-op 返回 (0, 0)：remove+add 作用于同一格的
    /// 净效应恒为零，而分格公式在过载边界处不具备该线性性，故显式短路。
    fn change_terms(&self, e: EntityId, new: ValueId) -> (i64, i64) {
        let old = self.value[e];
        if old == new {
            return (0, 0);
        }
        let ent = &self.model.entities[e];
        let mut dh: i64 = 0;
        for (r, res) in self.model.resources.iter().enumerate() {
            let d = ent.demands[r];
            let load = &self.load[r];
            let cap = &res.capacity;
            if new != UNASSIGNED {
                let vi = new as usize;
                dh += (load[vi] + d - cap[vi]).max(0) - (load[vi] - cap[vi]).max(0);
            }
            if old != UNASSIGNED {
                let vi = old as usize;
                dh += (load[vi] - d - cap[vi]).max(0) - (load[vi] - cap[vi]).max(0);
            }
        }
        if new != UNASSIGNED {
            for &c in &ent.conflicts {
                if self.value[c] == new {
                    dh += 1;
                }
            }
        }
        if old != UNASSIGNED {
            for &c in &ent.conflicts {
                if self.value[c] == old {
                    dh -= 1;
                }
            }
        }
        let ds = if new != UNASSIGNED {
            self.model.cost[e][new as usize]
        } else {
            0
        } - if old != UNASSIGNED {
            self.model.cost[e][old as usize]
        } else {
            0
        };
        (dh, ds)
    }

    /// (Δhard, Δsoft)：`a`、`b` 互换（va ≠ vb；va == vb 为 no-op，返回 (0, 0)）。
    ///
    /// Δhard_load 对两格 {va, vb} 各做 −da+db / −db+da；Δhard_conf =
    /// Σ_{c ∈ conf(a)\{b}} ([vc==vb] − [vc==va]) + Σ_{c ∈ conf(b)\{a}} ([vc==va] − [vc==vb])
    /// （a、b 互为冲突时该项为 0，因 va ≠ vb）；
    /// Δsoft = cost[a][vb] + cost[b][va] − cost[a][va] − cost[b][vb]。
    ///
    /// 前置条件：`a`、`b` 均已分配（由 pub 接口断言）。
    fn swap_terms(&self, a: EntityId, b: EntityId) -> (i64, i64) {
        let (va, vb) = (self.value[a], self.value[b]);
        if va == vb {
            return (0, 0);
        }
        let (ea, eb) = (&self.model.entities[a], &self.model.entities[b]);
        let mut dh: i64 = 0;
        for (r, res) in self.model.resources.iter().enumerate() {
            let (da, db) = (ea.demands[r], eb.demands[r]);
            let load = &self.load[r];
            let cap = &res.capacity;
            let (ia, ib) = (va as usize, vb as usize);
            dh += (load[ia] - da + db - cap[ia]).max(0) - (load[ia] - cap[ia]).max(0);
            dh += (load[ib] - db + da - cap[ib]).max(0) - (load[ib] - cap[ib]).max(0);
        }
        for &c in &ea.conflicts {
            if c == b {
                continue;
            }
            let vc = self.value[c];
            if vc == vb {
                dh += 1;
            } else if vc == va {
                dh -= 1;
            }
        }
        for &c in &eb.conflicts {
            if c == a {
                continue;
            }
            let vc = self.value[c];
            if vc == va {
                dh += 1;
            } else if vc == vb {
                dh -= 1;
            }
        }
        let ds = self.model.cost[a][vb as usize]
            + self.model.cost[b][va as usize]
            - self.model.cost[a][va as usize]
            - self.model.cost[b][vb as usize];
        (dh, ds)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{PlanningEntity, Resource};
    use crate::rng::Rng;

    /// 手造小模型（不经 parse，直接程序化构建）：
    /// m=4 值、2 资源、6 实体、冲突链 + 三角、非平凡代价与需求。
    fn toy_model() -> Arc<ProblemModel> {
        let cap_sentinel = i64::MAX / 2;
        let resources = vec![
            Resource {
                name: "cpu".to_string(),
                capacity: vec![2, 2, 3, cap_sentinel],
            },
            Resource {
                name: "mem".to_string(),
                capacity: vec![3, 1, 2, 2],
            },
        ];
        let dom = |d: &[u32]| d.to_vec();
        let ent = |dom: Vec<u32>, demands: Vec<i64>, conflicts: Vec<usize>| PlanningEntity {
            domain: dom,
            demands,
            conflicts,
        };
        let entities = vec![
            ent(dom(&[0, 1, 2, 3]), vec![2, 1], vec![1, 2]),      // 0
            ent(dom(&[0, 1, 3]), vec![1, 2], vec![0, 3]),          // 1
            ent(dom(&[1, 2, 3]), vec![3, 1], vec![0, 3, 4]),       // 2
            ent(dom(&[0, 2, 3]), vec![1, 1], vec![1, 2]),          // 3
            ent(dom(&[0, 1, 2, 3]), vec![2, 3], vec![2, 5]),       // 4
            ent(dom(&[0, 2, 3]), vec![1, 1], vec![4]),             // 5
        ];
        let cost = vec![
            vec![4, 1, 7, 0],
            vec![2, 9, 0, 3],
            vec![0, 5, 6, 2],
            vec![8, 0, 1, 4],
            vec![3, 2, 8, 5],
            vec![6, 0, 4, 7],
        ];
        let model = ProblemModel {
            name: "toy".to_string(),
            value_count: 4,
            resources,
            entities,
            cost,
        };
        model.validate().expect("toy model must be valid");
        Arc::new(model)
    }

    #[test]
    fn fresh_scorer_is_zero_and_total_agrees_on_empty() {
        let model = toy_model();
        let sol = Solution::empty(model.clone());
        let calc = ScoreCalculator::new(model.clone());
        let inc = IncrementalScorer::new(model);
        assert_eq!(inc.score(), HardSoftScore::ZERO);
        assert_eq!(calc.total_score(&sol), HardSoftScore::ZERO);
        assert!(calc.total_score(&sol).is_feasible());
    }

    #[test]
    fn hand_computed_total_score() {
        let model = toy_model();
        let calc = ScoreCalculator::new(model.clone());
        // 赋值 [3, 0, 3, 2, 0, 2]
        let sol = Solution::from_assignment(model.clone(), &[3, 0, 3, 2, 0, 2]).expect("in domain");
        // 负载（cpu, mem）：值 0: e1(1,2)+e4(2,3) = (3,5)；值 2: e3(1,1)+e5(1,1) = (2,2)；
        // 值 3: e0(2,1)+e2(3,1) = (5,2)
        // cpu caps [2,2,3,SENTINEL]: 值0 → max(0,3-2)=1；值2 → max(0,2-3)=0；值3 → 哨兵 → 0 ⇒ 1
        // mem caps [3,1,2,2]: 值0 → max(0,5-3)=2；值2 → 0；值3 → 0 ⇒ 2
        // 冲突边: (0,2): 3==3 → 1；其余边均不同值 ⇒ 1
        // soft = cost[0][3]+cost[1][0]+cost[2][3]+cost[3][2]+cost[4][0]+cost[5][2]
        //      = 0+2+2+1+3+4 = 12
        assert_eq!(calc.total_score(&sol), HardSoftScore::new(1 + 2 + 1, 12));
        // 部分解（e5 未分配）：其负载/代价/冲突边全部消失
        let sol = Solution::from_assignment(model, &[3, 0, 3, 2, 0, UNASSIGNED]).expect("in domain");
        // cpu: 值0 = 3 → 1；值2 = e3 = 1 → 0；值3 → 哨兵 → 0 ⇒ 1
        // mem: 值0 = 5 → 2；其余 0 ⇒ 2；冲突 (0,2) 仍命中 ⇒ 1
        // soft = 12 - cost[5][2] = 8
        assert_eq!(calc.total_score(&sol), HardSoftScore::new(1 + 2 + 1, 8));
    }

    #[test]
    fn reset_returns_total_and_shadow_matches() {
        let model = toy_model();
        let calc = ScoreCalculator::new(model.clone());
        let mut inc = IncrementalScorer::new(model.clone());
        let sol = Solution::from_assignment(model, &[1, 3, 1, 0, 2, 3]).expect("in domain");
        let s = inc.reset(&sol);
        assert_eq!(s, calc.total_score(&sol));
        assert_eq!(inc.score(), s);
        for e in 0..sol.model().entity_count() {
            assert_eq!(inc.value(e), sol.value(e));
        }
    }

    #[test]
    fn preview_swap_on_unassigned_panics() {
        let model = toy_model();
        let mut inc = IncrementalScorer::new(model);
        inc.apply_change(0, 1); // e0 已分配，其余未分配
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = inc.preview_swap(0, 1); // e1 未分配
        }));
        assert!(result.is_err(), "preview_swap on unassigned must panic");
        // 消息逐字匹配契约（assert! 字面量消息为 &str payload）
        let payload = result.err().expect("panic payload");
        let msg = payload
            .downcast_ref::<String>()
            .map(|s| s.as_str())
            .or_else(|| payload.downcast_ref::<&str>().copied());
        assert_eq!(
            msg,
            Some("preview_swap on unassigned entity"),
            "panic message must match contract, got {msg:?}"
        );
        // apply_swap 同语义
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            inc.apply_swap(0, 2);
        }));
        assert!(result.is_err(), "apply_swap on unassigned must panic");
    }

    #[test]
    #[should_panic(expected = "preview_swap on unassigned entity")]
    fn preview_swap_unassigned_panic_attribute_form() {
        let model = toy_model();
        let inc = IncrementalScorer::new(model);
        let _ = inc.preview_swap(0, 1); // 全未分配 ⇒ 违反前置条件
    }

    /// 随机游走一致性：每步随机 change / UNASSIGN / swap，断言
    /// (1) preview == 假想解的全量重算；(2) apply 后影子分 == 全量重算；
    /// (3) 影子赋值与解逐位一致。含 UNASSIGNED 中间态。
    #[test]
    fn random_walk_incremental_equals_full_recompute() {
        let model = toy_model();
        let n = model.entity_count();
        let calc = ScoreCalculator::new(model.clone());
        let mut sol = Solution::empty(model.clone());
        let mut inc = IncrementalScorer::new(model.clone());
        let mut rng = Rng::new(2026_1007);
        assert_eq!(inc.score(), calc.total_score(&sol));

        let hypothetical = |sol: &Solution, e: usize, v: ValueId, calc: &ScoreCalculator| {
            let mut tmp = sol.clone();
            tmp.set(e, v);
            calc.total_score(&tmp)
        };

        for step in 0..6000 {
            let pick = rng.gen_range(0, 10);
            if pick < 6 {
                // change（含随机未分配化）
                let e = rng.gen_range(0, n);
                let v = if rng.gen_ratio(0.9) {
                    let dom = sol.model().domain(e);
                    dom[rng.gen_range(0, dom.len())]
                } else {
                    UNASSIGNED
                };
                let expected = hypothetical(&sol, e, v, &calc);
                assert_eq!(
                    inc.preview_change(e, v),
                    expected,
                    "step {step}: preview_change({e}, {v}) mismatch"
                );
                sol.set(e, v);
                inc.apply_change(e, v);
            } else {
                // swap：取两个已分配实体（不足则跳过该步）
                let assigned: Vec<usize> = (0..n).filter(|&e| sol.value(e) != UNASSIGNED).collect();
                if assigned.len() < 2 {
                    continue;
                }
                let a = assigned[rng.gen_range(0, assigned.len())];
                let b = assigned[rng.gen_range(0, assigned.len())];
                if a == b || sol.value(a) == sol.value(b) {
                    continue;
                }
                let (va, vb) = (sol.value(a), sol.value(b));
                let mut tmp = sol.clone();
                tmp.set(a, vb);
                tmp.set(b, va);
                let expected = calc.total_score(&tmp);
                assert_eq!(
                    inc.preview_swap(a, b),
                    expected,
                    "step {step}: preview_swap({a}, {b}) mismatch"
                );
                inc.apply_swap(a, b);
                sol.set(a, vb);
                sol.set(b, va);
            }
            assert_eq!(inc.score(), calc.total_score(&sol), "step {step}: shadow score drift");
            for e in 0..n {
                assert_eq!(inc.value(e), sol.value(e), "step {step}: shadow value drift at {e}");
            }
        }
    }

    /// 专门覆盖 swap 的"互为冲突邻居"情形：拆成两次 change 会算错，联合计算必须正确。
    #[test]
    fn swap_between_mutual_conflict_neighbors() {
        let model = toy_model();
        let calc = ScoreCalculator::new(model.clone());
        let mut sol = Solution::empty(model.clone());
        let mut inc = IncrementalScorer::new(model.clone());
        // e0 与 e1 互为冲突；同置 1（产生 1 条违反），再与其他实体交织
        inc.apply_change(0, 1);
        sol.set(0, 1);
        inc.apply_change(1, 1);
        sol.set(1, 1);
        inc.apply_change(2, 3);
        sol.set(2, 3);
        assert_eq!(inc.score(), calc.total_score(&sol));
        // e1 <-> e2（域约束满足：e2 的 3? e1=1,e2=3 → swap 后 e1=3,e2=1，均在域内）
        let expected = {
            let mut tmp = sol.clone();
            tmp.set(1, 3);
            tmp.set(2, 1);
            calc.total_score(&tmp)
        };
        assert_eq!(inc.preview_swap(1, 2), expected);
        // e0 <-> e1 互为冲突（同值分支已被上一步换走；此处 va≠vb）
        inc.apply_swap(1, 2);
        sol.set(1, 3);
        sol.set(2, 1);
        assert_eq!(inc.score(), calc.total_score(&sol));
        // 构造互为冲突且 va≠vb 的直接 swap：e0=1 与 e1=3
        let expected = {
            let mut tmp = sol.clone();
            tmp.set(0, 3);
            tmp.set(1, 1);
            calc.total_score(&tmp)
        };
        assert_eq!(inc.preview_swap(0, 1), expected);
        inc.apply_swap(0, 1);
        sol.set(0, 3);
        sol.set(1, 1);
        assert_eq!(inc.score(), calc.total_score(&sol));
    }

    /// change 到相同值 = no-op（preview 与 apply 都恒等于当前分）。
    #[test]
    fn change_to_same_value_is_noop() {
        let model = toy_model();
        let mut sol = Solution::empty(model.clone());
        let mut inc = IncrementalScorer::new(model.clone());
        let calc = ScoreCalculator::new(model);
        inc.apply_change(3, 2);
        sol.set(3, 2);
        let before = inc.score();
        assert_eq!(inc.preview_change(3, 2), before);
        inc.apply_change(3, 2);
        assert_eq!(inc.score(), before);
        assert_eq!(inc.score(), calc.total_score(&sol));
        // UNASSIGNED -> UNASSIGNED 同样 no-op
        assert_eq!(inc.preview_change(5, UNASSIGNED), before);
        inc.apply_change(5, UNASSIGNED);
        assert_eq!(inc.score(), before);
    }
}
