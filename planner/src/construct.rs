//! 构造启发式：First Fit / Best Fit（确定性：同模型同启发式 ⇒ 同解，不使用随机数）。

use std::sync::Arc;

use crate::model::{PlannerError, ProblemModel, Solution, ValueId};
use crate::score::IncrementalScorer;

/// 构造启发式种类。
///
/// 语义（normative，平局规则决定确定性）：逐实体 e（按下标 0..n）考察其值域（升序），
/// 令 dh(v) = 该赋值新增的硬违反（= `preview_change(e, v).hard − 当前 hard`），
/// key(v) = `(dh(v), cost[e][v])`：
/// - [`FirstFit`](ConstructionHeuristic::FirstFit)：取值域中第一个 dh(v) == 0 的值；
///   若无，取 key 字典序最小者（平局取值域序最前）。
/// - [`BestFit`](ConstructionHeuristic::BestFit)：直接取 key 字典序最小者
///   （平局取值域序最前）。
///
/// 选定后立即 `IncrementalScorer::apply_change` 提交再处理下一实体。构造过程允许中间解
/// 带硬违反（如 cloud-small 的 FirstFit 起点为 hard=6，由局部搜索修复）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstructionHeuristic {
    /// First Fit：第一个不新增硬违反的值，否则回退 (dh, cost) 字典序最小。
    FirstFit,
    /// Best Fit：直接取 (dh, cost) 字典序最小。
    BestFit,
}

/// 构造完整初始解。
///
/// 先 `model.validate()`：模型非法（如空值域）返回 `Err(PlannerError::Model(..))`——
/// 不对空迭代器取首元素而 panic；合法模型恒成功，返回 `is_complete() == true` 的解。
/// 内部用 [`crate::score::IncrementalScorer`] 的 `preview_change`/`apply_change` 驱动
/// （与求解器共用同一条增量路径，天然一致）。
pub fn construct(
    model: &Arc<ProblemModel>,
    heuristic: ConstructionHeuristic,
) -> Result<Solution, PlannerError> {
    model.validate()?;
    let n = model.entity_count();
    let mut sol = Solution::empty(model.clone());
    let mut inc = IncrementalScorer::new(model.clone());
    for e in 0..n {
        let cur_hard = inc.score().hard;
        let mut first_feasible: Option<ValueId> = None; // FF：第一个 dh == 0 的值
        let mut best: Option<(ValueId, (i64, i64))> = None; // key 字典序最小（严格小于才替换 → 平局取前）
        for &v in model.domain(e) {
            let dh = inc.preview_change(e, v).hard - cur_hard;
            let key = (dh, model.cost[e][v as usize]);
            if first_feasible.is_none() && dh == 0 {
                first_feasible = Some(v);
            }
            if best.map_or(true, |(_, bk)| key < bk) {
                best = Some((v, key));
            }
        }
        let chosen = if heuristic == ConstructionHeuristic::FirstFit {
            first_feasible.or(best.map(|(v, _)| v))
        } else {
            best.map(|(v, _)| v)
        };
        // validate() 已保证值域非空；此处仍以 Err 兜底而非 unwrap（双保险路径）
        let Some(v) = chosen else {
            return Err(PlannerError::Model(format!("entity {e} has empty domain")));
        };
        sol.set(e, v);
        inc.apply_change(e, v);
    }
    Ok(sol)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{HardSoftScore, PlanningEntity, Resource};
    use crate::score::ScoreCalculator;

    /// FF/BF 选择序列锁定用的小模型：
    /// m=2 值、3 实体两两冲突（三角）、cost 全 1（平局）、无资源。
    fn triangle_model() -> Arc<ProblemModel> {
        let model = ProblemModel {
            name: "triangle".to_string(),
            value_count: 2,
            resources: vec![],
            entities: vec![
                PlanningEntity { domain: vec![0, 1], demands: vec![], conflicts: vec![1, 2] },
                PlanningEntity { domain: vec![0, 1], demands: vec![], conflicts: vec![0, 2] },
                PlanningEntity { domain: vec![0, 1], demands: vec![], conflicts: vec![0, 1] },
            ],
            cost: vec![vec![1, 1], vec![1, 1], vec![1, 1]],
        };
        model.validate().expect("triangle model must be valid");
        Arc::new(model)
    }

    /// 带资源与差异代价的小模型：m=3、cap 糖语义手工重建、域受限。
    fn resource_model() -> Arc<ProblemModel> {
        let model = ProblemModel {
            name: "res".to_string(),
            value_count: 3,
            resources: vec![Resource {
                name: "cap".to_string(),
                capacity: vec![1, 2, i64::MAX / 2],
            }],
            entities: vec![
                PlanningEntity { domain: vec![0, 1, 2], demands: vec![1], conflicts: vec![] },
                PlanningEntity { domain: vec![0, 1, 2], demands: vec![1], conflicts: vec![] },
                PlanningEntity { domain: vec![0, 2], demands: vec![1], conflicts: vec![] },
            ],
            cost: vec![vec![1, 5, 2], vec![4, 1, 6], vec![7, 0, 3]],
        };
        model.validate().expect("resource model must be valid");
        Arc::new(model)
    }

    #[test]
    fn ff_tie_breaks_to_first_and_falls_back() {
        // 三角冲突 + 全 1 代价：e0 取首个 dh==0 的值 0；e1 取 1；e2 无 dh==0 值 →
        // FF 回退 key 最小（dh=1 平局、cost 平局 → 值域序最前 = 0）
        let model = triangle_model();
        let sol = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        assert_eq!(sol.assignment(), &[0, 1, 0]);
        // 中间允许硬违反：e2 与 e0 同值 → hard = 1（三角不可能完全可行，2 值 3 实体）
        let calc = ScoreCalculator::new(model.clone());
        assert_eq!(calc.total_score(&sol), HardSoftScore::new(1, 3));
    }

    #[test]
    fn bf_picks_min_key_sequence() {
        // BF 对 e0/e1/e2 均直接取 key 最小：e0 → (0,1) 平局取前 = 0；
        // e1 → dh(0)=1/cost4, dh(1)=0/cost1, dh(2)=0/cost6 → 值 1；
        // e2 → dh(0)=1/cost7, dh(2)=0/cost3 → 值 2（注意 e2 域 = [0,2]）
        let model = resource_model();
        let sol = construct(&model, ConstructionHeuristic::BestFit).expect("construct");
        assert_eq!(sol.assignment(), &[0, 1, 2]);
        let calc = ScoreCalculator::new(model.clone());
        // soft = 1 + 1 + 3 = 5；负载：值0 = 1 ≤ 1、值1 = 1 ≤ 2、值2 = 1（哨兵）→ hard 0
        assert_eq!(calc.total_score(&sol), HardSoftScore::new(0, 5));
    }

    #[test]
    fn ff_takes_first_feasible_while_bf_takes_cheapest() {
        // 单实体、无约束、cost = [5, 1]：FF 取首个 dh==0 的值 0（soft=5），
        // BF 取 key 最小的值 1（soft=1）——两启发式在此明确分叉
        let model = Arc::new(ProblemModel {
            name: "cheap".to_string(),
            value_count: 2,
            resources: vec![],
            entities: vec![PlanningEntity {
                domain: vec![0, 1],
                demands: vec![],
                conflicts: vec![],
            }],
            cost: vec![vec![5, 1]],
        });
        model.validate().expect("valid");
        let ff = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        let bf = construct(&model, ConstructionHeuristic::BestFit).expect("construct");
        assert_eq!(ff.assignment(), &[0]);
        assert_eq!(bf.assignment(), &[1]);
        let calc = ScoreCalculator::new(model);
        assert_eq!(calc.total_score(&ff), HardSoftScore::new(0, 5));
        assert_eq!(calc.total_score(&bf), HardSoftScore::new(0, 1));
    }

    #[test]
    fn ff_skips_overloaded_first_value() {
        // 资源约束下 FF 的"首个可行"逐实体推进：
        // e0 → 首个 dh==0 = 0；e1 → 值 0 会令 load[0]=2>1（dh=1），首个 dh==0 是 1；
        // e2（域 [0,2]）→ 值 0 超（dh=1）、值 2 → dh=0 → 取 2
        let model = resource_model();
        let sol = construct(&model, ConstructionHeuristic::FirstFit).expect("construct");
        assert_eq!(sol.assignment(), &[0, 1, 2]);
        let calc = ScoreCalculator::new(model.clone());
        assert_eq!(calc.total_score(&sol), HardSoftScore::new(0, 5));
    }

    #[test]
    fn deterministic_same_model_same_assignment() {
        for model in [triangle_model(), resource_model()] {
            for h in [ConstructionHeuristic::FirstFit, ConstructionHeuristic::BestFit] {
                let a = construct(&model, h).expect("construct");
                let b = construct(&model, h).expect("construct");
                assert_eq!(a.assignment(), b.assignment(), "construct must be deterministic");
            }
        }
    }

    #[test]
    fn invalid_model_returns_err_not_panic() {
        // 空值域模型（程序化构建，绕过 parse）：construct 必须 Err 而非 panic
        let model = Arc::new(ProblemModel {
            name: "bad".to_string(),
            value_count: 2,
            resources: vec![],
            entities: vec![PlanningEntity {
                domain: vec![], // 空
                demands: vec![],
                conflicts: vec![],
            }],
            cost: vec![vec![0, 0]],
        });
        match construct(&model, ConstructionHeuristic::FirstFit) {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("empty domain"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }
        match construct(&model, ConstructionHeuristic::BestFit) {
            Err(PlannerError::Model(_)) => {}
            other => panic!("expected Model error, got {other:?}"),
        }
    }

    #[test]
    fn result_is_complete_and_in_domain() {
        for model in [triangle_model(), resource_model()] {
            for h in [ConstructionHeuristic::FirstFit, ConstructionHeuristic::BestFit] {
                let sol = construct(&model, h).expect("construct");
                assert!(sol.is_complete(), "constructed solution must be complete");
                for e in 0..model.entity_count() {
                    assert!(
                        model.in_domain(e, sol.value(e)),
                        "entity {e} value {} must be in domain",
                        sol.value(e)
                    );
                }
            }
        }
    }

    #[test]
    fn incremental_final_score_equals_total() {
        for model in [triangle_model(), resource_model()] {
            for h in [ConstructionHeuristic::FirstFit, ConstructionHeuristic::BestFit] {
                let sol = construct(&model, h).expect("construct");
                let calc = ScoreCalculator::new(model.clone());
                let mut inc = IncrementalScorer::new(model.clone());
                let s = inc.reset(&sol);
                assert_eq!(s, calc.total_score(&sol), "reset must equal total (h={h:?})");
            }
        }
    }

    /// 内嵌文本算例（经 format::parse 全链路），FF 与 BF 构造分人工推演如下：
    ///
    /// 模型：m=3；资源 "cap"（CAP 0 1 ⇒ 容量 [1, 2, i64::MAX/2]）；全员 demand 1；
    /// 冲突 0–2；cost: e0=[2,9,4]，e1=[7,6,1]（COSTROW 全 0 后被 COST 覆盖），e2=[0,0,0]；
    /// 域：e0/e1 全值，e2=[0,2]（FORBID 2 1）。
    ///
    /// **FF**：e0 全空 → 首个 dh==0 取 v0（cost 2）；
    /// e1：dh(v0)=1（load 1+1>cap 1），dh(v1)=0，dh(v2)=0 → 首个可行 v1（cost 6）；
    /// e2：dh(v0)=2（超载 1 + 与 e0 冲突 1），dh(v2)=0 → v2。
    /// ⇒ 赋值 [0,1,2]，hard=0，soft = 2+6+0 = **8**。
    ///
    /// **BF**：e0 key 最小 (0,2)@v0；e1 keys：(1,7)@v0、(0,6)@v1、(0,1)@v2 → v2（cost 1）；
    /// e2 keys：(2,0)@v0、(0,0)@v2 → v2。
    /// ⇒ 赋值 [0,2,2]，hard=0，soft = 2+1+0 = **3**。
    ///
    /// （两启发式明确分叉：FF 的"首个可行"v1 不及其 key 最小 v2 便宜。）
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

    #[test]
    fn embedded_text_instance_ff_bf_scores() {
        let model = Arc::new(crate::format::parse(HAND_INSTANCE).expect("parse hand instance"));
        let calc = ScoreCalculator::new(model.clone());

        let ff = construct(&model, ConstructionHeuristic::FirstFit).expect("construct ff");
        assert_eq!(ff.assignment(), &[0, 1, 2]);
        assert_eq!(calc.total_score(&ff), HardSoftScore::new(0, 8));

        let bf = construct(&model, ConstructionHeuristic::BestFit).expect("construct bf");
        assert_eq!(bf.assignment(), &[0, 2, 2]);
        assert_eq!(calc.total_score(&bf), HardSoftScore::new(0, 3));

        // 同一文本同启发式重复构造逐位一致（确定性）
        let ff2 = construct(&model, ConstructionHeuristic::FirstFit).expect("construct ff again");
        assert_eq!(ff.assignment(), ff2.assignment());
    }

    #[test]
    fn heuristic_enum_traits() {
        // 契约仅 derive PartialEq/Eq（无 Ord）；Copy/Clone/Debug 可用
        assert_ne!(ConstructionHeuristic::FirstFit, ConstructionHeuristic::BestFit);
        let copy = ConstructionHeuristic::FirstFit;
        assert_eq!(copy, ConstructionHeuristic::FirstFit);
        assert_eq!(format!("{:?}", copy), "FirstFit");
        let cloned = copy.clone();
        assert_eq!(cloned, copy);
    }
}
