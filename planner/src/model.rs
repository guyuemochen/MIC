//! 域模型：规划问题（[`ProblemModel`]）、规划实体（[`PlanningEntity`]）、规划变量值域、
//! 工作解（[`Solution`]）与硬软评分（[`HardSoftScore`]）。
//!
//! OptaPlanner 概念映射：`ProblemModel` ↔ `@PlanningSolution` 的静态部分；
//! `PlanningEntity` ↔ `@PlanningEntity`；`domain` ↔ `@PlanningVariable` 的值域；
//! `Solution` ↔ 工作解（未分配实体记 [`UNASSIGNED`]，对标 nullable planning variable）。
//!
//! 评分语义：`hard = Σ_{r,v} max(0, load[r][v] − cap[r][v])`
//! `+ |{(e,c) ∈ 冲突边 : val(e) = val(c) ≠ UNASSIGNED}|`（可行 ⟺ hard == 0）；
//! `soft = Σ_e cost[e][val(e)]`；未分配实体不计任何项；比较为 `(hard, soft)` 字典序、
//! 越小越好（与 OptaPlanner HardSoftScore 同构，方向相反：彼处最大化负惩罚）。

use std::sync::Arc;

/// 规划实体 id（0 起始，等于 [`ProblemModel::entities`] 的下标）。
pub type EntityId = usize;

/// 规划变量取值 id（0 起始）。保留 [`UNASSIGNED`] 表示"未分配"。
pub type ValueId = u32;

/// 未分配标记。未分配实体不占负载、不计代价、不参与冲突命中。
pub const UNASSIGNED: ValueId = u32::MAX;

/// 规模硬界（validate 与 format::parse 双层强制；求和在 i128 中判定）。
const MAX_VALUES: usize = 1 << 31; // 2^31
const MAX_ENTITIES: usize = 1 << 20; // 2^20
const MAX_RESOURCES: usize = 64;
const MAX_COST_OR_CAP: i128 = 1 << 40; // 单值 cost / cap
const MAX_DEMAND: i128 = 1 << 31; // 单值 demand
const MAX_TOTAL: i128 = 1 << 40; // Σcost / 每资源 Σdemand

/// CAP 糖创建 "cap" 资源时，未被 `CAP v c` 指明值上的默认容量哨兵。
/// 它是唯一的"超 2^40 合法容量"（解析层保证用户书写的容量 ≤ 2^40；
/// 负载算术安全性由 Σdemand ≤ 2^40 硬界保证，见模块级约定）。
pub(crate) const CAP_SENTINEL: i64 = i64::MAX / 2;

/// 硬软评分：`hard` = 硬约束违反总惩罚（0 = 可行），`soft` = 软惩罚总值；越小越好。
///
/// 比较为 `(hard, soft)` 字典序（derive 字段序即所需语义：hard 优先、soft 次之）。
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct HardSoftScore {
    pub hard: i64,
    pub soft: i64,
}

impl HardSoftScore {
    /// 全零评分（空解 / 全未分配解的分数）。
    pub const ZERO: Self = Self { hard: 0, soft: 0 };

    /// 构造评分。
    pub fn new(hard: i64, soft: i64) -> Self {
        Self { hard, soft }
    }

    /// 是否可行（`self.hard == 0`）。
    pub fn is_feasible(&self) -> bool {
        self.hard == 0
    }
}

/// 资源：每个可能取值上的容量（负载类硬约束）。
#[derive(Debug, Clone)]
pub struct Resource {
    pub name: String,
    /// 每值容量，非负；`capacity.len() == ProblemModel::value_count`。
    /// （文本格式中 `CAP` 糖产生的 "cap" 资源，未指明值的容量为 `i64::MAX / 2`。）
    pub capacity: Vec<i64>,
}

/// 规划实体：一个规划变量 + 其值域与静态属性。
#[derive(Debug, Clone)]
pub struct PlanningEntity {
    /// 允许取值：升序、去重、非空（建模期保证；文本格式中 = 全值 [0,m) 减 FORBID）。
    pub domain: Vec<ValueId>,
    /// 对各资源（按 [`ProblemModel::resources`] 声明顺序）的需求数；
    /// 非负；`demands.len() == resources.len()`。
    pub demands: Vec<i64>,
    /// 冲突邻接表（与本实体不得同值的实体）；对称、去重、无自环（建模期保证）。
    pub conflicts: Vec<EntityId>,
}

/// 规划问题模型（不可变；经 [`crate::format::parse`] 或程序化构建，使用前须过 [`ProblemModel::validate`]）。
#[derive(Debug, Clone)]
pub struct ProblemModel {
    pub name: String,
    /// 可能取值个数 m（值 id ∈ [0, m)，且 m ≤ 2^31 保证值 id 与 [`UNASSIGNED`] 不冲突）。
    pub value_count: usize,
    pub resources: Vec<Resource>,
    pub entities: Vec<PlanningEntity>,
    /// `cost[e][v]`：实体 e 取值 v 的软惩罚；非负、缺省 0；`cost.len() == entities.len()`，
    /// `cost[e].len() == value_count`。
    pub cost: Vec<Vec<i64>>,
}

impl ProblemModel {
    /// 实体个数 n。
    pub fn entity_count(&self) -> usize {
        self.entities.len()
    }

    /// 可能取值个数 m。
    pub fn value_count(&self) -> usize {
        self.value_count
    }

    /// 实体 `e` 的值域（升序切片）。
    pub fn domain(&self, e: EntityId) -> &[ValueId] {
        &self.entities[e].domain
    }

    /// `v` 是否在实体 `e` 的值域内（[`UNASSIGNED`] 恒为 false）。
    pub fn in_domain(&self, e: EntityId, v: ValueId) -> bool {
        if v == UNASSIGNED {
            return false;
        }
        self.entities[e].domain.binary_search(&v).is_ok()
    }

    /// 结构校验（值域非空 / 行长度一致 / 冲突对称且无自环）+ 规模硬界校验：
    /// `m ≤ 2^31`、`n ≤ 2^20`、资源数 R ≤ 64、单值 cost/cap ≤ 2^40、单值 demand ≤ 2^31、
    /// `Σ_{e,v} cost ≤ 2^40`、每资源 `Σ_e demand ≤ 2^40`（求和在 i128 中判定）。
    /// 违例返回 `Err(PlannerError::Model(..))`（不 panic；这些界保证评分全程 i64 安全，
    /// 并使 SA 标量 `hard·W + soft` 只需 i128）。
    ///
    /// 容量例外：[`CAP_SENTINEL`]（= i64::MAX/2，CAP 糖的内部默认值）是唯一允许超过
    /// 2^40 的容量值；解析层另行保证用户书写的容量本身 ≤ 2^40。
    pub fn validate(&self) -> Result<(), PlannerError> {
        // ---- 便宜的规模界先行（避免在超界模型上做 O(n·m) 工作） ----
        if self.value_count == 0 {
            return Err(PlannerError::Model("value_count must be >= 1".to_string()));
        }
        if self.value_count > MAX_VALUES {
            return Err(PlannerError::Model(format!(
                "value_count {} exceeds 2^31",
                self.value_count
            )));
        }
        if self.entities.is_empty() {
            return Err(PlannerError::Model("entity_count must be >= 1".to_string()));
        }
        if self.entities.len() > MAX_ENTITIES {
            return Err(PlannerError::Model(format!(
                "entity_count {} exceeds 2^20",
                self.entities.len()
            )));
        }
        if self.resources.len() > MAX_RESOURCES {
            return Err(PlannerError::Model(format!(
                "resource count {} exceeds 64",
                self.resources.len()
            )));
        }

        // ---- 资源结构 ----
        for res in &self.resources {
            if res.capacity.len() != self.value_count {
                return Err(PlannerError::Model(format!(
                    "resource '{}': capacity length {} != value_count {}",
                    res.name,
                    res.capacity.len(),
                    self.value_count
                )));
            }
            for (v, &c) in res.capacity.iter().enumerate() {
                if c < 0 {
                    return Err(PlannerError::Model(format!(
                        "resource '{}': negative capacity {} at value {v}",
                        res.name, c
                    )));
                }
                if c as i128 > MAX_COST_OR_CAP && c != CAP_SENTINEL {
                    return Err(PlannerError::Model(format!(
                        "resource '{}': capacity {} at value {v} exceeds 2^40",
                        res.name, c
                    )));
                }
            }
        }

        // ---- cost 矩阵（行长度、单值界、总和界） ----
        if self.cost.len() != self.entities.len() {
            return Err(PlannerError::Model(format!(
                "cost row count {} != entity count {}",
                self.cost.len(),
                self.entities.len()
            )));
        }
        let mut cost_total: i128 = 0;
        for (e, row) in self.cost.iter().enumerate() {
            if row.len() != self.value_count {
                return Err(PlannerError::Model(format!(
                    "entity {e}: cost row length {} != value_count {}",
                    row.len(),
                    self.value_count
                )));
            }
            for (v, &c) in row.iter().enumerate() {
                if c < 0 {
                    return Err(PlannerError::Model(format!(
                        "entity {e}: negative cost {c} at value {v}"
                    )));
                }
                if c as i128 > MAX_COST_OR_CAP {
                    return Err(PlannerError::Model(format!(
                        "entity {e}: cost {c} at value {v} exceeds 2^40"
                    )));
                }
                cost_total += c as i128;
            }
        }
        if cost_total > MAX_TOTAL {
            return Err(PlannerError::Model(format!(
                "total cost {cost_total} exceeds 2^40"
            )));
        }

        // ---- 实体结构（值域 / 需求 / 冲突邻接） ----
        for (e, ent) in self.entities.iter().enumerate() {
            if ent.domain.is_empty() {
                return Err(PlannerError::Model(format!("entity {e}: empty domain")));
            }
            for w in ent.domain.windows(2) {
                if w[0] >= w[1] {
                    return Err(PlannerError::Model(format!(
                        "entity {e}: domain not strictly ascending ({} after {})",
                        w[1], w[0]
                    )));
                }
            }
            for &v in &ent.domain {
                if v as usize >= self.value_count {
                    return Err(PlannerError::Model(format!(
                        "entity {e}: domain value {v} >= value_count {}",
                        self.value_count
                    )));
                }
            }
            if ent.demands.len() != self.resources.len() {
                return Err(PlannerError::Model(format!(
                    "entity {e}: demands length {} != resource count {}",
                    ent.demands.len(),
                    self.resources.len()
                )));
            }
            for (r, &d) in ent.demands.iter().enumerate() {
                if d < 0 {
                    return Err(PlannerError::Model(format!(
                        "entity {e}: negative demand {d} on resource {r}"
                    )));
                }
                if d as i128 > MAX_DEMAND {
                    return Err(PlannerError::Model(format!(
                        "entity {e}: demand {d} on resource {r} exceeds 2^31"
                    )));
                }
            }
            for &c in &ent.conflicts {
                if c == e {
                    return Err(PlannerError::Model(format!(
                        "entity {e}: self-loop conflict"
                    )));
                }
                if c >= self.entities.len() {
                    return Err(PlannerError::Model(format!(
                        "entity {e}: conflict id {c} out of range"
                    )));
                }
            }
            let mut sorted_conflicts = ent.conflicts.clone();
            sorted_conflicts.sort_unstable();
            if sorted_conflicts.windows(2).any(|w| w[0] == w[1]) {
                return Err(PlannerError::Model(format!(
                    "entity {e}: duplicate conflict entry"
                )));
            }
            for &c in &ent.conflicts {
                if !self.entities[c].conflicts.contains(&e) {
                    return Err(PlannerError::Model(format!(
                        "conflict {c} <-> {e} not symmetric"
                    )));
                }
            }
        }

        // ---- 每资源需求总和（i128） ----
        for (r, res) in self.resources.iter().enumerate() {
            let total: i128 = self.entities.iter().map(|ent| ent.demands[r] as i128).sum();
            if total > MAX_TOTAL {
                return Err(PlannerError::Model(format!(
                    "resource '{}': total demand {total} exceeds 2^40",
                    res.name
                )));
            }
        }
        Ok(())
    }
}

/// 工作解：模型 + 每实体当前取值（[`UNASSIGNED`] 表示未分配）。
///
/// crate 内唯一赋值入口是 [`Solution::set`]（`pub(crate)`），因此 lib 外部无法绕过
/// [`crate::moves::Move::apply`] 直接变更解——这是"影子状态与解永不失步"的结构保证。
#[derive(Debug, Clone)]
pub struct Solution {
    model: Arc<ProblemModel>,
    assignment: Box<[ValueId]>,
}

impl Solution {
    /// 空解：全部实体 [`UNASSIGNED`]，评分 [`HardSoftScore::ZERO`]。
    pub fn empty(model: Arc<ProblemModel>) -> Solution {
        Solution {
            assignment: vec![UNASSIGNED; model.entity_count()].into_boxed_slice(),
            model,
        }
    }

    /// 由赋值向量构造解；校验长度与值域（每项须为域内值或 [`UNASSIGNED`]），违例返回 `Err`。
    pub fn from_assignment(
        model: Arc<ProblemModel>,
        assignment: &[ValueId],
    ) -> Result<Solution, PlannerError> {
        if assignment.len() != model.entity_count() {
            return Err(PlannerError::Invalid(format!(
                "assignment length {} != entity count {}",
                assignment.len(),
                model.entity_count()
            )));
        }
        for (e, &v) in assignment.iter().enumerate() {
            if v != UNASSIGNED && !model.in_domain(e, v) {
                return Err(PlannerError::Invalid(format!(
                    "value {v} not in domain of entity {e}"
                )));
            }
        }
        Ok(Solution {
            model,
            assignment: assignment.to_vec().into_boxed_slice(),
        })
    }

    /// 所属模型。
    pub fn model(&self) -> &ProblemModel {
        &self.model
    }

    /// 实体 `e` 的当前取值；未分配返回 [`UNASSIGNED`]。
    pub fn value(&self, e: EntityId) -> ValueId {
        self.assignment[e]
    }

    /// 完整赋值向量（下标 = 实体 id）。
    pub fn assignment(&self) -> &[ValueId] {
        &self.assignment
    }

    /// 是否全部实体已分配（无 [`UNASSIGNED`]）。
    pub fn is_complete(&self) -> bool {
        self.assignment.iter().all(|&v| v != UNASSIGNED)
    }

    /// crate 内唯一赋值入口（lib 外不可调用）；仅由 [`crate::construct::construct`] 与
    /// [`crate::moves::Move::apply`] 调用，且必须与 [`crate::score::IncrementalScorer`]
    /// 的 apply_* 在同一函数内成对执行。
    pub(crate) fn set(&mut self, e: EntityId, v: ValueId) {
        self.assignment[e] = v;
    }
}

/// 统一错误类型：解析（带 1 起始行号）、模型校验、运行期非法参数、IO。
/// 用户输入路径不引发 panic——所有可预期错误都以本类型返回。
#[derive(Debug)]
pub enum PlannerError {
    /// 文本算例解析错误（语法 / 指令顺序 / id 越界 / 数值超界）；`line` 为 1 起始行号。
    Parse { line: usize, message: String },
    /// 模型结构 / 规模硬界校验失败（如值域为空、Σcost 超界）。
    Model(String),
    /// 运行期非法参数（如 local_search 配置越界：history_len = 0 等）。
    Invalid(String),
    /// 文件 / IO 错误。
    Io(std::io::Error),
}

impl std::fmt::Display for PlannerError {
    /// 每个 variant 的可读消息（英文；格式如 `parse error at line 3: ...`）。
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PlannerError::Parse { line, message } => {
                write!(f, "parse error at line {line}: {message}")
            }
            PlannerError::Model(message) => write!(f, "invalid model: {message}"),
            PlannerError::Invalid(message) => write!(f, "invalid argument: {message}"),
            PlannerError::Io(err) => write!(f, "io error: {err}"),
        }
    }
}

impl std::error::Error for PlannerError {}

#[cfg(test)]
mod tests {
    use super::*;

    /// 构造最小合法模型：1 资源（容量 [3,3]）、2 实体全值域、需求 0、无冲突、代价全 0。
    fn minimal_model() -> ProblemModel {
        ProblemModel {
            name: "minimal".to_string(),
            value_count: 2,
            resources: vec![Resource {
                name: "r0".to_string(),
                capacity: vec![3, 3],
            }],
            entities: vec![
                PlanningEntity {
                    domain: vec![0, 1],
                    demands: vec![0],
                    conflicts: vec![],
                },
                PlanningEntity {
                    domain: vec![0, 1],
                    demands: vec![0],
                    conflicts: vec![],
                },
            ],
            cost: vec![vec![0, 0], vec![0, 0]],
        }
    }

    #[test]
    fn score_lexicographic_order() {
        // (hard, soft) 字典序：hard 优先于 soft
        assert!(HardSoftScore::new(0, 5) < HardSoftScore::new(1, 0));
        assert!(HardSoftScore::new(0, i64::MAX) < HardSoftScore::new(1, i64::MIN));
        assert!(HardSoftScore::new(2, 10) < HardSoftScore::new(2, 11));
        assert_eq!(HardSoftScore::new(3, 4), HardSoftScore::new(3, 4));
        assert!(HardSoftScore::ZERO < HardSoftScore::new(0, 1));
        assert!(HardSoftScore::new(-1, 0) < HardSoftScore::ZERO, "负 hard 亦按字典序（构造器不生成，但序必须一致）");
        assert_eq!(HardSoftScore::default(), HardSoftScore::ZERO);
    }

    #[test]
    fn score_feasibility() {
        assert!(HardSoftScore::ZERO.is_feasible());
        assert!(HardSoftScore::new(0, 100).is_feasible());
        assert!(!HardSoftScore::new(1, 0).is_feasible());
        assert!(!HardSoftScore::new(-1, 0).is_feasible());
    }

    #[test]
    fn accessors_and_domain_queries() {
        let m = minimal_model();
        assert_eq!(m.entity_count(), 2);
        assert_eq!(m.value_count(), 2);
        assert_eq!(m.domain(1), &[0, 1]);
        assert!(m.in_domain(0, 0));
        assert!(m.in_domain(1, 1));
        assert!(!m.in_domain(0, 2));
        assert!(!m.in_domain(0, UNASSIGNED), "UNASSIGNED is never in a domain");
        assert!(m.validate().is_ok(), "minimal model must validate");
    }

    #[test]
    fn validate_rejects_scale_bound_violations() {
        // m 超 2^31（行长度故意不匹配：validate 必须先拦 value_count 本身）
        let mut m = minimal_model();
        m.value_count = (1usize << 31) + 1;
        match m.validate() {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("value_count"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }

        // n 超 2^20
        let mut m = minimal_model();
        m.entities = vec![
            PlanningEntity { domain: vec![0, 1], demands: vec![0], conflicts: vec![] };
            (1usize << 20) + 1
        ];
        m.cost = vec![vec![0, 0]; m.entities.len()];
        match m.validate() {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("2^20"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }

        // 资源数 65
        let mut m = minimal_model();
        m.resources = (0..65)
            .map(|i| Resource { name: format!("r{i}"), capacity: vec![3, 3] })
            .collect();
        for ent in &mut m.entities {
            ent.demands = vec![0; 65];
        }
        match m.validate() {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("exceeds 64"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }

        // 单值 cost 超 2^40
        let mut m = minimal_model();
        m.cost[0][0] = (1i64 << 40) + 1;
        match m.validate() {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("exceeds 2^40"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }

        // 单值 cap 超 2^40（且非 CAP 哨兵）
        let mut m = minimal_model();
        m.resources[0].capacity[0] = (1i64 << 40) + 7;
        match m.validate() {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("exceeds 2^40"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }
        // CAP 哨兵（i64::MAX/2）是唯一合法的超 2^40 容量
        m.resources[0].capacity[0] = i64::MAX / 2;
        assert!(m.validate().is_ok());

        // 单值 demand 超 2^31
        let mut m = minimal_model();
        m.entities[0].demands[0] = (1i64 << 31) + 1;
        match m.validate() {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("exceeds 2^31"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }

        // Σcost 超 2^40（两处恰为 2^40，单值合法、总和非法）
        let mut m = minimal_model();
        m.cost[0][0] = 1i64 << 40;
        m.cost[1][1] = 1i64 << 40;
        match m.validate() {
            Err(PlannerError::Model(msg)) => {
                assert!(msg.contains("total cost") && msg.contains("exceeds"), "got: {msg}")
            }
            other => panic!("expected Model error, got {other:?}"),
        }

        // 每资源 Σdemand 超 2^40（单值均 ≤ 2^31：513 实体 × 2^31 > 2^40）
        let mut m = minimal_model();
        m.entities = (0..513)
            .map(|_| PlanningEntity { domain: vec![0, 1], demands: vec![1i64 << 31], conflicts: vec![] })
            .collect();
        m.cost = vec![vec![0, 0]; m.entities.len()];
        match m.validate() {
            Err(PlannerError::Model(msg)) => {
                assert!(msg.contains("total demand") && msg.contains("exceeds"), "got: {msg}")
            }
            other => panic!("expected Model error, got {other:?}"),
        }
    }

    #[test]
    fn validate_rejects_structural_violations() {
        // 空域
        let mut m = minimal_model();
        m.entities[0].domain = vec![];
        match m.validate() {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("empty domain"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }

        // 域值越界 / 非升序
        let mut m = minimal_model();
        m.entities[0].domain = vec![0, 2];
        match m.validate() {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("value_count"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }
        let mut m = minimal_model();
        m.entities[0].domain = vec![1, 1];
        match m.validate() {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("ascending"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }

        // demands / capacity / cost 行长不匹配
        let mut m = minimal_model();
        m.entities[0].demands = vec![];
        assert!(matches!(m.validate(), Err(PlannerError::Model(_))));
        let mut m = minimal_model();
        m.resources[0].capacity = vec![3];
        assert!(matches!(m.validate(), Err(PlannerError::Model(_))));
        let mut m = minimal_model();
        m.cost[1] = vec![0, 0, 0];
        assert!(matches!(m.validate(), Err(PlannerError::Model(_))));

        // 自环冲突
        let mut m = minimal_model();
        m.entities[0].conflicts = vec![0];
        match m.validate() {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("self-loop"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }

        // 非对称冲突
        let mut m = minimal_model();
        m.entities[0].conflicts = vec![1];
        match m.validate() {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("symmetric"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }

        // 重复冲突条目
        let mut m = minimal_model();
        m.entities[0].conflicts = vec![1, 1];
        m.entities[1].conflicts = vec![0];
        match m.validate() {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("duplicate"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }

        // 冲突 id 越界
        let mut m = minimal_model();
        m.entities[0].conflicts = vec![9];
        match m.validate() {
            Err(PlannerError::Model(msg)) => assert!(msg.contains("out of range"), "got: {msg}"),
            other => panic!("expected Model error, got {other:?}"),
        }

        // 负值
        let mut m = minimal_model();
        m.cost[0][0] = -1;
        assert!(matches!(m.validate(), Err(PlannerError::Model(_))));
        let mut m = minimal_model();
        m.entities[0].demands[0] = -1;
        assert!(matches!(m.validate(), Err(PlannerError::Model(_))));
        let mut m = minimal_model();
        m.resources[0].capacity[0] = -5;
        assert!(matches!(m.validate(), Err(PlannerError::Model(_))));
    }

    #[test]
    fn solution_empty_and_from_assignment() {
        let model = std::sync::Arc::new(minimal_model());
        let sol = Solution::empty(model.clone());
        assert_eq!(sol.value(0), UNASSIGNED);
        assert_eq!(sol.value(1), UNASSIGNED);
        assert!(!sol.is_complete());
        assert_eq!(sol.assignment(), &[UNASSIGNED, UNASSIGNED]);
        // 模型访问与归属
        assert_eq!(sol.model().name, "minimal");

        let sol = Solution::from_assignment(model.clone(), &[0, 1]).expect("valid assignment");
        assert_eq!(sol.value(0), 0);
        assert_eq!(sol.value(1), 1);
        assert!(sol.is_complete());

        // 长度不符
        match Solution::from_assignment(model.clone(), &[0]) {
            Err(PlannerError::Invalid(_)) => {}
            other => panic!("expected Invalid, got {other:?}"),
        }
        // 域外值
        match Solution::from_assignment(model.clone(), &[0, 2]) {
            Err(PlannerError::Invalid(msg)) => assert!(msg.contains("not in domain"), "got: {msg}"),
            other => panic!("expected Invalid, got {other:?}"),
        }
        // UNASSIGNED 允许（部分解）
        let partial = Solution::from_assignment(model.clone(), &[1, UNASSIGNED]).expect("partial");
        assert!(!partial.is_complete());
        assert_eq!(partial.value(1), UNASSIGNED);
    }

    #[test]
    fn error_display_messages() {
        let e = PlannerError::Parse { line: 3, message: "unknown directive 'X'".to_string() };
        assert_eq!(e.to_string(), "parse error at line 3: unknown directive 'X'");
        let e = PlannerError::Model("empty domain".to_string());
        assert_eq!(e.to_string(), "invalid model: empty domain");
        let e = PlannerError::Invalid("bad arg".to_string());
        assert_eq!(e.to_string(), "invalid argument: bad arg");
        let e = PlannerError::Io(std::io::Error::new(std::io::ErrorKind::NotFound, "no file"));
        assert!(e.to_string().starts_with("io error: "));
        // 实现 std::error::Error（source / display 可用）
        let _: &dyn std::error::Error = &e;
    }
}
