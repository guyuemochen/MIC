//! 业务表达式按定义顺序组成无环图。全量计算从工作解独立重建各项，
//! 搜索时由变动的叶节点向使用它的表达式和具名规则传播。
//! 未赋值的标量为未知；算术和比较传播未知，聚合只统计已知项。
//! 构造期间尚未确定的罚项和目标暂记零。条件只使用当前命中的分支，未命中分支
//! 上的数据越界或算术溢出随节点缓存，只有影响到具名规则或目标时才返回错误。

use std::collections::{BTreeMap, BTreeSet, VecDeque};

mod collisions;

use super::sequence::{self, Edit, Schedule};
use super::{Expr, ListSpec, Model, State};
use collisions::PairCounts;

/// 各规则的罚分按声明顺序返回，软分还包括目标表达式。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Breakdown {
    pub hard: i64,
    pub soft: i64,
    pub penalties: Vec<i64>,
}

fn checked(value: i128) -> Result<i64, String> {
    i64::try_from(value).map_err(|_| "评分数值超出 i64 范围".to_owned())
}

fn inputs(expr: &Expr) -> Vec<usize> {
    use Expr::*;
    match expr {
        Add(a, b)
        | Sub(a, b)
        | Mul(a, b)
        | Eq(a, b)
        | Ne(a, b)
        | Lt(a, b)
        | Le(a, b)
        | Gt(a, b)
        | Ge(a, b) => vec![*a, *b],
        If(c, a, b) => vec![*c, *a, *b],
        Sum(xs) | Max(xs) => xs.clone(),
        Lookup(_, idx) => vec![*idx],
        _ => vec![],
    }
}

fn unique<T: Ord>(items: &[T]) -> bool {
    items.iter().collect::<BTreeSet<_>>().len() == items.len()
}

/// 检查建模接口产生的索引、静态值域、维度和工艺前置关系。
/// 所有列表组的所有者总数最多为 2²⁰，避免一个计数在空工作解中触发无界分配。
pub fn validate(model: &Model) -> Result<(), String> {
    for (i, domain) in model.domains.iter().enumerate() {
        if domain.is_empty() || !unique(domain) {
            return Err(format!("变量 {i} 的值域必须非空且不能重复"));
        }
    }
    let mut remaining_owners = 1usize << 20;
    for (g, spec) in model.lists.iter().enumerate() {
        remaining_owners = remaining_owners
            .checked_sub(spec.owners)
            .ok_or("列表所有者总数超过 2²⁰")?;
        let n = spec.eligible.len();
        if spec.owners == 0
            || (!spec.durations.is_empty() && spec.durations.len() != n)
            || spec.durations.iter().any(|&d| d < 0)
        {
            return Err(format!("列表组 {g} 的所有者数量或工时表无效"));
        }
        for owners in &spec.eligible {
            if owners.is_empty() || owners.iter().any(|&o| o >= spec.owners) || !unique(owners) {
                return Err(format!("列表组 {g} 含空、重复或越界的候选所有者"));
            }
        }
        if spec.precedences.iter().any(|&(a, b)| a >= n || b >= n) || !unique(&spec.precedences) {
            return Err(format!("列表组 {g} 的工艺前置关系无效"));
        }
        let mut next = vec![Vec::new(); n];
        let mut deg = vec![0; n];
        for &(a, b) in &spec.precedences {
            next[a].push(b);
            deg[b] += 1;
        }
        let mut queue: VecDeque<_> = (0..n).filter(|&i| deg[i] == 0).collect();
        let mut count = 0;
        while let Some(i) = queue.pop_front() {
            count += 1;
            for &j in &next[i] {
                deg[j] -= 1;
                if deg[j] == 0 {
                    queue.push_back(j);
                }
            }
        }
        if count != n {
            return Err(format!("列表组 {g} 的工艺前置关系存在环"));
        }
    }
    for (i, node) in model.nodes.iter().enumerate() {
        if inputs(node).iter().any(|&j| j >= i) {
            return Err(format!("表达式 {i} 引用了未定义的表达式"));
        }
        use Expr::*;
        match node {
            Var(v) if *v >= model.domains.len() => {
                return Err(format!("表达式 {i} 引用了未知变量"))
            }
            Load {
                group,
                owner,
                demands,
            } => {
                let s = model
                    .lists
                    .get(*group)
                    .ok_or_else(|| format!("表达式 {i} 引用了未知列表组"))?;
                if *owner >= s.owners || demands.len() != s.eligible.len() {
                    return Err(format!("表达式 {i} 的负载表或所有者无效"));
                }
            }
            Distance {
                group,
                owner,
                matrix,
            } => {
                let s = model
                    .lists
                    .get(*group)
                    .ok_or_else(|| format!("表达式 {i} 引用了未知列表组"))?;
                let size = s
                    .eligible
                    .len()
                    .checked_add(1)
                    .and_then(|n| n.checked_mul(n));
                if *owner >= s.owners || size != Some(matrix.len()) {
                    return Err(format!("表达式 {i} 的距离表或所有者无效"));
                }
            }
            Adjacent {
                group,
                owner,
                before,
                after,
            } => {
                let s = model
                    .lists
                    .get(*group)
                    .ok_or_else(|| format!("表达式 {i} 引用了未知列表组"))?;
                if *owner >= s.owners || *before >= s.eligible.len() || *after >= s.eligible.len() {
                    return Err(format!("表达式 {i} 的相邻关系索引无效"));
                }
            }
            Finish { group, item } => {
                if model
                    .lists
                    .get(*group)
                    .map_or(true, |s| *item >= s.eligible.len())
                {
                    return Err(format!("表达式 {i} 的工序索引无效"));
                }
            }
            _ => {}
        }
    }
    let mut names = BTreeSet::new();
    for r in &model.rules {
        if r.name.is_empty()
            || !names.insert(&r.name)
            || r.weight < 0
            || (r.hard && r.weight == 0)
            || r.expr >= model.nodes.len()
        {
            return Err(format!("约束名称、权重或表达式无效：{}", r.name));
        }
    }
    if model.objective.is_some_and(|e| e >= model.nodes.len()) {
        return Err("目标引用了未知表达式".into());
    }
    Ok(())
}

fn check_state(model: &Model, state: &State) -> Result<(), String> {
    if state.values.len() != model.domains.len() || state.lists.len() != model.lists.len() {
        return Err("工作解维度不符".into());
    }
    for (x, domain) in state.values.iter().zip(&model.domains) {
        if x.is_some_and(|v| !domain.contains(&v)) {
            return Err("工作解变量取值不在值域内".into());
        }
    }
    for (routes, spec) in state.lists.iter().zip(&model.lists) {
        if routes.len() != spec.owners {
            return Err("工作解列表维度不符".into());
        }
        let mut seen = vec![false; spec.eligible.len()];
        for (o, route) in routes.iter().enumerate() {
            for &i in route {
                if i >= seen.len() || seen[i] || !spec.eligible[i].contains(&o) {
                    return Err("工作解列表包含重复、越界或不合格的元素".into());
                }
                seen[i] = true;
            }
        }
    }
    Ok(())
}

/// 将静态前置和同机器相邻边都计入，按拓扑顺序计算最长完工时间。
/// 这里直接读取工作解，独立于搜索时维护的前驱和后继缓存。
fn finish_times(spec: &ListSpec, routes: &[Vec<usize>]) -> Result<Vec<Option<i128>>, String> {
    let n = spec.eligible.len();
    let mut present = vec![false; n];
    for &i in routes.iter().flatten() {
        present[i] = true;
    }
    let mut next = vec![Vec::new(); n];
    let mut degree = vec![0; n];
    let mut start = vec![0i128; n];
    let mut finish = vec![None; n];
    let machine_edges = routes
        .iter()
        .flat_map(|r| r.windows(2).map(|p| (p[0], p[1])));
    for (a, b) in spec.precedences.iter().copied().chain(machine_edges) {
        if present[a] && present[b] {
            next[a].push(b);
            degree[b] += 1;
        }
    }
    let mut queue: VecDeque<_> = (0..n).filter(|&i| present[i] && degree[i] == 0).collect();
    while let Some(i) = queue.pop_front() {
        let t = start[i] + i128::from(spec.durations.get(i).copied().unwrap_or(0));
        finish[i] = Some(t);
        for &j in &next[i] {
            start[j] = start[j].max(t);
            degree[j] -= 1;
            if degree[j] == 0 {
                queue.push_back(j);
            }
        }
    }
    if present.iter().zip(&finish).any(|(p, f)| *p && f.is_none()) {
        return Err(sequence::CYCLE.into());
    }
    Ok(finish)
}

type Value = Option<i64>;
type Faults = BTreeMap<usize, String>;

/// 只编译可从规则或目标到达的表达式。条件的两支都可以在未来的移动中被选中。
fn used_nodes(model: &Model) -> Vec<bool> {
    let mut used = vec![false; model.nodes.len()];
    for i in model.rules.iter().map(|r| r.expr).chain(model.objective) {
        used[i] = true;
    }
    for i in (0..used.len()).rev() {
        if used[i] {
            for j in inputs(&model.nodes[i]) {
                used[j] = true;
            }
        }
    }
    used
}

/// 求值错误属于产生它的节点。条件仅传播选择分支上的错误，未赋值的条件暂不选择。
fn inherited_fault(node: &Expr, values: &[Value], faults: &Faults) -> Option<String> {
    if faults.is_empty() {
        return None;
    }
    match node {
        Expr::If(c, a, b) => faults
            .get(c)
            .or_else(|| values[*c].and_then(|v| faults.get(if v != 0 { a } else { b })))
            .cloned(),
        Expr::Sum(xs) | Expr::Max(xs) => xs.iter().find_map(|i| faults.get(i)).cloned(),
        _ => inputs(node)
            .into_iter()
            .find_map(|i| faults.get(&i))
            .cloned(),
    }
}

/// 共用纯算术定义。列表叶节点由调用者提供，全量和增量计算有各自的数据来源。
fn evaluate(node: &Expr, values: &[Value], state: &State) -> Result<Value, String> {
    use Expr::*;
    match node {
        Const(v) => Ok(Some(*v)),
        Var(v) => Ok(state.values[*v]),
        If(c, a, b) => Ok(values[*c].and_then(|c| if c != 0 { values[*a] } else { values[*b] })),
        Sum(xs) => Ok(Some(checked(
            xs.iter().filter_map(|&i| values[i]).map(i128::from).sum(),
        )?)),
        Max(xs) => Ok(Some(
            xs.iter().filter_map(|&i| values[i]).max().unwrap_or(0),
        )),
        Lookup(table, i) => match values[*i] {
            None => Ok(None),
            Some(k) => usize::try_from(k)
                .ok()
                .and_then(|j| table.get(j))
                .copied()
                .map(Some)
                .ok_or_else(|| format!("查表下标 {k} 超出长度 {}", table.len())),
        },
        Add(a, b)
        | Sub(a, b)
        | Mul(a, b)
        | Eq(a, b)
        | Ne(a, b)
        | Lt(a, b)
        | Le(a, b)
        | Gt(a, b)
        | Ge(a, b) => {
            let (Some(a), Some(b)) = (values[*a], values[*b]) else {
                return Ok(None);
            };
            let v = match node {
                Add(..) => checked(i128::from(a) + i128::from(b))?,
                Sub(..) => checked(i128::from(a) - i128::from(b))?,
                Mul(..) => checked(i128::from(a) * i128::from(b))?,
                Eq(..) => i64::from(a == b),
                Ne(..) => i64::from(a != b),
                Lt(..) => i64::from(a < b),
                Le(..) => i64::from(a <= b),
                Gt(..) => i64::from(a > b),
                Ge(..) => i64::from(a >= b),
                _ => unreachable!(),
            };
            Ok(Some(v))
        }
        Load { .. } | Distance { .. } | Adjacent { .. } | Finish { .. } => {
            unreachable!("列表表达式由调用方直接计算")
        }
    }
}

fn penalty(value: Value, weight: i64) -> Result<i64, String> {
    let v = value.unwrap_or(0);
    if v < 0 {
        return Err("罚量必须非负".into());
    }
    checked(i128::from(v) * i128::from(weight))
}

fn breakdown(model: &Model, values: &[Value], faults: &Faults) -> Result<Breakdown, String> {
    for i in model.rules.iter().map(|r| r.expr).chain(model.objective) {
        if let Some(e) = faults.get(&i) {
            return Err(e.clone());
        }
    }
    let mut hard = 0i128;
    let mut soft = i128::from(model.objective.and_then(|i| values[i]).unwrap_or(0));
    let mut penalties = Vec::with_capacity(model.rules.len());
    for r in &model.rules {
        let p = penalty(values[r.expr], r.weight)?;
        if r.hard {
            hard += i128::from(p);
        } else {
            soft += i128::from(p);
        }
        penalties.push(p);
    }
    Ok(Breakdown {
        hard: checked(hard)?,
        soft: checked(soft)?,
        penalties,
    })
}

/// 全量路线分数通过工作解中的完整序列独立计算。在某个路线叶节点原值超过
/// i64 的罕见情况下，下一次变化也用它恢复数值，不会拿未知值当作零。
fn route_value(node: &Expr, model: &Model, state: &State) -> Result<Value, String> {
    let total = match node {
        Expr::Load {
            group,
            owner,
            demands,
        } => state.lists[*group][*owner]
            .iter()
            .map(|&j| i128::from(demands[j]))
            .sum(),
        Expr::Distance {
            group,
            owner,
            matrix,
        } => {
            let n = model.lists[*group].eligible.len();
            let mut prev = n;
            let mut total = 0i128;
            let route = &state.lists[*group][*owner];
            for &j in route {
                total += i128::from(matrix[prev * (n + 1) + j]);
                prev = j;
            }
            if !route.is_empty() {
                total += i128::from(matrix[prev * (n + 1) + n]);
            }
            total
        }
        Expr::Adjacent {
            group,
            owner,
            before,
            after,
        } => i128::from(
            state.lists[*group][*owner]
                .windows(2)
                .any(|w| w == [*before, *after]),
        ),
        _ => unreachable!("仅用于路线分数"),
    };
    checked(total).map(Some)
}

fn initial_values(
    model: &Model,
    state: &State,
    used: &[bool],
) -> Result<(Vec<Value>, Faults), String> {
    // 配送路线不需要工序图；有静态前置边时，即使不读取时间也要检查成环。
    let mut timed: Vec<_> = model
        .lists
        .iter()
        .map(|spec| !spec.precedences.is_empty())
        .collect();
    for (i, node) in model.nodes.iter().enumerate() {
        if let Expr::Finish { group, .. } = node {
            timed[*group] |= used[i];
        }
    }
    let finishes = model
        .lists
        .iter()
        .zip(&state.lists)
        .zip(timed)
        .map(|((spec, routes), needed)| {
            if needed {
                finish_times(spec, routes)
            } else {
                Ok(vec![])
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut values = vec![None; model.nodes.len()];
    let mut faults = Faults::new();
    for (i, node) in model.nodes.iter().enumerate() {
        if !used[i] {
            continue;
        }
        let result = if let Some(error) = inherited_fault(node, &values, &faults) {
            Err(error)
        } else {
            match node {
                Expr::Load { .. } | Expr::Distance { .. } | Expr::Adjacent { .. } => {
                    route_value(node, model, state)
                }
                Expr::Finish { group, item } => finishes[*group][*item].map(checked).transpose(),
                _ => evaluate(node, &values, state),
            }
        };
        match result {
            Ok(v) => values[i] = v,
            Err(e) => {
                faults.insert(i, e);
            }
        }
    }
    Ok((values, faults))
}

pub fn full(model: &Model, state: &State) -> Result<Breakdown, String> {
    validate(model)?;
    check_state(model, state)?;
    let (values, faults) = initial_values(model, state, &used_nodes(model))?;
    breakdown(model, &values, &faults)
}

/// 对求和保留精确的小计，对最大值记录频次。一个变量多次出现在同一聚合中时，
/// 反向索引也保留相同次数，更新重复输入不会漏算。
enum Aggregate {
    Sum(i128),
    Max(BTreeMap<i64, usize>),
    Pairs(Box<PairCounts>),
}

impl Aggregate {
    fn new(expr: &Expr, values: &[Value]) -> Option<Self> {
        match expr {
            Expr::Sum(xs) => Some(Self::Sum(
                xs.iter().filter_map(|&i| values[i]).map(i128::from).sum(),
            )),
            Expr::Max(xs) => {
                let mut counts = BTreeMap::new();
                for v in xs.iter().filter_map(|&i| values[i]) {
                    *counts.entry(v).or_default() += 1;
                }
                Some(Self::Max(counts))
            }
            _ => None,
        }
    }

    fn replace(&mut self, group: usize, old: Value, new: Value) {
        match self {
            Self::Pairs(pairs) => pairs.replace(group, old, new),
            Self::Sum(total) => {
                *total += i128::from(new.unwrap_or(0)) - i128::from(old.unwrap_or(0))
            }
            Self::Max(counts) => {
                if let Some(v) = old {
                    let count = counts.get_mut(&v).expect("旧聚合值应有正频次");
                    *count -= 1;
                    if *count == 0 {
                        counts.remove(&v);
                    }
                }
                if let Some(v) = new {
                    *counts.entry(v).or_default() += 1;
                }
            }
        }
    }

    fn value(&self) -> Result<Value, String> {
        match self {
            Self::Pairs(pairs) => checked(pairs.total()).map(Some),
            Self::Sum(total) => Ok(Some(checked(*total)?)),
            Self::Max(counts) => Ok(Some(counts.last_key_value().map_or(0, |(&v, _)| v))),
        }
    }
}

#[derive(Default)]
struct ListInputs {
    owners: BTreeMap<usize, Vec<usize>>,
    finishes: BTreeMap<usize, Vec<usize>>,
}

/// 每次工作解变更只通知已改的变量和列表编辑。表达式对此模型只编译一次。
pub(super) struct Scorer {
    values: Vec<Value>,
    faults: Faults,
    rules: Vec<Vec<usize>>,
    parents: Vec<Vec<(usize, usize)>>,
    variables: Vec<Vec<usize>>,
    lists: Vec<ListInputs>,
    schedules: Vec<Option<Schedule>>,
    aggregates: Vec<Option<Aggregate>>,
    scores: Breakdown,
}

impl Scorer {
    pub fn new(model: &Model, state: &State) -> Result<Self, String> {
        validate(model)?;
        check_state(model, state)?;
        let (values, mut faults) = initial_values(model, state, &used_nodes(model))?;
        let scores = breakdown(model, &values, &faults)?;
        let n = model.nodes.len();
        let mut parents = vec![vec![]; n];
        let mut rules = vec![vec![]; n];
        let mut used = vec![false; n];
        for i in model.rules.iter().map(|r| r.expr).chain(model.objective) {
            used[i] = true;
        }
        let mut aggregates = (0..n).map(|_| None).collect::<Vec<_>>();
        // 从规则向下编译。被冲突桶替代的比较若还有其它用处，会经普通依赖保留。
        for i in (0..n).rev() {
            if !used[i] {
                continue;
            }
            let expr = &model.nodes[i];
            let mut subscribe = |child: usize, group| {
                used[child] = true;
                parents[child].push((i, group));
            };
            if let Some(pairs) = PairCounts::compile(model, expr, &values, &mut subscribe) {
                aggregates[i] = Some(Aggregate::Pairs(Box::new(pairs)));
            } else {
                for child in inputs(expr) {
                    subscribe(child, 0);
                }
                aggregates[i] = Aggregate::new(expr, &values);
            }
        }
        faults.retain(|&i, _| used[i]);
        let mut variables = vec![vec![]; model.domains.len()];
        let mut lists: Vec<_> = model.lists.iter().map(|_| ListInputs::default()).collect();
        for (i, node) in model.nodes.iter().enumerate() {
            if !used[i] {
                continue;
            }
            match node {
                Expr::Var(v) => variables[*v].push(i),
                Expr::Load { group, owner, .. }
                | Expr::Distance { group, owner, .. }
                | Expr::Adjacent { group, owner, .. } => {
                    lists[*group].owners.entry(*owner).or_default().push(i);
                }
                Expr::Finish { group, item } => {
                    lists[*group].finishes.entry(*item).or_default().push(i)
                }
                _ => {}
            }
        }
        for (i, r) in model.rules.iter().enumerate() {
            rules[r.expr].push(i);
        }
        let mut schedules = Vec::with_capacity(model.lists.len());
        for (g, s) in model.lists.iter().enumerate() {
            schedules.push(
                if s.precedences.is_empty() && lists[g].finishes.is_empty() {
                    None
                } else {
                    Some(Schedule::new(s, &state.lists[g])?)
                },
            );
        }
        Ok(Self {
            scores,
            values,
            faults,
            rules,
            parents,
            variables,
            lists,
            schedules,
            aggregates,
        })
    }

    pub fn score(&self) -> (i64, i64) {
        (self.scores.hard, self.scores.soft)
    }

    /// 工作解须已应用这次移动；一次候选可改若干变量和同一列表组中的多个位置。
    /// 在写入评分缓存之前先检查机器顺序是否成环，拒绝这一候选时评分原样保留。
    pub fn update(
        &mut self,
        model: &Model,
        state: &State,
        vars: &[usize],
        edits: &[Edit],
    ) -> Result<(), String> {
        let mut dirty = BTreeSet::<usize>::new();
        let mut leaves = BTreeMap::<usize, Result<Value, String>>::new();
        for &v in vars {
            let ids = self
                .variables
                .get(v)
                .ok_or_else(|| "更新引用了未知变量".to_owned())?;
            dirty.extend(ids);
        }
        if let Some(first) = edits.first() {
            let g = first.group;
            if g >= self.lists.len() || edits.iter().any(|e| e.group != g) {
                return Err("一次列表更新必须属于同一已知组".into());
            }
            // 跨列表、同列表组合移动可能多次触及同一个叶节点，先合并所有增量，
            // 不让某个中途状态的上溢掩盖可表示的最终结果。
            let mut deltas = BTreeMap::<usize, i128>::new();
            for e in edits {
                if let Some(ids) = self.lists[g].owners.get(&e.owner) {
                    for &i in ids {
                        let delta = match &model.nodes[i] {
                            Expr::Load { demands, .. } => sequence::load_delta(e, demands)?,
                            Expr::Distance { matrix, .. } => {
                                sequence::distance_delta(e, matrix, model.lists[g].eligible.len())?
                            }
                            Expr::Adjacent { before, after, .. } => {
                                sequence::adjacent_delta(e, *before, *after)
                            }
                            _ => unreachable!(),
                        };
                        *deltas.entry(i).or_default() += delta;
                    }
                }
            }
            for (i, d) in deltas {
                if d != 0 {
                    let value = match self.values[i] {
                        Some(old) => checked(i128::from(old) + d).map(Some),
                        None => route_value(&model.nodes[i], model, state),
                    };
                    leaves.insert(i, value);
                    dirty.insert(i);
                }
            }
            if let Some(sched) = &mut self.schedules[g] {
                for item in sched.update(&model.lists[g], edits)? {
                    if let Some(ids) = self.lists[g].finishes.get(&item) {
                        for &i in ids {
                            dirty.insert(i);
                            leaves.insert(i, sched.finish(item).map(checked).transpose());
                        }
                    }
                }
            }
        }
        let mut changes = Vec::new();
        let mut penalties = Vec::new();
        let result = (|| -> Result<(i64, i64), String> {
            let mut hard = i128::from(self.scores.hard);
            let mut soft = i128::from(self.scores.soft);
            // 父表达式总在子表达式之后定义；按编号出队即可在同一步合并多个子项的变化。
            while let Some(i) = dirty.pop_first() {
                let old = self.values[i];
                let fault = if self.faults.is_empty() {
                    None
                } else if let Some(Aggregate::Pairs(pairs)) = &self.aggregates[i] {
                    pairs.fault(&self.faults)
                } else {
                    inherited_fault(&model.nodes[i], &self.values, &self.faults)
                };
                let result = if let Some(e) = fault {
                    Err(e)
                } else if let Some(v) = leaves.remove(&i) {
                    v
                } else if let Some(a) = &self.aggregates[i] {
                    a.value()
                } else {
                    evaluate(&model.nodes[i], &self.values, state)
                };
                let (new, error) = match result {
                    Ok(v) => (v, None),
                    Err(e) => (None, Some(e)),
                };
                if old == new && self.faults.get(&i) == error.as_ref() {
                    continue;
                }
                if !self.rules[i].is_empty() || model.objective == Some(i) {
                    if let Some(e) = &error {
                        return Err(e.clone());
                    }
                }
                let old_error = if let Some(e) = error {
                    self.faults.insert(i, e)
                } else {
                    self.faults.remove(&i)
                };
                changes.push((i, old, old_error));
                self.values[i] = new;
                for &(p, group) in &self.parents[i] {
                    if let Some(a) = &mut self.aggregates[p] {
                        a.replace(group, old, new);
                    }
                    dirty.insert(p);
                }
                for &j in &self.rules[i] {
                    let rule = &model.rules[j];
                    let val = penalty(new, rule.weight)?;
                    let delta = i128::from(val) - i128::from(self.scores.penalties[j]);
                    if rule.hard {
                        hard += delta;
                    } else {
                        soft += delta;
                    }
                    penalties.push((j, val));
                }
                if model.objective == Some(i) {
                    soft += i128::from(new.unwrap_or(0)) - i128::from(old.unwrap_or(0));
                }
            }
            Ok((checked(hard)?, checked(soft)?))
        })();
        match result {
            Ok((hard, soft)) => {
                self.scores.hard = hard;
                self.scores.soft = soft;
                for (i, val) in penalties {
                    self.scores.penalties[i] = val;
                }
            }
            Err(error) => {
                // 随候选记录触及的节点；错误只撤销这一小部分，原状态仍可继续评分。
                for (i, old, error) in changes.into_iter().rev() {
                    if let Some(e) = error {
                        self.faults.insert(i, e);
                    } else {
                        self.faults.remove(&i);
                    }
                    let new = std::mem::replace(&mut self.values[i], old);
                    for &(p, group) in &self.parents[i] {
                        if let Some(a) = &mut self.aggregates[p] {
                            a.replace(group, new, old);
                        }
                    }
                }
                if let Some(g) = edits.first().map(|e| e.group) {
                    if let Some(schedule) = &mut self.schedules[g] {
                        let inverse = edits.iter().rev().map(Edit::inverse).collect::<Vec<_>>();
                        schedule.update(&model.lists[g], &inverse)?;
                    }
                }
                return Err(error);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
