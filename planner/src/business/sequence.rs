//! 有序列表的局部编辑、配送路段差分与工序完工时间传播。
//! 每个编辑保存操作当时的两侧邻居，同一批编辑按生成顺序解释。

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use super::{ListSpec, State};

pub const CYCLE: &str = "工序顺序成环";

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Edit {
    pub group: usize,
    pub owner: usize,
    pub pos: usize,
    pub removed: Vec<usize>,
    pub added: Vec<usize>,
    pub left: Option<usize>,
    pub right: Option<usize>,
}

impl Edit {
    pub fn inverse(&self) -> Self {
        Self {
            group: self.group,
            owner: self.owner,
            pos: self.pos,
            removed: self.added.clone(),
            added: self.removed.clone(),
            left: self.left,
            right: self.right,
        }
    }
}

/// 替换连续区段。位置和长度检查在改动列表之前完成。
pub fn splice(
    state: &mut State,
    group: usize,
    owner: usize,
    pos: usize,
    count: usize,
    added: &[usize],
) -> Result<Edit, String> {
    let route = state
        .lists
        .get_mut(group)
        .and_then(|lists| lists.get_mut(owner))
        .ok_or("列表编号越界")?;
    if pos > route.len() || count > route.len() - pos {
        return Err("列表区段越界".into());
    }
    let left = pos.checked_sub(1).map(|i| route[i]);
    let right = route.get(pos + count).copied();
    let removed = route
        .splice(pos..pos + count, added.iter().copied())
        .collect();
    Ok(Edit {
        group,
        owner,
        pos,
        removed,
        added: added.to_vec(),
        left,
        right,
    })
}

/// 撤销由本模块刚刚应用的一整批编辑，按相反顺序恢复，故原始位置仍有效。
pub fn undo(state: &mut State, edits: &[Edit]) {
    for edit in edits.iter().rev() {
        let route = &mut state.lists[edit.group][edit.owner];
        route.splice(
            edit.pos..edit.pos + edit.added.len(),
            edit.removed.iter().copied(),
        );
    }
}

/// 将一段任务保持原顺序地移到目标位置，目标位置在移除区段之后计数。
pub fn relocate(
    state: &mut State,
    group: usize,
    src: usize,
    start: usize,
    count: usize,
    dst: usize,
    pos: usize,
) -> Result<Vec<Edit>, String> {
    let lists = state.lists.get(group).ok_or("列表组编号越界")?;
    let from = lists.get(src).ok_or("来源列表编号越界")?;
    let to = lists.get(dst).ok_or("目标列表编号越界")?;
    if count == 0 || start > from.len() || count > from.len() - start {
        return Err("搬移区段越界或为空".into());
    }
    let remaining = if src == dst {
        to.len() - count
    } else {
        to.len()
    };
    if pos > remaining {
        return Err("搬入位置越界".into());
    }
    if src == dst && pos == start {
        return Ok(vec![]);
    }
    let removed = splice(state, group, src, start, count, &[])?;
    let inserted = splice(state, group, dst, pos, 0, &removed.removed)?;
    Ok(vec![removed, inserted])
}

/// 同列表交换从较大位置开始移除，再从较小位置开始插回。中间状态也没有重复成员。
pub fn swap(
    state: &mut State,
    group: usize,
    a: usize,
    i: usize,
    b: usize,
    j: usize,
) -> Result<Vec<Edit>, String> {
    let lists = state.lists.get(group).ok_or("列表组编号越界")?;
    let first = *lists.get(a).and_then(|r| r.get(i)).ok_or("交换位置越界")?;
    let second = *lists.get(b).and_then(|r| r.get(j)).ok_or("交换位置越界")?;
    if a == b && i == j {
        return Ok(vec![]);
    }
    let ((lo, p, x), (hi, q, y)) = if a == b && i > j {
        ((b, j, second), (a, i, first))
    } else {
        ((a, i, first), (b, j, second))
    };
    Ok(vec![
        splice(state, group, hi, q, 1, &[])?,
        splice(state, group, lo, p, 1, &[])?,
        splice(state, group, lo, p, 0, &[y])?,
        splice(state, group, hi, q, 0, &[x])?,
    ])
}

fn demand_sum(items: &[usize], demands: &[i64]) -> Result<i128, String> {
    items.iter().try_fold(0i128, |sum, &item| {
        Ok(sum + i128::from(*demands.get(item).ok_or("载重项目编号越界")?))
    })
}

#[cfg(test)]
fn load(route: &[usize], demands: &[i64]) -> Result<i64, String> {
    i64::try_from(demand_sum(route, demands)?).map_err(|_| "载重合计溢出".into())
}

pub fn load_delta(edit: &Edit, demands: &[i64]) -> Result<i128, String> {
    Ok(demand_sum(&edit.added, demands)? - demand_sum(&edit.removed, demands)?)
}

/// 仅两边都不存在且区段为空时，车辆没有出库，即使仓库对角元素非零费用也为零。
fn path_cost(
    left: Option<usize>,
    items: &[usize],
    right: Option<usize>,
    matrix: &[i64],
    n: usize,
) -> Result<i128, String> {
    let width = n.checked_add(1).ok_or("距离矩阵过大")?;
    if width.checked_mul(width) != Some(matrix.len()) {
        return Err("距离矩阵维数错误".into());
    }
    if left.is_none() && items.is_empty() && right.is_none() {
        return Ok(0);
    }
    let mut prev = left.unwrap_or(n);
    if left.is_some_and(|x| x >= n) || right.is_some_and(|x| x >= n) {
        return Err("路段编号越界".into());
    }
    let mut sum = 0i128;
    for &item in items {
        if item >= n {
            return Err("路段编号越界".into());
        }
        sum += i128::from(matrix[prev * width + item]);
        prev = item;
    }
    Ok(sum + i128::from(matrix[prev * width + right.unwrap_or(n)]))
}

#[cfg(test)]
fn distance(route: &[usize], matrix: &[i64], n: usize) -> Result<i64, String> {
    i64::try_from(path_cost(None, route, None, matrix, n)?).map_err(|_| "路线距离合计溢出".into())
}

pub fn distance_delta(edit: &Edit, matrix: &[i64], n: usize) -> Result<i128, String> {
    let new = path_cost(edit.left, &edit.added, edit.right, matrix, n)?;
    let old = path_cost(edit.left, &edit.removed, edit.right, matrix, n)?;
    Ok(new - old)
}

/// 只统计业务元素之间的有向相邻边。空区段将两侧邻居直接相连。
pub fn adjacent_delta(edit: &Edit, before: usize, after: usize) -> i128 {
    let matches = |items: &[usize]| {
        let mut previous = edit.left;
        let mut count = 0i128;
        for next in items.iter().copied().map(Some).chain([edit.right]) {
            if previous == Some(before) && next == Some(after) {
                count += 1;
            }
            previous = next;
        }
        count
    };
    matches(&edit.added) - matches(&edit.removed)
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Link {
    owner: Option<usize>,
    prev: Option<usize>,
    next: Option<usize>,
}

/// 静态边为工艺依赖，动态边为同机相邻关系；每道工序只保存一个动态前驱与后继。
#[derive(Debug, Clone)]
pub struct Schedule {
    pred: Vec<Vec<usize>>,
    succ: Vec<Vec<usize>>,
    links: Vec<Link>,
    finish: Vec<Option<i128>>,
}

impl Schedule {
    pub fn new(spec: &ListSpec, lists: &[Vec<usize>]) -> Result<Self, String> {
        let n = spec.eligible.len();
        if lists.len() != spec.owners {
            return Err("列表个数错误".into());
        }
        if (!spec.durations.is_empty() && spec.durations.len() != n)
            || spec.durations.iter().any(|&d| d < 0)
        {
            return Err("工序时长无效".into());
        }
        let mut s = Self {
            pred: vec![vec![]; n],
            succ: vec![vec![]; n],
            links: vec![Link::default(); n],
            finish: vec![None; n],
        };
        for &(x, y) in &spec.precedences {
            if x >= n || y >= n {
                return Err("工艺依赖编号越界".into());
            }
            s.pred[y].push(x);
            s.succ[x].push(y);
        }
        let mut assigned = BTreeSet::new();
        for (owner, list) in lists.iter().enumerate() {
            for (pos, &x) in list.iter().enumerate() {
                if x >= n || !spec.eligible[x].contains(&owner) || !assigned.insert(x) {
                    return Err("列表含有重复、越界或不合格的成员".into());
                }
                s.links[x] = Link {
                    owner: Some(owner),
                    prev: pos.checked_sub(1).map(|p| list[p]),
                    next: list.get(pos + 1).copied(),
                };
            }
        }
        s.propagate(spec, &assigned)?;
        Ok(s)
    }

    pub fn finish(&self, item: usize) -> Option<i128> {
        self.finish.get(item).copied().flatten()
    }

    /// 失败时同时保留原邻接关系和完工时间；调用方只需撤销工作解中的编辑。
    pub fn update(&mut self, spec: &ListSpec, edits: &[Edit]) -> Result<Vec<usize>, String> {
        let mut journal = BTreeMap::new();
        let result = (|| {
            for edit in edits {
                self.apply_edit(spec, edit, &mut journal)?;
            }
            let seeds = journal
                .iter()
                .filter_map(|(&x, old)| {
                    let new = self.links[x];
                    (old.prev != new.prev || old.owner.is_some() != new.owner.is_some())
                        .then_some(x)
                })
                .collect();
            self.propagate(spec, &seeds)
        })();
        if result.is_err() {
            for (x, old) in journal {
                self.links[x] = old;
            }
        }
        result
    }

    fn apply_edit(
        &mut self,
        spec: &ListSpec,
        edit: &Edit,
        journal: &mut BTreeMap<usize, Link>,
    ) -> Result<(), String> {
        if edit.owner >= spec.owners {
            return Err("列表编号越界".into());
        }
        for &x in edit
            .left
            .iter()
            .chain(&edit.right)
            .chain(&edit.removed)
            .chain(&edit.added)
        {
            let link = *self.links.get(x).ok_or("列表成员编号越界")?;
            journal.entry(x).or_insert(link);
        }
        for x in edit.left.iter().chain(&edit.right).chain(&edit.removed) {
            if self.links[*x].owner != Some(edit.owner) {
                return Err("编辑的来源列表不匹配".into());
            }
        }
        // 编辑记录可以跨试探/撤销传递。拒绝过期的邻居信息，避免静默破坏拓扑缓存。
        let mut left = edit.left;
        for right in edit.removed.iter().copied().map(Some).chain([edit.right]) {
            if left.is_some_and(|x| self.links[x].next != right)
                || right.is_some_and(|x| self.links[x].prev != left)
            {
                return Err("列表邻接关系与编辑记录不一致".into());
            }
            left = right;
        }
        for &x in &edit.removed {
            self.links[x] = Link::default();
        }
        for (p, &x) in edit.added.iter().enumerate() {
            if self.links[x].owner.is_some() || !spec.eligible[x].contains(&edit.owner) {
                return Err("重复分配或列表不在候选集合中".into());
            }
            self.links[x] = Link {
                owner: Some(edit.owner),
                prev: p.checked_sub(1).map(|i| edit.added[i]).or(edit.left),
                next: edit.added.get(p + 1).copied().or(edit.right),
            };
        }
        if let Some(x) = edit.left {
            self.links[x].next = edit.added.first().copied().or(edit.right);
        }
        if let Some(x) = edit.right {
            self.links[x].prev = edit.added.last().copied().or(edit.left);
        }
        Ok(())
    }

    /// Kahn 拓扑遍历仅处理前驱改变的工序及其后继。边界前驱直接读取原缓存。
    /// 先计算到局部暂存表，确认整个受影响子图无环后再一次提交完工时间。
    fn propagate(
        &mut self,
        spec: &ListSpec,
        seeds: &BTreeSet<usize>,
    ) -> Result<Vec<usize>, String> {
        let mut affected = seeds.clone();
        let mut stack: Vec<_> = seeds.iter().copied().collect();
        while let Some(x) = stack.pop() {
            for &y in self.succ[x].iter().chain(self.links[x].next.iter()) {
                if self.links[y].owner.is_some() && affected.insert(y) {
                    stack.push(y);
                }
            }
        }

        let mut ready = VecDeque::new();
        let mut remaining = BTreeMap::new();
        for &x in &affected {
            let count = if self.links[x].owner.is_none() {
                0
            } else {
                self.pred[x]
                    .iter()
                    .chain(self.links[x].prev.iter())
                    .filter(|y| affected.contains(y) && self.links[**y].owner.is_some())
                    .count()
            };
            if count == 0 {
                ready.push_back(x);
            }
            remaining.insert(x, count);
        }
        let mut computed = BTreeMap::new();
        while let Some(x) = ready.pop_front() {
            if self.links[x].owner.is_none() {
                computed.insert(x, None);
                continue;
            }
            let mut start = 0i128;
            for &p in self.pred[x].iter().chain(self.links[x].prev.iter()) {
                if self.links[p].owner.is_some() {
                    let end = computed
                        .get(&p)
                        .copied()
                        .unwrap_or(self.finish[p])
                        .ok_or("工序前驱缺少完工时间")?;
                    start = start.max(end);
                }
            }
            // 一条路径最多包含 usize::MAX 道工序，每道时长不超过 i64::MAX。
            let duration = i128::from(spec.durations.get(x).copied().unwrap_or(0));
            let end = start + duration;
            computed.insert(x, Some(end));
            for &y in self.succ[x].iter().chain(self.links[x].next.iter()) {
                if self.links[y].owner.is_some() {
                    let count = remaining.get_mut(&y).ok_or("工序后继缺少拓扑计数")?;
                    *count -= 1;
                    if *count == 0 {
                        ready.push_back(y);
                    }
                }
            }
        }
        if computed.len() != affected.len() {
            return Err(CYCLE.into());
        }
        let mut changed = Vec::new();
        for (x, end) in computed {
            if self.finish[x] != end {
                self.finish[x] = end;
                changed.push(x);
            }
        }
        Ok(changed)
    }
}

#[cfg(test)]
fn acyclic(spec: &ListSpec, lists: &[Vec<usize>]) -> Result<bool, String> {
    match Schedule::new(spec, lists) {
        Ok(_) => Ok(true),
        Err(e) if e == CYCLE => Ok(false),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests;
