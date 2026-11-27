//! 在非负装箱列中保留低检验数候选，反复撤掉部分已选列，在释放的容量内重新组合。
//! 邻域搜索中的收益界仅用于挑选可行候选；调用者还原后在完整模型验收。
use std::time::Instant;

pub(super) struct Item {
    pub column: usize,
    pub profit: f64,
    pub rows: Vec<(usize, f64)>,
}

struct Choice<'a> {
    item: &'a Item,
    root: f64,
    cost: f64,
    priority: f64,
}

struct Search<'a> {
    choices: Vec<Choice<'a>>,
    load: Vec<f64>,
    taken: Vec<bool>,
    path: Vec<usize>,
    completion: Vec<usize>,
    completion_value: f64,
    nodes: usize,
    used_nodes: u64,
    deadline: Instant,
}

impl Search<'_> {
    fn fits(&self, j: usize) -> bool {
        self.choices[j]
            .item
            .rows
            .iter()
            .all(|&(i, w)| self.load[i] + w <= 1.)
    }

    fn change(&mut self, j: usize, add: bool) {
        self.taken[j] = add;
        for &(i, w) in &self.choices[j].item.rows {
            if add {
                self.load[i] += w;
            } else {
                self.load[i] -= w;
            }
        }
    }

    // 同一资源中每项超过剩余容量一半时至多选一项。把每列归到一个这样的组，
    // 各组的最高收益相加即可估计余下改进空间。搜索节点封顶只影响候选质量。
    fn complete(&mut self, available: &[usize], value: f64) {
        if self.nodes == 0 || (self.nodes % 32 == 0 && Instant::now() >= self.deadline) {
            self.nodes = 0;
            return;
        }
        self.nodes -= 1;
        self.used_nodes += 1;
        if value > self.completion_value {
            self.completion_value = value;
            self.completion.clone_from(&self.path);
        }
        if available.is_empty() || self.path.len() >= 64 {
            return;
        }
        let mut counts = vec![0usize; self.load.len()];
        for &j in available {
            for &(i, w) in &self.choices[j].item.rows {
                if 2. * w > 1. - self.load[i] {
                    counts[i] += 1;
                }
            }
        }
        let mut top = vec![0.0_f64; counts.len()];
        let mut cap = value;
        let mut most = 0;
        let mut pivot = None;
        for &j in available {
            let item = self.choices[j].item;
            let mut through = None;
            let mut degree = 0;
            for &(i, w) in &item.rows {
                if 2. * w > 1. - self.load[i] && counts[i] > degree {
                    degree = counts[i];
                    through = Some(i);
                }
            }
            if let Some(i) = through {
                top[i] = top[i].max(item.profit);
                if degree > most {
                    most = degree;
                    pivot = Some(i);
                }
            } else {
                cap += item.profit;
            }
        }
        cap += top.iter().sum::<f64>();
        if cap <= self.completion_value {
            return;
        }

        let (mut alternatives, mut others): (Vec<_>, Vec<_>) =
            available.iter().copied().partition(|&j| {
                pivot.is_some_and(|i| {
                    self.choices[j]
                        .item
                        .rows
                        .iter()
                        .any(|&(r, w)| r == i && 2. * w > 1. - self.load[i])
                })
            });
        if alternatives.is_empty() {
            alternatives.push(others.pop().unwrap());
        }
        alternatives.sort_unstable_by(|&a, &b| {
            self.choices[b]
                .item
                .profit
                .total_cmp(&self.choices[a].item.profit)
                .then(a.cmp(&b))
        });
        for j in alternatives {
            if self.nodes == 0 {
                return;
            }
            self.change(j, true);
            self.path.push(j);
            let next: Vec<_> = others.iter().copied().filter(|&k| self.fits(k)).collect();
            self.complete(&next, value + self.choices[j].item.profit);
            self.path.pop();
            self.change(j, false);
        }
        self.complete(&others, value);
    }
}

// 固定种子与列号破平局，使相同工作预算下的搜索顺序可复现。
fn sample(state: &mut u64) -> f64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    (*state >> 11) as f64 / (1u64 << 53) as f64
}

/// 每项的正权重以对应资源容量为单位。返回选择的原工作列号，每列最多取一次。
pub(super) fn selected(
    resources: usize,
    items: &[Item],
    root: &[f64],
    reduced_costs: &[f64],
    deadline: Instant,
    node_limit: u64,
) -> (Vec<usize>, u64) {
    if resources == 0
        || resources > 4096
        || items.is_empty()
        || node_limit == 0
        || Instant::now() >= deadline
    {
        return (Vec::new(), 0);
    }
    // 用典型列收益定标，整体放大目标函数时仍采用相同的检验数筛选和搜索温度。
    let mut gains: Vec<_> = items.iter().map(|it| it.profit).collect();
    let mid = gains.len() / 2;
    let scale = *gains.select_nth_unstable_by(mid, |a, b| a.total_cmp(b)).1;
    if !scale.is_finite() || scale <= 0. {
        return (Vec::new(), 0);
    }
    let mut choices: Vec<_> = items
        .iter()
        .filter_map(|item| {
            let x = *root.get(item.column)?;
            let cost = reduced_costs.get(item.column)?.max(0.) / scale;
            (cost < 0.06 || x > 1e-7).then(|| Choice {
                item,
                root: x,
                cost,
                priority: (item.profit / scale).ln() - 0.4 * (item.rows.len().max(1) as f64).ln(),
            })
        })
        .collect();
    // 限制一次局部搜索的排序和枚举代价。先保留根支撑，再按检验数和原列号截断。
    const MAX_CHOICES: usize = 3072;
    if choices.len() > MAX_CHOICES {
        choices.select_nth_unstable_by(MAX_CHOICES, |a, b| {
            (b.root > 1e-7)
                .cmp(&(a.root > 1e-7))
                .then(a.cost.total_cmp(&b.cost))
                .then(a.item.column.cmp(&b.item.column))
        });
        choices.truncate(MAX_CHOICES);
        choices.sort_unstable_by_key(|it| it.item.column);
    }
    if choices.is_empty() {
        return (Vec::new(), 0);
    }
    let n = choices.len();
    let mut s = Search {
        choices,
        load: vec![0.; resources],
        taken: vec![false; n],
        path: Vec::new(),
        completion: Vec::new(),
        completion_value: 0.,
        nodes: 0,
        used_nodes: 0,
        deadline,
    };
    let mut random = 1;
    let mut best = Vec::new();
    let mut current = Vec::new();
    let mut best_value = 0.;
    let mut current_value = 0.;
    let mut scores = vec![0.; n];
    let mut order: Vec<usize> = (0..n).collect();
    for trial in 0..2048 {
        if Instant::now() >= deadline || s.used_nodes >= node_limit {
            break;
        }
        s.load.fill(0.);
        s.taken.fill(false);
        let destroy = if trial % 13 == 0 {
            1.
        } else if trial % 3 == 0 {
            0.28
        } else {
            0.12
        };
        let mut point = Vec::new();
        for &j in if trial % 5 == 0 { &best } else { &current } {
            if sample(&mut random) > destroy {
                point.push(j);
                s.change(j, true);
            }
        }
        if trial > 4 && trial % 3 != 0 {
            let available: Vec<_> = (0..n).filter(|&j| !s.taken[j] && s.fits(j)).collect();
            s.completion_value = 0.;
            s.completion.clear();
            s.nodes = (node_limit - s.used_nodes).min(500) as usize;
            s.complete(&available, 0.);
            point.extend(&s.completion);
        } else {
            let temperature = (0.25 + 1.75 * sample(&mut random)) / 100.;
            let noise = sample(&mut random) * 0.7;
            let bias = sample(&mut random) * 0.4;
            for (j, choice) in s.choices.iter().enumerate() {
                scores[j] = choice.priority
                    + noise * sample(&mut random).max(1e-6).ln()
                    + bias * choice.root
                    - choice.cost / temperature;
            }
            order.sort_unstable_by(|&a, &b| scores[b].total_cmp(&scores[a]).then(a.cmp(&b)));
            for &j in &order {
                if !s.taken[j] && s.fits(j) {
                    s.change(j, true);
                    point.push(j);
                }
            }
        }
        let value: f64 = point.iter().map(|&j| s.choices[j].item.profit).sum();
        if value > best_value {
            best_value = value;
            best.clone_from(&point);
        }
        if value >= current_value
            || trial % 31 == 0
            || sample(&mut random)
                < 0.06 * ((value - current_value) / (scale * 0.035)).min(0.).exp()
        {
            current = point;
            current_value = value;
        } else if trial % 13 == 0 {
            current.clone_from(&best);
            current_value = best_value;
        }
    }
    (
        best.into_iter().map(|j| s.choices[j].item.column).collect(),
        s.used_nodes,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn enumeration_obeys_shared_node_budget_and_deadline() {
        let items = [
            Item {
                column: 0,
                profit: 3.,
                rows: vec![(0, 1.)],
            },
            Item {
                column: 1,
                profit: 1.,
                rows: vec![(0, 1.)],
            },
        ];
        let root = [0.; 2];
        let costs = [0.; 2];
        let end = Instant::now() + Duration::from_secs(1);
        assert_eq!(selected(1, &items, &root, &costs, end, 0), (vec![], 0));
        assert_eq!(
            selected(1, &items, &root, &costs, Instant::now(), 10),
            (vec![], 0),
        );
        let (choice, nodes) = selected(1, &items, &root, &costs, end, 1);
        assert_eq!(choice, vec![0]);
        assert_eq!(nodes, 1);
    }
}
