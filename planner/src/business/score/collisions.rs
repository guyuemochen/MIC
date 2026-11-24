//! 同一组内每对键都比较一次时，值相同的 n 个键恰好贡献 n*(n-1)/2。
//! 键是原表达式节点，未知值不入桶，与逐对比较后求和的构造期语义一致。

use std::collections::{BTreeMap, HashMap};

use super::{Expr, Model, Value};

pub(super) struct PairCounts {
    counts: Vec<BTreeMap<i64, usize>>,
    total: i128,
    // 按每个键在原求和中首次比较的顺序保留，多个键同时出错时仍返回原来的错误。
    keys: Vec<usize>,
}

impl PairCounts {
    pub(super) fn compile(
        model: &Model,
        expr: &Expr,
        values: &[Value],
        mut subscribe: impl FnMut(usize, usize),
    ) -> Option<Self> {
        let Expr::Sum(terms) = expr else { return None };
        if terms.len() < 3 {
            return None;
        }
        let mut ids = HashMap::new();
        let mut keys = Vec::new();
        let mut neighbors = Vec::<Vec<usize>>::new();
        let mut intern = |key| {
            *ids.entry(key).or_insert_with(|| {
                let index = keys.len();
                keys.push(key);
                index
            })
        };
        for &term in terms {
            let Expr::Eq(a, b) = model.nodes[term] else {
                return None;
            };
            if a == b {
                return None;
            }
            let (a, b) = (intern(a), intern(b));
            if neighbors.len() <= a.max(b) {
                neighbors.resize_with(a.max(b) + 1, Vec::new);
            }
            neighbors[a].push(b);
            neighbors[b].push(a);
        }
        // 互不相交的单条比较没有可复用的计数，保留普通聚合。
        if neighbors.iter().all(|edges| edges.len() == 1) {
            return None;
        }
        let mut groups = vec![usize::MAX; keys.len()];
        let mut num_groups = 0;
        let mut component = Vec::new();
        for start in 0..keys.len() {
            if groups[start] != usize::MAX {
                continue;
            }
            groups[start] = num_groups;
            component.clear();
            component.push(start);
            let mut next = 0;
            while next < component.len() {
                let key = component[next];
                next += 1;
                let edges = &mut neighbors[key];
                edges.sort_unstable();
                if edges.windows(2).any(|pair| pair[0] == pair[1]) {
                    return None;
                }
                for &other in edges.iter() {
                    if groups[other] == usize::MAX {
                        groups[other] = num_groups;
                        component.push(other);
                    }
                }
            }
            // 排除了自环和重边后，每个顶点都须与本分量其余所有顶点相连。
            if component
                .iter()
                .any(|&key| neighbors[key].len() != component.len() - 1)
            {
                return None;
            }
            num_groups += 1;
        }
        let mut result = Self {
            counts: (0..num_groups).map(|_| BTreeMap::new()).collect(),
            total: 0,
            keys,
        };
        for (i, group) in groups.into_iter().enumerate() {
            let key = result.keys[i];
            subscribe(key, group);
            result.replace(group, None, values[key]);
        }
        Some(result)
    }

    pub(super) fn replace(&mut self, group: usize, old: Value, new: Value) {
        let counts = &mut self.counts[group];
        if let Some(v) = old {
            let n = counts.get_mut(&v).expect("旧冲突桶应有正频次");
            *n -= 1;
            self.total -= *n as i128;
            if *n == 0 {
                counts.remove(&v);
            }
        }
        if let Some(v) = new {
            let n = counts.entry(v).or_default();
            self.total += *n as i128;
            *n += 1;
        }
    }

    pub(super) fn total(&self) -> i128 {
        self.total
    }

    pub(super) fn fault(&self, faults: &super::Faults) -> Option<String> {
        self.keys.iter().find_map(|key| faults.get(key)).cloned()
    }
}
