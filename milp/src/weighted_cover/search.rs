//! 非负行价格给出集合覆盖的拉格朗日下界。各列可独立决定选或不选，
//! 故在一个固定价格下，强制改变某列选择的额外费用为其约化费用的绝对值。

use std::cmp::Reverse;
use std::time::Instant;

const SCALE: i64 = 1 << 20;
const MAX_PRICE: f64 = 1_000_000.;
const SMOOTHING: [(f64, usize); 4] = [(1.5, 3), (0.3, 4), (0.06, 6), (0.012, 8)];

const EXCLUDED: i8 = -1;
const FREE: i8 = 0;
const SELECTED: i8 = 1;

struct Node {
    rows: Vec<Vec<usize>>,
    // 当前列在原模型中的位置。
    columns: Vec<usize>,
    costs: Vec<i64>,
    prices: Vec<f64>,
}

struct Search<'a> {
    end: Option<Instant>,
    limit: u64,
    nodes: u64,
    root_bound: i64,
    upper: i64,
    incumbent: &'a mut Vec<bool>,
    update: &'a mut dyn FnMut(&[bool]),
}

impl Node {
    // 浮点运算只用于挑选行价格。令每列的松弛值为 max(0, -r/宽度)，
    // 其中 r 是约化费用。对单行调整价格，直到该行的松弛覆盖量为一。
    fn improve(
        &mut self,
        width: f64,
        passes: usize,
        residual: &mut [f64],
        end: Option<Instant>,
    ) -> bool {
        let mut candidates = Vec::new();
        for t in 0..passes {
            if end.is_some_and(|d| Instant::now() >= d) {
                return false;
            }
            for k in 0..self.rows.len() {
                let i = if t % 2 == 0 {
                    k
                } else {
                    self.rows.len() - 1 - k
                };
                let row = &self.rows[i];
                if self.prices[i] == 0. && row.iter().any(|&j| residual[j] <= -width) {
                    continue;
                }
                let first = row
                    .iter()
                    .map(|&j| residual[j])
                    .fold(f64::INFINITY, f64::min);
                // 其他列若比最小值大至少一个宽度，则最便宜列独自贡献一单位。
                candidates.clear();
                for &j in row {
                    if residual[j] < first + width {
                        candidates.push(residual[j] - first);
                    }
                }
                candidates.sort_unstable_by(f64::total_cmp);
                let mut water = width;
                let mut count = 1;
                for &val in candidates.iter().skip(1) {
                    if val * count as f64 >= water {
                        break;
                    }
                    water += val;
                    count += 1;
                }
                let delta = (first + water / count as f64).max(-self.prices[i]);
                if delta == 0. {
                    continue;
                }
                self.prices[i] += delta;
                for &j in row {
                    residual[j] -= delta;
                }
            }
        }
        true
    }

    // 对行价格取整，从整数系数计算 fixed + Σy + Σmin(0,c-Aᵀy)。
    // 入口最多 1024 行和列、10 万个非零项。各列原费用与每个行价格均 ≤ 10⁶，
    // 乘 SCALE 后，绝对值总和 ≤ (100_000 + 3 × 1024) × 10⁶ × SCALE < 2⁵⁷。
    fn certificate(&self, fixed: i64) -> (i64, Vec<i64>) {
        let mut cost: Vec<_> = self.costs.iter().map(|c| c * SCALE).collect();
        let mut value = fixed * SCALE;
        let mut prices = Vec::with_capacity(self.rows.len());
        for (r, &p) in self.rows.iter().zip(&self.prices) {
            let y = (p.clamp(0., MAX_PRICE) * SCALE as f64).round() as i64;
            prices.push(y);
            value += y;
            for &j in r {
                cost[j] -= y;
            }
        }
        // 逐行把最小的正约化费用降到零，或令第二小的负约化费用升到零。
        // 更新后行价格仍非负；一行至多提高到所含列的原费用，整数证书仍有效。
        for (row, price) in self.rows.iter().zip(prices) {
            let mut first = i64::MAX;
            let mut second = i64::MAX;
            for &j in row {
                let c = cost[j];
                if c < first {
                    second = first;
                    first = c;
                } else {
                    second = second.min(c);
                }
            }
            let delta = if first > 0 {
                first
            } else {
                second.min(0).max(-price)
            };
            if delta != 0 {
                value += delta;
                for &j in row {
                    cost[j] -= delta;
                }
            }
        }
        (value + cost.iter().map(|&x| x.min(0)).sum::<i64>(), cost)
    }

    fn child(&self, decisions: &[i8], solution: &mut [bool], fixed: &mut i64) -> Self {
        let mut map = vec![0; self.columns.len()];
        let mut columns = Vec::new();
        let mut costs = Vec::new();
        for (j, &x) in decisions.iter().enumerate() {
            match x {
                SELECTED => {
                    solution[self.columns[j]] = true;
                    *fixed += self.costs[j];
                }
                FREE => {
                    map[j] = columns.len();
                    columns.push(self.columns[j]);
                    costs.push(self.costs[j]);
                }
                _ => {}
            }
        }
        let mut rows = Vec::new();
        let mut prices = Vec::new();
        for (i, r) in self.rows.iter().enumerate() {
            if r.iter().any(|&j| decisions[j] == SELECTED) {
                continue;
            }
            rows.push(
                r.iter()
                    .filter(|&&j| decisions[j] == FREE)
                    .map(|&j| map[j])
                    .collect(),
            );
            prices.push(self.prices[i]);
        }
        Self {
            rows,
            prices,
            costs,
            columns,
        }
    }
}

impl Search<'_> {
    fn visit(
        &mut self,
        mut node: Node,
        mut solution: Vec<bool>,
        mut fixed: i64,
        root: bool,
    ) -> bool {
        if self.nodes >= self.limit || self.end.is_some_and(|d| Instant::now() >= d) {
            return false;
        }
        self.nodes += 1;
        let mut reuse_prices = false;
        loop {
            if self.end.is_some_and(|d| Instant::now() >= d) {
                return false;
            }
            if fixed >= self.upper || node.rows.iter().any(Vec::is_empty) {
                return true;
            }
            if node.rows.is_empty() {
                self.upper = fixed;
                (self.update)(&solution);
                *self.incumbent = solution;
                return true;
            }
            let mut residual: Vec<_> = node.costs.iter().map(|&c| c as f64).collect();
            for (r, &p) in node.rows.iter().zip(&node.prices) {
                for &j in r {
                    residual[j] -= p;
                }
            }
            let target = (self.upper - 1) * SCALE;
            let mut decisions = vec![FREE; node.columns.len()];
            for r in &node.rows {
                if r.len() == 1 {
                    decisions[r[0]] = SELECTED;
                }
            }
            for (width, passes) in SMOOTHING {
                if !reuse_prices && !node.improve(width, passes, &mut residual, self.end) {
                    return false;
                }
                let (value, reduced) = node.certificate(fixed);
                if root {
                    self.root_bound = self.root_bound.max((value + SCALE - 1).div_euclid(SCALE));
                }
                if value > target {
                    return true;
                }
                for (j, &c) in reduced.iter().enumerate() {
                    if value + c.abs() > target {
                        let force = if c > 0 { EXCLUDED } else { SELECTED };
                        if decisions[j] == -force {
                            return true;
                        }
                        decisions[j] = force;
                    }
                }
                let changed = decisions.iter().filter(|&&a| a != FREE).count();
                if reuse_prices || changed * 4 > node.columns.len() {
                    break;
                }
            }
            let changed = decisions.iter().filter(|&&d| d != FREE).count();
            if changed > 0 {
                // 只固定少数列时先沿用行价格，避免在几乎相同的子问题上重复松弛。
                reuse_prices = changed * 8 < node.columns.len();
                node = node.child(&decisions, &mut solution, &mut fixed);
                continue;
            }
            // 先分支候选列最少的行；同样长时，优先选择能覆盖更多其余行的列所在行。
            let mut degree = vec![0usize; node.columns.len()];
            for row in &node.rows {
                for &j in row {
                    degree[j] += 1;
                }
            }
            let row = node
                .rows
                .iter()
                .min_by_key(|r| {
                    (
                        r.len(),
                        Reverse(r.iter().map(|&j| degree[j]).sum::<usize>()),
                    )
                })
                .unwrap();
            let mut order = row.clone();
            order.sort_unstable_by(|&a, &b| residual[a].total_cmp(&residual[b]).then(a.cmp(&b)));
            // 每个子树选择该行中第一个被选的列；此前选项全部置零，分支不遗漏也不重复。
            for j in order {
                decisions[j] = SELECTED;
                let mut x = solution.clone();
                let mut cost = fixed;
                let child = node.child(&decisions, &mut x, &mut cost);
                if !self.visit(child, x, cost, false) {
                    return false;
                }
                decisions[j] = EXCLUDED;
            }
            return true;
        }
    }
}

// 若 k 费用不高于 j 且覆盖 j 的每一行，则选择 j 可替换为选择 k。
// 费用相同且覆盖也相同的列保留最小下标。按这个顺序，连续替换必终止。
fn dominated(rows: &[Vec<usize>], cols: &[Vec<usize>], costs: &[i64]) -> Vec<i8> {
    let words = rows.len().div_ceil(64);
    let mut masks = vec![vec![0u64; words]; cols.len()];
    for (mask, column) in masks.iter_mut().zip(cols) {
        for &i in column {
            mask[i / 64] |= 1 << (i % 64);
        }
    }
    let mut decisions = vec![FREE; cols.len()];
    for (j, column) in cols.iter().enumerate() {
        let Some(&i) = column.iter().min_by_key(|&&i| rows[i].len()) else {
            decisions[j] = EXCLUDED;
            continue;
        };
        let rank = |k: usize| (costs[k], Reverse(cols[k].len()), k);
        if rows[i].iter().any(|&k| {
            rank(k) < rank(j)
                && cols[k].len() >= column.len()
                && masks[j].iter().zip(&masks[k]).all(|(a, b)| a & !b == 0)
        }) {
            decisions[j] = EXCLUDED;
        }
    }
    decisions
}

pub(super) fn search(
    rows: &[Vec<usize>],
    cols: &[Vec<usize>],
    costs: &[i64],
    incumbent: &mut Vec<bool>,
    deadline: Option<Instant>,
    limit: u64,
    update: &mut dyn FnMut(&[bool]),
) -> (bool, u64, i64) {
    if limit == 0 || deadline.is_some_and(|d| Instant::now() >= d) {
        return (false, 0, 0);
    }
    let upper = costs
        .iter()
        .zip(incumbent.iter())
        .filter(|(_, v)| **v)
        .map(|(c, _)| c)
        .sum();
    let mut s = Search {
        end: deadline,
        limit,
        nodes: 0,
        upper,
        root_bound: 0,
        incumbent,
        update,
    };
    let n = costs.len();
    let node = Node {
        rows: rows.to_vec(),
        columns: (0..n).collect(),
        costs: costs.to_vec(),
        prices: vec![0.; rows.len()],
    };
    let mut x = vec![false; n];
    let mut fixed = 0;
    let decisions = dominated(rows, cols, costs);
    let node = if decisions.contains(&EXCLUDED) {
        node.child(&decisions, &mut x, &mut fixed)
    } else {
        node
    };
    let done = s.visit(node, x, fixed, true);
    // 根上的固定针对费用严格小于当时上界的解。即便固定后的模型为空，
    // 所有未考虑的解费用仍不低于当时上界；对外下界以当前可行费用截断。
    (
        done,
        s.nodes,
        if done {
            s.upper
        } else {
            s.root_bound.min(s.upper)
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_certificate_bounds_each_opposing_fixing() {
        let node = Node {
            rows: vec![vec![0, 1], vec![0, 2], vec![1, 2, 3], vec![3, 4]],
            costs: vec![9, 4, 3, 1_000_000, 5],
            columns: (0..5).collect(),
            prices: vec![1.25, 500_000.23, 1_000_000., 2.5],
        };
        let fixed = 17;
        let (bound, reduced) = node.certificate(fixed);
        for mask in 0..1 << node.columns.len() {
            let picked = |j| mask & (1usize << j) != 0usize;
            if !node.rows.iter().all(|r| r.iter().any(|&j| picked(j))) {
                continue;
            }
            let total = fixed
                + node
                    .costs
                    .iter()
                    .enumerate()
                    .filter_map(|(j, &c)| picked(j).then_some(c))
                    .sum::<i64>();
            // 非负约化费用列倾向取零，负约化费用列倾向取一。
            let contrary: i64 = reduced
                .iter()
                .enumerate()
                .filter_map(|(j, &c)| (picked(j) != (c < 0)).then_some(c.abs()))
                .sum();
            assert!(bound + contrary <= total * SCALE);
        }
    }

    #[test]
    fn dominance_keeps_one_representative_and_checks_all_rows() {
        let rows = vec![vec![0, 1, 2, 3], vec![0, 1, 2, 4], vec![2]];
        let columns = vec![
            vec![0, 1],
            vec![0, 1],
            vec![0, 1, 2],
            vec![0],
            vec![1],
            vec![],
        ];
        // 列 2 是超集但较贵；列 4 费用便宜，覆盖范围较小，二者都必须保留。
        assert_eq!(
            dominated(&rows, &columns, &[3, 3, 4, 3, 1, 1]),
            vec![FREE, EXCLUDED, FREE, EXCLUDED, FREE, EXCLUDED]
        );
    }
}
