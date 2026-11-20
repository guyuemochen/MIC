//! 用覆盖费用的拉格朗日松弛计算下界。浮点迭代只挑选行乘子；每次发布前，
//! 将乘子换成非负定点数，在原 0/1 盒上用整数重算证书。
use super::expired;
use std::collections::BTreeSet;
use std::time::Instant;

const SCALE: i64 = 1 << 22;

struct Row {
    cols: Vec<usize>,
    rhs: i64,
    price: f32,
    step: f32,
}

impl Row {
    fn new(cols: Vec<usize>, rhs: i64) -> Self {
        let step = 0.3 / cols.len() as f32;
        Self {
            cols,
            rhs,
            price: 0.,
            step,
        }
    }
}

fn certificate(rows: &[Row], n: usize) -> usize {
    let mut reduced = vec![SCALE; n];
    let mut value = 0;
    for r in rows {
        // 乘子截在 [0, 2] 仍给出有效下界。即使所有 20 万原始项和新割同时取最大值，
        // 各列剩余价格之和与 rhs 贡献也远在 i64 范围内。
        let y = (r.price.clamp(0., 2.) * SCALE as f32).round() as i64;
        value += y * r.rhs;
        for &j in &r.cols {
            reduced[j] -= y;
        }
    }
    // min_{x∈{0,1}} (1-Σ y_i)x = min(0,1-Σ y_i)。原变量代价为 1，最优费用为整数。
    value += reduced.into_iter().map(|c| c.min(0)).sum::<i64>();
    ((value.max(0) + SCALE - 1) / SCALE) as usize
}

// 三条原覆盖行若没有一列能同时覆盖三行，并集必选至少两列。
// 以下阈值、索引截断和工作限额仅用于寻找被松弛点违反的并集。
fn separate(
    rows: &mut Vec<Row>,
    base: usize,
    x: &[f32],
    seen: &mut BTreeSet<Vec<usize>>,
    end: Option<Instant>,
) {
    let n = x.len();
    let mut posts = vec![Vec::new(); n];
    let mut tight = vec![false; base];
    for (i, r) in rows[..base].iter().enumerate() {
        if r.cols.iter().map(|&j| x[j]).sum::<f32>() < 1.25 {
            tight[i] = true;
            for &j in &r.cols {
                if x[j] > 0.001 {
                    posts[j].push(i);
                }
            }
        }
    }
    let mut neighbors = vec![Vec::new(); base];
    for i in 0..base {
        if i % 32 == 0 && expired(end) {
            return;
        }
        if !tight[i] {
            continue;
        }
        let a = &mut neighbors[i];
        for &j in &rows[i].cols {
            if posts[j].len() <= 128 && x[j] > 0.001 {
                a.extend(posts[j].iter().copied().filter(|&k| k > i));
            }
        }
        a.sort_unstable();
        a.dedup();
        a.truncate(128);
    }
    let mut marked = vec![usize::MAX; base];
    let mut count = vec![0u8; n];
    let mut union = Vec::new();
    let mut budget = 4_000_000usize;
    for i in 0..base {
        if i % 8 == 0 && expired(end) {
            return;
        }
        for &j in &neighbors[i] {
            marked[j] = i;
        }
        for &j in &neighbors[i] {
            for &k in &neighbors[j] {
                if marked[k] != i {
                    continue;
                }
                let mut activity = 0.;
                let mut common = false;
                for &t in &[i, j, k] {
                    budget = budget.saturating_sub(rows[t].cols.len());
                    for &col in &rows[t].cols {
                        if count[col] == 0 {
                            union.push(col);
                            activity += x[col];
                        }
                        count[col] += 1;
                        common |= count[col] == 3;
                    }
                }
                if !common && activity < 1.98 {
                    union.sort_unstable();
                    if seen.insert(union.clone()) {
                        rows.push(Row::new(union.clone(), 2));
                    }
                }
                for &col in &union {
                    count[col] = 0;
                }
                union.clear();
                if budget == 0 || seen.len() >= 1024 {
                    return;
                }
            }
        }
    }
}

pub(super) fn lower_bound(source: &[Vec<usize>], n: usize, end: Option<Instant>) -> usize {
    if source.is_empty() || expired(end) {
        return 0;
    }
    let base = source.len();
    let mut rows: Vec<_> = source.iter().map(|r| Row::new(r.clone(), 1)).collect();
    let mut x = vec![0.3f32; n];
    let mut extrapolated = x.clone();
    let mut cost = vec![0f32; n];
    let mut seen = BTreeSet::new();
    let mut best = 0;
    let mut work = 100_000_000usize;
    for stage in 0..12 {
        if expired(end) {
            break;
        }
        let mut lengths = vec![0usize; n];
        let mut nonzeros = 0;
        for r in &rows {
            nonzeros += r.cols.len();
            for &j in &r.cols {
                lengths[j] += 1;
            }
        }
        let step: Vec<_> = lengths
            .iter()
            .map(|&d| 1. / (0.3 * d.max(1) as f32))
            .collect();
        for it in 0..if stage == 0 { 160 } else { 200 } {
            if work < nonzeros || (it % 8 == 0 && expired(end)) {
                return best.max(certificate(&rows, n));
            }
            work -= nonzeros;
            cost.fill(-1.);
            for r in &mut rows {
                let mut gap = r.rhs as f32;
                for &j in &r.cols {
                    gap -= extrapolated[j];
                }
                r.price = (r.price + r.step * gap).max(0.);
                for &j in &r.cols {
                    cost[j] += r.price;
                }
            }
            for j in 0..n {
                let new = (x[j] + step[j] * cost[j]).clamp(0., 1.);
                extrapolated[j] = new + new - x[j];
                x[j] = new;
            }
        }
        best = best.max(certificate(&rows, n));
        if (stage < 4 || stage % 3 == 0) && seen.len() < 1024 {
            separate(&mut rows, base, &x, &mut seen, end);
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intersecting_pairs_require_two_and_boxes_correct_overpricing() {
        let src = vec![vec![0, 1], vec![1, 2], vec![0, 2]];
        let mut rows: Vec<_> = src.iter().map(|r| Row::new(r.clone(), 1)).collect();
        let mut seen = BTreeSet::new();
        separate(&mut rows, 3, &[0.5; 3], &mut seen, None);
        assert_eq!(rows.len(), 4);
        assert_eq!(rows[3].rhs, 2);
        assert_eq!(rows[3].cols, vec![0, 1, 2]);
        for mask in 0..8 {
            let activity =
                |r: &Row| r.cols.iter().filter(|&&j| mask & (1 << j) != 0).count() as i64;
            if rows[..3].iter().all(|r| activity(r) >= 1) {
                assert!(activity(&rows[3]) >= 2);
            }
        }
        for (r, y) in rows.iter_mut().zip([1.9, 0.25, 1.2, 0.6]) {
            r.price = y;
        }
        assert!(certificate(&rows, 3) <= 2);
        assert_eq!(lower_bound(&src, 3, None), 2);
        // 三行共享列 0，取它即可全覆盖，不能强加至少取两列的割。
        let star = vec![vec![0, 1], vec![0, 2], vec![0, 3]];
        assert_eq!(lower_bound(&star, 4, None), 1);
    }
}
