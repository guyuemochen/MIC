//! 大型等价费用集合覆盖的可行解。单选行先固定，覆盖范围被同价列包含的列可以替换。
//! 贪心解建好后，每次少选一列并交换集合，反复漏掉的行逐渐加权。只回传原模型可行点。

use crate::model::{Problem, RelOp, Sense, VarType};
use std::time::Instant;

#[path = "cover_repair/bound.rs"]
mod bound;

fn expired(end: Option<Instant>) -> bool {
    end.is_some_and(|t| Instant::now() >= t)
}

fn incidence(rows: &[Vec<usize>], n: usize) -> Vec<Vec<usize>> {
    let mut cols = vec![Vec::new(); n];
    for (i, row) in rows.iter().enumerate() {
        for &j in row {
            cols[j].push(i);
        }
    }
    cols
}

struct Cover {
    rows: Vec<Vec<usize>>,
    fixed: Vec<bool>,
}

impl Cover {
    fn detect(p: &Problem) -> Option<Self> {
        let n = p.vars.len();
        if !(256..=4096).contains(&n)
            || !(1024..=16384).contains(&p.constraints.len())
            || p.sense != Sense::Minimize
            || p.vars
                .iter()
                .any(|v| v.var_type != VarType::Binary || v.lb != 0. || v.ub != 1.)
            || p.obj.terms.len() != n
            || p.obj
                .terms
                .iter()
                .enumerate()
                .any(|(k, &(j, c))| j != k || c != 1.)
        {
            return None;
        }
        let mut rows = Vec::with_capacity(p.constraints.len());
        let mut remaining = 200_000usize;
        for r in &p.constraints {
            if !matches!((r.op, r.rhs), (RelOp::Ge, 1.) | (RelOp::Le, -1.))
                || r.expr.terms.is_empty()
                || r.expr.terms.len() > 128
                || r.expr.terms.iter().any(|&(_, a)| a != r.rhs)
            {
                return None;
            }
            remaining = remaining.checked_sub(r.expr.terms.len())?;
            rows.push(r.expr.terms.iter().map(|&(j, _)| j).collect());
        }
        Some(Self {
            rows,
            fixed: vec![false; n],
        })
    }

    fn reduce(&mut self, end: Option<Instant>) {
        for _ in 0..12 {
            if expired(end) {
                return;
            }
            for r in &self.rows {
                if r.len() == 1 {
                    self.fixed[r[0]] = true;
                }
            }
            self.rows.retain(|r| !r.iter().any(|&j| self.fixed[j]));
            let cols = incidence(&self.rows, self.fixed.len());
            let mut drop = vec![false; cols.len()];
            for (j, col) in cols.iter().enumerate() {
                if j % 16 == 0 && expired(end) {
                    return;
                }
                let Some(&i) = col.iter().min_by_key(|&&i| self.rows[i].len()) else {
                    continue;
                };
                drop[j] = self.rows[i].iter().any(|&k| {
                    (cols[k].len() > col.len() || (cols[k].len() == col.len() && k < j))
                        && col.iter().all(|v| cols[k].binary_search(v).is_ok())
                });
            }
            if !drop.contains(&true) {
                return;
            }
            for r in &mut self.rows {
                r.retain(|&j| !drop[j]);
            }
        }
    }
}

struct Search {
    rows: Vec<Vec<usize>>,
    cols: Vec<Vec<usize>>,
    count: Vec<usize>,
    // 选中列号的异或。在覆盖数恰为 1 时给出唯一覆盖列。
    owners: Vec<usize>,
    price: Vec<i64>,
    gain: Vec<i64>,
    loss: Vec<i64>,
    chosen: Vec<usize>,
    position: Vec<usize>,
    missed: Vec<usize>,
    missed_pos: Vec<usize>,
}

impl Search {
    fn new(rows: Vec<Vec<usize>>, n: usize) -> Self {
        let m = rows.len();
        let cols = incidence(&rows, n);
        let gain = cols.iter().map(|c| c.len() as i64).collect();
        Self {
            rows,
            cols,
            count: vec![0; m],
            owners: vec![0; m],
            price: vec![1; m],
            gain,
            loss: vec![0; n],
            chosen: vec![],
            position: vec![usize::MAX; n],
            missed: (0..m).collect(),
            missed_pos: (0..m).collect(),
        }
    }

    fn select(&mut self, j: usize) {
        self.position[j] = self.chosen.len();
        self.chosen.push(j);
        for &i in &self.cols[j] {
            if self.count[i] == 0 {
                let pos = self.missed_pos[i];
                self.missed.swap_remove(pos);
                if pos < self.missed.len() {
                    self.missed_pos[self.missed[pos]] = pos;
                }
                for &k in &self.rows[i] {
                    self.gain[k] -= self.price[i];
                }
                self.loss[j] += self.price[i];
            } else if self.count[i] == 1 {
                self.loss[self.owners[i]] -= self.price[i];
            }
            self.owners[i] ^= j;
            self.count[i] += 1;
        }
    }

    fn remove(&mut self, j: usize) {
        let pos = self.position[j];
        self.chosen.swap_remove(pos);
        if pos < self.chosen.len() {
            self.position[self.chosen[pos]] = pos;
        }
        self.position[j] = usize::MAX;
        self.loss[j] = 0;
        for &i in &self.cols[j] {
            self.owners[i] ^= j;
            self.count[i] -= 1;
            if self.count[i] == 0 {
                self.missed_pos[i] = self.missed.len();
                self.missed.push(i);
                for &k in &self.rows[i] {
                    self.gain[k] += self.price[i];
                }
            } else if self.count[i] == 1 {
                self.loss[self.owners[i]] += self.price[i];
            }
        }
    }

    fn run(
        &mut self,
        p: &Problem,
        fixed: &[bool],
        end: Option<Instant>,
        accept: &mut dyn FnMut(&[f64]),
    ) -> Option<Vec<f64>> {
        let n = fixed.len();
        while !self.missed.is_empty() {
            if expired(end) {
                return None;
            }
            let j = (0..n).min_by_key(|&j| (-self.gain[j], j))?;
            if self.gain[j] == 0 {
                return None;
            }
            self.select(j);
        }
        for j in (0..n).rev() {
            if self.position[j] != usize::MAX && self.loss[j] == 0 {
                self.remove(j);
            }
        }
        let mut best_len = usize::MAX;
        let mut best = None;
        let mut age = vec![0usize; n];
        let mut last = usize::MAX;
        let mut rng = 1u64;
        // 按已选列扫描和受影响的关联项计工作量，大模型同样遵守启发式预算。
        let mut work = 150_000_000usize;
        for step in 1..=350_000 {
            if expired(end) {
                break;
            }
            let full = self.missed.is_empty();
            if full && self.chosen.len() < best_len {
                let x: Vec<_> = (0..n)
                    .map(|j| f64::from(fixed[j] || self.position[j] != usize::MAX))
                    .collect();
                if p.constraints.iter().any(|r| r.expr.eval(&x) * r.rhs < 1.) {
                    break;
                }
                accept(&x);
                best = Some(x);
                best_len = self.chosen.len();
                if best_len <= 1 {
                    break;
                }
            }
            if work == 0 {
                break;
            }
            for _ in 0..if full { 2 } else { 1 } {
                let Some(j) = self
                    .chosen
                    .iter()
                    .copied()
                    .filter(|&j| j != last)
                    .min_by_key(|&j| (self.loss[j], age[j], j))
                else {
                    break;
                };
                work = work.saturating_sub(self.chosen.len() + self.cols[j].len());
                self.remove(j);
                age[j] = step;
                last = j;
            }
            if self.missed.is_empty() {
                continue;
            }
            // 固定序列抽取漏行，避免反复修同一个局部矛盾。
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let i = self.missed[(rng >> 32) as usize % self.missed.len()];
            let j = self.rows[i]
                .iter()
                .copied()
                .filter(|&j| j != last)
                .min_by_key(|&j| (-self.gain[j], age[j], j))
                .or_else(|| self.rows[i].first().copied())?;
            self.select(j);
            age[j] = step;
            last = j;
            work = work.saturating_sub(self.cols[j].len());
            for &i in &self.missed {
                self.price[i] += 1;
                work = work.saturating_sub(self.rows[i].len());
                for &k in &self.rows[i] {
                    self.gain[k] += 1;
                }
            }
        }
        best
    }
}

pub(crate) fn search(
    p: &Problem,
    deadline: Option<Instant>,
    accept: &mut dyn FnMut(&[f64]),
) -> Option<(Vec<f64>, Option<f64>)> {
    if expired(deadline) {
        return None;
    }
    let mut model = Cover::detect(p)?;
    model.reduce(deadline);
    let mut state = Search::new(model.rows, p.vars.len());
    let x = state.run(p, &model.fixed, deadline, accept)?;
    // 固定列只来自单选行与等价替换，贪心和换列的选择不会进入下界。
    let fixed = model.fixed.iter().filter(|&&f| f).count();
    let count = fixed + bound::lower_bound(&state.rows, p.vars.len(), deadline);
    let total = p.obj_offset + count as f64;
    // 非整数偏移的加法向外舍入；整数安全范围内保持精确的目标尺度。
    let lb = if p.obj_offset.fract() == 0. && p.obj_offset.abs() <= (1u64 << 52) as f64 {
        total
    } else {
        total.next_down()
    };
    Some((x, lb.is_finite().then_some(lb)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check_scores(s: &Search) {
        let mut gain = vec![0; s.cols.len()];
        let mut loss = gain.clone();
        let mut missing = Vec::new();
        for (i, r) in s.rows.iter().enumerate() {
            let chosen: Vec<_> = r
                .iter()
                .copied()
                .filter(|&j| s.position[j] != usize::MAX)
                .collect();
            assert_eq!(s.count[i], chosen.len());
            assert_eq!(s.owners[i], chosen.iter().fold(0, |v, j| v ^ j));
            match chosen.as_slice() {
                [] => {
                    missing.push(i);
                    for &j in r {
                        gain[j] += s.price[i];
                    }
                }
                &[j] => loss[j] += s.price[i],
                _ => {}
            }
        }
        assert_eq!(s.gain, gain);
        assert_eq!(s.loss, loss);
        let mut recorded = s.missed.clone();
        recorded.sort_unstable();
        assert_eq!(recorded, missing);
        for (k, &j) in s.chosen.iter().enumerate() {
            assert_eq!(s.position[j], k);
        }
        for (k, &i) in s.missed.iter().enumerate() {
            assert_eq!(s.missed_pos[i], k);
        }
    }

    #[test]
    fn gain_loss_and_uncovered_rows_match_full_recount() {
        let mut s = Search::new(vec![vec![0, 2, 3], vec![0, 1], vec![1, 2], vec![2, 3]], 4);
        for (i, &w) in [3, 8, 17, 9].iter().enumerate() {
            s.price[i] = w;
            for &j in &s.rows[i] {
                s.gain[j] += w - 1;
            }
        }
        for j in [0, 3, 2, 0, 1, 3, 2, 1] {
            if s.position[j] == usize::MAX {
                s.select(j);
            } else {
                s.remove(j);
            }
            check_scores(&s);
        }
    }

    #[test]
    fn mandatory_and_dominated_columns_restore_original_cover() {
        let mut p = Problem::new("覆盖", Sense::Minimize);
        for j in 0..256 {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        p.set_objective((0..256).map(|j| (j, 1.)).collect(), -0.5)
            .unwrap();
        let rows = [
            vec![(0, 1.)],
            vec![(0, 1.), (4, 1.)],
            vec![(1, 1.), (2, 1.)],
            vec![(1, 1.), (3, 1.)],
        ];
        for i in 0..1024 {
            p.add_constraint(None, rows[i % 4].clone(), RelOp::Ge, 1.)
                .unwrap();
        }
        let mut events = Vec::new();
        let (result, lower) = search(&p, None, &mut |x| events.push(x.to_vec())).unwrap();
        assert!(lower.is_some_and(|b| b <= 1.5));
        assert_eq!(events, vec![result.clone()]);
        assert_eq!(result[0..4], [1., 1., 0., 0.]);
        assert!(p.constraints.iter().all(|r| r.expr.eval(&result) >= r.rhs));
        assert_eq!(p.obj_offset + p.obj.eval(&result), 1.5);
        assert!(search(&p, Some(Instant::now()), &mut |_| panic!(
            "到时不应发布新解"
        ))
        .is_none());
        p.obj.terms[1].1 = 0.75;
        assert!(Cover::detect(&p).is_none());
        p.obj.terms[1].1 = 1.;
        p.constraints[0].rhs = 2.;
        assert!(Cover::detect(&p).is_none());
    }
}
