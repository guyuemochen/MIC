//! 按每单位费用新增的覆盖行数补齐可行解，删去冗余列，再局部拆开重建。
//! 扰动只影响尝试选列的顺序。每次完整覆盖后按原整数费用比较。

use super::Cover;
use std::time::Instant;

fn next_random(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

impl Cover {
    fn complete(&self, chosen: &mut [bool], price: &[i64], expensive_first: &[usize]) {
        let mut coverage = vec![0u16; self.rows.len()];
        for (j, &used) in chosen.iter().enumerate() {
            if used {
                for &i in &self.cols[j] {
                    coverage[i] += 1;
                }
            }
        }
        let mut missing = 0;
        let mut gain = vec![0i64; chosen.len()];
        for (i, row) in self.rows.iter().enumerate() {
            if coverage[i] == 0 {
                missing += 1;
                for &j in row {
                    gain[j] += 1;
                }
            }
        }
        while missing > 0 {
            let mut best = 0;
            for j in 1..gain.len() {
                if gain[j] * price[best] > gain[best] * price[j] {
                    best = j;
                }
            }
            chosen[best] = true;
            for &i in &self.cols[best] {
                if coverage[i] == 0 {
                    missing -= 1;
                    for &j in &self.rows[i] {
                        gain[j] -= 1;
                    }
                }
                coverage[i] += 1;
            }
        }
        for &j in expensive_first {
            if chosen[j] && self.cols[j].iter().all(|&i| coverage[i] > 1) {
                chosen[j] = false;
                for &i in &self.cols[j] {
                    coverage[i] -= 1;
                }
            }
        }
    }

    pub(super) fn seed(
        &self,
        deadline: Option<Instant>,
        update: &mut dyn FnMut(&[bool]),
    ) -> Vec<bool> {
        let mut order: Vec<_> = (0..self.costs.len()).collect();
        order.sort_by_key(|&j| (std::cmp::Reverse(self.costs[j]), j));
        let mut best = vec![false; self.costs.len()];
        self.complete(&mut best, &self.costs, &order);
        update(&best);
        let mut value = self.value(&best);
        let mut state = 1u64;
        let mut price = self.costs.clone();
        let mut last_improvement = 0;
        for round in 0..512 {
            if round - last_improvement >= 192 || super::expired(deadline) {
                break;
            }
            let mut chosen = best.clone();
            let keep = (round % 4) as u64;
            for j in 0..chosen.len() {
                let r = next_random(&mut state);
                chosen[j] &= r % 5 < keep;
                price[j] = self.costs[j] * (64 + ((r >> 8) % 64) as i64);
            }
            self.complete(&mut chosen, &price, &order);
            let cost = self.value(&chosen);
            if cost < value {
                last_improvement = round;
                value = cost;
                best = chosen;
                update(&best);
            }
        }
        best
    }
}
