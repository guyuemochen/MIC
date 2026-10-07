//! 确定性伪随机数生成器（契约逐位钉死，实现必须逐字一致，不得改动任何常数或消耗顺序）。
//!
//! 算法：SplitMix64 种子扩展 + xoshiro256\*\* 主循环；[`Rng::gen_range`] 用标准 Lemire
//! 拒绝采样；[`Rng::gen_ratio`] 用 53 位均匀判定。决策路径禁止依赖 HashMap 迭代序。

/// xoshiro256\*\* 状态（由 [`Rng::new`] 以 SplitMix64 扩展种子产生）。
#[derive(Debug, Clone)]
pub struct Rng {
    s: [u64; 4],
}

/// SplitMix64 单步（种子扩展专用，常数与混合序为 normative）：
/// `*x += 0x9E3779B97F4A7C15`；`z = *x`；`z = (z ^ (z>>30)) * 0xBF58476D1CE4E5B9`；
/// `z = (z ^ (z>>27)) * 0x94D049BB133111EB`；返回 `z ^ (z>>31)`，并经 `x` 传出更新后的状态。
fn splitmix(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9E3779B97F4A7C15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

impl Rng {
    /// 以 `seed` 构造：SplitMix64 连续产出 4 个状态字。
    ///
    /// splitmix(x)：`*x = x.wrapping_add(0x9E3779B97F4A7C15); z = *x;`
    /// `z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);`
    /// `z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31)`。
    pub fn new(seed: u64) -> Rng {
        let mut x = seed;
        let mut s = [0u64; 4];
        for slot in &mut s {
            *slot = splitmix(&mut x);
        }
        Rng { s }
    }

    /// xoshiro256\*\* 主循环：
    /// `r = rotl(s[1].wrapping_mul(5), 7).wrapping_mul(9)`；`t = s[1] << 17`；
    /// `s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = s[3].rotate_left(45)`；返回 r。
    pub fn next_u64(&mut self) -> u64 {
        let r = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        r
    }

    /// 均匀采样 `[low, high)`。前置条件 `low < high`（违例 panic，内部调用点恒满足）。
    ///
    /// 标准 Lemire 拒绝采样（接受方向与取位钉死）：`len = high - low`；
    /// `m = (next_u64() as u128) * (len as u128)`；`l = m 低 64 位`；
    /// `t = len.wrapping_neg() % len`（即 2^64 mod len）；`l >= t` 接受、结果
    /// `low + (m >> 64) as usize`；`l < t` 整个乘法重抽。
    pub fn gen_range(&mut self, low: usize, high: usize) -> usize {
        assert!(low < high, "gen_range: low ({low}) must be less than high ({high})");
        let len = (high - low) as u64;
        let t = len.wrapping_neg() % len;
        loop {
            let m = u128::from(self.next_u64()) * u128::from(len);
            let l = m as u64;
            if l >= t {
                return low + (m >> 64) as usize;
            }
        }
    }

    /// 以概率 p 返回 true：`p <= 0.0` 恒 false；`p >= 1.0` 恒 true；否则
    /// `((next_u64() >> 11) as f64) < p * 9007199254740992.0`。
    /// [`crate::moves::MoveKind::Mixed`] 下 true ⇒ 生成 Change、false ⇒ 生成 Swap。
    pub fn gen_ratio(&mut self, p: f64) -> bool {
        if p <= 0.0 {
            return false;
        }
        if p >= 1.0 {
            return true;
        }
        ((self.next_u64() >> 11) as f64) < p * 9007199254740992.0
    }

    /// Fisher-Yates（向后）均匀洗牌。
    pub fn shuffle<T>(&mut self, xs: &mut [T]) {
        for i in (1..xs.len()).rev() {
            let j = self.gen_range(0, i + 1);
            xs.swap(i, j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Rng;

    /// 参考流：由本 crate 规范公式（SplitMix64 播种 + xoshiro256**）的独立 Python
    /// 转写离线计算（2026-10-07，Python 3.13），用于钉死逐位一致性。
    const SEED0_FIRST6: [u64; 6] = [
        0x99ec5f36cb75f2b4,
        0xbf6e1f784956452a,
        0x1a5f849d4933e6e0,
        0x6aa594f1262d2d2c,
        0xbba5ad4a1f842e59,
        0xffef8375d9ebcaca,
    ];
    const SEED42_FIRST6: [u64; 6] = [
        0x15780b2e0c2ec716,
        0x6104d9866d113a7e,
        0xae17533239e499a1,
        0xecb8ad4703b360a1,
        0xfde6dc7fe2ec5e64,
        0xc50da53101795238,
    ];
    /// seed 7 对 [0,1,...,9] 做规范 Fisher-Yates（向后）的结果（同上离线计算）。
    const SHUFFLE_SEED7_0TO9: [usize; 10] = [1, 8, 3, 0, 4, 5, 9, 6, 2, 7];

    #[test]
    fn same_seed_same_stream_bit_for_bit() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..1000 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
        // 克隆亦逐位一致
        let mut c = a.clone();
        assert_eq!(a.next_u64(), c.next_u64());
    }

    #[test]
    fn different_seeds_differ() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(43);
        let mut equal = 0usize;
        for _ in 0..64 {
            if a.next_u64() == b.next_u64() {
                equal += 1;
            }
        }
        assert!(equal < 8, "different seeds should give different streams (equal={equal})");
    }

    #[test]
    fn matches_offline_reference_vectors() {
        let mut r = Rng::new(0);
        for want in SEED0_FIRST6 {
            assert_eq!(r.next_u64(), want);
        }
        let mut r = Rng::new(42);
        for want in SEED42_FIRST6 {
            assert_eq!(r.next_u64(), want);
        }
    }

    #[test]
    fn gen_range_stays_in_domain() {
        let mut r = Rng::new(2026);
        for &(low, high) in &[(0usize, 1usize), (0, 2), (3, 17), (100, 1000), (0, 1 << 20)] {
            for _ in 0..10_000 {
                let v = r.gen_range(low, high);
                assert!(v >= low && v < high, "gen_range({low},{high}) = {v} out of range");
            }
        }
        // 退化为单值区间：恒返回 low（且 len=1 时 t=0，恒接受，不耗额外重抽路径出错）
        for _ in 0..100 {
            assert_eq!(r.gen_range(5, 6), 5);
        }
    }

    #[test]
    fn gen_range_invalid_interval_panics() {
        let mut r = Rng::new(1);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| r.gen_range(3, 3)));
        assert!(result.is_err(), "gen_range with low == high must panic");
    }

    #[test]
    fn gen_ratio_degenerate_probabilities() {
        let mut r = Rng::new(42);
        for _ in 0..8 {
            assert!(!r.gen_ratio(0.0));
        }
        for _ in 0..8 {
            assert!(!r.gen_ratio(-0.5));
        }
        for _ in 0..8 {
            assert!(r.gen_ratio(1.0));
        }
        for _ in 0..8 {
            assert!(r.gen_ratio(1.5));
        }
    }

    #[test]
    fn gen_ratio_half_is_about_half() {
        let mut r = Rng::new(42);
        let mut trues = 0usize;
        let n = 100_000;
        for _ in 0..n {
            if r.gen_ratio(0.5) {
                trues += 1;
            }
        }
        let ratio = trues as f64 / n as f64;
        assert!(ratio > 0.49 && ratio < 0.51, "gen_ratio(0.5) ratio {ratio} far from 0.5");
    }

    #[test]
    fn shuffle_is_permutation_and_matches_reference() {
        let mut r = Rng::new(7);
        let mut xs: [usize; 10] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        r.shuffle(&mut xs);
        assert_eq!(xs, SHUFFLE_SEED7_0TO9);

        let mut r = Rng::new(99);
        let mut ys: Vec<usize> = (0..64).collect();
        let orig = ys.clone();
        r.shuffle(&mut ys);
        assert_eq!(ys.len(), orig.len());
        let mut sorted = ys.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, orig, "shuffle must be a permutation");
    }
}
