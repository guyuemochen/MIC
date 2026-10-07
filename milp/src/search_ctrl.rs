//! 搜索控制层（纯逻辑：零锁、零线程、全部可单元测试）。
//!
//! 交付物：
//! - 进度快照 [`ProgressInfo`] 与其唯一格式化口径（进度行 / gap 判定共用，
//!   数值规则与 CLI objective 行一致）；
//! - gap 百分比 [`gap_percent`]（一律 clamp 到 ≥0，见 N1）与对偶界聚合
//!   [`aggregate_bound`]（空迭代 ⇒ None，无哨兵毒化）；
//! - 自适应启发式控制器 [`HeuristicController`]（『无 incumbent ⇄ stall ⇄
//!   近期命中率高』三条件 + 节流 + 时间预算保护）。
//!
//! 纪律：本模块不含锁、线程与真实计时——时间一律经 [`SearchSnapshot`]
//! 注入，测试可伪造；多线程编排见 [`crate::parallel_bb`]。

use std::cell::Cell;
use std::collections::VecDeque;
use std::sync::Arc;

use crate::model::Sense;

/// 相对 gap 的分母保护阈值：|incumbent| 低于该值时相对口径无意义。
const TINY_DENOM: f64 = 1e-10;

/// 进度快照（进度行与 gap 判定唯一数据源）。
#[derive(Debug, Clone, Copy)]
pub struct ProgressInfo {
    /// 目标方向。
    pub sense: Sense,
    /// 距求解开始的墙钟秒数。
    pub elapsed_secs: f64,
    /// 已处理节点数（根计 1）。
    pub nodes: u64,
    /// 当前 incumbent（原始 sense + obj_offset 口径；仅有限值，否则 None）。
    pub incumbent: Option<f64>,
    /// 全局对偶界（保守聚合；仅有限值，否则 None）。
    pub bound: Option<f64>,
}

/// gap 百分比（恒非负）。
///
/// Max 方向 `raw = (bound − inc)/|inc| × 100`、Min 镜像，返回
/// `Some(raw.max(0.0))`——db ≤ inc（incumbent 跳跃越过陈旧 est，或
/// est = 39.999999999999996 < inc = 40.0 的浮点噪声）时 gap = 0 而非负值
/// （实测 `{:.2}` 对 −2.5e-11 输出 `-0.00`、对 −7.69 输出 `-7.69`，负值会
/// 击穿 `gap=\d+\.\d+%` 正则与『gap_reached 行永不矛盾』契约）。
/// db ≤ inc 触发 0% gap 是诚实的：全部剩余子树 est ≤ inc ⇒ 最优 ≤
/// max(inc, db) = inc ⇒ incumbent 在 est 有效性内已达最优。
///
/// 返回 `None`（表示 inf）当且仅当：任一输入非有限，或 |inc| < 1e-10 且
/// 方向差不 ≈ 0（分母淹没，相对 gap 无意义）。
pub fn gap_percent(sense: Sense, incumbent: f64, bound: f64) -> Option<f64> {
    if !incumbent.is_finite() || !bound.is_finite() {
        return None;
    }
    let diff = match sense {
        Sense::Maximize => bound - incumbent,
        Sense::Minimize => incumbent - bound,
    };
    if incumbent.abs() < TINY_DENOM {
        // 分母 ≈ 0：差也 ≈ 0（含浮点噪声）⇒ 0%；否则相对口径无意义 ⇒ None。
        if diff.abs() <= TINY_DENOM {
            Some(0.0)
        } else {
            None
        }
    } else {
        Some((diff / incumbent.abs() * 100.0).max(0.0))
    }
}

/// 对偶界聚合：Max 取 est 最大值、Min 取最小值；空迭代 ⇒ None。
///
/// 调用方负责先 filter 掉 None/非有限槽位（空闲 worker 的 inflight 槽不
/// 参与聚合——无哨兵毒化）。
pub fn aggregate_bound(sense: Sense, ests: impl Iterator<Item = f64>) -> Option<f64> {
    ests.fold(None, |acc: Option<f64>, est| match (acc, sense) {
        (None, _) => Some(est),
        (Some(a), Sense::Maximize) => Some(a.max(est)),
        (Some(a), Sense::Minimize) => Some(a.min(est)),
    })
}

/// 与 CLI objective 行同规则的数值格式化：
/// `0 < |v| < 1e-6` 或 `|v| ≥ 1e15` 用 `{:e}`，否则 `{:.6}`。
///
/// 仅有限值属契约输入；非有限输入按 Display 输出（`NaN`/`inf`/`-inf`，
/// 防御性兜底，不 panic）。
pub fn format_number(v: f64) -> String {
    if !v.is_finite() {
        return format!("{v}");
    }
    let magnitude = v.abs();
    if (0.0 < magnitude && magnitude < 1e-6) || magnitude >= 1e15 {
        format!("{v:e}")
    } else {
        format!("{v:.6}")
    }
}

/// 进度行：
/// `progress t={:.2} nodes={} incumbent={}|none bound={}|none gap={:.2}%|inf`。
///
/// incumbent/bound 缺失（或非有限）打印 `none`；gap 按 [`gap_percent`]
/// （clamp 后恒非负）打印 `{:.2}%`，无法计算（缺任一侧或 inf 语义）打印
/// `inf`。
pub fn format_progress_line(info: &ProgressInfo) -> String {
    let incumbent = match info.incumbent {
        Some(v) if v.is_finite() => format_number(v),
        _ => "none".to_string(),
    };
    let bound = match info.bound {
        Some(v) if v.is_finite() => format_number(v),
        _ => "none".to_string(),
    };
    let gap = match gap_percent(
        info.sense,
        info.incumbent.unwrap_or(f64::NAN),
        info.bound.unwrap_or(f64::NAN),
    ) {
        Some(g) => format!("{g:.2}%"),
        None => "inf".to_string(),
    };
    format!(
        "progress t={:.2} nodes={} incumbent={} bound={} gap={}",
        info.elapsed_secs, info.nodes, incumbent, bound, gap
    )
}

/// 进度回调（Arc 共享多线程；回调内不得获取求解器任何锁、应快速返回）。
pub type ProgressCallback = Arc<dyn Fn(&ProgressInfo) + Send + Sync>;

/// 自适应控制器参数（默认 32/50/8/0.25/0.8）。
#[derive(Debug, Clone, Copy)]
pub struct HeuristicConfig {
    /// 两次 dive 最小节点间隔（节流）。
    pub min_interval_nodes: u64,
    /// 连续无改进阈值。
    pub stall_nodes: u64,
    /// 命中率统计窗口。
    pub hit_window: usize,
    /// 窗口命中率 ≥ 该值继续投资。
    pub keep_hit_rate: f64,
    /// elapsed > frac × budget 且有 incumbent ⇒ 停 dive（把尾段留给最优性
    /// 证明）。
    pub time_budget_frac: f64,
}

impl Default for HeuristicConfig {
    /// 默认值：32/50/8/0.25/0.8。
    fn default() -> Self {
        HeuristicConfig {
            min_interval_nodes: 32,
            stall_nodes: 50,
            hit_window: 8,
            keep_hit_rate: 0.25,
            time_budget_frac: 0.8,
        }
    }
}

/// 控制器决策。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HeuristicDecision {
    /// 运行一次 dive。
    Run,
    /// 跳过。
    Skip,
}

/// 决策输入（纯数据可伪造 ⇒ 测试不依赖真实计时）。
#[derive(Debug, Clone, Copy)]
pub struct SearchSnapshot {
    /// 当前已处理节点数。
    pub nodes: u64,
    /// 当前 incumbent（仅有限值语义）。
    pub incumbent: Option<f64>,
    /// 距上次改进的节点数。
    pub nodes_since_improvement: u64,
    /// 时间预算（秒）；None = 不限时。
    pub time_budget_secs: Option<f64>,
    /// 已耗时（秒）。
    pub elapsed_secs: f64,
}

/// 自适应启发式控制器。
///
/// auto 策略理由：dive 成本 ≈ dive_max_fixes 次 LP 重解；仅在『无 incumbent
/// ⇄ stall ⇄ 近期命中率高（命中律）』时运行，节流防过频；时间预算保护把
/// 尾段留给最优性证明。正确性：dive 只供应 incumbent，不改剪枝/分支语义
/// ⇒ auto 仍精确最优；默认阈值高于原 7 算例节点数（≤ 7）⇒ 原算例
/// auto ≡ on。
///
/// 节流锚点说明：`decide` 可观察行为纯（不改命中率/节流基准），但会把快照
/// 节点缓存到内部 `Cell`；`record` 以该缓存推进 `last_run_node`——即
/// 『授权某次 dive 的那次 decide 所在节点』成为下一次节流的基准，实现
/// 『两次 dive 最小节点间隔』。多线程下锁内串行调用，缓存可能被并发
/// decide 推进到更晚的节点（只会把下一次 dive 推得更迟，保守方向）。
pub struct HeuristicController {
    /// 参数。
    cfg: HeuristicConfig,
    /// 上次（根/dive）运行所在节点——节流基准。
    last_run_node: u64,
    /// 命中窗口（尾部 hit_window 条）。
    window: VecDeque<bool>,
    /// decide 最近一次看到的快照节点（record 的节流推进锚点）。
    last_seen_node: Cell<u64>,
}

impl HeuristicController {
    /// 以给定参数构建。
    pub fn new(cfg: HeuristicConfig) -> Self {
        HeuristicController {
            cfg,
            last_run_node: 0,
            window: VecDeque::new(),
            last_seen_node: Cell::new(0),
        }
    }

    /// 根阶段完成后调用：`last_run_node = root_nodes`（根启发式刚运行过）。
    pub fn notify_root_done(&mut self, root_nodes: u64) {
        self.last_run_node = root_nodes;
        self.last_seen_node.set(root_nodes);
    }

    /// 决策（可观察行为不改状态）：
    /// 1. 节流未到（距上次运行 < min_interval_nodes）⇒ Skip；
    /// 2. 时间预算保护（有预算且有 incumbent 且 elapsed > frac×budget）⇒ Skip；
    /// 3. 无 incumbent，或 stall 达标（nodes_since_improvement ≥ stall_nodes），
    ///    或窗口命中率 ≥ keep_hit_rate ⇒ Run；否则 Skip。
    pub fn decide(&self, snap: &SearchSnapshot) -> HeuristicDecision {
        self.last_seen_node.set(snap.nodes);
        if snap.nodes.saturating_sub(self.last_run_node) < self.cfg.min_interval_nodes {
            return HeuristicDecision::Skip;
        }
        if let Some(budget) = snap.time_budget_secs {
            if snap.incumbent.is_some() && snap.elapsed_secs > self.cfg.time_budget_frac * budget {
                return HeuristicDecision::Skip;
            }
        }
        if snap.incumbent.is_none() {
            return HeuristicDecision::Run;
        }
        if snap.nodes_since_improvement >= self.cfg.stall_nodes {
            return HeuristicDecision::Run;
        }
        if let Some(rate) = self.hit_rate() {
            if rate >= self.cfg.keep_hit_rate {
                return HeuristicDecision::Run;
            }
        }
        HeuristicDecision::Skip
    }

    /// 记录一次 dive 结果（`hit` = 是否产生新全局 incumbent）。
    pub fn record(&mut self, hit: bool) {
        self.last_run_node = self.last_seen_node.get();
        self.window.push_back(hit);
        while self.window.len() > self.cfg.hit_window {
            self.window.pop_front();
        }
    }

    /// 窗口命中率（按已有个数计；无记录 ⇒ None）。
    pub fn hit_rate(&self) -> Option<f64> {
        if self.window.is_empty() {
            return None;
        }
        let hits = self.window.iter().filter(|h| **h).count();
        Some(hits as f64 / self.window.len() as f64)
    }
}

// ===========================================================================
// 策略级自动降级（auto 模式 + 时间预算；「继续证最优」vs「提前返回最好解」）
// ===========================================================================

/// 策略级自动降级参数（默认 0.5 / 0.25 / 0.01）。
///
/// 语义：仅当同时满足 ①`--heuristic-mode auto`（自适应启用）与 ②给了时间
/// 预算（用户表达了速度优先的需求）时才评估——无预算时 auto 仍保证证明
/// 最优（契约不变）。预算消耗达 `start_frac` 之后，若 incumbent 已连续
/// `plateau_frac`×预算 无改进、且期间对偶界的相对改进不足
/// `bound_progress_rel`（双停滞：解侧与证明侧都不再前进），判定继续证明
/// 最优的边际收益过低，提前终止并返回当前最好解（终态诚实标注
/// LimitReached + gap 证据，与 `--mipgap` 同径，另附 `auto_degraded` 标记）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DegradePolicy {
    /// 预算消耗比例达到该值后才允许降级（先给精确搜索公平机会）。
    pub start_frac: f64,
    /// incumbent 无改进窗口占预算的比例。
    pub plateau_frac: f64,
    /// 窗口内对偶界相对改进阈值（朝 incumbent 方向的改进低于该比例视为停滞）。
    pub bound_progress_rel: f64,
}

impl Default for DegradePolicy {
    /// 默认值：0.5 / 0.25 / 0.01。
    fn default() -> Self {
        DegradePolicy {
            start_frac: 0.5,
            plateau_frac: 0.25,
            bound_progress_rel: 0.01,
        }
    }
}

/// 降级判定输入（纯数据可伪造 ⇒ 测试不依赖真实计时）。
#[derive(Debug, Clone, Copy)]
pub struct DegradeSnapshot {
    /// 目标方向（决定对偶界的改进朝向）。
    pub sense: Sense,
    /// 已耗时（秒）。
    pub elapsed_secs: f64,
    /// 时间预算（秒）；None = 无预算（永不降级）。
    pub budget_secs: Option<f64>,
    /// 是否已有 incumbent（无解时绝不降级——没有可返回的东西）。
    pub has_incumbent: bool,
    /// 最近一次 incumbent 严格改进的时刻（秒）。
    pub last_improvement_secs: f64,
    /// 当前前沿聚合对偶界（仅有限值语义；None ⇒ 本次不判）。
    pub bound: Option<f64>,
}

/// 策略级自动降级判定器（有状态：停滞监测窗口的界基准快照）。
///
/// 纪律与 [`HeuristicController`] 一致：无锁无线程，时间与界均由快照注入。
pub struct AutoDegrade {
    cfg: DegradePolicy,
    /// 窗口起点的对偶界基准（None = 尚未采样——首个合格快照只建立基准）。
    window_db0: Option<f64>,
}

impl AutoDegrade {
    /// 以给定参数构建。
    pub fn new(cfg: DegradePolicy) -> Self {
        AutoDegrade {
            cfg,
            window_db0: None,
        }
    }

    /// 参数只读访问（调用方做资格预检时使用）。
    pub fn cfg(&self) -> DegradePolicy {
        self.cfg
    }

    /// 资格预检（无状态，可静态调用）：预算存在（有限且 > 0）且已消耗
    /// ≥ `start_frac`×预算。调用方先经此过滤再做 O(前沿) 的界聚合，避免
    /// 逐节点无谓开销。
    pub fn eligible(cfg: DegradePolicy, elapsed_secs: f64, budget_secs: Option<f64>) -> bool {
        match budget_secs {
            Some(b) if b.is_finite() && b > 0.0 => elapsed_secs >= cfg.start_frac * b,
            _ => false,
        }
    }

    /// 判定是否降级（调用方先用 [`AutoDegrade::eligible`] 过滤）：
    /// 1. 资格不满足（无预算 / 未到 start_frac）⇒ false；
    /// 2. 无 incumbent 或 incumbent 改进尚新（< plateau_frac×预算）⇒ false；
    /// 3. 界缺失或非有限 ⇒ false（防御；聚合层已滤非有限）；
    /// 4. 首个合格快照只建立窗口基准 ⇒ false；
    /// 5. 自基准以来界的相对改进（朝 incumbent 方向）< bound_progress_rel
    ///    ⇒ true（双停滞成立）；否则刷新基准继续监测 ⇒ false。
    pub fn decide(&mut self, snap: &DegradeSnapshot) -> bool {
        if !Self::eligible(self.cfg, snap.elapsed_secs, snap.budget_secs) {
            return false;
        }
        let budget = snap.budget_secs.unwrap_or(f64::NAN);
        if !snap.has_incumbent {
            return false;
        }
        if snap.elapsed_secs - snap.last_improvement_secs < self.cfg.plateau_frac * budget {
            return false;
        }
        let Some(db) = snap.bound.filter(|v| v.is_finite()) else {
            return false;
        };
        let Some(db0) = self.window_db0 else {
            self.window_db0 = Some(db);
            return false;
        };
        // Max 方向对偶界向下逼近 incumbent（改进 = db0 − db）；Min 镜像。
        let progress = match snap.sense {
            Sense::Maximize => db0 - db,
            Sense::Minimize => db - db0,
        };
        if progress / (1.0 + db0.abs()) < self.cfg.bound_progress_rel {
            return true;
        }
        // 界有实质进展 ⇒ 以当前界刷新基准（证明侧仍在前进，不降级）。
        self.window_db0 = Some(db);
        false
    }
}

// ===========================================================================
// 单元测试（全部不依赖计时：时间经 SearchSnapshot 伪造）
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// 默认参数逐字段符合契约（32/50/8/0.25/0.8）。
    #[test]
    fn heuristic_config_default_values() {
        let c = HeuristicConfig::default();
        assert_eq!(c.min_interval_nodes, 32);
        assert_eq!(c.stall_nodes, 50);
        assert_eq!(c.hit_window, 8);
        assert_eq!(c.keep_hit_rate, 0.25);
        assert_eq!(c.time_budget_frac, 0.8);
    }

    /// 无记录 ⇒ hit_rate None；有记录按个数计。
    #[test]
    fn hit_rate_none_then_ratio() {
        let mut c = HeuristicController::new(HeuristicConfig::default());
        assert!(c.hit_rate().is_none(), "无记录 ⇒ None");
        c.record(true);
        c.record(false);
        c.record(true);
        assert!((c.hit_rate().unwrap() - 2.0 / 3.0).abs() < 1e-12);
    }

    /// 窗口容量：第 9 条记录挤出最旧的一条。
    #[test]
    fn window_capacity_evicts_oldest() {
        let mut c = HeuristicController::new(HeuristicConfig::default());
        c.record(true); // 将被挤出
        for _ in 0..8 {
            c.record(false);
        }
        assert_eq!(c.hit_rate(), Some(0.0), "8 条 miss 全部在窗内");
    }

    /// 节流 Skip：距根 dive 的节点间隔未到 min_interval_nodes。
    #[test]
    fn controller_throttle_skips() {
        let mut c = HeuristicController::new(HeuristicConfig::default());
        c.notify_root_done(1);
        let snap = SearchSnapshot {
            nodes: 32, // 32 − 1 = 31 < 32
            incumbent: None,
            nodes_since_improvement: 100,
            time_budget_secs: None,
            elapsed_secs: 0.0,
        };
        assert_eq!(c.decide(&snap), HeuristicDecision::Skip);
        // 恰达间隔（差 == min_interval_nodes）⇒ 不再节流。
        let snap = SearchSnapshot {
            nodes: 33,
            incumbent: None,
            nodes_since_improvement: 0,
            time_budget_secs: None,
            elapsed_secs: 0.0,
        };
        assert_eq!(c.decide(&snap), HeuristicDecision::Run);
    }

    /// 无 incumbent ⇒ Run（节流已过、无预算）。
    #[test]
    fn controller_no_incumbent_runs() {
        let mut c = HeuristicController::new(HeuristicConfig::default());
        c.notify_root_done(1);
        let snap = SearchSnapshot {
            nodes: 100,
            incumbent: None,
            nodes_since_improvement: 0,
            time_budget_secs: None,
            elapsed_secs: 1e9,
        };
        assert_eq!(c.decide(&snap), HeuristicDecision::Run);
    }

    /// stall 达标 ⇒ Run（有 incumbent、无命中记录）。
    #[test]
    fn controller_stall_runs() {
        let mut c = HeuristicController::new(HeuristicConfig::default());
        c.notify_root_done(1);
        let snap = SearchSnapshot {
            nodes: 100,
            incumbent: Some(10.0),
            nodes_since_improvement: 50,
            time_budget_secs: None,
            elapsed_secs: 0.0,
        };
        assert_eq!(c.decide(&snap), HeuristicDecision::Run);
        // 差一条不达标 ⇒ Skip（无命中记录时）。
        let snap = SearchSnapshot {
            nodes: 100,
            incumbent: Some(10.0),
            nodes_since_improvement: 49,
            time_budget_secs: None,
            elapsed_secs: 0.0,
        };
        assert_eq!(c.decide(&snap), HeuristicDecision::Skip);
    }

    /// 窗口命中率 ≥ keep_hit_rate ⇒ Run（默认 keep=0.25：命中率 0.5 ⇒ Run）；
    /// 镜像：keep=0.8 时 0.5 < 0.8 ⇒ Skip。
    #[test]
    fn controller_hit_rate_runs() {
        let mut c = HeuristicController::new(HeuristicConfig::default());
        c.notify_root_done(1);
        c.record(true);
        c.record(false);
        c.record(true);
        c.record(false); // 命中率 0.5
        let snap = SearchSnapshot {
            nodes: 100,
            incumbent: Some(10.0),
            nodes_since_improvement: 0,
            time_budget_secs: None,
            elapsed_secs: 0.0,
        };
        assert_eq!(c.hit_rate(), Some(0.5));
        assert_eq!(
            c.decide(&snap),
            HeuristicDecision::Run,
            "默认 keep=0.25：0.5 ≥ 0.25 ⇒ Run"
        );

        let cfg = HeuristicConfig {
            keep_hit_rate: 0.8,
            ..HeuristicConfig::default()
        };
        let mut d = HeuristicController::new(cfg);
        d.notify_root_done(1);
        d.record(true);
        d.record(false);
        d.record(true);
        d.record(false);
        assert_eq!(
            d.decide(&snap),
            HeuristicDecision::Skip,
            "keep=0.8：0.5 < 0.8 ⇒ Skip"
        );
    }

    /// 连 miss 8 次后 Skip（命中率 0 < keep，且 stall 未达、有 incumbent）。
    #[test]
    fn controller_eight_misses_skip() {
        let mut c = HeuristicController::new(HeuristicConfig::default());
        c.notify_root_done(1);
        // 先让一次 decide（stall 达标）授权 dive，随后 8 次 miss。
        let snap = SearchSnapshot {
            nodes: 60,
            incumbent: Some(10.0),
            nodes_since_improvement: 50,
            time_budget_secs: None,
            elapsed_secs: 0.0,
        };
        assert_eq!(c.decide(&snap), HeuristicDecision::Run);
        for _ in 0..8 {
            c.record(false);
        }
        let snap = SearchSnapshot {
            nodes: 100,
            incumbent: Some(10.0),
            nodes_since_improvement: 0,
            time_budget_secs: None,
            elapsed_secs: 0.0,
        };
        assert_eq!(c.hit_rate(), Some(0.0));
        assert_eq!(c.decide(&snap), HeuristicDecision::Skip);
    }

    /// record 之后节流以 dive 所在节点为基准（两次 dive 最小间隔）。
    #[test]
    fn record_advances_throttle_anchor() {
        let mut c = HeuristicController::new(HeuristicConfig::default());
        c.notify_root_done(1);
        let run = SearchSnapshot {
            nodes: 60,
            incumbent: Some(10.0),
            nodes_since_improvement: 50,
            time_budget_secs: None,
            elapsed_secs: 0.0,
        };
        assert_eq!(c.decide(&run), HeuristicDecision::Run);
        c.record(false);
        // dive 在节点 60 运行过：节点 61–91 均节流 Skip（即使 stall 达标）。
        let snap = SearchSnapshot {
            nodes: 91,
            incumbent: Some(10.0),
            nodes_since_improvement: 50,
            time_budget_secs: None,
            elapsed_secs: 0.0,
        };
        assert_eq!(c.decide(&snap), HeuristicDecision::Skip, "61 + 32 − 1 = 92 前均节流");
        let snap = SearchSnapshot {
            nodes: 92,
            incumbent: Some(10.0),
            nodes_since_improvement: 50,
            time_budget_secs: None,
            elapsed_secs: 0.0,
        };
        assert_eq!(c.decide(&snap), HeuristicDecision::Run, "距节点 60 恰 32 ⇒ 放行");
    }

    /// 时间预算保护双分支：
    /// (i) 有预算 + 有 incumbent + elapsed > frac×budget ⇒ Skip（即使 stall）；
    /// (ii) elapsed 恰未超 ⇒ 正常决策 Run；
    /// 附：无预算（None）永不因此 Skip。
    #[test]
    fn controller_budget_protection_branches() {
        let mut c = HeuristicController::new(HeuristicConfig::default());
        c.notify_root_done(1);
        let stall = SearchSnapshot {
            nodes: 100,
            incumbent: Some(10.0),
            nodes_since_improvement: 50,
            time_budget_secs: Some(100.0),
            elapsed_secs: 80.000001, // > 0.8 × 100
        };
        assert_eq!(c.decide(&stall), HeuristicDecision::Skip, "尾段预算保护");
        let not_yet = SearchSnapshot {
            nodes: 100,
            incumbent: Some(10.0),
            nodes_since_improvement: 50,
            time_budget_secs: Some(100.0),
            elapsed_secs: 80.0, // == frac×budget，未超
        };
        assert_eq!(c.decide(&not_yet), HeuristicDecision::Run);
        let no_budget = SearchSnapshot {
            nodes: 100,
            incumbent: Some(10.0),
            nodes_since_improvement: 50,
            time_budget_secs: None,
            elapsed_secs: 1e9,
        };
        assert_eq!(c.decide(&no_budget), HeuristicDecision::Run, "无预算不保护");
        // 有预算但无 incumbent：不保护（仍 Run）。
        let no_inc = SearchSnapshot {
            nodes: 100,
            incumbent: None,
            nodes_since_improvement: 0,
            time_budget_secs: Some(1.0),
            elapsed_secs: 1e9,
        };
        assert_eq!(c.decide(&no_inc), HeuristicDecision::Run);
    }

    /// gap_percent：Max/Min 镜像、正向值。
    #[test]
    fn gap_percent_sense_mirror() {
        let g = gap_percent(Sense::Maximize, 21.0, 21.5).unwrap();
        assert!((g - 100.0 * 0.5 / 21.0).abs() < 1e-12, "Max: (21.5−21)/21×100");
        let g = gap_percent(Sense::Minimize, 21.0, 20.5).unwrap();
        assert!((g - 100.0 * 0.5 / 21.0).abs() < 1e-12, "Min: (21−20.5)/21×100");
        assert_eq!(gap_percent(Sense::Minimize, 40.0, 39.5), Some(1.25));
    }

    /// gap_percent：None 语义（非有限输入；|inc| ≈ 0 且差不 ≈ 0 ⇒ None）。
    #[test]
    fn gap_percent_none_semantics() {
        assert!(gap_percent(Sense::Maximize, f64::NAN, 1.0).is_none());
        assert!(gap_percent(Sense::Maximize, 1.0, f64::NAN).is_none());
        assert!(gap_percent(Sense::Maximize, f64::INFINITY, 1.0).is_none());
        assert!(gap_percent(Sense::Minimize, 1.0, f64::NEG_INFINITY).is_none());
        // |inc| ≈ 0：差 ≈ 0 ⇒ Some(0.0)；差显著 ⇒ None（表示 inf）。
        assert_eq!(gap_percent(Sense::Maximize, 0.0, 0.0), Some(0.0));
        assert_eq!(gap_percent(Sense::Minimize, 0.0, -0.0), Some(0.0));
        assert!(gap_percent(Sense::Maximize, 0.0, 5.0).is_none());
        assert!(gap_percent(Sense::Minimize, 0.0, -5.0).is_none());
    }

    /// gap_percent：负值 clamp（N1）。
    #[test]
    fn gap_percent_negative_clamps_to_zero() {
        // 浮点噪声：est = 39.999999999999996 < inc = 40.0。
        assert_eq!(
            gap_percent(Sense::Maximize, 40.0, 39.999999999999996),
            Some(0.0)
        );
        // Min 镜像：bound = 40.000000000000004 > inc = 40.0。
        assert_eq!(
            gap_percent(Sense::Minimize, 40.0, 40.000000000000004),
            Some(0.0)
        );
        // 显著负差（raw = −7.69% 场景）：Max bound < inc、Min bound > inc。
        assert_eq!(gap_percent(Sense::Maximize, 40.0, 36.924), Some(0.0));
        assert_eq!(gap_percent(Sense::Minimize, 40.0, 43.076), Some(0.0));
        // clamp 后不再出现 -0.00。
        let g = gap_percent(Sense::Maximize, 40.0, 39.999999999999996).unwrap();
        assert_eq!(format!("{g:.2}"), "0.00");
    }

    /// aggregate_bound：空迭代 None；Max 取最大；Min 取最小。
    #[test]
    fn aggregate_bound_empty_and_directions() {
        assert_eq!(aggregate_bound(Sense::Maximize, std::iter::empty()), None);
        assert_eq!(aggregate_bound(Sense::Minimize, std::iter::empty()), None);
        let ests = [1.0, 3.0, 2.0];
        assert_eq!(
            aggregate_bound(Sense::Maximize, ests.iter().copied()),
            Some(3.0)
        );
        assert_eq!(
            aggregate_bound(Sense::Minimize, ests.iter().copied()),
            Some(1.0)
        );
        assert_eq!(aggregate_bound(Sense::Maximize, [7.5].into_iter()), Some(7.5));
    }

    /// format_number：与 CLI objective 行同规则（镜像 main.rs 既有向量）。
    #[test]
    fn format_number_branches() {
        assert_eq!(format_number(0.0), "0.000000");
        assert_eq!(format_number(13.0), "13.000000");
        assert_eq!(format_number(-2.5), "-2.500000");
        assert_eq!(format_number(1e-6), "0.000001", "1e-6 本身不走科学记数");
        assert_eq!(format_number(9.999999e-7), "9.999999e-7");
        assert_eq!(format_number(-1e-7), "-1e-7");
        assert_eq!(format_number(9.99e14), "999000000000000.000000");
        assert_eq!(format_number(1e15), "1e15");
        assert_eq!(format_number(-2.5e18), "-2.5e18");
    }

    /// format_progress_line：精确字符串（数值 / none / inf / 零 gap / 科学记数）。
    #[test]
    fn format_progress_line_exact_strings() {
        let base = ProgressInfo {
            sense: Sense::Maximize,
            elapsed_secs: 1.234,
            nodes: 5,
            incumbent: Some(21.0),
            bound: Some(21.5),
        };
        assert_eq!(
            format_progress_line(&base),
            "progress t=1.23 nodes=5 incumbent=21.000000 bound=21.500000 gap=2.38%"
        );
        // none / inf。
        let none = ProgressInfo {
            sense: Sense::Maximize,
            elapsed_secs: 0.0,
            nodes: 1,
            incumbent: None,
            bound: None,
        };
        assert_eq!(
            format_progress_line(&none),
            "progress t=0.00 nodes=1 incumbent=none bound=none gap=inf"
        );
        // 零 gap（inc == bound）。
        let zero = ProgressInfo {
            sense: Sense::Maximize,
            elapsed_secs: 2.5,
            nodes: 9,
            incumbent: Some(21.0),
            bound: Some(21.0),
        };
        assert_eq!(
            format_progress_line(&zero),
            "progress t=2.50 nodes=9 incumbent=21.000000 bound=21.000000 gap=0.00%"
        );
        // 极小值科学记数 + 精确 gap（1e-7 与 2e-7 均为 2 的幂次缩放，差恰为
        // 1e-7 的最近双精度表示 ⇒ gap 恰为 100）。
        let tiny = ProgressInfo {
            sense: Sense::Maximize,
            elapsed_secs: 0.05,
            nodes: 3,
            incumbent: Some(1e-7),
            bound: Some(2e-7),
        };
        assert_eq!(
            format_progress_line(&tiny),
            "progress t=0.05 nodes=3 incumbent=1e-7 bound=2e-7 gap=100.00%"
        );
        // Min 镜像负 gap 经 clamp 输出 0.00%（不是 -12.50%）。
        let neg = ProgressInfo {
            sense: Sense::Minimize,
            elapsed_secs: 0.1,
            nodes: 2,
            incumbent: Some(40.0),
            bound: Some(45.0),
        };
        assert_eq!(
            format_progress_line(&neg),
            "progress t=0.10 nodes=2 incumbent=40.000000 bound=45.000000 gap=0.00%"
        );
    }

    // —— 策略级自动降级（AutoDegrade；全部伪造时间/界，不依赖真实计时）——

    /// 合格快照（默认参数下 elapsed=60、budget=100 已过 start_frac=0.5，
    /// last_improvement=10 ⇒ 停滞 50 ≥ plateau 25；界任填）。
    fn degrade_snap(sense: Sense, bound: Option<f64>) -> DegradeSnapshot {
        DegradeSnapshot {
            sense,
            elapsed_secs: 60.0,
            budget_secs: Some(100.0),
            has_incumbent: true,
            last_improvement_secs: 10.0,
            bound,
        }
    }

    /// 无预算 / 未到 start_frac / 无 incumbent / 改进尚新 ⇒ 一律不降级。
    #[test]
    fn auto_degrade_guards() {
        let mut d = AutoDegrade::new(DegradePolicy::default());
        // 无预算。
        let mut s = degrade_snap(Sense::Maximize, Some(10.0));
        s.budget_secs = None;
        assert!(!d.decide(&s));
        // 预算非正。
        s.budget_secs = Some(0.0);
        assert!(!d.decide(&s));
        // 未到 start_frac（elapsed=49 < 50）。
        let mut s = degrade_snap(Sense::Maximize, Some(10.0));
        s.elapsed_secs = 49.0;
        assert!(!d.decide(&s));
        // 无 incumbent。
        let mut s = degrade_snap(Sense::Maximize, Some(10.0));
        s.has_incumbent = false;
        assert!(!d.decide(&s));
        // incumbent 改进尚新（停滞 10 < plateau 25）。
        let mut s = degrade_snap(Sense::Maximize, Some(10.0));
        s.last_improvement_secs = 55.0;
        assert!(!d.decide(&s));
        // 界缺失。
        assert!(!d.decide(&degrade_snap(Sense::Maximize, None)));
        // 资格预检口径一致。
        assert!(AutoDegrade::eligible(DegradePolicy::default(), 50.0, Some(100.0)));
        assert!(!AutoDegrade::eligible(DegradePolicy::default(), 49.9, Some(100.0)));
        assert!(!AutoDegrade::eligible(DegradePolicy::default(), 60.0, None));
    }

    /// 双停滞（界不动）⇒ 触发；界有实质进展 ⇒ 不触发且基准被刷新。
    #[test]
    fn auto_degrade_double_stall_triggers() {
        // Max：界自 10.0 起纹丝不动 ⇒ 相对改进 0 < 1% ⇒ 第二次判定触发。
        let mut d = AutoDegrade::new(DegradePolicy::default());
        assert!(!d.decide(&degrade_snap(Sense::Maximize, Some(10.0)))); // 建基准
        assert!(d.decide(&degrade_snap(Sense::Maximize, Some(10.0)))); // 双停滞
        // Min 镜像：界自 20.0 不动 ⇒ 同样触发。
        let mut d = AutoDegrade::new(DegradePolicy::default());
        assert!(!d.decide(&degrade_snap(Sense::Minimize, Some(20.0))));
        assert!(d.decide(&degrade_snap(Sense::Minimize, Some(20.0))));
        // 界实质进展（Max：10→8，相对改进 20% ≥ 1%）⇒ 不触发且基准刷新到 8。
        let mut d = AutoDegrade::new(DegradePolicy::default());
        assert!(!d.decide(&degrade_snap(Sense::Maximize, Some(10.0))));
        assert!(!d.decide(&degrade_snap(Sense::Maximize, Some(8.0))));
        // 基准已刷新：8→7.99（相对改进 0.125% < 1%）⇒ 触发。
        assert!(d.decide(&degrade_snap(Sense::Maximize, Some(7.99))));
        // Min 镜像：20→22（改进 10%）⇒ 不触发并刷新。
        let mut d = AutoDegrade::new(DegradePolicy::default());
        assert!(!d.decide(&degrade_snap(Sense::Minimize, Some(20.0))));
        assert!(!d.decide(&degrade_snap(Sense::Minimize, Some(22.0))));
    }

    /// 自定义参数（测试与注入用）：start=0 / plateau=0 ⇒ 资格即判。
    #[test]
    fn auto_degrade_custom_policy() {
        let cfg = DegradePolicy {
            start_frac: 0.0,
            plateau_frac: 0.0,
            ..DegradePolicy::default()
        };
        let mut d = AutoDegrade::new(cfg);
        let snap = DegradeSnapshot {
            sense: Sense::Maximize,
            elapsed_secs: 0.001,
            budget_secs: Some(10.0),
            has_incumbent: true,
            last_improvement_secs: 0.0,
            bound: Some(21.5),
        };
        assert!(!d.decide(&snap)); // 建基准
        assert!(d.decide(&snap)); // 界不动 ⇒ 触发
    }

    /// 默认参数契约：0.5 / 0.25 / 0.01。
    #[test]
    fn degrade_policy_default_values() {
        let c = DegradePolicy::default();
        assert!((c.start_frac - 0.5).abs() < 1e-12);
        assert!((c.plateau_frac - 0.25).abs() < 1e-12);
        assert!((c.bound_progress_rel - 0.01).abs() < 1e-12);
    }
}
