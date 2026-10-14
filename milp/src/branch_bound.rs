//! 深度优先分支定界（WP4 契约，本文件为 WP4 实现）。
//!
//! 架构：单一 [`crate::simplex::Simplex`] 实例 + Bounds 界栈（trail）+ 深度优先，
//! 节点零拷贝；下行分支先探索；回溯严格按 trail 逆序恢复旧界（恢复后由
//! resolve 的松弛失贴 snap 处理非基内点）。
//!
//! 节点处理（resolve 返回后逐条判定，顺序固定）：
//! (a) Infeasible ⇒ 剪枝回溯；
//! (b) Unbounded：根节点 → 整树终止返回 [`SolveStatus::Unbounded`]；非根节点
//!     （父有界时理论不发生，数值异常时可能出现）⇒ 立即整树终止返回
//!     `Err(MilpError::Numerical)`，不得静默继续；
//! (c) IterationLimit ⇒ 剪枝并置 limit_hit 标记；
//! (d) LP bound 不优于 incumbent 的有效剪枝界（[`effective_prune_bound`]：
//!     整数目标门控 [`crate::model::objective_is_integral`] 时取整数 cutoff
//!     与普通容差的更激进者——Max: max(inc+1−ε, inc+prune_abs+prune_rel·|inc|)、
//!     Min 镜像；非整数目标退化为普通容差 bound ≤ inc + prune_abs +
//!     prune_rel·|inc|）⇒ 剪枝；
//! (e) 整数性满足（全部 int_var |x−round| ≤ tol.integrality）⇒ round_and_check
//!     复核：通过且严格更优 ⇒ 更新 incumbent 并回溯；复核不通过（伪整数解）⇒
//!     不得剪枝、不得记 incumbent，必须按分支规则（MostFractional/Pseudocost）
//!     继续分支；
//! (f) 空界（tighten 返回 false）⇒ 剪枝。
//!
//! 根节点顺序：int_vars 空 ⇒ 单次 LP 直接返回；根 solve 非 Optimal ⇒ 不进割循环
//! 与启发式，直接按 (a)/(b)/(c) 终局；根 Optimal ⇒ （可选）cut_loop_recorded
//! ⇒ （可选）round_and_check + dive ⇒ DFS 主循环。割循环的终态分流：Infeasible
//! 保留整解终局（有效割切空根 LP 是整数不可行的合法证明）；**IterationLimit
//! 不作终态**——丢弃全部根割、重建无割 Simplex 恢复 Optimal 根态后继续启发式
//! 与分支（无割松弛更松的界不剪任何整数可行解；重建前迭代经 base_
//! lp_iterations 累计不丢，stats.cuts 如实清零）。
//!
//! 终态判定（与上限诚实性整合）：主循环因栈空自然退出或因 node_limit/
//! time_limit 提前退出后统一判定——limit_hit 为真或上限触发 ⇒ LimitReached
//! （附 incumbent，若有）；否则有 incumbent ⇒ Optimal、无 ⇒ Infeasible。
//!
//! # 实现说明（WP4）
//! - 节点记账：根 LP 求解与每个树节点的 resolve 各计 1 节点；割循环后取回
//!   根 LP 解、dive 后恢复根热状态属记账性 resolve，不计节点（其迭代数仍
//!   全部计入 `stats.lp_iterations`）。
//! - 树结构：显式 `stack`（待应用的分支决策 [`Pending`]）+ `trail`（已应用
//!   界修改的回放记录）。每个 `Pending` 记录应用前 trail 应回退到的长度
//!   `unwind_to`，保证兄弟节点在父节点界状态上施加自己的分支界；回退严格
//!   按 trail 逆序恢复旧值。
//! - 进展保证：子节点分支界相对当前界严格收缩且新区间非空时才压栈
//!   （[`tighten_child`]）。非分数 x̄ 的两个分支界中至多一侧无收缩，伪整数
//!   回退分支因此在最坏情形下也单调收缩，节点上限兜底。
//! - 上限检查位于栈弹出之后、resolve 之前：栈空 ⇒ 搜索自然耗尽（先于上限
//!   判定，完整搜索不因恰好用满预算而谎报 LimitReached）；弹出待处理决策后
//!   触发上限 ⇒ 搜索确实未完成 ⇒ 诚实置 limit_hit。
//! - 可选优化（整数目标 cutoff 剪枝）：目标项系数全为整数值且全落在
//!   Binary/Integer 变量、obj_offset 为整数值（[`crate::model::objective_is_integral`]）
//!   ⇒ 目标值恒在整数网格上，Max 下任何严格更优解 obj ≥ inc+1（Min 镜像
//!   ≤ inc−1）。此时 (d) 的有效剪枝界取整数 cutoff（inc+1−ε，ε=1e-6·(1+|inc|)
//!   工程余量，见 [`effective_prune_bound`] 契约）与用户容差
//!   prune_abs/prune_rel 的更激进者；同一门控下主循环弹出陈旧前沿
//!   （父节点 LP 目标 est 已不优于有效剪枝界）直接跳过——不 resolve、
//!   不计 nodes、不耗 node_limit（并行侧同口径：池弹出后、reserve 之前）。
//!   非整数目标（连续变量入目标 / 分数系数 / 分数 offset）机制性零影响。
//!
//! # 受控扩展（默认全关 ⇒ 与原实现逐字节一致）
//! [`branch_and_bound_controlled`] 经 [`SearchControl`]（threads / mipgap /
//! adaptive_heuristics / log_period_secs）与 [`SearchHooks`]（进度回调）注入
//! 新能力，既有 [`branch_and_bound`] 签名与行为完全不变：
//! - 顺序路径（threads ≤ 1）全部改动门控：`Pending.est`（压栈时父节点 LP
//!   目标）仅在整数目标门控下参与陈旧前沿弹出剪枝（见实现说明），否则为
//!   纯数据不参与决策；首行进度在根阶段结束的汇合点发射（含纯 LP
//!   捷径与根终局，先发射再返回）；主循环每轮顶部（pop 之前）仅在栈非空时
//!   ①发射周期行（bound = est 聚合，gap 经 [`crate::search_ctrl::gap_percent`]
//!   clamp ≥ 0）与 ②做 mipgap 判定（达标即粘滞跳出并记触发时刻证据）；
//!   自适应启发式在 branch_or_close 的分支路径上（(d) 未剪枝且非整数记录、
//!   选分支变量之前）由控制器决定是否补跑 dive。
//! - threads ≥ 2 组装 [`crate::parallel_bb::ParallelConfig`] 委托
//!   [`crate::parallel_bb::parallel_branch_and_bound`]，并按统一口径映射
//!   ParallelOutcome
//!   （根终局直返；NodeLimit/TimeLimit/Gap ⇒ LimitReached，Gap 附
//!   `gap_terminated = true` 与证据 bound；无 stop ⇒ incumbent?Optimal:Infeasible）。

use std::time::Instant;

use crate::cuts::{self, CutOptions};
use crate::error::{MilpError, Result};
use crate::heuristics::{self, Incumbent};
use crate::model::{objective_is_integral, Problem, RelOp, Sense, SolveStatus, Tolerances, VarType};
use crate::search_ctrl::{
    aggregate_bound, gap_percent, AutoDegrade, DegradePolicy, DegradeSnapshot, Determinism,
    HeuristicConfig, HeuristicController, HeuristicDecision, ProgressInfo, Pseudocosts,
    SearchSnapshot,
};
use crate::simplex::{Bounds, LpSolution, LpStatus};
use crate::simplex_sparse::{new_engine, LpEngine, LpEngineKind};

/// 分支变量选择规则。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BranchRule {
    /// 最分数变量（max min(f, 1−f)，并列取最小下标）。
    MostFractional,
    /// 伪代价分支（[`crate::search_ctrl::Pseudocosts`] 打分：历史增益 ×
    /// 分数度乘积；无观测时退化为 f⁻·f⁺ 乘积规则——对强退化 LP 的顶点
    /// 漂移比 most-fractional 鲁棒，大规模树形显著更小）。
    ///
    /// 【R3 证据裁决（定档矩阵 {mostfrac, pseudo×probes=0, pseudo×probes=8}
    /// × 5 算例 × 3 rep，t1 release，R2 落地后口径）】旧注释「作默认时
    /// marketsplit-4x30 倒退 ~11×」与 baseline.json 隔离实验（pseudocost
    /// −45%）的矛盾已裁决：**11× 倒退 = pseudo×probes=8 口径**（marketsplit
    /// 0.108→1.04s ≈ 9.6×；纯被动 pseudo 反而最优 −54%——第三轮的倒退
    /// 实为可靠性探测预算在该算例等式行结构上的白烧）；mknap-5x300 则
    /// pseudo×probes=8 最优（−25% 证最优墙钟）。故默认走类自适应
    /// （[`BbOptions::branch_rule_auto`] + [`detect_pseudocost_favorable`]）：
    /// 背包/多维背包类（全 Binary + 全正系数 Le 行 + 行数 2..=64，P3 收紧
    /// 单行回落）⇒ Pseudocost + probes=8，其余 ⇒ MostFractional；
    /// marketsplit（等式行）不命中类，探测白烧路径不进默认。
    /// 可靠性冷启动见 [`BbOptions::pc_reliability_probes`]。
    Pseudocost,
}

/// 节点选择规则（首版仅 DFS；best-bound 需按节点存基或冷重启，代价大）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeRule {
    /// 深度优先 + 界栈热启动 + 回溯恢复旧界。
    DepthFirst,
}

/// 分支定界参数。
#[derive(Debug, Clone, Copy)]
pub struct BbOptions {
    /// 分支变量规则。默认 Pseudocost（2026-10 起；可靠性冷启动由
    /// [`Self::pc_reliability_probes`] 门控）。
    pub branch_rule: BranchRule,
    /// 节点选择规则。默认 DepthFirst。
    pub node_rule: NodeRule,
    /// 节点上限（每次 resolve 计 1 node）。默认 2_000_000（2026-10 对偶
    /// 热修复使节点吞吐 ↑ ~30×，旧 10 万上限秒级打满）。
    pub node_limit: u64,
    /// 时间上限毫秒；0 = 不限。默认 0。
    pub time_limit_ms: u64,
    /// 是否启用根节点割循环。默认 true。
    pub root_cuts: bool,
    /// 是否启用启发式。默认 true。
    pub heuristics: bool,
    /// dive 修复变量数上限（深潜预算；对偶热修复下单次 resolve 成本极低，
    /// 大预算可控）。默认 4096。
    pub dive_max_fixes: usize,
    /// 剪枝绝对容差。默认 1e-9。
    pub prune_abs: f64,
    /// 剪枝相对容差。默认 1e-9。
    pub prune_rel: f64,
    /// Pseudocost 可靠性强分支冷启动的每节点探测预算（探测 resolve 次数
    /// 上限；强分支折叠为 Pseudocost 内置能力，不设独立规则枚举值）。
    /// 仅 `branch_rule = Pseudocost` 时生效：对 `n_obs < η` 的分数候选
    /// 做 down+up 双向探测（收紧→resolve→record→恢复）冷启动伪代价
    /// 统计。默认 8；0 = 纯被动（只靠子节点 resolve 被动记账）。预算按
    /// 节点/迭代计、绝不按墙钟（确定性；深度 > 8 的节点不探测；根节点
    /// 按探测日程全量冷启动，见 [`probe_reliability`]）。
    pub pc_reliability_probes: usize,
    /// 子集和 MITM 可行性启发式的窗口重试预算（R5-retry）：全 Binary +
    /// 全 Eq + 整系数类、且根阶段廉价启发式（round/dive）全部 miss 时在根
    /// 一次求解（n ≤ 42 全变量精确；n ≤ 84 窗口半枚举）。默认 8；0 = 关闭
    ///（CLI `--mitm=off` 供消融）。
    pub mitm_windows: usize,
    /// 全局 DFS 工作栈（R6-retry，轴3；仅并行 threads ≥ 2 生效）：true ⇒
    /// 子代直接压全局栈顶、worker 恒从栈顶弹——全局探索序对齐顺序 DFS，
    /// 抑制并行搜索膨胀（本地亲和模式下各 worker 独立 DFS 树的 incumbent
    /// 延迟发现是 t4 节点膨胀实测根因）。字段默认 false；有效值由
    /// [`Self::work_stack_auto`]（W1 按类选栈，多线程默认）在并行派发层
    /// 解析；CLI `--work-stack auto|local|global`。
    pub global_work_stack: bool,
    /// 类自适应工作栈开关（W1）：true 时 [`Self::global_work_stack`] 字段视为
    /// 「未指定」，并行派发层（[`parallel_entry`]，经
    /// [`resolve_work_stack_mode`]）按 [`detect_pseudocost_favorable`]（与 R3
    /// 分支规则类自适应**同一检测器**）解析——命中（全 Binary + 正系数 Le
    /// 行 + 行数 2..=64，即 mknap 多维背包族）⇒ global（P2 §8.2/8.3 实测：
    /// mknap-5x300 t4 节点中位 −14.3%、mknap-3x100 −26%~−45%）；未命中（单行
    /// 背包 / 等式行 marketsplit / Ge 行 setcover / 连续列 milp-mixed 类）⇒
    /// local（global 在 knapsack-500 回退 ×1.2-3.0、setcover-100x150
    /// +8%~+26%）。CLI `--work-stack auto|local|global`，多线程默认 auto；
    /// 显式 local/global（auto=false）覆盖优先且不过小算例调度门（保持栈
    /// 形态 A/B 可测量）。解析只发生在并行派发层：threads=1 顺序路径不读
    /// 该字段（t1 行为逐字节不变）。默认 true。
    pub work_stack_auto: bool,
    /// 树内局部割（W6）：深度 ≤ 8 的节点上按配额跑轻量局部割轮（复用根
    /// 割批序；每节点 ≤ LOCAL_CUT_PER_NODE 条、每 LOCAL_CUT_NODE_INTERVAL
    /// 个节点一次、全局预算封顶），子树作用域经 Simplex 快照恢复（冷审
    /// 调整 6 的机制 spike 结论）。仅顺序路径；并行 worker 不施加（窃取
    /// 交错会使局部割泄漏到子树外——共享 Simplex 无行作用域）。默认 true。
    pub local_cuts: bool,
    /// LP 引擎选择（S3）：默认 [`LpEngineKind::Tableau`]（本轮不改默认；
    /// CLI `--lp-engine tableau|revised` / 库字段两层接线）。`Revised` =
    /// 稀疏修正单纯形（simplex_sparse.rs，CSC+稀疏 LU+ETA+对偶热修复），
    /// **仅接线顺序路径**：threads ≥ 2 在 [`branch_and_bound_controlled`]
    /// 入口显式报错（并行 worker 直接构造稠密 Simplex——F5 冷审：静默
    /// 忽略旗标会形成文档外行为分叉）。
    pub lp_engine: LpEngineKind,
    /// 类自适应默认开关（R3）：true 时 `branch_rule` 字段视为「未指定」，
    /// 求解入口按问题结构类检测（[`detect_pseudocost_favorable`]）解析
    /// 有效规则——背包/多维背包类（全 Binary + 全正系数 Le 行 + 行数
    /// 2..=64，P3 收紧单行回落）⇒ Pseudocost + `pc_reliability_probes ≥ 8`，
    /// 否则 MostFractional；false 时按 `branch_rule` 字段逐字生效（显式覆盖
    /// 优先）。默认 true。
    pub branch_rule_auto: bool,
}

impl Default for BbOptions {
    /// 默认值：MostFractional / DepthFirst / 2_000_000 / 0 / true / true / 4096 /
    /// 1e-9 / 1e-9 / 8 / global_work_stack=false + work_stack_auto=true（W1：
    /// 多线程默认按类选栈）。（2026-10 复核后回退：第三轮曾默认 Pseudocost，
    /// 但 market-split/mknap-5x300 等类实测大幅倒退且无实测收益——全部
    /// 已验证性能数字均在 MostFractional 口径下取得；Pseudocost 经
    /// `BbOptions::branch_rule` / CLI `--branch-rule pseudocost` 可选。）
    fn default() -> Self {
        BbOptions {
            branch_rule: BranchRule::MostFractional,
            node_rule: NodeRule::DepthFirst,
            node_limit: 2_000_000,
            time_limit_ms: 0,
            root_cuts: true,
            heuristics: true,
            dive_max_fixes: 4096,
            prune_abs: 1e-9,
            prune_rel: 1e-9,
            pc_reliability_probes: 8,
            mitm_windows: 8,
            global_work_stack: false,
            work_stack_auto: true,
            branch_rule_auto: true,
            local_cuts: true,
            lp_engine: LpEngineKind::default(),
        }
    }
}

/// 搜索控制旋钮（全部默认 = 现行为；纯数据，Debug+Clone+Copy 安全）。
///
/// 经 [`branch_and_bound_controlled`] 注入；[`SearchControl::default`] 与
/// `hooks = None` 的组合下顺序路径与既有 [`branch_and_bound`] 逐字节一致
/// （所有新能力均门控，默认全关）。
#[derive(Debug, Clone, Copy)]
pub struct SearchControl {
    /// 线程数。默认 1；1 = 纯顺序路径（与现版本一致），≥2 组装
    /// [`crate::parallel_bb::ParallelConfig`] 委托并行引擎。
    pub threads: usize,
    /// MIP gap 百分比阈值。默认 None = 关闭；Some(0.0) 等价关闭；Some(p) 且
    /// p > 0 时剩余对偶界与 incumbent 的 gap（一律 clamp ≥ 0，见
    /// [`crate::search_ctrl::gap_percent`]）≤ p% 即提前终止，终态诚实标注
    /// LimitReached 且 `gap_terminated = true`。
    pub mipgap_percent: Option<f64>,
    /// 是否启用自适应启发式（树节点分支路径上按控制器决策补跑 dive）。
    /// 默认 false（现行为：仅根节点启发式）。
    pub adaptive_heuristics: bool,
    /// 周期进度行的间隔秒数。默认 None = 不输出周期行（首行仍随钩子发出）。
    pub log_period_secs: Option<f64>,
    /// 确定性模式（W3，默认 [`Determinism::Wall`] = 现状）：budget = 纯
    /// 预算制——全部墙钟判定改为等价节点数预算（详见 [`Determinism`]），
    /// 决策路径零墙钟读点、仅支持 threads=1。
    pub determinism: Determinism,
    /// 策略级自动降级参数（「继续证最优」vs「提前返回最好解」）。仅在
    /// `adaptive_heuristics = true` 且给了时间预算（`BbOptions::time_limit_ms
    /// != 0`）时评估；触发时与 mipgap 同径提前终止（LimitReached +
    /// gap_terminated），另附 `auto_degraded = true` 标记。默认
    /// [`DegradePolicy::default`]。
    pub degrade: DegradePolicy,
}

impl Default for SearchControl {
    /// 默认值：threads=1 / mipgap_percent=None / adaptive_heuristics=false /
    /// log_period_secs=None / determinism=Wall / degrade=DegradePolicy::default()
    ///（= 现版本行为，零变化）。
    fn default() -> Self {
        SearchControl {
            threads: 1,
            mipgap_percent: None,
            adaptive_heuristics: false,
            log_period_secs: None,
            determinism: Determinism::Wall,
            degrade: DegradePolicy::default(),
        }
    }
}

/// 进度回调注入（仅 derive(Clone)：回调闭包不可 Copy，derive(Copy) 触发 E0277）。
#[derive(Clone)]
pub struct SearchHooks {
    /// 进度回调（`Arc<dyn Fn(&ProgressInfo) + Send + Sync>`）。
    pub progress: crate::search_ctrl::ProgressCallback,
}

/// 分支定界统计。
#[derive(Debug, Clone, Copy, Default)]
pub struct BbStats {
    /// 处理节点数（每次 resolve 计 1）。
    pub nodes: u64,
    /// 累计单纯形迭代数。
    pub lp_iterations: u64,
    /// 施加割数。
    pub cuts: u64,
    /// 累计耗时毫秒。
    pub time_ms: u64,
}

/// 分支定界结果。
#[derive(Debug, Clone)]
pub struct BbResult {
    /// 求解终态。gap 提前终止时 = LimitReached（诚实：未证明最优）。
    pub status: SolveStatus,
    /// 最优目标（原始 sense + obj_offset 口径；LimitReached 时为已有 incumbent）。
    pub best_objective: Option<f64>,
    /// 最优解（与入参 problem 同维同序，未经 postsolve；整数变量已取整）。
    pub best_values: Option<Vec<f64>>,
    /// 统计信息。
    pub stats: BbStats,
    /// 求解结束时的对偶界（原始 sense + obj_offset 口径）。gap 终止 =
    /// 触发时刻证据 db（粘滞快照，不依赖后续栈态）；截断（上限触发）= 剩余
    /// 待处理子树的 est 聚合（栈空时可 None）；完备且有 incumbent =
    /// best_objective；Infeasible/Unbounded = None。
    pub bound: Option<f64>,
    /// 是否因达到 mipgap 阈值提前终止（true 时 status 必为 LimitReached，
    /// bound 为触发时刻证据）。
    pub gap_terminated: bool,
    /// 是否因策略级自动降级提前终止（auto 模式 + 时间预算 + 双停滞触发）。
    /// true ⇒ `gap_terminated` 必为 true、status 必为 LimitReached——与
    /// mipgap 同径的诚实提前终止，此字段区分触发原因。
    pub auto_degraded: bool,
}

/// 深度优先分支定界主入口（签名与行为完全不变；默认受控路径的薄委托）。
pub fn branch_and_bound(problem: &Problem, opts: &BbOptions, tol: Tolerances) -> Result<BbResult> {
    branch_and_bound_controlled(problem, opts, &SearchControl::default(), tol, None)
}

/// Pseudocost 类自适应默认的结构检测（R3；P3 收紧行数下限）：**全变量
/// Binary、全部约束为正系数 Le 行、2 ≤ 行数 ≤ 64** 同时成立——强结构特征
/// 覆盖 mknap-* 家族（多维背包）；setcover（Ge 行 + 混合符号目标语义）、
/// marketsplit（等式行）、milp-mixed（连续列 + 等式行）均不命中，回落
/// MostFractional。条件取强特征交集，宁缺勿滥（置信不足时回落 mostfrac，
/// 风险提示原文）。
///
/// 行数下限 2（P3，2026-10-13 同窗 A/B 6 rep）：单行背包类
/// （knapsack-500/2000/5000）在 Pseudocost+根探测下零收益纯付探测开销
/// （墙钟 +48%/持平/+26%，节点 −33%/−2.5%/−29% 全在 MostFractional 侧），
/// 而多行类收益真实（mknap-3x100 墙钟 −67%、节点 −73%）——下限收紧到 2
/// 保留多行收益、消除单行开销；objective 全程已知最优不变。
pub fn detect_pseudocost_favorable(problem: &Problem) -> bool {
    if problem.vars.is_empty() || problem.vars.iter().any(|v| v.var_type != VarType::Binary) {
        return false; // 非「全变量 Binary」
    }
    if problem.constraints.len() < 2 || problem.constraints.len() > 64 {
        return false; // 行数 < 2（单行背包类回落 mostfrac，见模块注释 P3）或 > 64
    }
    problem.constraints.iter().all(|con| {
        con.op == RelOp::Le && con.expr.terms.iter().all(|&(_, a)| a > 0.0)
    })
}

/// 求解入口解析类自适应默认（R3 分支规则）：`branch_rule_auto = true` 时按
/// [`detect_pseudocost_favorable`] 翻转有效规则（命中 ⇒ Pseudocost 且
/// `pc_reliability_probes ≥ 8`；未命中 ⇒ MostFractional），返回解析后的
/// 具体选项（`branch_rule_auto` 置 false，幂等）；显式指定
/// （`branch_rule_auto = false`）原样透传——覆盖优先。
///
/// W1 按类选栈（`work_stack_auto`）的解析**不在本函数**——工作栈是纯并行
/// 概念（threads=1 顺序路径不读该字段），解析连同小算例调度门一起在
/// [`parallel_entry`]（并行派发层）完成：auto ⇒ 类检测翻转 + 调度门；
/// 显式 local/global ⇒ 逐字生效且不进门（保持 P2 式栈形态 A/B 可测量）。
fn resolve_branch_rule_defaults(problem: &Problem, opts: &BbOptions) -> BbOptions {
    if !opts.branch_rule_auto {
        return *opts;
    }
    let mut resolved = *opts;
    resolved.branch_rule_auto = false;
    if detect_pseudocost_favorable(problem) {
        resolved.branch_rule = BranchRule::Pseudocost;
        resolved.pc_reliability_probes = resolved.pc_reliability_probes.max(8);
    } else {
        resolved.branch_rule = BranchRule::MostFractional;
    }
    resolved
}

/// W1 按类选栈解析（纯函数，仅 [`parallel_entry`] 调用；可单元级直接测试）：
/// `work_stack_auto = true`（多线程默认）⇒ 按 [`detect_pseudocost_favorable`]
///（与 R3 分支规则**同一检测器**）翻转 [`BbOptions::global_work_stack`]——
/// 命中（全 Binary + 正系数 Le 行 + 行数 2..=64，即 mknap 多维背包族）⇒
/// global（P2 §8.2/8.3 实测 mknap 族 t4 节点 −14%~−45%）；未命中（单行
/// 背包 / 等式行 marketsplit / Ge 行 setcover / 连续列 milp-mixed 类）⇒
/// local（global 在 knapsack-500 回退 ×1.2-3.0、setcover-100x150
/// +8%~+26%）。显式（auto=false）逐字透传（覆盖优先）。返回解析后选项与
/// 「是否处于 auto 决策链」（调度门只对 auto 链生效）。
fn resolve_work_stack_mode(problem: &Problem, opts: &BbOptions) -> (BbOptions, bool) {
    let mut resolved = *opts;
    let was_auto = resolved.work_stack_auto;
    if was_auto {
        resolved.work_stack_auto = false;
        resolved.global_work_stack = detect_pseudocost_favorable(problem);
    }
    (resolved, was_auto)
}

/// 受控入口：threads ≤ 1 走纯顺序路径；threads ≥ 2 组装
/// [`crate::parallel_bb::ParallelConfig`] 委托
/// [`crate::parallel_bb::parallel_branch_and_bound`] 并
/// 按统一口径映射 ParallelOutcome：
/// - 根阶段终局 ⇒ 直返；
/// - stop ∈ {NodeLimit, TimeLimit}（或顺序侧 IterationLimit 类截断）⇒
///   LimitReached（附 incumbent，若有）；
/// - stop = Gap ⇒ LimitReached 且 `gap_terminated = true`（bound/gap_percent
///   用触发时刻证据口径）；
/// - stop = None ⇒ 有 incumbent 则 Optimal、无则 Infeasible（粘滞 stop 优先
///   于完备判定）。
///
/// 默认（[`SearchControl::default`] + `hooks = None`）与 [`branch_and_bound`]
/// 逐字节一致。
pub fn branch_and_bound_controlled(
    problem: &Problem,
    opts: &BbOptions,
    ctrl: &SearchControl,
    tol: Tolerances,
    hooks: Option<&SearchHooks>,
) -> Result<BbResult> {
    // R3 类自适应默认（分支规则）：入口一次解析（顺序/并行两路共用；worker
    // 收到的即具体规则，确定性）。显式 branch_rule（auto=false）覆盖优先。
    // W1 按类选栈的解析在 parallel_entry（纯并行概念，t1 不读该字段）。
    let opts = &resolve_branch_rule_defaults(problem, opts);
    if ctrl.threads >= 2 {
        // S3/F5 防御：revised 引擎仅接线顺序路径（并行 worker 直接构造稠密
        // Simplex）；对直接 API 调用者显式报错而非静默回 tableau（静默忽略
        // 旗标会形成文档外行为分叉）。CLI 层已在解析层给出用法错误。
        if opts.lp_engine == LpEngineKind::Revised {
            return Err(MilpError::InvalidModel(
                "--lp-engine revised 仅支持 threads=1（S3 接线范围：并行引擎直接构造稠密 Simplex）"
                    .to_string(),
            ));
        }
        return parallel_entry(problem, opts, ctrl, tol, hooks);
    }
    sequential_bb(problem, opts, ctrl, tol, hooks)
}

/// 顺序分支定界（既有 [`branch_and_bound`] 的实现体 + 全部门控的受控扩展；
/// 门控默认全关 ⇒ 行为与原实现逐字节一致）。
fn sequential_bb(
    problem: &Problem,
    opts: &BbOptions,
    ctrl: &SearchControl,
    tol: Tolerances,
    hooks: Option<&SearchHooks>,
) -> Result<BbResult> {
    let start = Instant::now();
    let int_vars = problem.integer_vars();
    let sense = problem.sense;
    // 整数目标门控（cutoff 剪枝 + est 弹出剪枝共用；纯问题属性，算一次）。
    let obj_integral = objective_is_integral(problem);
    let mut stats = BbStats::default();
    // —— W3 确定性模式：budget = 纯预算制——time_limit 换算为节点预算，
    // 决策通道的进度量单位从墙钟秒改为节点数（progress_units 注入），节点
    // 边界截断改查节点预算；决策路径零墙钟读点（展示通道——进度行时刻/
    // stats.time_ms——仍读钟，不参与任何判定）。wall（默认）零变化。
    let budget_mode = ctrl.determinism == Determinism::Budget;
    // budget 模式的等价节点预算（time_limit_ms > 0 才有；0 = 不限 ⇒ None，
    // 降级与预算保护关闭，行为 ≡ 无时限）。
    let budget_node_cap = (budget_mode && opts.time_limit_ms > 0).then(|| budget_node_cap_of(opts));
    // 决策通道的预算量（单位与进度量一致）：budget = 节点数 / wall = 秒。
    let budget_progress_cap: Option<f64> = if budget_mode {
        budget_node_cap.map(|n| n as f64)
    } else if opts.time_limit_ms == 0 {
        None
    } else {
        Some(opts.time_limit_ms as f64 / 1000.0)
    };
    // 割丢弃回退（根割 IterationLimit 重建 Simplex）前累计的 LP 迭代数：
    // finish 以 base + 当前 Simplex 迭代数记账，重建不丢根 LP + 割循环迭代。
    let mut base_lp_iterations: u64 = 0;
    let mut simplex = new_engine(problem, tol, opts.lp_engine)?;
    // W2 内部墙钟时限：0（不限）= None ⇒ 零墙钟读点（条件读钟）；重建
    // Simplex（预割轮 / R8 回滚）同款重注入。W3：budget 模式不注入（决策
    // 路径零墙钟读点契约；时限语义已换算为节点预算）。
    let time_deadline = if budget_mode {
        None
    } else {
        bb_time_deadline(&start, opts)
    };
    simplex.set_time_deadline(time_deadline);
    let mut bounds = Bounds::from_problem(problem);
    // incumbent = (目标, 取值)；目标一律按取整后的值评估（WP3 复核口径）。
    let mut incumbent: Option<(f64, Vec<f64>)> = None;
    // —— 受控扩展状态（默认全关 ⇒ 零行为变化）——
    // 上行进度行的发射时刻（秒）；仅钩子启用时被写入。
    let mut last_log: Option<f64> = None;
    // gap 达标提前终止（粘滞）与触发时刻证据 db（与并行证据口径一致，不依赖
    // 后续栈态）。
    let mut gap_hit = false;
    let mut gap_evidence: Option<f64> = None;
    // 自适应启发式状态（adaptive_heuristics=false ⇒ 恒 None，零开销零行为）。
    let mut adaptive = ctrl
        .adaptive_heuristics
        .then(|| AdaptiveState::new(HeuristicConfig::default(), ctrl.degrade));

    // —— 根 LP（计 1 节点）——
    let mut sol = simplex.solve(&bounds);
    stats.nodes = 1;

    // 纯 LP 捷径：无整数变量 ⇒ 单次 LP 直接按其状态返回。
    if int_vars.is_empty() {
        let status = match sol.status {
            LpStatus::Optimal => SolveStatus::Optimal,
            LpStatus::Infeasible => SolveStatus::Infeasible,
            LpStatus::Unbounded => SolveStatus::Unbounded,
            // W2：纯 LP 单节点正是节点边界检查的盲区——根 LP 内部墙钟超限
            // ⇒ 诚实 LimitReached（不冒充 Optimal；纯 LP 无 incumbent 概念，
            // 不携带 best）。
            LpStatus::IterationLimit | LpStatus::TimeLimit => SolveStatus::LimitReached,
        };
        let lp_obj = sol.objective;
        let incumbent = if sol.status == LpStatus::Optimal {
            Some((sol.objective, sol.x))
        } else {
            None
        };
        // 完备（Optimal）⇒ bound = 最优值；其余终局无对偶界可报。
        let bound = if status == SolveStatus::Optimal {
            Some(lp_obj)
        } else {
            None
        };
        // 首行进度：根阶段终局汇合点之一（先发射再返回）。
        if let Some(h) = hooks {
            emit_progress(
                h,
                &start,
                stats.nodes,
                incumbent.as_ref().map(|(o, _)| *o),
                bound,
                sense,
            );
        }
        return Ok(finish(
            status,
            incumbent,
            &mut stats,
            &simplex,
            start,
            base_lp_iterations,
            bound,
            false,
            false,
        ));
    }

    // 根 solve 非 Optimal ⇒ 不进割循环与启发式（gomory_cuts 的前提是 Optimal
    // 表格），直接按 (a)/(b)/(c) 终局。
    if let Some(status) = root_terminal(classify_node(sol.status, true))? {
        // R8（③）终态 incumbent 提取：根 solve 以 IterationLimit 或 TimeLimit
        //（W2：中断点同为原始可行——比值测试每步保持 rhs 在界内、Phase-1
        // 已完成，x 非最优但可行）收场时，round_and_check（原始问题逐约束
        // + 整数性复核，不依赖表格语义）通过则作为 incumbent 随
        // LimitReached 携带返回：终态诚实不变（仍 limit_reached），但不再
        // 以零 incumbent 传播（setcover-300x600/400x800 类根 LP 上限锁的
        // 失败模式）。Infeasible/Unbounded 不提取（不可行信号不得吞）；提取
        // 失败（分数点取整违反约束）如实保持无 incumbent。
        let extracted = if status == SolveStatus::LimitReached
            && matches!(sol.status, LpStatus::IterationLimit | LpStatus::TimeLimit)
            && opts.heuristics
        {
            heuristics::round_and_check(problem, &bounds, &sol.x, &tol)
        } else {
            None
        };
        // 首行进度：根阶段终局汇合点（先发射再返回；此时尚无前沿栈 ⇒ bound=none）。
        if let Some(h) = hooks {
            emit_progress(h, &start, stats.nodes, incumbent.as_ref().map(|(o, _)| *o), None, sense);
        }
        let extracted_pair = extracted.map(|inc| (inc.objective, inc.values));
        return Ok(finish(
            status,
            extracted_pair,
            &mut stats,
            &simplex,
            start,
            base_lp_iterations,
            None,
            false,
            false,
        ));
    }

    // —— 根预割启发式（仅 root_cuts 开启时额外执行一轮）——
    // 割循环可能把根 LP 顶点从「单断点近整数」推向「多分数中面」（界更紧
    // 但取整/潜水质量下降——max 背包类：纯 LP 贪心顶点仅容量断点一个分数
    // 变量，dive 的断点回退链恰好在其上跨过；割后中面顶点上同一 dive 只得
    // 弱解）。割前顶点上先跑一轮 round_and_check + dive、与割后一轮取最
    // 优，两轮互为补充：incumbent 只增不减（maybe_update），割生成只读表格
    // 不读 incumbent，割序不受影响。关割路径顶点唯一（割前 ≡ 割后），
    // `root_cuts=false` 时不执行 ⇒ 行为逐字节不变。
    if opts.root_cuts && opts.heuristics {
        if let Some(inc) = heuristics::round_and_check(problem, &bounds, &sol.x, &tol) {
            maybe_update(sense, &mut incumbent, inc);
        }
        if let Some(inc) = heuristics::dive(
            problem,
            &mut simplex,
            &bounds,
            &sol.x,
            opts.dive_max_fixes,
            &tol,
        ) {
            maybe_update(sense, &mut incumbent, inc);
        }
        // 子集和 MITM（R5-retry）：全 Eq + 全 Binary 类的最后手段——
        // round/dive（及一切 LP 驱动机制）在该类零产出（marketsplit 大规格
        // 200 万节点零可行）；MITM 脱离 LP 信息源精确/窗口求解。仅当廉价
        // 启发式全部 miss（incumbent 为空）才启动，非该类恒零成本跳过。
        if opts.heuristics
            && opts.mitm_windows > 0
            && incumbent.is_none()
        {
            if let Some(inc) = heuristics::subset_sum_mitm(
                problem,
                &bounds,
                &sol.x,
                opts.mitm_windows,
                &tol,
            ) {
                maybe_update(sense, &mut incumbent, inc);
            }
        }
        // dive 扰动了 Simplex 内部基：割循环入口的防御性 resolve(bounds) 只
        // 恢复 Optimal 性，不恢复**同一个**最优基——重退化根 LP（如背包类）
        // 的备择最优基上 GMI 分离结果不同，割集会系统性偏离无预割路径
        // （实测 knapsack-2000：12→4 割、594→1628 节点；knapsack-500：
        // 360→1614 节点）。此处重建冷 Simplex + 重解根 LP（确定性求解器 ⇒
        // 与无预割路径逐字节一致的表格），使割循环恒在冷基上分离；incumbent
        // 已由 maybe_update 收取，不受重建影响。迭代记账先累计再重建
        // （与下方 PruneLimit 重建路径同款），不丢数。
        base_lp_iterations = base_lp_iterations.saturating_add(simplex.iterations());
        simplex = new_engine(problem, tol, opts.lp_engine)?;
        simplex.set_time_deadline(time_deadline);
        sol = simplex.solve(&bounds);
        if let Some(status) = root_terminal(classify_node(sol.status, true))? {
            // 首行进度：根阶段终局汇合点（先发射再返回）。
            if let Some(h) = hooks {
                emit_progress(h, &start, stats.nodes, incumbent.as_ref().map(|(o, _)| *o), None, sense);
            }
            return Ok(finish(
                status,
                None,
                &mut stats,
                &simplex,
                start,
                base_lp_iterations,
                None,
                false,
                false,
            ));
        }
    }

    // —— 根割循环（可选；割行全局保留至树结束）——
    if opts.root_cuts {
        let (applied, st, _recorded) = cuts::cut_loop_recorded(
            &mut simplex,
            &bounds,
            problem,
            &int_vars,
            CutOptions::default(),
            &tol,
            || new_engine(problem, tol, opts.lp_engine),
        )?;
        match classify_node(st, true) {
            // (c) 割后 resolve 触发迭代上限：割循环停在非 Optimal 态，不可作
            // 为树搜索的根态。不走 finish——丢弃全部根割，重建无割 Simplex
            // + solve 恢复 Optimal 根态（进割循环前根 LP 必为 Optimal，重建
            // 复现同一状态）后继续启发式与分支。正确性：无割原始松弛的可行
            // 域 ⊇ 含割版本（有效割不切任何整数可行解）⇒ 回退后的完备搜索
            // 覆盖全部整数可行解，最坏等同于关割行为（不误终态、不误剪）。
            NodeAction::PruneLimit => {
                // 记账完整：重建前 Simplex 的累计迭代（根 LP + 割循环）不丢。
                base_lp_iterations = base_lp_iterations.saturating_add(simplex.iterations());
                stats.cuts = 0; // 割已随 Simplex 一起丢弃，如实清零
                simplex = new_engine(problem, tol, opts.lp_engine)?;
                simplex.set_time_deadline(time_deadline);
                sol = simplex.solve(&bounds);
                if let Some(status) = root_terminal(classify_node(sol.status, true))? {
                    // 首行进度：根阶段终局汇合点（先发射再返回）。
                    if let Some(h) = hooks {
                        emit_progress(
                            h,
                            &start,
                            stats.nodes,
                            incumbent.as_ref().map(|(o, _)| *o),
                            None,
                            sense,
                        );
                    }
                    return Ok(finish(
                        status,
                        None,
                        &mut stats,
                        &simplex,
                        start,
                        base_lp_iterations,
                        None,
                        false,
                        false,
                    ));
                }
            }
            _ => {
                stats.cuts = applied;
                if let Some(status) = root_terminal(classify_node(st, true))? {
                    // 首行进度：根阶段终局汇合点（先发射再返回；此时尚无前沿栈 ⇒ bound=none）。
                    if let Some(h) = hooks {
                        emit_progress(
                            h,
                            &start,
                            stats.nodes,
                            incumbent.as_ref().map(|(o, _)| *o),
                            None,
                            sense,
                        );
                    }
                    return Ok(finish(
                        status,
                        None,
                        &mut stats,
                        &simplex,
                        start,
                        base_lp_iterations,
                        None,
                        false,
                        false,
                    ));
                }
                // 取回割后的根 LP 解（记账性 resolve，不计节点）。
                sol = simplex.resolve(&bounds);
                if let Some(status) = root_terminal(classify_node(sol.status, true))? {
                    // 首行进度：根阶段终局汇合点（先发射再返回；此时尚无前沿栈 ⇒ bound=none）。
                    if let Some(h) = hooks {
                        emit_progress(
                            h,
                            &start,
                            stats.nodes,
                            incumbent.as_ref().map(|(o, _)| *o),
                            None,
                            sense,
                        );
                    }
                    return Ok(finish(
                        status,
                        None,
                        &mut stats,
                        &simplex,
                        start,
                        base_lp_iterations,
                        None,
                        false,
                        false,
                    ));
                }
            }
        }
    }

    // —— 根启发式（可选）：round_and_check + dive + MITM 供应初始 incumbent ——
    if opts.heuristics {
        if let Some(inc) = heuristics::round_and_check(problem, &bounds, &sol.x, &tol) {
            maybe_update(sense, &mut incumbent, inc);
        }
        if let Some(inc) = heuristics::dive(
            problem,
            &mut simplex,
            &bounds,
            &sol.x,
            opts.dive_max_fixes,
            &tol,
        ) {
            maybe_update(sense, &mut incumbent, inc);
        }
        if opts.mitm_windows > 0 && incumbent.is_none() {
            if let Some(inc) = heuristics::subset_sum_mitm(
                problem,
                &bounds,
                &sol.x,
                opts.mitm_windows,
                &tol,
            ) {
                maybe_update(sense, &mut incumbent, inc);
            }
        }
        // dive 扰动了 Simplex 内部基：用真实界恢复根热状态（记账性 resolve）。
        sol = simplex.resolve(&bounds);
        if let Some(status) = root_terminal(classify_node(sol.status, true))? {
            // 首行进度：根阶段终局汇合点（先发射再返回；此时尚无前沿栈 ⇒ bound=none）。
            if let Some(h) = hooks {
                emit_progress(h, &start, stats.nodes, incumbent.as_ref().map(|(o, _)| *o), None, sense);
            }
            return Ok(finish(
                status,
                None,
                &mut stats,
                &simplex,
                start,
                base_lp_iterations,
                None,
                false,
                false,
            ));
        }
    }

    // 根启发式供应的 incumbent 记为根阶段末尾时刻的改进（降级停滞锚点；
    // W3：budget 模式锚点单位 = 节点数——根 = 1）。
    if let Some(state) = adaptive.as_mut() {
        if incumbent.is_some() {
            state.last_improve_progress = progress_units(&start, stats.nodes, budget_mode);
        }
    }

    // —— DFS 主循环状态 ——
    let mut stack: Vec<Pending> = Vec::new();
    let mut trail: Vec<TrailEntry> = Vec::new();
    let mut limit_hit = false;
    // 策略级自动降级触发标记（与 gap_hit 同径终止，报告侧另附原因）。
    let mut degraded = false;
    // 伪代价分支统计（仅 Pseudocost 规则读写；MostFractional 下为零开销
    // 空表）。
    let mut pc = Pseudocosts::new(problem.vars.len());
    // —— W6 树内局部割状态 ——
    // 作用域栈：(施加时 trail 长度, Simplex 快照, 快照时累计迭代)。DFS 弹出
    // unwind_to < 施加点的 Pending ⇒ 已退出该子树 ⇒ 恢复快照（局部行只在
    // 子树内有效——泄漏会切子树外整数可行解）。恢复前把迭代差额累计进
    // base（快照回退 iterations_total，分离/重解迭代不丢账）。
    let mut cut_scopes: Vec<(usize, Box<dyn LpEngine>, u64)> = Vec::new();
    let mut local_cut_budget = LOCAL_CUT_GLOBAL_BUDGET;
    let mut next_local_cut_node: u64 = LOCAL_CUT_NODE_INTERVAL;
    // W6 类门控：仅多维背包类（detect_pseudocost_favorable——与 R3/W1 同一
    // 检测器）施加局部割。实测依据：mknap 族收益真实（5x300 节点 −17%、
    // 3x100 −11%）；marketsplit-4x40 零压缩纯开销（节点 −1%、墙钟中位
    // 10.74→12.98s +21%，5 rep 区间不相交——28 条局部割无树压缩）；
    // 其余类树形/根闭合本来不触发。单行背包类（knapsack/hard-knapsack）
    // 不命中门 ⇒ 行为逐字节不变。
    let local_cuts_on = opts.local_cuts && detect_pseudocost_favorable(problem);

    // 根节点按 (d)/(e)/(f) 判定（无界等已终局，此处的 Unbounded 恒为根语义）。
    branch_or_close(
        problem,
        &int_vars,
        &mut bounds,
        &sol,
        tol,
        opts,
        obj_integral,
        &mut incumbent,
        &mut stack,
        &mut trail,
        &mut simplex,
        stats.nodes,
        progress_units(&start, stats.nodes, budget_mode),
        budget_progress_cap,
        &mut pc,
        adaptive.as_mut(),
        0,
    );

    // 根阶段完成（根启发式刚运行过）：控制器节流锚点 = 根节点数（=1）。
    if let Some(state) = adaptive.as_mut() {
        state.controller.notify_root_done(stats.nodes);
    }

    // —— 根阶段结束：首行进度发射点（钩子启用时无条件发射，口径与并行一致；
    // 先发射再进入主循环。栈空时若已有 incumbent，搜索事实上完备 ⇒ bound =
    // incumbent（gap=0），杜绝 incumbent=X bound=none gap=inf 矛盾行）——
    if let Some(h) = hooks {
        let inc_obj = incumbent.as_ref().map(|(o, _)| *o);
        let bound = frontier_bound(sense, &stack, inc_obj);
        emit_progress(h, &start, stats.nodes, inc_obj, bound, sense);
        last_log = Some(start.elapsed().as_secs_f64());
    }

    loop {
        // —— 受控扩展（默认全关 ⇒ 与原实现逐字节一致）。①/② 仅在栈非空时
        // 参与：空栈轮次直接落到下方 pop 的自然 break（Optimal/Infeasible），
        // 杜绝 incumbent=X bound=none gap=inf 矛盾周期行——
        if !stack.is_empty() && !gap_hit {
            // ① 周期进度行：log 启用且距上行 ≥ log_period_secs 才发射。
            if let (Some(h), Some(period)) = (hooks, ctrl.log_period_secs) {
                let now = start.elapsed().as_secs_f64();
                if period > 0.0 && last_log.is_some_and(|t| now - t >= period) {
                    let inc_obj = incumbent.as_ref().map(|(o, _)| *o);
                    let db = aggregate_bound(sense, stack.iter().map(|q| q.est));
                    emit_progress(h, &start, stats.nodes, inc_obj, db, sense);
                    last_log = Some(now);
                }
            }
            // ② mipgap：gap（经统一 clamp ≥ 0）≤ P ⇒ 粘滞置位并跳出，记下
            // 触发时刻证据 (db, inc)——与并行证据口径一致，不依赖后续栈态。
            if let Some(p) = ctrl.mipgap_percent {
                if p > 0.0 {
                    let inc_obj = incumbent.as_ref().map(|(o, _)| *o);
                    let db = aggregate_bound(sense, stack.iter().map(|q| q.est));
                    if let (Some(inc), Some(db)) = (inc_obj, db) {
                        if let Some(g) = gap_percent(sense, inc, db) {
                            if g <= p {
                                gap_hit = true;
                                gap_evidence = Some(db);
                                break;
                            }
                        }
                    }
                }
            }
            // ③ 策略级自动降级（仅 auto 模式且给了时间预算）：资格预检通过
            // 才聚合前沿界（避免逐节点 O(栈) 开销）；双停滞（incumbent 与界
            // 都不再前进）⇒ 与 ② 同径提前终止，证据口径一致。W3：budget
            // 模式的进度/预算/停滞锚点单位 = 节点数（progress_units 注入，
            // 零墙钟读点）。
            if let (Some(state), Some(budget)) = (adaptive.as_mut(), budget_progress_cap) {
                let now = progress_units(&start, stats.nodes, budget_mode);
                if AutoDegrade::eligible(state.degrade.cfg(), now, Some(budget)) {
                    let db = aggregate_bound(sense, stack.iter().map(|q| q.est));
                    let snap = DegradeSnapshot {
                        sense,
                        elapsed_secs: now,
                        budget_secs: Some(budget),
                        has_incumbent: incumbent.is_some(),
                        last_improvement_secs: state.last_improve_progress,
                        bound: db,
                    };
                    if state.degrade.decide(&snap) {
                        gap_hit = true;
                        degraded = true;
                        gap_evidence = db;
                        break;
                    }
                }
            }
        }
        // 栈空 ⇒ 搜索自然耗尽（先于上限判定，见模块实现说明）。
        let Some(p) = stack.pop() else {
            break;
        };
        // W6：退出局部割作用域——弹出的 Pending 属祖先层（unwind_to 更短）
        // ⇒ 其施加点在作用域栈顶的割全部出域，恢复快照（栈序 = 深度序）。
        while let Some(&(scope_trail, _, _)) = cut_scopes.last() {
            if scope_trail > p.unwind_to {
                let (_, snap, snap_iters) = cut_scopes.pop().expect("checked last");
                base_lp_iterations =
                    base_lp_iterations.saturating_add(simplex.iterations().saturating_sub(snap_iters));
                simplex = snap;
            } else {
                break;
            }
        }
        // est 弹出剪枝（仅整数目标启用，与 cutoff 同门控）：父节点 LP 目标
        // est 已不优于 incumbent 的有效剪枝界 ⇒ 该子树不可能再出严格更优解，
        // 直接丢弃——不 unwind/tighten、不 resolve、不计 stats.nodes、不耗
        // node_limit（est 是子树最优值的方向性乐观界，判定与 (d) 同口径）。
        if obj_integral {
            if let Some((obj, _)) = incumbent.as_ref() {
                if bound_no_better(sense, p.est, *obj, true, opts.prune_abs, opts.prune_rel) {
                    continue;
                }
            }
        }
        // 节点/时间上限：每次 resolve 前检查；0 = 不限。W3：budget 模式的
        // 时限语义 = 节点预算（确定性，零墙钟读点）；wall 模式照旧墙钟检查。
        if opts.node_limit != 0 && stats.nodes >= opts.node_limit {
            limit_hit = true;
            break;
        }
        if let Some(cap) = budget_node_cap {
            if stats.nodes >= cap {
                limit_hit = true;
                break;
            }
        } else if opts.time_limit_ms != 0 && time_exceeded(&start, opts.time_limit_ms) {
            limit_hit = true;
            break;
        }
        // 回退到父路径界状态后施加本分支界；空界/无收缩 ⇒ (f) 剪枝（不 resolve）。
        unwind(&mut trail, &mut bounds, p.unwind_to);
        if !tighten_child(&mut bounds, &mut trail, p.var, p.tighten_lb, p.bound) {
            continue;
        }
        let mut sol = simplex.resolve(&bounds);
        stats.nodes += 1;
        // 伪代价观测：子节点 LP 目标相对父 est 的方向化增益 ÷ 方向分数
        // 距离（非 Optimal 子节点不记——不可行=无穷增益语义留空，保守）。
        if opts.branch_rule == BranchRule::Pseudocost && sol.status == LpStatus::Optimal {
            let gain = match sense {
                Sense::Maximize => p.est - sol.objective,
                Sense::Minimize => sol.objective - p.est,
            };
            if gain.is_finite() && gain > 0.0 && p.frac_dir > tol.integrality {
                pc.record(p.var, p.tighten_lb, gain / p.frac_dir);
            }
        }

        // (a)/(b)/(c)：resolve 状态判定（顺序固定）。
        match classify_node(sol.status, false) {
            NodeAction::Prune => continue,
            NodeAction::PruneLimit => {
                limit_hit = true;
                continue;
            }
            // W2：节点内墙钟超限 ⇒ 与 PruneLimit 同处置（剪枝 + limit_hit
            // 粘滞），下一轮顶部的节点/时间上限检查立即停树 ⇒ 诚实
            // LimitReached。
            NodeAction::TruncateTime => {
                limit_hit = true;
                continue;
            }
            NodeAction::FailUnbounded => return Err(nonroot_unbounded_error()),
            NodeAction::TreeUnbounded => {
                return Err(MilpError::Internal(
                    "非根节点分类不可能返回根无界动作".to_string(),
                ));
            }
            NodeAction::Process => {}
        }

        // —— W6：树内局部割轮（配额制；仅顺序路径）——
        if local_cuts_on
            && p.depth >= 1
            && p.depth <= LOCAL_CUT_MAX_DEPTH
            && stats.nodes >= next_local_cut_node
            && local_cut_budget > 0
        {
            let quota = LOCAL_CUT_PER_NODE.min(local_cut_budget as usize);
            let batch = cuts::separate_local_cuts(
                &simplex,
                &bounds,
                problem,
                &int_vars,
                CutOptions::default(),
                &tol,
                quota,
            );
            if !batch.is_empty() {
                // 快照（施加前状态；子树出口恢复）。
                cut_scopes.push((trail.len(), simplex.snapshot(), simplex.iterations()));
                let mut applied: u64 = 0;
                let mut add_ok = true;
                for cut in &batch {
                    if simplex.add_row(&cut.terms, RelOp::Le, cut.rhs).is_err() {
                        add_ok = false;
                        break;
                    }
                    applied += 1;
                }
                let strengthened = if add_ok {
                    simplex.resolve(&bounds)
                } else {
                    // add_row 失败（防御）：视同分离失败，走下方防御回退。
                    LpSolution {
                        status: LpStatus::IterationLimit,
                        objective: f64::NAN,
                        x: Vec::new(),
                        duals: Vec::new(),
                    }
                };
                match strengthened.status {
                    LpStatus::Optimal => {
                        sol = strengthened;
                        stats.cuts += applied;
                        local_cut_budget -= applied;
                        next_local_cut_node = stats.nodes + LOCAL_CUT_NODE_INTERVAL;
                    }
                    LpStatus::Infeasible => {
                        // 局部割切空本子树 LP ⇒ 子树无整数可行解（割在子树内
                        // 有效）：节点关闭——恢复快照（弹出刚压的作用域）、
                        // 不分支、不计割。
                        let (_, snap, snap_iters) = cut_scopes.pop().expect("just pushed");
                        base_lp_iterations = base_lp_iterations
                            .saturating_add(simplex.iterations().saturating_sub(snap_iters));
                        simplex = snap;
                        continue;
                    }
                    _ => {
                        // IterationLimit/TimeLimit/Unbounded：防御性回退——
                        // 恢复快照，用原 sol 继续本节点判定（不冒充加强解）。
                        let (_, snap, snap_iters) = cut_scopes.pop().expect("just pushed");
                        base_lp_iterations = base_lp_iterations
                            .saturating_add(simplex.iterations().saturating_sub(snap_iters));
                        simplex = snap;
                    }
                }
            }
        }

        // (d)/(e)/(f) 与分支压栈。
        branch_or_close(
            problem,
            &int_vars,
            &mut bounds,
            &sol,
            tol,
            opts,
            obj_integral,
            &mut incumbent,
            &mut stack,
            &mut trail,
            &mut simplex,
            stats.nodes,
            progress_units(&start, stats.nodes, budget_mode),
            budget_progress_cap,
            &mut pc,
            adaptive.as_mut(),
            p.depth,
        );
    }

    // —— 终态判定（与上限诚实性整合，禁止字面分裂）：gap 达标与上限触发
    // 一律 LimitReached（诚实：未证明最优），粘滞标记优先于完备判定——
    let status = if limit_hit || gap_hit {
        SolveStatus::LimitReached
    } else if incumbent.is_some() {
        SolveStatus::Optimal
    } else {
        SolveStatus::Infeasible
    };
    // —— bound 口径：gap 终止 = 触发时刻证据 db；截断 = 剩余栈 est 聚合
    //（可 None）；完备有 incumbent = best_objective；Infeasible = None——
    let best_objective = incumbent.as_ref().map(|(o, _)| *o);
    let bound = if gap_hit {
        gap_evidence
    } else if limit_hit {
        aggregate_bound(sense, stack.iter().map(|q| q.est))
    } else if status == SolveStatus::Optimal {
        best_objective
    } else {
        None
    };
    Ok(finish(
        status,
        incumbent,
        &mut stats,
        &simplex,
        start,
        base_lp_iterations,
        bound,
        gap_hit,
        degraded,
    ))
}

// ===========================================================================
// 节点分类与判定（纯函数，便于单元级直接测试）
// ===========================================================================

/// 节点 resolve 状态 → 动作（模块内纯函数）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum NodeAction {
    /// (a) LP 不可行 ⇒ 剪枝回溯。
    Prune,
    /// (c) 迭代上限 ⇒ 剪枝并置 limit_hit 标记。
    PruneLimit,
    /// (c') 内部墙钟超限（W2，[`LpStatus::TimeLimit`]）⇒ 诚实截断：树节点
    /// 侧与 PruneLimit 同处置（剪枝 + limit_hit 粘滞，随后节点边界检查立即
    /// 停树）；**根侧刻意与 PruneLimit 分离**——割循环据此走「直接
    /// LimitReached 终态」而非 R8 回滚（回滚重建+重放是又一次整段根 LP，
    /// 会把超限点再拖过 deadline，见 cuts.rs 回滚语义注）。
    TruncateTime,
    /// (b) 根节点无界 ⇒ 整树终止返回 [`SolveStatus::Unbounded`]。
    TreeUnbounded,
    /// (b) 非根节点无界（父有界时理论不发生的数值异常）⇒ 整树终止
    /// `Err(MilpError::Numerical)`，不得静默继续。
    FailUnbounded,
    /// Optimal ⇒ 继续 (d)/(e)/(f) 判定。
    Process,
}

/// 节点 resolve 状态分类：(a) 不可行、(b) 无界（根/非根分流）、(c) 迭代
/// 上限 / (c') 内部墙钟超限（W2，独立动作，防与 R8 回滚混淆）。
fn classify_node(status: LpStatus, is_root: bool) -> NodeAction {
    match status {
        LpStatus::Infeasible => NodeAction::Prune,
        LpStatus::IterationLimit => NodeAction::PruneLimit,
        LpStatus::TimeLimit => NodeAction::TruncateTime,
        LpStatus::Unbounded => {
            if is_root {
                NodeAction::TreeUnbounded
            } else {
                NodeAction::FailUnbounded
            }
        }
        LpStatus::Optimal => NodeAction::Process,
    }
}

/// 非根节点 LP 无界的数值错误：子节点仅为父界的收紧，收紧不可能使有界 LP
/// 变无界；出现即数值异常，整树终止。
fn nonroot_unbounded_error() -> MilpError {
    MilpError::Numerical("node LP unbounded after bounded root".to_string())
}

/// 根节点动作 → 整树终态（`None` = 继续树搜索）。
fn root_terminal(action: NodeAction) -> Result<Option<SolveStatus>> {
    match action {
        NodeAction::Process => Ok(None),
        NodeAction::Prune => Ok(Some(SolveStatus::Infeasible)),
        NodeAction::PruneLimit => Ok(Some(SolveStatus::LimitReached)),
        // W2：墙钟超限 ⇒ 诚实 LimitReached 终态（无 R8 回滚，见
        // NodeAction::TruncateTime 文档）。
        NodeAction::TruncateTime => Ok(Some(SolveStatus::LimitReached)),
        NodeAction::TreeUnbounded => Ok(Some(SolveStatus::Unbounded)),
        NodeAction::FailUnbounded => Err(MilpError::Internal(
            "根节点分类不可能返回非根无界动作".to_string(),
        )),
    }
}

/// (e) 整数节点判定结果。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum NodeEval {
    /// 复核通过 ⇒ 调用方更新 incumbent（若严格更优）后回溯。
    Record,
    /// 非整数，或伪整数复核失败 ⇒ 继续分支。
    Branch,
}

/// (e) 纯判定：`integral`（整数性满足）且 `check_passed`（round_and_check
/// 复核通过）⇒ [`NodeEval::Record`]；否则 [`NodeEval::Branch`]。
///
/// 伪整数回退：复核不通过时不得剪枝、不得记 incumbent，必须继续分支——
/// 堵死“丢最优/可行问题误报 Infeasible”路径。
fn integer_node_eval(integral: bool, check_passed: bool) -> NodeEval {
    if integral && check_passed {
        NodeEval::Record
    } else {
        NodeEval::Branch
    }
}

/// (d) 有效剪枝界合成（pub(crate)：顺序/并行三调用点同口径的唯一实现，
/// 并行引擎原逐字副本已删除、改导入本函数）。
///
/// `integral = false`（[`crate::model::objective_is_integral`] 判否）⇒ 普通
/// 容差口径：Max = inc + prune_abs + prune_rel·|inc|、Min 镜像，与历史行为
/// 完全一致；`integral = true`（整数目标，目标值恒在整数网格）⇒ 与整数
/// cutoff 取更激进者：
/// - Max：`max(inc + 1 − ε, inc + prune_abs + prune_rel·|inc|)`，ε = 1e-6·(1+|inc|)；
/// - Min：`min(inc − 1 + ε, inc − prune_abs − prune_rel·|inc|)`。
/// 理由：整数 cutoff 是精确推理（Max 下任何严格更优解 obj ≥ inc+1），
/// prune_abs/prune_rel 是用户自担的启发放宽，二者取更激进者。
///
/// # ε 的工程余量论证（非形式证明）
/// incumbent 侧精确（整数系数 × 整数取值 + 整数 offset，|inc| < 2^53 内
/// f64 精确表示），风险仅在 bound 侧低估：表格漂移受 REFRESH_INTERVAL=200
/// 周期重算与终检 10·feasibility 夹回抑制，工作容差 1e-9 量级，ε 相对量级
/// 1e-6 高出 3 个数量级（适用域 |inc| ≲ 1e4 内低于整数网格 6 个数量级）。
/// 误剪需双条件合取：真实上确界落入宽 ε 窗口 [inc+1−ε, inc+1) 且计算误差
/// 同向超过窗口余量——工程上以高容差余量兑底，另有 bench DP 真值与验收
/// 钉死目标值兜底。
///
/// # 契约注记（|inc| ≥ 1e6 的保守退化）
/// |inc| ≥ 1e6 时 ε ≥ 1 ⇒ inc+1−ε ≤ inc，经 max()/min() 合成后退化为普通
/// prune 容差——有意保守（永不弱于现状、也不在该量级加强）；『ε 低于网格
/// 6 个数量级』论证的适用域为 |inc| ≲ 1e4（bench 全部算例满足）。
pub(crate) fn effective_prune_bound(
    sense: Sense,
    inc: f64,
    integral: bool,
    prune_abs: f64,
    prune_rel: f64,
) -> f64 {
    let slack = prune_abs + prune_rel * inc.abs();
    match sense {
        Sense::Maximize => {
            if integral {
                let eps = 1e-6 * (1.0 + inc.abs());
                (inc + 1.0 - eps).max(inc + slack)
            } else {
                inc + slack
            }
        }
        Sense::Minimize => {
            if integral {
                let eps = 1e-6 * (1.0 + inc.abs());
                (inc - 1.0 + eps).min(inc - slack)
            } else {
                inc - slack
            }
        }
    }
}

/// (d) 剪枝判定（按原始 sense 方向；pub(crate) 供并行引擎共享，杜绝分叉）：
/// Max: bound ≤ 有效剪枝界；Min: bound ≥ 有效剪枝界（[`effective_prune_bound`]）。
pub(crate) fn bound_no_better(
    sense: Sense,
    bound: f64,
    inc: f64,
    integral: bool,
    prune_abs: f64,
    prune_rel: f64,
) -> bool {
    let eff = effective_prune_bound(sense, inc, integral, prune_abs, prune_rel);
    match sense {
        Sense::Maximize => bound <= eff,
        Sense::Minimize => bound >= eff,
    }
}

/// 整数性判定（pub(crate) 供并行引擎共享）：全部 int_var 满足
/// |x−round(x)| ≤ integrality。
pub(crate) fn is_integral(int_vars: &[usize], x: &[f64], integrality: f64) -> bool {
    int_vars.iter().all(|&j| match x.get(j) {
        Some(v) => (v - v.round()).abs() <= integrality,
        None => false,
    })
}

/// MostFractional 分支变量选择：max min(f, 1−f)（f = x−⌊x⌋），并列取最小下标。
fn select_branch_var(int_vars: &[usize], x: &[f64]) -> Option<(usize, f64)> {
    let mut best: Option<(usize, f64, f64)> = None; // (var, 值, 得分)
    for &j in int_vars {
        let Some(&v) = x.get(j) else {
            continue;
        };
        if !v.is_finite() {
            continue;
        }
        let f = v - v.floor();
        let score = f.min(1.0 - f);
        let replace = match best {
            None => true,
            Some((_, _, s)) => score > s, // 严格大于 ⇒ 并列保留最小下标
        };
        if replace {
            best = Some((j, v, score));
        }
    }
    best.map(|(j, v, _)| (j, v))
}

// ===========================================================================
// 界栈（trail）与分支决策
// ===========================================================================

/// trail 记录：某次界收紧前的旧值（逆序回放恢复）。
struct TrailEntry {
    var: usize,
    is_lb: bool,
    old: f64,
}

/// 待应用的分支决策（stack 元素）。
struct Pending {
    var: usize,
    /// true = 上行分支（收紧 lb 至 ⌊x̄⌋+1）；false = 下行（收紧 ub 至 ⌊x̄⌋）。
    tighten_lb: bool,
    bound: f64,
    /// 应用本决策前 trail 应回到的长度（即父路径的 trail 状态）。
    unwind_to: usize,
    /// 压栈时父节点 LP 目标（该子树最优值的上界估计）。整数目标门控开启时
    /// 参与陈旧前沿弹出剪枝（est 不优于有效剪枝界 ⇒ 弹出即弃，见主循环）；
    /// 否则为纯数据，仅供受控口径（对偶界聚合
    /// [`crate::search_ctrl::aggregate_bound`] / gap 判定 / 进度行）读取；
    /// 伪代价分支也用它记录子节点增益。
    est: f64,
    /// 该方向距被违反界的分数距离（下行 = x_j−⌊x_j⌋；上行 = ⌈x_j⌉−x_j），
    /// 伪代价观测用。非伪代价规则下不读。
    frac_dir: f64,
    /// 节点深度（根 = 0；W6 局部割的深度门控）。
    depth: u32,
}

/// 回退 trail 至长度 `to`（严格逆序恢复旧界）。
fn unwind(trail: &mut Vec<TrailEntry>, bounds: &mut Bounds, to: usize) {
    while trail.len() > to {
        let Some(e) = trail.pop() else {
            break;
        };
        if e.is_lb {
            bounds.lb[e.var] = e.old;
        } else {
            bounds.ub[e.var] = e.old;
        }
    }
}

/// 施加子节点分支界：相对当前界严格收缩且新区间非空才应用并记账 trail，
/// 返回 true；空界（(f)）或无收缩（无进展）返回 false，调用方跳过该分支
/// （不 resolve、不计节点）。
fn tighten_child(
    bounds: &mut Bounds,
    trail: &mut Vec<TrailEntry>,
    var: usize,
    tighten_lb: bool,
    bound: f64,
) -> bool {
    if tighten_lb {
        if bound <= bounds.lb[var] || bound > bounds.ub[var] {
            return false;
        }
        trail.push(TrailEntry {
            var,
            is_lb: true,
            old: bounds.lb[var],
        });
        bounds.lb[var] = bound;
    } else {
        if bound >= bounds.ub[var] || bound < bounds.lb[var] {
            return false;
        }
        trail.push(TrailEntry {
            var,
            is_lb: false,
            old: bounds.ub[var],
        });
        bounds.ub[var] = bound;
    }
    true
}

// ===========================================================================
// Pseudocost 可靠性强分支冷启动（P1-4；pub(crate) 辅助供并行 worker 共用）
// ===========================================================================

/// 可靠性判定阈值 η：变量两方向均有 ≥η 条伪代价观测视为可靠，不再探测。
const PC_RELIABILITY_ETA: u32 = 4;
/// 探测深度门：深度超过该值的节点不再探测（冷启动收益集中在浅层）。
const PC_PROBE_MAX_DEPTH: usize = 8;
/// 根节点全量冷启动的候选数上限（R3 探测日程）：根预算 = 2×min(分数候选
/// 数, 该值) 次 resolve——强分支收益集中在浅层，根一次全量初始化即让
/// 全树统计可靠；上限防止大实例（数百分数变量）根探测白烧。
const PC_ROOT_PROBE_MAX_CANDIDATES: usize = 64;

/// P1-4 可靠性强分支冷启动：对 `n_obs < η` 的分数候选（按当前伪代价得分
/// 降序、并列最小下标——先探测规则认为最值得分支的变量）做 down+up 双向
/// 探测：收紧分支界 → resolve → 按主循环子节点同口径 `pc.record` → 恢复
/// 原界（探测不留任何状态）。
///
/// 预算与门控（**确定性，绝不读墙钟**——规避自适应/降级层的时间面）：
/// - **探测日程（R3）**：根节点（depth 0）全量冷启动——预算 =
///   2×min(分数候选数, [`PC_ROOT_PROBE_MAX_CANDIDATES`])（每候选 down+up
///   两次 resolve），一次初始化让全树统计可靠；深度 1..=8 维持每节点
///   预算 `opts.pc_reliability_probes`（默认 8）；深度 >
///   [`PC_PROBE_MAX_DEPTH`] 不探测；
/// - 一个候选耗 down+up 两次 resolve（方向对称），剩余预算 < 2 即停；
/// - 0 = 纯被动（调用方已按 branch_rule = Pseudocost 门控）；
/// - 探测 resolve 不计 nodes（记账性，同 dive 回取口径），其迭代数如实
///   累计在 Simplex 统计里。
///
/// 正确性论证：探测只读取候选子节点的 LP 目标并精确恢复原界——不剪枝、
/// 不更新 incumbent、不写 trail（界修改直接改写并逐条恢复，trail 长度
/// 不变）⇒ 本节点的 (d)/(e) 判定与分支决策仍在原 LP 解 `sol` 上进行，
/// 两分支域划分完备性不受影响，搜索完整性不变。非 Optimal 探测结果不记
/// （与主循环子节点被动观测同口径：不可行 = 无穷增益语义留空，保守）。
pub(crate) fn probe_reliability(
    sense: Sense,
    bounds: &mut Bounds,
    int_vars: &[usize],
    x: &[f64],
    est: f64,
    tol: Tolerances,
    opts: &BbOptions,
    pc: &mut Pseudocosts,
    simplex: &mut dyn LpEngine,
    depth: usize,
) {
    if opts.pc_reliability_probes == 0 || depth > PC_PROBE_MAX_DEPTH {
        return;
    }
    // 分数候选（与 Pseudocosts::select 同口径过滤）中的低观测者。
    let mut candidates: Vec<(usize, f64, f64)> = Vec::new(); // (var, LP 值, 得分)
    for &j in int_vars {
        let Some(&v) = x.get(j) else {
            continue;
        };
        if !v.is_finite() {
            continue;
        }
        let f = v - v.floor();
        if f <= tol.integrality || f >= 1.0 - tol.integrality {
            continue;
        }
        if pc.n_obs(j) >= PC_RELIABILITY_ETA {
            continue;
        }
        candidates.push((j, v, pc.score(j, f)));
    }
    candidates.sort_by(|a, b| {
        b.2.partial_cmp(&a.2)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.0.cmp(&b.0))
    });
    // 探测日程（R3）：根全量（候选数截到 PC_ROOT_PROBE_MAX_CANDIDATES），
    // 深度 1-8 维持每节点预算。
    let mut budget = if depth == 0 {
        2 * candidates.len().min(PC_ROOT_PROBE_MAX_CANDIDATES)
    } else {
        opts.pc_reliability_probes
    };
    for &(j, v, _) in &candidates {
        if budget < 2 {
            break; // 一个候选需 down+up 两次 resolve（方向对称成对探测）
        }
        let floor = v.floor();
        let frac = v - floor;
        // down：ub → ⌊x⌋（frac_dir = f）；up：lb → ⌊x⌋+1（frac_dir = 1−f）。
        budget -= probe_one(sense, bounds, simplex, j, false, floor, est, frac, pc, tol);
        budget -= probe_one(
            sense,
            bounds,
            simplex,
            j,
            true,
            floor + 1.0,
            est,
            1.0 - frac,
            pc,
            tol,
        );
    }
}

/// 单方向探测：收紧 → resolve →（Optimal 且严格正增益时）按方向记账 →
/// 恢复原界。返回实际发生的 resolve 次数（0 = 分支界相对当前界无收缩/
/// 区间空——与 [`tighten_child`] 同判定，跳过）。
fn probe_one(
    sense: Sense,
    bounds: &mut Bounds,
    simplex: &mut dyn LpEngine,
    j: usize,
    tighten_lb: bool,
    bound: f64,
    est: f64,
    frac_dir: f64,
    pc: &mut Pseudocosts,
    tol: Tolerances,
) -> usize {
    // 与 tighten_child 同判定的可应用性（严格收缩且新区间非空）。
    let applicable = if tighten_lb {
        bound > bounds.lb[j] && bound <= bounds.ub[j]
    } else {
        bound < bounds.ub[j] && bound >= bounds.lb[j]
    };
    if !applicable {
        return 0;
    }
    let old = if tighten_lb { bounds.lb[j] } else { bounds.ub[j] };
    if tighten_lb {
        bounds.lb[j] = bound;
    } else {
        bounds.ub[j] = bound;
    }
    let res = simplex.resolve(bounds);
    if res.status == LpStatus::Optimal {
        // 与主循环子节点被动观测完全同口径的方向化增益 ÷ 分数距离。
        let gain = match sense {
            Sense::Maximize => est - res.objective,
            Sense::Minimize => res.objective - est,
        };
        if gain.is_finite() && gain > 0.0 && frac_dir > tol.integrality {
            pc.record(j, tighten_lb, gain / frac_dir);
        }
    }
    // 精确恢复原界：探测不留状态（trail 未动，节点判定仍基于原 LP 解）。
    if tighten_lb {
        bounds.lb[j] = old;
    } else {
        bounds.ub[j] = old;
    }
    1
}

// ===========================================================================
// Optimal 节点处理与 incumbent 维护
// ===========================================================================

/// Optimal 节点的 (d)/(e) 判定与分支压栈：
/// - (d) 界不优于 incumbent 的有效剪枝界（整数目标门控 obj_integral 时含
///   整数 cutoff，见 [`effective_prune_bound`]）⇒ 节点关闭；
/// - (e) 整数性满足 ⇒ round_and_check 复核，通过 ⇒ 更新 incumbent（若严格
///   更优）后节点关闭；不通过（伪整数）⇒ 落入分支；
/// - （受控扩展，默认关闭）自适应启发式：(d) 未剪枝且本节点未按整数记录
///   （即落入分支路径）、选分支变量之前，控制器 decide = Run ⇒ 在本节点
///   补跑 dive（是否产生新全局 incumbent 作为 hit record）并用真实界
///   resolve 回取（记账性，不计节点）；
/// - （P1-4）Pseudocost 规则下分支选择前做可靠性冷启动探测
///   （[`probe_reliability`]，预算由 `pc_reliability_probes` 门控）；
/// - 分支：按规则选变量（MostFractional / Pseudocost），下行（x ≤ ⌊x̄⌋）
///   后压 ⇒ 先探索（P1-5 best-estimate 压栈序仅在两方向均有伪代价观测
///   且上行估计严格更优时翻转，见压栈处注释）。
#[allow(clippy::too_many_arguments)]
fn branch_or_close(
    problem: &Problem,
    int_vars: &[usize],
    bounds: &mut Bounds,
    sol: &LpSolution,
    tol: Tolerances,
    opts: &BbOptions,
    obj_integral: bool,
    incumbent: &mut Option<(f64, Vec<f64>)>,
    stack: &mut Vec<Pending>,
    trail: &mut Vec<TrailEntry>,
    simplex: &mut dyn LpEngine,
    nodes: u64,
    // 决策通道进度量（wall = 墙钟秒 / budget = 节点数，W3 单位契约）。
    elapsed_secs: f64,
    // 决策通道预算量（单位与 elapsed_secs 一致；None = 不限）。
    budget_progress_cap: Option<f64>,
    pc: &mut Pseudocosts,
    mut adaptive: Option<&mut AdaptiveState>,
    depth: u32,
) {
    // (d)：LP bound 不优于 incumbent ⇒ 剪枝（按原始 sense 方向比较）。
    if let Some((obj, _)) = incumbent.as_ref() {
        if bound_no_better(
            problem.sense,
            sol.objective,
            *obj,
            obj_integral,
            opts.prune_abs,
            opts.prune_rel,
        ) {
            return;
        }
    }
    // (e)：整数性满足 ⇒ round_and_check 复核（目标按取整后的值评估）。
    let integral = is_integral(int_vars, &sol.x, tol.integrality);
    let mut check_passed = false;
    let mut improved = false;
    if integral {
        if let Some(inc) = heuristics::round_and_check(problem, bounds, &sol.x, &tol) {
            improved = maybe_update(problem.sense, incumbent, inc);
            check_passed = true;
        }
    }
    if integer_node_eval(integral, check_passed) != NodeEval::Branch {
        // 复核通过 ⇒ 回溯（节点关闭）；incumbent 改进刷新自适应锚点。
        if let Some(state) = adaptive {
            if improved {
                state.last_improvement_node = nodes;
                state.last_improve_progress = elapsed_secs;
            }
        }
        return;
    }
    // 自适应启发式（仅分支路径、选分支变量之前；adaptive=None ⇒ 零行为）。
    if let Some(state) = adaptive.as_deref_mut() {
        let snap = SearchSnapshot {
            nodes,
            incumbent: incumbent.as_ref().map(|(o, _)| *o),
            nodes_since_improvement: nodes.saturating_sub(state.last_improvement_node),
            // W3 单位契约：budget 模式下 elapsed/预算均为节点数（调用方注入）。
            time_budget_secs: budget_progress_cap,
            elapsed_secs,
        };
        if state.controller.decide(&snap) == HeuristicDecision::Run {
            let hit = match heuristics::dive(
                problem,
                simplex,
                bounds,
                &sol.x,
                opts.dive_max_fixes,
                &tol,
            ) {
                Some(inc) => maybe_update(problem.sense, incumbent, inc),
                None => false,
            };
            state.controller.record(hit);
            if hit {
                state.last_improvement_node = nodes;
                state.last_improve_progress = elapsed_secs;
            }
            // dive 扰动了 Simplex 内部基：用本节点真实界 resolve 回取
            //（记账性，不计节点；分支决策仍基于传入的本节点 LP 解 sol）。
            let _ = simplex.resolve(bounds);
        }
    }
    // P1-4 可靠性强分支冷启动（仅 Pseudocost 且预算 > 0 时有实际动作）：
    // 分支选择前对低观测候选做 down+up 双向探测，为本地伪代价统计供应
    // 可靠初值。探测只读 LP 界并精确恢复（见 probe_reliability 正确性
    // 论证）——不剪枝、不记节点、不动 incumbent，搜索完整性不变。
    if opts.branch_rule == BranchRule::Pseudocost {
        probe_reliability(
            problem.sense,
            bounds,
            int_vars,
            &sol.x,
            sol.objective,
            tol,
            opts,
            pc,
            simplex,
            trail.len(),
        );
    }
    // 分支：按规则选变量（MostFractional / Pseudocost）；无可分支变量为
    // 内部不可达防御（不剪枝即丢弃节点）。
    let selected = match opts.branch_rule {
        BranchRule::MostFractional => select_branch_var(int_vars, &sol.x),
        BranchRule::Pseudocost => pc.select(int_vars, &sol.x, tol.integrality),
    };
    let Some((var, val)) = selected else {
        return;
    };
    let floor = val.floor();
    let frac = val - floor;
    let unwind_to = trail.len();
    // P1-5 best-estimate 子代压栈序：两方向均有伪代价观测（n_obs ≥ 1）且
    // 上行估计退化严格更小（est_gain = p 方向均值 × 分数距离；更小 ⇒ 该
    // 子代估计目标更优）时翻转压栈序——上行后压 ⇒ 先探索；任一方向无
    // 观测（p = 1 的中性默认与真实增益量级不可比，比较失真）或估计不
    // 严格更优 ⇒ 保持现序（上行先压、下行先探索）。压栈序不改变两分支
    // 域划分的完备性，DFS 全覆盖不变。
    let explore_up_first = pc.n_obs(var) > 0
        && pc.est_gain(var, true, 1.0 - frac) < pc.est_gain(var, false, frac);
    // est = 压栈时父节点 LP 目标；frac_dir = 该方向距界的分数距离
    // （伪代价观测用）。
    let up = Pending {
        depth: depth + 1,
        var,
        tighten_lb: true,
        bound: floor + 1.0,
        unwind_to,
        est: sol.objective,
        frac_dir: 1.0 - frac,
    };
    let down = Pending {
        depth: depth + 1,
        var,
        tighten_lb: false,
        bound: floor,
        unwind_to,
        est: sol.objective,
        frac_dir: frac,
    };
    if explore_up_first {
        stack.push(down);
        stack.push(up); // 后压 ⇒ 先弹出探索
    } else {
        stack.push(up); // 现序：上行先压（后弹出）、下行先探索
        stack.push(down);
    }
}

/// 顺序路径的自适应启发式状态：WP1 控制器 + incumbent 改进锚点（供
/// [`SearchSnapshot::nodes_since_improvement`] 与降级停滞判定使用）。
struct AdaptiveState {
    controller: HeuristicController,
    last_improvement_node: u64,
    /// 最近一次 incumbent 严格改进时的进度量（wall = 墙钟秒 / budget =
    /// 节点数，W3 单位契约——根启发式供应的解记根阶段末尾进度）。策略级
    /// 自动降级的双停滞信号之一。
    last_improve_progress: f64,
    /// 策略级自动降级判定器（窗口界基准有状态）。
    degrade: AutoDegrade,
}

impl AdaptiveState {
    fn new(cfg: HeuristicConfig, degrade_policy: DegradePolicy) -> Self {
        AdaptiveState {
            controller: HeuristicController::new(cfg),
            last_improvement_node: 0,
            last_improve_progress: 0.0,
            degrade: AutoDegrade::new(degrade_policy),
        }
    }
}

/// 严格更优才替换 incumbent（按原始 sense 方向；目标与取值均为 WP3 取整后
/// 口径，防报不可行值）。返回是否发生了替换（自适应控制器的 hit 语义）。
fn maybe_update(sense: Sense, incumbent: &mut Option<(f64, Vec<f64>)>, inc: Incumbent) -> bool {
    if !inc.objective.is_finite() {
        return false; // 防御：非有限目标不入 incumbent
    }
    let better = match incumbent.as_ref() {
        None => true,
        Some((obj, _)) => match sense {
            Sense::Maximize => inc.objective > *obj,
            Sense::Minimize => inc.objective < *obj,
        },
    };
    if better {
        *incumbent = Some((inc.objective, inc.values));
    }
    better
}

/// W6 树内局部割参数（配额制——分离成本 ≫ 节点 resolve 成本，配额与
/// 间隔小步标定；见 docs/parallel-profiling.md §15 的 A/B）。
/// 深度门：≤ 8（深节点 LP 已贴局部整数包，分离效率低）。
const LOCAL_CUT_MAX_DEPTH: u32 = 8;
/// 间隔：每 64 个 resolve 最多一次分离。
const LOCAL_CUT_NODE_INTERVAL: u64 = 64;
/// 每节点配额：≤ 2 条。
const LOCAL_CUT_PER_NODE: usize = 2;
/// 全局局部割预算（条）。
const LOCAL_CUT_GLOBAL_BUDGET: u64 = 2_000;

/// 时间上限检查（调用方保证 `limit_ms > 0`）。
fn time_exceeded(start: &Instant, limit_ms: u64) -> bool {
    start.elapsed().as_millis() >= limit_ms as u128
}

/// W2 内部墙钟时限合成（pub(crate) 供并行侧共用）：`time_limit_ms > 0` ⇒
/// 求解起点 + 预算的绝对 deadline（注入 Simplex，迭代循环每
/// [`crate::simplex`] TIME_CHECK_INTERVAL 次主元条件读钟，超限
/// `LpStatus::TimeLimit` 诚实截断）；`0 = 不限` ⇒ None（**零墙钟读点**，
/// 无时限路径 t1 确定性零扰动——条件读钟红线）。
pub(crate) fn bb_time_deadline(start: &Instant, opts: &BbOptions) -> Option<Instant> {
    (opts.time_limit_ms > 0)
        .then(|| *start + std::time::Duration::from_millis(opts.time_limit_ms))
}

/// W3 budget 模式的时限→节点预算换算率（节点/秒，**固定常量、机器无关**）。
/// 标定（ask 风险栏的两参照）：
/// - mknap-5x300 健康态 t1 吞吐 87,800 节点 / 12.8-13.7s ≈ 6.4-6.9k 节点/s
///   （本窗实测 6.7-6.9k、R8 窗 6.4-7.6k）——取 6,000 略保守：30s 名义
///   预算 = 180,000 节点 > 完备需求 87,800，且降级资格线（start_frac 0.5 ×
///   180k = 90k）也高于之 ⇒ budget-30s 恒完备最优 4762@87,800（恢复 R8
///   窗参照行为，消除 wall 模式的 optimal/auto_degraded 漂移）；
/// - marketsplit-4x30 参照：根闭合 1 节点即最优，任何预算下不受影响。
/// 语义边界：budget 模式以确定性优先于墙钟诚实——慢机上墙钟可超名义值
/// （换算率不随机器浮动，这正是逐字节可复现的来源）。
pub(crate) const BUDGET_NODES_PER_SEC: f64 = 6_000.0;

/// W3：budget 模式下 `time_limit_ms` 的等价节点预算（确定性换算；
/// `time_limit_ms = 0` = 不限 ⇒ None——降级与预算保护关闭，行为 ≡ 无时限）。
/// ceil + max(1)：任何正时限至少允许 1 个节点（根）。
fn budget_node_cap_of(opts: &BbOptions) -> u64 {
    debug_assert!(opts.time_limit_ms > 0);
    let secs = opts.time_limit_ms as f64 / 1000.0;
    ((secs * BUDGET_NODES_PER_SEC).ceil() as u64).max(1)
}

/// W3：决策通道的「进度量」——wall 模式 = 墙钟秒（现状）；budget 模式 =
/// 节点数（**不读钟**）。注入 SearchSnapshot/DegradeSnapshot 的进度字段与
/// 自适应锚点（单位契约见 search_ctrl 各快照文档）。
fn progress_units(start: &Instant, nodes: u64, budget_mode: bool) -> f64 {
    if budget_mode {
        nodes as f64
    } else {
        start.elapsed().as_secs_f64()
    }
}

/// 组装结果并补全统计（耗时与累计迭代 = `base_iters` + 当前 Simplex 迭代——
/// 割丢弃回退重建过 Simplex 时 base 计入重建前的根 LP + 割循环迭代）；
/// bound/gap_terminated 口径由调用方按终态类别给出（见 [`BbResult::bound`]
/// 字段说明）。
fn finish(
    status: SolveStatus,
    incumbent: Option<(f64, Vec<f64>)>,
    stats: &mut BbStats,
    simplex: &dyn LpEngine,
    start: Instant,
    base_iters: u64,
    bound: Option<f64>,
    gap_terminated: bool,
    auto_degraded: bool,
) -> BbResult {
    stats.lp_iterations = base_iters.saturating_add(simplex.iterations());
    stats.time_ms = u64::try_from(start.elapsed().as_millis()).unwrap_or(u64::MAX);
    let (best_objective, best_values) = match incumbent {
        Some((obj, values)) => (Some(obj), Some(values)),
        None => (None, None),
    };
    BbResult {
        status,
        best_objective,
        best_values,
        stats: *stats,
        bound,
        gap_terminated,
        auto_degraded,
    }
}

// ===========================================================================
// 受控口径的纯辅助（进度行与对偶界聚合；gap 统一走 search_ctrl::gap_percent）
// ===========================================================================

/// 首行进度的 bound 口径：栈非空 ⇒ 剩余聚合（每个开放子树的最优值不会优于
/// 其压栈时的父节点 LP 目标 est）；栈空且已有 incumbent ⇒ 搜索事实上完备，
/// bound = incumbent（gap=0）。
fn frontier_bound(sense: Sense, stack: &[Pending], incumbent: Option<f64>) -> Option<f64> {
    aggregate_bound(sense, stack.iter().map(|p| p.est)).or(incumbent)
}

/// 组装并发射一行进度（字段口径与并行引擎一致；gap 由展示层经
/// [`gap_percent`] 现算，缺任一侧或 inf 语义渲染为 inf）。
fn emit_progress(
    hooks: &SearchHooks,
    start: &Instant,
    nodes: u64,
    incumbent: Option<f64>,
    bound: Option<f64>,
    sense: Sense,
) {
    (hooks.progress)(&ProgressInfo {
        sense,
        elapsed_secs: start.elapsed().as_secs_f64(),
        nodes,
        incumbent,
        bound,
    });
}

// ===========================================================================
// 并行委托（threads ≥ 2；引擎本体见 crate::parallel_bb，本文件只做组装与映射）
// ===========================================================================

/// 小算例调度门探测预算（W1，P1 §5 M3 落地）：threads ≥ 2 时先用顺序路径
/// 以该节点预算探测，预算内出终态 ⇒ 门命中直返。512 的取值依据既有 t1
/// 节点数的两簇空档：**并行纯亏簇**全部 ≤ 360（milp-mixed-250i40c 53/85、
/// knapsack-500 360、knapsack-2000 38、knapsack-5000 37、setcover-100x150
/// 13-15、hard-knapsack-12 159），**并行有收益/大树簇**全部 ≥ 2,120
///（mknap-3x100 t1 2,120；mknap-5x300 87,800）——512 落在两簇之间。
const SMALL_TREE_PROBE_NODES: u64 = 512;

/// 调度门静态规模阈（LP 单位 = 变量数 × (行数 + 1)，与每节点 LP 成本同
/// 量级的代理）：门内最大 = lp-wide-200x300（60,200）与 setcover-100x150
///（15,150）；门外最小 = assign-40x40（129,600）与 setcover-300x600
///（180,600）——大 LP 类不进门，探测前缀的墙钟上限被钉在「512 个小 LP」
///（mknap-5x300 门未命中前缀实测 ≈ 0.2s / 30s 量级，见 W1 落档）。
const SMALL_TREE_LP_UNITS: usize = 65_536;

/// 调度门适用判定（纯函数，确定性：规模 + 节点预算，零墙钟读点；threads
/// ≥ 2 与 auto 决策链由调用方保证）。取证依据（P1 §3）：小树并行 = 纯
/// 开销——milp-mixed-250i40c t4 下 root 单线程段 ~42%、Σcv_wait 57.29%、
/// 3/4 worker 饥饿（某 worker 中位节点数 0）、lp_resolve 仅 23.50%；树
/// ≤ 85 节点填不满 2/4 线程。用户紧节点预算（≤ 探测预算）是对搜索形态的
/// 显式选择，不进门（保持既有并行口径）。
fn small_tree_gate_applies(problem: &Problem, opts: &BbOptions) -> bool {
    let units = problem.vars.len() * (problem.constraints.len() + 1);
    units <= SMALL_TREE_LP_UNITS
        && (opts.node_limit == 0 || opts.node_limit > SMALL_TREE_PROBE_NODES)
}

/// threads ≥ 2 的受控入口：**先过小算例调度门**（W1），再由 [`SearchControl`]
/// 组装 [`crate::parallel_bb::ParallelConfig`]（BbOptions 原样透传），委托
/// [`crate::parallel_bb::parallel_branch_and_bound`]，再按本文件统一口径把
/// [`crate::parallel_bb::ParallelOutcome`] 映射回 [`BbResult`]。
///
/// W1 按类选栈（`work_stack_auto`，多线程默认）在**本层**解析（工作栈是纯
/// 并行概念，threads=1 顺序路径不读该字段 ⇒ t1 行为逐字节不变）：auto ⇒
/// [`resolve_work_stack_mode`] 按类检测翻转（mknap 类 ⇒ global，其余 ⇒
/// local）；显式 local/global 逐字生效（覆盖优先）。
///
/// 调度门**只对 auto 决策链**生效：显式 `--work-stack local|global` 是用户
/// 点名并行引擎的栈形态（P2 式栈形态 A/B 的测量入口），不进门。门命中
/// （[`small_tree_gate_applies`] 且顺序探测在 [`SMALL_TREE_PROBE_NODES`]
/// 预算内出终态——Optimal/Infeasible/Unbounded，或 gap/降级/时限等与
/// threads=1 同口径的诚实截断）⇒ 直接返回：**threads ≥ 2 的输出与
/// threads=1 逐字节相同**（同一条确定性代码路径），免掉池协调/worker 冷
/// 启动开销；预算耗尽（LimitReached 且 nodes ≥ 预算且非 gap/降级）⇒ 树
/// 不小，丢弃探测、并行从零派发（既有行为 + 有界前缀浪费；进度行在门未
/// 命中时会重复发射首行——probe 与并行各一次，--log-secs 用户可见，纯
/// 外观）。
fn parallel_entry(
    problem: &Problem,
    opts: &BbOptions,
    ctrl: &SearchControl,
    tol: Tolerances,
    hooks: Option<&SearchHooks>,
) -> Result<BbResult> {
    // W3：budget 确定性模式仅支持 threads=1——并行侧无确定性契约（README
    // 口径：nodes/迭代/耗时/取值跨运行可变），节点预算语义也未在并行引擎
    // 定义；显式报错优于静默给出无契约的组合。
    if ctrl.determinism == Determinism::Budget {
        return Err(MilpError::Internal(
            "--determinism budget 仅支持 threads=1（t1 确定性契约：决策路径零墙钟读点、双跑逐字节相同）；并行侧无确定性契约".to_string(),
        ));
    }
    // —— W1 按类选栈解析（仅并行路径；t1 不读该字段）——
    let (opts, auto_stack) = resolve_work_stack_mode(problem, opts);
    let opts = &opts;
    // —— 小算例调度门（W1，仅 auto 决策链）：确定性判定 + 顺序探测 ——
    if auto_stack && small_tree_gate_applies(problem, opts) {
        let mut probe_opts = *opts;
        probe_opts.node_limit = SMALL_TREE_PROBE_NODES;
        let mut probe_ctrl = *ctrl;
        probe_ctrl.threads = 1;
        let probe = sequential_bb(problem, &probe_opts, &probe_ctrl, tol, hooks)?;
        let budget_exhausted = probe.status == SolveStatus::LimitReached
            && !probe.gap_terminated
            && !probe.auto_degraded
            && probe.stats.nodes >= SMALL_TREE_PROBE_NODES;
        if !budget_exhausted {
            return Ok(probe);
        }
        // 门未命中：树不小 ⇒ 探测作废（计数器不串），并行从零。
    }
    let config = crate::parallel_bb::ParallelConfig {
        threads: ctrl.threads,
        mipgap_percent: ctrl.mipgap_percent,
        adaptive_heuristics: ctrl.adaptive_heuristics,
        log_period_secs: ctrl.log_period_secs,
        degrade: ctrl.degrade,
        progress: hooks.map(|h| h.progress.clone()),
    };
    let outcome = crate::parallel_bb::parallel_branch_and_bound(problem, opts, &config, tol)?;
    Ok(map_parallel_outcome(outcome))
}

/// 并行终态 → [`BbResult`]（与顺序路径口径逐一对应）：
/// - 根阶段终局（root_terminal）⇒ 直返（bound 引擎已按口径给出：纯 LP
///   Optimal = Some(objective)，其余 None）；
/// - stop ∈ {NodeLimit, TimeLimit}（及其他截断类）⇒ LimitReached，bound =
///   引擎上报的剩余聚合（截断不可知时诚实 None）；
/// - stop = Gap ⇒ LimitReached 且 gap_terminated = true，bound = 触发时刻
///   证据快照（引擎以 CAS 胜者独占写入，绝不join后重聚合）；
/// - stop = None ⇒ incumbent? Optimal : Infeasible（粘滞 stop 优先于完备），
///   完备有 incumbent 时 bound = best_objective。
fn map_parallel_outcome(outcome: crate::parallel_bb::ParallelOutcome) -> BbResult {
    use crate::parallel_bb::ParallelStop;
    // 根阶段终局 ⇒ 直返。
    if let Some(status) = outcome.root_terminal {
        return BbResult {
            status,
            best_objective: outcome.best_objective,
            best_values: outcome.best_values,
            stats: outcome.stats,
            bound: outcome.bound,
            gap_terminated: false,
            auto_degraded: false,
        };
    }
    // 粘滞 stop 优先于完备判定。NodeLimit/TimeLimit（含迭代上限类截断，
    // 引擎侧统一映射为 NodeLimit 停止码）⇒ 诚实 LimitReached。
    let (status, gap_terminated, auto_degraded, bound) = match outcome.stop {
        Some(ParallelStop::Gap) => (SolveStatus::LimitReached, true, false, outcome.bound),
        Some(ParallelStop::AutoDegrade) => (SolveStatus::LimitReached, true, true, outcome.bound),
        Some(ParallelStop::NodeLimit | ParallelStop::TimeLimit) => {
            (SolveStatus::LimitReached, false, false, outcome.bound)
        }
        None => {
            if outcome.best_objective.is_some() {
                (SolveStatus::Optimal, false, false, outcome.bound)
            } else {
                (SolveStatus::Infeasible, false, false, None)
            }
        }
    };
    BbResult {
        status,
        best_objective: outcome.best_objective,
        best_values: outcome.best_values,
        stats: outcome.stats,
        bound,
        gap_terminated,
        auto_degraded,
    }
}

// ===========================================================================
// 单元测试
// ===========================================================================

#[cfg(test)]
mod tests {
    use crate::simplex::Simplex;
    use super::*;
    use crate::model::{RelOp, VarType};

    /// 默认容差。
    fn tol() -> Tolerances {
        Tolerances::default()
    }

    // ===================== R3：类自适应默认 ======================

    /// （R3）结构检测：全 Binary + 全正系数 Le 行 + 2 ≤ 行数 ≤ 64 命中
    /// （多维背包类）；连续列 / 等式行 / Ge 行 / 混合符号系数 / 单行 /
    /// 行数超限均回落。（P3 重钉：单行背包向量从命中改为排除——行数下限
    /// 1→2 收紧，单行类 Pseudocost 零收益纯付探测开销，见
    /// detect_pseudocost_favorable 文档注释的 A/B 数字；新增 2 行命中向量。）
    #[test]
    fn detect_pseudocost_favorable_classification() {
        // 命中：2 Binary + 2 条正系数 Le 行（多维背包）。
        let mut p = Problem::new("k", Sense::Maximize);
        p.add_var("a", VarType::Binary).unwrap();
        p.add_var("b", VarType::Binary).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 7.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 5.0), (1, 3.0)], RelOp::Le, 7.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 2.0), (1, 4.0)], RelOp::Le, 6.0)
            .unwrap();
        assert!(detect_pseudocost_favorable(&p));
        // 排除（P3）：单行背包（行数下限 2 收紧——零收益纯付探测开销）。
        let mut single = p.clone();
        single.constraints.truncate(1);
        assert!(!detect_pseudocost_favorable(&single));
        // 排除：Integer 变量（非全 Binary）。
        let mut q = Problem::new("int", Sense::Maximize);
        q.add_var("a", VarType::Integer).unwrap();
        q.set_var_bounds(0, 0.0, 5.0).unwrap();
        q.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        q.add_constraint(None, vec![(0, 2.0)], RelOp::Le, 7.0).unwrap();
        q.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 5.0).unwrap();
        assert!(!detect_pseudocost_favorable(&q));
        // 排除：连续列。
        let mut c = p.clone();
        c.vars.push(crate::model::Variable {
            name: "t".into(),
            var_type: VarType::Continuous,
            lb: 0.0,
            ub: 5.0,
        });
        assert!(!detect_pseudocost_favorable(&c));
        // 排除：Ge 行（setcover 类）。
        let mut g = p.clone();
        g.constraints[0].op = RelOp::Ge;
        assert!(!detect_pseudocost_favorable(&g));
        // 排除：Eq 行（marketsplit 类）。
        let mut e = p.clone();
        e.constraints[0].op = RelOp::Eq;
        assert!(!detect_pseudocost_favorable(&e));
        // 排除：负系数。
        let mut n = p.clone();
        n.constraints[0].expr.terms[0].1 = -5.0;
        assert!(!detect_pseudocost_favorable(&n));
        // 排除：行数 > 64。
        let mut many = Problem::new("many", Sense::Maximize);
        many.vars.push(crate::model::Variable {
            name: "a".into(),
            var_type: VarType::Binary,
            lb: 0.0,
            ub: 1.0,
        });
        many.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        for i in 0..65 {
            many.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 1.0 + i as f64)
                .unwrap();
        }
        assert!(!detect_pseudocost_favorable(&many));
        let mut few = many.clone();
        few.constraints.truncate(64);
        assert!(detect_pseudocost_favorable(&few));
    }

    /// （R3）入口解析：auto=true 命中类 ⇒ Pseudocost + probes ≥ 8；未命中 ⇒
    /// MostFractional；显式（auto=false）原样透传（覆盖优先）。
    ///（P3 重钉：命中载体从单行背包改为 2 行——行数下限收紧后单行不再
    /// 命中，见 detect_pseudocost_favorable 文档注释。W1：工作栈解析移至
    /// parallel_entry（纯并行概念），本函数只管分支规则。）
    #[test]
    fn resolve_branch_rule_defaults_paths() {
        let mut p = Problem::new("k", Sense::Maximize);
        p.add_var("a", VarType::Binary).unwrap();
        p.add_var("b", VarType::Binary).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 7.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 5.0), (1, 3.0)], RelOp::Le, 7.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 2.0), (1, 4.0)], RelOp::Le, 6.0)
            .unwrap();
        let mut q = p.clone();
        q.constraints[0].op = RelOp::Eq; // 非命中类
        // auto + 命中 ⇒ Pseudocost + probes 抬到 ≥8。
        let r = resolve_branch_rule_defaults(&p, &BbOptions::default());
        assert_eq!(r.branch_rule, BranchRule::Pseudocost);
        assert!(!r.branch_rule_auto); // 解析后幂等
        assert!(r.pc_reliability_probes >= 8);
        // auto + 未命中 ⇒ MostFractional。
        let r2 = resolve_branch_rule_defaults(&q, &BbOptions::default());
        assert_eq!(r2.branch_rule, BranchRule::MostFractional);
        assert!(!r2.branch_rule_auto);
        // 显式覆盖优先：pseudo+probes0 原样透传（纯被动不被抬升）。
        let mut ex = BbOptions::default();
        ex.branch_rule = BranchRule::Pseudocost;
        ex.pc_reliability_probes = 0;
        ex.branch_rule_auto = false;
        let r3 = resolve_branch_rule_defaults(&p, &ex);
        assert_eq!(r3.branch_rule, BranchRule::Pseudocost);
        assert_eq!(r3.pc_reliability_probes, 0); // 显式 0 不被抬升
        // 显式 mostfrac 在命中类上仍 mostfrac。
        let mut ex2 = BbOptions::default();
        ex2.branch_rule = BranchRule::MostFractional;
        ex2.branch_rule_auto = false;
        let r4 = resolve_branch_rule_defaults(&p, &ex2);
        assert_eq!(r4.branch_rule, BranchRule::MostFractional);
    }

    /// （W1）按类选栈解析：auto ⇒ 命中类（mknap 族）global、未命中类
    ///（单行/等式行/Ge/连续列）local；显式 local/global 逐字透传
    ///（覆盖优先）；解析后 auto 位清零（幂等）且返回 was_auto=true 供
    /// 调度门判定；显式时 was_auto=false（门不生效，保持栈形态可测量）。
    #[test]
    fn resolve_work_stack_mode_paths() {
        let mut p = Problem::new("k", Sense::Maximize);
        p.add_var("a", VarType::Binary).unwrap();
        p.add_var("b", VarType::Binary).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 7.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 5.0), (1, 3.0)], RelOp::Le, 7.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 2.0), (1, 4.0)], RelOp::Le, 6.0)
            .unwrap();
        let mut q = p.clone();
        q.constraints[0].op = RelOp::Eq; // 非命中类（marketsplit 显式分支）
        // auto + 命中 ⇒ global，was_auto = true。
        let (r, was_auto) = resolve_work_stack_mode(&p, &BbOptions::default());
        assert!(r.global_work_stack, "命中类（mknap 族）⇒ global（W1）");
        assert!(!r.work_stack_auto); // 解析后幂等
        assert!(was_auto);
        // auto + 未命中 ⇒ local。
        let (r2, was_auto2) = resolve_work_stack_mode(&q, &BbOptions::default());
        assert!(!r2.global_work_stack, "未命中类保持 local（W1）");
        assert!(!r2.work_stack_auto);
        assert!(was_auto2);
        // 显式 local 在命中类上仍 local；显式 global 在未命中类上仍 global
        //（覆盖优先），was_auto = false（调度门不生效）。
        let mut ex = BbOptions::default();
        ex.work_stack_auto = false;
        ex.global_work_stack = false;
        let (r3, w3) = resolve_work_stack_mode(&p, &ex);
        assert!(!r3.global_work_stack, "显式 local 在命中类上仍 local");
        assert!(!w3);
        let mut ex2 = BbOptions::default();
        ex2.work_stack_auto = false;
        ex2.global_work_stack = true;
        let (r4, w4) = resolve_work_stack_mode(&q, &ex2);
        assert!(r4.global_work_stack, "显式 global 在未命中类上仍 global");
        assert!(!w4);
        // 幂等：已解析选项二次解析逐字段不变。
        let (r5, w5) = resolve_work_stack_mode(&p, &r);
        assert_eq!(r5.global_work_stack, r.global_work_stack);
        assert!(!r5.work_stack_auto && !w5);
    }

    /// （W2）内部墙钟时限分类：TimeLimit ⇒ 独立动作 TruncateTime（根/非根
    /// 一致），根侧诚实映射 LimitReached（不进 R8 回滚——PruneLimit 专属）；
    /// 迭代上限仍为 PruneLimit（回归钉：两状态码不得混淆）。
    #[test]
    fn classify_node_time_limit_is_distinct_from_iteration_limit() {
        assert_eq!(classify_node(LpStatus::TimeLimit, true), NodeAction::TruncateTime);
        assert_eq!(classify_node(LpStatus::TimeLimit, false), NodeAction::TruncateTime);
        assert_eq!(classify_node(LpStatus::IterationLimit, true), NodeAction::PruneLimit);
        assert_eq!(classify_node(LpStatus::IterationLimit, false), NodeAction::PruneLimit);
        assert_eq!(
            root_terminal(NodeAction::TruncateTime).unwrap(),
            Some(SolveStatus::LimitReached)
        );
        assert_eq!(
            root_terminal(NodeAction::PruneLimit).unwrap(),
            Some(SolveStatus::LimitReached)
        );
    }

    /// （W2）时限开启 ⇒ 诚实截断终态：setcover-100x150（根 LP 数百 ms）以
    /// 1ms 预算求解必在时限内截断——status = LimitReached、不得冒充 Optimal、
    /// 求解快速返回（不把整段根 LP 跑满——节点边界检查盲区已消除）。
    #[test]
    fn tiny_time_limit_truncates_root_lp_honestly() {
        let path = format!(
            "{}/../bench/instances/setcover-100x150.milp",
            env!("CARGO_MANIFEST_DIR")
        );
        let problem = crate::instance::load_instance(std::path::Path::new(&path))
            .unwrap_or_else(|e| panic!("算例加载失败：{e}"));
        let opts = BbOptions {
            time_limit_ms: 1,
            ..BbOptions::default()
        };
        let t0 = std::time::Instant::now();
        let r = branch_and_bound(&problem, &opts, tol()).unwrap();
        let elapsed = t0.elapsed();
        assert_eq!(
            r.status,
            SolveStatus::LimitReached,
            "1ms 预算下根 LP（数百 ms）必须诚实截断"
        );
        assert!(
            elapsed < std::time::Duration::from_secs(5),
            "截断必须及时（实测 {elapsed:?}）——盲区消除"
        );
    }

    /// （W2）纯 LP 捷径的内部时限：lp-wide-200x300（全连续，根 LP ~15ms）
    /// 以 1ms 预算 ⇒ LimitReached 且不携带 objective（纯 LP 无 incumbent
    /// 概念）；无预算 ⇒ Optimal 不受影响（对照）。
    #[test]
    fn pure_lp_time_limit_truncates_without_fake_objective() {
        let path = format!(
            "{}/../bench/instances/lp-wide-200x300.milp",
            env!("CARGO_MANIFEST_DIR")
        );
        let problem = crate::instance::load_instance(std::path::Path::new(&path))
            .unwrap_or_else(|e| panic!("算例加载失败：{e}"));
        // 对照：无预算 ⇒ Optimal（既有行为）。
        let full = branch_and_bound(&problem, &BbOptions::default(), tol()).unwrap();
        assert_eq!(full.status, SolveStatus::Optimal);
        assert!(full.best_objective.is_some());
        // 1ms 预算 ⇒ 诚实 LimitReached、无 objective。
        let opts = BbOptions {
            time_limit_ms: 1,
            ..BbOptions::default()
        };
        let r = branch_and_bound(&problem, &opts, tol()).unwrap();
        assert_eq!(r.status, SolveStatus::LimitReached);
        assert!(
            r.best_objective.is_none(),
            "纯 LP 截断不携带 best（无 incumbent 概念，不冒充）"
        );
    }

    // ===================== W3：确定性模式（budget） ======================

    /// （W3）budget 模式双跑逐字节：mknap-3x100 + 100ms 时限（换算 600 节点
    /// 预算 < t1 完备需求 2,120 ⇒ 必截断）+ 自适应启发式（auto——wall 模式
    /// 下的漂移源）。两次求解的**全部行为字段**逐字节相同（status/objective/
    /// values/nodes/lp_iterations/cuts/bound/gap_terminated/auto_degraded），
    /// 仅 time_ms 计时字段不同。
    #[test]
    fn budget_mode_double_run_byte_identical() {
        let path = format!(
            "{}/../bench/instances/mknap-3x100.milp",
            env!("CARGO_MANIFEST_DIR")
        );
        let problem = crate::instance::load_instance(std::path::Path::new(&path))
            .unwrap_or_else(|e| panic!("算例加载失败：{e}"));
        let opts = BbOptions {
            time_limit_ms: 100, // budget 换算：ceil(0.1 × 6000) = 600 节点 < 2,120
            ..BbOptions::default()
        };
        let ctrl = SearchControl {
            adaptive_heuristics: true,
            determinism: Determinism::Budget,
            ..SearchControl::default()
        };
        let a = branch_and_bound_controlled(&problem, &opts, &ctrl, tol(), None).unwrap();
        let b = branch_and_bound_controlled(&problem, &opts, &ctrl, tol(), None).unwrap();
        assert_eq!(a.status, b.status, "终态逐字节");
        assert_eq!(a.best_objective, b.best_objective, "目标逐字节");
        assert_eq!(a.best_values, b.best_values, "取值逐字节");
        assert_eq!(a.stats.nodes, b.stats.nodes, "节点数逐字节");
        assert_eq!(a.stats.lp_iterations, b.stats.lp_iterations);
        assert_eq!(a.stats.cuts, b.stats.cuts);
        assert_eq!(a.bound, b.bound);
        assert_eq!(a.gap_terminated, b.gap_terminated);
        assert_eq!(a.auto_degraded, b.auto_degraded);
        // 截断形态确认：600 节点预算（完备需求 2,120）⇒ 非完备终态。
        assert_eq!(
            a.status,
            SolveStatus::LimitReached,
            "预算 600 < 完备需求 2,120 ⇒ 截断（实际 nodes={}）",
            a.stats.nodes
        );
        assert!(a.stats.nodes <= 600, "节点预算不超限（实际 {}）", a.stats.nodes);
    }

    /// （W3）budget 节点预算的确定性截断点：time_limit_ms=100 ⇒
    /// ceil(0.1 × BUDGET_NODES_PER_SEC) = 600 节点——截断恰在预算处（pop 后
    /// resolve 前检查），双跑恒 600。
    #[test]
    fn budget_node_cap_conversion() {
        assert_eq!(BUDGET_NODES_PER_SEC, 6_000.0);
        let opts = BbOptions {
            time_limit_ms: 100,
            ..BbOptions::default()
        };
        assert_eq!(budget_node_cap_of(&opts), 600);
        // 1ms ⇒ ceil(6) = 6；极大值饱和不 panic。
        let tiny = BbOptions {
            time_limit_ms: 1,
            ..BbOptions::default()
        };
        assert_eq!(budget_node_cap_of(&tiny), 6);
        let huge = BbOptions {
            time_limit_ms: u64::MAX,
            ..BbOptions::default()
        };
        assert!(budget_node_cap_of(&huge) > 0, "饱和转换不 panic");
    }

    /// （W3）budget + threads ≥ 2 显式报错（并行无确定性契约；不做静默
    /// 降级语义）。
    #[test]
    fn budget_mode_rejects_parallel() {
        let p = Problem::new("tiny", Sense::Maximize);
        let r = branch_and_bound_controlled(
            &p,
            &BbOptions::default(),
            &SearchControl {
                threads: 2,
                determinism: Determinism::Budget,
                ..SearchControl::default()
            },
            tol(),
            None,
        );
        assert!(matches!(r, Err(MilpError::Internal(_))), "budget+threads=2 应报错");
    }

    /// （W3）budget 无时限 ≡ wall 无时限（无预算 ⇒ 降级/预算保护关闭，
    /// 决策通道零差异）：同算例两模式行为字段全等——同时钉住「wall 默认
    /// 零变化」的结构性对照。
    #[test]
    fn budget_no_limit_equals_wall_no_limit() {
        let path = format!(
            "{}/../bench/instances/mknap-3x100.milp",
            env!("CARGO_MANIFEST_DIR")
        );
        let problem = crate::instance::load_instance(std::path::Path::new(&path))
            .unwrap_or_else(|e| panic!("算例加载失败：{e}"));
        let wall = branch_and_bound_controlled(
            &problem,
            &BbOptions::default(),
            &SearchControl::default(),
            tol(),
            None,
        )
        .unwrap();
        let budget = branch_and_bound_controlled(
            &problem,
            &BbOptions::default(),
            &SearchControl {
                determinism: Determinism::Budget,
                ..SearchControl::default()
            },
            tol(),
            None,
        )
        .unwrap();
        assert_eq!(wall.status, SolveStatus::Optimal);
        assert_eq!(budget.status, wall.status);
        assert_eq!(budget.best_objective, wall.best_objective);
        assert_eq!(budget.stats.nodes, wall.stats.nodes);
        assert_eq!(budget.stats.lp_iterations, wall.stats.lp_iterations);
        assert_eq!(budget.stats.cuts, wall.stats.cuts);
        assert_eq!(budget.bound, wall.bound);
    }

    /// （W6）树内局部割：mknap-3x100（多维背包类，类门控命中）默认开启 ⇒
    /// 触发（cuts > 根割基线 57）、双跑全部行为字段逐字节确定、目标恒同
    /// 已知最优 2331；`local_cuts = false` ⇒ 与关闭路径逐字节一致
    /// （57 cuts / 2,120 节点 = W6 前基线，回归钉）。
    #[test]
    fn local_cuts_fire_deterministically_and_toggle_restores() {
        let path = format!(
            "{}/../bench/instances/mknap-3x100.milp",
            env!("CARGO_MANIFEST_DIR")
        );
        let problem = crate::instance::load_instance(std::path::Path::new(&path))
            .unwrap_or_else(|e| panic!("算例加载失败：{e}"));
        let on1 = branch_and_bound(&problem, &BbOptions::default(), tol()).unwrap();
        let on2 = branch_and_bound(&problem, &BbOptions::default(), tol()).unwrap();
        assert_eq!(on1.status, SolveStatus::Optimal);
        assert_close(on1.best_objective.unwrap(), 2331.0);
        // 确定性：行为字段双跑逐字节。
        assert_eq!(on1.stats.nodes, on2.stats.nodes);
        assert_eq!(on1.stats.lp_iterations, on2.stats.lp_iterations);
        assert_eq!(on1.stats.cuts, on2.stats.cuts);
        assert_eq!(on1.bound, on2.bound);
        // 触发：局部割计入 cuts（根割基线 57 → 更高）。
        assert!(
            on1.stats.cuts > 57,
            "局部割应触发（cuts = {}，根割基线 57）",
            on1.stats.cuts
        );
        // 开关：关闭后无局部割（cuts = 根割基线 57）且目标恒同；对照
        // adaptive-on 口径（探针/CLI 默认）与 pre-W6 基线 2,120 节点逐位
        // 一致——证明 off 路径纯净。
        let off = branch_and_bound(
            &problem,
            &BbOptions {
                local_cuts: false,
                ..BbOptions::default()
            },
            tol(),
        )
        .unwrap();
        assert_eq!(off.stats.cuts, 57, "关闭后只剩根割");
        assert_close(off.best_objective.unwrap(), 2331.0);
        let off_auto = branch_and_bound_controlled(
            &problem,
            &BbOptions {
                local_cuts: false,
                ..BbOptions::default()
            },
            &SearchControl {
                adaptive_heuristics: true,
                ..SearchControl::default()
            },
            tol(),
            None,
        )
        .unwrap();
        assert_eq!(off_auto.stats.cuts, 57);
        // raw（无 presolve）+ adaptive-on 的 pre-W6 基线 = 2,106 节点 /
        // 24,072 迭代（HEAD 实测钉；探针口径的 2,120 是 presolve 重建后
        // 的问题——与本测试的 raw 路径无关）。
        assert_eq!(off_auto.stats.nodes, 2106, "adaptive-on off-cut = pre-W6 raw 基线");
        assert_eq!(off_auto.stats.lp_iterations, 24072);
        assert_close(off_auto.best_objective.unwrap(), 2331.0);
    }

    /// 关闭根割与根启发式的选项：隔离分支定界本体行为
    /// （整数节点复核 round_and_check 属跨模块联调路径，见相关用例）。
    fn plain_opts() -> BbOptions {
        BbOptions {
            root_cuts: false,
            heuristics: false,
            ..BbOptions::default()
        }
    }

    fn assert_close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 1e-7 * (1.0 + b.abs()),
            "数值不符: {a} vs {b}"
        );
    }

    /// BbOptions 默认值逐字段符合契约。
    #[test]
    fn bb_options_default_contract() {
        let o = BbOptions::default();
        assert_eq!(o.branch_rule, BranchRule::MostFractional);
        assert_eq!(o.node_rule, NodeRule::DepthFirst);
        assert_eq!(o.node_limit, 2_000_000);
        assert_eq!(o.time_limit_ms, 0);
        assert!(o.root_cuts);
        assert!(o.heuristics);
        assert_eq!(o.dive_max_fixes, 4096);
        assert_eq!(o.prune_abs, 1e-9);
        assert_eq!(o.prune_rel, 1e-9);
        assert_eq!(o.pc_reliability_probes, 8);
        assert!(!o.global_work_stack);
        assert!(o.work_stack_auto, "W1：多线程默认按类选栈（auto）");
    }

    // ===================== W1：小算例调度门 ======================

    /// （W1）调度门适用判定：LP 单位阈（65,536 = 2^16）与节点预算边界。
    /// 边界向量按 bench 实测单位数取值（lp-wide-200x300 = 60,200 门内、
    /// assign-40x40 = 129,600 门外；node_limit = 0 不限或 > 512 进门，
    /// 紧预算 ≤ 512 不进）。
    #[test]
    fn small_tree_gate_applies_boundaries() {
        // 60,000 变量单位问题（200 变量 × 300 行 = 60,200 ≤ 65,536）进门。
        let mut small_lp = Problem::new("small-lp", Sense::Maximize);
        small_lp.vars.push(crate::model::Variable {
            name: "z".into(),
            var_type: VarType::Continuous,
            lb: 0.0,
            ub: 1.0,
        });
        small_lp.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        for i in 0..300 {
            small_lp
                .add_constraint(None, vec![(0, 1.0)], RelOp::Le, 1.0 + i as f64)
                .unwrap();
        }
        // 129,600 单位问题（1,600 变量 × 80 行）不进门。
        let mut big_lp = Problem::new("big-lp", Sense::Maximize);
        for j in 0..1600 {
            big_lp.vars.push(crate::model::Variable {
                name: format!("b{j}"),
                var_type: VarType::Binary,
                lb: 0.0,
                ub: 1.0,
            });
        }
        let terms: Vec<(usize, f64)> = (0..1600).map(|j| (j, 1.0)).collect();
        big_lp.set_objective(terms.clone(), 0.0).unwrap();
        for _ in 0..80 {
            big_lp.add_constraint(None, terms.clone(), RelOp::Le, 800.0).unwrap();
        }
        let unlimited = BbOptions::default();
        assert!(small_tree_gate_applies(&small_lp, &unlimited));
        assert!(!small_tree_gate_applies(&big_lp, &unlimited));
        // node_limit：0（不限）与 > 512 进门；512 与更紧预算不进门
        //（紧预算 = 用户对搜索形态的显式选择，保持既有并行口径）。
        for (limit, expected) in [(0u64, true), (513, true), (512, false), (100, false)] {
            let opts = BbOptions {
                node_limit: limit,
                ..BbOptions::default()
            };
            assert_eq!(
                small_tree_gate_applies(&small_lp, &opts),
                expected,
                "node_limit = {limit}"
            );
        }
    }

    /// （W1）门命中 ⇒ threads=4 输出与 threads=1 **逐字节相同**（同一顺序
    /// 代码路径；仅 time_ms 计时字段不同）：hard-knapsack-12 t1 = 159 节点
    ///（P3 口径）< 512 探测预算 ⇒ 门命中。行为字段全比：status/objective/
    /// values/nodes/lp_iterations/cuts/bound/gap_terminated/auto_degraded。
    #[test]
    fn small_tree_gate_hit_returns_t1_identical_result() {
        let weights = [
            7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0, 31.0, 37.0, 41.0, 43.0, 47.0,
        ];
        let values = [
            11.0, 20.0, 32.0, 44.0, 55.0, 70.0, 88.0, 99.0, 111.0, 122.0, 133.0, 144.0,
        ];
        let mut p = Problem::new("hard-knapsack-12", Sense::Maximize);
        for i in 0..12 {
            p.add_var(&format!("x{i}"), VarType::Binary).unwrap();
        }
        p.set_objective((0..12).map(|i| (i, values[i])).collect(), 0.0)
            .unwrap();
        p.add_constraint(
            None,
            (0..12).map(|i| (i, weights[i])).collect(),
            RelOp::Le,
            172.0,
        )
        .unwrap();
        let t1 = branch_and_bound_controlled(
            &p,
            &BbOptions::default(),
            &SearchControl::default(),
            tol(),
            None,
        )
        .unwrap();
        assert_eq!(t1.status, SolveStatus::Optimal);
        assert_eq!(t1.stats.nodes, 159, "P3 口径 t1 节点数（159）");
        let t4 = branch_and_bound_controlled(
            &p,
            &BbOptions::default(),
            &SearchControl {
                threads: 4,
                ..SearchControl::default()
            },
            tol(),
            None,
        )
        .unwrap();
        assert_eq!(t4.status, t1.status);
        assert_eq!(t4.best_objective, t1.best_objective);
        assert_eq!(t4.best_values, t1.best_values, "取值逐字节（含并列破键序）");
        assert_eq!(t4.stats.nodes, t1.stats.nodes, "节点数逐字节（门命中）");
        assert_eq!(t4.stats.lp_iterations, t1.stats.lp_iterations);
        assert_eq!(t4.stats.cuts, t1.stats.cuts);
        assert_eq!(t4.bound, t1.bound);
        assert_eq!(t4.gap_terminated, t1.gap_terminated);
        assert_eq!(t4.auto_degraded, t1.auto_degraded);
    }

    /// （W1）门未命中 ⇒ 并行路径照常（mknap-3x100：t1 = 2,120 节点 > 512
    /// 探测预算 ⇒ 探测耗尽、丢弃、并行从零）；终态完备且最优恒同（2331.0，
    /// bench .opt 独立真值口径）。3 rep 覆盖窃取调度漂移。
    #[test]
    fn small_tree_gate_miss_dispatches_parallel() {
        let path = format!(
            "{}/../bench/instances/mknap-3x100.milp",
            env!("CARGO_MANIFEST_DIR")
        );
        let problem = crate::instance::load_instance(std::path::Path::new(&path))
            .unwrap_or_else(|e| panic!("算例加载失败：{e}"));
        for rep in 0..3 {
            let out = branch_and_bound_controlled(
                &problem,
                &BbOptions::default(),
                &SearchControl {
                    threads: 4,
                    ..SearchControl::default()
                },
                tol(),
                None,
            )
            .unwrap_or_else(|e| panic!("rep {rep} 求解失败：{e}"));
            assert_eq!(out.status, SolveStatus::Optimal, "rep {rep}");
            assert_close(out.best_objective.unwrap(), 2331.0);
        }
    }

    /// 根预割启发式轮对割循环**基透明**（2026-10-09 O1 修复的回归钉）：
    /// 预割 dive 扰动 Simplex 基后，修复前割循环在备择最优基上分离 GMI，
    /// 割集偏离无预割路径（knapsack-2000 实测 12→4 割、594→1628 节点）；
    /// 修复后预割轮以冷 Simplex 重建收尾 ⇒ `heuristics` 开关不再改变
    /// 根割施加数（incumbent 可能更优、节点数只减不增，但割集恒冷基）。
    /// knapsack 类重退化 LP（多备择最优基）对该扰动最敏感，选作钉子。
    #[test]
    fn precut_round_is_basis_transparent_to_cut_loop() {
        // 11 变量单行背包（随机族搜索定位的分离算例：修复前预割 dive 扰动
        // 基使割集 5→4——dive 恰好落在备择最优基上、GMI 分离随之改变）。
        let obj = [22.0, 33.0, 29.0, 13.0, 18.0, 33.0, 10.0, 12.0, 15.0, 27.0, 9.0];
        let weights = [21.0, 41.0, 35.0, 25.0, 39.0, 11.0, 42.0, 28.0, 20.0, 36.0, 12.0];
        let cap = 145.07999999999998_f64;
        let mut p = crate::model::Problem::new("knapsack-degen", crate::model::Sense::Maximize);
        let vars: Vec<usize> = (0..obj.len())
            .map(|i| p.add_var(&format!("x{i}"), VarType::Binary).unwrap())
            .collect();
        p.set_objective(
            vars.iter().zip(obj.iter()).map(|(&j, &c)| (j, c)).collect(),
            0.0,
        )
        .unwrap();
        p.add_constraint(
            None,
            vars.iter().zip(weights.iter()).map(|(&j, &w)| (j, w)).collect(),
            RelOp::Le,
            cap,
        )
        .unwrap();
        let with_heur = branch_and_bound(
            &p,
            &BbOptions {
                root_cuts: true,
                heuristics: true,
                ..BbOptions::default()
            },
            tol(),
        )
        .unwrap();
        let without_heur = branch_and_bound(
            &p,
            &BbOptions {
                root_cuts: true,
                heuristics: false,
                ..BbOptions::default()
            },
            tol(),
        )
        .unwrap();
        assert_eq!(
            with_heur.stats.cuts, without_heur.stats.cuts,
            "预割轮不得改变根割施加数（基透明）"
        );
        assert_eq!(with_heur.status, crate::model::SolveStatus::Optimal);
        assert_eq!(without_heur.status, crate::model::SolveStatus::Optimal);
        assert_close(
            with_heur.best_objective.unwrap(),
            without_heur.best_objective.unwrap(),
        );
        // 预割 incumbent 只增不减 ⇒ 节点数不高于无预割路径。
        assert!(
            with_heur.stats.nodes <= without_heur.stats.nodes,
            "预割 incumbent 只应增强剪枝: {} > {}",
            with_heur.stats.nodes,
            without_heur.stats.nodes
        );
    }

    /// 节点分类纯函数：(Unbounded, is_root=false) → FailUnbounded（整树终止
    /// Err(Numerical) 分支）；(Unbounded, is_root=true) → TreeUnbounded
    /// （→ SolveStatus::Unbounded 分支）。非根 Unbounded 数学上不可用算例
    /// 构造（父有界时收紧界不可能无界），故只对纯函数断言。
    #[test]
    fn classify_node_unbounded_root_vs_nonroot() {
        assert_eq!(
            classify_node(LpStatus::Unbounded, true),
            NodeAction::TreeUnbounded
        );
        assert_eq!(
            classify_node(LpStatus::Unbounded, false),
            NodeAction::FailUnbounded
        );
        // FailUnbounded 动作到错误的映射。
        assert!(
            matches!(nonroot_unbounded_error(), MilpError::Numerical(msg) if msg.contains("unbounded"))
        );
        // 其余映射。
        assert_eq!(
            classify_node(LpStatus::Infeasible, false),
            NodeAction::Prune
        );
        assert_eq!(classify_node(LpStatus::Infeasible, true), NodeAction::Prune);
        assert_eq!(
            classify_node(LpStatus::IterationLimit, true),
            NodeAction::PruneLimit
        );
        assert_eq!(classify_node(LpStatus::Optimal, false), NodeAction::Process);
        // 根动作 → 终态。
        assert_eq!(
            root_terminal(NodeAction::TreeUnbounded).unwrap(),
            Some(SolveStatus::Unbounded)
        );
        assert_eq!(
            root_terminal(NodeAction::Prune).unwrap(),
            Some(SolveStatus::Infeasible)
        );
        assert_eq!(
            root_terminal(NodeAction::PruneLimit).unwrap(),
            Some(SolveStatus::LimitReached)
        );
        assert_eq!(root_terminal(NodeAction::Process).unwrap(), None);
    }

    /// 伪整数复核失败回退判定（纯函数）：integral 且复核不通过 ⇒ 必须分支
    /// （不得剪枝、不得记 incumbent）；复核通过 ⇒ Record；非整数 ⇒ Branch。
    #[test]
    fn pseudo_integer_fallback_branches() {
        assert_eq!(integer_node_eval(true, false), NodeEval::Branch);
        assert_eq!(integer_node_eval(true, true), NodeEval::Record);
        assert_eq!(integer_node_eval(false, false), NodeEval::Branch);
        assert_eq!(integer_node_eval(false, true), NodeEval::Branch);
    }

    /// (d) 剪枝比较按原始 sense 方向（Max/Min 镜像；integral=false 退化口径）。
    #[test]
    fn bound_prune_respects_sense() {
        let (pa, pr) = (1e-9, 1e-9);
        // Max：bound ≤ inc + slack ⇒ 剪枝。
        assert!(bound_no_better(Sense::Maximize, 10.0, 10.0, false, pa, pr));
        assert!(!bound_no_better(Sense::Maximize, 10.5, 10.0, false, pa, pr));
        // Min：bound ≥ inc − slack ⇒ 剪枝。
        assert!(bound_no_better(Sense::Minimize, 2.0, 2.0, false, pa, pr));
        assert!(!bound_no_better(Sense::Minimize, 1.5, 2.0, false, pa, pr));
    }

    /// effective_prune_bound 合成契约（cutoff 纯函数）：
    /// - integral=false ⇒ 普通容差（与历史口径逐值一致）；
    /// - integral=true ⇒ Max 取 max(inc+1−ε, inc+slack)、Min 取
    ///   min(inc−1+ε, inc−slack)，ε=1e-6·(1+|inc|)；
    /// - 用户容差更激进时（大 prune_abs）合成退让到用户值（取更激进者）；
    /// - |inc| ≥ 1e6 ⇒ ε ≥ 1 ⇒ inc+1−ε ≤ inc，合成保守退化为普通容差
    ///   （有意行为：永不弱于现状、也不在该量级加强）。
    #[test]
    fn effective_prune_bound_composition() {
        let (pa, pr) = (1e-9, 1e-9);
        // integral=false：普通容差逐值（slack = pa + pr·|inc|）。
        assert_eq!(effective_prune_bound(Sense::Maximize, 10.0, false, pa, pr), 10.0 + (1e-9 + 1e-9 * 10.0));
        assert_eq!(effective_prune_bound(Sense::Minimize, 10.0, false, pa, pr), 10.0 - (1e-9 + 1e-9 * 10.0));
        // integral=true, Max: inc=21 ⇒ ε=2.2e-5, cutoff=21+1−ε > inc+slack。
        let eff = effective_prune_bound(Sense::Maximize, 21.0, true, pa, pr);
        assert!((eff - (22.0 - 2.2e-5)).abs() < 1e-12, "实际 {eff}");
        // 边界窗口：est ∈ (inc+slack, inc+1−ε] 剪枝、est > inc+1−ε 保留
        //（inc=21 ⇒ ε=2.2e-5，cutoff=21.999978；探测值 22−5e-5=21.99995
        // 在窗口内，已数值核对 ≤ cutoff）。
        assert!(bound_no_better(Sense::Maximize, 21.5, 21.0, true, pa, pr), "21.5 ≤ 21+1−ε ⇒ 剪");
        assert!(bound_no_better(Sense::Maximize, 22.0 - 5.0e-5, 21.0, true, pa, pr), "窗口内 ⇒ 剪");
        assert!(
            !bound_no_better(Sense::Maximize, 22.0, 21.0, true, pa, pr),
            "22.0 = inc+1 恰在网格上 ⇒ 不剪（可能存在并列为 inc+1 的更优解）"
        );
        // 同窗口在 integral=false 下全部保留（机制性对照）。
        assert!(!bound_no_better(Sense::Maximize, 21.5, 21.0, false, pa, pr));
        assert!(!bound_no_better(Sense::Maximize, 22.0 - 5.0e-5, 21.0, false, pa, pr));
        // Min 镜像：inc=30 ⇒ ε=3.1e-5。
        let eff = effective_prune_bound(Sense::Minimize, 30.0, true, pa, pr);
        assert!((eff - (29.0 + 3.1e-5)).abs() < 1e-12, "实际 {eff}");
        assert!(bound_no_better(Sense::Minimize, 29.5, 30.0, true, pa, pr));
        assert!(!bound_no_better(Sense::Minimize, 29.0, 30.0, true, pa, pr));
        // 用户容差更激进（prune_abs=2）：Max 合成 = inc+2+rel（> cutoff）。
        let eff = effective_prune_bound(Sense::Maximize, 21.0, true, 2.0, pr);
        assert!((eff - 23.0).abs() < 1e-6, "取更激进者（用户容差胜出），实际 {eff}");
        // 大 incumbent 保守退化：inc=1e6 ⇒ ε≈1.000001 > 1 ⇒ cutoff ≤ inc，
        // 合成退化为 inc+slack（普通容差；括号对齐函数内的求和结合序）。
        let eff = effective_prune_bound(Sense::Maximize, 1.0e6, true, pa, pr);
        assert_eq!(eff, 1.0e6 + (1e-9 + 1e-9 * 1.0e6), "实际 {eff}");
        let eff = effective_prune_bound(Sense::Minimize, 1.0e6, true, pa, pr);
        assert_eq!(eff, 1.0e6 - (1e-9 + 1e-9 * 1.0e6), "实际 {eff}");
    }

    /// MostFractional 选择与并列最小下标破键。
    #[test]
    fn most_fractional_selection_and_tie_break() {
        let iv = [0usize, 1, 2];
        // 得分 0.5 / 0.2 / 0.2 → 变量 0。
        assert_eq!(select_branch_var(&iv, &[1.5, 2.2, 0.8]), Some((0, 1.5)));
        // 1.25/3.75 二进制精确表示：得分恰为 0.25 / 0.25 / 0 的真并列
        // （1.2/3.8 等十进制值的分数性受 fl 噪声影响并非精确并列）。
        assert_eq!(select_branch_var(&iv, &[1.25, 3.75, 0.0]), Some((0, 1.25)));
        // 得分 0 / 0.3 / 0 → 变量 1。
        assert_eq!(select_branch_var(&iv, &[2.0, 1.3, 4.0]), Some((1, 1.3)));
        // 负值：floor(−1.5) = −2，f = 0.5。
        assert_eq!(select_branch_var(&iv, &[-1.5, 0.0, 0.0]), Some((0, -1.5)));
        assert_eq!(select_branch_var(&[], &[1.5]), None);
    }

    /// 空整数集（纯 LP）走捷径：单节点、无割，直接按 LP 状态返回。
    #[test]
    fn empty_int_vars_pure_lp_shortcut() {
        // max 3x+2y s.t. x+y<=4, x+3y<=6, [0,3]^2 → (3,1)，目标 11。
        let mut p = Problem::new("lp", Sense::Maximize);
        p.add_var("x", VarType::Continuous).unwrap();
        p.add_var("y", VarType::Continuous).unwrap();
        p.set_var_bounds(0, 0.0, 3.0).unwrap();
        p.set_var_bounds(1, 0.0, 3.0).unwrap();
        p.set_objective(vec![(0, 3.0), (1, 2.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 4.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 3.0)], RelOp::Le, 6.0)
            .unwrap();
        let r = branch_and_bound(&p, &BbOptions::default(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_close(r.best_objective.unwrap(), 11.0);
        let v = r.best_values.unwrap();
        assert_eq!(v.len(), 2);
        assert_close(v[0], 3.0);
        assert_close(v[1], 1.0);
        assert_eq!(r.stats.nodes, 1);
        assert_eq!(r.stats.cuts, 0);
    }

    /// （M5）整数不可行问题在根割把 LP 松弛割空时：终态 Infeasible（曾误报
    /// Err(Numerical)）——CLI 按既定映射输出 `status = infeasible`、退出码 2。
    /// 算例见 cuts::tests::cut_loop_infeasible_after_cuts_is_infeasible_terminal
    ///（穷举 0 个整数可行点）；关割同问少量节点内同样诚实报 Infeasible。
    #[test]
    fn integer_infeasible_proved_by_cuts_reports_infeasible() {
        let mut p = Problem::new("cut-infeasible", Sense::Maximize);
        p.add_var("x1", VarType::Integer).unwrap();
        p.add_var("x2", VarType::Integer).unwrap();
        let c3 = p.add_var("c", VarType::Continuous).unwrap();
        p.set_var_bounds(0, 0.0, 5.0).unwrap();
        p.set_var_bounds(1, 0.0, 2.0).unwrap();
        p.set_var_bounds(2, 0.0, 3.0).unwrap();
        p.set_objective(vec![(0, 6.0), (1, 2.0), (2, 5.0)], 0.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 2.0), (1, 5.0), (2, 2.0)], RelOp::Eq, 8.0)
            .unwrap();
        p.add_constraint(
            None,
            vec![(0, 3.0), (1, 3.0), (2, 0.5)],
            RelOp::Le,
            3.5,
        )
        .unwrap();
        let _ = c3;
        // 开割（默认）：割空 ⇒ 根终局 Infeasible，不外泄数值错误。
        let r = branch_and_bound(&p, &BbOptions::default(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Infeasible);
        assert!(r.best_objective.is_none());
        assert!(r.stats.cuts >= 1, "该算例根 LP 分数处必产割");
        // 关割对照：树搜索同样收敛到 Infeasible（少量节点内）。
        let mut plain = BbOptions::default();
        plain.root_cuts = false;
        let r2 = branch_and_bound(&p, &plain, tol()).unwrap();
        assert_eq!(r2.status, SolveStatus::Infeasible);
    }

    /// 根 Infeasible 透传：不触割循环与启发式（stats.cuts == 0）。
    #[test]
    fn root_infeasible_passthrough_no_cuts() {
        let mut p = Problem::new("inf", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Ge, 2.5)
            .unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 2.0)
            .unwrap();
        let r = branch_and_bound(&p, &BbOptions::default(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Infeasible);
        assert!(r.best_objective.is_none());
        assert!(r.best_values.is_none());
        assert_eq!(r.stats.nodes, 1);
        assert_eq!(r.stats.cuts, 0, "根非 Optimal 不得进入割循环");
    }

    /// 根 Unbounded 透传（整数变量存在而根 LP 松弛无界 ⇒ 按契约近似直接判
    /// Unbounded）：不触割循环与启发式。
    #[test]
    fn root_unbounded_passthrough_no_cuts() {
        let mut p = Problem::new("unb", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        let r = branch_and_bound(&p, &BbOptions::default(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Unbounded);
        assert!(r.best_objective.is_none());
        assert!(r.best_values.is_none());
        assert_eq!(r.stats.nodes, 1);
        assert_eq!(r.stats.cuts, 0, "根非 Optimal 不得进入割循环");
    }

    /// 极小 node_limit 驱动 limit_hit ⇒ LimitReached（附 incumbent，若有）：
    /// 绝不把截断搜索谎报为 Optimal/Infeasible。
    #[test]
    fn tiny_node_limit_reports_limit_reached() {
        let make = || {
            let mut p = Problem::new("lim", Sense::Maximize);
            p.add_var("x", VarType::Integer).unwrap();
            p.add_var("y", VarType::Integer).unwrap();
            p.set_objective(vec![(0, 1.0), (1, 1.0)], 0.0).unwrap();
            p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 3.5)
                .unwrap();
            p
        };
        // node_limit = 1：仅根节点（LP 最优 x=3.5 分数）。
        let opts = BbOptions {
            node_limit: 1,
            ..plain_opts()
        };
        let r = branch_and_bound(&make(), &opts, tol()).unwrap();
        assert_eq!(r.status, SolveStatus::LimitReached);
        assert_eq!(r.stats.nodes, 1);
        assert!(
            r.best_objective.is_none(),
            "未取得 incumbent 时不带目标返回"
        );
        // node_limit = 2：根 + 下行子节点（x=3, y=0.5 仍分数）后截断。
        let opts = BbOptions {
            node_limit: 2,
            ..plain_opts()
        };
        let r = branch_and_bound(&make(), &opts, tol()).unwrap();
        assert_eq!(r.status, SolveStatus::LimitReached);
        assert_eq!(r.stats.nodes, 2);
    }

    /// 整数变量界夹在分数区间内（无可行整数点）：两个分支界均空/无收缩 ⇒
    /// 不压栈、不 resolve，Infeasible 且仅 1 节点。
    #[test]
    fn infeasible_integer_band_no_children() {
        let mut p = Problem::new("band", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.4, 0.6).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Infeasible);
        assert_eq!(r.stats.nodes, 1, "空界分支不得 resolve");
    }

    /// 兄弟节点 LP 不可行剪枝 + trail 回溯正确性：根 2.6 分数，下行 [0,2] 与
    /// 上行 [3,5] 的 LP 均不可行 ⇒ 全树 3 节点后 Infeasible。
    #[test]
    fn sibling_lp_infeasible_prune_and_backtrack() {
        let mut p = Problem::new("sib", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 5.0).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Ge, 2.5)
            .unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 2.6)
            .unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Infeasible);
        assert_eq!(r.stats.nodes, 3, "根 + 两个 LP 不可行子节点");
        assert!(r.best_values.is_none());
    }

    /// （跨模块联调：整数节点复核经 WP3 heuristics::round_and_check）
    /// 单变量 MILP 全路径：max x s.t. 2x≤7 → 最优 3。整数目标（系数 1 ×
    /// Integer）⇒ est 弹出剪枝启用：下行分支取 incumbent 3 后，上行分支
    /// est=3.5 ≤ 有效剪枝界 3+1−ε ⇒ 跳过不 resolve——2 节点（原 3 节点，
    /// 上行子节点曾实际 resolve 后按 (d) 剪枝）。
    #[test]
    fn small_milp_optimal_integration() {
        let mut p = Problem::new("m1", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 2.0)], RelOp::Le, 7.0)
            .unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_close(r.best_objective.unwrap(), 3.0);
        let v = r.best_values.unwrap();
        assert_eq!(v.len(), 1);
        assert_close(v[0], 3.0);
        assert_eq!(r.stats.nodes, 2);
    }

    /// （跨模块联调）Min 方向镜像：min x s.t. x≥2.5 → 最优 3，3 节点
    /// （下行子节点 LP 不可行，上行取得整数 incumbent；两颗子节点弹出时
    /// 尚无 incumbent ⇒ est 弹出剪枝不触发，节点数与机制引入前一致）。
    #[test]
    fn minimize_sense_integration() {
        let mut p = Problem::new("min1", Sense::Minimize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Ge, 2.5)
            .unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_close(r.best_objective.unwrap(), 3.0);
        assert_close(r.best_values.unwrap()[0], 3.0);
        assert_eq!(r.stats.nodes, 3);
    }

    /// （跨模块联调）两层分支 + incumbent 两次更新 + 深层回溯：
    /// max 10a+7b s.t. 5a+3b≤7, a,b∈{0,1} → LP 松弛 15（a=0.8,b=1），
    /// 整数最优 10（a=1,b=0），共 5 节点（整数目标 cutoff/est 剪枝均未
    /// 触发：各弹出节点 est（15/14.667）均高于有效剪枝界 11−ε）。
    #[test]
    fn binary_knapsack_two_level_integration() {
        let mut p = Problem::new("k2", Sense::Maximize);
        p.add_var("a", VarType::Binary).unwrap();
        p.add_var("b", VarType::Binary).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 7.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 5.0), (1, 3.0)], RelOp::Le, 7.0)
            .unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_close(r.best_objective.unwrap(), 10.0);
        let v = r.best_values.unwrap();
        assert_close(v[0], 1.0);
        assert_close(v[1], 0.0);
        assert_eq!(r.stats.nodes, 5);
    }

    /// （跨模块联调）行范数大（‖row‖₁=1e6）、整数偏差≈5e-7 的近整数根节点：
    /// 取整后行违反 0，复核应通过 ⇒ 根节点直接收敛（1 节点、目标 3）。
    #[test]
    fn near_integer_root_accept_integration() {
        let mut p = Problem::new("near1", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0e6)], RelOp::Le, 3_000_000.5)
            .unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_eq!(r.stats.nodes, 1, "近整数根节点经复核通过后不再分支");
        assert_close(r.best_objective.unwrap(), 3.0);
        assert_close(r.best_values.unwrap()[0], 3.0);
    }

    /// （跨模块联调）行范数大、取整造成行违反的伪整数边缘算例：
    /// LP 解 2.9999995（偏差≈5e-7 ≤ integrality），取整为 3 使行违反 0.5。
    /// - WP3 阶梯复核（violation ≤ 10feas + integrality×‖row‖₁ ≈ 1.0）通过
    ///   ⇒ 根节点记录 incumbent=3；
    /// - 复核若拒绝 ⇒ 伪整数回退分支，收敛到真整数最优 2。
    /// 两条路径都必须诚实给出 Optimal —— 关键回归点：复核失败时绝不能因
    /// 剪枝而误报 Infeasible。
    #[test]
    fn near_integer_rounding_edge_integration() {
        let mut p = Problem::new("near2", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0e6)], RelOp::Le, 2_999_999.5)
            .unwrap();
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(
            r.status,
            SolveStatus::Optimal,
            "伪整数复核失败必须回退分支，不得误报 Infeasible"
        );
        let obj = r.best_objective.unwrap();
        assert!(
            obj == 2.0 || obj == 3.0,
            "阶梯语义两种结局：回退分支得 2；复核通过得 3。实际 {obj}"
        );
    }

    // =====================================================================
    // 受控入口（SearchControl / SearchHooks / gap / 进度行 / 并行委托）
    // =====================================================================

    /// SearchControl 默认值逐字段符合契约（= 现版本行为，零变化）。
    #[test]
    fn search_control_default_contract() {
        let c = SearchControl::default();
        assert_eq!(c.threads, 1, "默认 1 = 纯顺序路径");
        assert_eq!(c.mipgap_percent, None);
        assert!(!c.adaptive_heuristics, "默认关闭（现行为：仅根启发式）");
        assert_eq!(c.log_period_secs, None);
    }

    // —— P1-4 可靠性强分支冷启动（probe_reliability 行为/边界）——

    /// probe_reliability：预算 0 ⇒ 零动作（纯被动等价——观测与界均不动）；
    /// 预算 8（默认）⇒ 低观测分数候选获得双向观测（增益口径 = 方向化增益 ÷
    /// 分数距离）且界被**精确恢复**（探测不留状态 ⇒ 搜索完整性不受影响）；
    /// 深度 > 8 ⇒ 不探测。
    #[test]
    fn probe_reliability_budget_gating_and_state_restoration() {
        // binary_knapsack 问题：根 LP a=0.8（分数）、b=1（整数），目标 15。
        let mut p = Problem::new("k2", Sense::Maximize);
        p.add_var("a", VarType::Binary).unwrap();
        p.add_var("b", VarType::Binary).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 7.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 5.0), (1, 3.0)], RelOp::Le, 7.0)
            .unwrap();
        let t = tol();
        let int_vars = p.integer_vars();
        let mut simplex = Simplex::new(&p, t).unwrap();
        let mut bounds = Bounds::from_problem(&p);
        let sol = simplex.solve(&bounds);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.objective - 15.0).abs() < 1e-9);
        let snapshot = bounds.clone();

        // 预算 0：零观测、界零扰动（纯被动等价的机制性保证）。
        let mut pc0 = Pseudocosts::new(2);
        let opts0 = BbOptions {
            pc_reliability_probes: 0,
            ..BbOptions::default()
        };
        probe_reliability(
            Sense::Maximize,
            &mut bounds,
            &int_vars,
            &sol.x,
            sol.objective,
            t,
            &opts0,
            &mut pc0,
            &mut simplex,
            0,
        );
        assert_eq!(pc0.n_obs(0), 0, "预算 0 ⇒ 不探测");
        assert_eq!(pc0.n_obs(1), 0);
        assert_eq!(bounds.lb, snapshot.lb, "预算 0 ⇒ 界零扰动");
        assert_eq!(bounds.ub, snapshot.ub);

        // 预算 8：a 获得双向观测——下行探测（a≤0 ⇒ obj 7）增益 (15−7)/0.8
        // = 10；上行探测（a≥1 ⇒ obj 44/3）增益 (1/3)/0.2 = 5/3。界精确恢复。
        let mut pc = Pseudocosts::new(2);
        let opts = BbOptions::default();
        probe_reliability(
            Sense::Maximize,
            &mut bounds,
            &int_vars,
            &sol.x,
            sol.objective,
            t,
            &opts,
            &mut pc,
            &mut simplex,
            0,
        );
        assert!(pc.n_obs(0) >= 1, "分数变量 a 双向探测后必有观测");
        assert!(
            (pc.est_gain(0, false, 1.0) - 10.0).abs() < 1e-9,
            "下行均值 = 10，实际 {}",
            pc.est_gain(0, false, 1.0)
        );
        assert!(
            (pc.est_gain(0, true, 1.0) - 5.0 / 3.0).abs() < 1e-9,
            "上行均值 = 5/3，实际 {}",
            pc.est_gain(0, true, 1.0)
        );
        assert_eq!(pc.n_obs(1), 0, "整数变量 b 非候选");
        assert_eq!(bounds.lb, snapshot.lb, "探测后界精确恢复（下行）");
        assert_eq!(bounds.ub, snapshot.ub, "探测后界精确恢复（上行）");

        // 深度门：depth > 8 ⇒ 不探测。
        let mut pc2 = Pseudocosts::new(2);
        probe_reliability(
            Sense::Maximize,
            &mut bounds,
            &int_vars,
            &sol.x,
            sol.objective,
            t,
            &opts,
            &mut pc2,
            &mut simplex,
            9,
        );
        assert_eq!(pc2.n_obs(0), 0, "深度 9 > 8 ⇒ 不探测");
    }

    /// 受控默认（SearchControl::default + hooks=None）与既有
    /// branch_and_bound 结果全等（含 nodes 计数与 bound/gap_terminated 口径）。
    #[test]
    fn controlled_default_matches_branch_and_bound() {
        // 纯 LP。
        let mut lp = Problem::new("lp", Sense::Maximize);
        lp.add_var("x", VarType::Continuous).unwrap();
        lp.add_var("y", VarType::Continuous).unwrap();
        lp.set_var_bounds(0, 0.0, 3.0).unwrap();
        lp.set_var_bounds(1, 0.0, 3.0).unwrap();
        lp.set_objective(vec![(0, 3.0), (1, 2.0)], 0.0).unwrap();
        lp.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 4.0)
            .unwrap();
        lp.add_constraint(None, vec![(0, 1.0), (1, 3.0)], RelOp::Le, 6.0)
            .unwrap();
        // 整数两变量 + 极小 node_limit（LimitReached 路径）。
        let mut lim = Problem::new("lim", Sense::Maximize);
        lim.add_var("x", VarType::Integer).unwrap();
        lim.add_var("y", VarType::Integer).unwrap();
        lim.set_objective(vec![(0, 1.0), (1, 1.0)], 0.0).unwrap();
        lim.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 3.5)
            .unwrap();
        // 根不可行。
        let mut inf = Problem::new("inf", Sense::Maximize);
        inf.add_var("x", VarType::Integer).unwrap();
        inf.set_var_bounds(0, 0.0, 10.0).unwrap();
        inf.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        inf.add_constraint(None, vec![(0, 1.0)], RelOp::Ge, 2.5).unwrap();
        inf.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 2.0).unwrap();

        let cases: Vec<Problem> = vec![lp, lim, inf];
        for p in &cases {
            for opts in [plain_opts(), BbOptions::default()] {
                let a = branch_and_bound(p, &opts, tol()).unwrap();
                let b = branch_and_bound_controlled(
                    p,
                    &opts,
                    &SearchControl::default(),
                    tol(),
                    None,
                )
                .unwrap();
                assert_eq!(a.status, b.status, "{:?} {:?}", p.name, opts.root_cuts);
                assert_eq!(a.best_objective, b.best_objective);
                assert_eq!(a.best_values, b.best_values);
                assert_eq!(a.stats.nodes, b.stats.nodes, "nodes 计数必须全等");
                assert_eq!(a.stats.lp_iterations, b.stats.lp_iterations);
                assert_eq!(a.stats.cuts, b.stats.cuts);
                assert_eq!(a.bound, b.bound);
                assert_eq!(a.gap_terminated, b.gap_terminated);
                assert!(!b.gap_terminated, "默认路径绝不可 gap 终止");
            }
        }
        // LimitReached 子句：node_limit=2 时两路径同为 LimitReached 且 nodes=2。
        let opts = BbOptions {
            node_limit: 2,
            ..plain_opts()
        };
        let a = branch_and_bound(&cases[1], &opts, tol()).unwrap();
        let b = branch_and_bound_controlled(
            &cases[1],
            &opts,
            &SearchControl::default(),
            tol(),
            None,
        )
        .unwrap();
        assert_eq!(a.status, SolveStatus::LimitReached);
        assert_eq!(a.stats.nodes, b.stats.nodes);
        assert_eq!(a.stats.nodes, 2);
    }

    /// dive 问题：max 10x1+x2+0.5 s.t. 20x1+3x2≤43，整数 [0,10]。
    /// 根 LP 22.0（x1=2.15 分数），根启发式 dive 得 (2,1)=21.5，整数最优
    /// 21.5。offset=0.5 有意非整数化：整数目标 cutoff/est 弹出剪枝的门控
    /// （objective_is_integral）判否 ⇒ 该 fixture 家族保住 gap/降级/
    /// node-limit 路径覆盖（offset=0.0 时根 bound 21.5 ≤ 21+1−ε 会在 (d)
    /// 直接闭合根节点，终态从 LimitReached/gap 触发翻转为完备 Optimal）。
    fn dive_problem() -> Problem {
        let mut p = Problem::new("dive-gap", Sense::Maximize);
        p.add_var("x1", VarType::Integer).unwrap();
        p.add_var("x2", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_var_bounds(1, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 1.0)], 0.5).unwrap();
        p.add_constraint(None, vec![(0, 20.0), (1, 3.0)], RelOp::Le, 43.0)
            .unwrap();
        p
    }

    /// dive 问题的受控选项：关根割（割可能改变根 LP 界，破坏用例确定性），
    /// 保留根启发式（dive 供应 incumbent=21.5）。
    fn dive_opts() -> BbOptions {
        BbOptions {
            root_cuts: false,
            heuristics: true,
            ..BbOptions::default()
        }
    }

    /// 顺序确定性 gap 用例（dive 问题，offset=0.5 非整数目标 ⇒ cutoff 关闭）：
    /// - mipgap=Some(3.0)：主循环首轮（栈=[根的两个子节点，est=22.0]）即
    ///   gap=(22.0−21.5)/21.5×100≈2.326 ≤ 3 ⇒ LimitReached + gap_terminated，
    ///   best=21.5，bound=触发时刻证据 22.0；
    /// - mipgap=Some(2.0)：2.326 > 2 不触发，全树 3 节点收敛 Optimal 21.5。
    #[test]
    fn sequential_gap_termination_dive_problem() {
        let p = dive_problem();
        // Some(3.0) ⇒ gap 终止。
        let ctrl = SearchControl {
            mipgap_percent: Some(3.0),
            ..SearchControl::default()
        };
        let r = branch_and_bound_controlled(&p, &dive_opts(), &ctrl, tol(), None).unwrap();
        assert_eq!(r.status, SolveStatus::LimitReached, "gap 终止诚实标注截断");
        assert!(r.gap_terminated);
        assert_close(r.best_objective.unwrap(), 21.5);
        assert_close(r.bound.expect("gap 终止必有证据 bound"), 22.0);
        assert_eq!(r.stats.nodes, 1, "首轮即触发：仅根节点");
        // gap_percent 以触发时刻证据重算：Some(≈2.33) 且非负有限。
        let g = crate::search_ctrl::gap_percent(
            Sense::Maximize,
            r.best_objective.unwrap(),
            r.bound.unwrap(),
        )
        .expect("gap_terminated ⇒ gap_percent 必为 Some");
        assert!(g.is_finite() && g >= 0.0);
        assert!((g - 2.3255813953488373).abs() < 1e-9, "实际 {g}");

        // Some(2.0) ⇒ 不触发，完备最优。
        let ctrl = SearchControl {
            mipgap_percent: Some(2.0),
            ..SearchControl::default()
        };
        let r = branch_and_bound_controlled(&p, &dive_opts(), &ctrl, tol(), None).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert!(!r.gap_terminated);
        assert_close(r.best_objective.unwrap(), 21.5);
        assert_close(r.bound.unwrap(), 21.5);
        assert_eq!(r.stats.nodes, 3);
    }

    /// gap 负值防护（N1）：gap_percent 一律 clamp ≥ 0——db ≤ inc（incumbent
    /// 跳跃越过陈旧 est，或 est=39.999999999999996 < inc=40.0 的浮点噪声）时
    /// gap = 0 而非负值（负 gap 会击穿 'gap=\d+\.\d+%' 正则与
    /// 『gap_reached 行永不矛盾』契约）。单元级直接调用覆盖为主。
    #[test]
    fn gap_percent_clamps_negative_to_zero() {
        use crate::search_ctrl::gap_percent;
        // Max：db ≤ inc ⇒ Some(0.0)。
        assert_eq!(
            gap_percent(Sense::Maximize, 40.0, 39.999999999999996),
            Some(0.0),
            "浮点噪声方向的 db<inc 必须 clamp 为 0"
        );
        assert_eq!(gap_percent(Sense::Maximize, 21.0, 21.0), Some(0.0));
        // 微小负 raw（约 -2.5e-12）clamp 后恰为 0（不会打印 -0.00）。
        assert_eq!(
            gap_percent(Sense::Maximize, 40.0, 40.0 - 1.0e-11),
            Some(0.0)
        );
        // 正常正值：Max=(bound−inc)/|inc|×100；Min 镜像=(inc−bound)/|inc|×100。
        let g = gap_percent(Sense::Maximize, 21.0, 21.5).unwrap();
        assert!((g - 50.0 / 21.0).abs() < 1e-9, "Max 方向 (21.5-21)/21×100，实际 {g}");
        let g = gap_percent(Sense::Minimize, 30.0, 27.0).unwrap();
        assert!((g - 10.0).abs() < 1e-9, "Min 镜像 (30-27)/30×100，实际 {g}");
        // Min 方向 db ≥ inc（无可改进）⇒ 0。
        assert_eq!(gap_percent(Sense::Minimize, 27.0, 30.0), Some(0.0));
        // 非有限输入（缺 incumbent/bound 的表达方式）⇒ None（inf 语义）。
        assert_eq!(gap_percent(Sense::Maximize, f64::NAN, 21.5), None);
        assert_eq!(gap_percent(Sense::Maximize, 21.0, f64::NAN), None);
    }

    /// API 级顺序侧 clamp 场景：大阈值 mipgap 触发时 gap_percent 必为非负
    /// 有限 Some（触发即有 incumbent + 证据 bound，统一 clamp 保证不 inf/负）。
    #[test]
    fn sequential_gap_trigger_yields_nonnegative_finite_gap() {
        // binary_knapsack 两层问题：root LP 15 分数。默认（2026-10 复核后
        // 回退）MostFractional 固定下行先探索：根下行（a≤0，LP 14）先探索
        // 记 incumbent=7 后其 est 14 留栈；上行子树继续至 b≤0 记 incumbent
        // 提升；触发时刻证据口径 ⇒ gap=46.67% ≤ 50 触发。（Pseudocost +
        // best-estimate 序下为 50%，见 --branch-rule pseudocost。R3 类自适应
        // 默认会把本例（全 Binary + 正系数 Le 单行）翻到 Pseudocost——本测试
        // 钉的是 gap 机制而非默认规则，固定 MostFractional 显式覆盖。）
        let mut p = Problem::new("k2", Sense::Maximize);
        p.add_var("a", VarType::Binary).unwrap();
        p.add_var("b", VarType::Binary).unwrap();
        p.set_objective(vec![(0, 10.0), (1, 7.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 5.0), (1, 3.0)], RelOp::Le, 7.0)
            .unwrap();
        let ctrl = SearchControl {
            mipgap_percent: Some(50.0),
            ..SearchControl::default()
        };
        let mut opts = plain_opts();
        opts.branch_rule = BranchRule::MostFractional;
        opts.branch_rule_auto = false; // 显式覆盖（R3）
        let r = branch_and_bound_controlled(&p, &opts, &ctrl, tol(), None).unwrap();
        assert_eq!(r.status, SolveStatus::LimitReached);
        assert!(r.gap_terminated);
        let g = crate::search_ctrl::gap_percent(
            Sense::Maximize,
            r.best_objective.unwrap(),
            r.bound.unwrap(),
        )
        .expect("gap_terminated ⇒ Some");
        assert!(g.is_finite() && g >= 0.0, "gap 必须非负有限，实际 {g}");
        assert!((g - 140.0 / 3.0).abs() < 1e-6, "实际 {g}");
    }

    /// mipgap=Some(0.0) ≡ None（结果全等；Some(0.0)=关闭）。
    #[test]
    fn controlled_mipgap_zero_equals_none() {
        let p = dive_problem();
        let none_r =
            branch_and_bound_controlled(&p, &dive_opts(), &SearchControl::default(), tol(), None)
                .unwrap();
        let zero_r = branch_and_bound_controlled(
            &p,
            &dive_opts(),
            &SearchControl {
                mipgap_percent: Some(0.0),
                ..SearchControl::default()
            },
            tol(),
            None,
        )
        .unwrap();
        assert_eq!(none_r.status, zero_r.status);
        assert_eq!(none_r.best_objective, zero_r.best_objective);
        assert_eq!(none_r.best_values, zero_r.best_values);
        assert_eq!(none_r.stats.nodes, zero_r.stats.nodes);
        assert_eq!(none_r.stats.lp_iterations, zero_r.stats.lp_iterations);
        assert_eq!(none_r.bound, zero_r.bound);
        assert!(!zero_r.gap_terminated, "Some(0.0) 不得触发 gap 终止");
    }

    // —— 整数目标 cutoff + 陈旧前沿 est 弹出剪枝（本包新机制）——

    /// 塌缩载体：1 大件（w=7,v=10）+ n 个相同小件（w=2,v=1），容量 12。
    /// 最优 = 12（大件 + 2 小件；2^（n+1） 穷举可证）；LP 松弛恒 12.5
    ///（小件分数 0.5）——旧口径下 down 链每层仍 12.5 分数、每个 up 兄弟
    /// 各自再展开子树；cutoff 生效后全部 est=12.5 ≤ 有效剪枝界
    /// 12+1−ε 的陈旧前沿在弹出时剪除（est-skip）。
    /// offset 参数：0.0 = 整数目标（cutoff/est-skip 启用）；0.5 = 非整数
    /// 孪生——objective_is_integral 判否 ⇒ 门控关闭 ⇒ 与机制引入前的代码
    /// 路径逐字段一致（全部系数/界为二进有理数，+0.5 等距平移在 f64 中
    /// 精确，保持全部比较语义）。
    fn chain_knapsack(n_smalls: usize, offset: f64) -> Problem {
        let mut p = Problem::new("chain-knap", Sense::Maximize);
        p.add_var("b", VarType::Binary).unwrap();
        for i in 0..n_smalls {
            p.add_var(&format!("s{i}"), VarType::Binary).unwrap();
        }
        let mut terms = vec![(0usize, 10.0)];
        terms.extend((0..n_smalls).map(|i| (i + 1, 1.0)));
        p.set_objective(terms, offset).unwrap();
        let mut weights = vec![(0usize, 7.0)];
        weights.extend((0..n_smalls).map(|i| (i + 1, 2.0)));
        p.add_constraint(None, weights, RelOp::Le, 12.0).unwrap();
        p
    }

    /// 塌缩载体的 2^（n+1） 穷举真值（与 [`chain_knapsack`] 同参）。
    fn chain_knapsack_bruteforce(n_smalls: usize, offset: f64) -> f64 {
        let mut best = f64::NEG_INFINITY;
        for mask in 0u32..(1u32 << (n_smalls + 1)) {
            let (mut ww, mut vv) = (0.0f64, offset);
            for j in 0..=n_smalls {
                if mask & (1 << j) != 0 {
                    ww += if j == 0 { 7.0 } else { 2.0 };
                    vv += if j == 0 { 10.0 } else { 1.0 };
                }
            }
            if ww <= 12.0 + 1e-9 {
                best = best.max(vv);
            }
        }
        best
    }

    /// cutoff 生效小背包例：节点数相对非整数孪生（= 旧口径）塌缩 ≥10×，
    /// objective 与 2^（n+1） 穷举真值精确一致（整数算术，无浮点误差）。
    #[test]
    fn integral_cutoff_collapses_tree_at_least_10x() {
        let n = 8;
        let opt = chain_knapsack_bruteforce(n, 0.0);
        assert_eq!(opt, 12.0, "穷举真值（大件 + 2 小件）");
        // 整数目标（cutoff + est-skip 启用）。
        let r = branch_and_bound(&chain_knapsack(n, 0.0), &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_eq!(r.best_objective.unwrap(), opt, "目标值与穷举精确一致");
        // 非整数孪生（门控关闭 = 机制引入前行为）。
        let r2 = branch_and_bound(&chain_knapsack(n, 0.5), &plain_opts(), tol()).unwrap();
        assert_eq!(r2.status, SolveStatus::Optimal);
        assert_eq!(r2.best_objective.unwrap(), opt + 0.5, "孪生目标 = 真值 + offset");
        assert_eq!(
            r.best_values, r2.best_values,
            "两孪生最优取值一致（等距平移不改Argmax集）"
        );
        // 实测（2026-10-08，plain_opts）：cutoff 路径 7 节点（根 + 6 层
        // down 链恰至剩余小件整数填充），旧口径孪生 167 节点（up 兄弟各自
        // 展开）——塌缩 23.9×；钉关系断言（≥10×）而非具体值，避免与并发
        // 工作包的单纯形顶点破键微调耦合。
        assert!(
            10 * r.stats.nodes <= r2.stats.nodes,
            "节点数塌缩须 ≥10×：cutoff={} 旧口径={}",
            r.stats.nodes,
            r2.stats.nodes
        );
    }

    /// est 弹出剪枝的节点记账：陈旧前沿子树不 resolve——同 node_limit 下
    /// 整数目标侧完备 Optimal（nodes 恰为无陈旧开销的需求数），非整数
    /// 孪生侧（旧口径）在陈旧前沿上烧尽预算 ⇒ LimitReached 且恰耗满。
    #[test]
    fn est_skip_stale_frontier_not_resolved_under_node_limit() {
        let n = 8;
        // plain_opts（无割/无启发式）下两者皆确定性，预算内自测规模。
        let need_new =
            branch_and_bound(&chain_knapsack(n, 0.0), &plain_opts(), tol()).unwrap().stats.nodes;
        let need_old =
            branch_and_bound(&chain_knapsack(n, 0.5), &plain_opts(), tol()).unwrap().stats.nodes;
        assert!(need_new < need_old, "前置：cutoff 生效（{need_new} < {need_old}）");
        let limit = need_old - 1;
        let opts = BbOptions {
            node_limit: limit,
            ..plain_opts()
        };
        // 整数目标：陈旧前沿（est=12.5 ≤ cutoff 12+1−ε）弹出即弃——不
        // resolve、不计 nodes、不耗 node_limit ⇒ 同预算完备。
        let r = branch_and_bound(&chain_knapsack(n, 0.0), &opts, tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal, "陈旧子树未 resolve ⇒ 预算够用");
        assert_eq!(r.stats.nodes, need_new, "nodes 计数不含被跳过的陈旧子树");
        assert_eq!(r.best_objective.unwrap(), 12.0);
        // 非整数孪生（= 旧口径）：同预算烧在陈旧前沿的 resolve 上。
        let r2 = branch_and_bound(&chain_knapsack(n, 0.5), &opts, tol()).unwrap();
        assert_eq!(r2.status, SolveStatus::LimitReached);
        assert_eq!(r2.stats.nodes, limit, "旧口径恰耗满预算");
    }

    /// 非整数目标两型（分数 offset / 连续变量入目标）：门控判否 ⇒ 搜索
    /// 行为与机制引入前逐字段一致（机制性保证 + 确定性钉值复验）。
    #[test]
    fn nonintegral_objectives_keep_search_unchanged() {
        // ① 分数 offset 型：max 4x s.t. 8x≤18，x 整数 [0,3]。旧口径语义：
        // 根 LP x=2.25；down（x≤2）记 incumbent 8、up（x≥3）resolve 后
        // 不可行剪枝——3 节点 Optimal 8.25（含 offset 0.25…此处取 0.5）。
        let mut p = Problem::new("frac-offset", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 3.0).unwrap();
        p.set_objective(vec![(0, 4.0)], 0.5).unwrap();
        p.add_constraint(None, vec![(0, 8.0)], RelOp::Le, 18.0).unwrap();
        assert!(!crate::model::objective_is_integral(&p), "分数 offset ⇒ 门控关闭");
        let r = branch_and_bound(&p, &plain_opts(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_close(r.best_objective.unwrap(), 8.5);
        assert_eq!(r.stats.nodes, 3, "根 + down 记录 + up 不可行剪枝（旧口径）");
        // 对照：同问题 offset=0.0（整数目标）⇒ 门控开启，但 up 兄弟
        // est=9 = inc+1 恰在整数网格点上（可能存在并列为 8+1 的更优解）
        // ——cutoff 窗口 [inc+1−ε, inc+1) 不含网格点本身 ⇒ 不剪，仍 3 节点
        //（窗口边界语义与 effective_prune_bound_composition 的网格点断言
        // 一致；est 严格落入窗口内的塌缩证明见 chain_knapsack 用例）。
        let mut q = Problem::new("frac-offset-int", Sense::Maximize);
        q.add_var("x", VarType::Integer).unwrap();
        q.set_var_bounds(0, 0.0, 3.0).unwrap();
        q.set_objective(vec![(0, 4.0)], 0.0).unwrap();
        q.add_constraint(None, vec![(0, 8.0)], RelOp::Le, 18.0).unwrap();
        assert!(crate::model::objective_is_integral(&q));
        let r2 = branch_and_bound(&q, &plain_opts(), tol()).unwrap();
        assert_eq!(r2.stats.nodes, 3, "est=inc+1 恰在网格 ⇒ 不剪（保守窗口语义）");
        assert_close(r2.best_objective.unwrap(), 8.0);
        // ② 连续变量入目标型：max 4x + y s.t. 8x + 3y ≤ 18，x 整数 [0,3]、
        // y 连续 [0,2]。门控判否（连续入目标）⇒ 旧口径钉值：根 LP x=2.25
        // 分数；down x=2 ⇒ y=2/3 记 incumbent 8+2/3；up x=3 ⇒ 24>18 不可行
        // 剪枝——3 节点 Optimal 26/3。
        let mut c = Problem::new("cont-obj", Sense::Maximize);
        c.add_var("x", VarType::Integer).unwrap();
        c.add_var("y", VarType::Continuous).unwrap();
        c.set_var_bounds(0, 0.0, 3.0).unwrap();
        c.set_var_bounds(1, 0.0, 2.0).unwrap();
        c.set_objective(vec![(0, 4.0), (1, 1.0)], 0.0).unwrap();
        c.add_constraint(None, vec![(0, 8.0), (1, 3.0)], RelOp::Le, 18.0)
            .unwrap();
        assert!(!crate::model::objective_is_integral(&c), "连续入目标 ⇒ 门控关闭");
        let r3 = branch_and_bound(&c, &plain_opts(), tol()).unwrap();
        assert_eq!(r3.status, SolveStatus::Optimal);
        assert_close(r3.best_objective.unwrap(), 8.0 + 2.0 / 3.0);
        assert_eq!(r3.stats.nodes, 3, "根 + down 记录 + up 不可行剪枝（旧口径）");
    }

    /// bound 字段口径：完备 Optimal ⇒ best_objective；Infeasible ⇒ None；
    /// 截断（node_limit）⇒ 剩余栈 est 聚合（此处 = 根 LP 22.0）。
    #[test]
    fn controlled_bound_accounting() {
        // 完备 Optimal。
        let p = dive_problem();
        let r =
            branch_and_bound_controlled(&p, &dive_opts(), &SearchControl::default(), tol(), None)
                .unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_eq!(r.bound, Some(21.5));
        // Infeasible ⇒ None。
        let mut q = Problem::new("band", Sense::Maximize);
        q.add_var("x", VarType::Integer).unwrap();
        q.set_var_bounds(0, 0.4, 0.6).unwrap();
        q.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        let r = branch_and_bound_controlled(&q, &plain_opts(), &SearchControl::default(), tol(), None)
            .unwrap();
        assert_eq!(r.status, SolveStatus::Infeasible);
        assert_eq!(r.bound, None);
        // 截断：node_limit=1 ⇒ 根后栈非空即断 ⇒ bound=剩余聚合=根 est 22.0
        //（offset=0.5 非整数目标 ⇒ est 弹出剪枝关闭，弹出即触发上限）。
        let opts = BbOptions {
            node_limit: 1,
            ..dive_opts()
        };
        let r = branch_and_bound_controlled(
            &p,
            &opts,
            &SearchControl::default(),
            tol(),
            None,
        )
        .unwrap();
        assert_eq!(r.status, SolveStatus::LimitReached);
        assert!(!r.gap_terminated);
        assert_close(r.bound.unwrap(), 22.0);
        assert_close(r.best_objective.unwrap(), 21.5);
    }

    // —— hard-knapsack-12（N3 合成算例：并发强度与进度行载体）——

    /// hard-knapsack-12：12 变量 0-1 背包（权重/价值/容量见常量）。
    fn hard_knapsack_12() -> Problem {
        const W: [f64; 12] =
            [7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0, 31.0, 37.0, 41.0, 43.0, 47.0];
        const V: [f64; 12] =
            [11.0, 20.0, 32.0, 44.0, 55.0, 70.0, 88.0, 99.0, 111.0, 122.0, 133.0, 144.0];
        let mut p = Problem::new("hard-knapsack-12", Sense::Maximize);
        for j in 0..12 {
            p.add_var(&format!("x{j}"), VarType::Binary).unwrap();
        }
        p.set_objective((0..12).map(|j| (j, V[j])).collect(), 0.0)
            .unwrap();
        p.add_constraint(None, (0..12).map(|j| (j, W[j])).collect(), RelOp::Le, 172.0)
            .unwrap();
        p
    }

    /// 2^12 穷举背包最优（真值基准）。
    fn brute_force_knapsack(w: &[f64], v: &[f64], cap: f64) -> f64 {
        let n = w.len();
        assert_eq!(v.len(), n);
        let mut best = 0.0;
        for mask in 0u32..(1u32 << n) {
            let (mut ww, mut vv) = (0.0, 0.0);
            for j in 0..n {
                if mask & (1 << j) != 0 {
                    ww += w[j];
                    vv += v[j];
                }
            }
            if ww <= cap + 1e-9 && vv > best {
                best = vv;
            }
        }
        best
    }

    /// 顺序默认路径与 2^12 穷举一致（最优 523）。
    #[test]
    fn hard_knapsack_12_sequential_matches_bruteforce() {
        const W: [f64; 12] =
            [7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0, 31.0, 37.0, 41.0, 43.0, 47.0];
        const V: [f64; 12] =
            [11.0, 20.0, 32.0, 44.0, 55.0, 70.0, 88.0, 99.0, 111.0, 122.0, 133.0, 144.0];
        assert_eq!(brute_force_knapsack(&W, &V, 172.0), 523.0, "穷举真值须为 523");
        let p = hard_knapsack_12();
        let r = branch_and_bound(&p, &BbOptions::default(), tol()).unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_close(r.best_objective.unwrap(), 523.0);
    }

    // —— 进度行（首行契约 + 周期行仅栈非空）——

    /// 进度行快照（从 ProgressInfo 提取的纯数据，避免依赖 ProgressInfo: Clone）。
    #[derive(Debug, Clone, Copy)]
    struct Line {
        nodes: u64,
        incumbent: Option<f64>,
        bound: Option<f64>,
        gap: Option<f64>,
    }

    /// 收集器钩子：回调把每行快照压入共享 Vec（回调须 Send+Sync）。gap 由
    /// ProgressInfo 的 incumbent/bound 经 search_ctrl::gap_percent 现算（与
    /// format_progress_line 同口径）。
    fn collector() -> (std::sync::Arc<std::sync::Mutex<Vec<Line>>>, SearchHooks) {
        use std::sync::{Arc, Mutex};
        let lines = Arc::new(Mutex::new(Vec::new()));
        let sink = Arc::clone(&lines);
        let hooks = SearchHooks {
            progress: Arc::new(move |i: &crate::search_ctrl::ProgressInfo| {
                if let Ok(mut v) = sink.lock() {
                    v.push(Line {
                        nodes: i.nodes,
                        incumbent: i.incumbent,
                        bound: i.bound,
                        gap: crate::search_ctrl::gap_percent(
                            i.sense,
                            i.incumbent.unwrap_or(f64::NAN),
                            i.bound.unwrap_or(f64::NAN),
                        ),
                    });
                }
            }),
        };
        (lines, hooks)
    }

    /// hard-knapsack-12 + 极小 log_period：首行在根阶段结束（nodes=1）；
    /// 此后周期行仅在栈非空时出现 ⇒ 全部行 bound 均 Some（无
    /// 'incumbent=X bound=none gap=inf' 矛盾行）；行数 ≥ 2（首行后仍有行）。
    #[test]
    fn progress_lines_first_after_root_and_never_bound_none() {
        let p = hard_knapsack_12();
        let (lines, hooks) = collector();
        let ctrl = SearchControl {
            log_period_secs: Some(1e-9),
            ..SearchControl::default()
        };
        let r = branch_and_bound_controlled(&p, &BbOptions::default(), &ctrl, tol(), Some(&hooks))
            .unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        assert_close(r.best_objective.unwrap(), 523.0);
        let lines = lines.lock().unwrap();
        assert!(!lines.is_empty(), "必须收到首行进度");
        assert_eq!(
            lines[0].nodes, 1,
            "首行在根阶段结束的汇合点发射（nodes=1）"
        );
        assert!(
            lines.len() >= 2,
            "极小周期下首行之后至少一行周期行，实际 {} 行",
            lines.len()
        );
        for (i, l) in lines.iter().enumerate() {
            assert!(l.bound.is_some(), "第 {i} 行（nodes={}）bound 不得为 none", l.nodes);
        }
        // nodes 单调不减；行内容不钉数值。
        for w in lines.windows(2) {
            assert!(w[1].nodes >= w[0].nodes, "nodes 单调不减");
        }
    }

    /// 纯 LP 算例也必须收到首行：恰好一行、nodes=1、bound=incumbent=LP 最优、
    /// gap=Some(0)（完备搜索，clamp 后无 -0.00）。
    #[test]
    fn progress_first_line_pure_lp() {
        let mut p = Problem::new("lp", Sense::Maximize);
        p.add_var("x", VarType::Continuous).unwrap();
        p.add_var("y", VarType::Continuous).unwrap();
        p.set_var_bounds(0, 0.0, 3.0).unwrap();
        p.set_var_bounds(1, 0.0, 3.0).unwrap();
        p.set_objective(vec![(0, 3.0), (1, 2.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 4.0)
            .unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 3.0)], RelOp::Le, 6.0)
            .unwrap();
        let (lines, hooks) = collector();
        let ctrl = SearchControl {
            log_period_secs: Some(1e-9),
            ..SearchControl::default()
        };
        let r = branch_and_bound_controlled(&p, &BbOptions::default(), &ctrl, tol(), Some(&hooks))
            .unwrap();
        assert_eq!(r.status, SolveStatus::Optimal);
        let lines = lines.lock().unwrap();
        assert_eq!(lines.len(), 1, "纯 LP 只有根，恰一行");
        assert_eq!(lines[0].nodes, 1);
        assert_close(lines[0].incumbent.unwrap(), 11.0);
        assert_close(lines[0].bound.unwrap(), 11.0);
        assert_eq!(lines[0].gap, Some(0.0));
    }

    /// 根不可行终局：仍先发射首行再返回（nodes=1；incumbent/bound 均 none）。
    #[test]
    fn progress_first_line_root_infeasible() {
        let mut p = Problem::new("inf", Sense::Maximize);
        p.add_var("x", VarType::Integer).unwrap();
        p.set_var_bounds(0, 0.0, 10.0).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Ge, 2.5).unwrap();
        p.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 2.0).unwrap();
        let (lines, hooks) = collector();
        let r = branch_and_bound_controlled(
            &p,
            &BbOptions::default(),
            &SearchControl::default(),
            tol(),
            Some(&hooks),
        )
        .unwrap();
        assert_eq!(r.status, SolveStatus::Infeasible);
        let lines = lines.lock().unwrap();
        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0].nodes, 1);
        assert_eq!(lines[0].incumbent, None);
        assert_eq!(lines[0].bound, None);
    }

    // —— 并行委托与自适应（对接 crate::parallel_bb / crate::search_ctrl）——

    /// 5 个原验收算例（已知最优）。
    fn acceptance_problems() -> Vec<(String, Problem)> {
        use std::path::Path;
        let cases = [
            ("lp-prod-3", 13.0),
            ("knapsack-10", 17.0),
            ("assign-4x4", 10.0),
            ("milp-wolsey-2", 40.0),
            ("milp-mixed-5", 29.0),
        ];
        cases
            .iter()
            .map(|&(name, _)| {
                let path = format!("{}/instances/{name}.milp", env!("CARGO_MANIFEST_DIR"));
                let p = crate::instance::load_instance(Path::new(&path))
                    .unwrap_or_else(|e| panic!("算例 {name} 加载失败：{e}"));
                (name.to_string(), p)
            })
            .collect()
    }

    /// threads=4 于 5 个原验收算例：最优值与默认（threads=1）完全一致
    ///（池互斥弹出每节点恰处理一次 + incumbent 单调 + 粘滞 stop ⇒ 不丢最优）。
    ///（W1 测试调整，断言不变：显式 node_limit=512 旁路小算例调度门——这些
    /// 算例 t1 节点 ≪ 512，W1 起默认会被调度门接走走顺序路径，本测试的意图
    /// 是受控入口的**并行**正确性；紧节点预算（≤ 门预算）按门契约保持既有
    /// 并行口径，512 对这些算例的完备需求仍充裕。）
    #[test]
    fn threads4_matches_default_on_acceptance_instances() {
        for (name, p) in acceptance_problems() {
            let seq = branch_and_bound(&p, &BbOptions::default(), tol())
                .unwrap_or_else(|e| panic!("{name} 顺序求解失败：{e}"));
            let par = branch_and_bound_controlled(
                &p,
                &BbOptions {
                    node_limit: 512,
                    ..BbOptions::default()
                },
                &SearchControl {
                    threads: 4,
                    ..SearchControl::default()
                },
                tol(),
                None,
            )
            .unwrap_or_else(|e| panic!("{name} 并行求解失败：{e}"));
            assert_eq!(par.status, SolveStatus::Optimal, "算例 {name}");
            assert!(
                seq.status.is_optimal() && par.best_objective == seq.best_objective,
                "算例 {name}：并行最优值必须与顺序一致（seq={:?} par={:?}）",
                seq.best_objective,
                par.best_objective
            );
        }
    }

    /// threads=4 于 hard-knapsack-12（最优 523）重复 ≥100 次：每次断言
    /// objective=523 且 Optimal（丢唤醒/死锁/竞态在重复下必现；不依赖计时）。
    ///（W1 测试调整，断言不变：显式 node_limit=512 旁路调度门——t1 = 159
    /// 节点 < 512 会被门接走走顺序路径，而本测试的意图是并行引擎压测；
    /// 512 ≥ 完备需求 159，预算不截断搜索。）
    #[test]
    fn threads4_hard_knapsack_12_repeated_100() {
        let p = hard_knapsack_12();
        let ctrl = SearchControl {
            threads: 4,
            ..SearchControl::default()
        };
        let opts = BbOptions {
            node_limit: 512,
            ..BbOptions::default()
        };
        for i in 0..100 {
            let r = branch_and_bound_controlled(&p, &opts, &ctrl, tol(), None)
                .unwrap_or_else(|e| panic!("第 {i} 次并行运行失败：{e}"));
            assert_eq!(r.status, SolveStatus::Optimal, "第 {i} 次运行终态");
            assert_close(r.best_objective.unwrap(), 523.0);
        }
    }

    /// adaptive_heuristics=true 于 5 个原验收算例：仍 Optimal 且目标不变。
    #[test]
    fn adaptive_heuristics_keeps_acceptance_optima() {
        for (name, p) in acceptance_problems() {
            let base = branch_and_bound(&p, &BbOptions::default(), tol())
                .unwrap_or_else(|e| panic!("{name} 默认求解失败：{e}"));
            let ctrl = SearchControl {
                adaptive_heuristics: true,
                ..SearchControl::default()
            };
            let r = branch_and_bound_controlled(&p, &BbOptions::default(), &ctrl, tol(), None)
                .unwrap_or_else(|e| panic!("{name} 自适应求解失败：{e}"));
            assert_eq!(r.status, SolveStatus::Optimal, "算例 {name}");
            assert_close(r.best_objective.unwrap(), base.best_objective.unwrap());
            assert!(!r.gap_terminated);
        }
    }
}
