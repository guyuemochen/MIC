//! 句柄数据核心：`ModelBox`（模型 + 参数 + 最近一次求解结果 + 中毒标记）。
//!
//! `create` 成功返回 `Box::into_raw` 的裸指针；FFI 侧一切操作经该指针
//! 重建 `&mut ModelBox`。句柄非线程安全（禁止多线程并发使用同一句柄）。

use milp::branch_bound::BbStats;
use milp::model::{Problem, Sense, SolveStatus};

/// 启发式模式（FFI `HeuristicMode` 参数的内部表示；默认 `On`）。
///
/// 与 CLI 默认（auto）的差异：FFI 默认对齐 `SearchControl::default` 的库语义
/// （adaptive_heuristics=false），避免 auto + TimeLimit 意外自动降级为
/// LIMIT_REACHED（见 README 参数表）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum HMode {
    /// on：heuristics=true 且不自适应（库默认行为）。
    On,
    /// off：heuristics=false。
    Off,
    /// auto：heuristics=true 且 adaptive_heuristics=true。
    Auto,
}

/// FFI 参数镜像（Python `setParam` 经 FFI setter 逐项写入）。
///
/// 关断编码（0 是规范关断值，任何后设值可覆盖先设值）：
/// `TimeLimit`/`MIPGap`/`LogPeriodSecs` 的 0.0 ⇒ `None`（不限/关）。
#[derive(Debug, Clone, Copy)]
pub(crate) struct Params {
    /// 时间上限秒；None = 不限（0.0 关断编码）。
    pub time_limit_secs: Option<f64>,
    /// 线程数（1..=1024）；默认 1 = 纯顺序路径。
    pub threads: usize,
    /// MIP gap 百分比阈值；None = 关闭（0.0 关断编码）。
    pub mipgap_percent: Option<f64>,
    /// 启发式模式；默认 On。
    pub heuristic_mode: HMode,
    /// 进度行周期秒；None = 不注入进度钩子（0.0 关断编码）。
    pub log_period_secs: Option<f64>,
}

impl Default for Params {
    /// 默认：无时限 / threads=1 / 无 gap / HeuristicMode=on / 无日志。
    fn default() -> Self {
        Params {
            time_limit_secs: None,
            threads: 1,
            mipgap_percent: None,
            heuristic_mode: HMode::On,
            log_period_secs: None,
        }
    }
}

/// 最近一次成功求解的存档（每次 solve 整体覆盖）。
///
/// `values.len()` 是解维数唯一锚点（`n_solved`）：结构性修改（add_var 等）
/// 后仍以存档维数为准，当前 `problem.vars.len()` 只用于下标越界判定。
#[derive(Debug, Clone)]
pub(crate) struct StoredSolve {
    /// 求解终态。
    pub status: SolveStatus,
    /// 目标值（原 sense + obj_offset 口径）；无可行解时 None。
    pub objective: Option<f64>,
    /// 对偶界证据；无证据时 None。
    pub bound: Option<f64>,
    /// gap 触发时刻的百分比；未触发时 None。
    pub gap_percent: Option<f64>,
    /// 变量取值（原始变量空间）；无可行解时 None。
    pub values: Option<Vec<f64>>,
    /// 分支定界统计。
    pub stats: BbStats,
    /// 预处理化简量。
    pub presolve_reductions: usize,
    /// 是否因 mipgap 阈值提前终止。
    pub gap_terminated: bool,
    /// 是否因策略级自动降级提前终止。
    pub auto_degraded: bool,
}

/// FFI 句柄本体（`void*` 指向的对象）。
pub(crate) struct ModelBox {
    /// 不可变数据核心（milp::model::Problem）。
    pub problem: Problem,
    /// FFI 参数镜像。
    pub params: Params,
    /// 最近一次成功求解的存档；未求解为 None。
    pub result: Option<StoredSolve>,
    /// panic 中毒标记：置位后除 destroy 外的所有入口返回 7。
    pub poisoned: bool,
}

impl ModelBox {
    /// 创建句柄本体（`milp_ffi_model_create` 成功路径）。
    pub(crate) fn new(name: &str, sense: Sense) -> Self {
        ModelBox {
            problem: Problem::new(name, sense),
            params: Params::default(),
            result: None,
            poisoned: false,
        }
    }
}
