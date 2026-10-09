# milp

纯 Rust 标准库实现的混合整数线性规划（MILP）精确求解器：零第三方依赖、
无 `unsafe`、默认单线程确定性（所有并列一律按最小下标破键；`--threads ≥2`
并行下 objective/status 仍确定，统计量可变——见「并行架构与确定性」
一节），提供人类可读的 `.milp` 算例格式、命令行入口与库 API。

## 定位与能力边界

- **能力**：纯 LP（全连续）、0-1 规划、一般整数与混合类型 MILP 的**精确**
  求解——两阶段有界单纯形 + 深度优先分支定界 + 根节点 Gomory 混合整数割
  （GMI）+ 取整/潜水启发式 + 保守预处理，可选多线程并行与 MIP gap 提前
  终止。`threads = 1` 时求解完全确定，同一算例同一结果（`threads ≥ 2` 见
  「并行架构与确定性」一节）。
- **规模目标**：以 ≤500 结构变量、≤2000 行的正确性为设计目标。分支定界
  使用显式 dense 表格，内存 ≈ m×(n+m+a)×8B（2000 行 × 500 变量约 72MB
  量级）；运行时间随分支节点数增长，由节点/时间上限兜底。
- **诚实性铁律**：只有满足最优性终止条件才报 Optimal；触发迭代/节点/时间
  上限一律报 `LimitReached`（附已有 incumbent），绝不把截断结果谎报为最优。
- **边界**：变量下界必须有限（不支持 −∞ 下界；上界默认 +∞ 合法）；`.milp`
  文法无目标常数 token（`obj_offset` 恒为 0）；仅线性目标与线性约束。

## 构建与 CLI

纯标准库，`Cargo.toml` 无 `[dependencies]`、无 `build.rs`：

```text
cargo build --manifest-path milp/Cargo.toml --release
cargo test  --manifest-path milp/Cargo.toml
```

命令行用法（`--time-limit-secs N` 的秒数经 `N.saturating_mul(1000)` 饱和
映射为 `BbOptions.time_limit_ms`，极值不溢出、不 panic）：

```text
milp <算例路径> [--time-limit-secs N] [--threads N] [--log-secs F] [--mipgap P] [--heuristic-mode auto|on|off] [--branch-rule most_fractional|pseudocost] [--verbose]
```

无参时向 stderr 打印用法；各参数非法（见下表）一律 stderr 中文 `error: ...`
+ 退出码 1。

### 参数

| 参数 | 语义 | 默认值 | 校验（违例 ⇒ 退出码 1） |
|------|------|--------|--------------------------|
| `--time-limit-secs N` | 时间上限秒数（饱和映射毫秒） | 不限 | 非负整数；重复出现 |
| `--threads N` | 线程数；`1` = 纯顺序路径（输出与历史版本逐字节一致），≥2 走并行引擎 | `1` | 整数 `1..=1024`（0/越界/非整数报错）；重复出现 |
| `--log-secs F` | 进度行周期秒数（首行必发，见「进度日志」一节） | 不输出任何进度行 | 有限且 `F > 0`（0/负/NaN/inf 报错）；重复出现 |
| `--mipgap P` | MIP gap 提前终止阈值（百分比，见「MIP gap 提前终止」一节） | 关闭 | 有限且 `P ≥ 0`（`0` ≡ 关闭；负/NaN/inf 报错）；重复出现 |
| `--heuristic-mode auto\|on\|off` | 启发式模式（见「自适应启发式」一节） | `auto` | 未知 token（大小写敏感）报错；重复出现 |
| `--branch-rule most_fractional\|pseudocost` | 分支规则（默认 `most_fractional`；`pseudocost` 含可靠性冷启动，个别问题类可能更优/更差，可切换对比） | `most_fractional` | 未知 token / 重复出现报错 |
| `--verbose` | Optimal 时逐行输出变量取值 | 关 | — |

### 输出格式

- **Optimal**：stdout 依次输出 `status = optimal`、`nodes = N`（处理的
  分支定界节点数；数值取决于算例与实现）；`--verbose` 时再逐行输出
  `name = value` 的变量取值；**最后一行精确为** `objective = <数值>`。
  数值常规为 `{:.6}` 六位小数（如 `objective = 17.000000`）；当
  `0 < |v| < 1e-6` 或 `|v| ≥ 1e15` 时改用科学记数 `{:e}`（如 `1e-7`）。
- **Infeasible / Unbounded / LimitReached**：仅输出
  `status = infeasible|unbounded|limit_reached`，**没有 objective 行**；
  `LimitReached` 且已有 incumbent 时附一行 `best = <数值>`。
- **gap 终止**（`--mipgap P` 且 P>0 触发）：依次输出 `status = gap_reached`、
  `nodes = N`、`gap = <百分比>`（触发时刻证据 + 统一 clamp 双保证为**非负
  有限**，永不 inf/负）、`best = <数值>`（**无 objective 行**——未证明最优，
  绝不冒充）。
- **错误**：向 stderr 输出 `error: <MilpError>`（中文信息，解析错误附
  1 起始行号），退出码 1。

```text
$ cargo run --manifest-path milp/Cargo.toml --release -- milp/instances/knapsack-10.milp
status = optimal
nodes = 1                # 实测（threads=1 默认路径；节点数由求解器统计）
objective = 17.000000
```

### 退出码

| 退出码 | 含义 |
|--------|------|
| 0 | 最优（stdout 最后一行 `objective = <数值>`） |
| 1 | 错误/用法（加载、解析、模型非法、参数非法、数值异常等） |
| 2 | 不可行（Infeasible） |
| 3 | 无界（Unbounded） |
| 4 | 未证明最优的提前终止——`limit_reached`（时间/节点上限触发）**或** `gap_reached`（`--mipgap` 达标触发），均附已有 incumbent |

## 进度日志（--log-secs）

- **默认关闭**：不带 `--log-secs` 时 stdout 不出现任何 `progress` 行，
  保持既有输出契约逐字节稳定（脚本可解析不受扰动）；需要观测搜索过程时
  显式开启。
- **首行契约（不依赖计时）**：给出 `--log-secs` 后，首行进度在**根阶段
  结束的唯一汇合点**无条件发射——覆盖一切捷径（纯 LP 单节点直返、根
  Infeasible/Unbounded 终局、根割循环后终局均先发射再返回），先于任何
  周期行与 `status` 行。
- **周期行**：此后每 `F` 秒一行，**仅在搜索未结束时发射**；顺序路径以
  待处理栈非空为前提（空栈轮次不发射，杜绝 `incumbent=X bound=none
  gap=inf` 矛盾行）；并行协调线程接管节奏即发首周期行（此刻树搜索进行中、
  池/在处理槽非空，不会空聚合出 `bound=none`），此后每 `F` 秒一行，
  停止/完成标志置位后只退出不打印——亚毫秒搜索配极小 `F`（如 0.001）时
  不依赖 OS 定时精度也必有周期行。
- **行格式**：`progress t={耗时秒} nodes={已处理节点数} incumbent={当前最好解或none} bound={全局对偶界或none} gap={百分比或inf}`，例如：

```text
progress t=0.12 nodes=42 incumbent=511.000000 bound=526.870859 gap=3.11%
```

- **gap 恒非负**：展示与触发判定统一走 `search_ctrl::gap_percent`，raw
  值一律 clamp 到 ≥0——`db ≤ inc`（incumbent 跳跃越过陈旧估计，或浮点
  噪声如 est=39.999999999999996 < inc=40.0）时输出 `0.00%` 而非 `-0.00%`；
  缺任一侧或相对口径无意义时打印 `inf`。
- Optimal 终局仍以 `objective = <数值>` 最后一行收尾；进度行写 stdout
  失败（如管道断裂）静默忽略。

## MIP gap 提前终止（--mipgap）

- **`P = 0` 等价关闭**（与不给出该参数完全一致）；`P > 0` 时，剩余对偶界
  与 incumbent 的相对 gap ≤ P% 即提前终止。
- **终止行为**：库层终态为 `LimitReached` + `gap_terminated = true`（诚实：
  未证明最优），CLI 依次输出 `status = gap_reached`、`nodes = N`、
  `gap = <百分比>`、`best = <数值>`（无 objective 行），退出码 4；`best`
  为触发时刻的 incumbent（真实可行解，经安全网终检）。
- **gap 口径**：`(方向差)/|incumbent| × 100`（Max 用 `bound − inc`，Min
  镜像），其中 `bound` 为**触发时刻证据快照**——顺序路径在主循环顶部判定
  即断即停（剩余子树 est 保守聚合，粘滞、不依赖后续栈态）；并行路径由
  stop CAS 的胜者独占写入证据（绝不 join 后重聚合）。CLI 报告的 gap 以
  该证据重算，与触发判定同源。
- **clamp ≥ 0**：`db ≤ inc` 时 gap 记 0% 而非负值——此时全部剩余子树
  `est ≤ inc` ⇒ 最优 ≤ `max(inc, db)` = incumbent ⇒ incumbent 在 est 有效
  性内已达最优，0% 触发是诚实的；同时杜绝 `-0.00%` 之类的矛盾输出。
- **诚实性**：gap 终止绝不冒充 optimal（无 objective 行、退出码 4）；
  剩余子树的 est 均为父谱系 LP 目标（对其子树有效），聚合保守 ⇒ 触发
  只会推迟、不会虚假。

## 自适应启发式（--heuristic-mode）

- **三档语义**：`off` ⇒ 关闭启发式（`BbOptions.heuristics = false`）；
  `on` ⇒ 仅根节点启发式（取整复核 + dive，即历史版本行为）；
  `auto`（默认）⇒ 根节点启发式 + **自适应**树节点 dive。
- **根预割启发式轮**（2026-10）：割循环可能把根 LP 顶点从「单断点近整数」
  推向「多分数中面」（界更紧但取整/潜水质量塌缩——hard-knapsack-12 割后
  顶点上 dive 仅得 453）。root_cuts 开启时在**割前顶点**先跑一轮
  round_and_check + dive、与割后一轮取最优：incumbent 只增不减，割生成
  只读表格不受影响；关割路径顶点唯一，行为逐字节不变。
- **控制器策略**（`search_ctrl::HeuristicController`，参数
  256/400/8/0.25/0.8，2026-10 随 dive 深预算化重调）：『无 incumbent ⇄
  连续 400 节点无改进（stall）⇄ 近期 8 次 dive 命中率 ≥ 0.25』三条件任一
  满足才授权 dive，叠加 **256 节点最小间隔节流**与**时间预算保护**（有
  `--time-limit-secs` 且已耗 80% 预算并有 incumbent 时停投，把尾段留给
  最优性证明）。理由：对偶热修复下单次 LP 重解 ≈ 数次迭代，一次深潜
  （预算 4096）≈ 数百次热重解——期望回报（更快拿到好 incumbent 以加速
  剪枝）不够高时不投资。
- **正确性**：dive 只供应 incumbent，不改剪枝/分支语义 ⇒ `auto` 仍**精确
  最优**（目标与 `on`/`off` 一致）。
- **原 7 个验收算例**（节点数 ≤ 7）均低于节流/阈值 ⇒ `auto ≡ on`，CLI
  默认输出逐字节一致（由验收测试钉死）。
- **hard-knapsack-12**：184 节点 / objective 523（2026-10 启发式增强后
  的确定性基线——根 dive 现在能供应 incumbent；增强前根 dive 在断点链
  上必耗尽 8 次预算而失败。根预割启发式轮加入后根 incumbent 453→511 ⇒
  剪枝提前，235→184）。一般而言中途 dive 若改进 incumbent 可使
  节点数**变少**（剪枝更早），`auto` 的 nodes/迭代数不承诺与 `on` 相同
  （目标与终态恒同）。

### 策略级自动降级（auto + --time-limit-secs）

在「树节点上要不要投资 dive」的微观自适应之上，`auto` 还有一层**策略级**
自适应：自动判断「继续证明最优」还是「提前返回当前最好解」。

- **激活条件（缺一不可）**：`--heuristic-mode auto`（默认）**且**给了
  `--time-limit-secs N`。不给时间预算 = 没有表达速度优先的需求 ⇒ 契约
  红线：auto 仍**证明最优**，绝不降级。
- **触发判定**（`search_ctrl::AutoDegrade`，默认参数 0.5/0.25/0.01，
  经 `SearchControl::degrade` 可注入）：预算消耗 ≥ 50% 之后，若
  ① 已有 incumbent 且连续 ≥ 25%×预算 无改进，② 且同期对偶界的相对
  改进 < 1%（**双停滞**：解侧与证明侧都不再前进）⇒ 判定继续证明的边际
  收益过低，提前终止。界每有实质进展会刷新停滞基准——证明仍在收敛时
  绝不放弃。
- **输出契约**：触发时输出 `status = auto_degraded`、`nodes = N`、
  `gap = <%.2f>`（触发时刻证据）、`best = <数值>`，**无 objective 行**
  （未证明最优，绝不冒充），退出码 4（与 limit/gap 终止同级）；库层
  `ControlledSolveReport` 以 `auto_degraded = true` 区分于 `--mipgap`
  触发（`gap_terminated = true` 两者同置）。
- **与 `--mipgap` 的关系**：`--mipgap` 是手动的「gap 到多少就停」；
  自动降级是求解器自己判断「这预算内证不完了，别耗着」。两者可并用，
  谁先满足谁终止。

## 并行架构与确定性（--threads）

- **threads = 1（默认）**：纯顺序路径，输出与历史版本**逐字节一致**
  （所有并列一律按最小下标破键）。
- **threads ≥ 2**：调用线程先完成与顺序路径同构的根阶段（根 LP → 纯 LP
  捷径 → 根割循环（记录割清单）→ 根启发式），随后把根子节点（各自携带
  完整 `Bounds` 快照与父节点 LP 目标 est）压入共享**节点池**，经
  `std::thread::scope` 派生 N 个 worker 抢占式处理；worker 各持独立
  `Simplex` 实例，按序 `add_row` 重放根阶段记录的割清单，任意界跳转由
  snap/Phase-1 重建处理。
- **共享状态与同步原语**：节点池 `Mutex`（互斥弹出 ⇒ 每节点恰处理一次）；
  incumbent `RwLock`（单调改进，剪枝判定与顺序路径逐字相同 ⇒ 永不丢更优
  解）；**双 Condvar**——`cv_workers` 供 worker 等节点、`cv_progress` 供
  协调线程按 `--log-secs` 周期发射进度行，互不唤醒；**粘滞 stop** 原子量
  （CAS 一次性置位：节点/时间上限与 gap 达标，镜像顺序路径 limit_hit
  语义 ⇒ 截断绝不冒充完备）；节点 LP 迭代上限与顺序 (c) 分支同口径——
  剪枝该节点 + 置独立粘滞标记 `limit_hit`（不终止树搜索），探索完备后
  终态仍诚实 LimitReached；**gap 证据**为触发时刻快照（CAS 胜者独占
  写入，叶级短锁）。节点上限采用**预约制**（弹出后按预约序号判定，
  与顺序路径 `nodes >= node_limit` 逐节点对齐）。
- **借用交接（无 'static）**：`Shared` 为调用线程栈上局部结构体，worker
  闭包在 `thread::scope` 内直接借用 `&Problem`、`&割清单`、`&Shared`
  （问题与割的交接即经此借用传递，不装入 `Shared`、无需 `Arc<Shared>`）；
  协调线程（进度节奏）同样在 scope 内于调用线程执行；
  `ScopedJoinHandle::join` 的 `Err`（worker panic）映射
  `MilpError::Internal`，不外泄 panic。
- **锁纪律**：任一线程任一时刻**至多持一把锁**，绝不嵌套；gap 检查两段式
  （池锁内聚合 → 释放 → incumbent 读 → 释放 → 无锁计算 + CAS）；锁中毒
  `into_inner` 恢复并置 stop。
- **Send/Sync 契约**：仅 std、零外部依赖、无 unsafe；纯数据类型
  （`Problem`/`Bounds`/`Tolerances`/`BbOptions`/`Cut`/`ProgressInfo`/
  `SearchControl` 等）自动 `Send + Sync`；`Simplex` 仅 `Send`（每 worker
  独占，绝不共享）；进度回调为 `Arc<dyn Fn(&ProgressInfo) + Send + Sync>`。
- **内存提示**：每 worker 一份 dense 表格 ⇒ 内存 ≈ `threads × m×(n+m+a)×8B`
  （2000 行 × 500 变量约 72MB/线程量级），线程数请按内存预算选取。
- **确定性**：`threads = 1` 逐字节一致；`threads > 1` 时 **status 与最优
  objective 确定**（池互斥弹出每节点恰处理一次、子节点必入池 ⇒ 树覆盖
  无遗漏；incumbent 单调 + 剪枝与顺序路径逐字相同 ⇒ 永不丢最优；粘滞
  stop 保证截断绝不冒充 Optimal），而 **nodes/迭代数/耗时/取值（多最优
  时取哪个）跨运行可变**——验收测试只断言 objective/status 与粘滞证据
  字段。

## .milp 算例格式

文法（行式解析，`#` 起注释至行尾，空白分隔 token，大小写敏感）——
WP1 契约逐字引用：

```text
[name <算例名>]          # 可选；算例名=单个非空白 token，允许 '-''_''.'，
sense max|min            #   仅显示用途，不得含 '#'
vars
  <变量名> continuous|integer|binary   # 每行一个；变量名 [A-Za-z_][A-Za-z0-9_.]*，
                                        # 不重复、不得为保留字
obj                      # <系数> <变量名> 对，可跨多行直到下一节关键字；
                          #   数值为 Rust f64 语法（允许 -6、2.5、1e-2）
cons                     # 每行一条：[约束名:] <系数> <变量名> ... <=|>=|= <右端>
                          #   约束名=去尾 ':' 的非空白序列（允许 '-'）；
                          #   系数 0.0 的项解析层直接丢弃
[bounds                  # 可选；每行：<变量名> <=|>=|= <数值>，
  ...]                   #   叠加在类型默认界上取交集；Binary 只允许 [0,1] 内
end                      # 必须存在
```

- **保留字**（不得用作变量名，撞之报解析错误）：`name`、`sense`、`vars`、
  `obj`、`cons`、`bounds`、`end`、`continuous`、`integer`、`binary`、
  `max`、`min`。注意算例名与约束名允许含 `-`（WP1 文法定义，与变量名正则
  无关）。
- 变量默认界：Continuous/Integer 为 `[0, +∞)`、Binary 为 `[0, 1]`；文法不
  提供无穷字面量，默认上界 +∞ 以**省略** bounds 行表达。
- 确定性补充约定（均以 `Err(Parse{line, ..})` 报告，1 起始行号）：一切数值
  字面量必须解析为**有限** f64（`inf`/`NaN`/`1e999` 拒绝）；节头行
  （`vars`/`obj`/`cons`/`bounds`/`end`）必须独占一行；obj 节内、单条约束
  行内同一变量重复出现即报错（不静默归并）；bounds 与现有界交集为空或
  Binary 界值超出 [0,1] 即报错（整数变量非空但取整后为空的区间合法通过，
  交由预处理/单纯形按空界判不可行）；`end` 之后不得再出现任何非注释 token。

示例（`instances/milp-mixed-5.milp`，binary + 一般整数 + 连续混合）：

```text
name milp-mixed-5
sense max
vars
  y1 binary
  y2 binary
  z integer
  c continuous
obj
  10 y1 8 y2 4 z 1 c
cons
  cap: 6 y1 5 y2 3 z 1 c <= 20
  count: 1 y1 1 y2 1 z <= 4
  linkc: 1 c -6 y1 <= 0
  linkz: 1 z -2 y2 >= 0
bounds
  z <= 10
  c <= 10
end
```

## 算法说明

求解管线（`milp::solver`）：`validate` →（可选）presolve → 分支定界 →
postsolve 还原 → 安全网终检 → 报告。库 API 统一 `Result<T, MilpError>`
（五变体：Io/Parse/InvalidModel/Numerical/Internal）；用户输入路径禁止
`unwrap`/`expect`/`panic`。

### 两阶段有界变量单纯形（WP2）

- 显式 dense 表格 + 有界变量机制（非基变量位于 lb 或 ub）；分支定界只改
  `Bounds` 不改约束矩阵，热启动代价最低。内部统一 min（Maximize 目标系数
  取负，报告时还原并加 `obj_offset`）。
- **按行归一、行永不翻转**：构建期每行除以 `max|a_ij|` 后施用容差，杜绝
  除零/NaN；行不翻转保证松弛列符号只有两个合法来源（Le=+1、Ge=−1），对偶
  与解的还原无符号歧义。
- **带符号人工变量初始基**：以 x⁰ = 各变量取 lb 计算残差 r_i = b̂_i − â_i·x⁰，
  Le 行残差非负松弛列入基、Ge 行残差非正剩余列入基，否则加人工列
  （σ=±1、值为正）入基；Eq 行恒加人工列。人工变量绝不限于等式行——正右端
  ≥ 行（如 `x ≥ 3` 而 x 停在 lb=0）的剩余量恒为负，取反行救不了它；旧的
  “b_i<0 整行取反”规则已废除。
- Phase-1 目标 min Σ人工变量，最优值超容差 ⇒ Infeasible；Phase-2 用原目标。
  近零行（`max|a_ij| ≤ tol.zero`）不入表，按原尺度 rhs 判定不可行或永久剔除。
- **双侧定价合格性**：贴下界且 rc < −tol.optimality（增向入基）或贴上界且
  rc > +tol.optimality（减向入基）——贴上界列的减向入基是二进制/有限上界
  变量的常态操作，缺它会在非最优点误报 Optimal。
- **完整 Bland 反循环规则**：Dantzig 阶段（|rc| 得分最大、并列最小列下标）
  连续退化迭代 > 50 次或总迭代超过 3×(m+n) 后切入 Bland 模式（本次调用内
  永久生效）：入基取合格列最小下标；出基侧全部并列候选（被阻挡基变量、
  入基变量自身界翻转）一律取最小下标。比值测试同时考虑基变量上下界与
  入基变量自身界翻转；|主元| ≥ tol.pivot 才允许；入基列无阻挡 ⇒ Unbounded。
  检验数每次由当前表格重算；每 200 次主元用原始数据重算 B⁻¹b 抑制漂移；
  迭代上限 `max(10_000, 20×(m+n))` ⇒ IterationLimit（诚实上报）。

### 深度优先分支定界（WP4）

- 单一 `Simplex` 实例 + `Bounds` 界栈（trail）+ 深度优先，节点零拷贝；下行
  分支先探索；回溯严格按 trail 逆序恢复旧界，**非基双向 snap**（收紧越界
  snap / 恢复后严格内点 snap 到较近界、等距按 rc 符号选侧）处理后热启动。
- 节点处理顺序固定：LP Infeasible ⇒ 剪枝；根节点 LP Unbounded ⇒ 整树判
  Unbounded（见理论近似声明）；IterationLimit ⇒ 剪枝并记 limit_hit；LP bound
  不优于 incumbent ⇒ 剪枝；整数性满足 ⇒ 取整复核（伪整数解**不得**剪枝，
  必须继续分支）；空界 ⇒ 剪枝。
- 分支变量规则 MostFractional（max min(f, 1−f)，并列取最小下标）；可选
  [`BranchRule::Pseudocost`]（`BbOptions::branch_rule`，编程 API；无观测时
  退化为 f⁻·f⁺ 乘积规则——目前仅顺序路径生效）；剪枝条件
  `node_bound ≤ incumbent + 1e-9×(1+|incumbent|)`。
- **界变化后的对偶修复快路径**（2026-10，大规模主优化）：分支收紧界后若
  基变量越界，不再整基重建 + 完整 Phase-1（冷重解），而是有界变量对偶
  单纯形几次旋转恢复原始可行（基与成本未动 ⇒ 对偶可行性保持；出基变量
  落到其被违反的界）。前置条件不满足（Phase-1 欠账 / 基中含人工列 /
  对偶不可行 / 超预算 / 数值异常）一律回退冷路径，快路径只加速、绝不
  改变语义。实测每节点 LP 迭代 123 → ~2.5（约 46×），节点吞吐约 30×。
  一致性由 `tests/warm_cold_parity.rs` 差分钉死（界收紧序列上热解与冷解
  逐点对比状态与目标）。
- **深潜启发式**（2026-10）：`dive_max_fixes` 默认 8 → 4096（对偶热修复
  使单次 resolve 成本 ≈ 数次迭代，深潜变得便宜）；预算耗尽或修复不可行
  时退化为取整复核而非直接放弃；修复不可行先做**一步反向回退**（ceil
  不行改 floor 再试，每变量至多一次，防振荡）。根 dive 由此在背包类
  断点链上从必然失败变为稳定供应 incumbent。
- 节点上限默认 2_000_000（吞吐 ↑ ~30× 后旧 10 万上限秒级打满）。
- 终态诚实判定：栈空自然退出且有 incumbent ⇒ Optimal、无 ⇒ Infeasible；
  节点/时间上限触发 ⇒ LimitReached（附 incumbent）。
- **基准套件**：`bench/`（`gen_instances.py` 生成，种子确定）：强相关背包
  200/500（附 DP 独立核对最优）、多背包 3×100、指派 16×16、混合 MILP、
  宽 LP。2026-10 优化后（对偶修复 + 深潜 + 控制器重调 256/400）：可完成
  算例单线程 2.1–3.3×、4 线程 2.2–3.4× 加速，目标值与基线全部一致。

### 根节点 GMI 割（WP3）

只在根节点生成（对原问题全体整数可行解全局有效；分支子树内调用无效，
禁止）。取基变量为整数类型结构变量的行推导，割为 `Σ π_j·d_j ≥ 1`（π_j ≥ 0）。
**整数列前提**：整数类型结构列仅当所贴界为整数值（|bound−round| ≤
tol.zero）时才可用整数公式，否则必须按连续列公式处理，否则割可切掉整数
最优而报 Optimal。数值纪律：|c_j| < tol.zero 的项只允许按放宽割的方向
丢弃，严禁加严；violation ≥ min_violation 才收；动态范围超 1e9 整条丢弃。
割循环后根 LP 不可行是整数不可行问题的合法终局（割对全部整数可行解有效、
该集合为空时割可切空松弛区域）——终态如实映射 Infeasible（CLI
`status = infeasible`、退出码 2），不是数值错误。

### 预处理（WP3）

保守规则仅四类且可证明安全，每轮严格按 **R1→R2→R3→R4** 顺序施加，禁止
重排，最多 5 轮无化简即停：R1 整数变量界取整；R2 固定列（单点界）代入删列；
R3 空行判定/删行；R4 单系数行换算为变量界后删行。

### 启发式（WP3）

`round_and_check`（简单取整 + 可行性复核）与 `dive`（fix-and-dive，迭代
≤ 8 次固定最接近整数的非整变量后重解 LP）只供应 incumbent 以加速剪枝，
不影响正确性；复核容差与 integrality 联动（见下）。

### 容差表（集中定义于 `model::Tolerances`，库内禁止散落 magic number）

| 容差 | 默认值 | 用途 |
|------|--------|------|
| `feasibility` | 1e-9 | 原始/对偶可行性（按行归一后施用） |
| `optimality` | 1e-9 | 检验数最优性 |
| `pivot` | 1e-9 | 主元最小绝对值 |
| `integrality` | 1e-6 | 整数判定：\|x−round(x)\| ≤ 1e-6 |
| `zero` | 1e-12 | 系数清零阈值 |

- 剪枝容差：`node_bound ≤ incumbent + 1e-9×(1+|incumbent|)`。
- **取整复核容差联动公式**（容差阶梯一致性）：约束
  `violation ≤ 10×tol.feasibility + tol.integrality×‖row‖₁`（‖row‖₁ = 该行
  系数绝对值和）——保证被判“整数性满足”的节点解（每变量偏差 ≤ 1e-6）在行
  范数放大后仍能通过复核；复核失败的伪整数节点不得剪枝，必须继续分支。
- **安全网终检**（solver）：对原始约束
  `violation > 1e-6×(1+|rhs|+‖row‖₁)` 即返回 `Err(Numerical)` 而非输出错解；
  报告的整数变量值已取整且目标按该值重算（值与目标自洽）。

## 理论近似声明

- **整数变量存在而根 LP 松弛无界时直接判 Unbounded**：理论上原整数问题
  可能不可行（LP 无界方向可能不存在整数可行点），本契约接受该近似。
- **割平面仅在根节点生成**：以保证割对原问题全体整数可行解全局有效
  （分支子树内的割只对该子树有效，混用会破坏正确性）。

## 验收算例

`milp/instances/` 下 8 个算例（6 个数值型 + 2 个状态型）；期望值与状态
判定均已用独立穷举脚本实跑验证。验收测试见 `milp/tests/acceptance.rs`
（数值解析比较，容差 1e-6，不依赖具体小数位数；原 7 算例的 CLI 默认
输出另以逐字节基线钉死）。

| 算例 | 类型 | 期望结果 | 已知最优解 | 验证方式 |
|------|------|----------|------------|----------|
| `lp-prod-3` | 纯 LP（3 连续变量） | Optimal，objective = 13 | (x1,x2,x3) = (2,0,1) | 3×3 顶点枚举 |
| `knapsack-10` | 0-1 背包（5 变量） | Optimal，objective = 17 | {x2,x3,x5} | 2^5 穷举；LP 松弛 17.667 触发分支 |
| `assign-4x4` | 4×4 指派（16 binary） | Optimal，objective = 10 | W1→J3, W2→J1, W3→J4, W4→J2 | 4!=24 穷举；LP 松弛即整最优且高度退化（Bland 规则测试点） |
| `milp-wolsey-2` | 一般整数 MILP | Optimal，objective = 40 | (x1,x2) = (0,5) | 穷举；LP 松弛 41.25 于 (2.25,3.75) 触发分支 + GMI 机会；整数可行点 25 个（x2=0..5 时 x1 上限 6,5,4,3,1,0），GMI 公式对 25 点全部验证有效 |
| `milp-mixed-5` | binary+整数+连续混合 | Optimal，objective = 29 | (y1,y2,z,c) = (1,1,2,3) | 4 个二元组合 × z 枚举并解析求 c；根 LP 松弛即整数最优 |
| `hard-knapsack-12` | 0-1 背包（12 变量，权重/价值/容量 7..47 / 11..144 / 172） | Optimal，objective = 523 | 184 节点（threads=1；并行下节点数可变） | 2^12 = 4096 穷举；实测 184 节点 ≫ 原算例（≤7），作为 `--threads 2/4` 真实压满节点池、`--log-secs` 周期行与 `--heuristic-mode` 的并发/受控参数验收载体 |
| `infeasible-tiny` | 状态型 | Infeasible，退出码 2，无 objective 行 | — | x≥3 ∧ x≤1 空集；正右端 ≥ 行（x≥3 在 lb=0 处违反，走人工变量路径）端到端覆盖 |
| `unbounded-tiny` | 状态型 | Unbounded，退出码 3，无 objective 行 | — | max x，x≥3 且上界 +∞；默认上界 +∞ 的合法性由本例顺带验收 |

## 仓库布局

```text
milp/
├── src/
│   ├── main.rs          # CLI 入口（WP5 + 受控搜索参数：threads/log-secs/mipgap/heuristic-mode）
│   ├── error.rs         # MilpError 与 Result 别名（WP1）
│   ├── model.rs         # 建模 + 集中容差（WP1）
│   ├── instance.rs      # .milp 解析/序列化（WP1）
│   ├── simplex.rs       # 两阶段有界单纯形（WP2）
│   ├── presolve.rs      # 保守预处理与还原（WP3）
│   ├── cuts.rs          # 根节点 GMI 割（WP3）
│   ├── heuristics.rs    # 取整/潜水启发式（WP3）
│   ├── branch_bound.rs  # 深度优先分支定界 + 受控入口（WP4）
│   ├── solver.rs        # 统一编排、受控求解入口与安全网终检（WP4）
│   ├── search_ctrl.rs   # 搜索控制层：进度快照/gap 口径/自适应控制器（纯逻辑零锁）
│   └── parallel_bb.rs   # 并行分支定界引擎（节点池/incumbent/双 Condvar/粘滞 stop）
├── instances/           # 8 个验收算例（7 原算例 + hard-knapsack-12）
└── tests/acceptance.rs  # 验收测试（WP5 + 受控搜索参数扩展）
```

单元测试 `#[cfg(test)]` 就近放在各模块内；验收统一在 `tests/acceptance.rs`。
