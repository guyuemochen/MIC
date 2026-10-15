# milp

纯 Rust 标准库实现的混合整数线性规划（MILP）精确求解器：零第三方依赖、
无 `unsafe`、默认单线程确定性（所有并列一律按最小下标破键；`--threads ≥2`
并行下 objective/status 仍确定，统计量可变——见「并行架构与确定性」
一节），提供人类可读的 `.milp` 算例格式、命令行入口与库 API。

## 定位与能力边界

- **能力**：纯 LP（全连续）、0-1 规划、一般整数与混合类型 MILP 的**精确**
  求解——两阶段有界单纯形**双引擎**（稠密表格 / 稀疏修正单纯形，S4 起
  `--lp-engine auto` 类自适应默认；大规模启用稀疏簿记、退化停滞的 Bland
  释放）+
  深度优先分支定界（类自适应分支默认）+ 根节点 GMI/聚合 MIR 割（循环预算
  收紧 + 失败轮回滚）+ 取整/潜水/子集和 MITM 启发式 + 数值缩放与保守预
  处理，可选多线程并行（本地亲和或全局 DFS 工作栈）与 MIP gap 提前终止。
  `threads = 1` 时求解完全确定，同一算例同一结果（`threads ≥ 2` 见
  「并行架构与确定性」一节）。
- **规模目标**：以 ≤500 结构变量、≤2000 行的正确性为设计目标；2026-10
  大规模优化轮（R1–R8，W7 起阈值为 2^16）后，表格单元数超阈值的大算例（如 bench 的
  lp-sparse-4000x800、setcover-400x800）启用稀疏簿记路径以压低每迭代成本
（小算例逐字节走旧稠密路径，行为零扰动）；S2–S4 起再加一层引擎级类自适应
（`--lp-engine auto` 默认：表格 m×(n+m) ≥ 10^6 格走稀疏修正单纯形 revised
引擎，lp-sparse-4000x800 大幅提速、目标恒同——迭代序列 S6 起随 Devex
定价有意分叉，见「双引擎架构」一节）。分支定界使用显式 dense 表格，
  内存 ≈ m×(n+m+a)×8B（2000 行 × 500 变量约 72MB 量级）；运行时间随分支
  节点数增长，由节点/时间上限兜底。逐轮实测数据见
  `docs/large-scale-report.md`。
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
milp <算例路径> [--time-limit-secs N] [--threads N] [--log-secs F] [--mipgap P] [--heuristic-mode auto|on|off] [--branch-rule most_fractional|pseudocost] [--scaling=on|off] [--mitm=on|off] [--work-stack auto|local|global] [--determinism wall|budget] [--lp-engine auto|tableau|revised] [--verbose]
```

无参时向 stderr 打印用法；各参数非法（见下表）一律 stderr 中文 `error: ...`
+ 退出码 1。

### 参数

| 参数 | 语义 | 默认值 | 校验（违例 ⇒ 退出码 1） |
|------|------|--------|--------------------------|
| `--time-limit-secs N` | 时间上限秒数（饱和映射毫秒）。两级检查（W2）：节点边界 + **LP 迭代循环内部**（每 64 次主元条件读钟，超限诚实 `limit_reached`）——单节点纯 LP 不再越窗跑满；**时限关闭（默认）时求解路径零墙钟读点**，t1 确定性逐字节不变 | 不限 | 非负整数；重复出现 |
| `--threads N` | 线程数；`1` = 纯顺序路径（输出与历史版本逐字节一致），≥2 走并行引擎 | `1` | 整数 `1..=1024`（0/越界/非整数报错）；重复出现 |
| `--log-secs F` | 进度行周期秒数（首行必发，见「进度日志」一节） | 不输出任何进度行 | 有限且 `F > 0`（0/负/NaN/inf 报错）；重复出现 |
| `--mipgap P` | MIP gap 提前终止阈值（百分比，见「MIP gap 提前终止」一节） | 关闭 | 有限且 `P ≥ 0`（`0` ≡ 关闭；负/NaN/inf 报错）；重复出现 |
| `--heuristic-mode auto\|on\|off` | 启发式模式（见「自适应启发式」一节） | `auto` | 未知 token（大小写敏感）报错；重复出现 |
| `--branch-rule most_fractional\|pseudocost` | 分支规则；**缺省 = 类自适应**（R3：全 Binary + 全正系数 Le 行 + 行数 2..=64 的多维背包类自动用 `pseudocost` 含可靠性冷启动，否则 `most_fractional`；P3 收紧行数下限 2——单行背包类探测零收益纯开销，回落 `most_fractional`）；显式指定覆盖优先 | 类自适应（auto） | 未知 token / 重复出现报错 |
| `--scaling=on\|off` | 问题级数值缩放（R4：行列迭代几何均值，2 的幂系数；目标/解还原原尺度，输出契约不变） | `on` | 非 `on\|off` token / 重复出现报错 |
| `--mitm=on\|off` | 子集和 MITM 可行性启发式（R5：全 Eq + 全 Binary 等式类的最后手段，见「启发式」一节） | `on` | 非 `on\|off` token / 重复出现报错 |
| `--work-stack auto\|local\|global` | 并行工作栈模式（R6 接线，W1 默认按类选栈；仅 `--threads ≥ 2` 生效，threads=1 不受影响）：`auto` **多线程默认**——类检测命中多维背包类（全 Binary + 正系数 Le 行 + 行数 2..=64）⇒ `global`，单行背包/等式行/Ge 行/连续列类 ⇒ `local`；`local` 本地亲和 + 偷取；`global` 单一全局 DFS 栈，并行探索序对齐顺序 DFS | `auto` | 非 `auto\|local\|global` token / 重复出现报错 |
| `--determinism wall\|budget` | 确定性模式（W3）：`wall`（默认）= 现状（时限/自适应/降级按墙钟判定，长跑终态随机器状态漂移）；`budget` = 纯预算制——`--time-limit-secs` 换算为节点预算（6,000 节点/秒），全部墙钟判定改为等价节点数预算，决策路径零墙钟读点 ⇒ **同参数双跑全部行为字段逐字节相同**（仅 `--threads 1`，与 ≥2 组合报错；慢机上墙钟可超名义值——确定性优先） | `wall` | 非 `wall\|budget` token / 重复出现报错 |
| `--lp-engine auto\|tableau\|revised` | LP 引擎选择（S3 接线、S4 默认 auto，见「双引擎架构」一节）：`tableau` = 稠密表格两阶段单纯形（`simplex.rs`，历史默认路径）；`revised` = 稀疏修正单纯形（`simplex_sparse.rs`：CSC + 稀疏 LU + ETA + 对偶修复热启动，lp-sparse-4000x800 大幅提速、目标恒同——定价规模分档：m×(n+m) ≥ 10^6 格阈上 Devex（S6）、阈下 Dantzig 逐字节兼容旧版，迭代序列不承诺与稠密引擎逐次一致）；`auto`（默认）= 类自适应——表格规模 m×(n+m) ≥ 10^6 格 ⇒ revised，阈下留 tableau（小算例走旧路径防 codegen 扰动）；显式指定覆盖优先 | `auto` | 未知 token / 缺参 / 重复出现报错；`revised` × `--threads ≥2` 报错（S3 范围仅 t1；auto 解析为大算例 revised × ≥2 由库层按有效引擎拒绝，不静默分叉） |
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
- **hard-knapsack-12**：159 节点 / objective 523（P3 起的确定性基线——
  hard-knapsack-12 为**单行**背包，类检测行数下限收到 2 后不再命中，默认
  回落 MostFractional；目标 523 恒为已知最优；演进史 235→184（根
  预割启发式轮使根 incumbent 453→511、剪枝提前）→159（默认规则复核回退
  MostFractional）→196（R3 类自适应命中背包类翻 Pseudocost + 可靠性
  冷启动）→159（P3 行数下限 1→2，单行类回落；`--branch-rule
  pseudocost` 显式复现 196 旧径）。一般而言中途 dive 若改进 incumbent
  可使节点数**变少**（剪枝更早），`auto` 的 nodes/迭代数不承诺与 `on`
  相同（目标与终态恒同）。

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

### 确定性模式（--determinism budget，W3）

墙钟判定的代价是**终态漂移**：`wall` 模式（默认）下时限截断、自适应
预算保护与自动降级全部按 `start.elapsed()` 判定，45s 级长跑的节点数/
终态随机器功率态漂移——mknap-5x300 `--time-limit-secs 30` 同参数跨
运行在 `optimal 4762@87,800` 与 `auto_degraded 4754@39,533` 之间漂移
（同窗实测甚至出现 optimal@21.1s 与 degraded@15.0s 并存，见
`docs/parallel-profiling.md` §12）。

`--determinism budget`（仅 threads=1）把全部墙钟判定改为等价**节点数**
预算：

- `--time-limit-secs N` 换算为节点预算 `⌈N × 6,000⌉`（固定常量、机器
  无关；按 mknap-5x300 健康态 t1 吞吐 ~6.4-6.9k 节点/s 标定取 6,000——
  30s 名义预算 = 180k 节点 > 该例完备需求 87,800 且降级资格线 90k 亦
  高于之 ⇒ budget-30s 恒完备最优，消除漂移）；
- 自适应预算保护（elapsed > 0.8×预算）与自动降级（start_frac/plateau
  分数与停滞窗）的注入单位从秒改为节点数（`search_ctrl` 判定逻辑只做
  同单位无量纲比较，结构不变）；
- 节点边界截断改查节点预算；Simplex 内部 deadline 不注入（W2 条件读钟
  在 budget 模式零执行）⇒ **决策路径零墙钟读点**（进度行时刻与
  `stats.time_ms` 等展示/报告字段仍读钟，不参与任何判定）；
- 效果：**同参数双跑全部行为字段（status/objective/values/nodes/
  lp_iterations/cuts/bound/…）逐字节相同**；无时限时行为 ≡ wall 无时限
  （mknap-3x100 双模式全等实测）；
- 语义边界：确定性优先于墙钟诚实——慢机上墙钟可超名义值（换算率不随
  机器浮动正是可复现性的来源）；与 `--threads ≥ 2` 组合显式报错（并行
  无确定性契约）。

## 并行架构与确定性（--threads）

- **threads = 1（默认）**：纯顺序路径，输出与历史版本**逐字节一致**
  （所有并列一律按最小下标破键）。
- **threads ≥ 2**：调用线程先完成与顺序路径同构的根阶段（根 LP → 纯 LP
  捷径 → 根割循环（记录割清单）→ 根启发式），随后把根子节点（各自携带
  完整 `Bounds` 快照与父节点 LP 目标 est）压入共享**节点池**，经
  `std::thread::scope` 派生 N 个 worker 抢占式处理；worker 各持独立
  `Simplex` 实例，按序 `add_row` 重放根阶段记录的割清单，任意界跳转由
  snap/Phase-1 重建处理。
- **工作栈模式（大规模轮 R6 接线，W1 起默认按类选栈）**：`--work-stack
  local` = 本地亲和 + 偷取；`--work-stack global` = 单一全局 DFS steal
  栈——`complete_node` 把子代压全局栈顶、worker 恒弹栈顶（LIFO），并行
  探索序对齐顺序 DFS（谁弹到谁延续 t1 的下一节点），incumbent 时序/剪枝
  与 t1 同构，抑制并行搜索膨胀。实测 mknap-5x300 t4：local 29.3–33.3s
  （nodes 193k–251k）→ global 16.1s（nodes ≈100k），墙钟 −45%~−52%；
  t1 默认路径 15 算例行为逐位不变。诚实边界：t4 聚合吞吐仍被单锁池竞争与
  并行子树交错吃掉（t4/t1 ≈ 0.6–1.15×，本架构下无稳定并行加速），不承诺
  线性加速。**W1：`--work-stack auto`（多线程默认）**——global 是按问题
  类别分化的杠杆（mknap 族 t4 节点 −26%~−45%，但 knapsack-500 回退
  ×1.2-3.0、setcover-100x150 轻微回退），故入口按 `detect_pseudocost_
  favorable` 类检测（与 R3 分支规则同一检测器：全 Binary + 正系数 Le 行 +
  行数 2..=64）选栈：命中多维背包类 ⇒ global，单行背包/等式行/Ge 行/
  连续列类 ⇒ local；显式 `local`/`global` 覆盖优先（`threads=1` 顺序路径
  不读该字段，t1 行为逐字节不变）。
- **小算例调度门（W1，`--threads ≥ 2` 生效，仅 `--work-stack auto` 默认
  决策链）**：小树并行 = 纯开销（P1 插桩取证：milp-mixed t4 下 root 单线程
  段 ~42%、Σcv_wait 57.29%、3/4 worker 饥饿、树 ≤ 85 节点填不满线程）。
  threads ≥ 2 且问题规模低于静态阈（LP 单位 = 变量数 × (行数+1) ≤ 65,536）
  且未设紧节点预算（node_limit = 0 或 > 512）时，先以 512 节点预算跑
  **顺序路径**探测：预算内出终态 ⇒ 门命中直返——**threads ≥ 2 的输出与
  threads=1 逐字节相同**（同一条确定性代码路径），免掉池协调与 worker
  冷启动开销；预算耗尽 ⇒ 树不小，丢弃探测、并行从零派发（既有行为 +
  有界前缀浪费，mknap-5x300 实测 ≈ 0.2s/30s 量级）。512 预算落在既有
  t1 节点数的两簇空档（并行纯亏簇 ≤ 360：milp-mixed 53-85、knapsack-500
  360、knapsack-2000 38、knapsack-5000 37、setcover-100x150 13-15；
  并行有收益/大树簇 ≥ 2,120：mknap-3x100）。显式 `--work-stack
  local|global` 点名并行引擎的栈形态，**不进门**（P2 式栈形态 A/B 保持
  可测量，覆盖优先）。
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
- **稀疏混合表格（大规模轮 R1）**：表格单元数（m×ncols）超过
  `RC_SPARSE_CELL_THRESHOLD = 2^16 = 65536`（W7 自 2^18 下移——assign-40x40
  实测 +32.6% 区间不相交、其余类行为逐字节零扰动；lp-sparse/setcover-400x800
  在 2^18 时代已启用，门槛下移对它们是零操作）才启用每行非零列索引簿记
  （`row_nz`）：主元消元同步归并游标簿记（fill-in 插入 / 恰零惰性墓碑 /
  墓碑复活，每 64 主元压缩 + debug 全量重建比对），密度过高的行（nnz×4 ≥
  ncols）回退稠密标记；检验数批量散射（`#[inline(never)]` 外置函数，防
  codegen 扰动旧路径）、加行列插入索引右移、割分离稀疏读表
  （`tableau_row_terms`）均走稀疏路径。**阈值之下逐字节走旧稠密路径**
  （`do_pivot_legacy_inner` 与旧实现函数体经程序化比对字节一致），bench
  13 个既有算例全部留驻旧路径、行为零扰动。实测 setcover-400x800 同一
  24000 迭代序列每迭代成本 −6%（两轮独立复现）；lp-sparse-4000x800
  收益 0~−5%（不稳定）。
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
  连续退化迭代 > 50 次或总迭代超过 3×(m+n) 后切入 Bland 模式：入基取合格列
  最小下标；出基侧全部并列候选（被阻挡基变量、入基变量自身界翻转）一律取
  最小下标。**Bland 释放（大规模轮 R8）**：退化停滞触发的 Bland 在非退化步
  （目标真改善 ⇒ 已脱离触发退化区）且驻留满 8 步后回 Dantzig 定价——Bland
  的最小下标入基在宽表（ncols ≫ m）上单步进展极小，纯 Bland 爬行实测打满
  迭代上限（setcover 类根 LP）；每次调用最多释放 64 次，重入预算耗尽后
  Bland 永久驻留，无循环定理的终止路径恒在场。比值测试同时考虑基变量上下
  界与入基变量自身界翻转；|主元| ≥ tol.pivot 才允许；入基列无阻挡 ⇒
  Unbounded。检验数每次由当前表格重算；每 200 次主元用原始数据重算 B⁻¹b
  抑制漂移；迭代上限 `max(10_000, 20×(m+n))`（割循环内每轮 resolve 收紧为
  `max(3_000, 12×(m+n))`，见「根节点 GMI 割」一节）⇒ IterationLimit
  （诚实上报）。

### 双引擎架构（--lp-engine，大规模轮 S2–S4、S6）

LP 求解器有两个实现，是**同一 LP 语义的两个实现**（比值测试、反循环、
两阶段纪律的契约逐条写在 `simplex_sparse.rs` 模块头注释，改任一侧都要
检查另一侧是否需要对齐）：

- **`tableau`（`simplex.rs`）**：稠密表格两阶段有界单纯形（上节），历史
  默认路径；小表格式问题上每次调用底价更低。
- **`revised`（`simplex_sparse.rs`）**：稀疏修正单纯形——约束矩阵 CSC
  压缩列存储 + Gilbert-Peierls 左视稀疏 LU（阈值主元 |v| ≥ 0.1·max|v|，
  基三角形部分先行列序抑制填充）+ ETA 乘子形式基更新 + 周期整基重因子化
  并从原始数据重算 x_B（漂移双保险）+ 界变化后的对偶修复热启动（对偶
  可行性保持时数步旋转恢复原始可行，异常回退整基重建 + Phase-1）。7 项
  共享常量与 Bland/W4/R8 触发语义与稠密引擎逐款对齐（S5 契约测试钉值）。
  **定价规模分档（S6）**：构造期按原始 m×(n+m) ≥ 10^6 格（与 auto 定档
  同阈同维）定死——阈上 = revised 的 auto 默认域 ⇒ **Devex** 定价
  （Forrest–Goldfarb 参考框架近似最陡边：入基打分 rc²/w_j，权重随主元
  更新、任一越阈 1e6 即重置参考框架，solve/resolve 入口与 Phase-1→2
  切换亦重置）；阈下保持 **Dantzig**、与本引擎 Devex 化之前逐字节同
  入基序列（A/B 兼容锚——阈下中规模 B&B 类 setcover-300x600 上 Devex
  顶点轨迹劣化搜索形态，探针实测见 `simplex_sparse.rs` 模块头定价节）。
  稠密引擎恒 Dantzig。**两引擎对账红线 = 目标恒同**；迭代序列不承诺
  逐次一致（Dantzig 期的「逐次一致 7774」是同定价规则下的偶然副产品，
  Devex 域内有意分叉——多解顶点下 x 允许不同）。

两实现经 `LpEngine` trait（15 方法 = branch_bound/cuts/heuristics 全
消费面）统一接入 B&B（持有 `Box<dyn LpEngine>`；W6 局部割快照走 trait
`snapshot()`）。`auto`（S4 起默认）按表格规模定档：m×(n+m) ≥
`ENGINE_AUTO_CELL_THRESHOLD` = 10^6 格 ⇒ revised，阈下 tableau——阈值
落在 setcover-400x800（480k）与 lp-sparse-4000x800（3.84M）之间的整数
量级空带（「表格超百万格」一句话可解释）；显式 `tableau`/`revised`
覆盖优先。`revised` 仅支持 `--threads 1`（并行 worker 直接构造稠密
Simplex，旗标与 ≥2 组合显式报错，不静默分叉）。

实测（S2–S5 目标恒同为前提；迭代序列 Dantzig 期曾与稠密引擎逐次一致
7774 次，S6 Devex 起有意分叉——对账红线为目标恒同）：lp-sparse-4000x800
standalone 与全路径均约 3.5×（dense 32–35s → sparse 9.0–10.2s 中位）；
S6 Devex 定价后同窗 A/B（3 rep min）再降至 6.59s（旧 11.27s，0.58×）、
lp-wide-1000x2000 8.84→7.56s（0.85×），目标极差 0.00；默认 auto 生效后
该例 gate t1 从 30.3s limit_reached（W2 内部时限截断）翻为约 9s optimal
（gate 基线随轮重钉）。阈下类留 tableau 有应用层证据（acceptance 默认
输出逐字节基线不变）。诚实降档：setcover-300x600 类（270k 格、8% 密度
内核）revised 曾 1.37×（<2× 线，S4 同窗）——LU 填充瓶颈（重构 ~4.5ms/
次微基准留档），需 Markowitz 级列序工程，该类默认保持 tableau（S6 窗
同窗比值 2.83×、窗口间分歧，翻档未执行、留专门轮次定档）。

### 深度优先分支定界（WP4）

- 单一 `Simplex` 实例 + `Bounds` 界栈（trail）+ 深度优先，节点零拷贝；下行
  分支先探索；回溯严格按 trail 逆序恢复旧界，**非基双向 snap**（收紧越界
  snap / 恢复后严格内点 snap 到较近界、等距按 rc 符号选侧）处理后热启动。
- 节点处理顺序固定：LP Infeasible ⇒ 剪枝；根节点 LP Unbounded ⇒ 整树判
  Unbounded（见理论近似声明）；IterationLimit ⇒ 剪枝并记 limit_hit；LP bound
  不优于 incumbent ⇒ 剪枝；整数性满足 ⇒ 取整复核（伪整数解**不得**剪枝，
  必须继续分支）；空界 ⇒ 剪枝。
- **根 LP IterationLimit 的 incumbent 提取（大规模轮 R8）**：根 solve 以
  迭代上限收场时，Phase-2 中断点的基是原始可行的（比值测试每步保持 rhs
  在界内、Phase-1 已完成）——取整复核（原始问题逐约束 + 整数性，不依赖
  表格语义）通过则作为 incumbent 随 LimitReached 携带返回（终态仍诚实
  limit_reached，只供应不改语义）；Infeasible/Unbounded 不提取，提取失败
  如实保持无 incumbent（如 setcover-300x600/400x800 的分数点取整违反覆盖
  约束）。
- 分支变量规则默认**类自适应**（大规模轮 R3，`BbOptions::branch_rule_auto`
  默认 true）：`detect_pseudocost_favorable`（全变量 Binary ∧ 全部约束为
  正系数 Le 行 ∧ 行数 ∈ [2,64]——覆盖多维背包类；P3 收紧行数下限 2，
  单行背包类回落）命中 ⇒
  [`BranchRule::Pseudocost`] + 可靠性冷启动（根节点全量探测，预算
  2×min(分数候选数, 64)；深度 1-8 每节点探测预算 8；深度 > 8 不探测；
  探测 = 收紧界 → 重解 → 严格正增益才记账 → 精确恢复原界，确定性不读
  墙钟），否则 MostFractional（max min(f, 1−f)，并列取最小下标）；显式
  指定（CLI `--branch-rule` / `branch_rule_auto = false`）覆盖优先。实测
  mknap 类 nodes −44%~−56%（mknap-5x300 无时限 341446→163641）；单行
  背包类探测白烧已由 P3 收紧消除（行数下限 1→2：knapsack-500 墙钟
  −48%、knapsack-5000 −26%、knapsack-2000 持平，objective 恒最优），
  不命中类（setcover/marketsplit/milp-mixed）行为逐字节不变。剪枝条件
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
  节点/时间上限触发 ⇒ LimitReached（附 incumbent）。时间上限自 W2 起
  在 LP 迭代循环内部每 64 次主元条件读钟（`LpStatus::TimeLimit`，独立于
  迭代上限状态码——不触发 R8 割回滚）；时限关闭时零墙钟读点。
- **基准套件**：`bench/`（`gen_instances.py` 生成，种子确定）：强相关背包
  200/500/2000/5000（附 DP 独立核对最优 `.opt`）、多背包 3×100/5×300/
  10×200、指派 16×16/40×40、混合 MILP、宽 LP、稀疏 LP 4000×800、
  setcover 至 400×800、marketsplit 至 5×50，共 22 个 `.milp` 算例。计时
  矩阵分两集（`bench/run_timings.py --set gate|full`）：gate = 6 算例 ×
  30s 时限（单次实测约 2 分钟、设计上限 3 分钟，快速回归）；full = 18 算例 × threads {1,4} ×
  90s（中位数 JSON + 逐次原始记录）。2026-10 大规模轮（R1–R8 保留 + R9 空集
  回退）逐轮同窗交错 A/B 实测与数据见 `docs/large-scale-report.md`。

### 根节点 GMI 割（WP3）

只在根节点生成（对原问题全体整数可行解全局有效；分支子树内调用无效，
禁止）。取基变量为整数类型结构变量的行推导，割为 `Σ π_j·d_j ≥ 1`（π_j ≥ 0）。
**整数列前提**：整数类型结构列仅当所贴界为整数值（|bound−round| ≤
tol.zero）时才可用整数公式，否则必须按连续列公式处理，否则割可切掉整数
最优而报 Optimal。数值纪律：|c_j| < tol.zero 的项只允许按放宽割的方向
丢弃，严禁加严；violation ≥ min_violation 才收；动态范围超 1e9 整条丢弃。
S7 起全部割族另在施加边界统一过数值卫生闸（清洁化/密度闸/violation
复核，见「分离边界数值卫生闸」小节）。
割循环后根 LP 不可行是整数不可行问题的合法终局（割对全部整数可行解有效、
该集合为空时割可切空松弛区域）——终态如实映射 Infeasible（CLI
`status = infeasible`、退出码 2），不是数值错误。

#### 表格行聚合 MIR 割（大规模轮 R2）

根割循环的割族在 c-MIR、GMI 之后新增**聚合 MIR**（`cuts::aggregated_mir_cuts`）：
候选表格行（GMI 同口径）与其**非基松弛列的原始归一行**合并（μ = −α/σ 恰消去
该松弛列、并入其结构内容；行支撑含锚列自身时整行除锚系数重新归一，漏归一会
静默丢项产无效割——由 fuzz 抓出并钉住），全部列取单侧偏差（非基按所贴界、
基本结构列按较近界）后套与 GMI 相同的混合整数公式；**连续/松弛列 π ≠ 0
保留**（松弛列回代按 GMI 既有记账纪律）。每轮尝试上限 64（`AGG_MIR_MAX_TRIES`）、
配额 8 条（`AGG_MIR_MAX_PER_ROUND`，无配额时 setcover 类根割白烧 +244%），
批序 c-MIR → GMI → 聚合 MIR → 隐含界 → 覆盖统一截断。两个朴素构造（表格行
内松弛列全部回代 ⇒ 恒等式塌缩 0=0；表格行对 λ 直接合并产生双基列）经离线
证明否决。正确性由 50 种子随机小模型穷举 fuzz 钉死（3072 整数点 × 连续两端点，
零误切，随测试套件运行）。实测：mknap-5x300 无时限 nodes −42%
（591752→341446）、墙钟 −39%；knapsack-2000 nodes −94%（594→38）；
milp-wolsey-2 根割界直接收紧到整数最优 40 ⇒ 根节点即证最优（nodes 7→1）。

#### 树内局部割（大规模轮 W6）

割不再只在根：**多维背包类**（`detect_pseudocost_favorable`，与 R3/W1 同一
检测器）的顺序搜索在深度 ≤ 8 的节点上按配额跑轻量局部割轮（复用根割批序
c-MIR→GMI→隐含界→覆盖，批序产出统一过 S7 分离边界卫生闸后再按配额截断；
每 64 个 resolve 一次、每节点 ≤ 2 条、全局 2,000 条封顶）。局部割只作用于
该子树——施加前克隆 Simplex（derive(Clone) 全字段
纯数据，深拷贝即精确状态还原），DFS 退出子树时恢复（remove_row 在表格基
语义下不可行；并行 worker 不施加——窃取交错会泄漏到子树外切真最优）。
`BbOptions::local_cuts = false` 逐字节回到关割路径。实测 mknap-5x300 无时限
t1 节点 −17.0%（87,800→72,891）、迭代 −15.1%，目标恒同 4762；其余类经类
门控零扰动（无门控版 marketsplit-4x40 零压缩纯开销 +21% 墙钟，实测否决）。

#### 割循环预算与回滚（大规模轮 R8）

割循环每轮 resolve 是暖启动（对偶修复快路径，收敛轮实测 ~1.5k 迭代），
经 `set_iter_limit_hint` 收紧为 `max(3_000, 12×(m+n))`（resolve 后立即
复位，不泄漏）；新割行引发的退化爬行超限即走**回滚语义**——弃本轮增量割
（ledger/recorded/applied 截断回轮前）、冷重建 Simplex 重放已提交割（确定性
路径 ⇒ 恢复上轮末 Optimal 根态），已提交割的根界保留、树搜索继续，不再以
limit_reached 无解终态整树回退。配合单纯形的 Bland 释放（见「完整 Bland
反循环规则」），实测解开 setcover-100x150 的根割循环白烧：lp_iterations
18924→5246（−72%）、树 45→13 节点、割 0→23 条存活，目标 276 恒证。
Infeasible 信号不受回滚影响（Phase-1 判定在预算内第一时间如实返回；真实
不可行最迟由树搜索发现，不会被伪造）。

#### 分离边界数值卫生闸（大规模轮 S7）

割在**施加边界**统一过三道闸后才进表与记账（`cuts::sanitize_cut`，私有
纯函数）：施加点为根割循环 `cut_loop_recorded` 与树内局部割
`separate_local_cuts`——全部五族割（c-MIR/GMI/聚合 MIR/隐含界/覆盖）与
稠密/稀疏两引擎的唯一汇聚点，单点过闸即全覆盖。三道闸依次为：① 小系数
清洁化（|c_j| ≤ 1e-9 且该列界幅有限 ⇒ 丢项并按界幅抬高右端补偿——数学上
**只弱化**，清洁后可行域只扩不缩，整数可行点绝不被切；界幅无界的列原样
保留）；② 密度闸（清洁后 nnz > max(200, ⌈0.1·n⌉) 整条丢弃）；③ violation
复核（右端补偿吃掉分离量至 min_violation 以下的割丢弃）。过闸发生在每轮
配额截断（max_cuts_per_round / max_cuts）之前——ledger/recorded 记账与
worker 重放的即过闸后的行。丢割永远安全（只丢界强度）。动机与实测（Devex
基 160 条 74% 密度噪声割毒化全树的根因链；mknap-5x300 −25% 墙钟、
knapsack-2000/mknap 节点数变化与目标恒同对账）见
`docs/large-scale-report.md` S7 节。

### 预处理（WP3；W5 增 R5 等价列聚合）

保守规则五类且可证明安全，每轮严格按 **R1→R2→R3→R4→R5** 顺序施加，禁止
重排，最多 5 轮无化简即停：R1 整数变量界取整；R2 固定列（单点界）代入删列；
R3 空行判定/删行；R4 单系数行换算为变量界后删行（单例行消去的既定形态：
换算成界后点位由次轮 R2 完成变量替换，级联由轮循环驱动）；**R5（W5）等价
二元列聚合**——全部存活行中系数向量相同**且目标系数相等**的 m 个 Binary
列聚合为单个 Integer 变量 x'∈[0,m]（行/目标只看组内总和，与组内取值组合
双射 ⇒ 最优精确保留；postsolve 把 x'=t 还原为组内前 t 个成员置 1）。实测
knapsack-2000 聚合 1,596 列（2000→404）⇒ t1 墙钟 0.241→0.011s（−95%）、
lpi 16,787→2,247；knapsack-5000 聚合 4,590 列 ⇒ 1.606→0.021s（−99%）；
目标恒同 DP 真值。**红线记档（实测否决的规则）**：①不同目标系数的重复列
做支配固定（x_k←0）——最优解可同时取多个同权物品（x_j=x_k=1），交换论证
需要 x_j+x_k≤1 才成立，首版 A/B 背包族目标全面错误（3409→1415 等）已删；
②R6 多行隐含界传播 + R7 冗余行删除——全 bench 零树压缩（milp-mixed 节点
47→47/53→53）且 lp-sparse 确定性回退（墙钟中位 34.3→40.2s、lpi 7,774→
8,535，隐含连续上界扰动 LP 路径），实现并 A/B 后删除（数字见
docs/parallel-profiling.md §14）。

### 数值缩放（大规模轮 R4）

`scaling::scale_problem` 在求解编排最前端做问题级行列缩放（管线变为
scale → presolve → solve → postsolve → unscale；CLI `--scaling=off` 关闭，
关闭时逐字节复现旧路径）：迭代几何均值（≤ 8 轮 sweep），行因子
r_i = 1/√(max|a|·min|a|)（行非零上）吸附最近 **2 的幂**后整行乘上，列乘子
同式经代换 ŷ = x·s_j 施加；全部因子吸附 2^0 即停。代换保证**目标值不变**
（ĉ·ŷ = c·x）、正缩放保序（incumbent 比较合法）、线性双射保可行集（割在
缩放空间生成对原问题有效）；2 的幂乘除 f64 精确无舍入 ⇒ unscale 逐位还原，
终检与目标重算在**原问题**上做（输出契约不变）。**正确性红线（实测抓出）**：
列缩放仅 **Continuous** 列——Integer/Binary 列的整格双射 x∈Z ⟺ ŷ∈Z 对
2 的幂因子仅 s=1 成立（首次按宽口径落地被安全网终检在 milp-mixed-250i40c
上当场拒绝半整数解），Binary 的 0/1 语义同禁；全整数算例只有行缩放参与。
实测 mknap-5x300 无时限 nodes −46%（163641→87800）、墙钟 −41%
（32.7→19.2s 中位）；milp-mixed 连续列收益 lpi −18%/−21%；其余算例行为
逐位不变（knapsack DP 真值 3409/8975/35328/88526 在 on/off 双口径下精确）。

### 启发式（WP3）

`round_and_check`（简单取整 + 可行性复核）与 `dive`（fix-and-dive，迭代
≤ 8 次固定最接近整数的非整变量后重解 LP）只供应 incumbent 以加速剪枝，
不影响正确性；复核容差与 integrality 联动（见下）。

`subset_sum_mitm`（子集和 MITM，大规模轮 R5；CLI `--mitm=on|off`，默认
on）是 LP 驱动启发式全部失败后的**最后手段**，在根阶段 dive 之后、仅
incumbent 为空时启动。类门控（全部满足才启动，否则恒零成本跳过）：全变量
Binary ∧ 全约束 Eq ∧ 系数/右端可整数化（2 的幂分母有理数，|v| < 32000——
R4 缩放产生的 2^k 分母可精确恢复）∧ 行数 ≤ 8 ∧ n ≤ 84。机制（脱离 LP
信息源）：n ≤ 42 时全变量半枚举 MITM（2×2^(n/2)，u128 打包键 + 二分）
**精确**求解 A·x = b（planted 解必被找到）；n ≤ 84 时窗口半枚举（基点取
根 LP 解取整，|窗口| ≤ 42 的固定种子 LCG 随机子集上重解，默认 8 窗重试）。
产出经 `round_and_check` 原始界复核，只供应 incumbent 不影响正确性；固定
种子、不读墙钟（确定性）。实测解开 marketsplit 大规格「零可行解解锁死」：
marketsplit-4x40 从 200 万节点零 incumbent（25.2s）变为 854,952 节点证最优
19.0（9.9s）；marketsplit-5x50 拿到首个 incumbent 25.0（证最优未闭合）。

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
| `hard-knapsack-12` | 0-1 背包（12 变量，权重/价值/容量 7..47 / 11..144 / 172） | Optimal，objective = 523 | 159 节点（threads=1 默认——单行背包类 P3 起不命中类检测、回落 `most_fractional`；`--branch-rule pseudocost` 下为 196；并行下节点数可变） | 2^12 = 4096 穷举；实测节点数 ≫ 原算例（≤7），作为 `--threads 2/4` 真实压满节点池、`--log-secs` 周期行与 `--heuristic-mode` 的并发/受控参数验收载体 |
| `infeasible-tiny` | 状态型 | Infeasible，退出码 2，无 objective 行 | — | x≥3 ∧ x≤1 空集；正右端 ≥ 行（x≥3 在 lb=0 处违反，走人工变量路径）端到端覆盖 |
| `unbounded-tiny` | 状态型 | Unbounded，退出码 3，无 objective 行 | — | max x，x≥3 且上界 +∞；默认上界 +∞ 的合法性由本例顺带验收 |

## 大规模算例（bench/）

`bench/instances/` 共 22 个 `.milp` 算例（`bench/gen_instances.py` 生成，
种子确定），2026-10 大规模轮新增 9 个，其中进入计时矩阵的 5 个（t1 + 90s
探针标定）：

| 算例 | 规模 | t1 实测（探针标定） |
|------|------|---------------------|
| `lp-sparse-4000x800` | 稀疏纯 LP：800 行 × 4000 列，密度 8%（256,000 非零元） | **S4 起默认 auto ⇒ revised、S6 Devex 定价后同窗 A/B 约 6.6s optimal**（obj 5893.705757 恒同；迭代序列 S6 起与稠密引擎有意分叉，7774 为 Dantzig 期口径）——稀疏修正单纯形引擎主秀场（tableau 约 32–36s）；历史口径：tableau 35.8s optimal、W2 内部时限后 30s 名义时限下 30.1–30.3s 诚实 limit_reached——S4 默认翻 auto 后 30s 时限不再截断、S6 Devex 再降迭代数，gate 基线随轮重钉 |
| `setcover-400x800` | 集合覆盖：400 覆盖行 × 800 二元列，密度 0.08 | **W4 解锁根 LP**：曾 17.3s limit_reached（1 节点，24000 迭代上限锁死、零 incumbent）；W4 给总迭代 Bland 触发加停滞证据门后根 LP 首次可解——90s 预算下诚实终态携带 incumbent 129（auto_degraded/limit_reached 视降级时机），30s 预算仍为根内诚实截断 |
| `knapsack-5000` | 强相关 0/1 背包 5000 变量 | W5 前实测 1.6s optimal；W5 等价列聚合后 0.021s optimal（聚合 4,590/5,000 列），objective 88526 与独立 DP 真值（`knapsack-5000.opt`）一致 |
| `mknap-10x200` | 多维背包：200 二元变量 × 10 约束行 | 26.4s optimal，obj 1791（254,119 节点 / 3,132,767 LP 迭代） |
| `assign-40x40` | 指派 40×40（1600 二元变量 / 80 等式行） | 0.12s optimal，obj 176（1 节点——指派约束矩阵全单模，根 LP 即整数解） |

矩阵外另有拉伸/消融算例：`lp-wide-1000x2000`、`setcover-300x600`（W4
起根 LP 解锁：无时限 t1 optimal 144.0@1,867 节点 / 275,471 迭代，~60-97s；
曾为 18000 迭代上限锁死零 incumbent）、`marketsplit-4x40` /
`marketsplit-5x50`（R5 MITM 目标算例）。计时
矩阵：`python bench/run_timings.py --set gate`（6 算例 × 30s，约 3 分钟快速
回归）或 `--set full`（18 算例 × threads {1,4} × 90s）。逐轮改动与实测数据
见 [docs/large-scale-report.md](../docs/large-scale-report.md)。

外部交叉验证（S5）：`bench/python_bench.py` 把 micpython（经 milp-ffi DLL，
吃到 Auto 默认引擎）与 CBC/HiGHS/SCIP 同机对比——14 算例 micpython 全部
optimal、13 个有共识对的算例与 HiGHS/SCIP 全部一致（纯 LP 三例位级吻合），
MISMATCH 5 条均为既有外部离群（scip knapsack-5000、cbc Unbounded×3）零新增。
复现（conda 环境 `mic-bench`，创建命令见脚本头注释）：

```bash
conda run --no-capture-output -n mic-bench python -u bench/python_bench.py
# 单算例核验（仅打印不落盘）：python -u bench/python_bench.py --pair <算例名> micpython
```

## 仓库布局

```text
milp/
├── src/
│   ├── main.rs          # CLI 入口（WP5 + 受控搜索参数：threads/log-secs/mipgap/heuristic-mode/branch-rule/scaling/mitm/work-stack/determinism）
│   ├── error.rs         # MilpError 与 Result 别名（WP1）
│   ├── model.rs         # 建模 + 集中容差（WP1）
│   ├── instance.rs      # .milp 解析/序列化（WP1）
│   ├── simplex.rs       # 两阶段有界单纯形 + 稀疏混合表格簿记（WP2，大规模轮 R1）
│   ├── simplex_sparse.rs # 稀疏修正单纯形引擎（CSC + 稀疏 LU + ETA，大规模轮 S2-S5；--lp-engine revised/auto）
│   ├── simplex_sparse.rs # 稀疏修正单纯形引擎 + LpEngine 双引擎接口与类自适应默认（大规模轮 S2–S4）
│   ├── presolve.rs      # 保守预处理与还原（WP3）
│   ├── cuts.rs          # 根节点 GMI/聚合 MIR 割（WP3，大规模轮 R2 扩展）
│   ├── heuristics.rs    # 取整/潜水启发式 + 子集和 MITM（WP3，大规模轮 R5）
│   ├── scaling.rs       # 行列迭代几何均值数值缩放（大规模轮 R4）
│   ├── branch_bound.rs  # 深度优先分支定界 + 受控入口 + 类自适应分支默认（WP4）
│   ├── solver.rs        # 统一编排（含缩放接线）、受控求解入口与安全网终检（WP4）
│   ├── search_ctrl.rs   # 搜索控制层：进度快照/gap 口径/自适应控制器（纯逻辑零锁）
│   └── parallel_bb.rs   # 并行分支定界引擎（节点池/incumbent/双 Condvar/粘滞 stop/全局工作栈）
├── instances/           # 8 个验收算例（7 原算例 + hard-knapsack-12）
└── tests/acceptance.rs  # 验收测试（WP5 + 受控搜索参数扩展）
```

单元测试 `#[cfg(test)]` 就近放在各模块内；验收统一在 `tests/acceptance.rs`。
