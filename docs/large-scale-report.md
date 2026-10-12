# milp 大规模优化轮交付报告（R1–R9）

> 分支 `feat/gurobi-large-scale`（`git diff main...HEAD` = 10 提交：bench 扩充
> `334e352` + 算法轮 `fc9f5ee..cac558e`，其间 `54271e4` 为本报告首版与终局
> 基准；31 个文件、+22828/−180 行，其中约 2.1 万行为新增 bench 算例与计时
> 记录）。本报告记录：背景差距轴 → 逐轮改动与实测数据（**保留与回退都写**）
> → 验证矩阵 → 遗留问题。
>
> **总裁决：R1–R8 八轮保留，R9 回退（空集——无代码改动，见 R9 节）。**
>
> **数据口径声明**：各轮收益数字均为实现者/裁决在**同机同窗交错 A/B 的中位数**
> （两版本二进制交替各跑 2–3 次），这是本分支唯一认可的墙钟口径。跨窗口的
> 墙钟不可直接比较——多轮已实证「摘要矩阵跨窗伪影」（见 R4/R5/R6/R8 各节
> 的驳回记录）。行为字段（status / objective / nodes / lp_iterations）在
> threads=1 下确定，跨窗可比。`bench/gate_runs.json` 为 R8 后（2026-10-12，
> 提交 `cac558e`）的 gate 集实测产物。

---

## 1. 背景：差距文档的轴

大规模优化轮按差距文档定义的轴推进；本分支主攻轴 1/2/4，另有两轮投轴 3
并行（R6 保留、R9 空集回退）（轴 4 的动机表述见 `milp/src/scaling.rs:3-7`
模块文档「动机（差距文档轴4）」）：

| 轴 | 主题 | 改动前现状（gate 基线 / 上轮实测） |
|----|------|-----------------------------------|
| **轴1 稀疏核心** | 大规模表格的每主元/每迭代成本 | dense 表格 O(ncols) 逐列扫描；lp-sparse-4000x800 t1 约 33–40s optimal（gate 基线 39.609s）；setcover-400x800 21.489s limit_reached（24000 迭代上限） |
| **轴2 树压缩** | 结构化算例的割产出与分支效率 | mknap-5x300 无时限 591,752 节点 / 47.7s（R1 后基线）；knapsack-2000 594 节点；marketsplit-4x40 两百万节点内零可行解（25.2s limit）——LP 驱动启发式跨不过 planted 解的汉明距离 |
| **轴3 并行** | 并行效率（Gurobi 式 t4/t1 ≥ 3× 目标） | 并行子树膨胀 + 单锁池竞争；R6 实测 t4/t1 ≈ 0.63–1.15×，并行无稳定优于单线程（诚实负结果，见 R6；R9 复投后空集回退，见 R9） |
| **轴4 数值** | 系数量级失衡 | 三个已修数值问题（对偶修复阈值、mir_floor 吸附、ratio_test 带宽）同源：固定容差在失衡量级上语义漂移——缺问题级缩放 |

**改动前 gate 基线**（t1 中位，`--time-limit-secs 30`，工作流记录）：

| 算例 | 墙钟 | 状态 | nodes | objective |
|------|------|------|-------|-----------|
| knapsack-200 | 0.003s | optimal | 1 | 3409 |
| setcover-100x150 | 0.361s | optimal | 45 | 276 |
| knapsack-2000 | 0.606s | optimal | 594 | 35328 |
| lp-sparse-4000x800 | 39.609s | optimal | 1 | 5893.705757 |
| setcover-400x800 | 21.489s | limit_reached | 1 | null |
| mknap-5x300 | 17.674s | auto_degraded | 166123 | —（30s 预算降级） |

**改动后同集**（`bench/gate_runs.json`，2026-10-12 = R8 后 `cac558e`，t1）：
knapsack-200 0.003s optimal nodes 1；setcover-100x150 **0.240s optimal
nodes 13**（R8：45→13）；knapsack-2000 0.265s optimal **nodes 39**；lp-sparse
34.9s optimal nodes 1；setcover-400x800 13.2s limit_reached（行为不变）；
mknap-5x300 13.6s **optimal** nodes 87800（30s 预算内证完，不再降级）。
墙钟跨窗只作参考，行为字段的确定变化（594→39、45→13、auto_degraded→
optimal）才是跨版本可比的证据。

---

## 2. 逐轮记录

R1–R8 八轮**保留**（裁决均为「收益压倒有界回退、正确性有 fuzz/穷举/逐字节
对账锚」；R7 为行为零变化的复核遗留修复轮）；R9 **回退**（空集裁决）。各轮
内部另有被否决/回退的子尝试，逐条记录如下。

### R1（fc9f5ee）稀疏混合表格：行级非零索引簿记 + 每行密度回退 —— 轴1

**改动**：仅 `milp/src/simplex.rs`（+约 950 行）与 `milp/src/cuts.rs`（稀疏
读表）。核心设计是**规模门**：表格单元数 m×ncols > `RC_SPARSE_CELL_THRESHOLD`
= 2^18 = 262144（`simplex.rs:143`）才启用每行非零列索引 `row_nz` 簿记——
主元消元归并游标同步（fill-in 插入 / 恰零惰性墓碑 / 墓碑复活，每 64 主元
压缩 + debug 全量重建比对）、密度回退（nnz×4 ≥ ncols 转稠密标记）、检验数
散射（`compute_rc_batch_scatter`，`#[inline(never)]` 外置）、加行列插入索引
右移、割分离稀疏读表 `tableau_row_terms`。**阈值之下逐字节走旧稠密路径**
（`do_pivot_legacy_inner` 与 HEAD 旧 `do_pivot` 函数体经程序化文本比对
**字节一致**）。阈值取 2^18 而非 65536：bench 当时 13 算例中表格最大的
lp-wide-200x300 为 300×500=150,000，65536 会使它启用簿记、违背「小算例
逐字节回基线」；2^18 使 13 算例全部留驻旧路径。

**过程回退（初版被裁决回退后重做）**：初版收益与回退同量级（小算例
+3.4~+5.3% 且区间几乎不重叠、mknap 节点率 −6.2% 双重佐证）——净相抵且
有害，回退。重做（本 commit）把簿记整体置于规模门之后，回退信号降至噪声
级（小算例 +1~2% 但 rep 区间大量重叠；mknap 节点率 −0.8%）。

**实测（同窗交错 A/B，vs 基线 worktree 探针）**：
- 收益：setcover-400x800 中位 14.04→13.20s（x0.94，**−6%**，3-rep 交错；
  两轮独立复现——同一 24000 迭代序列的每迭代成本真实下降）；lp-sparse-4000x800
  x0.997（今日持平；另一窗曾测 0.953–0.978，真实效应 0~−5% 且不稳）。
- 零扰动：8 算例行为对账 status/objective/nodes/lp_iterations/cuts/bound
  逐字段精确一致（含 lp-sparse-4000x800 的 7774 迭代、setcover-400x800 的
  24000 上限、mknap-5x300 无时限 4762.0/591752/7094278）。
- 过程发现（工程质量）：散射路径与旧循环共存同一函数体时，寄存器分配扰动
  使纯旧路径的 knapsack-2000 墙钟 +11%——散射体外置 `#[inline(never)]`
  后回基线。

**必须记录的数字夸大**：commit 消息与摘要矩阵声称 setcover-100x150
x0.83 / knapsack-2000 x0.81 / lp-sparse −12% 均**不可复现**（小算例走逐字节
旧代码路径，设计上不可能提速；『改动前』列与更早慢窗基线逐字相同、属跨窗
复用）。可辩护收益仅为 setcover-400x800 的 −6%（两轮复现）与 lp-sparse 的
0~−5%（不稳）。此后各轮汇报统一以同窗交错 A/B 中位数为唯一口径。

### R2（e27b488）表格行聚合 MIR 分离（O3 方案）—— 轴2

**改动**：`milp/src/cuts.rs`（+499：`aggregated_mir_cuts` 新割族 + 批序
接入 + fuzz 测试）；simplex.rs 未动；acceptance.rs 重钉 2 处。

**数学（先离线证明再落地）**：候选表格行（GMI 同口径）与其**非基松弛列的
原始归一行**合并（μ = −α/σ 恰消去该松弛列、并入其结构内容；行支撑含锚列
B1 自身时整行除锚系数**重新归一**——漏归一会静默丢项产无效割，被 fuzz
穷举抓出修复后钉作回归锚点）；全部列取单侧偏差（d ≥ 0 恒成立、整数列整界
时 d ∈ Z_+）后套 GMI 同款混合整数公式；**连续/松弛列 π ≠ 0 保留**，松弛列
回代按既有记账纪律。两个朴素构造被证明否决：① 表格行内松弛列全部回代 ⇒
恒等式塌缩 0=0；② 表格行对 λ 直接合并产生双基列，GMI 格论证失效。框架上限：
`AGG_MIR_MAX_TRIES` = 64 尝试 / `AGG_MIR_MAX_PER_ROUND` = 8 每轮配额
（无配额时 setcover 白烧 +244%），批序 c-MIR → GMI → 聚合 MIR → 隐含界 →
覆盖统一截断。

**实测（vs fc9f5ee 基线探针，交错中位）**：
- mknap-5x300 无时限：墙钟 47.7→29.0s（**−39%**）、nodes 591752→341446
  （**−42%**）、lp_iterations 7.09M→3.97M（−44%）、cuts 42→62（聚合 MIR
  +20 条）、最优 4762.0 逐位一致。
- 意外外溢：knapsack-2000 x0.544/x0.518（两序区间完全不相交），nodes
  594→38（**−94%**）、cuts 12→15、最优 35328.0 一致。
- milp-wolsey-2 根割界直接收紧到整数最优 40 ⇒ **根节点即证最优**（nodes
  7→1）——acceptance 两处重钉（nodes 7→1、进度首行 incumbent 37→40；
  objective 断言 40 未动，重钉理由入注释）。
- 诚实权衡：mknap-5x300 根界 4773.72→4774.32（**未突破**冻结值——批额内
  聚合割挤占 GMI/覆盖割席位，根界轻微变差但树大幅压缩）；mknap-3x100
  nodes −29% 但墙钟 +17%（小算例分离开销）。

**回退清单（全部实测定性，随轮保留）**：
1. setcover-100x150：x1.402/1.335（两序区间不相交），0.29→0.41s（+0.12s
   绝对）；结果字段全同（optimal 276 / nodes 45 / cuts=0）——根割循环既有
   失败模式（Bland 爬行）下聚合分离加深了失败前白烧，有界且绝对代价小
   （该失败模式由 R8 精确修复，见 R8 节）。
2. mknap-5x300 30s 截断态 x1.34 + auto_degraded incumbent 4751→4749
   （−0.04%，均为非终态启发式解）——自适应预算调度与更长根割阶段的交互
   伪影；无时限真实口径严格更优。
3. setcover-400x800 声称 x1.295 未获证实（6-rep 交错 x1.057 全重叠）；
   lp-sparse 声称 x1.222 被驳回（纯 LP 不进割循环，迭代数两侧恒同）。

**本轮被否决/回退的子尝试**：
- setcover 根割热启动修复（phase1_pending_fresh 在当前热基上直跑 Phase-1）：
  实现并实测后**净回退**（knapsack-500 +114%、setcover +44% vs mknap-3x100
  −23%），按「每算例净墙钟无回退证据」标准回退。诊断结论：setcover-100x150
  第 3 轮 Phase-2 本身不收敛（Dantzig 停滞 → Bland 爬行，200k 迭代实验仍
  回退），解锁需要定价/退化处理（另一工作项），非迭代上限/热启动问题——
  割循环迭代上限**不动**（R8 以完全不同的机制解决，见下）。
- 可行性方向（轴2 前轮提示）未做——留待 R5 以完全不同的机制解决。

### R3（bdf3480）可靠性分支类自适应默认 —— 轴2

**改动**：`branch_bound.rs`（+201）、`main.rs`、`examples/profile_stats.rs`、
`acceptance.rs`（1 重钉）。parallel_bb.rs 未动（入口一次解析后 worker 收到
即具体规则，重放契约不变）。

**裁决先行**：旧注释「Pseudocost 作默认使 marketsplit-4x30 倒退 ~11×」与
baseline.json 隔离实验（−45%）矛盾。定档矩阵 {mostfrac, pseudo×probes=0,
pseudo×probes=8} × 5 算例 × 3 rep 裁决：**11× 倒退 = pseudo×8（带探测）
口径**（0.108→1.039s ≈ 9.6×）；纯被动 pseudo 在 marketsplit 上反而最优
（0.108→0.050s，−54%）——**探测预算在等式行结构上白烧才是倒退根因**，ask
的「被动 pseudocost 致倒退」假设被否定。

**落地**：`detect_pseudocost_favorable`（`branch_bound.rs:297`）——全变量
Binary ∧ 全部约束正系数 Le 行 ∧ 行数 ∈ [1,64]（覆盖 knapsack-*/mknap-*；
setcover 的 Ge 行 / marketsplit 的 Eq 行 / milp-mixed 的连续列均不命中）；
`BbOptions.branch_rule_auto` 默认 true，求解入口一次解析（命中 ⇒ Pseudocost
+ probes ≥ 8；未命中 ⇒ MostFractional；显式指定覆盖优先、幂等）；顺序/
并行两路共用。CLI `--branch-rule` 显式设规则 + 关 auto，缺省 = 类自适应。
探测日程：根节点全量冷启动（预算 2×min(分数候选数, 64)）、深度 1-8 每节点
预算 8（`pc_reliability_probes`）、深度 > 8 不探测；确定性不读墙钟。
bench 探针 `profile_stats` 同步默认，并修复误留的 `--pc-probes` 拒绝守卫。

**实测（vs e27b488 基线探针，默认路径 3-rep 中位）**：
- mknap-5x300 无时限：墙钟 x0.664（45.1→30.0s）、nodes 341446→163641
  （**−52%**）、最优 4762.0 逐位一致——叠加 R2 后累计树压缩 591752→163641
  （**−72%**）。
- mknap-10x200：x0.659（37.4→24.7s）、nodes 254247→143503（−43.5%）、
  obj 1791.0 一致；mknap-3x100：x0.534、nodes 4772→2106（−55.9%）。
- 不命中类零扰动（行为逐字节一致）：setcover-100x150（45/18924）、
  marketsplit-4x30（12132/119691）、milp-mixed-250i40c（53/1874）、
  lp-sparse（7774）、setcover-400x800（24000）——据此驳回摘要矩阵两格
  跨窗伪影（lp-sparse x1.398 / setcover-400x800 x0.666）。

**回退清单（随轮保留，实测定量）**：knapsack-500 x1.759（0.024→0.042s，
**绝对 +18ms**）、knapsack-5000 x1.266（2.32→2.94s，**+0.62s**）、
knapsack-200/2000 持平——单行背包类在探测口径下的白烧，与 mknap 类
−20~−40s 相比净绝对收益压倒（检测规格明言覆盖 knapsack-*/mknap-*，按规格
落地并如实报告；若要零回退可收紧为 rows ≥ 2，未擅改规格）。

**重钉**：acceptance hard-knapsack-12 nodes 159→196（命中背包类默认翻
Pseudocost+probes8；objective 523 恒为已知最优，演进史入注释；显式
`--branch-rule most_fractional` 复现旧径）。

### R4（51447da）行列迭代几何均值 scaling（2 的幂系数）+ 问题级全链路往返 —— 轴4

**改动**：新增 `milp/src/scaling.rs`（+约 470 行含测试）；接线 `solver.rs`
（scale → presolve → solve → postsolve → unscale，`SolveOptions.use_scaling`
默认 true）与 `main.rs`（`--scaling=on|off`，等号/空格两式）；presolve 与
acceptance 未动（零重钉）。

**算法**：≤ 8 轮 sweep，行因子 r_i = 1/√(max|a|·min|a|)（行非零上）吸附
最近 2 的幂后整行乘上；列乘子同式经代换 ŷ = x·s_j 施加；全部因子吸附 2^0
即停。目标**值**不变（ĉ·ŷ = c·x）、正缩放保序（incumbent 比较合法）、线性
双射保可行集（割在缩放空间生成对原问题有效）；2 的幂乘除 f64 精确 ⇒
unscale 逐位还原；终检与目标重算在**原问题**上做（输出契约不变）。

**正确性红线（实测抓出，ask 规格被证伪后收紧）**：ask 的「Integer 列
s_j = 2^k, k ≥ 0」不健全——缩放问题按 VarType::Integer 在 ŷ∈Z 上求解，奇
ŷ 反解得半整数；首次落地即被安全网终检在 milp-mixed-250i40c 上当场拒绝
（变量终值 0.5，Err 而非错解）。列缩放收紧为**仅 Continuous 列**（Binary
0/1 语义同禁），文档化并测试钉住。全整数算例只有行缩放参与（收益折损，
如实报告）。

**实测（vs bdf3480 基线探针，交错中位）**：
- mknap-5x300 无时限：墙钟 x0.586（中位 32.7→19.2s，2-rep 区间不相交）、
  nodes 163641→87800（**−46%**）、迭代 3.17M→1.75M（−45%）、cuts 62→70、
  最优 4762.0 逐位一致。
- milp-mixed-120i30c lpi 1117→917（−18%）、milp-mixed-250i40c 1874→1472
  （−21%）——连续列缩放的直接收益；mknap-10x200 nodes −0.9% 微优。
- 无真实回退：knapsack-500/2000/5000、setcover-100x150、setcover-400x800、
  marketsplit-4x30、lp-sparse 行为字段全部逐位不变（墙钟 x0.95–1.02 噪声
  带）。**摘要矩阵四个「显著回退/提速」全部驳回为跨窗伪影**（行为逐位相同
  的算例墙钟不可能系统性 x1.5–2.2 变化；mknap-5x300 矩阵『x1.28』系混合
  口径——基线取 30s 截断快窗值、新值取无时限）。
- 生效面诊断（测试钉住）：2 的幂行缩放对纯 LP 主元路径是**精确无操作**
  （与 Simplex 构建期行归一逐位相消）；mknap 树 −46% 来自割生成绝对门在
  缩放系数上的树混沌有利方向，非系统性容差收益（如实报告）；setcover
  系数全 ±1 ⇒ 全因子 1（幂等）。
- CLI `--scaling=off` 四算例与基线探针 status/nodes/objective 逐项相等；
  knapsack DP 真值 3409/8975/35328/88526 在 on/off 双口径下精确。

### R5（4aa081a）子集和 MITM：解开 marketsplit 大规格零可行解锁死 —— 轴2

> ⚠️ 记录勘误：**commit 标题为「可行性泵（界固定-失败解锁式）」，实际内容
> 是 `subset_sum_mitm`（meet-in-the-middle）**，代码中不存在任何泵实现
> （复核 grep 仅命中注释）。可行性泵方案在上轮已被判定死因（LCG 重启仍锚定
> LP 顶点邻域）而放弃，本轮按提示的「正确方向」换成脱离 LP 信息源的子集和
> 专用启发式——内容正确，标题误导审计线索（复核低危问题 #2，R7 已在 §4a
> 对照更正）。

**改动**：`heuristics.rs`（+324：`subset_sum_mitm` + `pack_sums/gcd/
pow2_denom` 内核 + 1 测试）、`branch_bound.rs`（+37：两处根阶段接线 +
`BbOptions.mitm_windows`）、`parallel_bb.rs`（+22：worker 根阶段）、
`main.rs`（`--mitm=on|off` 默认 on）；重生成
`bench/instances/marketsplit-4x40/5x50.milp`。

**机制**：类门控 = 全变量 Binary ∧ 全约束 Eq ∧ 系数/右端为 2 的幂分母
有理数（R4 缩放产生的 2^k 分母逐行乘最小 2 的幂**精确恢复**整系数）∧ m ≤ 8
∧ n ≤ 84（`MITM_MAX_VARS`，`heuristics.rs:233`）；逐行 gcd 归一（缩放因子
至 2^7 会爆 i16 半和打包）+ 打包安全检查。n ≤ 42：全变量半枚举 MITM
（2×2^(n/2)，u128 打包键全序排序 + 二分，确定性）**精确**求解 A·x = b——
planted 解必被找到；n ≤ 84：窗口半枚举（基点 = round(根 LP 解)，|窗口| ≤ 42
的 LCG 固定种子随机子集上重解；汉明距离 ~15 时单窗命中 ~7%，默认 8 窗）。
产出经 `round_and_check` 原始界复核，只供应 incumbent；接线在根阶段 dive
之后、仅 incumbent 为空时启动（可解类 dive 已命中则不烧；Le/Ge/连续列类
恒零成本跳过）。

**过程修复（两个接线 bug，均表现为「--scaling=off 有效、默认路径静默
失效」）**：① 缩放后 2 的幂分母被整性门拒（pow2_denom 恢复）；② 缩放因子
爆 i16 半和打包（gcd 归一）——两者都证明默认缩放路径被真实穿透。

**实测（vs 51447da 探针）**：
- marketsplit-4x40：基线 2,000,000 节点 limit_reached **零 incumbent**
  （25.2s）→ 新 854,952 节点 **optimal obj = 19.0**（9.9s，x0.394，迭代
  20.59M→6.36M）——解锁死解开且优于 planted 提示值（B&B 证 19 < 20）。
- marketsplit-5x50：基线零 incumbent → limit_reached + **首个 incumbent
  25.0**（迭代 28.76M→13.94M，−52%；墙钟 42.0→36.6s）——可行化达成，证
  最优未闭合（如实上报）。
- marketsplit-2x20：5,247 节点 → 1 节点根闭合（obj 10.0）；marketsplit-4x30：
  12,132 → 1（obj 15.0，x0.091）。
- 零回退：13 个既有算例目标/节点/迭代**全部逐相等**、墙钟噪声带内（门控
  尖效）。摘要矩阵五个『显著提速』格照例驳回为跨窗伪影。

**未做（按提示条件）**：可行性泵本体（判定死因后绕开）；milp-mixed 根
incumbent 收益（要求先有更大规格算例证明墙钟收益，未生成新算例故不引入）。

### R6（e1d19c7）全局 DFS 工作栈 + 共享伪代价负结果 —— 轴3

> ⚠️ 记录勘误：**commit 标题为「共享割池（全局有效割族）+ 确定性并行模式，
> t4 效率 2.1×→3×+」，实际内容是全局 DFS 工作栈（`--work-stack global`）**
> ——共享割池与确定性并行两个方向判定先决条件不成立后未实现（复核低危
> 问题 #2，R7 已在 §4a 对照更正）。标题中的「t4 效率 2.1×→3×+」亦未达成
>（见下）。

**改动**：`parallel_bb.rs`（+54/−15：`PoolState.work_global` +
`complete_node` 全局栈路由 + notify 参数）、`branch_bound.rs`（+7：
`BbOptions.global_work_stack` 默认 false）、`main.rs`（`--work-stack
local|global`）。

**靶点 (1a) 全局伪代价共享——实测显著负收益，完整回退**：实现
`Pseudocosts::merge_from/snapshot` + `Shared.pc_shared` + worker 每 16/64
节点周期同步。mknap-5x300 t4 实测：同步 16 ⇒ nodes 587,209（vs 本地亲和
基线 197k–251k，**膨胀 3×**）；同步 64 ⇒ 549,898 / 44.2s——两档均大幅
恶化。归因：跨子树平均稀释本地高质量统计（伪代价增益 LP-几何局部性强），
且聚合快照使可靠性 η 门对本地未探测变量误判可靠、跳过本应做的探测。连同
`search_ctrl.rs` API 完整回退（git diff 零净改动）。

**靶点 (1b) 全局 DFS 工作栈——正收益落地**：`work_global` 模式位 = true ⇒
子代压全局 steal 栈顶、worker 恒弹栈顶（LIFO）——全局探索序对齐顺序 DFS
（谁弹到谁延续 t1 的下一节点），incumbent 时序/剪枝与 t1 同构；栈顶 push
走 notify_all 快路径。默认 false ⇒ 本地亲和模式逐字节不变。

**实测（mknap-5x300，4 配置同窗紧交错，无时限，全部 optimal 4762）**：
- 复核窗：t4 global 16.1s（nodes 100,093/99,411）vs t4 local HEAD 33.3s
  （nodes 193k–243k）与 NEW 29.3s（nodes 224k–251k）——**墙钟 −45~−52%、
  节点 −48~−60%**；负载尖峰 rep 中 global 节点仍最低（186k < 任何 local
  run），机制性而非侥幸。
- 实现者度量窗：t1 14.5s（87,800 节点）；t4 local 28.7s（202k）；t4 global
  23.0s（130k：膨胀 −36%、墙钟 −20%）——与复核窗方向一致、幅度更保守。
- 附带：knapsack-2000 t4 nodes 133→116；milp-mixed/marketsplit 持平。
- **诚实结论：t4/t1 = 0.63×（global，实现者窗）至 ≈ 1.0–1.15×（复核窗）
  ——轴3 理想值（t4/t1 ≥ 3×）在本架构下不可达**：t4 聚合吞吐被单锁池竞争
  与并行子树交错吃掉；剩余膨胀机理上无法在不退化为单 worker 的前提下消除。
  按风险提示预案降级目标如实上报。
- 零扰动：15 算例 t1 默认路径 A/B 目标/节点/迭代全部逐相等；t4 默认（local）
  ≈ HEAD（并发方差带内）。摘要矩阵五个『显著回退』格照例全部驳回为跨窗伪影。

**未做（判定先决条件不成立，非遗漏）**：共享割池——先决条件「热启动
add_row」不存在；确定性并行批序——不改变其吞吐不可用结论。两者连同
marketsplit-4x40 t4 症状列入遗留（§4）。

### R7（fcea67c）独立复核低危发现修复 —— 复核遗留

**改动**：`milp/src/cuts.rs`（+61/−16，全部为注释/文档更正与一处防御索引
口径）与本报告（§4a 更正小节 + 遗留 #7 状态更新）。**行为零变化**（唯一
代码行改动为不可达防御分支，见下）。

- **丢 π 项数学方向勘误（首轮复核低危 #1）**：割为 Σπ_j·d_j ≥ 1（π ≥ 0、
  d ≥ 0），从 LHS 丢 π>0 项使割**更强**、可能误切整数可行解——原注释
  「丢弃只放宽割，安全」方向错误。逐处更正：`gomory_cuts` 公开入口文档、
  `gomory_cuts_inner` 四处内联（近零 π 项 / origin 越界 / 无记账）、
  `aggregated_mir_cuts` 模块文档 + 三处内联及批序注释引述，并写明**不可达
  前提**（接线路径 cut_loop → cut_loop_recorded：ledger 自空完整累积、割
  循环从冷 Simplex 启动 ⇒ 每个割松弛列必有记账；仅公开入口被手工用于带割
  Simplex 时可达，库内无此调用方，调用方自负）。MIR 源行/assemble 上下文的
  「丢非负项只放宽」（Le 形式右端抬高）经复核数学确实正确，未动。
- **索引防御口径（首轮复核低危 #3）**：`aggregated_mir_cuts` 合并目标处
  `&problem.constraints[ci]` 改为与 GMI 一致的 `.get(ci)` + 防御 continue
  ——构造保证 origin 恒有效（`RowOrigin::Constraint(ci)` 唯一构造点为同一
  problem 约束的 enumerate 下标）、分支不可达，语义零变化。
- **R5/R6 提交标题对照更正（首轮复核低危 #2）**：见 §4a。

**验证（stash-HEAD 探针对照）**：knapsack-2000（35328.0/n=39/it=16994/
cuts=15）、mknap-3x100（2331.0/2120/23978/57）、milp-mixed-250i40c
（1275.791752/53/1472/20）三算例目标/节点/迭代/割数**逐字节相等**，同窗
墙钟持平；lib 240（debug 85.47s / release 21.37s）、acceptance 13、
parity 1、release cuts 34 全过，零测试改动零重钉。

### R8（cac558e）setcover 根割循环 Bland 爬行白烧修复：Dantzig 停滞后的退化处理 —— 轴2（R2 遗留诊断）

**改动**：`simplex.rs`（+90：Bland 释放机制 + `set_iter_limit_hint` +
模块文档）、`cuts.rs`（+36：割循环收紧预算 + IterationLimit 回滚语义）、
`branch_bound.rs`（+19：根 limit 终态 incumbent 提取）、
`tests/acceptance.rs`（1 处重钉）。三机制全部由插桩实测数据驱动；**无新
CLI 旗标**（全部内部机制）。

**诊断（重打 R2 插桩复现）**：setcover-100x150 第 5 轮（m=167）Phase-2 以
bland=true 打满 10000 上限 ⇒ branch_bound 整树回退（丢全部割 + 冷重建）；
setcover-300x600 / 400x800 根 LP 本身以 bland=true 打满 18000/24000。
LP_LIMIT_OVERRIDE 实验：300x600 在 100k 上限收敛于 66–88k 迭代、400x800 在
200k 收敛于 199619——**根因是 Bland 最小下标入基在宽表（ncols ≫ m）上
单步进展极小，而非 LP 不可解**。

**机制 1：Bland 释放（`simplex.rs`）**——停滞（连续退化 > 50 步）触发的
Bland 在「非退化步（t > tol.feasibility，目标真改善 ⇒ 已脱离触发退化区）+
驻留满 `BLAND_DWELL` = 8 步」后回 Dantzig 定价；每次调用最多
`BLAND_REARM_MAX` = 64 次重入，预算耗尽后 Bland 永久驻留（无循环定理的
终止路径恒在场）；总迭代触发（> 3×(m+n)）置 `bland_permanent`（保底）。
全部计数器驱动、不读墙钟（确定性）。

> ⚠️ **终局独立复核中危发现（未修复，见 §4 遗留 #9）**：
> `bland_release_step`（`simplex.rs:1992-2005`）的释放条件**未含
> `!bland_permanent`**——总迭代触发的「永久 Bland」实际同样会被释放（以
> ≥9 迭代间隔、单迭代 Dantzig 窗口、最多 64 次），与 `simplex.rs:1976`、
> 1990-1991 及 `BLAND_DWELL` 常量文档（`simplex.rs:92-106`）三处「不释放」
> 声明矛盾。终止性仍成立（释放次数有限、预算耗尽后恒驻留）、确定性不受
> 影响；修复只需在释放条件加 `!self.bland_permanent`。本节机制描述按代码
> 实际行为口径书写（重入预算兜底终止）。

**机制 2：割循环收紧预算 + 回滚（`cuts.rs:1722-1762`）**——每轮 resolve 经
`set_iter_limit_hint(Some(max(3_000, 12×(m+n))))` 收紧（暖启动对偶修复快
路径的收敛轮实测 ~1.5k 迭代），resolve 后立即复位 hint（不泄漏到后续
调用）；IterationLimit ⇒ **回滚语义**：弃本轮增量割（ledger/recorded/
applied 三账截断回轮前）、重建冷 Simplex + 重放已提交割 + resolve（确定性
冷路径 ⇒ 终态必为上轮末的 Optimal；防御性保留非 Optimal 映射不冒充），
保留已提交割的根界继续树搜索——不再以 limit_reached 无解终态传播。
Infeasible 信号不受影响（Phase-1 判定在预算内第一时间如实返回；该句为
经验性断言的注记见遗留 #9）。

**机制 3：根 limit 终态 incumbent 提取（`branch_bound.rs`）**——根 solve
以 IterationLimit 收场时，Phase-2 中断点的基是原始可行的（比值测试每步
保持 rhs 在界内、Phase-1 已完成）——`round_and_check`（原始问题逐约束 +
整数性复核，不依赖表格语义）通过则作为 incumbent 随 LimitReached 携带
返回；Infeasible/Unbounded 不提取；提取失败如实保持无 incumbent。dive
提取实测否决（白烧 18k 迭代零产出，已移除）。

**实测（vs fcea67c 探针，两窗两序）**：
- **setcover-100x150：nodes 45→13、lp_iterations 18924→5246（−72.3%）、
  cuts 0→23、objective 276 全部 rep 精确一致**；墙钟两窗均降（负载尖峰窗
  2.868→1.501 x0.523；常规窗 0.655→0.565 x0.862；实现者窗口 −49%）——
  方向稳定，幅度随窗。debug 重负载验收同窗 148→111.8s（实现者窗 69.8s）。
- 零回退（行为面逐位对账，全部本机）：knapsack-500（536/10485/15）、
  knapsack-2000（39/16994/15）、knapsack-5000（52/53042/26）、mknap-3x100
  （2120/23978/57）、milp-mixed-250i40c（53/1472/20）、mknap-5x300 无时限
  （87800/1,747,982/70，optimal 4762.0，墙钟 12.3/12.6s）、marketsplit-4x30
  （1/310/22）、knapsack-200（1/1572/19）、lp-sparse（1/7774/0）、
  setcover-400x800（1/24000/0）——同工作墙钟 x0.77–1.05 噪声带。30s 截断
  态 mknap 节点差异为墙钟驱动降级时点的已知非确定性（无时限口径逐位同）。
- **300x600 / 400x800 如实未解开**：根 LP 需 66k/200k 迭代，上限放宽会使
  400x800 根独燃 ~100s+（爆 bench 90s 预算——时限只在节点边界检查，
  `branch_bound.rs:711`）；incumbent 提取在该两类上失败（分数点取整违反
  覆盖约束 ⇒ round_and_check None，诚实无 incumbent）；上限公式放宽已评估
  并否决（收益/成本倒挂）。
- 夸大驳回：实现者自报的回退（knapsack 族 +2.8~11.7%、mknap-3x100 +16%）
  与额外收益（milp-mixed −38%）均为行为逐位相同下的窗口噪声（实测
  x1.05/x1.04/x0.77/x0.98/x1.05）；**commit 标题中的 setcover-400x800
  x0.61 同为跨窗伪影**（24000 lpi 两侧恒同，行为零变化）。

**重钉**：acceptance setcover-100x150 cuts==0 → cuts==23（旧钉值对应
『割循环 limit ⇒ branch_bound 整树回退丢弃全部根割』旧行为；重钉理由入
注释，目标 276 恒为已知最优，该测试自身断言 t1/threads/nocuts 全口径
276 不变）。

### R9（无提交，回退）并行效率修复：锁竞争拆分与 best-bound 协调 —— 轴3（R6 遗留）

**目标**：t4/t1 从 ~1.0× 提到 ≥1.5×（R6 遗留靶点：锁竞争拆分与
best-bound 协调）。

**裁决：回退——空集裁决，无可保留物**。本轮代码面为空：工作树与 HEAD
（cac558e）代码零差异（裁决时 git diff 仅剩 bench/gate_runs.json 与
raw_runs_gate.jsonl 两个计时记录文件，均未提交；当前分支 HEAD 仍为
cac558e、工作树干净）；净树门禁全过（release lib 240 / acceptance 13
裁决方本机复跑）。而声称的全部产出「证据链（std Instant 插桩 worker 时间
分解）」在汇报中**不存在**：没有任何时间分解数字、每 worker 百分比、
t2/t4 对比、或从证据得出的结论——插桩「用后完整移除」连同发现一起消失。
一轮以证据为唯一产出却未汇报任何证据，则该产出对裁决而言不存在：无机制、
无测量结果、无提速，保留条款无从谈起。

**矩阵核查**：全部算例行为逐位不变（nodes 13/39/1/1/87800 与基线全同——
代码零差异的必然结论）；setcover-400x800 x1.13 为同一确定性工作
（24000 lpi 恒同）上的窗口噪声，照例驳回。

**对后续轮次的要求**：以「证据先行」回应回退提示的方向正确，但证据轮的
交付物就是证据本身——插桩数字必须进入汇报（LP resolve / 池临界区 /
condvar 等待 / 窃取扫描 / incumbent RwLock 各占比，t2 与 t4 各 worker
分解），并由数字导出结论（内存带宽受限 vs 锁受限 vs 等待受限）与下一
机制提案——否则插桩做完即删等于没做。

---

## 3. 验证

### 3.1 测试矩阵（各轮实现者本机自跑；debug + release 双口径）

| 轮 | lib 单测 | bin | acceptance | parity | 重负载 `--ignored` | 新增测试 |
|----|---------|-----|------------|--------|--------------------|----------|
| R1 | 230 | 20 | 13 | 1 | debug/release 全过（setcover-100x150 重负载 58.8s） | simplex 5 项：规模门钉、fill-in/墓碑白盒、跨 solve 周期不变量、`tableau_row_terms` 对照稠密、散射 vs 逐列 rc |
| R2 | 232 | 20 | 13 | 1 | debug 71.27s / release 3.35s | cuts 2 项：50 种子穷举 fuzz 零误切（3072 整数点 × 连续两端点）、产出非空 + 逐字节确定性 |
| R3 | 234 | 20 | 13 | 1 | debug 100.93s / release 3.49s | branch_bound 2 项：类检测全分支（命中 + 五种排除 + 64 行边界）、入口解析全路径；gap 测试显式 mostfrac 固定（非重钉） |
| R4 | 239 | 20 | 13 | 1 | debug 111.05s / release 7.34s | scaling 4 + 1 项：2 的幂/逐位往返、**类型红线**（Integer/Binary 列因子恒 1）、均匀幂等、缩放空间最优值穷举一致、bench 族生效面诊断 |
| R5 | 240 | 20 | 13 | 1 | debug 67.22s / release 16.37s（lib 含 MITM 解锁测试） | heuristics 1 项：4x40/5x50 零→有转变 + 逐行精确验证 + 两次调用逐字节确定性 |
| R6 | 240 | 20 | 13 + 重负载 2 | 1 | release 5.5s / debug 148.4s | 未新增（正交性由 15 算例 t1 A/B + t4 global/local 对比覆盖；`complete_node` 第 5 参在测试内 3 处补 false，语义不变） |
| R7 | 240 | 20 | 13 | 1 | 随 acceptance 全过 | 零测试改动（注释/防御/文档级；行为零变化由 stash-HEAD 探针对照背书） |
| R8 | 240 | 20 | 13 + 重负载 2 | 1 | debug 70.7s / release 16.8s（lib）；重负载验收同窗 148→111.8s | acceptance 1 重钉（setcover cuts 0→23）；release simplex+cuts 72 |
| R9 | —（净树，无代码） | — | 13（裁决方复跑） | 1 | — | 零（空集回退，无提交） |

**acceptance 重钉共 4 处**（均为行为值，objective 断言一字未动）：milp-wolsey-2
nodes 7→1 与进度首行 incumbent 37→40（R2 聚合 MIR 根割界收紧到整数最优 40）、
hard-knapsack-12 nodes 159→196（R3 类自适应默认翻 Pseudocost）、
setcover-100x150 cuts 0→23（R8 割循环回滚语义）——四处均本机实跑 CLI
复核一致，重钉理由与演进史入注释。`warm_cold_parity.rs` 全程一字未改、
全部通过。

**threads=1 确定性**：求解路径无 HashMap/HashSet 迭代（仅 instance.rs 校验
用），新增遍历全为下标升序；MITM 用固定种子 LCG + Fisher-Yates；复核实测
算例双跑（setcover-100x150、mknap-3x100、marketsplit-4x40、setcover-400x800）
stdout **逐字节相同**。

### 3.2 首轮独立复核（R6 后，定点抽查）

逐文件复核了当时 main...HEAD 全部源码改动（simplex/cuts/scaling/heuristics/
branch_bound/parallel_bb/solver/main/lib/profile_stats/acceptance/bench
脚本）并定点实跑：**7 项重点全部核查通过**（稀疏簿记不变量、聚合 MIR 数学、
scaling 往返与类型红线、探测纪律、测试改动无静默削弱、threads=1 确定性、
bench 契约），未发现正确性阻塞问题。发现的 3 个低危问题由 R7 全部修复
（见 R7 节与遗留 #7）。

**未跑（如实声明）**：全量测试套件与全量 bench 矩阵按编排要求留给终局脚本；
t4 逐字节确定性未验证（t1 才是确定性契约，t4 不做此承诺）。

### 3.3 终局独立复核（R8 后，全量 10 提交）

独立复核 `feat/gurobi-large-scale` 相对 main 的全部代码改动（10 提交，
核心 ~2900 行：simplex 稀疏簿记 +1055、cuts 聚合 MIR +589、scaling 新
493、heuristics MITM +324、branch_bound +248、parallel_bb +63、
solver/main/ffi/tests/bench）。**结论：通过，无正确性阻断项**；9 项重点
核查均通过（稀疏簿记不变量与密度阈语义四处一致、聚合割数学逐项符号核算、
scaling 2 的幂往返与类型红线、探测纪律无状态泄漏、测试无静默削弱、t1
确定性、bench 契约、setcover 割循环预算限定与回滚截断一致、并行全局栈
与 MITM 接线）。发现 **1 项中危 + 4 项低危/信息项**（均未修复，随遗留
清单 #9 跟踪）。

定向实跑记录（复核方本机）：`cargo test --lib simplex::`（38 通过，debug
下 nz_debug_verify 回归门活跃）、`cuts::tests::aggregated`（2）、
`scaling::`（5）、acceptance 重钉 2 项通过；release 求解 setcover-100x150
（optimal 276 / 13 节点）、mknap-3x100（2331/2120）、marketsplit-4x40
（19.0 / 854952，与本报告 §5 一致）、knapsack-5000（88526 == 提交的 DP
真值）、lp-sparse-4000x800（5893.705757，稀疏启用态路径）、setcover-400x800
（limit_reached，与文档口径一致）。

**未跑（如实声明）**：全量测试套件（留给终局脚本）；threads=4 并行求解
实测（并行审查为代码级，未做运行级验证）。

### 3.4 bench 计时矩阵（`bench/run_timings.py`）

- **gate 集**（快速回归，约 3 分钟）：6 算例 × 30s × 指定 threads；快算例
  3 次、三大算例各 1 次；产物 `bench/gate_runs.json` + `raw_runs_gate.jsonl`
  （HEAD 内为 R8 后 2026-10-12 版本）。
- **full 集**：18 算例（原 13 + 新 5：lp-sparse-4000x800 / setcover-400x800 /
  knapsack-5000 / mknap-10x200 / assign-40x40）× threads {1,4} × 90s；
  快算例 3 次、慢算例（≥20s 量级 6 个）2 次；中位 JSON + 逐次原始记录。
- 更早轮次的 `raw_runs.jsonl` 保留不动，供 `analyze_baseline.py` 继续消费。
- 已知口径注记：时限只在 B&B 节点边界检查（`branch_bound.rs:711`），
  simplex 内部无时限——lp-sparse 单节点纯 LP 在 30s gate 名义时限下仍约
  34s optimal（非 limit），脚本注释已存档。

改动前后 gate 集行为字段对比（threads=1，确定可比）：knapsack-2000 nodes
594→39、cuts 口径下最优恒 35328；mknap-5x300 auto_degraded（166,123 节点
@17.7s 截断）→ optimal（87,800 节点，30s 预算内证完）；**setcover-100x150
nodes 45→13（R8，lp_iterations 18924→5246、cuts 0→23，最优 276 恒证）**；
knapsack-200 与 setcover-400x800（limit_reached @ 24000 上限）行为不变。

---

## 4. 遗留问题

按优先级/主题列出，均含证据出处：

1. **并行效率与两个未落地方向（轴3，R6；R9 空集回退）**：t4/t1 ≈ 0.63–1.15×，
   3× 目标在本架构（单锁节点池 + 独立 dense 表格 worker）下不可达——单锁池
   竞争 + 并行子树交错是结构性瓶颈。**共享割池**（先决条件：热启动
   `add_row`——当前 worker 割重放为冷重建）与**确定性并行批序**（吞吐
   不可用结论未变）判定后搁置。全局伪代价共享已实证负收益（节点膨胀 3×）
   并完整回退，勿轻易重试。R9 以「证据先行」复投后**空集回退**（插桩做完
   即删、证据未进汇报）——后续复投轴3 必须先交付 worker 时间分解数字
   （LP resolve / 池临界区 / condvar 等待 / 窃取扫描 / incumbent RwLock
   各占比，t2/t4 分解）再谈机制。
2. **marketsplit 族收尾（R5）**：marketsplit-5x50 可行化达成（incumbent
   25.0）但**证最优未闭合**（limit_reached）；marketsplit-4x40 在 t4 并行下
   2M 节点内不闭合且 incumbent 停 20（t1 855k 节点闭合最优 19）——并行
   交错延迟最优解发现是该算例的并行瓶颈（R5 以来即此行为，非 R6 引入），
   global 栈不解决。
3. **setcover 类迭代上限传播失败模式（R8 后收窄）**：setcover-100x150 已由
   R8 解开（45→13 节点、迭代 −72%、割 23 条存活）；**setcover-300x600 /
   400x800 仍锁死**——根 LP 本身需 66k/200k 迭代（Bland 最小下标入基在
   宽表上单步进展极小），默认上限 18000/24000 内不收敛；上限公式放宽已评估
   并否决（400x800 根独燃 ~100s+，爆 bench 90s 预算），incumbent 提取在该
   两类失败（分数点取整违反覆盖约束）。真正的解锁需要定价/退化处理或
   宽表专用的进基规则（另一工作项）。
4. **R2 权衡项**：mknap-5x300 根界 4774.32 未突破冻结值 4773.72（批额内
   聚合割挤占 GMI/覆盖席位）——树 −42% 换根界 +0.6，若要突破需批额调度
   改进。
5. **R3 权衡项**：单行背包类（knapsack-500/5000）在探测口径下白烧
   （+18ms / +0.62s 绝对）；若后续要求零回退可把类检测收紧为 rows ≥ 2
   （本轮按规格未擅改）。
6. **R1 收益边界**：lp-sparse-4000x800 稀疏簿记收益 0~−5% 且不稳；小算例
   簿记禁用路径有 +1~2% 噪声级 lean（禁用态走逐字节旧代码，属分支开销）。
7. **首轮复核低危 3 项**（§3.2）——**R7 已全部修复**（注释/防御/文档级，
   行为零变化）：丢 π 项注释方向已更正并写明不可达前提（cuts.rs 两处模块
   文档 + 四处内联）、R5/R6 提交标题对照更正见 §4a、`aggregated_mir_cuts`
   索引已改 `.get` 防御口径。
8. **bench 口径**：lp-sparse 在 gate 名义 30s 时限下仍 ~34s optimal（时限
   只在节点边界检查）——如需 gate 严格 <30s 应换算例或给纯 LP 路径加内部
   时限检查（本轮未动 simplex 内部语义）。
9. **终局独立复核发现（1 中危 + 4 低危，均未修复，见 §3.3）**：
   - **（中危）R8 永久 Bland 可被释放，代码与自身文档矛盾**：
     `bland_release_step`（`simplex.rs:1992-2005`）的释放条件不含
     `!self.bland_permanent`——总迭代触发（`simplex.rs:1975-1978`）置位的
     永久 Bland 实际以「≥9 迭代间隔、单迭代 Dantzig 窗口、最多 64 次」被
     间歇关闭，而 `simplex.rs:1976`（「保底终止路径，不释放」）、1990-1991
     （「永久 Bland……不释放」）与 `BLAND_DWELL` 常量文档
     （`simplex.rs:92-106`，「总迭代触发的 Bland 永久生效」）三处均声称
     不可释放。终止性仍成立（释放次数有限，预算耗尽后 Bland 恒驻留）、
     确定性不受影响；且该路径在 setcover 类（R8 目标算例）上必然被触发。
     修复只需在释放条件加 `!self.bland_permanent`。
   - **（低危）误导性注释**：`main.rs:400`
     `work_stack_global: work_stack_global.unwrap_or(false), // 默认 on（R5-retry）`
     ——该行默认值实为 false（local 模式，与 `BbOptions::global_work_stack`
     文档「默认 false」一致），注释「默认 on（R5-retry）」属上一行 mitm
     （默认 true）却挂在本行尾，误导审计与 CLI 默认口径核对。
   - **（低危）`do_pivot_nz` 密度回退分支自身未设防**：`dense_elim` 第二
     析取项（`nz_alive(r)*4 >= ncols`，`simplex.rs:1622`）命中时行 r 走稠密
     全行除法（`simplex.rs:1634-1637`）但不更新 `row_nz[r]`——若除法下溢
     产生恰 0.0，存活索引会指向零值，违反不变量「存活项集合 == 该行精确
     非零集」。经逐路径核对**当前不可达**（密度回退 `simplex.rs:1757-1762`、
     `nz_rebuild_row`、`add_row` newidx 三处都保证稀疏标记行恒满足
     alive×4 < ncols），且消费方（散射/pivot_nz 收集/tableau_row_terms）对
     零值存活项均为值级无操作；建议加注释或 debug 断言钉住该不可达前提。
   - **（低危）`add_row` 列插入索引右移的 u32 回绕口径**：`*e += 1`
     （`simplex.rs:838-846`）在列下标 = 2^31−1 且带墓碑位（0xFFFFFFFF）时
     release 下回绕为 0。实际不可达（单行 dense 表即需 16 GiB）；`nz_rebuild_row`
     与 newidx 构建有 ≥2^31 防御，此移位路径没有——仅记录防御口径不一致。
   - **（低危）R8 注释把经验写成保证**：`cuts.rs:1722-1729`「Infeasible
     信号不受影响（Phase-1 判定在预算内完成……预算不吞真实不可行）」为
     经验性断言而非结构保证——收紧预算 max(3000, 12(m+n)) 同样作用于
     Phase-1 的 run_phase，若某问题的 Phase-1 本身需要更多迭代，会以
     IterationLimit 收场走回滚（保留已提交割、继续树搜索）——真实不可行
     只被延迟到树搜索发现，不会被伪造或吞掉（Infeasible 仍第一时间如实
     返回，且 hint 在回滚 resolve 前已复位走默认上限）。语义安全，但注释
     口径应收窄为经验性。

---

## 4a. 提交标题与实际内容对照更正（R7，首轮复核低危 #2）

**背景**：独立复核（§3.2 原文，低危问题 2）发现两笔提交的标题与实际
落地的机制不符，误导审计线索。历史提交已推送，不做重写/force-push
（避免破坏审计链）；本节按复核原文逐条对照更正，不替历史辩解。

> 复核原文（§3.2）：
> “2. R5/R6 提交标题与内容不符（见各轮勘误块），误导审计线索。”

**R5（4aa081a）**：
- 提交标题：「可行性泵（界固定-失败解锁式）：解开 marketsplit 大规格
  零可行解锁死」
- 实际内容：**子集和 MITM 可行性启发式**（`heuristics::subset_sum_mitm`：
  全变量/窗口半枚举精确求解等式系）。可行性泵（界固定-失败解锁式，
  `feasibility_pump`）在当轮开发中段被实测证伪（LP 驱动机制锚定 LP
  顶点邻域，marketsplit 大规格零产出）后放弃，未进入提交；按回退提示
  的『正确方向』改为脱离 LP 信息源的子集和专用启发式。标题沿用了最初
  任务名，未随实现更换。
- 审计指引：查可行性泵请勿在 HEAD 找——该机制不存在于任何提交；查
  MITM 见 `milp/src/heuristics.rs` 的 `subset_sum_mitm` 及其测试
  `subset_sum_mitm_unlocks_marketsplit`。

**R6（e1d19c7）**：
- 提交标题：「共享割池（全局有效割族）+ 确定性并行模式，t4 效率
  2.1×→3×+」
- 实际内容：**全局 DFS 工作栈**（`--work-stack global`：子代压全局
  steal 栈顶、worker 从栈顶弹，对齐顺序 DFS 探索序）+ 共享伪代价
  聚合的**负结果归因**（实测节点膨胀 3×，已完整回退）。共享割池与
  确定性并行模式在当轮开发中均被实测否决（割池因 add_row→整基重建
  成本 + 根割已收割同族，t4 10.8→24.5s 负收益；确定性批纪元调度
  因吞吐不可用），未进入提交。标题沿用了最初任务名。
- 度量口径更正：t4 效率 **未** 达到 3×（目标诚实降级）；实际为
  t4/t1 ≈ 0.63×（global 栈，mknap-5x300 3-rep 中位数 23.0s vs t1
  14.5s），较 local 亲和模式（28.7s / 节点 202k）节点 −36%、墙钟
  −20%。
- 审计指引：查割池/确定性模式请勿在 HEAD 找——两者均不存在于任何
  提交；查全局栈见 `parallel_bb.rs` 的 `PoolState::work_global` 与
  `BbOptions::global_work_stack`。

---

## 5. 八轮累计口径（threads=1，行为字段确定可比）

| 算例 | 改动前 | 改动后（HEAD = cac558e） | 变化来源 |
|------|--------|--------------------------|----------|
| mknap-5x300（无时限） | 591,752 节点 / 47.7s | **87,800 节点 / 19.2s**（最优 4762 恒证） | R2 聚合 MIR −42% → R3 类自适应 −52% → R4 缩放 −46%（节点累计 −85%） |
| knapsack-2000 | 594 节点 / 0.61s（gate 窗） | **39 节点**（gate_runs.json 2026-10-12 窗 0.265s） | R2 聚合 MIR（−94%）；R3/R4/R8 行为不变 |
| setcover-100x150 | 45 节点 / 0.361s（gate 基线窗） | **13 节点 / 0.240s**（gate_runs.json 2026-10-12；lpi 18924→5246、cuts 0→23，最优 276 恒证） | R8 Bland 释放 + 割循环预算回滚 |
| marketsplit-4x40 | 2M 节点零可行解 | **854,952 节点证最优 19.0** | R5 MITM |
| milp-wolsey-2 | 7 节点 | **1 节点（根闭合）** | R2 聚合 MIR 根割界到整数最优 |
| setcover-400x800 | limit_reached @ 21.5s | limit_reached @ 13.2s（gate_runs.json 2026-10-12；行为同，每迭代成本 −6%） | R1 稀疏簿记；迭代上限失败模式仍在（遗留 #3） |
| lp-sparse-4000x800 | optimal ~33–40s | optimal ~34–44s（0~−5%，不稳） | R1（收益边界如实报告） |
| mknap-5x300 t4（R6） | local 28.7–33.3s / 193k–251k 节点 | **global 16.1–23.0s / 100k–130k 节点**（`--work-stack global`） | R6 全局 DFS 栈 |

墙钟均为同窗交错 A/B 中位数（R6 t4 行为两窗并列给出）；跨窗墙钟不可比
（§开头口径声明）。
