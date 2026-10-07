# milp

纯 Rust 标准库实现的混合整数线性规划（MILP）精确求解器：零第三方依赖、
无 `unsafe`、单线程确定性（所有并列一律按最小下标破键），提供人类可读的
`.milp` 算例格式、命令行入口与库 API。

## 定位与能力边界

- **能力**：纯 LP（全连续）、0-1 规划、一般整数与混合类型 MILP 的**精确**
  求解——两阶段有界单纯形 + 深度优先分支定界 + 根节点 Gomory 混合整数割
  （GMI）+ 取整/潜水启发式 + 保守预处理。求解完全确定，同一算例同一结果。
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
milp <算例路径> [--time-limit-secs N] [--verbose]
```

无参时向 stderr 打印用法；`N` 非数字或为负报错（退出码 1）。

### 输出格式

- **Optimal**：stdout 依次输出 `status = optimal`、`nodes = N`（处理的
  分支定界节点数；数值取决于算例与实现）；`--verbose` 时再逐行输出
  `name = value` 的变量取值；**最后一行精确为** `objective = <数值>`。
  数值常规为 `{:.6}` 六位小数（如 `objective = 17.000000`）；当
  `0 < |v| < 1e-6` 或 `|v| ≥ 1e15` 时改用科学记数 `{:e}`（如 `1e-7`）。
- **Infeasible / Unbounded / LimitReached**：仅输出
  `status = infeasible|unbounded|limit_reached`，**没有 objective 行**；
  `LimitReached` 且已有 incumbent 时附一行 `best = <数值>`。
- **错误**：向 stderr 输出 `error: <MilpError>`（中文信息，解析错误附
  1 起始行号），退出码 1。

```text
$ cargo run --manifest-path milp/Cargo.toml --release -- milp/instances/knapsack-10.milp
status = optimal
nodes = 3                # 示意；实际节点数由求解器统计
objective = 17.000000
```

### 退出码

| 退出码 | 含义 |
|--------|------|
| 0 | 最优（stdout 最后一行 `objective = <数值>`） |
| 1 | 错误/用法（加载、解析、模型非法、数值异常等） |
| 2 | 不可行（Infeasible） |
| 3 | 无界（Unbounded） |
| 4 | 上限触发（LimitReached，附已有 incumbent） |

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
- 分支变量规则 MostFractional（max min(f, 1−f)，并列取最小下标）；剪枝条件
  `node_bound ≤ incumbent + 1e-9×(1+|incumbent|)`。
- 终态诚实判定：栈空自然退出且有 incumbent ⇒ Optimal、无 ⇒ Infeasible；
  节点/时间上限触发 ⇒ LimitReached（附 incumbent）。

### 根节点 GMI 割（WP3）

只在根节点生成（对原问题全体整数可行解全局有效；分支子树内调用无效，
禁止）。取基变量为整数类型结构变量的行推导，割为 `Σ π_j·d_j ≥ 1`（π_j ≥ 0）。
**整数列前提**：整数类型结构列仅当所贴界为整数值（|bound−round| ≤
tol.zero）时才可用整数公式，否则必须按连续列公式处理，否则割可切掉整数
最优而错报 Optimal。数值纪律：|c_j| < tol.zero 的项只允许按放宽割的方向
丢弃，严禁加严；violation ≥ min_violation 才收；动态范围超 1e9 整条丢弃。

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

`milp/instances/` 下 7 个算例（5 个数值型 + 2 个状态型）；期望值与状态
判定均已用独立穷举脚本实跑验证。验收测试见 `milp/tests/acceptance.rs`
（数值解析比较，容差 1e-6，不依赖具体小数位数）。

| 算例 | 类型 | 期望结果 | 已知最优解 | 验证方式 |
|------|------|----------|------------|----------|
| `lp-prod-3` | 纯 LP（3 连续变量） | Optimal，objective = 13 | (x1,x2,x3) = (2,0,1) | 3×3 顶点枚举 |
| `knapsack-10` | 0-1 背包（5 变量） | Optimal，objective = 17 | {x2,x3,x5} | 2^5 穷举；LP 松弛 17.667 触发分支 |
| `assign-4x4` | 4×4 指派（16 binary） | Optimal，objective = 10 | W1→J3, W2→J1, W3→J4, W4→J2 | 4!=24 穷举；LP 松弛即整最优且高度退化（Bland 规则测试点） |
| `milp-wolsey-2` | 一般整数 MILP | Optimal，objective = 40 | (x1,x2) = (0,5) | 穷举；LP 松弛 41.25 于 (2.25,3.75) 触发分支 + GMI 机会；整数可行点 25 个（x2=0..5 时 x1 上限 6,5,4,3,1,0），GMI 公式对 25 点全部验证有效 |
| `milp-mixed-5` | binary+整数+连续混合 | Optimal，objective = 29 | (y1,y2,z,c) = (1,1,2,3) | 4 个二元组合 × z 枚举并解析求 c；根 LP 松弛即整数最优 |
| `infeasible-tiny` | 状态型 | Infeasible，退出码 2，无 objective 行 | — | x≥3 ∧ x≤1 空集；正右端 ≥ 行（x≥3 在 lb=0 处违反，走人工变量路径）端到端覆盖 |
| `unbounded-tiny` | 状态型 | Unbounded，退出码 3，无 objective 行 | — | max x，x≥3 且上界 +∞；默认上界 +∞ 的合法性由本例顺带验收 |

## 仓库布局

```text
milp/
├── src/
│   ├── main.rs          # CLI 入口（WP5）
│   ├── error.rs         # MilpError 与 Result 别名（WP1）
│   ├── model.rs         # 建模 + 集中容差（WP1）
│   ├── instance.rs      # .milp 解析/序列化（WP1）
│   ├── simplex.rs       # 两阶段有界单纯形（WP2）
│   ├── presolve.rs      # 保守预处理与还原（WP3）
│   ├── cuts.rs          # 根节点 GMI 割（WP3）
│   ├── heuristics.rs    # 取整/潜水启发式（WP3）
│   ├── branch_bound.rs  # 深度优先分支定界（WP4）
│   └── solver.rs        # 统一编排与安全网终检（WP4）
├── instances/           # 7 个验收算例
└── tests/acceptance.rs  # 验收测试（WP5）
```

单元测试 `#[cfg(test)]` 就近放在各模块内；验收统一在 `tests/acceptance.rs`。
