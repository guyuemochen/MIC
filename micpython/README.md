# micpython — milp 求解器的 gurobipy 风格 Python 绑定

零第三方依赖（ctypes + 标准库），Python 3.8+。目标求解器为本地 FFI 库
**milp-ffi**（`milp-ffi/src/api.rs` 的 23 个 C ABI 符号），底层是 `milp`
crate 的 `solve_controlled` 精确分支定界引擎。`import micpython` 无副作用
（惰性装载 DLL）。

## 安装与前置构建

前置：Rust 工具链（构建 milp-ffi）；在工作区根执行：

```bash
cargo build --manifest-path milp-ffi/Cargo.toml --release
```

产物为 `milp-ffi/target/release/milp_ffi.dll`（Linux 为 `libmilp_ffi.so`，
macOS 为 `libmilp_ffi.dylib`；已被根 `.gitignore` 的 `target/` 规则忽略）。

**仓库根即包根**（推荐用法）：`micpython` 是路径导入包，不需要安装——
在工作区根运行 `python`（仓库根在 `sys.path` 上）即可 `import micpython`；
独立脚本请照 `examples/knapsack.py` 的方式自行落路径：

```python
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # 指向仓库根
```

可选：`pip install -e micpython`（`micpython/pyproject.toml` 仅文档级元数据
声明——name、requires-python、packages，零依赖、零构建钩子，**不纳入门禁**，
安装路径映射未经验证，首选上述路径导入方式）。

## DLL 装载规则

装载发生在首次构造 `Model()` 或显式 `micpython.load()` 时（模块级惰性单例）：

1. 环境变量 `MIC_FFI_LIB` 非空 → 视为**库文件路径**，须存在且可加载；
2. 否则按默认候选探测：以 micpython 包目录上一级为工作区根，依次尝试
   `root/milp-ffi/target/release/` 下的 `milp_ffi.dll` → `libmilp_ffi.so`
   → `libmilp_ffi.dylib`（依次 `exists()`）；
3. 全部失败 → 抛 `FfiLoadError`：消息列出全部候选路径与各自失败原因，
   并附提示『先在工作区根执行 `cargo build --manifest-path milp-ffi/Cargo.toml
   --release`，或设 `MIC_FFI_LIB` 指向已构建库』；
4. 装载成功后校验 `milp_ffi_abi_version() == 1`，再逐一校验并绑定全部
   23 个符号（S1 起含批量 `milp_ffi_add_vars_batch`，纯增符号向后兼容；
   缺失或版本不符同样抛 `FfiLoadError`）。

`micpython._ffi._reset_cache()` 仅供测试强制重载。

## 快速上手（与 gurobipy 对照）

0-1 背包示例（改编自 `examples/knapsack.py`，等价
`milp/instances/knapsack-10.milp`，最优 17）。**自包含、零本地辅助函数**：
按上节构建 FFI 库后，把下面的代码保存为工作区根下的 `knapsack_demo.py`，
在工作区根执行 `python knapsack_demo.py` 即可运行：

```python
from micpython import BINARY, MAXIMIZE, Model, Status, quicksum

m = Model("knapsack")

x = m.addVars(5, vtype=BINARY, name="x")
values = [7, 5, 9, 6, 3]
weights = [4, 3, 5, 6, 2]
capacity = 10

m.setObjective(quicksum(values[i] * x[i] for i in range(5)), MAXIMIZE)
m.addConstr(quicksum(weights[i] * x[i] for i in range(5)) <= capacity)

status = m.optimize()
if status != Status.OPTIMAL:
    print("求解未达最优：status = %d（期望 OPTIMAL = %d）"
          % (status, Status.OPTIMAL))
    raise SystemExit(1)

for v in x:
    print(f"{v.VarName} = {v.X:.0f}")
print(f"objective = {m.ObjVal:.6f}")   # 17.000000（定点 6 位小数）

m.dispose()
```

与 gurobipy 的写法逐行对应（`Model`/`addVars`/`setObjective`/`addConstr`/
`optimize`/`v.X`/`m.ObjVal`）。末行输出与 milp CLI 的 objective 行同格式
（定点 6 位小数）；完整示例 `examples/knapsack.py` 另含 `format_number`
（优先经 ABI `milp_ffi_format_number`，与 CLI 逐字节一致，库不可用时退回
脚本自带的纯 Python 等价实现）与脚本自理的导入路径——本节从简省去。

## API 参考

### 常量

| 常量 | 值 | 说明 |
| --- | --- | --- |
| `BINARY` / `INTEGER` / `CONTINUOUS` | `'B'` / `'I'` / `'C'` | 变量类型（`vtype`） |
| `MAXIMIZE` / `MINIMIZE` | `-1` / `1` | 目标方向（gurobipy 同值） |
| `INF` | `float('inf')` | 无穷上界 / TimeLimit『不限』 |
| `GRB.*` | — | gurobipy 兼容别名命名空间（`GRB.BINARY` 等） |

### Status 码（gurobipy 数值风格）

| 属性 | 值 | 语义 |
| --- | --- | --- |
| `Status.NOT_SOLVED` | 0 | 尚未 optimize |
| `Status.OPTIMAL` | 2 | 已证明最优 |
| `Status.INFEASIBLE` | 3 | 不可行 |
| `Status.UNBOUNDED` | 5 | 无界 |
| `Status.LIMIT_REACHED`（别名 `TIME_LIMIT`） | 9 | 触发上限（时间/gap/自动降级） |

### Model

| 成员 | 说明 |
| --- | --- |
| `Model(name='Model', sense=MINIMIZE)` | 创建模型（触发 DLL 装载） |
| `addVar(lb=0.0, ub=None, obj=0.0, vtype=CONTINUOUS, name='')` | 加变量；`ub=None` 按类型默认（C/I→+INF，B→1.0）；`name=''` 自动生成唯一名（避开保留字） |
| `addVars(count, ..., name='x')` | 批量加变量（S1 起单次 FFI 跨界），第 i 个名为 `f"{name}{i}"`，返回 `list[Var]`；任一元素非法整批拒绝、零突变 |
| `addConstr(temp, name='')` | 落模一条 `<=/>=/==` 比较产生的 TempConstr |
| `setObjective(expr, sense=None)` | 设目标；`expr` 为 Var/LinExpr/数值，常数入目标偏移；`sense=None` 保持现方向 |
| `setParam(name, value)` / `getParam(name)` | 参数读写（见下节） |
| `optimize()` | 求解，返回 Status；可重复调用（结构性修改后重解） |
| `dispose()` | 显式释放句柄（幂等）；支持 `with Model(...) as m:` |
| `getVarByName(name)` | 命中返回 Var，未命中返回 `None` |
| 属性 | `ModelName` `NumVars` `NumConstrs` `Status` `ObjVal` `SolCount` `Runtime` `NodeCount` `Bound` `MIPGapPercent` `GapTerminated` `AutoDegraded` |

**setObjective 合并约定**（与 gurobipy 的『完全替换』不同）：最终目标项 =
表达式项 ∪ {未被表达式覆盖的 `addVar(obj=...)` 系数项}；`addVar` 记录的
obj 系数在每次 `setObjective` 时自动并入未被显式覆盖的列。

### Var / LinExpr / quicksum

- `Var` 属性：`VarName` `lb` `ub` `VType` `Obj`（建模镜像，只读）、`X`
  （最近一次 optimize 的解值；未求解/无解 → `NoSolutionError`）；
- `Var`/`LinExpr` 支持 `+`、`-`、`*`（双侧数乘）、`-x`，以及 `==/<=/>=`
  （对常量或另一表达式，双表达式移项合并常数）→ `TempConstr`；`Var*Var`
  等二次项不在 milp 能力边界内（`TypeError`）；
- `LinExpr.getConstant()` / `getValue()`（按最近一次解求值）/ `len()`
  （变量项数）；同下标系数归并、归并为零的项丢弃（对齐 milp `LinearExpr`
  不变量）；
- `quicksum(terms)` 顺序累加 iterable；空序列 → `LinExpr(0)`；实现为共享
  项字典的就地累加（O(n)；S1 修复，旧 O(n²) 算子链在 5000 项时约 1.3s）；
- 混用不同 Model 的变量 → `ArgumentError`（表达式/约束/目标构造时即拒绝）。

### 异常层次

`MicError(Exception)` 为基类，携带 `.message` / `.code`（`.code` 为 FFI
错误码；纯 Python 侧守卫错误 code=0）：

| 异常 | FFI 码 | 语义 |
| --- | --- | --- |
| `FfiLoadError` | — | DLL 装载失败（路径缺失/不可加载/符号缺失/ABI 不符） |
| `ArgumentError` | 1/8/9 | 空指针（公共 API 不可达）、UTF-8 编码失败（公共 API 不可达）、参数值域/类型非法；含内嵌 NUL 的名字/字符串值（Python 侧守卫，code=0）。**亦为 `ValueError` 子类**，可被 `except ValueError` 捕获 |
| `ModelError` | 2 | 模型非法（变量重名、撞保留字、lb>ub、Binary 界越界等） |
| `NumericalError` | 3 | 数值异常（终检失败等） |
| `ParseError` | 4 | 解析错误（FFI 层保留映射完整性） |
| `IoError` | 5 | IO 错误（FFI 层保留映射完整性） |
| `InternalError` | 6 | 内部不变量被破坏 |
| `PanicError` | 7 | Rust panic 被 catch_unwind 捕获；句柄中毒（仅 destroy 可用） |
| `UnknownParameterError` | 10 | 参数名不存在 |
| `NoSolutionError` | 11 | 未求解或本次求解无解时取解/统计/界 |

## 参数

| 参数 | 类型 | 值域 | 默认 | 透传 |
| --- | --- | --- | --- | --- |
| `TimeLimit` | float（秒） | ≥0；`INF` 或 0 = 不限 | `INF` | `milp_ffi_set_param_real` → `BbOptions.time_limit_ms`（秒→毫秒四舍五入、正值下限 1ms：`0.5`→500ms；FFI 0=不限） |
| `Threads` | int | 1..=1024 | 1 | `milp_ffi_set_param_int` → `SearchControl.threads`（1=纯顺序，≥2 并行引擎） |
| `MIPGap` | float（**小数**） | ≥0；0 = 关闭 | 0.0 | 小数 ×100 → 百分比经 FFI → `SearchControl.mipgap_percent` |
| `HeuristicMode` | str | `auto`/`on`/`off` | `'on'` | `milp_ffi_set_param_str`（off ⇒ heuristics=false；on/auto ⇒ true；auto 另置 adaptive_heuristics=true） |
| `LogPeriodSecs` | float（秒） | ≥0；0 = 关 | 0.0 | `>0` 时注入 stdout 进度钩子（镜像 milp CLI `--log-secs`） |

- **MIPGap 单位换算**：`setParam('MIPGap', 0.20)` 即 FFI 的 20.0 百分比
  （CLI `--mipgap 20`）。达阈值提前终止时 `optimize()` 返回
  `LIMIT_REACHED`、`GapTerminated` 为 True、`ObjVal`/`X` 携带 incumbent、
  `MIPGapPercent` 为触发时刻的百分比——终态诚实标注，绝不冒充 OPTIMAL。
- **关断编码（恒调 FFI，无跳过分支）**：
  - `TimeLimit` 的 `INF`（Python 语义『不限』）翻译为 **0.0** 传入（FFI
    拒绝非有限值，FFI 0 = 不限）。先 `setParam('TimeLimit', 5.0)` 再
    `setParam('TimeLimit', INF)` ⇒ FFI 收到 0.0 ⇒ `time_limit_secs=None`
    真正解除限速，无『残留 5 秒照旧限速』的失真；
  - `LogPeriodSecs` 关断直接传 `0.0`（FFI 0 = 关，钩子真正移除）；
  - `MIPGap` 关断传 `0.0`（FFI 0 = `mipgap_percent=None`）。
- **HeuristicMode 默认 `'on'`**（与 milp CLI 的默认 `auto` 不同）：对齐
  `SearchControl::default` 的库语义，避免 `auto` + `TimeLimit` 组合在
  双停滞时意外自动降级为 `LIMIT_REACHED`。需要自适应行为请显式
  `setParam('HeuristicMode', 'auto')`。
- **getParam 返回镜像原值且被拒值不更新镜像**：`setParam` 流程 = 路由表
  类型强制 → INF/单位翻译 → **恒调 FFI setter** → rc==0 才写本地镜像
  `_params[name] = 原始 Python 值`。被 FFI 拒绝的值（如 `Threads=0` →
  `ArgumentError`）不污染镜像；`getParam('MIPGap')` 精确回读 setParam 时
  的小数（0.07 回读 0.07，不做 ÷100 逆运算）。初始镜像：
  `{TimeLimit: INF, Threads: 1, MIPGap: 0.0, HeuristicMode: 'on',
  LogPeriodSecs: 0.0}`。`getParam` 只读本地镜像、不触碰 FFI。

## 生命周期

- **dispose 守卫与幂等短路**：Model 的一切方法（addVar/addConstr/
  setParam/optimize/属性读取，含 `Var.X`）仅经私有 `_call` 触碰 FFI——
  `dispose()` 后使用抛 `MicError('Model 已 dispose，禁止再使用')`（Python
  侧守卫，不会对悬垂句柄发起 FFI 调用、不是段错误）。`dispose()` 第一行
  先短路：`_handle is None` 时**静默返回不抛**（二次 dispose 幂等）。
- **`__del__` 幂等**：仅 `_handle` 非 None 时调 destroy 再置 None——与
  `dispose` 互斥，永不多调；即使误调 `destroy(NULL)` 也是安全幂等的。
  支持 `with Model() as m:`（退出即 dispose）。
- **solve 后结构性修改**：可以继续 addVar/addConstr 并再次 `optimize()`
  （每次 solve 整体覆盖上次结果）。注意两次同名 `addVars` 的第二次会抛
  `ModelError`（FFI/milp 拒绝重名变量，文档化行为）。
- **X/ObjVal 反映最近一次 optimize**：solve 后新增的变量在再次 solve 前
  取 `X` → `NoSolutionError`（解维数锚定最近一次求解）；未求解或
  Infeasible/Unbounded/LIMIT_REACHED-无-incumbent 时读取
  `X/ObjVal/Runtime/NodeCount` → `NoSolutionError`。
- 句柄非线程安全：禁止多线程并发使用同一 Model（ctypes 在 CPython 下自动
  释放 GIL，长求解不阻塞解释器）。

## 能力边界

- **lb 必须有限**：本库不支持 −∞ 下界（`addVar(lb=-INF)` →
  `ArgumentError`）；`ub` 允许有限数或 +INF；
- **名字/字符串值不得含内嵌 NUL**（`'\x00'`）：模型名、变量名（含
  `addVars` 的 name 前缀）、约束名、str 参数值（如 `HeuristicMode`）
  抛 `ArgumentError`（亦为 `ValueError`）——FFI 的 C 字符串会在首个
  NUL 处截断而 Python 镜像保留全名，两侧将静默分叉，入口统一拦截；
- **仅线性**：`Var*Var` / `LinExpr*LinExpr` → `TypeError`（无二次表达式）；
- **LIMIT_REACHED 诚实语义**：只有满足最优性终止条件才返回 OPTIMAL；触发
  时间/gap/自动降级上限一律 LIMIT_REACHED（有 incumbent 时解可用，
  `GapTerminated`/`AutoDegraded` 区分原因）；
- 求解选项固定 `SolveOptions::default()`（root_cuts=true、use_presolve=
  true），FFI 不暴露关闭入口；容差、分支规则等旋钮亦不透出。

## 门禁命令（工作区根）

```bash
# 1. 构建 FFI 库（含 Rust 侧测试）
cargo build --manifest-path milp-ffi/Cargo.toml --release

# 2. Python 测试套件（依赖 micpython/__init__.py 与 micpython/tests/__init__.py）
python -m unittest discover -s micpython/tests -t . -v

# 3. 示例脚本（exit 0 且最后一行恰为 objective = 17.000000）
python micpython/examples/knapsack.py
```

回归红线：`milp/` 与 `planner/` 零改动。

## 目录结构

```
micpython/
├── __init__.py        # 公共 API 导出（惰性装载，import 无副作用）
├── _ffi.py            # ctypes FFI 绑定层（23 符号签名、错误映射、DLL 装载）
├── core.py            # gurobipy 风格建模 API（Model/Var/LinExpr/quicksum）
├── pyproject.toml     # 文档级打包元数据（不作门禁、零依赖零构建钩子）
├── README.md          # 本文档
├── examples/
│   └── knapsack.py    # gurobipy 风格示例（门禁：末行 objective = 17.000000）
└── tests/
    ├── __init__.py    # unittest 包发现所需
    ├── test_loader.py      # DLL 装载规则 / 23 符号 / 错误映射完整性
    ├── test_expression.py  # LinExpr / quicksum / TempConstr 语义
    ├── test_model.py       # 建模计数 / 命名 / 界校验 / 目标方向
    ├── test_solve.py       # 求解终态 / 解读取 / dispose 守卫
    └── test_params.py      # 参数三铁律 / 关断编码 / gap 向量 / 非法值
```
