# milp-ffi —— milp 求解器的 C ABI 绑定层

在 **`milp/` 零改动** 前提下，向外部（Python `micpython` 包 / 任意 C 兼容语言）暴露
22 个 `milp_ffi_*` C 符号。`milp-ffi` 是工作区根下的独立 crate（根无 Cargo.toml，
path 依赖 `../milp`），仅依赖标准库与 `milp`，允许 `unsafe`（FFI 固有）。

- 源码：`src/lib.rs`（模块声明）、`src/api.rs`（22 个入口 + panic 兜底）、
  `src/error.rs`（错误码 + 粘滞错误消息）、`src/handle.rs`（`ModelBox` 句柄本体）；
- 测试：`src/api.rs` 内嵌单元测试（guard 兜底）+ `tests/abi.rs` 集成测试
  （仅经导出符号，不触碰 `pub(crate)`）。

## 构建与装载

```text
# 构建（工作区根执行；产物被根 .gitignore 的 target/ 规则忽略）
cargo build --manifest-path milp-ffi/Cargo.toml --release
#   Windows : milp-ffi/target/release/milp_ffi.dll
#   Linux   : milp-ffi/target/release/libmilp_ffi.so
#   macOS   : milp-ffi/target/release/libmilp_ffi.dylib

# 测试
cargo test --manifest-path milp-ffi/Cargo.toml
```

Python 侧装载规则（`micpython` 包实现，此处为契约）：

1. 环境变量 `MIC_FFI_LIB` 非空 → 视为库文件路径（须存在且可加载）；
2. 默认：以 micpython 包目录上一级为工作区根，依次探测
   `root/milp-ffi/target/release/milp_ffi.dll` → `libmilp_ffi.so` → `libmilp_ffi.dylib`；
3. 全失败 → `FfiLoadError`（列全部候选路径 + 提示先执行上述 cargo build 或设 `MIC_FFI_LIB`）；
4. 装载后校验 `milp_ffi_abi_version() == 1`。

## 函数总表（22 个符号）

返回值除注明外均为 `int32_t` 错误码（见错误码表）；`0` = 成功。
ctypes 对应列见文末类型对齐表。

| # | C 签名 | 语义要点 |
|---|--------|----------|
| 1 | `uint32_t milp_ffi_abi_version(void)` | 恒 `1` |
| 2 | `int64_t milp_ffi_last_error(char *buf, uint64_t cap)` | 粘滞错误消息拷出；返回所需字节数（不含 NUL）；`buf=NULL`/`cap=0` 只报长度；否则写 `min(cap-1,len)` 字节并补 `'\0'`（截断不报错） |
| 3 | `int64_t milp_ffi_format_number(double v, char *buf, uint64_t cap)` | 内容逐字节等于 `milp::search_ctrl::format_number`（与 CLI objective 行同格式，含 `0<\|v\|<1e-6` 与 `\|v\|≥1e15` 科学记数分支）；同一缓冲契约；panic（理论不可达）返回 `-1` 并置 last_error |
| 4 | `void* milp_ffi_model_create(const char *name, int32_t sense)` | `sense`：1=MINIMIZE, 2=MAXIMIZE（非法→NULL+last_error(9)）；`name` 可 NULL（=空名），非 UTF-8→NULL(8)；成功返回非空句柄 |
| 5 | `void milp_ffi_model_destroy(void *model)` | Box drop；NULL 安全幂等；**中毒句柄仍可销毁**；destroy 后句柄悬垂（FFI 无法检测，由 Python 守卫兜住） |
| 6 | `int32_t milp_ffi_model_num_vars(void *model, uint32_t *out)` | 变量数 |
| 7 | `int32_t milp_ffi_model_num_constraints(void *model, uint32_t *out)` | 约束数 |
| 8 | `int32_t milp_ffi_add_var(void *model, const char *name, int32_t vtype, double lb, double ub, uint32_t *index_out)` | `vtype`：1=CONTINUOUS, 2=INTEGER, 3=BINARY；`Problem::add_var`（默认界）+ `set_var_bounds(lb,ub)`（lb/ub 恒传入，Python 层已填类型默认）；`ub=+∞` 合法、`lb` 非有限→9；重名/保留字/lb>ub→2；name 为必填（NULL→1），非 UTF-8→8；设置界失败时回滚（变量不入模） |
| 9 | `int32_t milp_ffi_set_var_bounds(void *model, uint32_t index, double lb, double ub)` | 界校验同 add_var；`index ≥ 当前变量数`→9；lb>ub / Binary 越出 [0,1]→2 |
| 10 | `int32_t milp_ffi_add_constraint(void *model, const char *name_or_null, const uint32_t *var_idx, const double *coeffs, uint32_t count, int32_t relop, double rhs, uint32_t *index_out)` | `relop`：1=LE, 2=GE, 3=EQ；`rhs` 须有限（→9）；`name_or_null` 可 NULL（未命名，非 UTF-8→8）；`count=0` 允许两数组指针 NULL（空行照常入模）；系数非有限→9；下标越界→2 |
| 11 | `int32_t milp_ffi_set_objective(void *model, const uint32_t *var_idx, const double *coeffs, uint32_t count, int32_t sense, double offset)` | `sense`：0=保持现状, 1=MINIMIZE, 2=MAXIMIZE（**仅在目标写入成功后**写 `problem.sense`——下标越界→2 失败时模型零突变，对齐 add_var 回滚纪律）；`offset` 即 `obj_offset`；`count=0`（空目标）合法；offset/系数非有限→9 |
| 12 | `int32_t milp_ffi_set_param_real(void *model, const char *name, double value)` | real 型参数（见参数表）；非有限/负→9；类型不符→9；名字不存在→10 |
| 13 | `int32_t milp_ffi_set_param_int(void *model, const char *name, int64_t value)` | int 型参数（`Threads` 1..=1024，越界→9） |
| 14 | `int32_t milp_ffi_set_param_str(void *model, const char *name, const char *value)` | str 型参数（`HeuristicMode`：`auto\|on\|off`，大小写敏感，未知值→9） |
| 15 | `int32_t milp_ffi_solve(void *model)` | 选项组装见下文；`Err` 按 `MilpError` 变体映射 2-6 且**不覆盖旧结果**；`Ok`（含 Infeasible/Unbounded/LimitReached）存档并返回 0，状态经 getter 读取 |
| 16 | `int32_t milp_ffi_get_status(void *model, int32_t *out)` | 0=NOT_SOLVED, 1=OPTIMAL, 2=INFEASIBLE, 3=UNBOUNDED, 4=LIMIT_REACHED |
| 17 | `int32_t milp_ffi_get_objective(void *model, double *out)` | 无解（未求解/无可行解）→11 |
| 18 | `int32_t milp_ffi_get_bound(void *model, double *out)` | 无证据→11 |
| 19 | `int32_t milp_ffi_get_gap_percent(void *model, double *out)` | gap 未触发→11 |
| 20 | `int32_t milp_ffi_get_solution(void *model, double *out, uint32_t cap, uint32_t *written_or_null)` | 维数锚定 `n_solved`（见下文）：无解→11；`cap<n_solved`→9 且 `*written_or_null`（可 NULL）写 `n_solved`；否则恰拷贝 `n_solved` 个 f64 并写 written；`n_solved>0` 时 out NULL→1 |
| 21 | `int32_t milp_ffi_get_var_value(void *model, uint32_t index, double *out)` | `index ≥ 当前 problem.vars.len()`→9；否则 `index ≥ n_solved`→11；否则返回取值 |
| 22 | `int32_t milp_ffi_get_stat(void *model, int32_t stat_id, int64_t *out)` | `stat_id`：1=nodes, 2=lp_iterations, 3=cuts, 4=time_ms, 5=presolve_reductions, 6=gap_terminated(0/1), 7=auto_degraded(0/1)；未知 id→9；未求解→11 |

入参 `const char*` 一律 `CStr`→UTF-8 校验（失败→8）；C 契约要求 NUL 结尾。
枚举 i32 未知值→9（消息列合法值）。

## 错误码表（`int32_t` 返回值）

| 码 | 常量 | 含义 | Python 异常 |
|----|------|------|-------------|
| 0 | MILP_FFI_OK | 成功 | — |
| 1 | ERR_NULL | 句柄/必填指针/输出指针为 NULL | ArgumentError（公共 API 不可达） |
| 2 | ERR_INVALID_MODEL | ← `MilpError::InvalidModel`（重名、保留字、lb>ub 等） | ModelError |
| 3 | ERR_NUMERICAL | ← `MilpError::Numerical` | NumericalError |
| 4 | ERR_PARSE | ← `MilpError::Parse`（FFI 层保留映射完整性） | ParseError |
| 5 | ERR_IO | ← `MilpError::Io`（FFI 层保留映射完整性） | IoError |
| 6 | ERR_INTERNAL | ← `MilpError::Internal`（含并行 worker panic 的库内映射） | InternalError |
| 7 | ERR_PANIC | Rust panic 被 catch_unwind 捕获；该句柄自此**中毒** | PanicError |
| 8 | ERR_UTF8 | 入参字符串非 UTF-8 | ArgumentError（公共 API 不可达） |
| 9 | ERR_ARG | 枚举越界/数组计数不符/参数值域或类型非法/lb 非有限等 FFI 层校验 | ArgumentError |
| 10 | ERR_UNKNOWN_PARAM | 参数名不存在 | UnknownParameterError |
| 11 | ERR_NO_SOLUTION | 未求解或本次求解无解时取 objective/values/单变量值/统计/gap/bound | NoSolutionError |

## 消息缓冲契约（归属与释放）

- **所有缓冲归调用方所有**，FFI 侧绝不跨边界 malloc/free；
- 错误消息经线程局部**粘滞**存储传出：非 0 返回的调用覆盖旧消息，
  **成功调用不清除**；
- `milp_ffi_last_error(buf, cap)` 返回所需字节数（不含 NUL）：
  - `buf == NULL` 或 `cap == 0`：只返回长度，不写任何字节；
  - 否则写 `min(cap-1, len)` 字节并补 `'\0'`；截断不报错，可按返回值重试；
- 消息为 UTF-8（中文）；`milp_ffi_format_number` 同一缓冲规则。

## 句柄生命周期

- `model_create` 成功返回 `Box::into_raw(ModelBox)` 的非空 `void*`；失败 NULL + last_error；
- **句柄非线程安全**：禁止多线程并发使用同一句柄（错误消息也因此按线程隔离）；
- 所有带句柄函数：NULL→1；中毒→7（除 destroy）；结构性修改后可重复 solve，
  每次 solve **整体覆盖**上次结果（X/ObjVal 反映最近一次 optimize）；
- panic 兜底：每个入口经 `catch_unwind` 包装，panic 转错误码 7
  （`Rust panic: <消息>`）并将句柄置中毒——其后除 destroy 外一律返回 7。
  milp 并行引擎的 worker panic 已在库内映射 Internal（parallel_bb.rs 先例），
  FFI 层 catch_unwind 是最后兜底；
- `model_destroy`：NULL 直接返回（幂等）；中毒句柄仍可销毁；悬垂指针由
  Python 侧 dispose 守卫兜住。

### 解维数锚定（n_solved）

`get_solution` 的维数以最近一次求解存档的 `values.len() == n_solved` 为准，
**而非当前 `problem.vars.len()`**：

- 无解→11；`cap < n_solved`→9（`written_or_null` 写 n_solved）；否则恰拷贝 n_solved 个 f64；
- `get_var_value`：`index ≥ 当前 vars.len()`→9；否则 `index ≥ n_solved`→11
  （覆盖「solve 后又 add_var 再取新变量值」——不 panic、不中毒）；否则返回。

## 参数表与关断编码

| FFI 名 | 类型 | 值域 | 默认 | 透传目标 |
|--------|------|------|------|----------|
| TimeLimit | real 秒 | ≥0（**0=不限**） | 0 | `BbOptions.time_limit_ms`（秒→毫秒**四舍五入、正值下限 1ms**：`0.5`→500ms，`0<限时<0.0005s` 夹到 1ms，绝不舍回 0；`time_limit_ms=0` 即不限） |
| Threads | int | 1..=1024 | 1 | `SearchControl.threads`（1=纯顺序，≥2 并行引擎） |
| MIPGap | real 百分比 | ≥0（**0=关闭**） | 0 | `SearchControl.mipgap_percent`（None/Some(p>0)） |
| HeuristicMode | str | `auto\|on\|off`（大小写敏感） | **`on`** | off ⇒ `bb.heuristics=false`；on/auto ⇒ true；auto 另置 `adaptive_heuristics=true` |
| LogPeriodSecs | real 秒 | ≥0（**0=关闭**，>0 开启） | 0 | `SearchControl.log_period_secs` + 仅开启时注入 stdout 进度钩子（镜像 milp CLI `--log-secs`） |

- **0 是规范关断编码，任何后设值可覆盖先设值**：先 `TimeLimit=5.0` 再
  `TimeLimit=0.0` ⇒ 真正解除限速（不残留 5 秒限制）；LogPeriodSecs / MIPGap 同理；
- 非有限或负值一律 ERR_ARG（9）；Python 侧负责把 `float('inf')`（Python 语义
  「不限」）翻译为 0.0 传入——FFI 直收 INF 必须拒（→9）；
- 名字存在但类型入口不符（如 `set_param_int("TimeLimit")`）→9；名字不存在→10；
- **HeuristicMode 默认 `on` 而非 CLI 的 `auto`**：对齐 `SearchControl::default`
  库语义（adaptive_heuristics=false），并避免 auto + TimeLimit 意外触发策略级
  自动降级（LIMIT_REACHED）。需要自适应行为请显式 `set_param_str("HeuristicMode","auto")`。

### solve 的选项组装（固定，无关闭入口）

- `SolveOptions`：`Tolerances::default` + `BbOptions::default`（root_cuts=true、
  use_presolve=true）——仅按 HeuristicMode 覆盖 `heuristics`、按 TimeLimit 覆盖
  `time_limit_ms`（`secs_to_limit_ms`：秒→毫秒四舍五入 + 正值下限 1ms + 饱和，
  见参数表；与 CLI 的整秒 `u64×1000` 截断有意不同——FFI 收 f64 秒）；
- `SearchControl`：threads / mipgap_percent / adaptive_heuristics=(mode==Auto) /
  log_period_secs / degrade=default；
- 映射纪律照搬 milp CLI（main.rs 的 solve_options_from / search_control_from）。

## ctypes 类型对齐表

| C | Rust | ctypes |
|----|------|--------|
| `void*` 句柄 | `*mut ModelBox` | `c_void_p` |
| `int32_t` | `i32` | `c_int32` |
| `uint32_t` | `u32` | `c_uint32` |
| `int64_t` | `i64` | `c_int64` |
| `uint64_t` | `u64` | `c_uint64` |
| `double` | `f64` | `c_double` |
| `const char*` 入参 | `*const c_char` | `c_char_p`（`str.encode('utf-8')`） |
| `char*` 出参缓冲 | `*mut c_char` | `ctypes.create_string_buffer(n)` |
| `const uint32_t*` | `*const u32` | `(c_uint32 * n)` |
| `const double*` | `*const f64` | `(c_double * n)` |
| `double*` 出参 | `*mut f64` | `POINTER(c_double)` 或 `(c_double * n)` |

全部 22 个符号在 Python 侧显式声明 `argtypes`/`restype`（禁默认转换）；
f64 特殊值 ±inf 按位透传（`ub=+∞` 合法；`lb=-inf` 由 FFI 拒绝→9）。

## 验证命令（工作区根执行）

```text
cargo build --manifest-path milp-ffi/Cargo.toml --release   # 产出 milp_ffi.dll
cargo test  --manifest-path milp-ffi/Cargo.toml             # 单元 + 集成测试全绿
```

回归红线：`milp/` 与 `planner/` 零改动。
