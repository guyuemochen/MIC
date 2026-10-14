//! C ABI 入口：23 个 `#[no_mangle] pub extern "C"` 符号（总表见 README.md）。
//! （S1 增 `milp_ffi_add_vars_batch`：纯增符号，既有 22 个签名与语义零改动，
//! ABI 版本恒 1 不提升——装载方按符号表逐名校验，旧 DLL 配新绑定在装载期
//! 即报『缺少必需符号』，不存在静默半新半旧。）
//!
//! 统一纪律：
//! - 全部入口经 [`guard`]（`catch_unwind(AssertUnwindSafe(..))`）兜底，panic
//!   永不跨 FFI 边界：捕获后置粘滞错误（码 7，`Rust panic: <消息>`）并将该
//!   句柄 `poisoned = true`（其后除 destroy 外一律返回 7）；
//! - 带句柄入口先做 NULL（→1）/ 中毒（→7）检查；`destroy` 例外——中毒句柄
//!   仍可销毁，NULL 入参安全幂等；
//! - 入参 `const char*` 一律 `CStr`→UTF-8 校验（失败→8）；C 契约要求 NUL 结尾；
//! - 枚举 i32 未知值→9（消息列合法值）；u32 计数先转 usize 再建 Vec 防 OOM；
//! - 入参数组仅在 count>0 时解引用（count=0 允许 NULL，空表达式/空目标合法）。
//!
//! CLI 映射纪律照搬 milp/src/main.rs:294-327（HeuristicMode off ⇒
//! bb.heuristics=false、on/auto ⇒ true、auto 另置 adaptive_heuristics=true），
//! 仅 TimeLimit 秒→ms 换算有意偏离 CLI 的整秒截断：FFI 收 f64 秒，按
//! [`secs_to_limit_ms`] 四舍五入到最近毫秒（CLI 的 `u64 秒 ×1000` 对小数
//! 秒静默截断，0<TimeLimit<1 会变 0=不限——B1 已修复并钉死单测）。

use std::any::Any;
use std::io::Write;
use std::os::raw::{c_char, c_void};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;
use std::sync::Arc;

use milp::branch_bound::{BbOptions, SearchControl, SearchHooks};
use milp::error::MilpError;
use milp::model::{RelOp, Sense, SolveStatus, Tolerances, VarType};
use milp::search_ctrl::{format_number, format_progress_line, ProgressInfo};
use milp::solver::{solve_controlled, SolveOptions};

use crate::error::{
    copy_last_error, copy_out, set_error, MILP_FFI_OK, ERR_ARG, ERR_INTERNAL, ERR_INVALID_MODEL,
    ERR_IO, ERR_NO_SOLUTION, ERR_NULL, ERR_NUMERICAL, ERR_PANIC, ERR_PARSE, ERR_UNKNOWN_PARAM,
    ERR_UTF8,
};
use crate::handle::{HMode, ModelBox};

/// 消费面全部来自 milp 公开 API（`milp/` 零改动红线）：
/// Problem 建模（model.rs:186-301）、solve_controlled（solver.rs:96-166）、
/// SearchControl/BbOptions（branch_bound.rs:88-169）、format_number
/// （search_ctrl.rs:87）、进度钩子范式（main.rs:334-341）。

// ===========================================================================
// panic 兜底与公共校验
// ===========================================================================

/// 从 panic 载荷提取消息（`&str`/`String` 两类，其余给占位说明）。
fn panic_message(payload: Box<dyn Any + Send>) -> String {
    if let Some(s) = payload.downcast_ref::<&'static str>() {
        (*s).to_string()
    } else if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else {
        "（无消息的 panic 载荷）".to_string()
    }
}

/// 统一 panic 兜底包装：body 捕获 panic 后转错误码 7 并中毒句柄（若有）。
///
/// `model` 为 `Some(裸指针)` 时（各带句柄入口），panic 后置
/// `(*ptr).poisoned = true`；`None` 用于无句柄入口。milp 并行引擎的
/// worker panic 已在库内映射 `Internal`（parallel_bb.rs:1144 先例），
/// 本函数是最后兜底。
pub(crate) fn guard(model: Option<*mut ModelBox>, body: impl FnOnce() -> i32) -> i32 {
    match catch_unwind(AssertUnwindSafe(body)) {
        Ok(rc) => rc,
        Err(payload) => {
            if let Some(ptr) = model {
                unsafe {
                    (*ptr).poisoned = true;
                }
            }
            set_error(ERR_PANIC, &format!("Rust panic: {}", panic_message(payload)))
        }
    }
}

/// 句柄检查：NULL→1，中毒→7；通过返回 `*mut ModelBox`。
fn check_handle(model: *mut c_void) -> Result<*mut ModelBox, i32> {
    if model.is_null() {
        return Err(set_error(ERR_NULL, "句柄为 NULL"));
    }
    let m = model as *mut ModelBox;
    if unsafe { (*m).poisoned } {
        return Err(set_error(
            ERR_PANIC,
            "句柄已因 Rust panic 中毒，除 destroy 外的所有操作均被拒绝",
        ));
    }
    Ok(m)
}

/// 入参字符串 CStr→UTF-8 校验（NULL→1，非 UTF-8→8）。
/// 返回引用的生命周期由调用方语境约束（C 内存跨界即刻使用）。
fn utf8<'a>(ptr: *const c_char) -> Result<&'a str, i32> {
    if ptr.is_null() {
        return Err(set_error(ERR_NULL, "字符串入参为 NULL"));
    }
    let bytes = unsafe { std::ffi::CStr::from_ptr(ptr).to_bytes() };
    match std::str::from_utf8(bytes) {
        Ok(s) => Ok(s),
        Err(_) => Err(set_error(ERR_UTF8, "字符串入参不是合法 UTF-8")),
    }
}

/// 读 `count` 个元素的入参数组：count=0 允许 NULL（返回空切片），
/// count>0 且指针 NULL→1。
fn read_array<T>(ptr: *const T, count: usize) -> Result<&'static [T], i32> {
    if count == 0 {
        return Ok(&[]);
    }
    if ptr.is_null() {
        return Err(set_error(ERR_NULL, "数组入参指针为 NULL 但 count > 0"));
    }
    Ok(unsafe { slice::from_raw_parts(ptr, count) })
}

/// `MilpError` 变体→错误码（2-6），消息进粘滞存储。
fn map_milp_error(err: MilpError) -> i32 {
    let code = match &err {
        MilpError::Io(_) => ERR_IO,
        MilpError::Parse { .. } => ERR_PARSE,
        MilpError::InvalidModel(_) => ERR_INVALID_MODEL,
        MilpError::Numerical(_) => ERR_NUMERICAL,
        MilpError::Internal(_) => ERR_INTERNAL,
    };
    set_error(code, &err.to_string())
}

/// SolveStatus→FFI 紧凑状态码（0=NOT_SOLVED 由『未求解』另行给出）。
fn status_code(status: SolveStatus) -> i32 {
    match status {
        SolveStatus::Optimal => 1,
        SolveStatus::Infeasible => 2,
        SolveStatus::Unbounded => 3,
        SolveStatus::LimitReached => 4,
    }
}

/// TimeLimit 秒→毫秒换算（B1 契约，solve 与单测共用）：
///
/// - **四舍五入到最近毫秒**（`(secs*1000).round()`）——不再截断：旧实现
///   `(secs as u64).saturating_mul(1000)` 对秒截断，`0 < TimeLimit < 1.0`
///   静默变 0=不限，任意小数向下取整到整秒（0.5s → 0ms 而非 500ms）；
/// - **正值换算结果下限 1ms**：`0` 是 `BbOptions.time_limit_ms` 的『不限』
///   保留编码（branch_bound.rs:97-98），绝不能由舍入产生——`0 <
///   TimeLimit < 0.0005s` 舍入为 0 时夹到 1ms（有限时限就仍是时限）；
/// - 饱和：f64→u64 的 `as` 转换本身饱和（含上溢/+inf → `u64::MAX`）；
///   NaN（setter 校验下不可达）防御性给 1ms 而非 0（不静默变不限）。
///
/// 输入契约：经 `milp_ffi_set_param_real` 写入的值必为有限正秒（非有限/
/// 负→9，`0.0`→`None` 不经本函数）。
pub(crate) fn secs_to_limit_ms(secs: f64) -> u64 {
    let ms = (secs * 1000.0).round();
    if ms.is_nan() || ms < 1.0 {
        1
    } else {
        ms as u64 // f64→u64 饱和转换（+inf/上溢 → u64::MAX）
    }
}

// ===========================================================================
// 基础符号（1-3）
// ===========================================================================

/// ABI 版本；恒 1。
#[no_mangle]
pub extern "C" fn milp_ffi_abi_version() -> u32 {
    1
}

/// 粘滞错误消息拷出（缓冲契约：返回所需字节数不含 NUL；
/// buf=NULL/cap=0 只报长度；否则写 min(cap-1,len) 字节并补 NUL）。
#[no_mangle]
pub extern "C" fn milp_ffi_last_error(buf: *mut c_char, cap: u64) -> i64 {
    match catch_unwind(AssertUnwindSafe(|| copy_last_error(buf, cap))) {
        Ok(n) => n,
        Err(payload) => {
            set_error(ERR_PANIC, &format!("Rust panic: {}", panic_message(payload)));
            -1
        }
    }
}

/// 数值格式化 = `milp::search_ctrl::format_number`（与 CLI objective 行
/// 逐字节一致，含 1e-6/1e15 科学记数分支，search_ctrl.rs:87-97）。
/// 同一缓冲契约；panic（理论不可达）返回 -1 并置 last_error。
#[no_mangle]
pub extern "C" fn milp_ffi_format_number(v: f64, buf: *mut c_char, cap: u64) -> i64 {
    match catch_unwind(AssertUnwindSafe(|| format_number(v))) {
        Ok(s) => copy_out(&s, buf, cap),
        Err(payload) => {
            set_error(ERR_PANIC, &format!("Rust panic: {}", panic_message(payload)));
            -1
        }
    }
}

// ===========================================================================
// 句柄生命周期与建模（4-12）
// ===========================================================================

/// 创建模型句柄。sense：1=MINIMIZE，2=MAXIMIZE（非法→NULL+last_error）；
/// name 可 NULL（= 空名），非 UTF-8→NULL。成功返回 `Box::into_raw` 的非空指针。
#[no_mangle]
pub extern "C" fn milp_ffi_model_create(name: *const c_char, sense: i32) -> *mut c_void {
    let sense_val = match sense {
        1 => Sense::Minimize,
        2 => Sense::Maximize,
        _ => {
            set_error(ERR_ARG, "sense 非法（合法值：1=MINIMIZE, 2=MAXIMIZE）");
            return std::ptr::null_mut();
        }
    };
    let name_str: &str = if name.is_null() {
        ""
    } else {
        match utf8(name) {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(), // utf8 已置 last_error
        }
    };
    // 创建本体只做堆分配，理论不 panic；仍包 catch_unwind 统一纪律。
    match catch_unwind(AssertUnwindSafe(|| {
        Box::into_raw(Box::new(ModelBox::new(name_str, sense_val))) as *mut c_void
    })) {
        Ok(ptr) => ptr,
        Err(payload) => {
            set_error(ERR_PANIC, &format!("Rust panic: {}", panic_message(payload)));
            std::ptr::null_mut()
        }
    }
}

/// 销毁句柄（Box drop）。NULL 安全幂等；中毒句柄仍可销毁；无返回值。
#[no_mangle]
pub extern "C" fn milp_ffi_model_destroy(model: *mut c_void) {
    if model.is_null() {
        return;
    }
    let _ = catch_unwind(AssertUnwindSafe(|| {
        // drop 不可 panic（纯数据），catch_unwind 仅防御性兜底。
        drop(unsafe { Box::from_raw(model as *mut ModelBox) });
    }));
}

/// 变量数。
#[no_mangle]
pub extern "C" fn milp_ffi_model_num_vars(model: *mut c_void, out: *mut u32) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    if out.is_null() {
        return set_error(ERR_NULL, "输出指针 out 为 NULL");
    }
    guard(Some(m), || {
        let n = unsafe { (*m).problem.num_vars() };
        unsafe { *out = n as u32 };
        MILP_FFI_OK
    })
}

/// 约束数。
#[no_mangle]
pub extern "C" fn milp_ffi_model_num_constraints(model: *mut c_void, out: *mut u32) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    if out.is_null() {
        return set_error(ERR_NULL, "输出指针 out 为 NULL");
    }
    guard(Some(m), || {
        let n = unsafe { (*m).problem.num_constraints() };
        unsafe { *out = n as u32 };
        MILP_FFI_OK
    })
}

/// 添加变量：vtype 1=CONTINUOUS,2=INTEGER,3=BINARY；lb/ub 恒传入
/// （Python 层已填类型默认界）。`Problem::add_var`（默认界）+
/// `set_var_bounds(lb,ub)`；后者失败时回滚（弹出刚加的列）保持原子性。
/// lb 非有限→9（本库不支持 -∞ 下界）；ub 须为有限数或 +∞；重名/保留字→2。
#[no_mangle]
pub extern "C" fn milp_ffi_add_var(
    model: *mut c_void,
    name: *const c_char,
    vtype: i32,
    lb: f64,
    ub: f64,
    index_out: *mut u32,
) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    if index_out.is_null() {
        return set_error(ERR_NULL, "输出指针 index_out 为 NULL");
    }
    let name_str = match utf8(name) {
        Ok(s) => s,
        Err(rc) => return rc, // 变量名为必填：NULL→1，非 UTF-8→8
    };
    let var_type = match vtype {
        1 => VarType::Continuous,
        2 => VarType::Integer,
        3 => VarType::Binary,
        _ => {
            return set_error(
                ERR_ARG,
                "vtype 非法（合法值：1=CONTINUOUS, 2=INTEGER, 3=BINARY）",
            )
        }
    };
    if !lb.is_finite() {
        return set_error(ERR_ARG, "变量下界 lb 必须为有限数（本库不支持 -∞ 下界）");
    }
    if !(ub.is_finite() || ub == f64::INFINITY) {
        return set_error(ERR_ARG, "变量上界 ub 必须为有限数或 +∞");
    }
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        match mb.problem.add_var(name_str, var_type) {
            Ok(index) => match mb.problem.set_var_bounds(index, lb, ub) {
                Ok(()) => {
                    unsafe { *index_out = index as u32 };
                    MILP_FFI_OK
                }
                Err(err) => {
                    // 回滚：vars 为 pub 字段，弹出刚压入的列，模型恢复原状。
                    mb.problem.vars.pop();
                    map_milp_error(err)
                }
            },
            Err(err) => map_milp_error(err),
        }
    })
}

/// 批量添加变量（S1）：单次跨界调用建 `count` 个变量；起始下标经
/// `first_index_out` 返回，连续分配（第 i 个新变量下标 = first + i）。
///
/// ABI：
/// - `names`：`count` 个 `const char*` 的数组。数组指针 NULL ⇒ FFI 自动
///   命名 `x{k}`（k 从 0 自增，跳过与现有变量/批内已用名的重名；`x\d+`
///   形如不可能撞 .milp 保留字表——12 个保留字均为纯字母词）——仅为
///   C 消费方便利。micpython `addVars` 语义要求 `f"{prefix}{i}"` 命名与
///   重名拒绝（两次同名调用必须保持报错而非静默改名），恒显式传 names；
/// - `vtypes`/`lbs`/`ubs`：各 `count` 个元素，编码同 `milp_ffi_add_var`；
/// - `count=0`：四个数组指针均可 NULL，合法 no-op，`first_index_out` 写
///   当前变量数（下一个变量将取的下标）。
///
/// 纪律（B2 validate-then-mutate 的批量版）：**任一元素非法 → 整批拒绝、
/// 零突变**。两阶段实现：
/// 1. 预检（不触碰模型）：vtype 枚举（→9）、lb 非有限（→9）、ub 非有限/
///    +∞（→9）、name 元素 NULL（→1）或非 UTF-8（→8）、批内重名（→2）、
///    撞现有变量名（→2）——错误消息带元素下标；
/// 2. 入模：逐个 `Problem::add_var` + `set_var_bounds`（保留字 / lb>ub /
///    Binary 界越出 [0,1] 等模型语义的单一事实源仍在此层，→2）。预检后
///    理论全部成功；防御性兜底：任一步失败则 `vars.truncate(原长度)` 回滚
///    （本函数不触碰 obj/cons/sense，截断即精确还原）后映射错误码。
#[no_mangle]
pub extern "C" fn milp_ffi_add_vars_batch(
    model: *mut c_void,
    count: u32,
    names: *const *const c_char,
    vtypes: *const i32,
    lbs: *const f64,
    ubs: *const f64,
    first_index_out: *mut u32,
) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    if first_index_out.is_null() {
        return set_error(ERR_NULL, "输出指针 first_index_out 为 NULL");
    }
    let n = count as usize; // u32 计数前置 usize 转换
    let vt = match read_array(vtypes, n) {
        Ok(s) => s,
        Err(rc) => return rc,
    };
    let lb_arr = match read_array(lbs, n) {
        Ok(s) => s,
        Err(rc) => return rc,
    };
    let ub_arr = match read_array(ubs, n) {
        Ok(s) => s,
        Err(rc) => return rc,
    };
    // 预检一（纯 FFI 层值域）：枚举与数值界，逐元素、零模型接触。
    let mut var_types = Vec::with_capacity(n);
    for (i, &code) in vt.iter().enumerate() {
        match code {
            1 => var_types.push(VarType::Continuous),
            2 => var_types.push(VarType::Integer),
            3 => var_types.push(VarType::Binary),
            _ => {
                return set_error(
                    ERR_ARG,
                    &format!("第 {i} 个 vtype 非法（合法值：1=CONTINUOUS, 2=INTEGER, 3=BINARY）"),
                )
            }
        }
    }
    for (i, &lb) in lb_arr.iter().enumerate() {
        if !lb.is_finite() {
            return set_error(
                ERR_ARG,
                &format!("第 {i} 个变量下界 lb 必须为有限数（本库不支持 -∞ 下界）"),
            );
        }
    }
    for (i, &ub) in ub_arr.iter().enumerate() {
        if !(ub.is_finite() || ub == f64::INFINITY) {
            return set_error(ERR_ARG, &format!("第 {i} 个变量上界 ub 必须为有限数或 +∞"));
        }
    }
    // 预检二（名字）：names 非 NULL 时逐元素 NULL（→1）/UTF-8（→8）校验。
    let mut name_strs: Vec<String> = Vec::with_capacity(n);
    if !names.is_null() {
        let ptrs = match read_array(names, n) {
            Ok(s) => s,
            Err(rc) => return rc,
        };
        for (i, &p) in ptrs.iter().enumerate() {
            if p.is_null() {
                return set_error(ERR_NULL, &format!("names[{i}] 元素为 NULL（count>0 时每个名字都是必填）"));
            }
            let bytes = unsafe { std::ffi::CStr::from_ptr(p).to_bytes() };
            match std::str::from_utf8(bytes) {
                Ok(s) => name_strs.push(s.to_string()),
                Err(_) => {
                    return set_error(ERR_UTF8, &format!("names[{i}] 不是合法 UTF-8"));
                }
            }
        }
    }
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        if n == 0 {
            // 空批 no-op：写『下一个下标』，模型零接触。
            unsafe { *first_index_out = mb.problem.vars.len() as u32 };
            return MILP_FFI_OK;
        }
        // 预检三（重名）：现有名集合 + 批内增量插入，批内重复/撞现有均 →2。
        let mut used: std::collections::HashSet<String> = mb
            .problem
            .vars
            .iter()
            .map(|v| v.name.clone())
            .collect();
        if names.is_null() {
            // 自动命名 x{k}：跳过与现有/批内已用名重名（x\d+ 不撞保留字表）。
            let mut k = 0usize;
            for _ in 0..n {
                loop {
                    let candidate = format!("x{k}");
                    k += 1;
                    if !used.contains(&candidate) {
                        used.insert(candidate.clone());
                        name_strs.push(candidate);
                        break;
                    }
                }
            }
        } else {
            for (i, name) in name_strs.iter().enumerate() {
                if !used.insert(name.clone()) {
                    return set_error(
                        ERR_INVALID_MODEL,
                        &format!("变量名 `{name}` 重复（批量第 {i} 个，与批内更早元素或现有变量重名）"),
                    );
                }
            }
        }
        // 入模 + 防御性回滚：预检后理论全部成功，仍兜底 Problem 层语义拒绝
        // （保留字/lb>ub/Binary 界）——失败即 vars.truncate(原长度) 精确还原。
        let orig_len = mb.problem.vars.len();
        for i in 0..n {
            match mb.problem.add_var(&name_strs[i], var_types[i]) {
                Ok(index) => {
                    if let Err(err) = mb
                        .problem
                        .set_var_bounds(index, lb_arr[i], ub_arr[i])
                    {
                        mb.problem.vars.truncate(orig_len);
                        return map_milp_error(err);
                    }
                }
                Err(err) => {
                    mb.problem.vars.truncate(orig_len);
                    return map_milp_error(err);
                }
            }
        }
        unsafe { *first_index_out = orig_len as u32 };
        MILP_FFI_OK
    })
}

/// 设置变量上下界：lb 非有限→9、ub 须有限或 +∞、下标越界→9；
/// lb>ub / Binary 越出 [0,1] 等模型语义错误→2（Problem::set_var_bounds）。
#[no_mangle]
pub extern "C" fn milp_ffi_set_var_bounds(
    model: *mut c_void,
    index: u32,
    lb: f64,
    ub: f64,
) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    if !lb.is_finite() {
        return set_error(ERR_ARG, "变量下界 lb 必须为有限数（本库不支持 -∞ 下界）");
    }
    if !(ub.is_finite() || ub == f64::INFINITY) {
        return set_error(ERR_ARG, "变量上界 ub 必须为有限数或 +∞");
    }
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        let var = index as usize;
        if var >= mb.problem.vars.len() {
            return set_error(
                ERR_ARG,
                &format!(
                    "变量下标 {var} 越界（当前变量数 {}）",
                    mb.problem.vars.len()
                ),
            );
        }
        match mb.problem.set_var_bounds(var, lb, ub) {
            Ok(()) => MILP_FFI_OK,
            Err(err) => map_milp_error(err),
        }
    })
}

/// 添加约束：relop 1=LE,2=GE,3=EQ；rhs 须有限；name 可 NULL（未命名）；
/// count=0 允许两个数组指针为 NULL（空行照常入模型，model.rs:592-619 先例）；
/// 系数非有限→9；下标越界→2（Problem::add_constraint）。
#[no_mangle]
pub extern "C" fn milp_ffi_add_constraint(
    model: *mut c_void,
    name_or_null: *const c_char,
    var_idx: *const u32,
    coeffs: *const f64,
    count: u32,
    relop: i32,
    rhs: f64,
    index_out: *mut u32,
) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    if index_out.is_null() {
        return set_error(ERR_NULL, "输出指针 index_out 为 NULL");
    }
    let cname: Option<&str> = if name_or_null.is_null() {
        None
    } else {
        match utf8(name_or_null) {
            Ok(s) => Some(s),
            Err(rc) => return rc,
        }
    };
    let op = match relop {
        1 => RelOp::Le,
        2 => RelOp::Ge,
        3 => RelOp::Eq,
        _ => {
            return set_error(
                ERR_ARG,
                "relop 非法（合法值：1=LE, 2=GE, 3=EQ）",
            )
        }
    };
    if !rhs.is_finite() {
        return set_error(ERR_ARG, "约束右端 rhs 必须为有限数");
    }
    let n = count as usize; // u32 计数前置 usize 转换
    let idx = match read_array(var_idx, n) {
        Ok(s) => s,
        Err(rc) => return rc,
    };
    let coef = match read_array(coeffs, n) {
        Ok(s) => s,
        Err(rc) => return rc,
    };
    for (i, &c) in coef.iter().enumerate() {
        if !c.is_finite() {
            return set_error(ERR_ARG, &format!("第 {i} 项系数非有限数"));
        }
    }
    guard(Some(m), || {
        let terms: Vec<(usize, f64)> = idx
            .iter()
            .zip(coef.iter())
            .map(|(&v, &c)| (v as usize, c))
            .collect();
        let mb = unsafe { &mut *m };
        match mb.problem.add_constraint(cname, terms, op, rhs) {
            Ok(index) => {
                unsafe { *index_out = index as u32 };
                MILP_FFI_OK
            }
            Err(err) => map_milp_error(err),
        }
    })
}

/// 设置目标：sense 0=保持现状,1=MINIMIZE,2=MAXIMIZE（仅当目标项全部
/// 合法、写入成功后才写 `problem.sense`——下标越界返回 Err 时模型零突变，
/// 对齐 add_var 的回滚纪律，B2）；offset 即 obj_offset（model.rs:252）；
/// count=0（空目标）合法；offset/系数非有限→9；下标越界→2
/// （`Problem::set_objective`，写 sense 前先校验）。
#[no_mangle]
pub extern "C" fn milp_ffi_set_objective(
    model: *mut c_void,
    var_idx: *const u32,
    coeffs: *const f64,
    count: u32,
    sense: i32,
    offset: f64,
) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    let new_sense = match sense {
        0 => None,
        1 => Some(Sense::Minimize),
        2 => Some(Sense::Maximize),
        _ => {
            return set_error(
                ERR_ARG,
                "sense 非法（合法值：0=保持, 1=MINIMIZE, 2=MAXIMIZE）",
            )
        }
    };
    if !offset.is_finite() {
        return set_error(ERR_ARG, "目标常数 offset 必须为有限数");
    }
    let n = count as usize;
    let idx = match read_array(var_idx, n) {
        Ok(s) => s,
        Err(rc) => return rc,
    };
    let coef = match read_array(coeffs, n) {
        Ok(s) => s,
        Err(rc) => return rc,
    };
    for (i, &c) in coef.iter().enumerate() {
        if !c.is_finite() {
            return set_error(ERR_ARG, &format!("目标第 {i} 项系数非有限数"));
        }
    }
    guard(Some(m), || {
        let terms: Vec<(usize, f64)> = idx
            .iter()
            .zip(coef.iter())
            .map(|(&v, &c)| (v as usize, c))
            .collect();
        let mb = unsafe { &mut *m };
        // B2：先经 set_objective 校验并写入目标（其内部先全量校验下标再
        // 突变，失败时零副作用，model.rs:252-268），成功后才写 sense——
        // 旧实现在校验前先翻 sense，越界返回 Err 时模型半突变。
        // set_objective 不读 self.sense，故写入顺序无行为差。
        match mb.problem.set_objective(terms, offset) {
            Ok(()) => {
                if let Some(s) = new_sense {
                    mb.problem.sense = s;
                }
                MILP_FFI_OK
            }
            Err(err) => map_milp_error(err),
        }
    })
}

// ===========================================================================
// 参数（13-15）
// ===========================================================================

/// real 型参数：TimeLimit（秒，≥0，0=不限→None；solve 时经
/// [`secs_to_limit_ms`] 四舍五入到毫秒）/ MIPGap（百分比，≥0，
/// 0=关→None）/ LogPeriodSecs（≥0，0=关→None，>0→Some）。
/// 非有限或负→9；名字存在但类型不符→9；名字不存在→10。
#[no_mangle]
pub extern "C" fn milp_ffi_set_param_real(model: *mut c_void, name: *const c_char, value: f64) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    let name_str = match utf8(name) {
        Ok(s) => s,
        Err(rc) => return rc,
    };
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        match name_str {
            "TimeLimit" => {
                if !value.is_finite() || value < 0.0 {
                    return set_error(ERR_ARG, "TimeLimit 须为 ≥0 的有限秒数（0=不限）");
                }
                // 0.0 是规范关断编码（-0.0 == 0.0 一并归 None）。
                mb.params.time_limit_secs = if value == 0.0 { None } else { Some(value) };
                MILP_FFI_OK
            }
            "MIPGap" => {
                if !value.is_finite() || value < 0.0 {
                    return set_error(ERR_ARG, "MIPGap 须为 ≥0 的有限百分比（0=关闭）");
                }
                mb.params.mipgap_percent = if value == 0.0 { None } else { Some(value) };
                MILP_FFI_OK
            }
            "LogPeriodSecs" => {
                if !value.is_finite() || value < 0.0 {
                    return set_error(ERR_ARG, "LogPeriodSecs 须为 ≥0 的有限秒数（0=关闭）");
                }
                mb.params.log_period_secs = if value == 0.0 { None } else { Some(value) };
                MILP_FFI_OK
            }
            "Threads" | "HeuristicMode" => set_error(
                ERR_ARG,
                &format!("参数 `{name_str}` 不是 real 类型入口（请用对应类型的 setter）"),
            ),
            _ => set_error(
                ERR_UNKNOWN_PARAM,
                &format!(
                    "未知参数名 `{name_str}`（合法参数：TimeLimit/Threads/MIPGap/HeuristicMode/LogPeriodSecs）"
                ),
            ),
        }
    })
}

/// int 型参数：Threads（1..=1024）；越界→9；类型不符→9；名字不存在→10。
#[no_mangle]
pub extern "C" fn milp_ffi_set_param_int(model: *mut c_void, name: *const c_char, value: i64) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    let name_str = match utf8(name) {
        Ok(s) => s,
        Err(rc) => return rc,
    };
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        match name_str {
            "Threads" => {
                if !(1..=1024).contains(&value) {
                    return set_error(ERR_ARG, "Threads 须在 1..=1024 内（1=纯顺序路径）");
                }
                mb.params.threads = value as usize;
                MILP_FFI_OK
            }
            "TimeLimit" | "MIPGap" | "LogPeriodSecs" | "HeuristicMode" => set_error(
                ERR_ARG,
                &format!("参数 `{name_str}` 不是 int 类型入口（请用对应类型的 setter）"),
            ),
            _ => set_error(
                ERR_UNKNOWN_PARAM,
                &format!(
                    "未知参数名 `{name_str}`（合法参数：TimeLimit/Threads/MIPGap/HeuristicMode/LogPeriodSecs）"
                ),
            ),
        }
    })
}

/// str 型参数：HeuristicMode（auto|on|off，大小写敏感）；未知值→9；
/// 类型不符→9；名字不存在→10。
#[no_mangle]
pub extern "C" fn milp_ffi_set_param_str(
    model: *mut c_void,
    name: *const c_char,
    value: *const c_char,
) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    let name_str = match utf8(name) {
        Ok(s) => s,
        Err(rc) => return rc,
    };
    let value_str = match utf8(value) {
        Ok(s) => s,
        Err(rc) => return rc,
    };
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        match name_str {
            "HeuristicMode" => match value_str {
                "auto" => mb.params.heuristic_mode = HMode::Auto,
                "on" => mb.params.heuristic_mode = HMode::On,
                "off" => mb.params.heuristic_mode = HMode::Off,
                _ => {
                    return set_error(
                        ERR_ARG,
                        "HeuristicMode 非法（合法值：auto|on|off，大小写敏感）",
                    )
                }
            },
            "TimeLimit" | "Threads" | "MIPGap" | "LogPeriodSecs" => {
                return set_error(
                    ERR_ARG,
                    &format!("参数 `{name_str}` 不是 str 类型入口（请用对应类型的 setter）"),
                )
            }
            _ => {
                return set_error(
                    ERR_UNKNOWN_PARAM,
                    &format!(
                        "未知参数名 `{name_str}`（合法参数：TimeLimit/Threads/MIPGap/HeuristicMode/LogPeriodSecs）"
                    ),
                )
            }
        }
        MILP_FFI_OK
    })
}

// ===========================================================================
// 求解与结果读取（16-23）
// ===========================================================================

/// 求解。选项组装纪律（镜像 main.rs:294-327）：
/// - `SolveOptions`：Tolerances::default + BbOptions::default（FFI 固定
///   root_cuts=true、use_presolve=true，无关闭入口），仅按 HeuristicMode
///   覆盖 `heuristics`、按 TimeLimit 覆盖 `time_limit_ms`（秒→毫秒按
///   [`secs_to_limit_ms`] 四舍五入 + 下限 1ms，0=不限）；
/// - `SearchControl`：threads / mipgap_percent / adaptive_heuristics=
///   (mode==Auto) / log_period_secs / degrade=default；
/// - LogPeriodSecs 开启时注入 hooks（stdout 进度行，镜像 main.rs:334-341，
///   写失败静默）。
///
/// Err 按 MilpError 变体映射 2-6 且**不覆盖**旧 StoredSolve；
/// Ok（含 Infeasible/Unbounded/LimitReached）存 StoredSolve 并返回 0，
/// 状态经 getter 读取。
#[no_mangle]
pub extern "C" fn milp_ffi_solve(model: *mut c_void) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        let params = mb.params; // Copy
        let time_limit_ms = params
            .time_limit_secs
            .map_or(0u64, secs_to_limit_ms);
        let opts = SolveOptions {
            tol: Tolerances::default(),
            bb: BbOptions {
                heuristics: !matches!(params.heuristic_mode, HMode::Off),
                time_limit_ms,
                ..BbOptions::default()
            },
            use_presolve: true,
            use_scaling: true, // R4 数值缩放：与 CLI/库默认同源（目标与解还原原尺度）
        };
        let ctrl = SearchControl {
            threads: params.threads,
            mipgap_percent: params.mipgap_percent,
            adaptive_heuristics: matches!(params.heuristic_mode, HMode::Auto),
            log_period_secs: params.log_period_secs,
            ..SearchControl::default() // degrade: default
        };
        // 进度钩子仅在 LogPeriodSecs 开启时注入（与 CLI --log-secs 同门控）。
        let hooks = params.log_period_secs.map(|_| SearchHooks {
            progress: Arc::new(|info: &ProgressInfo| {
                let mut stdout = std::io::stdout().lock();
                let _ = writeln!(stdout, "{}", format_progress_line(info));
            }),
        });
        match solve_controlled(&mb.problem, &opts, &ctrl, hooks.as_ref()) {
            Ok(rep) => {
                mb.result = Some(crate::handle::StoredSolve {
                    status: rep.report.status,
                    objective: rep.report.objective,
                    bound: rep.bound,
                    gap_percent: rep.gap_percent,
                    values: rep.report.values,
                    stats: rep.report.stats,
                    presolve_reductions: rep.report.presolve_reductions,
                    gap_terminated: rep.gap_terminated,
                    auto_degraded: rep.auto_degraded,
                });
                MILP_FFI_OK
            }
            Err(err) => map_milp_error(err), // 失败不覆盖旧 StoredSolve
        }
    })
}

/// 求解终态：0=NOT_SOLVED,1=OPTIMAL,2=INFEASIBLE,3=UNBOUNDED,4=LIMIT_REACHED。
#[no_mangle]
pub extern "C" fn milp_ffi_get_status(model: *mut c_void, out: *mut i32) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    if out.is_null() {
        return set_error(ERR_NULL, "输出指针 out 为 NULL");
    }
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        let code = match &mb.result {
            None => 0, // NOT_SOLVED
            Some(stored) => status_code(stored.status),
        };
        unsafe { *out = code };
        MILP_FFI_OK
    })
}

/// 目标值；未求解或本次求解无解→11。
#[no_mangle]
pub extern "C" fn milp_ffi_get_objective(model: *mut c_void, out: *mut f64) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    if out.is_null() {
        return set_error(ERR_NULL, "输出指针 out 为 NULL");
    }
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        match mb.result.as_ref().and_then(|stored| stored.objective) {
            Some(v) => {
                unsafe { *out = v };
                MILP_FFI_OK
            }
            None => set_error(ERR_NO_SOLUTION, "无可用目标值（未求解，或本次求解无可行解）"),
        }
    })
}

/// 对偶界；未求解或无证据→11。
#[no_mangle]
pub extern "C" fn milp_ffi_get_bound(model: *mut c_void, out: *mut f64) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    if out.is_null() {
        return set_error(ERR_NULL, "输出指针 out 为 NULL");
    }
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        match mb.result.as_ref().and_then(|stored| stored.bound) {
            Some(v) => {
                unsafe { *out = v };
                MILP_FFI_OK
            }
            None => set_error(ERR_NO_SOLUTION, "无可用对偶界证据（未求解，或本次求解无证据）"),
        }
    })
}

/// gap 百分比；未求解或 gap 未触发→11。
#[no_mangle]
pub extern "C" fn milp_ffi_get_gap_percent(model: *mut c_void, out: *mut f64) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    if out.is_null() {
        return set_error(ERR_NULL, "输出指针 out 为 NULL");
    }
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        match mb.result.as_ref().and_then(|stored| stored.gap_percent) {
            Some(v) => {
                unsafe { *out = v };
                MILP_FFI_OK
            }
            None => set_error(ERR_NO_SOLUTION, "无可用 gap（未求解，或本次求解未因 gap 终止）"),
        }
    })
}

/// 解向量拷出，维数锚定 StoredSolve.values.len()==n_solved（而非当前
/// problem.vars.len()）：未求解/无解→11；cap<n_solved→9 且
/// `*written_or_null`（可 NULL）写 n_solved；否则恰拷贝 n_solved 个 f64
/// 并写 written_or_null。n_solved>0 时 out 为 NULL→1。
#[no_mangle]
pub extern "C" fn milp_ffi_get_solution(
    model: *mut c_void,
    out: *mut f64,
    cap: u32,
    written_or_null: *mut u32,
) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        let values = match mb.result.as_ref().and_then(|stored| stored.values.as_ref()) {
            Some(v) => v,
            None => {
                return set_error(
                    ERR_NO_SOLUTION,
                    "无可用解向量（未求解，或本次求解无可行解）",
                )
            }
        };
        let n_solved = values.len();
        if (cap as usize) < n_solved {
            if !written_or_null.is_null() {
                unsafe { *written_or_null = n_solved as u32 };
            }
            return set_error(
                ERR_ARG,
                &format!("缓冲区容量 {cap} 小于解维数 {n_solved}"),
            );
        }
        if n_solved > 0 && out.is_null() {
            return set_error(ERR_NULL, "输出指针 out 为 NULL");
        }
        if n_solved > 0 {
            unsafe { std::ptr::copy_nonoverlapping(values.as_ptr(), out, n_solved) };
        }
        if !written_or_null.is_null() {
            unsafe { *written_or_null = n_solved as u32 };
        }
        MILP_FFI_OK
    })
}

/// 单变量取值：index ≥ 当前 problem.vars.len()→9；否则 index ≥ n_solved→11
/// （覆盖『solve 后又 add_var 再取新变量值』）；否则写回取值。
#[no_mangle]
pub extern "C" fn milp_ffi_get_var_value(model: *mut c_void, index: u32, out: *mut f64) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    if out.is_null() {
        return set_error(ERR_NULL, "输出指针 out 为 NULL");
    }
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        let var = index as usize;
        if var >= mb.problem.vars.len() {
            return set_error(
                ERR_ARG,
                &format!(
                    "变量下标 {var} 越界（当前变量数 {}）",
                    mb.problem.vars.len()
                ),
            );
        }
        let values = match mb.result.as_ref().and_then(|stored| stored.values.as_ref()) {
            Some(v) => v,
            None => {
                return set_error(ERR_NO_SOLUTION, "无可用解（未求解，或本次求解无可行解）")
            }
        };
        if var >= values.len() {
            return set_error(ERR_NO_SOLUTION, "该变量不在最近一次求解的解向量内（解维数锚定为 n_solved）");
        }
        unsafe { *out = values[var] };
        MILP_FFI_OK
    })
}

/// 求解统计。stat_id：1=nodes,2=lp_iterations,3=cuts,4=time_ms,
/// 5=presolve_reductions,6=gap_terminated(0/1),7=auto_degraded(0/1)。
/// 未知 id→9；未求解→11。
#[no_mangle]
pub extern "C" fn milp_ffi_get_stat(model: *mut c_void, stat_id: i32, out: *mut i64) -> i32 {
    let m = match check_handle(model) {
        Ok(m) => m,
        Err(rc) => return rc,
    };
    if out.is_null() {
        return set_error(ERR_NULL, "输出指针 out 为 NULL");
    }
    if !(1..=7).contains(&stat_id) {
        return set_error(ERR_ARG, "stat_id 非法（合法值：1-7）");
    }
    guard(Some(m), || {
        let mb = unsafe { &mut *m };
        let stored = match mb.result.as_ref() {
            Some(s) => s,
            None => return set_error(ERR_NO_SOLUTION, "尚未求解（无统计可读）"),
        };
        let value: i64 = match stat_id {
            1 => stored.stats.nodes as i64,
            2 => stored.stats.lp_iterations as i64,
            3 => stored.stats.cuts as i64,
            4 => stored.stats.time_ms as i64,
            5 => stored.presolve_reductions as i64,
            6 => stored.gap_terminated as i64,
            7 => stored.auto_degraded as i64,
            _ => unreachable!("stat_id 已在 1..=7 校验"),
        };
        unsafe { *out = value };
        MILP_FFI_OK
    })
}

// ===========================================================================
// 单元测试：guard 兜底（guard 为 pub(crate)，同 crate 可直接访问，
// 不新增导出符号；集成测试见 tests/abi.rs，不触碰 guard）。
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::copy_last_error;

    /// guard 捕获 panic 并返回 7（要求：panic 永不跨边界外泄）。
    #[test]
    fn guard_catches_panic_returns_7() {
        let rc = guard(None, || panic!("boom"));
        assert_eq!(rc, ERR_PANIC);
    }

    /// panic 消息进粘滞 last_error 且带 "Rust panic:" 前缀。
    #[test]
    fn guard_panic_sets_sticky_error_message() {
        let _ = guard(None, || panic!("带消息的崩溃"));
        let mut buf = [0 as c_char; 512];
        let need = copy_last_error(buf.as_mut_ptr(), 512);
        assert!(need > 0, "粘滞消息应非空");
        let bytes = unsafe {
            std::slice::from_raw_parts(buf.as_ptr() as *const u8, need as usize)
        };
        let text = std::str::from_utf8(bytes).unwrap();
        assert!(text.starts_with("Rust panic:"), "实际：{text}");
        assert!(text.contains("带消息的崩溃"), "实际：{text}");
    }

    /// 带句柄的 panic 会中毒句柄；中毒句柄经 check_handle 返回 7。
    #[test]
    fn guard_poisons_model_handle() {
        let boxed = Box::new(ModelBox::new("t", Sense::Minimize));
        let ptr = Box::into_raw(boxed);
        let rc = guard(Some(ptr), || panic!("poison test"));
        assert_eq!(rc, ERR_PANIC);
        assert!(unsafe { (*ptr).poisoned }, "句柄应已中毒");
        assert_eq!(check_handle(ptr as *mut c_void), Err(ERR_PANIC));
        // 中毒句柄仍可销毁（drop 不因 poisoned 失败）。
        unsafe { drop(Box::from_raw(ptr)) };
    }

    /// guard 正常路径透传返回值。
    #[test]
    fn guard_passes_through_ok() {
        assert_eq!(guard(None, || MILP_FFI_OK), MILP_FFI_OK);
    }

    /// B1 回归：TimeLimit 秒→毫秒四舍五入（旧实现按秒截断）。
    ///
    /// 旧实现 `(secs as u64).saturating_mul(1000)` 的失真实证：0.5s→0ms
    /// （=不限！）、2.999s→2000ms、1e-9s→0ms（=不限）。
    #[test]
    fn secs_to_limit_ms_rounds_not_truncates() {
        // 核心钉死向量：0.5 秒必须换算为 500ms（曾静默变 0=不限）。
        assert_eq!(secs_to_limit_ms(0.5), 500, "0.5s → 500ms（B1）");
        assert_eq!(secs_to_limit_ms(1.0), 1000, "1.0s → 1000ms（原即精确）");
        assert_eq!(secs_to_limit_ms(0.999), 999, "亚秒小数不再丢成 0");
        assert_eq!(secs_to_limit_ms(2.999), 2999, "曾截断为 2000ms");
        assert_eq!(secs_to_limit_ms(2.5), 2500);
        // 正值下限 1ms：0 是『不限』保留编码，绝不由舍入产生。
        assert_eq!(secs_to_limit_ms(0.0001), 1, "0.1ms → 夹到 1ms");
        assert_eq!(secs_to_limit_ms(f64::MIN_POSITIVE), 1);
        // 饱和：f64→u64 饱和转换（上溢 → u64::MAX）。
        assert_eq!(secs_to_limit_ms(1e300), u64::MAX);
        // NaN 防御（setter 校验下不可达）：给 1ms 而非 0（不静默变不限）。
        assert_eq!(secs_to_limit_ms(f64::NAN), 1);
    }
}
