//! FFI 错误码常量与线程局部粘滞错误消息（缓冲契约见 README）。
//!
//! 消息为 UTF-8（中文）。粘滞语义：非 0 返回的调用覆盖旧消息，成功调用
//! 不清除。所有缓冲归调用方所有——本模块只做只读拷出，绝不跨界分配。

use std::cell::RefCell;
use std::os::raw::c_char;

/// 成功。
pub(crate) const MILP_FFI_OK: i32 = 0;
/// 句柄 / 必填指针 / 输出指针为 NULL。
pub(crate) const ERR_NULL: i32 = 1;
/// 模型非法（← `MilpError::InvalidModel`：重名、保留字、lb>ub 等）。
pub(crate) const ERR_INVALID_MODEL: i32 = 2;
/// 数值异常（← `MilpError::Numerical`）。
pub(crate) const ERR_NUMERICAL: i32 = 3;
/// 解析错误（← `MilpError::Parse`；FFI 层通常无此路径，保留映射完整性）。
pub(crate) const ERR_PARSE: i32 = 4;
/// IO 错误（← `MilpError::Io`；FFI 层通常无此路径，保留映射完整性）。
pub(crate) const ERR_IO: i32 = 5;
/// 内部错误（← `MilpError::Internal`；并行 worker panic 库内已映射至此）。
pub(crate) const ERR_INTERNAL: i32 = 6;
/// Rust panic 被 catch_unwind 捕获；该句柄自此『中毒』。
pub(crate) const ERR_PANIC: i32 = 7;
/// 入参字符串非 UTF-8。
pub(crate) const ERR_UTF8: i32 = 8;
/// FFI 层参数校验失败（枚举越界 / 值域或类型非法 / 数组计数不符等）。
pub(crate) const ERR_ARG: i32 = 9;
/// 参数名不存在。
pub(crate) const ERR_UNKNOWN_PARAM: i32 = 10;
/// 未求解或本次求解无解时取 objective/values/单变量值/统计/gap/bound。
pub(crate) const ERR_NO_SOLUTION: i32 = 11;

// 线程局部粘滞错误消息（句柄非线程安全 ⇒ 消息按线程隔离即够）。
thread_local! {
    static LAST_ERROR: RefCell<Option<String>> = const { RefCell::new(None) };
}

/// 记录错误消息并原样返回错误码（便于 `Err(set_error(..))` 一行式使用）。
pub(crate) fn set_error(code: i32, message: &str) -> i32 {
    LAST_ERROR.with(|cell| {
        *cell.borrow_mut() = Some(message.to_string());
    });
    code
}

/// 将粘滞错误消息按缓冲契约拷出（`milp_ffi_last_error` 的实现体）。
///
/// 无消息时按空串处理（所需长度 0）。返回所需字节数（不含 NUL）；
/// `buf` 为 NULL 或 `cap` 为 0 时只返回长度不写任何字节；否则写
/// `min(cap-1, len)` 字节并补 `'\0'`（截断不报错，可按返回值重试）。
pub(crate) fn copy_last_error(buf: *mut c_char, cap: u64) -> i64 {
    LAST_ERROR.with(|cell| {
        let message = cell.borrow().clone().unwrap_or_default();
        copy_out(&message, buf, cap)
    })
}

/// 缓冲契约拷出（`milp_ffi_last_error` / `milp_ffi_format_number` 共用）。
pub(crate) fn copy_out(message: &str, buf: *mut c_char, cap: u64) -> i64 {
    let bytes = message.as_bytes();
    let needed = bytes.len() as i64;
    if buf.is_null() || cap == 0 {
        return needed;
    }
    let cap = cap as usize;
    let n = bytes.len().min(cap - 1);
    unsafe {
        std::ptr::copy_nonoverlapping(bytes.as_ptr(), buf as *mut u8, n);
        *buf.add(n) = 0; // 终止 NUL（截断时也保证 NUL 结尾）
    }
    needed
}
