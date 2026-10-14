//! milp-ffi —— [`milp`] 求解器的 panic 安全 C ABI 绑定层。
//!
//! 目标：`milp/` 零改动前提下，向 Python ctypes（micpython 包）暴露 23 个
//! `milp_ffi_*` C 符号（符号总表见 `README.md`；S1 增批量 add_vars_batch，
//! 纯增符号向后兼容，ABI 版本恒 1）。
//!
//! 模块划分：
//! - [`api`]    ：全部 `#[no_mangle] pub extern "C"` 入口 + 统一 [`api::guard`]
//!                （catch_unwind 兜底）+ 内嵌单元测试；
//! - `error`    ：错误码常量与线程局部粘滞 `last_error`（UTF-8，缓冲契约）；
//! - `handle`   ：`ModelBox`（Problem + Params + StoredSolve + poisoned 标记）。
//!
//! 纪律：
//! - 跨 FFI 边界绝不 panic：所有入口经 catch_unwind 兜底，panic 转错误码 7
//!   并将句柄置『中毒』（其后除 destroy 外一律返回 7）；
//! - 所有缓冲归调用方所有，FFI 侧绝不跨边界 malloc/free；
//! - milp 并行引擎 worker panic 已在库内映射 `MilpError::Internal`
//!   （parallel_bb.rs:1144 先例），本层 catch_unwind 是最后兜底。

pub mod api;
mod error;
mod handle;
