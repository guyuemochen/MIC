//! 统一错误类型（WP1 契约）。
//!
//! 库代码统一 `Result<T, MilpError>`；用户输入路径禁止 unwrap/expect/panic；
//! 解析错误必须携带 1 起始行号。

use std::fmt;

/// 库统一错误类型（五变体）。
#[derive(Debug)]
pub enum MilpError {
    /// 文件读写错误。
    Io(std::io::Error),
    /// .milp 解析错误；`line` 为 1 起始行号。
    Parse {
        /// 1 起始的出错行号。
        line: usize,
        /// 人类可读的出错说明（中文）。
        message: String,
    },
    /// 模型非法（重复变量名、撞保留字、lb>ub、非有限系数/界/右端等）。
    InvalidModel(String),
    /// 数值异常（主元过滤后无候选、NaN/Inf、终检失败等）。
    Numerical(String),
    /// 内部不变量被破坏（不应发生）。
    Internal(String),
}

impl fmt::Display for MilpError {
    /// 人类可读中文信息。
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MilpError::Io(e) => write!(f, "IO 错误：{e}"),
            MilpError::Parse { line, message } => {
                write!(f, "解析错误（第 {line} 行）：{message}")
            }
            MilpError::InvalidModel(message) => write!(f, "模型非法：{message}"),
            MilpError::Numerical(message) => write!(f, "数值异常：{message}"),
            MilpError::Internal(message) => write!(f, "内部错误（不应发生）：{message}"),
        }
    }
}

impl std::error::Error for MilpError {
    /// 仅 `Io` 变体返回底层 `std::io::Error`，其余返回 `None`。
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            MilpError::Io(e) => Some(e),
            MilpError::Parse { .. }
            | MilpError::InvalidModel(_)
            | MilpError::Numerical(_)
            | MilpError::Internal(_) => None,
        }
    }
}

impl From<std::io::Error> for MilpError {
    fn from(e: std::io::Error) -> Self {
        MilpError::Io(e)
    }
}

/// 库统一 `Result` 别名。
pub type Result<T> = std::result::Result<T, MilpError>;

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error;

    /// Display 输出携带行号与出错说明。
    #[test]
    fn display_parse_contains_line_and_message() {
        let e = MilpError::Parse {
            line: 7,
            message: "未知变量名 `z9`".to_string(),
        };
        let s = e.to_string();
        assert!(s.contains('7'), "Display 应含行号：{s}");
        assert!(s.contains("z9"), "Display 应含出错说明：{s}");
    }

    /// 五个变体的 Display 均为人类可读中文。
    #[test]
    fn display_all_variants() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "no such file");
        assert!(MilpError::Io(io_err).to_string().contains("IO"));
        assert!(MilpError::InvalidModel("x".into())
            .to_string()
            .contains("模型非法"));
        assert!(MilpError::Numerical("x".into())
            .to_string()
            .contains("数值异常"));
        assert!(MilpError::Internal("x".into())
            .to_string()
            .contains("内部错误"));
    }

    /// source() 仅 Io 变体返回底层错误。
    #[test]
    fn source_only_for_io_variant() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "no such file");
        assert!(MilpError::Io(io_err).source().is_some());
        assert!(MilpError::Parse {
            line: 1,
            message: "m".into()
        }
        .source()
        .is_none());
        assert!(MilpError::InvalidModel("m".into()).source().is_none());
        assert!(MilpError::Numerical("m".into()).source().is_none());
        assert!(MilpError::Internal("m".into()).source().is_none());
    }

    /// `?` 语境下 std::io::Error 自动包装为 MilpError::Io。
    #[test]
    fn from_io_error_wraps() {
        fn fallible() -> Result<()> {
            Err(std::io::Error::new(
                std::io::ErrorKind::PermissionDenied,
                "denied",
            ))?;
            Ok(())
        }
        assert!(matches!(fallible(), Err(MilpError::Io(_))));
    }
}
