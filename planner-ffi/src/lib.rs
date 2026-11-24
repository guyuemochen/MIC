//! 通用规划器的本地接口。数字句柄在全局表中登记，释放不需要解引用调用方指针。
//! 求解时只锁住当前句柄，不阻塞独立模型。边界外不传递 Rust 对象。
use std::collections::HashMap;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::{Arc, Mutex, OnceLock};

#[derive(Default)]
struct Slot {
    output: String,
    failed: bool,
}

impl Slot {
    fn fail(&mut self, message: String) -> i32 {
        if !self.failed {
            self.output = message;
            self.failed = true;
        }
        1
    }
}

#[derive(Default)]
struct Handles {
    next: u64,
    slots: HashMap<u64, Arc<Mutex<Slot>>>,
}

fn registry() -> &'static Mutex<Handles> {
    static TABLE: OnceLock<Mutex<Handles>> = OnceLock::new();
    TABLE.get_or_init(Mutex::default)
}

fn with_slot<T>(id: u64, missing: T, f: impl FnOnce(&mut Slot) -> T) -> T {
    let slot = registry()
        .lock()
        .unwrap_or_else(|e| e.into_inner())
        .slots
        .get(&id)
        .cloned();
    match slot {
        None => missing,
        Some(slot) => f(&mut slot.lock().unwrap_or_else(|e| e.into_inner())),
    }
}

#[no_mangle]
pub extern "C" fn planner_abi_version() -> u32 {
    1
}

#[no_mangle]
pub extern "C" fn planner_create() -> u64 {
    catch_unwind(|| {
        let mut table = registry().lock().unwrap_or_else(|e| e.into_inner());
        let Some(id) = table.next.checked_add(1) else {
            return 0;
        };
        table.slots.insert(id, Arc::default());
        table.next = id;
        id
    })
    .unwrap_or(0)
}

#[no_mangle]
pub extern "C" fn planner_destroy(id: u64) {
    let _ = catch_unwind(|| {
        registry()
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .slots
            .remove(&id);
    });
}

/// # Safety
/// 非空输入指针须覆盖 len 个可读字节，并在本次调用期间保持有效。
#[no_mangle]
pub unsafe extern "C" fn planner_run(id: u64, data: *const u8, len: usize) -> i32 {
    catch_unwind(|| {
        with_slot(id, 2, |slot| {
            if slot.failed {
                return 1;
            }
            if len > isize::MAX as usize || (len > 0 && data.is_null()) {
                return slot.fail("输入缓冲区无效".into());
            }
            let result = catch_unwind(AssertUnwindSafe(|| {
                let bytes = if len == 0 {
                    &[]
                } else {
                    unsafe { std::slice::from_raw_parts(data, len) }
                };
                let input =
                    std::str::from_utf8(bytes).map_err(|_| "模型必须为 UTF-8".to_string())?;
                planner::business::wire::run(input)
            }));
            match result {
                Ok(Ok(output)) => {
                    slot.output = output;
                    0
                }
                Ok(Err(error)) => slot.fail(error),
                Err(_) => slot.fail("Rust 求解器内部 panic".into()),
            }
        })
    })
    .unwrap_or(2)
}

#[no_mangle]
pub extern "C" fn planner_size(id: u64) -> usize {
    catch_unwind(|| with_slot(id, usize::MAX, |s| s.output.len())).unwrap_or(usize::MAX)
}

/// # Safety
/// 非空输出指针须覆盖 capacity 个可写字节，调用期间由调用方独占。
#[no_mangle]
pub unsafe extern "C" fn planner_read(id: u64, out: *mut u8, capacity: usize) -> i32 {
    catch_unwind(|| {
        with_slot(id, 2, |slot| {
            let len = slot.output.len();
            if len == 0 {
                return 0;
            }
            if capacity < len || out.is_null() {
                return slot.fail("输出缓冲区无效或容量不足".into());
            }
            unsafe {
                std::ptr::copy_nonoverlapping(slot.output.as_ptr(), out, len);
            }
            0
        })
    })
    .unwrap_or(2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn independent_handles_and_caller_buffers() {
        let input = b"V 3 8\nE var 0\nO 0\nC 7 0 -1 1 tabu bf";
        let threads: Vec<_> = (0..4)
            .map(|_| {
                std::thread::spawn(move || {
                    let id = planner_create();
                    assert_ne!(id, 0);
                    let run = || unsafe { planner_run(id, input.as_ptr(), input.len()) };
                    for _ in 0..2 {
                        assert_eq!(run(), 0);
                        let mut bytes = vec![0; planner_size(id)];
                        assert_eq!(
                            unsafe { planner_read(id, bytes.as_mut_ptr(), bytes.len()) },
                            0
                        );
                        assert!(std::str::from_utf8(&bytes).unwrap().contains("\nV 3\n"));
                    }
                    assert_ne!(unsafe { planner_read(id, std::ptr::null_mut(), 0) }, 0);
                    assert_ne!(run(), 0);
                    planner_destroy(id);
                    assert_eq!(run(), 2);
                })
            })
            .collect();
        for worker in threads {
            worker.join().unwrap();
        }
    }

    #[test]
    fn error_is_sticky_and_destroy_is_idempotent() {
        let id = planner_create();
        let run = |bytes: &[u8]| unsafe { planner_run(id, bytes.as_ptr(), bytes.len()) };
        assert_ne!(run(&[255]), 0);
        let mut bytes = vec![0; planner_size(id)];
        assert_eq!(
            unsafe { planner_read(id, bytes.as_mut_ptr(), bytes.len()) },
            0
        );
        assert!(String::from_utf8(bytes.clone()).unwrap().contains("UTF-8"));
        assert_ne!(run(b"C 42 3 -1 1 tabu ff"), 0);
        assert_eq!(planner_size(id), bytes.len());
        planner_destroy(id);
        planner_destroy(id);
        assert_eq!(planner_size(id), usize::MAX);
    }
}
