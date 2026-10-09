//! milp-ffi C ABI 集成测试：仅经导出符号访问（不触碰 pub(crate) 内部）。
//!
//! 覆盖（对应 README 契约）：
//! ① 全流程（dive 21 / knapsack-10 17 / LP+offset 111）；
//! ② gap 向量（hard-knapsack-12，阈值 20%，事实由 milp CLI 实跑钉死：
//!    --mipgap 20 → gap_reached best=511 gap=3.11%；默认 → optimal 523）；
//! ③ 错误路径全码；
//! ④ 参数覆盖序列（0=关断编码钉死：TimeLimit/LogPeriodSecs 后设 0.0 真正关断）；
//! ⑤ last_error 截断/长度/粘滞；
//! ⑥ format_number 与 CLI objective 行逐字节一致（search_ctrl.rs:710-718 向量）；
//! ⑦ solve 后 get_solution 维数==n_solved、get_stat 全 id；
//! ⑧ 维数锚定（solve→add_var 后仍读旧维数，新变量取值→11 不中毒）；
//! ⑨ 重复 solve 覆盖。

use std::ffi::CString;
use std::os::raw::{c_char, c_void};
use std::ptr;

use milp_ffi::api::*;

const OK: i32 = 0;
const ERR_NULL: i32 = 1;
const ERR_INVALID_MODEL: i32 = 2;
const ERR_ARG: i32 = 9;
const ERR_UNKNOWN_PARAM: i32 = 10;
const ERR_NO_SOLUTION: i32 = 11;

// ===========================================================================
// 辅助
// ===========================================================================

fn cs(s: &str) -> CString {
    CString::new(s).expect("测试字符串不含内部 NUL")
}

/// 读取当前线程粘滞错误消息全文（两段式：先长度后内容）。
fn last_error_full() -> String {
    let need = milp_ffi_last_error(ptr::null_mut(), 0);
    assert!(need >= 0);
    let mut buf = vec![0u8; need as usize + 1];
    let need2 = milp_ffi_last_error(buf.as_mut_ptr() as *mut c_char, buf.len() as u64);
    assert_eq!(need, need2, "cap 足够时返回值不变");
    String::from_utf8(buf[..need as usize].to_vec()).expect("消息为 UTF-8")
}

fn expect(rc: i32, want: i32) {
    assert_eq!(
        rc, want,
        "期望返回码 {want}，实际 {rc}，last_error = {}",
        last_error_full()
    );
}

/// 创建无名模型句柄（sense：1=MIN, 2=MAX）。
fn create(sense: i32) -> *mut c_void {
    let m = milp_ffi_model_create(ptr::null(), sense);
    assert!(!m.is_null(), "create 应返回非空句柄");
    m
}

fn add_var(m: *mut c_void, name: &str, vtype: i32, lb: f64, ub: f64) -> u32 {
    let mut idx = 0u32;
    expect(
        milp_ffi_add_var(m, cs(name).as_ptr(), vtype, lb, ub, &mut idx),
        OK,
    );
    idx
}

#[allow(clippy::too_many_arguments)]
fn add_constraint(
    m: *mut c_void,
    name: Option<&str>,
    idx: &[u32],
    coeff: &[f64],
    relop: i32,
    rhs: f64,
) -> u32 {
    assert_eq!(idx.len(), coeff.len(), "测试辅助：数组等长");
    let mut out = 0u32;
    // CString 必须绑定变量：跨语句持有指针，避免悬垂。
    let name_c = name.map(cs);
    let name_ptr = name_c
        .as_ref()
        .map(|c| c.as_ptr())
        .unwrap_or_else(ptr::null);
    expect(
        milp_ffi_add_constraint(
            m,
            name_ptr,
            idx.as_ptr(),
            coeff.as_ptr(),
            idx.len() as u32,
            relop,
            rhs,
            &mut out,
        ),
        OK,
    );
    out
}

fn set_objective(m: *mut c_void, idx: &[u32], coeff: &[f64], sense: i32, offset: f64) {
    assert_eq!(idx.len(), coeff.len(), "测试辅助：数组等长");
    expect(
        milp_ffi_set_objective(
            m,
            idx.as_ptr(),
            coeff.as_ptr(),
            idx.len() as u32,
            sense,
            offset,
        ),
        OK,
    );
}

fn set_param_real(m: *mut c_void, name: &str, value: f64) -> i32 {
    milp_ffi_set_param_real(m, cs(name).as_ptr(), value)
}

fn set_param_int(m: *mut c_void, name: &str, value: i64) -> i32 {
    milp_ffi_set_param_int(m, cs(name).as_ptr(), value)
}

fn set_param_str(m: *mut c_void, name: &str, value: &str) -> i32 {
    milp_ffi_set_param_str(m, cs(name).as_ptr(), cs(value).as_ptr())
}

fn solve(m: *mut c_void) {
    expect(milp_ffi_solve(m), OK);
}

fn status(m: *mut c_void) -> i32 {
    let mut s = 0;
    expect(milp_ffi_get_status(m, &mut s), OK);
    s
}

fn objective(m: *mut c_void) -> f64 {
    let mut v = 0.0;
    expect(milp_ffi_get_objective(m, &mut v), OK);
    v
}

fn bound(m: *mut c_void) -> f64 {
    let mut v = 0.0;
    expect(milp_ffi_get_bound(m, &mut v), OK);
    v
}

fn gap_percent(m: *mut c_void) -> f64 {
    let mut v = 0.0;
    expect(milp_ffi_get_gap_percent(m, &mut v), OK);
    v
}

fn stat(m: *mut c_void, id: i32) -> i64 {
    let mut v = 0;
    expect(milp_ffi_get_stat(m, id, &mut v), OK);
    v
}

/// 取解向量：cap=0 探测维数（应报 9 并回写 n），再足量取回。
fn solution(m: *mut c_void) -> Vec<f64> {
    let mut n = 0u32;
    expect(
        milp_ffi_get_solution(m, ptr::null_mut(), 0, &mut n),
        ERR_ARG,
    );
    let len = n as usize;
    let mut buf = vec![0.0f64; len];
    let mut written = 0u32;
    expect(
        milp_ffi_get_solution(m, buf.as_mut_ptr(), len as u32, &mut written),
        OK,
    );
    assert_eq!(written as usize, len, "written == n_solved");
    buf
}

fn var_value(m: *mut c_void, index: u32) -> f64 {
    let mut v = 0.0;
    expect(milp_ffi_get_var_value(m, index, &mut v), OK);
    v
}

fn num_vars(m: *mut c_void) -> u32 {
    let mut n = 0;
    expect(milp_ffi_model_num_vars(m, &mut n), OK);
    n
}

/// dive 背包：max 10x1+x2, 20x1+3x2≤43, [0,10]² integer → OPTIMAL 21。
fn build_dive() -> *mut c_void {
    let m = create(2);
    add_var(m, "x1", 2, 0.0, 10.0);
    add_var(m, "x2", 2, 0.0, 10.0);
    add_constraint(m, None, &[0, 1], &[20.0, 3.0], 1, 43.0);
    set_objective(m, &[0, 1], &[10.0, 1.0], 0, 0.0);
    m
}

/// hard-knapsack-12（milp/instances/hard-knapsack-12.milp 等价）：
/// 12 个 BINARY，obj=[11,20,32,44,55,70,88,99,111,122,133,144]，
/// 重量 [7,11,13,17,19,23,29,31,37,41,43,47] ≤ 172 → 默认 OPTIMAL 523；
/// MIPGap=20% → LIMIT_REACHED best=511 gap≈3.11%。
fn build_hard12() -> *mut c_void {
    let obj = [
        11.0, 20.0, 32.0, 44.0, 55.0, 70.0, 88.0, 99.0, 111.0, 122.0, 133.0, 144.0,
    ];
    let weights = [
        7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0, 31.0, 37.0, 41.0, 43.0, 47.0,
    ];
    let m = create(2);
    for i in 0..12 {
        add_var(m, &format!("x{i}"), 3, 0.0, 1.0);
    }
    let idx: Vec<u32> = (0..12).collect();
    set_objective(m, &idx, &obj, 2, 0.0);
    add_constraint(m, Some("cap"), &idx, &weights, 1, 172.0);
    m
}

// ===========================================================================
// ① 全流程
// ===========================================================================

/// dive 背包：OPTIMAL 21，解 (2,1)，bound=21，gap 未触发→11。
#[test]
fn full_flow_dive_optimal_21() {
    let m = build_dive();
    assert_eq!(num_vars(m), 2);
    let mut nc = 0u32;
    expect(milp_ffi_model_num_constraints(m, &mut nc), OK);
    assert_eq!(nc, 1);

    solve(m);
    assert_eq!(status(m), 1, "OPTIMAL");
    assert!(
        (objective(m) - 21.0).abs() < 1e-9,
        "objective = {}",
        objective(m)
    );
    let x = solution(m);
    assert_eq!(x.len(), 2);
    assert!((x[0] - 2.0).abs() < 1e-9, "x1 = {}", x[0]);
    assert!((x[1] - 1.0).abs() < 1e-9, "x2 = {}", x[1]);
    assert!((bound(m) - 21.0).abs() < 1e-6, "完备搜索 bound == incumbent");
    let mut g = 0.0;
    expect(milp_ffi_get_gap_percent(m, &mut g), ERR_NO_SOLUTION);
    assert!((var_value(m, 0) - 2.0).abs() < 1e-9);
    assert!((var_value(m, 1) - 1.0).abs() < 1e-9);
    milp_ffi_model_destroy(m);
}

/// knapsack-10（milp/instances/knapsack-10.milp 等价）：OPTIMAL 17。
#[test]
fn knapsack_10_optimal_17() {
    let m = create(2);
    for i in 1..=5 {
        add_var(m, &format!("x{i}"), 3, 0.0, 1.0);
    }
    set_objective(m, &[0, 1, 2, 3, 4], &[7.0, 5.0, 9.0, 6.0, 3.0], 2, 0.0);
    add_constraint(
        m,
        Some("cap"),
        &[0, 1, 2, 3, 4],
        &[4.0, 3.0, 5.0, 6.0, 2.0],
        1,
        10.0,
    );
    solve(m);
    assert_eq!(status(m), 1);
    assert!((objective(m) - 17.0).abs() < 1e-9, "objective = {}", objective(m));
    milp_ffi_model_destroy(m);
}

/// LP + offset：max 3x+2y, x+y≤4, x+3y≤6, [0,3]², offset=100 → 111。
#[test]
fn lp_with_offset_111() {
    let m = create(2);
    add_var(m, "x", 1, 0.0, 3.0);
    add_var(m, "y", 1, 0.0, 3.0);
    add_constraint(m, None, &[0, 1], &[1.0, 1.0], 1, 4.0);
    add_constraint(m, None, &[0, 1], &[1.0, 3.0], 1, 6.0);
    set_objective(m, &[0, 1], &[3.0, 2.0], 0, 100.0);
    solve(m);
    assert_eq!(status(m), 1);
    assert!((objective(m) - 111.0).abs() < 1e-9, "objective = {}", objective(m));
    let x = solution(m);
    assert!((x[0] - 3.0).abs() < 1e-9 && (x[1] - 1.0).abs() < 1e-9);
    milp_ffi_model_destroy(m);
}

// ===========================================================================
// ② gap 向量
// ===========================================================================

/// MIPGap=20% → LIMIT_REACHED(4)、gap_terminated=1、objective≈511、
/// gap_percent ∈ (2.0, 4.0)；新句柄缺省 → OPTIMAL 523、gap_terminated=0。
/// （2026-10 割族增强后根 LP 终顶点漂移：旧 GMI 终顶点取整复核恰为 512；
/// 新割终顶点中面化，根 incumbent 改由割前顶点 dive 供应 = 511（CLI 实跑
/// 佐证 gap=3.11% / best=511.000000）。gap 区间断言 (2,4) 不变。）
#[test]
fn gap_vector_hard_knapsack_12() {
    let m = build_hard12();
    expect(set_param_real(m, "MIPGap", 20.0), OK);
    solve(m);
    assert_eq!(status(m), 4, "LIMIT_REACHED（gap 提前终止）");
    assert_eq!(stat(m, 6), 1, "gap_terminated");
    assert!(
        (objective(m) - 511.0).abs() <= 1e-6,
        "objective = {}",
        objective(m)
    );
    let g = gap_percent(m);
    assert!(g.is_finite() && g > 2.0 && g < 4.0, "gap = {g}");
    let b = bound(m);
    assert!(b.is_finite() && b > 511.0, "bound 证据 = {b}");
    milp_ffi_model_destroy(m);

    // 同一模型、新句柄缺省参数：完整证明最优。
    let m2 = build_hard12();
    solve(m2);
    assert_eq!(status(m2), 1, "OPTIMAL");
    assert_eq!(stat(m2, 6), 0, "无 gap 终止");
    assert!((objective(m2) - 523.0).abs() < 1e-9, "objective = {}", objective(m2));
    let mut g2 = 0.0;
    expect(milp_ffi_get_gap_percent(m2, &mut g2), ERR_NO_SOLUTION);
    milp_ffi_model_destroy(m2);
}

// ===========================================================================
// ③ 错误路径
// ===========================================================================

#[test]
fn error_paths_all_codes() {
    // NULL 句柄 → 1。
    let mut s: i32 = 0;
    expect(milp_ffi_get_status(ptr::null_mut(), &mut s), ERR_NULL);
    // NULL 输出指针 → 1。
    let m = create(1);
    expect(milp_ffi_get_status(m, ptr::null_mut()), ERR_NULL);
    // create 未知 sense → NULL。
    assert!(milp_ffi_model_create(ptr::null(), 7).is_null());
    // destroy NULL 幂等安全。
    milp_ffi_model_destroy(ptr::null_mut());

    let mut idx = 0u32;
    // 未知枚举 → 9。
    expect(
        milp_ffi_add_var(m, cs("x").as_ptr(), 7, 0.0, 1.0, &mut idx),
        ERR_ARG,
    );
    expect(
        milp_ffi_add_var(m, cs("x").as_ptr(), 1, f64::NEG_INFINITY, 1.0, &mut idx),
        ERR_ARG,
    );
    expect(
        milp_ffi_add_var(m, cs("x").as_ptr(), 1, 0.0, f64::NEG_INFINITY, &mut idx),
        ERR_ARG,
    );
    expect(
        milp_ffi_add_var(m, b"\xff\xfe\xfd\0".as_ptr() as *const c_char, 1, 0.0, 1.0, &mut idx),
        8, // ERR_UTF8
    );
    expect(
        milp_ffi_add_constraint(
            m,
            ptr::null(),
            [0u32].as_ptr(),
            [1.0f64].as_ptr(),
            1,
            9,
            1.0,
            &mut idx,
        ),
        ERR_ARG,
    );
    expect(
        milp_ffi_add_constraint(
            m,
            ptr::null(),
            [0u32].as_ptr(),
            [1.0f64].as_ptr(),
            1,
            1,
            f64::INFINITY,
            &mut idx,
        ),
        ERR_ARG,
    );
    // 建模语义 → 2。
    expect(milp_ffi_add_var(m, cs("x").as_ptr(), 1, 0.0, 1.0, &mut idx), OK);
    expect(
        milp_ffi_add_var(m, cs("x").as_ptr(), 1, 0.0, 1.0, &mut idx),
        ERR_INVALID_MODEL,
    );
    expect(
        milp_ffi_add_var(m, cs("max").as_ptr(), 1, 0.0, 1.0, &mut idx),
        ERR_INVALID_MODEL,
    );
    expect(
        milp_ffi_add_var(m, cs("y").as_ptr(), 1, 5.0, 4.0, &mut idx),
        ERR_INVALID_MODEL,
    );
    // 参数值域 → 9；类型不符 → 9；未知参数名 → 10。
    expect(set_param_real(m, "TimeLimit", -1.0), ERR_ARG);
    expect(set_param_real(m, "TimeLimit", f64::INFINITY), ERR_ARG);
    expect(set_param_real(m, "MIPGap", -1.0), ERR_ARG);
    expect(set_param_int(m, "Threads", 0), ERR_ARG);
    expect(set_param_int(m, "Threads", 1025), ERR_ARG);
    expect(set_param_str(m, "HeuristicMode", "fast"), ERR_ARG);
    expect(set_param_real(m, "Nope", 1.0), ERR_UNKNOWN_PARAM);
    expect(set_param_int(m, "TimeLimit", 5), ERR_ARG); // 存在但类型不符
    // 未求解 → 11；get_status 未求解 → 0。
    let mut f = 0.0f64;
    expect(milp_ffi_get_objective(m, &mut f), ERR_NO_SOLUTION);
    expect(milp_ffi_get_bound(m, &mut f), ERR_NO_SOLUTION);
    expect(milp_ffi_get_gap_percent(m, &mut f), ERR_NO_SOLUTION);
    let mut i64v = 0i64;
    expect(milp_ffi_get_stat(m, 1, &mut i64v), ERR_NO_SOLUTION);
    expect(milp_ffi_get_var_value(m, 0, &mut f), ERR_NO_SOLUTION);
    expect(milp_ffi_get_solution(m, ptr::null_mut(), 0, ptr::null_mut()), ERR_NO_SOLUTION);
    expect(milp_ffi_get_status(m, &mut s), OK);
    assert_eq!(s, 0, "NOT_SOLVED");
    milp_ffi_model_destroy(m);
}

// ===========================================================================
// ④ 参数覆盖序列（0=关断编码钉死）
// ===========================================================================

/// TimeLimit 5.0 后设 0.0 ⇒ 真正解除限速：硬算例完整求解 OPTIMAL 523。
/// LogPeriodSecs 2.5 后设 0.0 ⇒ 关闭钩子（仅断言 solve 返回 0）。
#[test]
fn zero_overrides_previous_time_limit_and_log() {
    let m = build_hard12();
    expect(set_param_real(m, "TimeLimit", 5.0), OK);
    expect(set_param_real(m, "TimeLimit", 0.0), OK);
    solve(m);
    assert_eq!(status(m), 1, "0.0 必须解除限速（完整证明最优）");
    assert!((objective(m) - 523.0).abs() < 1e-9, "objective = {}", objective(m));
    assert_eq!(stat(m, 7), 0, "非 auto 模式无自动降级");

    expect(set_param_real(m, "LogPeriodSecs", 2.5), OK);
    expect(set_param_real(m, "LogPeriodSecs", 0.0), OK);
    solve(m); // 无 stdout 断言，仅断言返回 0（helper expect）
    assert_eq!(status(m), 1);
    milp_ffi_model_destroy(m);
}

// ===========================================================================
// ⑤ last_error：截断、长度返回、粘滞
// ===========================================================================

#[test]
fn last_error_truncation_length_and_stickiness() {
    let m = create(1);
    expect(set_param_real(m, "Nope", 1.0), ERR_UNKNOWN_PARAM);
    let full = last_error_full();
    assert!(!full.is_empty());

    let need = milp_ffi_last_error(ptr::null_mut(), 0);
    assert_eq!(need as usize, full.as_bytes().len(), "返回所需字节数");

    // 截断：cap=4 → 写 3 字节 + NUL，返回值仍为全长。
    let mut small = [0u8; 4];
    let need2 = milp_ffi_last_error(small.as_mut_ptr() as *mut c_char, 4);
    assert_eq!(need2, need);
    assert_eq!(&small[..3], &full.as_bytes()[..3]);
    assert_eq!(small[3], 0, "截断也保证 NUL 结尾");

    // 成功调用不清除粘滞消息。
    assert_eq!(milp_ffi_abi_version(), 1u32);
    assert_eq!(
        milp_ffi_last_error(ptr::null_mut(), 0),
        need,
        "成功调用不得清除粘滞错误"
    );
    milp_ffi_model_destroy(m);
}

// ===========================================================================
// ⑥ format_number：与 CLI objective 行逐字节一致
// ===========================================================================

/// 向量逐字节断言（移植 milp/src/search_ctrl.rs:710-718 的
/// format_number_branches 测试向量）。
#[test]
fn format_number_matches_cli_vectors() {
    let vectors: &[(f64, &str)] = &[
        (0.0, "0.000000"),
        (13.0, "13.000000"),
        (-2.5, "-2.500000"),
        (1e-6, "0.000001"),
        (9.999999e-7, "9.999999e-7"),
        (-1e-7, "-1e-7"),
        (9.99e14, "999000000000000.000000"),
        (1e15, "1e15"),
        (-2.5e18, "-2.5e18"),
    ];
    for &(v, want) in vectors {
        let mut buf = [0u8; 64];
        let need = milp_ffi_format_number(v, buf.as_mut_ptr() as *mut c_char, 64);
        assert!(need > 0 && (need as usize) < 64);
        assert_eq!(&buf[..need as usize], want.as_bytes(), "v = {v}");
        assert_eq!(buf[need as usize], 0, "NUL 结尾");
    }
    // 长度探测：buf=NULL/cap=0 只报长度。
    let need = milp_ffi_format_number(1e15, ptr::null_mut(), 0);
    assert_eq!(need as usize, "1e15".len());
    // 非有限值防御性兜底（Display：NaN/inf/-inf），不 panic。
    for (v, want) in [(f64::NAN, "NaN"), (f64::INFINITY, "inf"), (f64::NEG_INFINITY, "-inf")] {
        let mut buf = [0u8; 32];
        let need = milp_ffi_format_number(v, buf.as_mut_ptr() as *mut c_char, 32);
        assert_eq!(&buf[..need as usize], want.as_bytes(), "v = {v}");
    }
}

// ===========================================================================
// ⑦ solve 后：解维数与全量统计
// ===========================================================================

#[test]
fn solution_dimension_and_all_stats_after_solve() {
    let m = build_dive();
    solve(m);
    let x = solution(m);
    assert_eq!(x.len(), num_vars(m) as usize, "n_solved == 当前变量数");

    // 容量不足 → 9 且 written 回写 n_solved。
    let mut written = 0u32;
    let mut buf1 = [0.0f64; 1];
    expect(
        milp_ffi_get_solution(m, buf1.as_mut_ptr(), 1, &mut written),
        ERR_ARG,
    );
    assert_eq!(written as usize, x.len(), "written_or_null 回写所需维数");

    for id in 1..=7 {
        let v = stat(m, id);
        assert!(v >= 0, "stat {id} = {v} 应非负");
    }
    assert!(stat(m, 1) >= 1, "nodes ≥ 1");
    assert_eq!(stat(m, 6), 0, "未触发 gap");
    assert_eq!(stat(m, 7), 0, "未自动降级");
    // 未知 stat_id → 9。
    let mut v: i64 = 0;
    expect(milp_ffi_get_stat(m, 0, &mut v), ERR_ARG);
    expect(milp_ffi_get_stat(m, 8, &mut v), ERR_ARG);
    milp_ffi_model_destroy(m);
}

// ===========================================================================
// ⑧ 维数锚定：solve → add_var → 旧维数仍可读、新变量→11 且不中毒
// ===========================================================================

#[test]
fn dimension_anchored_to_stored_solve() {
    let m = build_dive();
    solve(m);
    add_var(m, "z", 2, 0.0, 5.0); // 结构性修改：当前 vars.len()=3

    let x = solution(m);
    assert_eq!(x.len(), 2, "get_solution 仍返回旧维数 n_solved");
    assert!((x[0] - 2.0).abs() < 1e-9);

    // 新变量（index=2 < vars.len()=3 但 ≥ n_solved=2）→ 11，不 panic。
    let mut v = 0.0f64;
    expect(milp_ffi_get_var_value(m, 2, &mut v), ERR_NO_SOLUTION);
    // 句柄未中毒：其余入口照常工作。
    assert_eq!(status(m), 1);
    assert!((var_value(m, 0) - 2.0).abs() < 1e-9);
    // 下标越界（≥ 当前 vars.len()=3）→ 9。
    expect(milp_ffi_get_var_value(m, 3, &mut v), ERR_ARG);
    milp_ffi_model_destroy(m);
}

// ===========================================================================
// ⑨ 重复 solve 覆盖
// ===========================================================================

#[test]
fn repeated_solve_overwrites_result() {
    let m = build_dive();
    solve(m);
    assert!((objective(m) - 21.0).abs() < 1e-9);
    // 加约束 x1 ≤ 1 后重解：max 10x1+x2 → 10 + 7 = 17。
    add_constraint(m, Some("cap_x1"), &[0], &[1.0], 1, 1.0);
    solve(m);
    assert_eq!(status(m), 1);
    assert!(
        (objective(m) - 17.0).abs() < 1e-9,
        "结果应为最近一次 solve（objective = {}）",
        objective(m)
    );
    let x = solution(m);
    assert!((x[0] - 1.0).abs() < 1e-9 && (x[1] - 7.0).abs() < 1e-9);
    milp_ffi_model_destroy(m);
}

// ===========================================================================
// 补充：建模 API 与参数 happy path、并行路径
// ===========================================================================

/// set_var_bounds 生效与校验（lb>ub→2、越界→9、lb=-∞→9）。
#[test]
fn set_var_bounds_updates_and_validates() {
    let m = create(2);
    let i = add_var(m, "x", 1, 0.0, 10.0);
    expect(milp_ffi_set_var_bounds(m, i, 1.0, 5.0), OK);
    expect(milp_ffi_set_var_bounds(m, i, 5.0, 1.0), ERR_INVALID_MODEL);
    expect(milp_ffi_set_var_bounds(m, 99, 0.0, 1.0), ERR_ARG);
    expect(
        milp_ffi_set_var_bounds(m, i, f64::NEG_INFINITY, 1.0),
        ERR_ARG,
    );
    // ub=+∞ 合法。
    expect(milp_ffi_set_var_bounds(m, i, 1.0, f64::INFINITY), OK);
    expect(milp_ffi_set_var_bounds(m, i, 1.0, 5.0), OK);
    set_objective(m, &[i], &[1.0], 2, 0.0);
    solve(m);
    assert!((objective(m) - 5.0).abs() < 1e-9);
    milp_ffi_model_destroy(m);
}

/// 参数 setter happy path：全部参数合法值均返回 0。
#[test]
fn param_setters_happy_paths() {
    let m = create(2);
    expect(set_param_real(m, "TimeLimit", 0.0), OK);
    expect(set_param_real(m, "MIPGap", 0.0), OK);
    expect(set_param_real(m, "LogPeriodSecs", 0.0), OK);
    expect(set_param_int(m, "Threads", 1), OK);
    expect(set_param_int(m, "Threads", 1024), OK);
    for v in ["auto", "on", "off"] {
        expect(set_param_str(m, "HeuristicMode", v), OK);
    }
    milp_ffi_model_destroy(m);
}

/// 并行路径（Threads=2）经 SearchControl 透传：dive 仍精确最优 21。
#[test]
fn threads_two_parallel_path_still_exact() {
    let m = build_dive();
    expect(set_param_int(m, "Threads", 2), OK);
    solve(m);
    assert_eq!(status(m), 1);
    assert!((objective(m) - 21.0).abs() < 1e-9);
    milp_ffi_model_destroy(m);
}

/// 命名 create 与带名约束：name 非 NULL 正常入模。
#[test]
fn named_create_and_named_constraint() {
    let m = milp_ffi_model_create(cs("my_model").as_ptr(), 1);
    assert!(!m.is_null());
    let i = add_var(m, "x", 1, 0.0, f64::INFINITY);
    add_constraint(m, Some("c0"), &[i], &[2.0], 3, 7.0);
    set_objective(m, &[i], &[1.0], 1, 0.0);
    solve(m);
    assert_eq!(status(m), 1);
    assert!((objective(m) - 3.5).abs() < 1e-9);
    milp_ffi_model_destroy(m);
}

/// B2 回归：set_objective 越界失败（→2）时 sense 不得被半突变。
///
/// 场景（ctypes 实证同款）：MINIMIZE 模型 x∈[0,1]，先以 sense=2(MAXIMIZE)
/// 传越界下标 7 → 返回 2；旧实现在校验前已写 problem.sense=Max。随后以
/// sense=0（保持）合法设目标 min x 并求解——sense 若被翻为 MAX 则 x 取
/// 上界 1.0，保持 MIN 则取下界 0.0。
#[test]
fn set_objective_failure_leaves_sense_untouched() {
    let m = create(1); // MINIMIZE
    let i = add_var(m, "x", 1, 0.0, 1.0);
    let bad_idx = [7u32];
    let coef = [1.0f64];
    expect(
        milp_ffi_set_objective(
            m,
            bad_idx.as_ptr(),
            coef.as_ptr(),
            bad_idx.len() as u32,
            2, // MAXIMIZE：若先写 sense，失败后模型半突变
            0.0,
        ),
        ERR_INVALID_MODEL,
    );
    // 合法重设目标（sense=0 保持现状）并求解。
    set_objective(m, &[i], &[1.0], 0, 0.0);
    solve(m);
    assert_eq!(status(m), 1);
    assert!(
        (objective(m) - 0.0).abs() < 1e-9,
        "MINIMIZE 必须保持：objective = {}（1.0 即 sense 已被翻为 MAXIMIZE）",
        objective(m)
    );
    milp_ffi_model_destroy(m);
}
