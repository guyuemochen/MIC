#!/usr/bin/env python3
"""离线核对第一阶段记录的向量、界与约束。仅使用 Python 标准库。"""
import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import sys


ENCODING = "ieee754_binary64_hex_bits"
SLOTS = ("entry", "first_candidate", "first_failure", "before_cleanup")
FAILURE_KINDS = {"failure", "failed_recompute", "partial_refresh"}
TOLERANCES = ("feasibility", "optimality", "pivot", "zero", "integrality")
HEX64 = re.compile(r"[0-9a-fA-F]{16}\Z")


class InvalidCapture(ValueError):
    pass


def require(condition, path, reason):
    if not condition:
        raise InvalidCapture(f"{path}: {reason}")


def obj(value, path):
    require(isinstance(value, dict), path, "应为 JSON 对象")
    return value


def array(value, path, length=None):
    require(isinstance(value, list), path, "应为 JSON 数组")
    require(length is None or len(value) == length, path, f"应有 {length} 项，实际 {len(value)} 项")
    return value


def integer(value, path, limit=None):
    require(type(value) is int and value >= 0, path, "应为非负整数")
    require(limit is None or value < limit, path, f"下标 {value} 超出范围 0..{limit}")
    return value


def boolean(value, path):
    require(type(value) is bool, path, "应为布尔值")
    return value


@dataclass(frozen=True)
class Binary64:
    bits: str
    value: float


def decode(value, path):
    require(isinstance(value, str) and HEX64.fullmatch(value), path, "应为 16 位十六进制 binary64 位串")
    return Binary64(value, struct.unpack(">d", bytes.fromhex(value))[0])


def vector(raw, path, length):
    return [decode(value, f"{path}[{i}]") for i, value in enumerate(array(raw, path, length))]


def finite(number, path, allow_upper_inf=False):
    value = number.value
    require(math.isfinite(value) or (allow_upper_inf and value == math.inf), path,
            f"不允许非有限值 {number.bits}")
    return value


def relation(value, path):
    require(value in ("Le", "Ge", "Eq"), path, "未知行关系")
    return value


def terms(raw, path, n):
    result = []
    for k, pair in enumerate(array(raw, path)):
        at = f"{path}[{k}]"
        col, coeff = array(pair, at, 2)
        result.append((integer(col, at, n), finite(decode(coeff, at), at)))
    return result


def bounds(low, high, path):
    lo = finite(low, path)
    hi = finite(high, path, allow_upper_inf=True)
    require(lo <= hi, path, "上下界颠倒")
    return lo, hi


@dataclass
class Row:
    terms: list
    rhs: float
    op: str


@dataclass
class InputModel:
    columns: list
    rows: list


def read_model(raw, path):
    model = obj(raw, path)
    variables = array(model.get("variables"), f"{path}.variables")
    n = len(variables)
    offset_path = f"{path}.offset"
    finite(decode(model.get("offset"), offset_path), offset_path)
    terms(model.get("objective"), f"{path}.objective", n)
    columns = []
    for j, raw_var in enumerate(variables):
        at = f"{path}.variables[{j}]"
        var = obj(raw_var, at)
        columns.append(bounds(decode(var.get("lower"), at), decode(var.get("upper"), at), at))
    rows = []
    for i, raw_row in enumerate(array(model.get("rows"), f"{path}.rows")):
        at = f"{path}.rows[{i}]"
        row = obj(raw_row, at)
        rows.append(Row(terms(row.get("terms"), f"{at}.terms", n),
                        finite(decode(row.get("rhs"), f"{at}.rhs"), f"{at}.rhs"),
                        relation(row.get("op"), f"{at}.op")))
    return InputModel(columns, rows)


@dataclass
class Matrix:
    n: int
    ncs: int
    rows: list
    scales: list
    slack: list
    artificial: list
    constraints: list
    cuts: list
    sources: list
    artificial_by_row: dict

    @property
    def column_count(self):
        return self.ncs + len(self.artificial)


def read_matrix(raw):
    path = "normalized_matrix"
    matrix = obj(raw, path)
    n = integer(matrix.get("structural_columns"), f"{path}.structural_columns")
    ncs = integer(matrix.get("static_columns"), f"{path}.static_columns")
    require(n <= ncs, path, "结构列多于静态列")
    row_terms = array(matrix.get("row_terms"), f"{path}.row_terms")
    m = len(row_terms)
    rhs = vector(matrix.get("row_rhs"), f"{path}.row_rhs", m)
    ops = array(matrix.get("row_op"), f"{path}.row_op", m)
    scales = vector(matrix.get("row_scale"), f"{path}.row_scale", m)
    rows = []
    for i in range(m):
        rows.append(Row(terms(row_terms[i], f"{path}.row_terms[{i}]", n),
                        finite(rhs[i], f"{path}.row_rhs[{i}]"), relation(ops[i], f"{path}.row_op[{i}]")))
    scales = [finite(scale, f"{path}.row_scale[{i}]") for i, scale in enumerate(scales)]
    require(all(scale > 0 for scale in scales), path, "归一尺度必须为正")
    slack = array(matrix.get("row_slack"), f"{path}.row_slack", m)
    used_slack = set()
    for i, (col, row) in enumerate(zip(slack, rows)):
        at = f"{path}.row_slack[{i}]"
        if row.op == "Eq":
            require(col is None, at, "等式行没有松弛列")
        else:
            col = integer(col, at, ncs)
            require(col >= n and col not in used_slack, at, "松弛列无效或重复")
            used_slack.add(col)
    require(len(used_slack) == ncs - n, path, "静态列段存在无来源松弛列")
    artificial = []
    by_row = {}
    for t, pair in enumerate(array(matrix.get("artificial_columns"), f"{path}.artificial_columns")):
        at = f"{path}.artificial_columns[{t}]"
        i, sigma = array(pair, at, 2)
        i = integer(i, at, m)
        sigma = finite(decode(sigma, at), at)
        require(abs(sigma) == 1 and i not in by_row, at, "人工列系数须为 ±1，且同一行至多一列")
        by_row[i] = (ncs + t, sigma)
        artificial.append((i, sigma))
    sources = [None] * m
    constraints = array(matrix.get("constraint_row"), f"{path}.constraint_row")
    cuts = array(matrix.get("cut_row"), f"{path}.cut_row")
    for kind, mapping in (("constraint", constraints), ("cut", cuts)):
        for j, i in enumerate(mapping):
            if i is None:
                continue
            i = integer(i, f"{path}.{kind}_row[{j}]", m)
            require(sources[i] is None, path, f"归一行 {i} 的来源重复")
            require(kind != "cut" or scales[i] == 1.0, path, f"割行 {i} 的归一尺度应为 1")
            sources[i] = {"kind": kind, "index": j}
    require(all(source is not None for source in sources), path, "有归一行缺少约束或割来源")
    return Matrix(n, ncs, rows, scales, slack, artificial, constraints, cuts, sources, by_row)


def check_input_mapping(model, matrix, zero):
    require(len(model.columns) == matrix.n, "engine_input_lp", "结构列数与归一矩阵不符")
    require(len(model.rows) == len(matrix.constraints), "engine_input_lp", "行数与约束来源映射不符")
    scales = []
    mismatches = []
    for ci, (row, i) in enumerate(zip(model.rows, matrix.constraints)):
        g = max((abs(a) for _, a in row.terms), default=0.0)
        scales.append(g)
        if i is None:
            if g > zero:
                mismatches.append({"constraint": ci, "reason": "被省略的输入行超过近零阈值"})
            continue
        mapped = matrix.rows[i]
        if g <= zero:
            mismatches.append({"constraint": ci, "normalized_row": i, "reason": "近零输入行被映射入表"})
            continue
        expected = sorted((j, a / g) for j, a in row.terms)
        actual = sorted(mapped.terms)
        if (row.op != mapped.op or expected != actual or row.rhs / g != mapped.rhs
                or 1.0 / g != matrix.scales[i]):
            mismatches.append({"constraint": ci, "normalized_row": i, "reason": "行关系、列系数、右端或尺度与归一结果不符"})
    return scales, mismatches


def metric(value):
    """超出 binary64 的有限结果用精确分数表示。"""
    if value is None:
        return None
    try:
        number = float(value)
    except OverflowError:
        return {"overflow_binary64": True, "exact_fraction": str(value)}
    if not math.isfinite(number):
        return {"nonfinite": True, "bits": struct.pack(">d", number).hex()}
    if number == 0.0 and value != 0:
        return {"underflow_binary64": True, "exact_fraction": str(value)}
    return number


def safe_fsum(values):
    try:
        return math.fsum(values)
    except (OverflowError, ValueError):
        return math.nan


@dataclass
class Residual:
    value: object
    denominator: object
    exact: bool
    arithmetic_issue: bool

    def backward_error(self):
        if self.denominator == 0:
            return 0.0 if self.value == 0 else None
        return abs(self.value) / self.denominator


def residual(coefficients, x, rhs, threshold):
    # 保留 RHS 作为同一次求和的一项。严重抵消、临界与溢出才做精确复核。
    products = [a * x[j] for j, a in coefficients]
    value = safe_fsum([*products, -rhs])
    denominator = safe_fsum([*(abs(p) for p in products), abs(rhs)])
    overflow = not (math.isfinite(value) and math.isfinite(denominator))
    # 包含乘法舍入和次正规下溢的保守筛选带；它只选择复核行。
    roundoff = 8 * sys.float_info.epsilon * denominator + len(products) * math.ulp(0.0)
    near = False
    cancellation = False
    if not overflow:
        magnitude = Fraction.from_float(abs(value)) if isinstance(threshold, Fraction) else abs(value)
        near = abs(magnitude - threshold) <= roundoff
        cancellation = (denominator > 0 and abs(value) <= denominator * 1e-10
                        and roundoff >= threshold / 1000)
    if not (overflow or near or cancellation):
        return Residual(value, denominator, False, False)
    exact_products = [Fraction.from_float(a) * Fraction.from_float(x[j]) for j, a in coefficients]
    exact_rhs = Fraction.from_float(rhs)
    return Residual(sum(exact_products, -exact_rhs),
                    sum(map(abs, exact_products), abs(exact_rhs)), True, overflow)


def violation(value, op):
    if op == "Eq":
        return abs(value)
    return max(0.0, value if op == "Le" else -value)


def measure(entries, tolerance):
    # 输入顺序用于并列破键，行号和来源始终和对应指标一起返回。
    worst = max(entries, key=lambda entry: entry[1], default=None)
    bad = sum(value > tolerance for _, value in entries)
    return {"passed": bad == 0, "violations": bad,
            "max_violation": metric(worst[1]) if worst else 0.0,
            "worst": worst[0] if worst else None}


def bound_violation(value, lo, hi):
    if lo <= value <= hi:
        return 0.0
    left, right = (lo, value) if value < lo else (value, hi)
    delta = left - right
    return delta if math.isfinite(delta) else Fraction.from_float(left) - Fraction.from_float(right)


def bounds_report(values, limits, tolerance):
    entries = [({"column": j, "value": number.bits}, bound_violation(number.value, *limit))
               for j, (number, limit) in enumerate(zip(values, limits))]
    return measure(entries, tolerance)


def read_summary(raw, path):
    raw = obj(raw, path)
    require(isinstance(raw.get("kind"), str), path, "摘要种类必须为字符串")
    for field in ("state_revision", "iterations", "call_iterations", "eta_count", "artificial_basic",
                  "artificial_nonzero", "nonfinite"):
        integer(raw.get(field), f"{path}.{field}")
    require(type(raw.get("phase")) is int and raw["phase"] == 1, path, "摘要的阶段必须为 1")
    boolean(raw.get("bland"), f"{path}.bland")
    for field in ("positive_sum", "negative_sum", "minimum", "maximum", "basic_bound_violation"):
        if field in ("minimum", "maximum") and field in raw and raw[field] is None:
            continue
        decode(raw.get(field), f"{path}.{field}")
    pricing = obj(raw.get("pricing"), f"{path}.pricing")
    require(pricing.get("kind") in ("", "optimal", "candidate", "failed"), path, "定价结果种类无效")
    for field in ("column", "evaluated_revision"):
        require(field in pricing, path, f"缺少定价字段 {field}")
        if pricing[field] is not None:
            integer(pricing[field], f"{path}.pricing.{field}")
    require("direction" in pricing, path, "缺少定价方向")
    require(pricing["evaluated_revision"] is None or pricing["evaluated_revision"] <= raw["state_revision"],
            path, "定价世代晚于当前状态世代")
    if pricing["kind"] == "candidate":
        require(pricing["column"] is not None and pricing["evaluated_revision"] is not None,
                path, "定价候选缺少列号或状态世代")
        direction = finite(decode(pricing["direction"], f"{path}.pricing.direction"), path)
        require(abs(direction) == 1.0, path, "定价方向必须为 ±1")
    else:
        require(pricing["column"] is None and pricing["direction"] is None,
                path, "无定价候选却记录了列号或方向")
    return raw


def original_sum(numbers, positive):
    total = 0.0
    for number in numbers:
        if (number.value >= 0) == positive:
            total += number.value
    return total


def artificial_report(values, basic_artificial, summary, tolerance):
    finite_values = [number.value for number in values if math.isfinite(number.value)]
    bad = [{"index": i, "bits": number.bits} for i, number in enumerate(values)
           if not math.isfinite(number.value)]
    finite_basic = [number.value for number in basic_artificial if math.isfinite(number.value)]
    sums = {}
    sum_overflow = []
    if not bad:
        for name, components in (("positive_sum", [v for v in finite_values if v > 0]),
                                 ("negative_sum", [v for v in finite_values if v < 0]),
                                 ("absolute_sum", [abs(v) for v in finite_values])):
            value = safe_fsum(components)
            if not math.isfinite(value):
                sum_overflow.append(name)
                value = sum(map(Fraction.from_float, components), Fraction())
            elif name == "positive_sum" and abs(value - tolerance) <= math.ulp(value):
                value = sum(map(Fraction.from_float, components), Fraction())
            sums[name] = value
    checks = {"artificial_basic": len(basic_artificial),
              "artificial_nonzero": sum(number.value != 0.0 for number in basic_artificial)}
    if len(finite_basic) == len(basic_artificial):
        checks.update(minimum=min(finite_basic, default=None), maximum=max(finite_basic, default=None),
                      positive_sum=original_sum(basic_artificial, True),
                      negative_sum=original_sum(basic_artificial, False))
    differences = []
    for name, expected in checks.items():
        saved = summary[name]
        saved_value = decode(saved, name).value if isinstance(saved, str) else saved
        if saved_value != expected and not (isinstance(expected, float) and math.isnan(expected)
                                            and isinstance(saved_value, float) and math.isnan(saved_value)):
            differences.append({"field": name, "saved": saved, "recomputed": metric(expected)})
    near_zero = not bad and sums["positive_sum"] <= tolerance and all(v >= -tolerance for v in finite_values)
    return {**{name: metric(sums.get(name)) for name in ("positive_sum", "negative_sum", "absolute_sum")},
            "minimum": metric(min(finite_values, default=None)) if not bad else None,
            "maximum": metric(max(finite_values, default=None)) if not bad else None,
            "nonfinite": bad, "sum_overflow": sum_overflow, "near_zero": near_zero,
            "each_within_tolerance": not bad and all(abs(v) <= tolerance for v in finite_values),
            "summary_comparable": len(finite_basic) == len(basic_artificial),
            "saved_summary_bits": {name: summary[name] for name in
                                   ("positive_sum", "negative_sum", "minimum", "maximum")},
            "summary_differences": differences}


def row_checks(matrix, model, input_scales, values, tolerance):
    x = [number.value for number in values]
    normalized = []
    augmented = []
    originals = []
    original_raw = []
    skipped = []
    augmented_rows = []
    backward_errors = []
    original_rows = []
    exact_count = 0
    arithmetic_issues = []

    def calculate(coefficients, rhs, threshold, location):
        nonlocal exact_count
        result = residual(coefficients, x, rhs, threshold)
        exact_count += result.exact
        if result.arithmetic_issue:
            arithmetic_issues.append(location)
        return result

    for i, row in enumerate(matrix.rows):
        location = {"row": i, "source": matrix.sources[i]}
        base = calculate(row.terms, row.rhs, tolerance, {"layer": "normalized", **location})
        normalized.append((location, violation(base.value, row.op)))
        extras = []
        slack = matrix.slack[i]
        if slack is not None:
            extras.append((slack, 1.0 if row.op == "Le" else -1.0))
        if i in matrix.artificial_by_row:
            extras.append(matrix.artificial_by_row[i])
        aug = (calculate(row.terms + extras, row.rhs, tolerance, {"layer": "augmented", **location})
               if extras else base)
        augmented.append((location, abs(aug.value)))
        backward_error = aug.backward_error()
        if backward_error is not None:
            backward_errors.append(backward_error)
        augmented_rows.append({**location, "absolute_residual": metric(abs(aug.value)),
                               "relative_backward_error": metric(backward_error),
                               "exact_recheck": aug.exact})

    for ci, (row, g, i) in enumerate(zip(model.rows, input_scales, matrix.constraints)):
        location = {"constraint": ci, "normalized_row": i}
        threshold = tolerance if i is None else Fraction.from_float(g) * Fraction.from_float(tolerance)
        result = calculate(row.terms, row.rhs, threshold, {"layer": "engine_input", **location})
        raw = violation(result.value, row.op)
        original_raw.append((location, raw))
        scaled = None
        if i is None:
            skipped.append((location, raw))
            passed = raw <= tolerance
        else:
            # 在原尺度精确比较，避免分母或阈值的双精度溢出、下溢改变判定。
            scaled = Fraction(raw) / Fraction.from_float(g) if g > 0 else None
            if scaled is not None:
                originals.append((location, scaled))
            passed = scaled is not None and raw <= threshold
        original_rows.append({**location, "raw_violation": metric(raw), "scaled_violation": metric(scaled),
                              "passed": passed, "exact_recheck": result.exact})

    raw_worst = max(original_raw, key=lambda item: item[1], default=None)
    scaled_report = measure(originals, tolerance)
    scaled_report["violations"] = sum(not row["passed"] for row in original_rows if row["normalized_row"] is not None)
    scaled_report["passed"] = scaled_report["violations"] == 0
    augmented_report = measure(augmented, tolerance)
    augmented_report["max_relative_backward_error"] = metric(max(backward_errors, default=0.0))
    augmented_report["rows"] = augmented_rows
    return {"normalized": measure(normalized, tolerance),
            "cuts": measure([entry for entry in normalized if entry[0]["source"]["kind"] == "cut"], tolerance),
            "augmented": augmented_report,
            "engine_input": {"raw_max_violation": metric(raw_worst[1]) if raw_worst else 0.0,
                             "raw_worst": raw_worst[0] if raw_worst else None,
                             "scaled": scaled_report, "dropped_near_zero": measure(skipped, tolerance),
                             "rows": original_rows},
            "exact_recheck_count": exact_count, "binary64_arithmetic_issues": arithmetic_issues}


def inspect_snapshot(raw, matrix, model, input_scales, tolerance):
    kind = raw["kind"]
    summary = read_summary(raw.get("summary"), f"{kind}.summary")
    m = len(matrix.rows)
    n = matrix.n
    total = matrix.column_count
    basis = [integer(j, f"{kind}.basis[{k}]", total)
             for k, j in enumerate(array(raw.get("basis"), f"{kind}.basis", m))]
    require(len(set(basis)) == len(basis), kind, "基列重复")
    basic = vector(raw.get("basic_values"), f"{kind}.basic_values", m)
    lower = vector(raw.get("lower"), f"{kind}.lower", n)
    upper = vector(raw.get("upper"), f"{kind}.upper", n)
    saved_x = vector(raw.get("x_unclipped"), f"{kind}.x_unclipped", n)
    at_upper = [boolean(value, f"{kind}.at_upper[{j}]")
                for j, value in enumerate(array(raw.get("at_upper"), f"{kind}.at_upper", total))]
    broken = boolean(raw.get("factorization_broken"), f"{kind}.factorization_broken")
    if summary["pricing"]["column"] is not None:
        integer(summary["pricing"]["column"], f"{kind}.summary.pricing.column", matrix.ncs)
    require("failed_recomputed_values" in raw and "failed_refresh_row_value" in raw, kind, "缺少失败重算字段")
    failed_values = raw["failed_recomputed_values"]
    if failed_values is not None:
        failed_values = [v.bits for v in vector(failed_values, f"{kind}.failed_recomputed_values", m)]
    failed_refresh = raw["failed_refresh_row_value"]
    if failed_refresh is not None:
        i, value = array(failed_refresh, f"{kind}.failed_refresh_row_value", 2)
        failed_refresh = [integer(i, f"{kind}.failed_refresh_row_value", m),
                          decode(value, f"{kind}.failed_refresh_row_value").bits]
    require(kind != "failed_recompute" or failed_values is not None, kind, "缺少失败重算的新向量")
    require(kind != "partial_refresh" or failed_refresh is not None, kind, "缺少刷新失败的行和值")
    zero = Binary64("0000000000000000", 0.0)
    inf = Binary64("7ff0000000000000", math.inf)
    values = [upper[j] if at_upper[j] else lower[j] for j in range(n)]
    values.extend(inf if at_upper[j] else zero for j in range(n, total))
    for j, value in zip(basis, basic):
        values[j] = value
    mismatch = [{"column": j, "saved": saved.bits, "reconstructed": value.bits}
                for j, (saved, value) in enumerate(zip(saved_x, values))
                if saved.bits.lower() != value.bits.lower()]
    illegal = [{"column": j, "bits": number.bits} for j, number in enumerate(values)
               if not math.isfinite(number.value)]
    bound_issues = []
    current = []
    for j, (lo, hi) in enumerate(zip(lower, upper)):
        try:
            current.append(bounds(lo, hi, f"{kind}.bounds[{j}]"))
        except InvalidCapture as error:
            bound_issues.append({"column": j, "lower": lo.bits, "upper": hi.bits, "reason": str(error)})
    artificial = artificial_report(values[matrix.ncs:],
                                   [value for j, value in zip(basis, basic) if j >= matrix.ncs],
                                   summary, tolerance)
    bad_basic = sum(not math.isfinite(value.value) for value in basic)
    if bad_basic != summary["nonfinite"]:
        artificial["summary_differences"].append({"field": "nonfinite", "saved": summary["nonfinite"],
                                                  "recomputed": bad_basic})
    saved_basic_violation = decode(summary["basic_bound_violation"], "basic_bound_violation")
    result = {"kind": kind, "call_iterations": summary["call_iterations"],
              "iterations": summary["iterations"], "state_revision": summary["state_revision"],
              "factorization_broken": broken, "failed_recomputed_values": failed_values,
              "failed_refresh_row_value": failed_refresh, "x_unclipped_bits": [number.bits for number in saved_x],
              "x_bit_mismatches": mismatch, "nonfinite_reconstructed_columns": illegal,
              "invalid_bounds": bound_issues, "artificial": artificial, "issues": [],
              "saved_basic_bound_violation_bits": saved_basic_violation.bits,
              "reconstructed_column_bits": [number.bits for number in values],
              "current_structural_bounds": None, "current_augmented_bounds": None,
              "construction_bounds": None, "rows": None}
    if illegal:
        result["issues"].append("重建列含非有限值，包括贴在 +∞ 上界的非基列")
    if mismatch:
        result["issues"].append("重建结构向量和保存向量的原始位不一致")
    if bound_issues:
        result["issues"].append("当前界非法或非有限")
    if artificial["summary_differences"]:
        result["issues"].append("人工量或基本值非有限计数与摘要不符")
    if not math.isfinite(saved_basic_violation.value):
        result["issues"].append("保存的基本变量违界汇总非有限")
    if kind in ("partial_refresh", "failed_recompute"):
        result["issues"].append("失败重算包含部分更新或独立的新失败向量")
    if not (illegal or bound_issues):
        result["current_structural_bounds"] = bounds_report(values[:n], current, tolerance)
        result["current_augmented_bounds"] = bounds_report(values, current + [(0.0, math.inf)] * (total - n), tolerance)
        result["construction_bounds"] = bounds_report(values[:n], model.columns, tolerance)
        result["construction_bounds"]["tightened_columns"] = sum(
            lo > old_lo or hi < old_hi for (lo, hi), (old_lo, old_hi) in zip(current, model.columns))
        result["construction_bounds"]["relaxed_columns"] = sum(
            lo < old_lo or hi > old_hi for (lo, hi), (old_lo, old_hi) in zip(current, model.columns))
        result["rows"] = row_checks(matrix, model, input_scales, values, tolerance)
    return result


def classify(snapshot, missing):
    if missing or snapshot["issues"] or snapshot["rows"] is None:
        return "insufficient_evidence"
    rows = snapshot["rows"]
    if snapshot["artificial"]["sum_overflow"] or rows["binary64_arithmetic_issues"]:
        return "insufficient_evidence"
    if not snapshot["artificial"]["near_zero"]:
        return "artificial_not_near_zero"
    if not (snapshot["current_structural_bounds"]["passed"] and snapshot["construction_bounds"]["passed"]
            and rows["normalized"]["passed"]
            and rows["engine_input"]["scaled"]["passed"] and rows["engine_input"]["dropped_near_zero"]["passed"]):
        return "near_zero_but_primal_infeasible"
    if (snapshot["factorization_broken"] or not snapshot["current_augmented_bounds"]["passed"]
            or not rows["augmented"]["passed"]):
        return "primal_pass_but_augmented_inconsistent"
    return "primal_and_augmented_pass"


def read_slots(capture):
    require(capture.get("slot_order") == list(SLOTS), "slot_order", "槽名或顺序与 schema 1 不符")
    omitted = array(capture.get("slot_omitted"), "slot_omitted", len(SLOTS))
    attempts = array(capture.get("slot_attempt_iterations"), "slot_attempt_iterations", len(SLOTS))
    result = {}
    for i, slot in enumerate(SLOTS):
        skipped = boolean(omitted[i], f"slot_omitted[{i}]")
        attempt = attempts[i]
        if attempt is not None:
            integer(attempt, f"slot_attempt_iterations[{i}]")
        require(not skipped or attempt is not None, slot, "已省略却缺少尝试迭代")
        result[slot] = {"attempt_iteration": attempt, "state": "omitted" if skipped else
                        "not_triggered" if attempt is None else "captured"}
    return result


def inspect(document):
    report = {"schema": 1, "diagnostic_only": True, "candidate_conclusion": "insufficient_evidence",
              "issues": [], "snapshots": [], "capture": None,
              "external_input_coordinates_used": False}
    try:
        document = obj(document, "probe")
        require(document.get("float_encoding") == ENCODING, "probe.float_encoding", "不支持的浮点编码")
        result = {}
        if "status" in document:
            require(isinstance(document["status"], str), "probe.status", "终态必须为字符串")
            result["status"] = document["status"]
        for key in ("nodes", "lp_iterations"):
            if key in document:
                result[key] = integer(document[key], f"probe.{key}")
        if "solver_error" in document:
            result["solver_error"] = boolean(document["solver_error"], "probe.solver_error")
        report["probe_result"] = result
        capture = obj(document.get("capture"), "capture")
        require(type(capture.get("schema")) is int and capture["schema"] == 1
                and capture.get("float_encoding") == ENCODING, "capture", "不支持的 schema 或浮点编码")
        require(capture.get("engine") in ("Dense", "Revised"), "capture.engine", "未知引擎")
        metadata = {name: capture.get(name) for name in ("engine", "terminal", "first_failure_detail", "terminal_detail")}
        require(isinstance(metadata["terminal"], str), "capture.terminal", "终态必须为字符串")
        for key in ("first_failure_detail", "terminal_detail"):
            require(key in capture and (metadata[key] is None or isinstance(metadata[key], str)), key, "故障原因须为字符串或 null")
        for name in ("source_id", "events_omitted", "objects_omitted", "omitted_estimated_bytes", "estimated_bytes", "byte_limit"):
            metadata[name] = integer(capture.get(name), name)
        for name in ("source_unavailable", "artificial_mapping_unavailable"):
            metadata[name] = boolean(capture.get(name), name)
            if metadata[name]:
                report["issues"].append(f"{name} 为 true")
        report["capture"] = metadata
        require(metadata["estimated_bytes"] <= metadata["byte_limit"], "capture", "采集估算量超过预算")
        slots = read_slots(capture)
        metadata["slots"] = slots
        if metadata["events_omitted"] or metadata["objects_omitted"] or any(s["state"] == "omitted" for s in slots.values()):
            report["issues"].append("采集达到预算，部分证据已省略")
        tol = obj(capture.get("tolerances"), "tolerances")
        tolerances = {key: finite(decode(tol.get(key), f"tolerances.{key}"), f"tolerances.{key}")
                      for key in TOLERANCES}
        require(all(value >= 0 for value in tolerances.values()), "tolerances", "容差不能为负")
        metadata["tolerances_bits"] = {key: tol[key] for key in TOLERANCES}
        raw_events = array(capture.get("events"), "events")
        require(len(raw_events) <= 4096, "events", "超过 schema 1 事件上限")
        events = [read_summary(event, f"events[{i}]") for i, event in enumerate(raw_events)]
        matrix = read_matrix(capture.get("normalized_matrix"))
        skipped_cuts = [i for i, row in enumerate(matrix.cuts) if row is None]
        metadata["skipped_cut_indices"] = skipped_cuts
        if skipped_cuts:
            report["issues"].append("未入表割的原始内容未采集")
        model = read_model(capture.get("engine_input_lp"), "engine_input_lp")
        input_scales, map_mismatches = check_input_mapping(model, matrix, tolerances["zero"])
        if map_mismatches:
            report["issues"].append("引擎输入和归一矩阵的来源映射不一致")
        metadata["mapping_mismatches"] = map_mismatches
        for event in events:
            column = event["pricing"]["column"]
            if column is not None:
                integer(column, "events.pricing.column", matrix.ncs)
        seen = set()
        for i, raw in enumerate(array(capture.get("snapshots"), "snapshots")):
            try:
                raw = obj(raw, f"snapshots[{i}]")
                kind = raw.get("kind")
                require(isinstance(kind, str), f"snapshots[{i}].kind", "快照种类必须为字符串")
                slot = "first_failure" if kind in FAILURE_KINDS else kind
                require(slot in slots and slot not in seen and slots[slot]["state"] == "captured",
                        f"snapshots[{i}]", "快照缺槽、重复或状态与槽记录不符")
                seen.add(slot)
                snapshot = inspect_snapshot(raw, matrix, model, input_scales, tolerances["feasibility"])
                require(snapshot["call_iterations"] == slots[slot]["attempt_iteration"], slot, "快照迭代与槽记录不符")
            except InvalidCapture as error:
                report["issues"].append(str(error))
                continue
            report["snapshots"].append(snapshot)
        for slot, state in slots.items():
            if state["state"] == "captured" and slot not in seen:
                report["issues"].append(f"{slot} 槽记录已采集，但找不到快照")
        for snapshot in report["snapshots"]:
            snapshot["conclusion"] = classify(snapshot, report["issues"])
            if snapshot["kind"] == "first_candidate":
                report["candidate_conclusion"] = snapshot["conclusion"]
                pricing_events = [event for event in events if event.get("kind") == "pricing_at_first_candidate"
                                  and event["state_revision"] == snapshot["state_revision"]
                                  and event["call_iterations"] == snapshot["call_iterations"]
                                  and event["pricing"]["evaluated_revision"] == snapshot["state_revision"]]
                report["same_state_candidate_pricing"] = [event["pricing"] for event in pricing_events]
                report["feasible_candidate_with_entering_column"] = (
                    snapshot["conclusion"] == "primal_and_augmented_pass"
                    and any(event["pricing"]["kind"] == "candidate" for event in pricing_events))
    except InvalidCapture as error:
        report["issues"].append(str(error))
    return report


def no_duplicate_keys(pairs):
    obj = {}
    for key, value in pairs:
        require(key not in obj, "JSON", f"键重复 {key}")
        obj[key] = value
    return obj


def no_json_constant(value):
    raise InvalidCapture(f"JSON: 非标准常量 {value}，浮点证据须用十六进制位串")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture", type=Path, help="phase1_progress_probe 生成的 JSON")
    parser.add_argument("--input", type=Path, help="可选原始输入文件，仅记录 SHA-256")
    parser.add_argument("--solver-commit", help="采集程序的求解器提交")
    parser.add_argument("--build-options", help="采集程序的构建选项")
    parser.add_argument("--probe-args", help="原样记录采集命令参数")
    args = parser.parse_args()
    try:
        raw = args.capture.read_bytes()
        provenance = {"capture_sha256": hashlib.sha256(raw).hexdigest(), "solver_commit": args.solver_commit,
                      "build_options": args.build_options, "probe_args": args.probe_args,
                      "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest() if args.input else None}
        try:
            document = json.loads(raw, object_pairs_hook=no_duplicate_keys, parse_constant=no_json_constant)
        except (InvalidCapture, ValueError) as error:
            report = {"schema": 1, "diagnostic_only": True, "candidate_conclusion": "insufficient_evidence",
                      "issues": [str(error)], "snapshots": []}
        else:
            report = inspect(document)
        report["provenance"] = provenance
        json.dump(report, sys.stdout, ensure_ascii=False, indent=2, allow_nan=False)
        print()
    except OSError as error:
        parser.exit(2, f"读取输入失败：{error}\n")


if __name__ == "__main__":
    main()
