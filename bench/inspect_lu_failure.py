#!/usr/bin/env python3
"""只读复核LU捕获的换基日志与残差，不重解LP，不推断模型可行性。"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def norm(values):
    if any(v is None or not math.isfinite(v) for v in values):
        return None
    return max(map(abs, values), default=0.0)


def residual(basis, columns, alpha, rhs, rows):
    if len(alpha) != len(basis) or norm(alpha) is None:
        return None
    product, row_sums = [0.0] * rows, [0.0] * rows
    for column, value in zip(basis, alpha):
        for row, coefficient in columns[column]:
            product[row] += coefficient * value
            row_sums[row] += abs(coefficient)
    error = norm([a - b for a, b in zip(product, rhs)])
    a_norm, x_norm, rhs_norm = norm(row_sums), norm(alpha), norm(rhs)
    if any(value is None for value in (a_norm, x_norm, rhs_norm)):
        return None
    denominator = a_norm * x_norm + rhs_norm
    if error is None or not math.isfinite(denominator):
        return None
    return error / denominator if denominator else error


def inspect(document):
    failure = document.get("failure")
    if failure is None:
        return {"captured": False}
    rows = failure["rows"]
    columns = {column["id"]: column["entries"] for column in failure["columns"]}
    basis = list(failure["basis_last_refactor"])
    replay_ok = len(basis) == rows
    for pivot in failure["pivots"]:
        row = pivot["row"]
        if not replay_ok or basis[row] != pivot["leaving"]:
            replay_ok = False
            break
        basis[row] = pivot["entering"]
    replay_ok = replay_ok and basis == failure["basis_after"]
    failed = failure["failed_column"]
    summary = {
        "captured": True, "phase": failure["phase"], "reason": failure["reason"],
        "rows": rows, "iterations": failure["iterations"], "eta_count": failure["eta_count"],
        "failed_factor_column": failed, "max_candidate": failure["max_candidate"],
        "basis_column_at_failure": (failure["basis_after"][failure["factor_order"][failed]]
                                    if failed is not None else None),
        "basis_history_replays": replay_ok,
        "fresh_checks": [{key: value for key, value in c.items() if key != "alpha"}
                         for c in document.get("fresh_checks", [])],
    }
    history = document.get("history_checks", [])
    failed_history = [c["after_pivots"] for c in history if c["check"]["error"] is not None]
    summary["history_factorizations"] = len(history)
    summary["history_failed_after_pivots"] = failed_history
    if failed_history:
        first = min(failed_history)
        summary["first_history_failure_after_pivots"] = first
        if first > 0:
            summary["transition_into_first_failure"] = failure["pivots"][first - 1]
            previous = next((c["check"] for c in history if c["after_pivots"] == first - 1), None)
            summary["fresh_check_before_first_failure"] = previous
    if failure["pivots"]:
        last = failure["pivots"][-1]
        before = list(failure["basis_after"])
        assert before[last["row"]] == last["entering"]
        before[last["row"]] = last["leaving"]
        rhs = [0.0] * rows
        for row, value in columns[last["entering"]]:
            rhs[row] += value
        original = failure["last_pivot_alpha"]
        summary["last_pivot"] = last
        summary["eta_alpha_relative_residual"] = (
            residual(before, columns, original, rhs, rows) if original is not None else None)
        fresh = next((c for c in document.get("fresh_checks", [])
                      if c["name"] == "before_last_pivot"), None)
        if fresh and fresh["alpha"] is not None and original is not None:
            summary["alpha_max_absolute_difference"] = norm(
                [a - b for a, b in zip(original, fresh["alpha"])])
    return summary


def scipy_checks(document):
    """显式可选：只在benchmark环境使用独立SuperLU，不影响Rust求解器。"""
    import numpy as np
    import scipy
    from scipy.sparse import csc_matrix
    from scipy.sparse.linalg import splu

    failure = document["failure"]
    history = document.get("history_checks", [])
    bad = [c["after_pivots"] for c in history if c["check"]["error"] is not None]
    if not bad or min(bad) == 0:
        return {"scipy_version": scipy.__version__, "reason": "缺少可复核的首次失败转移"}
    first = min(bad)
    rows = failure["rows"]
    columns = {c["id"]: c["entries"] for c in failure["columns"]}
    pivot = failure["pivots"][first - 1]
    results = []
    for after in (first - 1, first):
        basis = list(failure["basis_last_refactor"])
        for event in failure["pivots"][:after]:
            assert basis[event["row"]] == event["leaving"]
            basis[event["row"]] = event["entering"]
        values, indices, starts = [], [], [0]
        for column in basis:
            for row, value in columns[column]:
                values.append(value)
                indices.append(row)
            starts.append(len(values))
        matrix = csc_matrix((np.array(values), np.array(indices), np.array(starts)), shape=(rows, rows))
        for order in ("NATURAL", "COLAMD"):
            result = {"after_pivots": after, "column_order": order, "diag_pivot_thresh": 1.0}
            try:
                lu = splu(matrix, permc_spec=order, diag_pivot_thresh=1.0)
                diagonal = np.abs(lu.U.diagonal())
                result["min_abs_u_diagonal"] = float(diagonal.min())
                result["max_abs_u"] = float(np.abs(lu.U.data).max())
                result["max_abs_basis"] = float(np.abs(matrix.data).max())
                if after == first - 1:
                    rhs = np.zeros(rows)
                    for row, value in columns[pivot["entering"]]:
                        rhs[row] += value
                    alpha = lu.solve(rhs)
                    result["alpha_r"] = float(alpha[pivot["row"]])
                    result["relative_residual"] = residual(basis, columns, alpha.tolist(), rhs.tolist(), rows)
            except RuntimeError as error:
                result["error"] = str(error)
            results.append(result)
    return {"scipy_version": scipy.__version__, "results": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture", type=Path)
    parser.add_argument("--scipy", action="store_true", help="用临时benchmark环境中的SciPy独立复核首次失败转移")
    args = parser.parse_args()
    raw = args.capture.read_bytes()
    document = json.loads(raw)
    result = inspect(document)
    if args.scipy:
        result["independent_superlu"] = scipy_checks(document)
    result["source_sha256"] = hashlib.sha256(raw).hexdigest()
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
