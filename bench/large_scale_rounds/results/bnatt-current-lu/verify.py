#!/usr/bin/env python3
"""按原捕获的二进制浮点系数验证基矩阵秩，不求解原模型。"""

import argparse
from fractions import Fraction
import gzip
import hashlib
import heapq
import json
from pathlib import Path


ROW_WEIGHTS = {
    27: -2, 102: 1, 106: -1, 496: -4, 508: 2, 810: -2, 818: -2,
    1734: -4, 1740: -4, 1744: -4, 1990: 1, 1992: -1, 2228: 2,
    2998: 2, 3002: 2, 3676: 2, 3678: 2, 4068: 2, 4076: 2,
}
PRIME = 2**31 - 1


def require(condition, message):
    if not condition:
        raise ValueError(message)


def product(entries):
    return sum((Fraction(coefficient) * ROW_WEIGHTS[row]
                for row, coefficient in entries if row in ROW_WEIGHTS), Fraction(0))


def nonsingular_mod_prime(basis, columns):
    """分母均为二的幂，用有限域消元验证原有理行列式非零。"""
    rows = [{} for _ in basis]
    active_rows = [set() for _ in basis]
    modular_values = {}
    for j, column in enumerate(basis):
        for row, coefficient in columns[column]:
            if coefficient not in modular_values:
                num, den = coefficient.as_integer_ratio()
                modular_values[coefficient] = num * pow(den, -1, PRIME) % PRIME
            value = modular_values[coefficient]
            if value:
                rows[row][j] = value
                active_rows[j].add(row)

    ready = [(len(indices), j) for j, indices in enumerate(active_rows)]
    heapq.heapify(ready)
    finished = set()
    while ready:
        count, j = heapq.heappop(ready)
        target_rows = active_rows[j]
        if j in finished or len(target_rows) != count:
            continue
        if not target_rows:
            return False
        pivot_row = min(target_rows, key=lambda r: (len(rows[r]), r))
        pivot = rows[pivot_row]
        divisor = pow(pivot[j], -1, PRIME)
        target_rows.remove(pivot_row)
        finished.add(j)
        for other in pivot.keys() - {j}:
            active_rows[other].remove(pivot_row)
            heapq.heappush(ready, (len(active_rows[other]), other))
        for r in sorted(target_rows):
            row = rows[r]
            multiplier = row.pop(j) * divisor % PRIME
            for other in pivot.keys() - {j}:
                old = row.get(other, 0)
                value = (old - multiplier * pivot[other]) % PRIME
                if value:
                    row[other] = value
                    if not old:
                        active_rows[other].add(r)
                elif old:
                    del row[other]
                    active_rows[other].remove(r)
                if bool(old) != bool(value):
                    heapq.heappush(ready, (len(active_rows[other]), other))
    return len(finished) == len(basis)


def verify(path, check_prior):
    raw = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
    capture = json.loads(raw)
    failure = capture["failure"]
    require((failure["rows"], failure["phase"], failure["iterations"],
             failure["eta_count"], failure["failed_column"], failure["max_candidate"])
            == (5058, "phase1", 3904, 64, 4950, 0), "输入并非当前 3904 步现场")
    columns = {col["id"]: col["entries"] for col in failure["columns"]}
    require(len(columns) == len(failure["columns"]), "重复基列来源")
    require(all(len({r for r, _ in entries}) == len(entries)
                and all(0 <= r < failure["rows"] for r, _ in entries)
                for entries in columns.values()), "异常行索引")
    basis = failure["basis_last_refactor"].copy()
    require(len(basis) == failure["rows"] and len(set(basis)) == len(basis), "初始基无效")
    require(len(failure["pivots"]) == 64, "换基日志缺失")
    prior_count = 0
    leaving_product = None
    for i, pivot in enumerate(failure["pivots"]):
        row, leaving, entering = (pivot[key] for key in ("row", "leaving", "entering"))
        require(basis[row] == leaving and entering not in basis, f"第 {i + 1} 条记录不能重放")
        if check_prior and i <= 24:
            require(nonsingular_mod_prime(basis, columns), f"第 {i} 次换基后未证基非奇异")
            prior_count += 1
        if i == 24:
            require((pivot["iterations_before"], row, leaving, entering, pivot["step"])
                    == (3865, 1906, 14, 8, 0), "第 25 次换基与冻结转移不一致")
            leaving_product = product(columns[leaving])
            require(leaving_product == Fraction(7205759403792785, 18014398509481984),
                    "换出列点积不匹配")
        basis[row] = entering
        if i == 24:
            require(all(product(columns[j]) == 0 for j in basis), "换基后左零向量验证失败")
    require(basis == failure["basis_after"], "最终基未能重放")
    print(json.dumps({
        "capture_sha256": hashlib.sha256(raw).hexdigest(),
        "replayed_pivots": len(failure["pivots"]),
        "prior_bases_proven_nonsingular": prior_count if check_prior else None,
        "prime": PRIME if check_prior else None,
        "row_weights": ROW_WEIGHTS,
        "departing_column_exact_product": str(leaving_product),
        "exact_zero_products_after_transition": len(basis),
        "exact_entering_pivot": 0,
    }, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture", type=Path, nargs="?",
                        default=Path(__file__).with_name("capture.json.gz"))
    parser.add_argument("--check-prior", action="store_true", help="另用模素数消元验证此前 25 个基")
    args = parser.parse_args()
    verify(args.capture, args.check_prior)
