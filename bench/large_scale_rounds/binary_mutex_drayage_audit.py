#!/usr/bin/env python3
"""在旧版缩减模型上独立核查两例 drayage 的短行互斥活动度。"""

import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import math
from pathlib import Path
from time import process_time


def load(path):
    rows, bounds, kinds = [], {}, {}
    section = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line in ('vars', 'obj', 'cons', 'bounds', 'end'):
            section = line
        elif section == 'vars':
            name, kind = line.split()
            kinds[name] = kind
            bounds[name] = (0.0, 1.0 if kind == 'binary' else math.inf)
        elif section == 'bounds':
            name, op, value = line.split()
            lo, hi = bounds[name]
            bounds[name] = (float(value), hi) if op == '>=' else (lo, float(value))
        elif section == 'cons':
            name, expression = line.split(':', 1)
            *lhs, op, rhs = expression.split()
            terms = [(lhs[i + 1], float(lhs[i])) for i in range(0, len(lhs), 2)]
            rows.append((name, terms, op, float(rhs)))
    return rows, bounds, kinds


def eligible_arcs(rows, bounds, kinds):
    """在报告的浮点模型上用精确有理活动度独立计算可固定弧。"""
    groups = defaultdict(set)
    count = 0
    for name, terms, op, rhs in rows:
        if op == '=' and rhs > 0 and all(a == rhs and kinds[j] == 'binary' for j, a in terms):
            count += 1
            for j, _ in terms:
                groups[j].add(name)
    assert count == 210 and all(len(g) == 2 for g in groups.values())

    independent, with_mutex = {}, {}
    checked = 0
    for name, terms, op, rhs in rows:
        binary = [(j, a) for j, a in terms if kinds[j] == 'binary']
        if len(terms) != 4 or len(binary) != 2:
            continue
        assert op == '>=', (name, op)
        checked += 1
        continuous = [(j, a) for j, a in terms if kinds[j] == 'continuous']
        maximum_bounds = [bounds[j][a > 0] for j, a in continuous]
        if not all(map(math.isfinite, maximum_bounds)):
            continue
        maximum = sum(Fraction(a) * Fraction(bound)
                      for (_, a), bound in zip(continuous, maximum_bounds))
        residual = Fraction(rhs) - maximum
        (x, a), (y, b) = binary
        assert bounds[x] == bounds[y] == (0.0, 1.0)
        assert groups[x] & groups[y], (name, x, y)
        # 本行两个二元都取零时已有一组连续变量可满足本行。假设本行
        # 之外的变量为 1 最多把其中某些二元置零，无法单独制造冲突。
        assert residual <= 0, (name, residual)
        for candidate, weight, other_weight in [(x, a, b), (y, b, a)]:
            if residual > Fraction(weight) + max(Fraction(other_weight), 0):
                independent[candidate] = name, residual - Fraction(weight) - max(Fraction(other_weight), 0)
            if residual > Fraction(weight):
                with_mutex[candidate] = name, residual - Fraction(weight)
    assert checked == 4160
    assert independent.keys() <= with_mutex.keys()
    return independent, with_mutex


def r6_single_row_fixed(rows, original, kinds):
    """模拟历史逐行盒传播，在旧版缩减模型上达到不动点。"""
    bounds = original.copy()
    inequalities = []
    for _, terms, op, rhs in rows:
        if op in ('>=', '='):
            inequalities.append((terms, rhs))
        if op in ('<=', '='):
            inequalities.append(([(j, -a) for j, a in terms], -rhs))
    guard = 1e-6
    for _ in range(20):
        changed = False
        for terms, rhs in inequalities:
            values = [a * bounds[j][a > 0] for j, a in terms]
            assert all(map(math.isfinite, values))
            total = math.fsum(values)
            assert total >= rhs - guard
            for (j, a), previous in zip(terms, values):
                implied = (rhs - guard - (total - previous)) / a
                lo, hi = bounds[j]
                if a > 0:
                    lo = max(lo, math.ceil(implied) if kinds[j] == 'binary' else implied)
                else:
                    hi = min(hi, math.floor(implied) if kinds[j] == 'binary' else implied)
                assert lo <= hi
                best = lo, hi
                old_lo, old_hi = bounds[j]
                if lo > old_lo + guard or hi < old_hi - guard:
                    changed = True
                    bounds[j] = best
        if not changed:
            return {j for j, (lo, hi) in bounds.items() if kinds[j] == 'binary' and lo == hi == 0}
    raise ValueError('两例的逐行界传播未在预期轮次内收敛')


def digest(names):
    return hashlib.sha256('\n'.join(sorted(names)).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('raw', type=Path)
    parser.add_argument('reduced_base', type=Path)
    parser.add_argument('base_fixes', type=Path)
    parser.add_argument('candidate_fixes', type=Path)
    args = parser.parse_args()

    raw = load(args.raw)
    reduced = load(args.reduced_base)
    assert eligible_arcs(*raw) == ({}, {})
    assert not args.base_fixes.read_text().strip()
    audit_start = process_time()
    independent, with_mutex = eligible_arcs(*reduced)
    elapsed_audit = process_time() - audit_start
    r6_start = process_time()
    r6_fixed = r6_single_row_fixed(*reduced)
    elapsed_r6 = process_time() - r6_start
    assert r6_fixed == independent.keys()
    actual = set()
    original_columns = list(raw[2])
    for line in args.candidate_fixes.read_text().splitlines():
        name, column, value = line.split('\t')
        assert original_columns[int(column)] == name
        assert float(value) == 0 and name not in actual
        actual.add(name)
    assert with_mutex.keys() == actual, (len(with_mutex), len(actual))
    extra = actual - r6_fixed
    print(f'{args.raw.stem} 旧版固定 0，逐行推断 {len(r6_fixed)}，互斥固定 {len(actual)}，增量 {len(extra)}')
    print(f'固定集 SHA-256 {digest(actual)}，增量集 SHA-256 {digest(extra)}')
    print(f'最小增量违反量 {float(min(gap for j, (_, gap) in with_mutex.items() if j in extra)):.12g}')
    print(f'独立有理数核查 {elapsed_audit * 1000:.3f} ms，逐行传播模拟 {elapsed_r6 * 1000:.3f} ms')


if __name__ == '__main__':
    main()
