#!/usr/bin/env python3
"""核查首次根 LP 报告点与逐列取整点在原始 drayage 约束下的活动度。"""

import argparse
from collections import Counter
from fractions import Fraction
import json
import math
from pathlib import Path

from binary_mutex_drayage_audit import load


def row_type(terms, op, kinds):
    binary = sum(kinds[j] == 'binary' for j, _ in terms)
    if op == '=' and binary == len(terms):
        return 'exactly_one'
    if len(terms) == 4 and binary == 2:
        return 'time_order'
    if len(terms) == 1:
        return 'single_bound'
    return 'other'


def examine(rows, bounds, kinds, point):
    by_type = Counter()
    violations = []
    maximum_bound = 0.0
    for j, value in point.items():
        low, high = bounds[j]
        maximum_bound = max(maximum_bound, low - value, value - high)
    for name, terms, op, rhs in rows:
        terms_at_point = [a * point[j] for j, a in terms]
        activity = math.fsum(terms_at_point)
        margin = 1e-6 + 1e-9 * max(abs(rhs), math.fsum(map(abs, terms_at_point)))
        difference = {'>=': rhs - activity, '<=': activity - rhs, '=': abs(activity - rhs)}[op]
        if difference <= margin:
            continue
        family = row_type(terms, op, kinds)
        by_type[family] += 1
        contribution = [
            dict(column=j, coefficient=a, value=point[j])
            for j, a in terms if point[j] != 0
        ]
        violations.append((difference, name, family, activity, rhs, contribution))
    violations.sort(key=lambda item: (-item[0], item[1]))
    return dict(violating_rows=len(violations), by_type=by_type, max_bound_violation=maximum_bound,
                worst=[dict(violation=v, row=name, family=family, activity=lhs, rhs=rhs, nonzeros=terms)
                       for v, name, family, lhs, rhs, terms in violations[:8]])


def retime(rows, bounds, kinds, rounded):
    """固定已取整的弧，仅检查原模型剩余时间变量的精确差分约束。"""
    origin = '时间原点'
    nodes = [origin] + [j for j, kind in kinds.items() if kind == 'continuous']
    edges = []
    for j in nodes[1:]:
        lo, hi = bounds[j]
        if math.isfinite(lo):
            edges.append((origin, j, Fraction(lo), f'{j} 下界'))
        if math.isfinite(hi):
            edges.append((j, origin, -Fraction(hi), f'{j} 上界'))
    for name, terms, op, rhs in rows:
        residual = Fraction(rhs) - sum(
            Fraction(a) for j, a in terms if kinds[j] == 'binary' and rounded[j] == 1.0)
        times = [(j, Fraction(a)) for j, a in terms if kinds[j] == 'continuous']
        signs = [1] if op == '>=' else [-1] if op == '<=' else [1, -1]
        for sign in signs:
            if not times:
                if sign * residual > 0:
                    return dict(status='infeasible_binary_row', row=name, violation=str(sign * residual))
            elif len(times) == 1:
                j, a = times[0]
                bound = residual / a
                edge = (origin, j, bound) if sign * a > 0 else (j, origin, -bound)
                edges.append((*edge, name))
            else:
                assert len(times) == 2 and times[0][1] == -times[1][1], (name, times)
                pos, neg = sorted(times, key=lambda term: -sign * term[1])
                edges.append((neg[0], pos[0], residual / pos[1], name))
    distance = dict.fromkeys(nodes, Fraction(0))
    predecessor = {}
    for _ in nodes:
        updated = None
        for src, dst, delta, name in edges:
            proposed = distance[src] + delta
            if proposed > distance[dst]:
                distance[dst] = proposed
                predecessor[dst] = (src, dst, delta, name)
                updated = dst
        if updated is None:
            rescheduled = rounded | {j: float(distance[j] - distance[origin]) for j in nodes[1:]}
            checked = examine(rows, bounds, kinds, rescheduled)
            assert checked['violating_rows'] == 0 and checked['max_bound_violation'] <= 1e-6
            return dict(status='feasible', checked=checked)
    for _ in nodes:
        updated = predecessor[updated][0]
    first = updated
    cycle = []
    while True:
        edge = predecessor[updated]
        cycle.append(edge)
        updated = edge[0]
        if updated == first:
            break
    total = sum(edge[2] for edge in cycle)
    assert total > 0
    return dict(status='infeasible_positive_cycle', total_gap=float(total), total_gap_exact=str(total),
                cycle=[dict(source=src, target=dst, gap=float(gap), gap_exact=str(gap), row=name)
                       for src, dst, gap, name in reversed(cycle)])


def analyze(raw, values):
    rows, bounds, kinds = load(raw)
    assert set(kinds.values()) <= {'continuous', 'binary'}
    columns = list(kinds)
    point = {}
    for line in values.read_text().splitlines():
        j, name, value = line.split('\t')
        assert name not in point and columns[int(j)] == name and math.isfinite(float(value))
        point[name] = float(value)
    assert point.keys() == kinds.keys()
    rounded = point.copy()
    deviations = []
    for j, value in point.items():
        if kinds[j] == 'binary':
            rounded[j] = float(value >= 0.5)
            deviations.append((abs(value - rounded[j]), j, value))
    deviations.sort(key=lambda item: (-item[0], item[1]))
    fractional = [(delta, j, value) for delta, j, value in deviations if delta > 1e-6]
    fraction_names = {j for _, j, _ in fractional}
    assignments = [
        [j for j, _ in terms]
        for _, terms, op, rhs in rows
        if op == '=' and rhs == 1 and all(kinds[j] == 'binary' and a == 1 for j, a in terms)
    ]
    assert len(assignments) == 210
    groups_per_fraction_count = Counter(sum(j in fraction_names for j in group) for group in assignments)
    rounded_selected_per_group = Counter(sum(rounded[j] == 1.0 for j in group) for group in assignments)
    raw_point = examine(rows, bounds, kinds, point)
    rounded_point = examine(rows, bounds, kinds, rounded)
    return dict(binary_columns=len(deviations), fractional=len(fractional), bins={
                    '1e-6_to_1e-3': sum(delta <= 1e-3 for delta, _, _ in fractional),
                    '1e-3_to_0.1': sum(1e-3 < delta <= 0.1 for delta, _, _ in fractional),
                    'above_0.1': sum(delta > 0.1 for delta, _, _ in fractional),
                },
                max_integrality_deviation=deviations[0][0],
                largest_fractional=[dict(deviation=delta, column=j, value=value)
                                    for delta, j, value in fractional[:12]],
                groups_per_fraction_count=dict(sorted(groups_per_fraction_count.items())),
                rounded_selected_per_group=dict(sorted(rounded_selected_per_group.items())),
                original_lp_point=raw_point, rounded_point=rounded_point,
                retimed_rounded_arcs=retime(rows, bounds, kinds, rounded))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('raw', type=Path)
    parser.add_argument('values', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    evidence = analyze(args.raw, args.values)
    with args.output.open('x', encoding='utf-8') as output:
        json.dump(evidence, output, ensure_ascii=False, indent=2)
        output.write('\n')


if __name__ == '__main__':
    main()
