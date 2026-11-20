"""从已独立核验的两例运输解中还原路线与原模型时刻。"""
import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path


def decode(model, values):
    variables = model['variables']
    names = [v['name'] for v in variables]
    clock = {j for j, v in enumerate(variables) if v['kind'] == 'continuous'}
    lower = {j: variables[j]['lb'] for j in clock}
    upper = {j: variables[j]['ub'] for j in clock}
    skipped = {}
    links = []
    memberships = defaultdict(set)
    for i, row in enumerate(model['rows']):
        terms = row['terms']
        if row['op'] == '=':
            for j, _ in terms:
                memberships[j].add(i)
            continue
        sign = 1 if row['op'] == '>=' else -1
        if len(terms) == 1 and terms[0][0] in clock:
            j, v = terms[0]
            limit = row['rhs'] / v
            if v * sign > 0:
                lower[j] = max(lower[j], limit)
            else:
                upper[j] = min(upper[j], limit) if upper[j] is not None else limit
        if len(terms) != 4 or sum(j in clock for j, _ in terms) != 2:
            continue
        earlier, = (j for j, c in terms if j in clock and sign * c == -1)
        later, = (j for j, c in terms if j in clock and sign * c == 1)
        diagonal, = (j for j, c in terms if j not in clock and sign * c > 0)
        arc, = (j for j, c in terms if j not in clock and sign * c < 0)
        assert skipped.setdefault(earlier, diagonal) == diagonal
        links.append((earlier, later, arc))
    assert set(skipped) == clock
    chosen = {j for j, v in enumerate(values) if j not in clock and v > 0.5}
    selected_in_group = {}
    for j in chosen:
        for group in memberships[j]:
            assert group not in selected_in_group
            selected_in_group[group] = j
    outgoing, incoming = {}, {}
    successors, prev = {}, set()
    for earlier, later, arc in links:
        src, = memberships[skipped[earlier]] & memberships[arc]
        dst, = memberships[skipped[later]] & memberships[arc]
        assert outgoing.setdefault(earlier, src) == src
        assert incoming.setdefault(later, dst) == dst
        if arc in chosen:
            assert later not in prev and earlier not in successors
            successors[earlier] = later
            prev.add(later)
    served = {j for j in clock if skipped[j] not in chosen}
    routes, seen = [], set()
    for j in sorted(served - prev):
        stops = []
        while True:
            assert j not in seen and j in served
            seen.add(j)
            assert lower[j] <= values[j] <= upper[j]
            stops.append(dict(clock=names[j], window=[lower[j], upper[j]], time=values[j],
                              inbound=names[selected_in_group[incoming[j]]],
                              outbound=names[selected_in_group[outgoing[j]]]))
            if j not in successors:
                break
            j = successors[j]
        routes.append(stops)
    assert seen == served
    return dict(routes=routes, tasks_served=len(seen), tasks_skipped=[names[j] for j in sorted(clock - served)])


def check_exact_point(model, values):
    """用原始双精度系数表示的精确分数核对 MIC 的首解和最终解。"""
    x = list(map(Fraction, values))
    assert len(x) == len(model['variables'])
    for value, var in zip(x, model['variables'], strict=True):
        assert value >= var['lb'] and (var['ub'] is None or value <= var['ub'])
        assert var['kind'] == 'continuous' or value.denominator == 1
    for row in model['rows']:
        lhs = sum(Fraction(a) * x[j] for j, a in row['terms'] if x[j])
        rhs = Fraction(row['rhs'])
        if row['op'] == '<=':
            assert lhs <= rhs
        elif row['op'] == '>=':
            assert lhs >= rhs
        else:
            assert row['op'] == '=' and lhs == rhs


if __name__ == '__main__':
    path = Path(sys.argv[1])
    solver = sys.argv[2]
    results = {}
    for row in map(json.loads, path.open()):
        if row['rep'] == 0 and row['solver'] == solver and row['status'] == 'optimal':
            model = json.loads((path.parent / row['model']).read_bytes())
            if solver == 'mic':
                check_exact_point(model, row['incumbent_events'][0]['values'])
                check_exact_point(model, row['values'])
            results[row['instance']] = dict(objective=row['objective'], model_sha256=row['model_sha256'],
                                            **decode(model, row['values']))
    assert len(results) == 2
    print(json.dumps(results, ensure_ascii=False, indent=2, allow_nan=False))
