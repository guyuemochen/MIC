#!/usr/bin/env python3
"""从原始 MILP 文本逐项验解，并与等价模型、完全枚举及 Gurobi 独立对账。"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess
import sys
from tempfile import TemporaryDirectory

ROUNDS = Path(__file__).resolve().parents[3]
ROOT = ROUNDS.parents[1]
sys.path.insert(0, str(ROUNDS))
from summarize_all_timing import parse_milp

LIMIT = 2_000_000
SOURCE = ROOT / 'bench/instances/setcover-400x800.milp'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def milp(costs, rows, flipped=()):
    flipped = set(flipped)
    lines = ['name audit', 'sense min', 'vars']
    lines.extend(f'z{j} binary' for j in range(len(costs)))
    lines.extend(['obj', ' '.join(f'{c} z{j}' for j, c in enumerate(costs)), 'cons'])
    for i, row in enumerate(rows):
        unit, op = (-1, '<=') if i in flipped else (1, '>=')
        terms = ' '.join(f'{unit} z{j}' for j in row)
        lines.append(f'{terms} {op} {unit}')
    return ('\n'.join([*lines, 'end']) + '\n').encode()


def witness(model, values, expected):
    assert values is not None and len(values) == len(model['variables'])
    assert all(math.isfinite(x) and abs(x - round(x)) <= 1e-6 and round(x) in (0, 1) for x in values)
    x = [round(v) for v in values]
    for i, row in enumerate(model['rows']):
        lhs = sum(int(a) * x[j] for j, a in row['terms'])
        assert (lhs >= row['rhs'] if row['op'] == '>=' else lhs <= row['rhs']), (i, lhs, row)
    value = sum(int(c) * x[j] for j, c in model['objective'])
    assert abs(expected - value) <= 1e-6, (expected, value)
    return [j for j, bit in enumerate(x) if bit]


def check(row, model, optimum, limit):
    assert row['status'] in ('optimal', 'limit_reached') and not row.get('callback_errors')
    assert 0 <= row['nodes'] <= limit
    assert row['bound'] is not None and math.isfinite(row['bound']) and row['bound'] <= optimum + 1e-6
    if row['status'] == 'optimal':
        assert abs(row['objective'] - optimum) <= 1e-6 and abs(row['bound'] - optimum) <= 1e-6
    else:
        assert row['objective'] >= optimum - 1e-6
    selected = witness(model, row['values'], row['objective'])
    events = row['incumbent_events']
    for event in events:
        witness(model, event['values'], event['objective'])
        assert event['objective'] >= optimum - 1e-6
    return dict(status=row['status'], objective=row['objective'], bound=row['bound'], nodes=row['nodes'],
                selected=selected, verified_incumbents=[event['objective'] for event in events])


def mic(binary, path, limit):
    command = [binary, 'timed-solve', path, '--threads', '1', '--time-limit-secs', '60',
               '--node-limit', str(limit), '--early-stop', 'off']
    row = json.loads(subprocess.check_output(command, timeout=90))
    assert row['lp_iterations'] == 0 and not row['auto_degraded'] and not row['gap_terminated']
    return row


def gurobi(library, temp, model):
    path = temp / 'model.json'
    path.write_text(json.dumps(model))
    result = subprocess.check_output([sys.executable, ROUNDS / 'compare_all_timing.py',
                                      '_worker', library, path, '60'], text=True, timeout=90)
    row = json.loads(result.strip().splitlines()[-1])
    assert row['status'] == 'optimal' and row['version'] == '10.0.3'
    return row


def exhaustive(seed):
    """用十二列作真值。新增列只覆盖某个原列覆盖行的子集，且费用不低于原列。"""
    rng = random.Random(seed)
    k, m, n = 12, 32, 256
    masks = [rng.getrandbits(m) & rng.getrandbits(m) for _ in range(k)]
    for i in range(m):
        masks[rng.randrange(k)] |= 1 << i
    costs = [rng.randint(1, 10) * (90_000 if seed % 3 == 0 else 1) for _ in range(k)]
    full = (1 << m) - 1
    cheapest = sum(costs)
    feasible = 0
    for assignment in range(1 << k):
        covered = value = 0
        for j in range(k):
            if (assignment >> j) & 1:
                covered |= masks[j]
                value += costs[j]
        if covered == full:
            feasible += 1
            cheapest = min(cheapest, value)
    for j in range(k, n):
        original = rng.randrange(k)
        mask, cost = masks[original], costs[original]
        if j % 3 == 0:
            masks.append(mask)
            costs.append(cost)
        else:
            masks.append(mask & rng.getrandbits(m))
            costs.append(1_000_000 if j % 3 == 1 else cost + rng.randint(0, 99))
    rows = [[j for j, mask in enumerate(masks) if mask & (1 << i)] for i in range(m)]
    return milp(costs, rows, range(0, m, 3)), cheapest, feasible


def random_cover(seed):
    rng = random.Random(seed)
    rows = [[j for j in range(800) if rng.random() < 0.08] for _ in range(400)]
    for row in rows:
        if not row:
            row.append(rng.randrange(800))
    costs = [rng.randint(1, 100) for _ in range(800)]
    return milp(costs, rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary', type=Path)
    parser.add_argument('--gurobi', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve()
    report = dict(binary_sha256=digest(binary.read_bytes()), gurobi_version='10.0.3',
                  gurobi_library_sha256=digest(args.gurobi.read_bytes()), cases={})
    with TemporaryDirectory() as tmp:
        temp = Path(tmp)

        def run(name, payload, expected=None, gurobi_check=False, limits=(), boundary=False):
            model = parse_milp(payload)
            assert model['sense'] == 'min' and model['offset'] == 0
            assert all(v['kind'] == 'binary' and v['lb'] == 0 and v['ub'] == 1 for v in model['variables'])
            assert all(c.is_integer() and c > 0 for _, c in model['objective'])
            for row in model['rows']:
                assert (row['op'], row['rhs']) in (('>=', 1), ('<=', -1))
                assert all(a == row['rhs'] for _, a in row['terms'])
            path = temp / (name + '.milp')
            path.write_bytes(payload)
            entry = dict(milp_sha256=digest(payload), variables=len(model['variables']), rows=len(model['rows']))
            if gurobi_check:
                grb = gurobi(args.gurobi, temp, model)
                if expected is None:
                    expected = round(grb['objective'])
                entry['gurobi'] = check(grb, model, expected, LIMIT)
            assert expected is not None
            entry['optimum'] = expected
            result = mic(binary, path, LIMIT)
            entry['mic'] = check(result, model, expected, LIMIT)
            assert result['status'] == 'optimal'
            if boundary:
                assert result['nodes'] > 1
                limits = sorted(set([*limits, result['nodes'] - 1, result['nodes']]))
            if limits:
                entry['limited'] = {}
                for cap in limits:
                    row = mic(binary, path, cap)
                    summary = check(row, model, expected, cap)
                    if cap < result['nodes']:
                        assert row['status'] == 'limit_reached' and row['nodes'] == cap
                    else:
                        assert row['status'] == 'optimal'
                    entry['limited'][cap] = summary
            report['cases'][name] = entry
            print(name, expected, result['nodes'], flush=True)
            return entry

        payload = SOURCE.read_bytes()
        native = parse_milp(payload)
        original = run('original', payload, gurobi_check=True, limits=(1,), boundary=True)
        opt = original['optimum']
        costs = [int(c) for _, c in native['objective']]
        rows = [[j for j, _ in row['terms']] for row in native['rows']]
        rng = random.Random(21837)
        cols = list(range(len(costs)))
        rng.shuffle(cols)
        inverse = {old: new for new, old in enumerate(cols)}
        shuffled = [[inverse[j] for j in row] for row in rows]
        rng.shuffle(shuffled)
        permuted = run('permuted', milp([costs[j] for j in cols], shuffled, range(1, len(rows), 2)), opt)
        selected = set(cols[j] for j in permuted['mic']['selected'])
        witness(native, [int(j in selected) for j in range(len(costs))], opt)
        forced = run('scaled-and-forced', milp([7*c for c in costs] + [19], rows + [[len(costs)]]),
                     opt * 7 + 19)
        selected = set(forced['mic']['selected'])
        assert len(costs) in selected
        witness(native, [int(j in selected) for j in range(len(costs))], opt)
        for seed in range(4101, 4113):
            payload, expected, feasible = exhaustive(seed)
            entry = run(f'enumerated-{seed}', payload, expected, limits=(1,))
            entry.update(base_columns=12, assignments=4096, feasible_assignments=feasible)
        for seed in (7341, 7342):
            run(f'random-{seed}', random_cover(seed), gurobi_check=True)
    print(Counter(e['mic']['status'] for e in report['cases'].values()))
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')


if __name__ == '__main__':
    main()
