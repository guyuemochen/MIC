#!/usr/bin/env python3
"""生成带闭环的二层网络，枚举所有主变量，与 MIC 正式入口逐一对账。"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import random
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[4]


def generate(seed):
    rng = random.Random(21000 + seed)
    n = 4 + seed % 6
    aux = n
    rows = []
    functions = []
    fixed_zero = None

    def gate(output, inputs, conjunction):
        k = len(inputs)
        signed = [(j, 1 if positive else -1) for j, positive in inputs]
        neg = sum(a < 0 for _, a in signed)
        factor = k if conjunction else 1
        rows.append(([(output, factor)] + [(j, -a) for j, a in signed], '<=', neg))
        factor = 1 if conjunction else k
        rows.append(([(output, -factor)] + signed, '<=', (k - 1 if conjunction else 0) - neg))

    for y in range(n):
        deps = tuple(rng.sample(range(n), 3))
        outer = []
        patterns = []
        for bits in itertools.product((0, 1), repeat=3):
            slot = aux
            aux += 1
            if rng.randrange(2):
                gate(slot, list(zip(deps, bits)), True)
                outer.append((slot, True))
                patterns.append(bits)
            else:
                value = rng.randrange(2)
                rows.append(([(slot, 1)], '=', value))
                outer.append((slot, value == 0))
                if value == 0:
                    fixed_zero = slot
        gate(y, outer, False)
        functions.append((deps, set(patterns)))
    obj_col = fixed_zero if seed % 8 == 0 and fixed_zero is not None else rng.randrange(n)
    coeff = 0 if seed % 7 == 0 else rng.choice([-7, -2, 1, 5])
    sense = 'min' if seed % 2 == 0 else 'max'
    perm = list(range(aux))
    rng.shuffle(perm)
    rng.shuffle(rows)
    objective = [] if coeff == 0 else [f'{coeff} v{obj_col}']
    lines = ['name small-network', 'sense ' + sense, 'vars', *(f'v{j} binary' for j in perm),
             'obj', *objective, 'cons']
    for terms, op, rhs in rows:
        lines.append(' '.join(f'{a} v{j}' for j, a in terms) + f' {op} {rhs}')
    lines.append('end')
    costs = []
    for bits in itertools.product((0, 1), repeat=n):
        if all(bits[j] == (tuple(bits[k] for k in deps) in patterns)
               for j, (deps, patterns) in enumerate(functions)):
            costs.append(coeff * (bits[obj_col] if obj_col < n else 0))
    return '\n'.join(lines) + '\n', perm, rows, obj_col, coeff, sense, costs


def check_one(probe, path, cap, perm, rows, obj_col, coeff, sense, costs):
    command = [probe, 'timed-solve', path, '--threads', '1', '--time-limit-secs', '10',
               '--node-limit', str(cap), '--early-stop', 'off']
    result = json.loads(subprocess.check_output(command))
    optimum = (min if sense == 'min' else max)(costs) if costs else None
    assert result['nodes'] <= cap
    if cap > 1:
        assert result['status'] == ('infeasible' if not costs else 'optimal')
    for point in [result, *result['incumbent_events']]:
        if point['values'] is None:
            assert point['objective'] is None
            continue
        x = dict(zip(perm, point['values'], strict=True))
        assert all(v in (0, 1) for v in x.values())
        for terms, op, rhs in rows:
            act = sum(a * x[j] for j, a in terms)
            assert (act == rhs if op == '=' else act <= rhs)
        assert point['objective'] == coeff * x[obj_col] and point['objective'] in costs
    if result['status'] == 'optimal':
        assert result['objective'] == result['bound'] == optimum
    elif result['status'] == 'infeasible':
        assert not costs and result['objective'] is None and result['bound'] is None
    else:
        assert result['status'] == 'limit_reached'
        if optimum is not None and result['bound'] is not None:
            assert (result['bound'] <= optimum if sense == 'min' else result['bound'] >= optimum)
    return {k: result[k] for k in ('status', 'nodes', 'lp_iterations', 'objective', 'bound')} | {
        'node_limit': cap, 'callbacks': len(result['incumbent_events'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    probe = ROOT / 'milp/target/release/examples/round_probe'
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'test.milp'
        for seed in range(80):
            source, perm, rows, obj_col, coeff, sense, goals = generate(seed)
            path.write_text(source)
            runs = [check_one(probe, path, cap, perm, rows, obj_col, coeff, sense, goals)
                    for cap in (1, 2_000_000)]
            results.append(dict(seed=seed + 21000, sense=sense,
                                sha256=hashlib.sha256(source.encode()).hexdigest(),
                                feasible_assignments=len(goals), runs=runs))
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('已核对', len(results), '个网络，', len(results) * 2, '次运行')


if __name__ == '__main__':
    main()
