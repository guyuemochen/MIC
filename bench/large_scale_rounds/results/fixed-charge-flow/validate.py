#!/usr/bin/env python3
"""对冻结原例、等价变形和可独立枚举的小规模配送模型复核可行解及界。"""
import argparse
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROUNDS = HERE.parents[1]
FROZEN = ROUNDS / 'results/main-gurobi10-100-60s'
sys.path.insert(0, str(ROUNDS))
from check_results import verify_solution
from summarize_all_timing import audit_truth, load_case, parse_milp, public_truths, read_bytes, read_events


def check_native(mps, model):
    """独立读取该原例使用的 MPS 段，并核对两个求解器读入的全部系数、界和类型。"""
    kinds, columns, costs, upper = {}, {}, {}, {}
    by_row, rhs, integer = defaultdict(list), {}, set()
    section, markers, sections = None, False, []
    for line in gzip.decompress(mps).decode('ascii').splitlines():
        if not line or line.startswith('*'):
            continue
        fields = line.split()
        if not line[0].isspace():
            section = fields[0]
            sections.append(section)
        elif section == 'ROWS':
            kind, row = fields
            assert row not in kinds
            kinds[row] = kind
        elif section == 'COLUMNS':
            name = fields[0]
            if fields[1] == "'MARKER'":
                assert fields[2] == ("'INTEND'" if markers else "'INTORG'")
                markers = not markers
                continue
            j = columns.setdefault(name, len(columns))
            if markers:
                integer.add(j)
            for r, number in zip(fields[1::2], fields[2::2], strict=True):
                if kinds[r] == 'N':
                    costs[j] = int(number)
                else:
                    by_row[r].append([j, int(number)])
        elif section == 'RHS':
            for r, value in zip(fields[1::2], fields[2::2], strict=True):
                rhs[r] = int(value)
        elif section == 'BOUNDS':
            kind, _, name, value = fields
            assert kind == 'UP'
            upper[columns[name]] = int(value)
        else:
            raise ValueError(section)
    assert sections == ['NAME','ROWS','COLUMNS','RHS','BOUNDS','ENDATA'] and not markers
    objective, = [r for r, k in kinds.items() if k == 'N']
    assert rhs.get(objective, 0) == model['offset'] == 0 and model['sense'] == 'min'
    assert model['objective'] == [[j, v] for j, v in sorted(costs.items()) if v]
    assert model['variables'] == [
        dict(name=f'x{j}', lb=0, ub=upper.get(j), kind='binary' if j in integer else 'continuous')
        for j in range(len(columns))]
    assert all(upper.get(j) == 1 for j in integer)
    ops = {'L': '<=', 'G': '>=', 'E': '='}
    assert model['rows'] == [dict(terms=sorted(by_row[r]), op=ops[k], rhs=rhs.get(r, 0))
                             for r, k in kinds.items() if k != 'N']


def solve(binary, path, sec, nodes):
    return json.loads(subprocess.check_output([
        binary, 'timed-solve', path, '--time-limit-secs', str(sec), '--node-limit', str(nodes),
        '--threads', '1', '--early-stop', 'off',
    ], timeout=sec + 45))


def check(model, row, optimum):
    errors = []
    objective = verify_solution(model, row, errors.append)
    first, goals, rejected = read_events(model, row, errors.append)
    witnesses = goals + ([objective] if objective is not None else [])
    audit_truth(model, dict(status='optimal', objective=optimum), row, witnesses, errors.append)
    assert not (rejected or errors or row.get('callback_errors')), (errors, rejected)
    summary = {k: row[k] for k in ('status', 'objective', 'bound', 'elapsed_secs', 'nodes',
                                      'lp_iterations', 'presolve_reductions')}
    summary.update(first=first, checked_events=len(goals))
    if objective is not None and row['bound'] is not None and abs(objective) > 0:
        summary['gap_percent'] = 100 * abs(objective - row['bound']) / abs(objective)
    if objective is not None:
        summary['best_at_secs'] = next((event['elapsed_secs'] for event in row['incumbent_events']
                                        if abs(event['objective'] - objective) <= 1e-6), None)
    return summary


def encode(model, seed):
    """打乱所有行列、更名、等价缩放行，将一部分二元列改用有界整数表示。"""
    rng = random.Random(seed)
    columns = list(range(len(model['variables'])))
    rows = list(model['rows'])
    rng.shuffle(columns)
    rng.shuffle(rows)
    names = {j: f'v{i}' for i, j in enumerate(columns)}
    declarations, bounds = [], []
    for j in columns:
        v = model['variables'][j]
        kind = 'integer' if v['kind'] == 'binary' and j % 2 == 0 else v['kind']
        declarations.append(f'{names[j]} {kind}')
        if v['lb'] != 0:
            bounds.append(f"{names[j]} >= {v['lb']:.17g}")
        if v['ub'] is not None:
            bounds.append(f"{names[j]} <= {v['ub']:.17g}")
    def terms(pairs, scale):
        return ' '.join(f'{a*scale:.17g} {names[j]}' for j, a in pairs)
    out = ['name transformed', 'sense '+model['sense'], 'vars', *declarations,
           'obj', terms(model['objective'], 1), 'cons']
    for r in rows:
        scale = rng.choice((0.5, -0.25, 2, -4))
        op = r['op'] if scale > 0 else {'<=': '>=', '>=': '<=', '=': '='}[r['op']]
        out.append(f"{terms(r['terms'], scale)} {op} {r['rhs']*scale:.17g}")
    out += ['bounds', *bounds, 'end']
    return ('\n'.join(out)+'\n').encode()


def shipment(seed):
    """固定开关时是整容量、整需求的网络流，线性最优可取整数，逐期枚举发货量。"""
    rng = random.Random(seed)
    m = dict(sense='min' if seed % 2 else 'max', variables=[], objective=[], rows=[])
    sign = 1 if m['sense'] == 'min' else -1
    def var(kind, cost=0):
        j = len(m['variables'])
        m['variables'].append(dict(kind=kind, lb=0, ub=1 if kind == 'binary' else None))
        if cost:
            m['objective'].append((j, cost * sign))
        return j
    def add(terms, op, rhs):
        m['rows'].append(dict(terms=terms, op=op, rhs=rhs))
    previous = [None, None]
    costs = {(0, 0): 0}
    for t in range(2 + seed % 2):
        cap, fee = rng.randint(1, 6), rng.randint(2, 15)
        demand = [rng.randint(0, 5) for _ in previous]
        penalty = [rng.randint(1, 12) for _ in previous]
        local_fee = [rng.randint(1, 10) for _ in previous]
        paid = [var('binary', c) for c in local_fee]
        exclusive = (t + seed) % 2 == 0
        flow = [var('continuous') for _ in previous]
        for i in range(2):
            remaining = var('continuous', penalty[i])
            terms = [(flow[i], 1), (remaining, 1)]
            if previous[i] is not None:
                terms.append((previous[i], -1))
            previous[i] = remaining
            add(terms, '=', demand[i])
            add([(flow[i], 1), (paid[i], -cap)], '<=', 0)
            add([(flow[i], 1), (var('binary'), -cap)], '<=', 0)
        add([(flow[0], 1), (flow[1], 1), (var('binary', fee), -cap)], '<=', 0)
        if exclusive:
            add([(paid[0], 1), (paid[1], 1)], '<=', 1)
        next_costs = {}
        for state, prev_cost in costs.items():
            a, b = [old + new for old, new in zip(state, demand)]
            for x in range(min(a, cap) + 1):
                for y in range(min(b, cap - x) + 1):
                    if exclusive and x and y:
                        continue
                    left = (a-x, b-y)
                    cost = (prev_cost + fee*bool(x+y) + local_fee[0]*bool(x)
                            + local_fee[1]*bool(y) + sum(i*j for i, j in zip(left, penalty)))
                    next_costs[left] = min(cost, next_costs.get(left, cost))
        costs = next_costs
    return encode(m, seed), sign * min(costs.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary', type=Path)
    parser.add_argument('runs', type=Path, help='含 solver 标签和原始结果 row 的 gzip JSONL')
    parser.add_argument('output', type=Path)
    parser.add_argument('--baseline', type=Path, help='重跑时提供 74275dc 的 round_probe')
    parser.add_argument('--gurobi', type=Path, help='重跑时提供 Gurobi 10.0.3 动态库')
    args = parser.parse_args()
    binary = args.binary.resolve()
    entry, = [c for c in json.loads((FROZEN/'manifest.json').read_bytes())['cases']
              if c['instance'] == 'neos-4954672-berkel']
    errors = []
    model, truth = load_case(FROZEN, entry, public_truths(), errors.append)
    assert not errors and truth['mps_rehashed'] and truth['status'] == 'optimal', errors
    check_native((FROZEN/entry['mps']).read_bytes(), model)
    if args.baseline or args.gurobi:
        if not (args.baseline and args.gurobi) or args.runs.exists():
            parser.error('重跑需要两种程序和一个尚未存在的原始记录路径')
        with tempfile.TemporaryDirectory() as scratch:
            path = Path(scratch) / 'original.milp'
            path.write_bytes(read_bytes(FROZEN/entry['input']))
            gurobi = json.loads(subprocess.check_output([sys.executable, ROUNDS/'compare_all_timing.py',
                '_worker', args.gurobi, FROZEN/entry['model'], '60'], timeout=120))
            rows = [('gurobi', args.gurobi, gurobi), ('baseline', args.baseline,
                     solve(args.baseline.resolve(), path, 60, 2000000)),
                    ('new', binary, solve(binary, path, 60, 2000000))]
            with gzip.open(args.runs, 'xt') as stream:
                for name, program, row in rows:
                    print(json.dumps(dict(solver=name, row=row,
                        sha256=hashlib.sha256(program.read_bytes()).hexdigest())), file=stream)
    report = dict(binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                  milp_sha256=entry['milp_sha256'], reference=truth, original={}, small=[],
                  baseline_commit='74275dc', seconds=60, threads=1, node_limit=2000000)
    with gzip.open(args.runs, 'rt') as stream:
        for run in map(json.loads, stream):
            summary = check(model, run['row'], truth['objective'])
            summary['program_sha256'] = run['sha256']
            if run['solver'] == 'new':
                assert run['sha256'] == report['binary_sha256']
            if run['solver'] == 'gurobi':
                summary['version'] = run['row']['version']
                summary['parameters'] = run['row']['parameters']
            report['original'][run['solver']] = summary
    assert set(report['original']) == {'baseline', 'new', 'gurobi'}
    with tempfile.TemporaryDirectory() as scratch:
        path = Path(scratch) / 'variant.milp'
        payload = encode(model, 3763)
        path.write_bytes(payload)
        permuted = solve(binary, path, 4, 2000000)
        report['permuted'] = check(parse_milp(payload), permuted, truth['objective'])
        columns = list(range(len(model['variables'])))
        random.Random(3763).shuffle(columns)
        for point in [permuted, *permuted['incumbent_events']]:
            if point['values'] is not None:
                point['values'] = [x for _, x in sorted(zip(columns, point['values'], strict=True))]
        check(model, permuted, truth['objective'])
        report['permuted']['sha256'] = hashlib.sha256(payload).hexdigest()
        for seed in range(371, 395):
            payload, truth = shipment(seed)
            path.write_bytes(payload)
            m = parse_milp(payload)
            row = dict(seed=seed, optimum=truth, runs=[])
            for nodes in (1, 2000000):
                result = check(m, solve(binary, path, 5, nodes), truth)
                if nodes > 1:
                    assert result['status'] == 'optimal', (seed, result)
                row['runs'].append(result)
            report['small'].append(row)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n')


if __name__ == '__main__':
    main()
