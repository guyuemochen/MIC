#!/usr/bin/env python3
"""直接从 s1234 的 MPS 恢复原行、以整数验解，并独立检查重排和穷举模型。"""
import argparse
from collections import defaultdict
import gzip
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
import random
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROUNDS = HERE.parents[2]
FROZEN = ROUNDS / 'results/main-gurobi10-100-60s'
sys.path.insert(0, str(ROUNDS))
from summarize_all_timing import parse_milp


def native():
    entry, = [x for x in json.loads((FROZEN / 'manifest.json').read_bytes())['cases'] if x['instance'] == 's1234']
    mps = gzip.decompress((FROZEN / entry['mps']).read_bytes())
    milp = gzip.decompress((FROZEN / entry['input']).read_bytes())
    assert sha256(mps).hexdigest() == entry['mps_sha256']
    assert sha256(milp).hexdigest() == entry['milp_sha256']
    sections, kinds, by_row, rhs, columns, binary = [], {}, defaultdict(dict), {}, {}, set()
    integers = False
    for line in mps.decode('ascii').splitlines():
        if not line or line.startswith('*'):
            continue
        f = line.split()
        if not line[0].isspace():
            sections.append(f[0])
        elif sections[-1] == 'ROWS':
            kind, name = f
            assert name not in kinds and kind in ('N', 'L', 'E', 'G')
            kinds[name] = kind
        elif sections[-1] == 'COLUMNS':
            if f[1] == "'MARKER'":
                assert f[2] == ("'INTEND'" if integers else "'INTORG'")
                integers = not integers
                continue
            assert integers
            j = columns.setdefault(f[0], len(columns))
            for row, value in zip(f[1::2], f[2::2], strict=True):
                assert row in kinds and j not in by_row[row]
                by_row[row][j] = int(value)
        elif sections[-1] == 'RHS':
            assert f[0] == 'RHS1'
            for row, value in zip(f[1::2], f[2::2], strict=True):
                assert row in kinds and row not in rhs
                rhs[row] = int(value)
        elif sections[-1] == 'BOUNDS':
            kind, vector, name = f
            assert (kind, vector) == ('BV', 'BND1') and name in columns and name not in binary
            binary.add(name)
        else:
            raise ValueError(f'意外的 MPS 段 {sections[-1]}')
    assert sections == ['NAME', 'ROWS', 'COLUMNS', 'RHS', 'BOUNDS', 'ENDATA'] and not integers
    assert binary == columns.keys()
    obj, = [row for row, k in kinds.items() if k == 'N']
    assert rhs.get(obj, 0) == 0
    costs = [by_row[obj].get(i, 0) for i in range(len(columns))]
    ops = {'E': '=', 'L': '<=', 'G': '>='}
    rows = [(sorted(by_row[r].items()), ops[k], rhs.get(r, 0)) for r, k in kinds.items() if k != 'N']
    parsed = parse_milp(milp)
    assert parsed['sense'] == 'min' and parsed['offset'] == 0
    assert parsed['variables'] == [dict(kind='binary', name=f'x{i}', lb=0, ub=1) for i in range(len(columns))]
    assert parsed['objective'] == [[i, c] for i, c in enumerate(costs) if c]
    assert parsed['rows'] == [dict(terms=[[j, a] for j, a in terms], op=op, rhs=b) for terms, op, b in rows]
    compiled = gzip.decompress((FROZEN / entry['model']).read_bytes())
    assert sha256(compiled).hexdigest() == entry['model_sha256']
    assert parsed == json.loads(compiled)
    return costs, rows, milp, {key: entry[key] for key in ('mps_sha256', 'milp_sha256', 'model_sha256')}


def holds(activity, op, rhs):
    return activity == rhs if op == '=' else activity >= rhs if op == '>=' else activity <= rhs


def witness(costs, rows, vector, claimed):
    assert len(vector) == len(costs) and all(isfinite(x) and abs(x - round(x)) <= 1e-6 and round(x) in (0, 1) for x in vector)
    bits = [round(x) for x in vector]
    for i, (terms, op, rhs) in enumerate(rows):
        activity = sum(a * bits[j] for j, a in terms)
        assert holds(activity, op, rhs), (i, activity, op, rhs)
        raw = sum(a * vector[j] for j, a in terms)
        violation = abs(raw - rhs) if op == '=' else rhs - raw if op == '>=' else raw - rhs
        assert violation <= 1e-6, (i, raw, op, rhs)
    obj = sum(c * x for c, x in zip(costs, bits, strict=True))
    assert isfinite(claimed) and abs(claimed - obj) <= 1e-6
    return bits


def activation_bound(costs, rows):
    """若覆盖行中的 j 取一，则只看同时含 j 的短行也有费用下界，剩余变量费用非负。"""
    long, = [r for r in rows if len(r[0]) > 10]
    terms, op, rhs = long
    assert (op, rhs) in (('>=', 1), ('<=', -1)) and all(a == rhs for _, a in terms) and min(costs) >= 0
    posting = defaultdict(list)
    for r in rows:
        if r is not long:
            for j, _ in r[0]:
                posting[j].append(r)
    cached, anchor_bounds = {}, []
    for j, _ in terms:
        local = posting[j]
        cols = sorted({j} | {col for t, _, _ in local for col, _ in t})
        assert len(cols) <= 12
        template = (cols.index(j), tuple(costs[c] for c in cols),
                    tuple((tuple((cols.index(c), a) for c, a in t), op, rhs) for t, op, rhs in local))
        if template not in cached:
            pivot, price, local_rows = template
            feasible = []
            for mask in range(1 << len(cols)):
                if (mask >> pivot) & 1 and all(holds(sum(a * ((mask >> c) & 1) for c, a in t), op, rhs)
                                                for t, op, rhs in local_rows):
                    feasible.append(sum(a * ((mask >> c) & 1) for c, a in enumerate(price)))
            cached[template] = min(feasible, default=None)
        anchor_bounds.append(cached[template])
    return dict(anchors=len(anchor_bounds), impossible_anchors=anchor_bounds.count(None),
                unique_local_tables=len(cached), lower_bound=min((v for v in anchor_bounds if v is not None), default=None))


def write_problem(path, costs, rows, names):
    def terms(t):
        return ' '.join(f'{a} {names[j]}' for j, a in t)
    with path.open('w') as f:
        print('sense min\nvars', file=f)
        for name in names:
            print(name, 'binary', file=f)
        print('obj', terms((j, c) for j, c in enumerate(costs) if c), 'cons', sep='\n', file=f)
        for i, (t, op, rhs) in enumerate(rows):
            print(f'r{i}:', terms(t), op, rhs, file=f)
        print('end', file=f)


def solve(binary, model, node_limit):
    command = [str(binary), 'timed-solve', str(model), '--threads', '1', '--node-limit', str(node_limit),
               '--time-limit-secs', '5', '--early-stop', 'off']
    return json.loads(subprocess.run(command, check=True, capture_output=True, timeout=12).stdout)


def audit(costs, rows, report, reorder=None):
    def check(record):
        vector = record['values']
        if reorder is not None:
            vector = [vector[i] for i in reorder]
        return witness(costs, rows, vector, record['objective'])

    assert not report.get('callback_errors') and report.get('raw_returncode', 0) == 0
    events = report['incumbent_events']
    times = [e['elapsed_secs'] for e in events]
    assert all(0 <= t <= report['elapsed_secs'] + 0.01 for t in times) and times == sorted(times)
    for event in events:
        check(event)
    bits = None
    if report['objective'] is not None:
        bits = check(report)
        assert all(report['objective'] <= e['objective'] + 1e-6 for e in events)
    else:
        assert not events and report['values'] is None
    bound, status, obj = report['bound'], report['status'], report['objective']
    if bound is not None:
        assert isfinite(bound) and (obj is None or bound <= obj + 1e-6)
    if status == 'optimal':
        assert obj is not None and bound is not None and abs(obj - bound) <= 1e-6
    elif status == 'infeasible':
        assert obj is None
    else:
        assert status == 'limit_reached'
    return dict(status=status, objective=obj, bound=bound, nodes=report['nodes'],
                callback_objectives=[event['objective'] for event in events],
                final_ones=sum(bits) if bits is not None else None)


def permuted(binary, tmp, costs, rows):
    rng = random.Random(817)
    permutation = list(range(len(costs)))
    rng.shuffle(permutation)
    to_new = [None] * len(costs)
    for i, j in enumerate(permutation):
        to_new[j] = i
    order = list(range(len(rows)))
    rng.shuffle(order)
    shuffled = [([(to_new[j], a) for j, a in rows[i][0]], rows[i][1], rows[i][2]) for i in order]
    path = tmp / 'renamed.milp'
    write_problem(path, [costs[j] for j in permutation], shuffled, [f'b{i}' for i in permutation])
    raw = solve(binary, path, 256)
    result = audit(costs, rows, raw, to_new)
    assert result['objective'] is not None, '重排模型在 256 节点内未给出可核验候选'
    result['model_sha256'] = sha256(path.read_bytes()).hexdigest()
    return result


def small_cases():
    # 每两个触发列组成一对。只选一列须支付共享的失配费；同选一对可免交该费用。
    for k in range(4):
        rng = random.Random(19 + k)
        costs = [rng.randrange(1, 5) for _ in range(11)] + [25, rng.randrange(1, 5)]
        rows = [([(j, 1) for j in range(11)], '>=', 1), ([(10, 1)], '=', 0)]
        for a in range(0, 10, 2):
            for i, j in ((a, a + 1), (a + 1, a)):
                rows.append(([(i, 1), (j, -1), (11, -1)], '<=', 0))
        rows.append(([(2, 1), (4, 1), (12, -1)], '<=', 1))
        if k == 1:
            rows.append(([(11, 1)], '=', 0))
        if k == 2:
            rows.extend([([(j, 1)], '=', 0) for j in range(10)])
        if k == 3:
            # 等价翻转长行方向，取第 0 列时还须支付最后一列的费用。
            rows[0] = ([(j, -1) for j in range(11)], '<=', -1)
            rows.append(([(0, 1), (12, -1)], '<=', 0))
        yield ['costly_singleton', 'no_singleton', 'infeasible', 'negative_anchor_row'][k], costs, rows
    # 每个触发列在两条不同的短行里都强制同一费用为 9 的辅助列取一。
    # 两行各自给的条件下界都是 10，合并时重复计入这笔费用会超过真值。
    costs = [1] * 11 + [9, 0, 0]
    rows = [([(j, 1) for j in range(11)], '>=', 1), ([(10, 1)], '=', 0)]
    for j in range(10):
        for aux in (12, 13):
            rows.append(([(j, 2), (11, -2), (aux, 1)], '<=', 1))
    yield 'shared_cost', costs, rows


def enumerate_small(binary, tmp):
    records = []
    for name, costs, rows in small_cases():
        solutions, singles = [], []
        for mask in range(1 << len(costs)):
            bits = [(mask >> j) & 1 for j in range(len(costs))]
            if all(holds(sum(a * bits[j] for j, a in t), op, rhs) for t, op, rhs in rows):
                obj = sum(c * x for c, x in zip(costs, bits, strict=True))
                solutions.append(obj)
                if sum(bits[:11]) == 1:
                    singles.append(obj)
        opt = min(solutions, default=None)
        singleton = min(singles, default=None)
        independent = activation_bound(costs, rows)
        assert independent['impossible_anchors'] > 0
        if opt is not None:
            assert independent['lower_bound'] is not None and independent['lower_bound'] <= opt
        path = tmp / (name + '.milp')
        write_problem(path, costs, rows, [f'b{i}' for i in range(len(costs))])
        reports = {}
        for limit in (1, 20000):
            r = audit(costs, rows, solve(binary, path, limit))
            if limit == 20000:
                assert (r['status'], r['objective']) == (('infeasible', None) if opt is None else ('optimal', opt)), (name, r)
            if opt is not None:
                assert r['status'] != 'infeasible' and (r['bound'] is None or r['bound'] <= opt + 1e-6), (name, opt, r)
                assert r['objective'] is None or r['objective'] >= opt
            else:
                assert r['objective'] is None
            reports[str(limit)] = r
        if name == 'no_singleton':
            assert opt is not None and singleton is None
        elif name == 'shared_cost':
            assert opt == singleton == independent['lower_bound'] == 10
        elif name != 'infeasible':
            assert opt < singleton
        records.append(dict(case=name, assignments=1 << len(costs), feasible=len(solutions),
                            optimum=opt, best_singleton=singleton, independent_bound=independent,
                            runs_by_node_limit=reports))
    return records


def read_reports(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        for line in f:
            r = json.loads(line)
            if r['instance'] == 's1234':
                yield r


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--runs', type=Path, default=HERE.parent / 'runs.jsonl.gz')
    parser.add_argument('--output', type=Path, default=HERE / 'results.json')
    args = parser.parse_args()
    costs, rows, _, hashes = native()
    result = dict(input_hashes=hashes, runs_sha256=sha256(args.runs.read_bytes()).hexdigest(),
                  probe_sha256=sha256(args.binary.read_bytes()).hexdigest(),
                  binary_columns=len(costs), original_rows=len(rows),
                  independent_bound=activation_bound(costs, rows))
    jobs = []
    for source, path in [('previous', FROZEN / 'timings.jsonl.gz'), ('current', args.runs)]:
        source_jobs = list(read_reports(path))
        solvers = {'mic', 'gurobi'} if source == 'previous' else {'old', 'mic', 'gurobi'}
        assert {r['solver'] for r in source_jobs} == solvers, (path, solvers)
        keys = {(r['solver'], r['rep']) for r in source_jobs}
        assert len(keys) == len(source_jobs)
        for r in source_jobs:
            assert all(r.get(k, val) == val for k, val in hashes.items())
            jobs.append(dict(dataset=source, solver=r['solver'], rep=r['rep'], **audit(costs, rows, r)))
    result['original_runs'] = jobs
    with TemporaryDirectory() as d:
        result['permutation'] = permuted(args.binary.resolve(), Path(d), costs, rows)
        result['exhaustive'] = enumerate_small(args.binary.resolve(), Path(d))
    checked = [*jobs, result['permutation']]
    upper = min(r['objective'] for r in checked if r['objective'] is not None)
    reference = json.loads((FROZEN / 'verification.json').read_bytes())['extra_official_truth']['s1234']
    assert reference['status'] == 'optimal' and result['independent_bound']['lower_bound'] <= reference['objective'] <= upper
    assert all(r['status'] != 'infeasible' for r in checked)
    assert all(r['bound'] is None or r['bound'] <= reference['objective'] + 1e-6 for r in checked)
    result['published_reference'] = reference
    result['best_verified_feasible'] = upper
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(f'已核实 {len(jobs)} 次原模型运行、1 个等价重排与 {len(result["exhaustive"])} 个穷举模型。')


if __name__ == '__main__':
    main()
