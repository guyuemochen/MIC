#!/usr/bin/env python3
"""生成等价行变体，以原始输入验解，独立核对旧 MIC、新 MIC 和 Gurobi 的解与界。"""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import random
import subprocess
import sys

ROUNDS = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROUNDS), str(ROUNDS.parent / 'auto_validation')]
from generate import covering, encode, packing
from summarize_all_timing import parse_milp, read_events
from check_results import verify_solution

SECONDS = 12
NODES = 2_000_000
FLIP = {'<=': '>=', '>=': '<=', '=': '='}
PARALLEL = {'mip_redundant', 'redundant_lp', 'near_parallel'}
ONE_NODE = {'cover_817391_scaled', 'cover_905731_integer', 'pack_medium_905749_scaled', 'strict_cover'}


def transform(payload, factors, single=False):
    """直接改原文，一条约束可单独改变符号、系数和右端。"""
    header, tail = payload.decode().split('cons\n')
    body, end = tail.split('end\n')
    lines = []
    for i, line in enumerate(body.splitlines()):
        factor = factors[i % len(factors)]
        if (single and i) or factor == 1:
            lines.append(line)
            continue
        label, *terms, op, rhs = line.split()
        tokens = [label]
        for a, col in zip(terms[::2], terms[1::2], strict=True):
            tokens.extend((f'{float(a) * factor:g}', col))
        lines.append('  ' + ' '.join([*tokens, FLIP[op] if factor < 0 else op, f'{float(rhs) * factor:g}']))
    return (header + 'cons\n' + '\n'.join(lines) + '\nend\n' + end).encode()


def small(kinds, rows, cost=(3, 5), ub=12):
    names = ['x', 'y']
    body = ['name audit', 'sense max', 'vars']
    body += [f'{n} {k}' for n, k in zip(names, kinds, strict=True)]
    body += ['obj', f'{cost[0]} x {cost[1]} y', 'cons']
    for a, b, op, rhs in rows:
        body.append(f'{a} x {b} y {op} {rhs}')
    body += ['bounds', f'x <= {ub}', f'y <= {ub}', 'end']
    return ('\n'.join(body) + '\n').encode()


def redundant_lp():
    rng = random.Random(904217)
    n, m = 80, 100
    cols = [f'z{j}' for j in range(n)]
    weights = [rng.randrange(1, 21) for _ in cols]
    rows = []
    for _ in range(m):
        js = rng.sample(range(n), 8)
        a = [rng.randrange(1, 9) for j in js]
        cap = sum(a) // 2
        for k, extra in [(1, 0), (1, 1), (-4, 0), (3, 2), (.5, 0), (-5, 3)]:
            terms = ' '.join(f'{v*k:g} {cols[j]}' for j, v in zip(js, a, strict=True))
            rows.append(f"{terms} {'<=' if k > 0 else '>='} {(cap+extra)*k:g}")
    rng.shuffle(rows)
    lines = ['name redundant_lp', 'sense max', 'vars', *(f'{z} continuous' for z in cols), 'obj',
             ' '.join(f'{cost} {col}' for col, cost in zip(cols, weights, strict=True)), 'cons',
             *rows, 'bounds', *(f'{z} <= 1' for z in cols), 'end']
    return ('\n'.join(lines) + '\n').encode()


def cases():
    scales = (3, -.5, .25, -5, 2)
    for seed in (817391, 905731):
        label = f'cover_{seed}'
        # 首个种子使用上一轮反例相同的模型名，保留其未改写原文。
        payload = covering('fresh_c_0' if seed == 817391 else 'fresh', 137, 81, seed)
        yield label, label, payload
        yield label + '_twice', label, transform(payload, (2,), single=True)
        yield label + '_scaled', label, transform(payload, scales)
        if seed == 905731:
            integer_bounds = 'bounds\n' + ''.join(f'x{j} <= 1\n' for j in range(1,138))
            integer_payload = transform(payload, scales).decode().replace(' binary\n',' integer\n')
            yield label + '_integer', label, integer_payload.replace('end\n', integer_bounds+'end\n').encode()
    # 一行要求两列同时取一，其余行只要求至少一列，最小费用分别是 3、1。
    base = [([(0,1), (1,1)], '>=', 1)] * 32
    costs = [1,2] + [10]*62
    yield 'strict_cover_base', 'strict_cover_base', encode('audit','min',64,costs,base)
    yield 'strict_cover', 'strict_cover', encode('audit','min',64,costs,[([(0,2),(1,2)],'>=',2.5), *base[1:]])
    for seed in (905743, 905747):
        label = f'pack_{seed}'
        payload = packing('fresh', 194, 9, 1, 28, False, seed)
        yield label, label, payload
        yield label + '_scaled', label, transform(payload, scales)
    label = 'pack_medium_905749'
    payload = packing('fresh', 194, 9, 1, 6, False, 905749)
    yield label, label, payload
    yield label + '_scaled', label, transform(payload, scales)
    for kind, types in [('lp', ('continuous', 'continuous')), ('mip', ('integer', 'continuous')),
                        ('ip', ('integer', 'integer'))]:
        base = [(1, 2, '<=', 9), (3, 1, '<=', 11)]
        extra = [(-.5, -1, '>=', -4.5), (3, 6, '<=', 30), (-15, -5, '>=', -55),
                 (1.5, .5, '<=', 6)]
        yield kind, kind, small(types, base)
        yield kind + '_redundant', kind, small(types, list(reversed(extra + base)))
    types = ('continuous', 'continuous')
    eqs = [(1, 2, '=', 9), (3, 1, '<=', 11)]
    yield 'equalities', 'equalities', small(types, eqs)
    yield 'equalities_redundant', 'equalities', small(types, [(-3, -6, '=', -27), (2, 4, '>=', 16),
                                                             (.5, 1, '<=', 4.5), *eqs])
    yield 'incompatible', 'incompatible', small(types, [(2, 4, '<=', 14), (-3, -6, '<=', -24)])
    yield 'near_parallel', 'near_parallel', small(
        types, [(1, 1, '=', 32768), (1, 1 + 2**-21, '<=', 32768 + 2**-7)], (0, 1), 65536)
    # 支撑列相同且有系数变号的非平行行，也要保留。每个模型的真值独立枚举。
    for k in range(9):
        rng = random.Random(271907+k)
        kinds = [('continuous','continuous'), ('integer','continuous'), ('integer','integer')][k%3]
        tx, ty = rng.randrange(2, 9), rng.randrange(2, 9)
        rows = []
        for _ in range(5):
            a, b = rng.choice((-3,-2,1,2,3)), rng.choice((-3,-2,1,2,3))
            rhs = a*tx+b*ty+rng.randrange(4)
            rows.extend(((a,b,'<=',rhs), (-3*a,-3*b,'>=',-3*rhs), (.25*a,.25*b,'<=',.25*(rhs+1))))
        rng.shuffle(rows)
        key = f'rational_{k}'
        yield key, key, small(kinds, rows, (rng.randrange(1,9),rng.randrange(1,9)))
    yield 'redundant_lp', 'redundant_lp', redundant_lp()


def oracle(model):
    """两变量有界模型。LP 枚举顶点，IP 枚举整数，混合模型枚举整数切片端点。"""
    if len(model['variables']) != 2:
        return None
    variables = model['variables']
    upper = [int(v['ub']) for v in variables]
    rows = [(Fraction(1), 0, upper[0]), (-1, 0, 0), (0, 1, upper[1]), (0, -1, 0)]
    for r in model['rows']:
        c = dict(r['terms'])
        a, b, h = Fraction(c.get(0, 0)), Fraction(c.get(1, 0)), Fraction(r['rhs'])
        if r['op'] in ('<=', '='):
            rows.append((a, b, h))
        if r['op'] in ('>=', '='):
            rows.append((-a, -b, -h))
    # 保证每个顶点都用有理数，不经过二进制浮点除法。
    rows = [tuple(map(Fraction, r)) for r in rows]
    if variables[1]['kind'] == 'integer':
        points = ((x, y) for x in range(upper[0]+1) for y in range(upper[1]+1))
    elif variables[0]['kind'] == 'integer':
        points = ((x, (h-a*x)/b) for x in range(upper[0]+1) for a,b,h in rows if b)
    else:
        points = (((h*d-b*k)/(a*d-b*c), (a*k-h*c)/(a*d-b*c))
                  for (a,b,h),(c,d,k) in combinations(rows, 2) if a*d != b*c)
    cost = dict(model['objective'])
    best = max((Fraction(cost.get(0, 0))*x + Fraction(cost.get(1, 0))*y for x,y in points
                if all(a*x+b*y <= h for a,b,h in rows)), default=None)
    return ('infeasible', None) if best is None else ('optimal', float(best))


def scale_sweep(binary, scratch, role, binary_hash):
    """跨越 2⁻²⁸ 到 2²⁸，对照有理数真值，所有计算都使用实际写出的原行。"""
    counts = {'optimal': 0, 'infeasible': 0}
    callbacks = 0
    inputs = sha256()
    kinds = [('integer', 'integer'), ('continuous', 'continuous'), ('integer', 'continuous')]
    for k in range(140):
        rng = random.Random(837859 + k)
        kind = rng.choice(kinds)
        x, y = rng.randrange(10), rng.randrange(10)
        rows = []
        for _ in range(4):
            a = rng.choice((-3, -2, -1, 1, 2, 3))
            b = rng.choice((-9, -3, -1, 1, 3, 9))
            limit = a*x + b*y + rng.randrange(-1, 5)
            first = rng.choice((-1, 1)) * 2.0 ** rng.randrange(-28, 29)
            second = rng.choice((.5, 2, 3, -4, -5))
            loose = rng.choice((0, 0, 1))
            for factor, rhs in ((first, limit), (second, limit + loose)):
                op = '<=' if factor > 0 else '>='
                rows.append((factor*a, factor*b, op, factor*rhs))
        rng.shuffle(rows)
        data = small(kind, rows, (rng.randrange(1, 8), rng.randrange(1, 8)))
        source = scratch / f'scale_{k}.milp'
        source.write_bytes(data)
        model = parse_milp(data)
        status, expected = oracle(model)
        cmd = [binary, 'timed-solve', source, '--threads', '1', '--time-limit-secs', str(SECONDS),
               '--node-limit', str(NODES), '--early-stop', 'off']
        run = json.loads(subprocess.check_output(cmd, timeout=SECONDS+20))
        value = audit(model, run)
        assert run['status'] == status, (k, run['status'], status)
        if expected is not None:
            assert abs(value-expected) < 1e-6 and abs(run['bound']-expected) < 1e-6, (k, value, expected)
        counts[status] += 1
        inputs.update(sha256(data).digest())
        callbacks += len(run['incumbent_events'])
    record = dict(binary_sha256=binary_hash, ordered_inputs_sha256=inputs.hexdigest(),
                  counts=counts, callback_vectors=callbacks)
    (scratch / f'{role}-sweep.json').write_text(json.dumps(record))


def run(role, binary, scratch):
    binary = binary.resolve()
    digest = sha256(binary.read_bytes()).hexdigest()
    for name, group, payload in cases():
        source, exported = scratch / (name + '.milp'), scratch / (name + '.model.json')
        source.write_bytes(payload)
        model = parse_milp(payload)
        exported.write_text(json.dumps(model))
        if role == 'gurobi':
            cmd = [sys.executable, ROUNDS/'compare_all_timing.py', '_worker', binary, exported, str(SECONDS)]
        else:
            got = json.loads(subprocess.check_output([binary, 'export', source]))
            assert all(model[k] == v for k,v in got.items()), name
            cmd = [binary, 'timed-solve', source, '--threads', '1', '--time-limit-secs', str(SECONDS),
                   '--node-limit', str(NODES), '--early-stop', 'off']
        raw = subprocess.check_output(cmd, text=True, timeout=SECONDS+20).strip().splitlines()[-1]
        row = json.loads(raw)
        checked = audit(model, row)
        (scratch/f'{role}-{name}.json').write_text(json.dumps(dict(
            input_sha256=sha256(payload).hexdigest(), binary_sha256=digest, row=row)))
        if role == 'new':
            extra_checks(binary, source, name, model, scratch, digest)
        print(role, name, row['status'], checked, row['lp_iterations'], round(row['elapsed_secs'], 6), flush=True)
    if role in ('old', 'new'):
        scale_sweep(binary, scratch, role, digest)


def extra_checks(binary, source, name, model, scratch, digest):
    """多线程走公开无回调入口，节点限额走单线程带回调入口。"""
    if name in PARALLEL:
        flags, threads, limit = 'solve', 2, NODES
    elif name in ONE_NODE:
        flags, threads, limit = 'timed-solve', 1, 1
    else:
        return
    cmd = [binary, flags, source, '--threads', str(threads), '--time-limit-secs', str(SECONDS),
           '--node-limit', str(limit), '--early-stop', 'off']
    r = json.loads(subprocess.check_output(cmd, timeout=SECONDS+20))
    audit(model, r)
    (scratch/f'extra-{name}.json').write_text(json.dumps(dict(
        binary_sha256=digest, input_sha256=sha256(source.read_bytes()).hexdigest(), threads=threads,
        limit=limit, row=r)))


def audit(model, row):
    errors = []
    checked = verify_solution(model, row, errors.append)
    _, _, rejected = read_events(model, row, errors.append)
    assert not errors and not rejected and not row.get('callback_errors'), (errors, rejected)
    assert row['status'] in ('infeasible', 'optimal', 'limit_reached')
    # 全二元模型再作无容差复核：所有系数乘 4 均为整数；对披露的每个向量取整再验。
    if all(v['kind'] in ('integer','binary') and (v['lb'],v['ub']) == (0,1) for v in model['variables']):
        records = ([row] if row['values'] is not None else []) + row.get('incumbent_events', [])
        for rec in records:
            x = [round(value) for value in rec['values']]
            assert all(v in (0,1) for v in x)
            for r in model['rows']:
                lhs = sum(int(4*a)*x[j] for j,a in r['terms'])
                rhs = int(4*r['rhs'])
                assert all(4*a == int(4*a) for j,a in r['terms']) and rhs == 4*r['rhs']
                assert {'<=':lhs <= rhs, '>=':lhs >= rhs, '=':lhs == rhs}[r['op']]
            assert abs(sum(a*x[j] for j,a in model['objective'])-rec['objective']) < 1e-6
    return checked


def summarize(scratch, output):
    results, hashes, witnesses, events = {}, {}, 0, 0
    truths = {}
    for name, group, payload in cases():
        model, digest = parse_milp(payload), sha256(payload).hexdigest()
        records = {}
        for role in ('old','new','gurobi'):
            item = json.loads((scratch/f'{role}-{name}.json').read_text())
            assert item['input_sha256'] == digest
            assert hashes.setdefault(role, item['binary_sha256']) == item['binary_sha256']
            records[role] = item['row']
        reference = records['gurobi']
        assert reference['status'] in ('infeasible','optimal') and reference['version'] == '10.0.3'
        assert reference['model_roundtrip_passed']
        params = reference['parameters']
        assert (params['Threads'], params['TimeLimit'], params['NodeLimit'], params['MIPGap'],
                params['MIPGapAbs']) == (1, SECONDS, NODES, 0., 0.)
        truth = (reference['status'], reference['objective'])
        exact = oracle(model)
        expected = exact or truths.get(group, truth)
        assert expected[0] == truth[0] and (truth[1] is None or abs(expected[1] - truth[1]) < 1e-6)
        truths.setdefault(group, truth)
        entry = dict(group=group, model_sha256=digest, rational_oracle=exact, runs={})
        for role, row in records.items():
            obj = audit(model, row)
            status, bound = row['status'], row['bound']
            assert status == 'limit_reached' or status == truth[0], (role, name, status, truth)
            if status == 'optimal':
                assert abs(obj - truth[1]) < 1e-6 and abs(bound-truth[1]) < 1e-6
            if truth[0] == 'infeasible':
                assert obj is None
            elif bound is not None:
                assert (bound - truth[1]) * (1 if model['sense']=='min' else -1) <= 1e-6
            count = len(row['incumbent_events'])
            witnesses += obj is not None
            events += count
            entry['runs'][role] = dict(status=status, objective=obj, bound=bound, callback_count=count,
                elapsed_secs=row['elapsed_secs'], nodes=row['nodes'], lp_iterations=row['lp_iterations'],
                reductions=row.get('presolve_reductions'))
        extra = scratch/f'extra-{name}.json'
        if name in ONE_NODE or name in PARALLEL:
            item = json.loads(extra.read_text())
            assert digest == item['input_sha256'] and item['binary_sha256'] == hashes['new']
            r = item['row']
            obj = audit(model, r)
            assert truth[0] == 'optimal' and r['status'] in ('optimal','limit_reached')
            sign = 1 if model['sense'] == 'min' else -1
            bound = r['bound']
            assert bound is None or sign*(bound-truth[1]) <= 1e-6
            assert obj is None or sign*(truth[1]-obj) <= 1e-6
            assert r['status'] != 'optimal' or abs(obj-truth[1]) < 1e-6
            assert item['threads'] != 2 or r['status'] == 'optimal'
            assert item['limit'] != 1 or r['nodes'] <= 1
            entry['additional'] = dict(threads=item['threads'], node_limit=item['limit'], status=r['status'],
                objective=obj, bound=bound, nodes=r['nodes'], reductions=r.get('presolve_reductions'),
                callback_count=len(r.get('incumbent_events',[])))
        results[name] = entry
    sweep = {role: json.loads((scratch / f'{role}-sweep.json').read_bytes()) for role in ('old', 'new')}
    for role, item in sweep.items():
        assert item['binary_sha256'] == hashes[role] and item['counts'] == sweep['old']['counts']
        assert item['ordered_inputs_sha256'] == sweep['old']['ordered_inputs_sha256']
    report = dict(exact_scale_sweep=sweep, gurobi_version='10.0.3', baseline_commit='0b47f2a', threads=1, repetitions=1, seconds=SECONDS, nodes=NODES, binary_sha256=hashes,
                  final_vectors=witnesses, callback_vectors=events, cases=results)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('old','new','gurobi','summarize'))
    p.add_argument('--binary', type=Path)
    p.add_argument('--scratch', type=Path, required=True)
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    args.scratch.mkdir(parents=True, exist_ok=True)
    if args.action == 'summarize':
        assert args.output
        summarize(args.scratch, args.output)
    else:
        assert args.binary
        run(args.action, args.binary, args.scratch)


if __name__ == '__main__':
    main()
