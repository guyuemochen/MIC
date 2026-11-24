#!/usr/bin/env python3
"""从默认单线程入口比较引擎。目标真值由生成时的原始可行点和独立费用界确定。"""
import argparse
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile

ROUNDS = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROUNDS))
from check_results import verify_solution
from summarize_all_timing import audit_truth, parse_milp, read_events


def encode(label, sense, kinds, costs, rows, bounds, rng):
    names = [f'z{i}' for i in rng.sample(range(len(kinds)), len(kinds))]

    def terms(pairs):
        return ' '.join(f'{a} {names[j]}' for j, a in pairs if a)

    columns = list(range(len(kinds)))
    rng.shuffle(columns)
    rng.shuffle(rows)
    lines = [f'name {label}', f'sense {sense}', 'vars', *(f'{names[j]} {kinds[j]}' for j in columns),
             'obj', terms(enumerate(costs)), 'cons',
             *(f'r{i}: {terms(row)} {op} {rhs}' for i, (row, op, rhs) in enumerate(rows)),
             'bounds', *(f'{names[j]} <= {b}' for j, b in bounds), 'end']
    return ('\n'.join(lines)+'\n').encode()


def matching(size, degree, mixed, seed):
    """完美匹配。每条边至少付两端势的和，预留的对角匹配恰付该价。"""
    rng = random.Random(seed)
    u, v = [[rng.randrange(1, 30) for _ in range(size)] for _ in range(2)]
    rows = [[] for _ in range(2*size)]
    kinds, costs = [], []
    for i in range(size):
        left = [k for k in range(size) if k != i]
        for k in [i, *rng.sample(left, degree - 1)]:
            j = len(costs)
            costs.append(u[i]+v[k] + (rng.randint(1, 30) if k != i else 0))
            kinds.append('binary' if mixed and rng.randrange(3) == 0 else 'continuous')
            rows[i].append((j, 1))
            rows[size+k].append((j, 1))
    key = f'matching_{size}_{degree}_{int(mixed)}_{seed}'
    return key, encode(key, 'min', kinds, costs, [(r, '=', 1) for r in rows], [], rng), sum(u)+sum(v)


def grid(side, mixed, seed):
    """有向格网络。一单位流从源到汇，弧价为势差加非负罚金，保留一条零罚金路径。"""
    rng = random.Random(seed)
    n = side*side
    potential = [rng.randrange(-30, 30) for _ in range(n)]
    rows = [[] for _ in range(n)]
    path = list(range(side))+list(range(2*side-1, n, side))
    spine = set(zip(path, path[1:]))
    costs, kinds = [], []
    for a in range(n):
        out = []
        if a % side < side-1:
            out.append(a+1)
        if a < n-side:
            out.append(a+side)
        for b in out:
            j = len(costs)
            costs.append(potential[b]-potential[a] + (0 if (a,b) in spine else rng.randrange(1,30)))
            kinds.append('binary' if mixed and rng.randrange(3) == 0 else 'continuous')
            rows[a].append((j, -1))
            rows[b].append((j, 1))
    key = f'grid_{side}_{int(mixed)}_{seed}'
    return key, encode(key, 'min', kinds, costs, [(r, '=', (i==n-1)-(i==0))
                       for i,r in enumerate(rows)], [], rng), potential[-1]-potential[0]


def capacity(n, width, banded, seed):
    """所有行在一个内点同时取等，目标等于正权约束之和，给出可达上界。"""
    rng = random.Random(seed)
    point = [rng.randrange(2, 9) for _ in range(n)]
    costs, rows = [0]*n, []
    for i in range(n):
        if banded:
            neighbors = [(i+k)%n for k in range(1, 2*width)]
        else:
            neighbors = [j for j in range(n) if i != j]
        selected = rng.sample(neighbors, width-1)
        row = [(j, rng.randrange(1, 4)) for j in selected]
        row.append((i, sum(a for _, a in row)+1))
        weight = rng.randrange(1, 12)
        for j,a in row:
            costs[j] += a*weight
        rows.append((row, '<=', sum(point[j]*a for j,a in row)))
    key = f'capacity_{n}_{width}_{int(banded)}_{seed}'
    return key, encode(key, 'max', ['continuous']*n, costs, rows, [], rng), sum(a*x for a,x in zip(costs, point))


def conflict_chain(n, width, seed):
    """带空闲收益的相邻作业冲突。空闲列等于 1 减去两端选择量，余下按最近 width 位动态规划。"""
    rng = random.Random(seed)
    weights = [rng.randint(1, 8) for _ in range(n)]
    costs = weights.copy()
    rows, kinds, masks = [], ['binary']*n, [0]*n
    for i in range(n):
        for k in range(1, min(width, i)+1):
            if k > 1 and rng.randrange(3) == 0:
                continue
            j, slack, premium = i-k, len(kinds), rng.randint(1, 3)
            kinds.append('continuous')
            costs.append(premium)
            costs[j] += premium
            costs[i] += premium
            rows.append(([(j, 1), (i, 1), (slack, 1)], '<=', 1))
            masks[i] |= 1 << (k-1)
    dp = {0: 0}
    for i, weight in enumerate(weights):
        scores = {}
        for bits, score in dp.items():
            shift = (bits << 1) & ((1 << width) - 1)
            scores[shift] = max(scores.get(shift, 0), score)
            if bits & masks[i] == 0:
                chosen = shift | 1
                scores[chosen] = max(scores.get(chosen, 0), score + weight)
        dp = scores
    key = f'conflict_{n}_{width}_{seed}'
    return key, encode(key, 'max', kinds, costs, rows, [], rng), sum(costs[n:])+max(dp.values())


def cases():
    yield matching(110, 5, True, 9317)
    yield grid(16, False, 2707)
    yield capacity(360, 4, True, 6323)
    yield capacity(360, 4, False, 6323)
    yield matching(160, 4, False, 2273)
    yield grid(19, True, 4517)
    yield capacity(240, 8, True, 6089)
    yield capacity(240, 8, False, 6089)
    yield capacity(240, 60, False, 6089)
    yield conflict_chain(64, 4, 7193)
    yield conflict_chain(112, 4, 6011)


def run(args):
    rows = {}
    models = cases()
    if args.density_boundary:
        models = (capacity(240, width, False, 6101) for width in (14, 15, 16))
    with tempfile.TemporaryDirectory() as scratch:
        for name, payload, optimum in models:
            path = Path(scratch) / f'{name}.milp'
            path.write_bytes(payload)
            model = parse_milp(payload)
            reduced = json.loads(subprocess.check_output([args.binary,'export-presolved',path]))
            m, n = len(reduced['rows']), len(reduced['variables'])
            result = {'sha256': hashlib.sha256(payload).hexdigest(), 'optimum': optimum,
                      'presolved_rows': m, 'presolved_cols': n,
                      'presolved_nnz': sum(len(r['terms']) for r in reduced['rows']), 'runs': []}
            for rep in range(args.reps):
                engines = args.engines if rep%2 == 0 else reversed(args.engines)
                for engine in engines:
                    cmd = [args.binary, 'timed-solve', path, '--time-limit-secs', str(args.seconds),
                           '--lp-engine', engine, '--node-limit', '2000000', '--threads', '1',
                           '--early-stop', 'off']
                    r = json.loads(subprocess.check_output(cmd, timeout=args.seconds+40))
                    errors = []
                    value = verify_solution(model, r, errors.append)
                    _, events, rejected = read_events(model, r, errors.append)
                    audit_truth(model, dict(status='optimal', objective=optimum), r,
                                events+([value] if value is not None else []), errors.append)
                    assert not (errors or rejected or r.get('callback_errors')), (name, engine, errors, rejected)
                    result['runs'].append(dict(engine=engine, rep=rep,
                        **{k:r[k] for k in ['status','objective','bound','elapsed_secs','lp_iterations','nodes']}))
            rows[name] = result
            print(name, [(x['engine'], x['status'], round(x['elapsed_secs'],4)) for x in result['runs']],flush=True)
            args.output.write_text(json.dumps(dict(probe_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(),
                seconds=args.seconds,threads=1,node_limit=2000000,models=rows), indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--seconds', default=12, type=int)
    parser.add_argument('--reps', default=3, type=int)
    parser.add_argument('--density-boundary', action='store_true',
                        help='仅测每行 14、15、16 项的 240 列容量矩阵，覆盖密度阈值两侧')
    parser.add_argument('--engines', nargs='+', choices=['auto','revised','tableau'], default=['tableau','revised'])
    run(parser.parse_args())
