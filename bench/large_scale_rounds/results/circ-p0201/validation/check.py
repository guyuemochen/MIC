#!/usr/bin/env python3
"""独立核对生产求解器输出。使用官方 MPS、等价变换和小模型完全枚举。"""
import argparse
from collections import Counter
from dataclasses import dataclass
import gzip
from hashlib import sha256
import itertools
import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROUNDS = HERE.parents[2]
sys.path.insert(0, str(ROUNDS))
from summarize_all_timing import parse_milp

ARCHIVE = ROUNDS / 'results/all-gurobi10-remaining-60s'
NODE_CAP = 2_000_000


@dataclass
class Model:
    sense: str
    costs: list
    rows: list  # 每项为（有序列号与系数，关系，右端）。
    names: list

    def text(self):
        def expression(terms):
            return ' '.join(f'{a} {self.names[j]}' for j, a in terms)
        lines = ['name audit', 'sense ' + self.sense, 'vars', *(name + ' binary' for name in self.names),
                 'obj', expression([(j, a) for j, a in enumerate(self.costs) if a]), 'cons']
        for terms, op, rhs in self.rows:
            lines.append(f'{expression(terms)} {op} {rhs}')
        return '\n'.join(lines + ['end', ''])

    def feasible(self, bits):
        for terms, op, rhs in self.rows:
            lhs = sum(a * bits[j] for j, a in terms)
            if ((op == '<=' and lhs > rhs) or (op == '>=' and lhs < rhs)
                    or (op == '=' and lhs != rhs)):
                return False
        return True

    def objective(self, bits):
        return sum(c * bit for c, bit in zip(self.costs, bits, strict=True))


def native():
    meta = json.loads((ROUNDS / 'results/public-sources.json').read_bytes())
    source, = [x for x in meta['instances'] if x['instance'] == 'p0201']
    mps = gzip.decompress((ARCHIVE / 'mps/p0201.mps.gz').read_bytes())
    mic = gzip.decompress((ARCHIVE / 'inputs/p0201.milp.gz').read_bytes())
    assert sha256(mps).hexdigest() == source['mps_sha256']
    assert sha256(mic).hexdigest() == source['milp_sha256']
    optimum = int(re.search(rb'^\*BEST SOLN:\s*(\d+) \(opt\)', mps, re.M)[1])

    # 该输入没有目标方向、常数或多 RHS 集。按原文件出现顺序恢复行列；
    # RHS 中没出现的行右端为 0，整数标记范围内的列默认下界为 0。
    kinds, rows, rhs, columns, upper = {}, {}, {}, {}, {}
    sections = []
    integer = False
    for line in mps.decode('ascii').splitlines():
        if not line or line.startswith('*'):
            continue
        fields = line.split()
        if not line[0].isspace():
            sections.append(fields[0])
        elif sections[-1] == 'ROWS':
            kind, name = fields
            assert name not in kinds and kind in ('N', 'L')
            kinds[name], rows[name] = kind, {}
        elif sections[-1] == 'COLUMNS':
            if fields[1] == "'MARKER'":
                assert fields[2] == ("'INTEND'" if integer else "'INTORG'")
                integer = not integer
                continue
            assert integer
            j = columns.setdefault(fields[0], len(columns))
            for name, coeff in zip(fields[1::2], fields[2::2], strict=True):
                assert j not in rows[name]
                rows[name][j] = int(coeff)
        elif sections[-1] == 'RHS':
            assert fields[0] == 'RHS'
            for name, value in zip(fields[1::2], fields[2::2], strict=True):
                assert name in kinds and name not in rhs
                rhs[name] = int(value)
        elif sections[-1] == 'BOUNDS':
            bound, vector, name, value = fields
            assert (bound, vector, value) == ('UP', 'ONE', '1') and name not in upper
            upper[name] = 1
        else:
            raise ValueError(f'未知 MPS 语段 {sections[-1]}')
    assert sections == ['NAME', 'ROWS', 'COLUMNS', 'RHS', 'BOUNDS', 'ENDATA'] and not integer
    assert upper.keys() == columns.keys()
    objective, = [name for name, op in kinds.items() if op == 'N']
    assert objective not in rhs
    model = Model('min', [rows[objective].get(j, 0) for j in range(len(columns))],
                  [(sorted(rows[name].items()), '<=', rhs.get(name, 0))
                   for name, kind in kinds.items() if kind == 'L'], list(columns))

    imported = parse_milp(mic)
    assert imported['sense'] == 'min' and imported['offset'] == 0
    assert len(imported['variables']) == len(model.names) == 201 and len(model.rows) == 133
    for j, var in enumerate(imported['variables']):
        assert var == dict(name=f'x{j}', kind='binary', lb=0, ub=1)
    assert imported['objective'] == [[j, a] for j, a in enumerate(model.costs) if a]
    assert imported['rows'] == [dict(terms=[list(x) for x in t], op=op, rhs=b) for t, op, b in model.rows]
    return model, optimum, mic, source


def transform(model, seed, objective_sign=1, shift_group=(), shift=0):
    rng = random.Random(seed)
    order = list(range(len(model.names)))
    rng.shuffle(order)
    inverse = {old: new for new, old in enumerate(order)}
    rows = []
    for k, (terms, op, rhs) in enumerate(model.rows):
        sign = -1 if k % 3 == 0 else 1
        if sign < 0:
            op = {'=': '=', '>=': '<=', '<=': '>='}[op]
        rows.append((sorted((inverse[j], sign * a) for j, a in terms), op, rhs * sign))
    rng.shuffle(rows)
    costs = [(model.costs[old] + (shift if old in shift_group else 0)) * objective_sign for old in order]
    sense = model.sense if objective_sign == 1 else {'max': 'min', 'min': 'max'}[model.sense]
    return Model(sense, costs, rows, [f'u{i}' for i in range(len(order))]), order


def small(seed):
    # 分组是定义域的一部分，每组恰有一位为 1。两个额外位可分别取 0 或 1。
    # 固定一组可行点后生成嵌套累计约束，费用包含两个符号。
    rng = random.Random(seed)
    groups = [list(range(3 * k, 3 * k + 3)) for k in range(5)]
    n = 3 * len(groups) + 2
    names = [f'x{i}' for i in range(n)]
    rows = []
    for g in groups:
        if (seed % 3) == 0:
            rows.append(([(j, 1) for j in g], '=', 1))
        else:
            rows += [([(j, s) for j in g], '<=', s) for s in (1, -1)]
    chosen = {rng.choice(g) for g in groups} | {n - 2}
    for k, g in enumerate(groups):
        t = [(j, -1 - (j + seed) % 3) for gg in groups[:k+1] for j in gg[:2]]
        t += [(n - 2, -1), (n - 1, -2)]
        rows.append((t, '<=', sum(a for j, a in t if j in chosen) + rng.randrange(3)))
        if k:
            t = [(groups[k-1][0], -1), (g[2], 1)]
            rows.append((t, '<=', sum(a for j, a in t if j in chosen)))
    if seed % 8 in (0, 5):
        t = [(g[1], 1) for g in groups] + [(n - 1, 2)]
        rows += [(t, '<=', 3), (t, '>=', 4)]
    model = Model(('min', 'max')[seed % 2], [rng.randint(-9, 8) for _ in names], rows, names)

    # 无状态压缩：逐项检查上述分组定义域的每个可能赋值。
    points = []
    for take in itertools.product(*groups):
        for a, b in itertools.product((0, 1), repeat=2):
            bits = [0] * n
            for j in take:
                bits[j] = 1
            bits[-2:] = [a, b]
            if model.feasible(bits):
                points.append(model.objective(bits))
    best = (min if model.sense == 'min' else max)(points) if points else None
    return model, best, len(points), 3**len(groups) * 4


def build_direct_probe(directory):
    # 库中只追加一个公开的检验入口，内部识别、搜索、验解和计数逻辑逐字沿用当前源码。
    crate = ROUNDS.parents[1] / 'milp'
    for entry in ['Cargo.toml', 'Cargo.lock']:
        shutil.copy2(crate / entry, directory / entry)
    for entry in ['src', 'examples']:
        shutil.copytree(crate / entry, directory / entry)
    with (directory / 'src/lib.rs').open('a') as fh:
        fh.write('''
pub fn direct_grouped_audit(path: &std::path::Path, limit: u64) {
    let problem = instance::load_instance(path).unwrap();
    let result = grouped_dp::solve(&problem, Default::default(), None, limit, &mut |_| {}).unwrap();
    let (status, bound, stats) = result.expect("该算例应当命中分组动态规划");
    println!("{status:?} {} {}", stats.nodes, bound.map_or("null".into(), |b| b.to_string()));
}
''')
    (directory / 'examples/direct_grouped_audit.rs').write_text('''
fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().unwrap();
    let limit = args.next().unwrap().parse().unwrap();
    milp::direct_grouped_audit(std::path::Path::new(&path), limit);
}
''')
    target = directory / 'target'
    subprocess.run(['cargo', 'build', '--offline', '--locked', '--jobs', '1', '--target-dir', str(target),
                    '--example', 'direct_grouped_audit'], cwd=directory, check=True)
    return target / 'debug/examples/direct_grouped_audit'


def check_run(binary, direct, path, model, truth, limit, expect):
    command = [str(binary), 'timed-solve', str(path), '--threads', '1', '--node-limit', str(limit),
               '--time-limit-secs', '5', '--early-stop', 'off']
    record = json.loads(subprocess.check_output(command, timeout=15))
    status, bound = record['status'], record['bound']
    assert status == expect, (path, status, expect)
    direct_status, nodes, direct_bound = subprocess.check_output(
        [str(direct), str(path), str(limit)], text=True, timeout=15).split()
    assert {'Optimal': 'optimal', 'LimitReached': 'limit_reached',
            'Infeasible': 'infeasible'}[direct_status] == status
    assert int(nodes) == record['nodes'] <= limit and json.loads(direct_bound) == bound
    assert record['lp_iterations'] == 0
    assert status in ('optimal', 'infeasible', 'limit_reached')
    assert status != 'infeasible' or truth is None
    if truth is not None and bound is not None:
        assert (bound <= truth + 1e-6 if model.sense == 'min' else bound >= truth - 1e-6), (truth, bound)
    assert not record.get('callback_errors')
    events = record['incumbent_events']
    times = [event['elapsed_secs'] for event in events]
    assert times == sorted(times) and all(0 <= t <= record['elapsed_secs'] for t in times)
    values = [event['objective'] for event in events]
    assert sorted(values, reverse=model.sense == 'min') == values
    checked = []
    for event in [*record['incumbent_events'], record]:
        if event['values'] is None:
            assert event is record and event['objective'] is None
            continue
        assert len(event['values']) == len(model.names)
        bits = [round(x) for x in event['values']]
        assert all(v in (0, 1) and abs(x-v) <= 1e-6 for x, v in zip(event['values'], bits))
        assert model.feasible(bits), (path.name, '返回向量不满足原约束')
        value = model.objective(bits)
        assert truth is not None and abs(value - event['objective']) <= 1e-6
        assert value >= truth if model.sense == 'min' else value <= truth
        checked.append(''.join(map(str, bits)))
    if status == 'optimal':
        assert values and values[-1] == truth and bound == truth == record['objective']
    elif status == 'limit_reached':
        assert bound is not None
    return dict(status=status, objective=record['objective'], bound=bound, nodes=record['nodes'],
                lp_iterations=record['lp_iterations'], callback_bits=checked[:len(events)],
                final_bits=checked[-1] if record['values'] is not None else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary', type=Path, help='候选版本的 release round_probe')
    parser.add_argument('--output', type=Path)
    opt = parser.parse_args()
    binary = opt.binary.resolve()
    model, optimum, mic, source = native()
    assert optimum == 7615
    sample_counts = Counter()
    results = []
    with TemporaryDirectory(prefix='mic-grouped-audit-') as tmp:
        direct = build_direct_probe(Path(tmp))
        path = Path(tmp) / 'input.milp'
        path.write_bytes(mic)
        real = check_run(binary, direct, path, model, optimum, NODE_CAP, 'optimal')
        selected = [name for name, bit in zip(model.names, real['final_bits'], strict=True) if bit == '1']
        results.append(dict(case='p0201', original_selected_columns=selected, **real))
        for lim in (1, max(1, real['nodes'] - 1), real['nodes']):
            expect = 'optimal' if lim == real['nodes'] else 'limit_reached'
            r = check_run(binary, direct, path, model, optimum, lim, expect)
            assert r['nodes'] <= lim, (r['nodes'], lim)
            results.append(dict(case='p0201-limit', limit=lim, **r))
        # 原模型中的一对相反行严格要求整组之和为 1，整组平移费用只会使每个解等额平移。
        ones = {(tuple(t), rhs) for t, op, rhs in model.rows if op == '<='}
        group = next([j for j, _ in t] for t, op, b in model.rows if op == '<=' and b == 1 and all(a == 1 for _, a in t)
                     and (tuple((j, -a) for j, a in t), -1) in ones)
        for seed, sign, shift in [(103, 1, 0), (412, -1, -503)]:
            variant, order = transform(model, seed, sign, group, shift)
            path.write_text(variant.text())
            r = check_run(binary, direct, path, variant, sign*(optimum+shift), NODE_CAP, 'optimal')
            original = [0]*len(order)
            for old, bit in zip(order, r['final_bits']):
                original[old] = int(bit)
            assert model.feasible(original) and model.objective(original) == optimum
            results.append(dict(case='排列与费用变换', seed=seed, sign=sign, shift=shift,
                                original_selected_columns=[model.names[j] for j, b in enumerate(original) if b], **r))
        for seed in range(24):
            instance, truth, feasible_count, assignments = small(seed)
            permuted, _ = transform(instance, 500 + seed)
            path.write_text(permuted.text())
            full = check_run(binary, direct, path, permuted, truth, NODE_CAP, 'infeasible' if truth is None else 'optimal')
            limited = check_run(binary, direct, path, permuted, truth, 1, 'limit_reached')
            sample_counts['assignments'] += assignments
            sample_counts['feasible_assignments'] += feasible_count
            sample_counts['infeasible' if truth is None else 'optimal'] += 1
            sample_counts[instance.sense + ('_infeasible' if truth is None else '_optimal')] += 1
            sample_counts['one_node_' + limited['status']] += 1
            results.append(dict(case='小型分组', seed=seed, sense=instance.sense, truth=truth,
                                feasible_count=feasible_count, full=full, limited=limited))
    result = dict(source=source['url'], original_mps_sha256=source['mps_sha256'], milp_sha256=source['milp_sha256'],
                  binary_sha256=sha256(binary.read_bytes()).hexdigest(), verified_original_variables=len(model.names),
                  direct_grouped_source_sha256=sha256((ROUNDS.parents[1]/'milp/src/grouped_dp.rs').read_bytes()).hexdigest(),
                  direct_confirmed_runs=sum(2 if 'limited' in row else 1 for row in results),
                  verified_original_rows=len(model.rows), official_optimum=optimum, small=sample_counts, runs=results)
    data = json.dumps(result, ensure_ascii=False, indent=2) + '\n'
    if opt.output:
        opt.output.write_text(data)
    else:
        print(data, end='')


if __name__ == '__main__':
    main()
