#!/usr/bin/env python3
"""核验真实规模变体的原始向量。可用 --probe 指定 round_probe 再从正式入口复跑。"""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import random
import subprocess
from tempfile import TemporaryDirectory

from verify import HERE, FROZEN, digest, original_network, read_cnf, require, row_truth, satisfies, uncompressed


def check_case(case, source, n, core, temp, probe):
    start = source.index
    variables = source[start('vars') + 1:start('obj')]
    rows = source[start('cons') + 1:start('bounds')]
    bounds = source[start('bounds'):]
    target = case.get('objective_variable', 0)
    changed = case.get('original_rows', [])
    if changed:
        for before, after in zip(changed, case['replacement_rows'], strict=True):
            require(rows.count(before) == 1, '找不到待替换的唯一原行')
            rows[rows.index(before)] = after
    elif 'permutation_seed' in case:
        rng = random.Random(case['permutation_seed'])
        rng.shuffle(variables)
        rng.shuffle(rows)
    else:
        rows.reverse()
    data = ('\n'.join([f"name {case['name']}", f"sense {case['sense']}", 'vars', *variables,
                       'obj', f'1 x{target}', 'cons', *rows, *bounds]) + '\n').encode()
    require(digest(data) == case['raw_sha256'], '变体输入摘要不同')
    constraints = []
    for line in rows:
        _, *terms, op, bound = line.split()
        cols = [int(tok[1:]) for tok in terms[1::2]]
        predicate, coeff, rhs = row_truth(op, tuple(map(Fraction, terms[::2])), Fraction(bound))
        constraints.append((cols, predicate, coeff, rhs))

    def check_point(x, obj):
        require(len(x) == n and all(v in (0, 1) for v in x) and x[target] == obj, '变量或目标不同')
        require(all(op(sum(a * x[j] for j, a in zip(cols, coeff)), rhs)
                    for cols, op, coeff, rhs in constraints), '点违反变体中的原行')

    if changed:
        # 当可行见证已达到二元目标的天然端点，就足以独立确定最优值。
        optimum = 0 if case['sense'] == 'min' else 1
        require(case['status'] == 'optimal' and case['objective'] == optimum == case['bound'],
                '证优记录未闭合')
        assignment = dict(zip(core, map(int, case['planted_primary_bits']), strict=True))
        for cols, op, coeff, rhs in constraints:
            if len(cols) == 1:
                j, = cols
                choices = [b for b in (0, 1) if op(coeff[0] * b, rhs)]
                assignment[j], = choices
        pairs = {}
        for cols, op, coeff, rhs in constraints:
            if len(cols) == 4:
                pairs.setdefault(tuple(cols), []).append((op, coeff, rhs))
        for cols, pair in pairs.items():
            j, = set(cols) - set(core)
            assignment[j], = [b for b in (0, 1) if satisfies(pair, [b if k == j else assignment[k] for k in cols])]
        check_point([assignment[j] for j in range(n)], optimum)
        require(len(case['original_order_vectors']) == 1 + case['checked_callbacks'], '记录缺少原向量')
    else:
        require(case['status'] == 'infeasible' and case['objective'] is None and not case['original_order_vectors'],
                '等价重排的状态或向量错误')
    for point in case['original_order_vectors']:
        packed = bytes.fromhex(point['bits_le'])
        require(len(packed) == (n + 7) // 8 and int.from_bytes(packed, 'little') >> n == 0,
                '向量长度或高位错误')
        x = [(packed[j // 8] >> (j % 8)) & 1 for j in range(n)]
        check_point(x, point['objective'])
    if probe:
        path = temp / (case['name'] + '.milp')
        path.write_bytes(data)
        run = json.loads(subprocess.check_output([probe, 'timed-solve', path, '--threads', '1',
                         '--time-limit-secs', '60', '--node-limit', '2000000', '--early-stop', 'off']))
        require(all(run[k] == case[k] for k in ('status', 'objective', 'bound', 'nodes', 'lp_iterations')),
                '正式入口运行结果与冻结变体不同')
        order = [int(v.split()[0][1:]) for v in variables]
        for pt in [*run['incumbent_events'], run]:
            if pt['objective'] is not None:
                raw = dict(zip(order, pt['values'], strict=True))
                check_point([raw[j] for j in range(n)], pt['objective'])
    return dict(name=case['name'], status=case['status'], checked_events=case['checked_callbacks'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path)
    args = parser.parse_args()
    n, _, manifest = original_network()
    archive = json.loads((HERE / 'real-variants.json').read_bytes())
    require(manifest['milp_sha256'] == archive['source_sha256'], '冻结原模型不同')
    source = uncompressed(FROZEN / manifest['input']).decode().splitlines()
    mapping, _ = read_cnf(uncompressed(HERE / 'network.cnf.gz'))
    # 归档见证按 DIMACS 升序；原列顺序只用于回填最终点。
    core = [mapping[i] for i in range(1, len(mapping) + 1)]
    with TemporaryDirectory() as tmp:
        for case in archive['cases']:
            print(json.dumps(check_case(case, source, n, core, Path(tmp), args.probe), ensure_ascii=False))


if __name__ == '__main__':
    main()
