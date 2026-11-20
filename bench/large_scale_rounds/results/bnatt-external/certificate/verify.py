#!/usr/bin/env python3
"""从 MIPLIB 的 bnatt400 原始系数重放最优性证明，仅需 Python 标准库与 rustc。"""

import argparse
from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path
import runpy
import subprocess
from tempfile import TemporaryDirectory


HERE = Path(__file__).resolve().parent
NATIVE = HERE.parent / 'native'
PREVIOUS = HERE.parents[1] / 'bnatt-boolean'
shared = runpy.run_path(str(PREVIOUS / 'verify.py'))
require = shared['require']
uncompressed = shared['uncompressed']
digest = shared['digest']
row_truth = shared['row_truth']
OPS = shared['OPS']
MPS_SHA256 = '50565a2a2e80dd9563ce7b5d06e3053d986acf72d57ab787d97d1d10ad519ef8'


def read_mps(raw):
    """按 MPS 原始顺序保留行、列与十进制系数。这里只接收该公开输入使用的语法。"""
    require(digest(raw) == MPS_SHA256, '官方 bnatt400 的未压缩 MPS 摘要不符')
    names, kinds, terms, rhs, upper = {}, {}, defaultdict(dict), {}, {}
    sections = []
    integer = False
    for line in raw.decode('ascii').splitlines():
        fields = line.split()
        if not fields or line.startswith('*'):
            continue
        if not line[0].isspace():
            sections.append(fields[0])
            continue
        section = sections[-1]
        if section == 'ROWS':
            kind, row = fields
            require(row not in kinds and kind in ('N', 'E', 'G', 'L'), 'MPS 行类型或名称不符')
            kinds[row] = kind
        elif section == 'COLUMNS':
            name, row, value = fields
            if row == "'MARKER'":
                require(value == ("'INTEND'" if integer else "'INTORG'"), '整数标记异常')
                integer = not integer
            else:
                require(integer and row in kinds, '原列非整数或引用不存在的行')
                j = names.setdefault(name, len(names))
                require(j not in terms[row], '原行含重复的列')
                terms[row][j] = Fraction(value)
        elif section == 'RHS':
            vector, row, value = fields
            require(vector == 'rhs' and row in kinds and row not in rhs, '右端项异常')
            rhs[row] = Fraction(value)
        elif section == 'BOUNDS':
            kind, vector, name, value = fields
            require(kind == 'UP' and vector == 'bnd' and name not in upper, '变量界异常')
            upper[name] = Fraction(value)
        else:
            raise ValueError(f'无法解释 {section} 段中的 {line}')
    # 无 OBJSENSE 时 MPS 默认最小化。整数变量的默认下界为零，此文件逐列给出上界一。
    require(sections == ['NAME', 'ROWS', 'COLUMNS', 'RHS', 'BOUNDS', 'ENDATA'] and not integer,
            'MPS 段不完整或存在未经验证的模型语义')
    require(set(upper) == set(names) and all(b == 1 for b in upper.values()), '并非所有列都是二元变量')
    objective, = [r for r in kinds if kinds[r] == 'N']
    require(rhs.get(objective, 0) == 0 and terms[objective] == {names['x1']: 1}, '原目标不是 min x1')
    senses = {'E': '=', 'G': '>=', 'L': '<='}
    rows = [(r, senses[k], sorted(terms[r].items()), rhs.get(r, Fraction(0)))
            for r, k in kinds.items() if r != objective]
    return names, rows


def verify_milp(raw, names, rows):
    """把实际 MIC 输入逐行与 MPS 对照，建立 x编号 到 MPS 原列的映射。"""
    lines = iter(line.split() for line in raw.decode().splitlines() if line.strip())
    require([next(lines), next(lines), next(lines)] == [['name', 'bnatt400'], ['sense', 'min'], ['vars']],
            'MIC 模型头不符')
    for j in range(len(names)):
        require(next(lines) == [f'x{j}', 'binary'], f'MIC 第 {j} 列的名称或类型不符')
    require([next(lines), next(lines), next(lines)] == [['obj'], ['1', f'x{names["x1"]}'], ['cons']],
            'MIC 模型目标发生变化')
    groups = defaultdict(list)
    for i, (name, op, terms, bound) in enumerate(rows):
        suffix = {'>=': 'l', '<=': 'u', '=': ''}[op]
        label, *body, sense, rhs = next(lines)
        require(label == f'r{i}{suffix}:' and sense == op and float(rhs) == float(bound),
                f'MIC 与 MPS 的第 {name} 行不一致')
        require(len(body) == 2 * len(terms) and
                all(x == f'x{j}' and float(a) == float(c)
                    for x, a, (j, c) in zip(body[1::2], body[::2], terms)), f'第 {name} 行的列或系数变化')
        # row_truth 另逐一对比原十进制与实际读取的二进制系数在该行所有二元赋值上的真值。
        groups[tuple(j for j, _ in terms)].append(row_truth(op, tuple(c for _, c in terms), bound))
    require(next(lines) == ['bounds'], 'MIC 输入缺少变量界')
    for j in range(len(names)):
        require([next(lines), next(lines)] == [[f'x{j}', '>=', '0'], [f'x{j}', '<=', '1']], 'MIC 变量界异常')
    require(list(lines) == [['end']], 'MIC 模型末尾不符')
    return groups


def verify_witness(raw, names, rows):
    result = json.loads(raw, parse_float=Fraction)
    values = result['values']
    require(len(values) == len(names) and all(v in (0, 1) for v in values), '见证长度或二元取值无效')
    bits = list(map(int, values))
    for name, op, terms, bound in rows:
        activity = sum(c * bits[j] for j, c in terms)
        require(OPS[op](activity, bound), f'见证违反官方 MPS 的第 {name} 行')
    value = bits[names['x1']]
    require(result['status'] == 'optimal' and result['objective'] == value, 'MIC 状态或报告值与见证不符')
    return value, bits


def verify_better(cnf, proof, mapping, clauses, j):
    # 目标仅有一个二元列，见证为 1 时，唯一更优取值为零，即一个负文字。
    lit, = [k for k, original in mapping.items() if original == j]
    assumptions = [line.split() for line in proof.splitlines() if line.split()[:2] == [b'c', b'assumption']]
    require(assumptions == [[b'c', b'assumption', str(-lit).encode()]], 'RUP 假设没有覆盖原目标列等于零')
    with TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        formula = tmp / 'formula.cnf'
        formula.write_bytes(cnf)
        steps = tmp / 'steps.rup'
        steps.write_bytes(proof)
        checker = tmp / 'verify_rup'
        subprocess.run(['rustc', '-O', '--edition=2021', PREVIOUS / 'verify_rup.rs', '-o', checker], check=True)
        output = subprocess.check_output([checker, str(len(mapping)), formula, steps, str(-lit)], text=True)
    initial, checked = map(int, output.split())
    require(initial == len(clauses) and checked > 0, 'RUP 校验未覆盖完整输入或缺少空子句')
    return dict(proof_sha256=digest(proof), assumption=-lit, verified_steps_including_empty=checked)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mps', type=Path, default=NATIVE / 'bnatt400.mps.gz')
    parser.add_argument('--milp', type=Path, default=NATIVE / 'bnatt400.milp.gz')
    parser.add_argument('--result', type=Path, default=NATIVE / 'current.jsonl.gz')
    parser.add_argument('--cnf', type=Path, default=NATIVE / 'network.cnf.gz')
    parser.add_argument('--proof', type=Path, default=NATIVE / 'better.rup.gz')
    args = parser.parse_args()
    mps, milp, result = map(uncompressed, (args.mps, args.milp, args.result))
    names, rows = read_mps(mps)
    groups = verify_milp(milp, names, rows)
    optimum, bits = verify_witness(result, names, rows)
    report = dict(instance='bnatt400', mps_sha256=digest(mps), milp_sha256=digest(milp),
                  witness_bits_sha256=digest(bytes(bits)), original_variables=len(names), original_rows=len(rows),
                  mps_objective='min x1', mic_objective_column=names['x1'], verified_objective=optimum)
    if optimum == 1:
        cnf = uncompressed(args.cnf)
        mapping, clauses = shared['read_cnf'](cnf)
        report['projection'] = shared['verify_projection'](len(names), groups, mapping, clauses)
        require(all(any(bool(bits[mapping[abs(lit)]]) == (lit > 0) for lit in c) for c in clauses),
                '原可行见证违反某条 CNF 子句')
        report['cnf_sha256'] = digest(cnf)
        report['better_value_excluded'] = verify_better(cnf, uncompressed(args.proof), mapping, clauses, names['x1'])
    else:
        report['lower_bound'] = '原目标变量的下界为零'
    report['conclusion'] = 'optimal'
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
