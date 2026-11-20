#!/usr/bin/env python3
"""从冻结的 MILP 原文核对布尔编码，并逐步重放不可行证明。运行 python3 verify.py。"""

import argparse
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
import gzip
from hashlib import sha256
from itertools import product
import json
from math import lcm
import operator
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory


HERE = Path(__file__).resolve().parent
FROZEN = HERE.parent / 'all-gurobi10-remaining-60s'
OPS = {'<=': operator.le, '=': operator.eq, '>=': operator.ge}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def uncompressed(path):
    return gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes()


def digest(data):
    return sha256(data).hexdigest()


@lru_cache(None)
def row_truth(op, coefficients, rhs):
    def integer_row(a, b):
        scale = lcm(b.denominator, *(v.denominator for v in a))
        return OPS[op], tuple(int(v * scale) for v in a), int(b * scale)

    text_row = integer_row(coefficients, rhs)
    # 原文写出的十进制与 JSON 的最短表示可不同。核对加载后的二进制浮点
    # 数，再按两套系数的精确分数逐一核对原行上所有可能的二元取值。
    binary_row = integer_row(tuple(Fraction(float(v)) for v in coefficients), Fraction(float(rhs)))
    for bits in product((0, 1), repeat=len(coefficients)):
        require(satisfies([text_row], bits) == satisfies([binary_row], bits),
                '十进制原问题与求解器二进制系数具有不同真值')
    return text_row


def original_network():
    manifest = json.loads((FROZEN / 'manifest.json').read_bytes())
    case, = [c for c in manifest['cases'] if c['instance'] == 'bnatt500']
    text = uncompressed(FROZEN / case['input'])
    require(digest(text) == case['milp_sha256'], '原始 MILP 与冻结摘要不符')
    exported = uncompressed(FROZEN / case['model'])
    require(digest(exported) == case['model_sha256'], '导出模型与冻结摘要不符')
    m = json.loads(exported, parse_float=Fraction)
    tokens = [line.split() for line in text.decode().splitlines() if line.strip()]
    require(tokens[:3] == [['name', 'bnatt500'], ['sense', 'min'], ['vars']], '实例头不符')
    obj_at = tokens.index(['obj'])
    cons_at = tokens.index(['cons'])
    box_at = tokens.index(['bounds'])
    require(tokens[-1] == ['end'], '缺少实例结束标记')
    names = []
    for name, kind in tokens[3:obj_at]:
        require(kind == 'binary', '出现非二元变量')
        names.append(name)
    ids = {name: j for j, name in enumerate(names)}
    require(len(ids) == len(names) == len(m['variables']), '变量重复或导出缺失')
    for var, name in zip(m['variables'], names):
        require(var == dict(name=name, kind='binary', lb=0, ub=1), '导出变量界不同')
    # 本冻结文件将每个变量的两端界都显式列出；改变界会改变本证明的适用输入。
    limits = {(name, op, Fraction(v)) for name, op, v in tokens[box_at + 1:-1]}
    require(limits == {(name, op, value) for name in names for op, value in (('>=', 0), ('<=', 1))},
            '二元变量上下界发生变化')
    require(m['offset'] == 0 and m['sense'] == 'min', '导出目标不符')
    cost, variable = tokens[obj_at + 1]
    require(obj_at + 2 == cons_at and m['objective'] == [[ids[variable], Fraction(cost)]],
            '导出目标与原文不符')

    groups = defaultdict(list)
    raw_rows = tokens[cons_at + 1:box_at]
    require(len(raw_rows) == len(m['rows']) == case['rows'], '原行数量不符')
    for raw, exported_row in zip(raw_rows, m['rows']):
        label, *coeffs, op, bound = raw
        require(label.endswith(':') and len(coeffs) % 2 == 0 and op in OPS, '无效的原始行')
        terms = sorted((ids[name], Fraction(a))
                       for name, a in zip(coeffs[1::2], coeffs[::2]))
        cols = tuple(j for j, _ in terms)
        rhs = Fraction(bound)
        require(len(set(cols)) == len(cols), '冻结行出现重复变量')
        require([(j, float(a)) for j, a in exported_row['terms']] == [(j, float(a)) for j, a in terms]
                and float(exported_row['rhs']) == float(rhs) and exported_row['op'] == op,
                '原行与导出模型的二进制系数不同')
        groups[cols].append(row_truth(op, tuple(a for _, a in terms), rhs))
    return len(ids), groups, case


def read_cnf(data):
    clauses, mapping, header = [], {}, None
    for line in data.decode().splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[:2] == ['c', 'var']:
            _, _, literal, column = fields
            k, j = int(literal), int(column)
            require(k not in mapping, 'DIMACS 变量映射重复')
            mapping[k] = j
        elif fields[0] == 'c':
            continue
        elif fields[0] == 'p':
            require(header is None and len(fields) == 4 and fields[1] == 'cnf', 'DIMACS 头不符')
            header = tuple(map(int, fields[2:]))
        else:
            require(fields[-1] == '0' and header is not None, '子句缺少头或终止符')
            c = tuple(map(int, fields[:-1]))
            require(all(1 <= abs(lit) <= header[0] for lit in c), '文字越界')
            clauses.append(c)
    require(header == (len(mapping), len(clauses)) and set(mapping) == set(range(1, len(mapping) + 1)),
            'DIMACS 数量或映射不符')
    require(len(set(mapping.values())) == len(mapping), '多个布尔变量对应同一原列')
    return mapping, clauses


def satisfies(rows, values):
    return all(op(sum(a * x for a, x in zip(coeff, values)), rhs) for op, coeff, rhs in rows)


def verify_projection(n, groups, mapping, clauses):
    # 以生产映射指定的主列作为自变量。用原始十进制系数的精确真值判定恢复
    # 私有列，避免复用生产端的门型识别、浮点比较或化简方法。
    to_lit = {j: k for k, j in mapping.items()}
    main = set(to_lit)
    fixed, auxiliary, outputs = {}, {}, {}
    for cols, rows in groups.items():
        if len(cols) == 1:
            j, = cols
            allowed = [v for v in (0, 1) if satisfies(rows, (v,))]
            require(len(allowed) == 1, '一元行不唯一确定变量')
            fixed[j], = allowed
        elif len(cols) == 4:
            aux, = set(cols) - main
            require(aux not in auxiliary, '同一私有列有多个定义')
            inputs = tuple(j for j in cols if j != aux)
            table = {}
            for bits in product((0, 1), repeat=len(inputs)):
                values = dict(zip(inputs, bits))
                possible = [v for v in (0, 1) if satisfies(rows, [v if j == aux else values[j] for j in cols])]
                require(len(possible) == 1, '私有列不能唯一还原')
                table[bits], = possible
            auxiliary[aux] = inputs, table
        elif len(cols) == 9:
            output, = set(cols) & main
            require(output not in outputs, '同一输出被多个主方程约束')
            outputs[output] = cols, rows
        else:
            raise ValueError('冻结输入出现未知类型的行')
    require(main == set(outputs) and len(main) + len(fixed) + len(auxiliary) == n
            and main | set(fixed) | set(auxiliary) == set(range(n)), '原变量没有得到互斥且完整的覆盖')

    # 第一方向：每个生产子句至少被一个原输出关系蕴含。第二方向：每个
    # 违反原关系的局部赋值，至少被一个作用域内的生产子句排除。
    # 多个原输出的作用域可以重合，同一子句可用于多个局部证明。
    postings = defaultdict(set)
    for i, c in enumerate(clauses):
        for lit in c:
            postings[abs(lit)].add(i)
    entailed = set()
    checked = rejected = 0
    for output, (cols, rows) in outputs.items():
        deps = {output}
        for j in cols:
            if j in auxiliary:
                deps.update(auxiliary[j][0])
        scope = tuple(sorted(deps))
        positions = {to_lit[j]: i for i, j in enumerate(scope)}
        related = set().union(*(postings[k] for k in positions))
        local = {i: [(positions[abs(lit)], lit > 0) for lit in clauses[i]]
                 for i in related if all(abs(lit) in positions for lit in clauses[i])}
        excluded_by_valid = set()
        for bits in product((0, 1), repeat=len(scope)):
            values = dict(zip(scope, bits))
            for j in cols:
                if j in auxiliary:
                    inputs, table = auxiliary[j]
                    values[j] = table[tuple(values[x] for x in inputs)]
                elif j in fixed:
                    values[j] = fixed[j]
            legal = satisfies(rows, [values[j] for j in cols])
            bad = {i for i, c in local.items() if not any(bits[k] == v for k, v in c)}
            if legal:
                excluded_by_valid.update(bad)
            else:
                require(bad, f'CNF 放过违反原输出列 {output} 的赋值 {bits}')
                rejected += 1
            checked += 1
        entailed.update(set(local) - excluded_by_valid)
    require(len(entailed) == len(clauses), '存在未由原模型蕴含的 CNF 子句')
    return dict(primary_variables=len(main), fixed_variables=len(fixed), restored_variables=len(auxiliary),
                private_input_tuples=sum(len(table) for _, table in auxiliary.values()),
                checked_output_tuples=checked, forbidden_output_tuples=rejected, entailed_clauses=len(entailed))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cnf', type=Path, default=HERE / 'network.cnf.gz')
    parser.add_argument('--proof', type=Path, action='append', help='全式或互补两个假设的单位传播证明')
    args = parser.parse_args()
    cnf = uncompressed(args.cnf)
    mapping, clauses = read_cnf(cnf)
    n, groups, source = original_network()
    report = verify_projection(n, groups, mapping, clauses)
    proofs = args.proof or [HERE / 'full.rup.gz']
    reports = []
    with TemporaryDirectory() as tmp:
        base = Path(tmp)
        formula = base / 'formula.cnf'
        formula.write_bytes(cnf)
        exe = base / 'verify_rup'
        subprocess.run(['rustc', '-O', '--edition=2021', str(HERE / 'verify_rup.rs'), '-o', str(exe)],
                       check=True)
        for proof in proofs:
            data = uncompressed(proof)
            assumptions = [int(line.split()[2]) for line in data.splitlines()
                           if line.split()[:2] == [b'c', b'assumption'] and len(line.split()) == 3]
            require(len(assumptions) == 1, '每条证明必须声明一个假设（0 表示无）')
            assumption, = assumptions
            path = base / 'steps.rup'
            path.write_bytes(data)
            out = subprocess.check_output([exe, str(len(mapping)), formula, path, str(assumption)], text=True)
            initial, steps = map(int, out.split())
            require(initial == len(clauses), '原始子句计数不一致')
            reports.append(dict(proof=proof.name, proof_sha256=digest(data), assumption=assumption,
                                verified_steps_including_empty=steps))
    assumptions = [r['assumption'] for r in reports]
    require(assumptions == [0] or len(assumptions) == 2 and assumptions[0] and sum(assumptions) == 0,
            '证明仅覆盖了部分取值')
    report.update(milp_sha256=source['milp_sha256'], model_sha256=source['model_sha256'],
                  original_rows=source['rows'], original_variables=n, cnf_sha256=digest(cnf),
                  proofs=reports, conclusion='infeasible')
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
