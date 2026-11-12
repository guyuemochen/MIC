#!/usr/bin/env python3
"""从原模型复核可行性和目标，交叉检查已证最优值和对偶界。"""
import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import sys


ABS_TOL = 1e-6
REL_TOL = 1e-9


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def tolerance(a, b):
    if a == 0 or b == 0:
        return ABS_TOL
    return ABS_TOL + REL_TOL * max(abs(a), abs(b))


def close(a, b):
    return abs(a - b) <= tolerance(a, b)


def exceeds(a, b):
    return a - b > tolerance(a, b)


def verify_solution(model, row, error):
    values = row['values']
    objective = row['objective']
    status = row['status']
    bound = row['bound']
    for field, value in (('objective', objective), ('bound', bound)):
        if value is not None and not finite(value):
            error(f'{field} 非有限数')
            return None
    if row['auto_degraded'] and (status != 'limit_reached' or not row['gap_terminated']):
        error('自动提前停止状态与 gap 标记不一致')
    if status == 'error':
        error(f'运行失败 {row.get("error", row["raw_status"])}')
    if status == 'optimal':
        if objective is None or bound is None:
            error('最优状态缺目标或界')
        elif not close(objective, bound):
            error(f'最优状态的目标 {objective} 与界 {bound} 不一致')
    if values is None:
        if objective is not None or status == 'optimal':
            error('缺少解向量')
        return None
    if status not in ('optimal', 'limit_reached'):
        error(f'终态 {status} 不应附带可行解')
    variables = model['variables']
    if len(values) != len(variables) or not all(finite(v) for v in values):
        error(f'解向量维度不符或含非有限值（期望 {len(variables)} 项）')
        return None

    for var, x in zip(variables, values):
        lo, hi = var['lb'], var['ub']
        if exceeds(lo, x) or (hi is not None and exceeds(x, hi)):
            error(f'变量 {var["name"]}={x} 超出界 [{lo}, {hi}]')
        if var['kind'] != 'continuous' and abs(x - round(x)) > ABS_TOL:
            error(f'整数变量 {var["name"]}={x} 不满足整数性')
    for i, constraint in enumerate(model['rows']):
        products = [values[j] * coefficient for j, coefficient in constraint['terms']]
        lhs, rhs = math.fsum(products), constraint['rhs']
        allowance = ABS_TOL + REL_TOL * max(abs(rhs), math.fsum(map(abs, products)))
        violation = {'<=': lhs - rhs, '>=': rhs - lhs, '=': abs(lhs - rhs)}[constraint['op']]
        if violation > allowance:
            error(f'行 {i} 约束违反 {violation:.9g}，容差 {allowance:.9g}')
    checked = math.fsum([model['offset']] + [values[j] * c for j, c in model['objective']])
    if objective is None or not close(checked, objective):
        error(f'重算目标 {checked} 与记录目标 {objective} 不符')
    return checked


def verify_bound(model, row, witness, description, error):
    bound = row['bound']
    if bound is None or not finite(bound):
        return
    incorrect = exceeds(bound, witness) if model['sense'] == 'min' else exceeds(witness, bound)
    if incorrect:
        error(f'界 {bound} 排除了{description} {witness}')


def inspect_results(inputs):
    errors = []
    runs = solutions_checked = mic_with_external_optimum = mic_bounds_checked = 0
    models = {}
    groups = defaultdict(list)
    by_instance = {}
    for source in inputs:
        with source.open(encoding='utf-8') as stream:
            for index, line in enumerate(stream, 1):
                row = json.loads(line)
                prefix = f'{source}:{index} {row["instance"]}/{row["variant"]}/{row["rep"]}'

                def error(message, prefix=prefix):
                    errors.append(f'{prefix} {message}')

                path = source.parent / row['model']
                if path not in models:
                    payload = path.read_bytes()
                    models[path] = hashlib.sha256(payload).hexdigest(), json.loads(payload)
                digest, model = models[path]
                runs += 1
                if row['schema'] != 1 or digest != row['model_sha256']:
                    error('结果版本不支持或原模型哈希不符')
                    continue
                if by_instance.setdefault(row['instance'], digest) != digest:
                    error('同名算例的原模型哈希不同')
                before = len(errors)
                objective = verify_solution(model, row, error)
                if objective is not None and len(errors) == before:
                    solutions_checked += 1
                else:
                    objective = None
                groups[(row['instance'], digest)].append((model, row, objective, error))

    for members in groups.values():
        external = [objective for _, row, objective, _ in members
                    if row['solver'] in ('highs', 'scip') and row['status'] == 'optimal'
                    and objective is not None]
        optimal = [objective for _, row, objective, _ in members
                   if row['status'] == 'optimal' and objective is not None]
        feasible = [objective for _, _, objective, _ in members if objective is not None]
        for model, row, objective, error in members:
            if feasible:
                best = min(feasible) if model['sense'] == 'min' else max(feasible)
                if row['status'] == 'infeasible':
                    error('不可行终态与已复核可行解矛盾')
                if row['solver'] != 'mic' or not external:
                    verify_bound(model, row, best, '已复核可行目标', error)
            if row['status'] == 'unbounded' and optimal:
                error('无界终态与已证最优状态矛盾')
            if row['status'] == 'optimal' and objective is not None and optimal:
                if any(not close(objective, other) for other in optimal):
                    error(f'已证最优目标与同模型其他运行不一致 {optimal}')
            if row['solver'] == 'mic' and external:
                mic_with_external_optimum += 1
                if finite(row['bound']):
                    mic_bounds_checked += 1
                for optimum in set(external):
                    verify_bound(model, row, optimum, '外部已证最优值', error)
                    if objective is not None:
                        improves = (exceeds(optimum, objective) if model['sense'] == 'min'
                                    else exceeds(objective, optimum))
                        if improves:
                            error(f'MIC 可行目标 {objective} 优于外部已证最优值 {optimum}')
    return dict(passed=not errors and runs > 0, runs=runs, instances=len(groups),
                solutions_checked=solutions_checked,
                mic_with_external_optimum=mic_with_external_optimum,
                mic_bounds_checked=mic_bounds_checked,
                errors=errors)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', nargs='+', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = inspect_results(args.inputs)
    rendered = json.dumps(result, ensure_ascii=False, indent=2) + '\n'
    if args.output:
        with args.output.open('x', encoding='utf-8') as target:
            target.write(rendered)
    print(rendered, end='')
    return not result['passed']


if __name__ == '__main__':
    sys.exit(main())
