#!/usr/bin/env python3
"""不运行求解器，从原始输入独立复核 MIC 与 Gurobi 10 的单遍计时记录。"""

import argparse
from collections import Counter, defaultdict
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from statistics import geometric_mean, median

from check_results import close, exceeds, finite, load_public_truth, verify_bound, verify_solution


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / 'results'
SOLVERS = ('mic', 'gurobi')
SOLVER_NAMES = {'mic': 'MIC', 'gurobi': 'Gurobi 10'}
CONCLUSIONS = {'optimal': '证优', 'infeasible': '不可行', 'unbounded': '无界'}
STATUSES = set(CONCLUSIONS) | {'limit_reached', 'infeasible_or_unbounded', 'error'}
STATUS_TEXT = dict(CONCLUSIONS, limit_reached='触限', infeasible_or_unbounded='不可行或无界', error='失败')
HASH_KEYS = ('model_sha256', 'milp_sha256', 'mps_sha256')
RUN_KEYS = ('threads', 'time_limit_secs', 'node_limit')
MODEL_KEYS = ('name', 'sense', 'offset', 'variables', 'objective', 'rows')
ROW_KEYS = ('schema', 'case_id', 'instance', 'solver', 'status', 'raw_status', 'rep', 'order',
            'model', 'objective', 'bound', 'values', 'elapsed_secs', 'version', *HASH_KEYS[:2], *RUN_KEYS)
MICROSECONDS = 0.001


def read_bytes(path):
    return gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes()


def parse_milp(payload):
    """按 MILP 文本重建导出字段，逐行对照双方使用的模型。"""
    lines = [tokens for line in payload.decode().splitlines()
             if (tokens := line.partition('#')[0].split())]
    name = lines.pop(0)[1] if lines[0][0] == 'name' else ''
    keyword, sense = lines.pop(0)
    if keyword != 'sense' or sense not in ('min', 'max'):
        raise ValueError('原始 .milp 缺少目标方向')
    sections = {}
    for tokens in lines:
        if tokens[0] in ('vars', 'obj', 'cons', 'bounds', 'end'):
            section, = tokens
            if section in sections:
                raise ValueError(f'.milp 节重复 {section}')
            sections[section] = []
        else:
            sections[section].append(tokens)
    if list(sections) not in (['vars', 'obj', 'cons', 'end'], ['vars', 'obj', 'cons', 'bounds', 'end']) or sections['end']:
        raise ValueError('.milp 节顺序错误或 end 后存在内容')

    def number(token):
        value = float(token)
        if not finite(value):
            raise ValueError(f'.milp 非有限系数 {token}')
        return value

    variables = [dict(name=key, kind=kind, lb=0.0, ub=1.0 if kind == 'binary' else None)
                 for key, kind in sections['vars']]
    columns = {var['name']: j for j, var in enumerate(variables)}
    if len(columns) != len(variables) or any(
        var['kind'] not in ('continuous', 'integer', 'binary') for var in variables
    ):
        raise ValueError('.milp 变量重名或类型无效')
    for key, op, token in sections.get('bounds', []):
        var, bound = variables[columns[key]], number(token)
        if op not in ('<=', '>=', '=') or (var['kind'] == 'binary' and not 0 <= bound <= 1):
            raise ValueError(f'.milp 变量界无效 {key} {op} {bound}')
        if op in ('>=', '='):
            var['lb'] = max(var['lb'], bound)
        if op in ('<=', '='):
            var['ub'] = bound if var['ub'] is None else min(var['ub'], bound)
        if var['ub'] is not None and var['lb'] > var['ub']:
            raise ValueError(f'.milp 变量界为空 {key}')

    def terms(tokens):
        entries = {columns[key]: number(token) for token, key in zip(tokens[::2], tokens[1::2], strict=True)}
        if len(entries) * 2 != len(tokens):
            raise ValueError('.milp 线性表达式内变量重复')
        return [[j, value] for j, value in sorted(entries.items()) if value != 0]

    rows = []
    for tokens in sections['cons']:
        body = tokens[1:] if tokens[0].endswith(':') else tokens
        *coefficients, op, rhs = body
        if op not in ('<=', '>=', '='):
            raise ValueError(f'.milp 约束运算符无效 {op}')
        rows.append(dict(op=op, rhs=number(rhs), terms=terms(coefficients)))
    objective = terms([token for line in sections['obj'] for token in line])
    return dict(name=name, sense=sense, offset=0.0, variables=variables, objective=objective, rows=rows)


def compare_model(model, original, error):
    for key in MODEL_KEYS:
        if key not in ('rows', 'variables'):
            if model[key] != original[key]:
                error(f'原始 .milp 的 {key} 与导出不同')
            continue
        if len(model[key]) != len(original[key]):
            error(f'原始 .milp 的 {key} 项数与导出不同')
            continue
        for index, (exported, parsed) in enumerate(zip(model[key], original[key]), 1):
            if exported != parsed:
                error(f'原始 .milp 的 {key} 第 {index} 项与导出不同')
                break


def public_truths():
    return [load_public_truth(RESULTS / filename)[0]
            for filename in ('public-sources.json', 'miplib2017-sources.json')]


def case_truth(base, case, public, error):
    name = case['instance']
    if any(name in entries for entries in public):
        legacy, modern = public
        entry = modern.get(name, legacy.get(name))
        for key in ('milp_sha256', 'mps_sha256'):
            if not entry.get(key) or case.get(key) != entry[key]:
                error(f'公开来源清单的 {key} 与冻结值不同')
                return None
        # MIPLIB 3.0 将官方目标直接写在原 MPS 的 BEST SOLN 注释中。
        source = ROOT / case['source_path']
        path = base / case['mps'] if case.get('mps') else source.with_suffix('.mps')
        if path.exists():
            payload = read_bytes(path)
            if sha256(payload).hexdigest() != entry['mps_sha256']:
                error('公开原始 MPS 的 SHA-256 与冻结值不同')
                return None
            if name in legacy:
                match = re.search(rb'(?m)^\*BEST SOLN:\s*(\S+)\s+\(opt\)', payload[:4096])
                if match is None or not finite(value := float(match[1])):
                    error('MIPLIB 3.0 原始 MPS 未给出有效的官方最优值')
                    return None
                return dict(source='MIPLIB 3.0 原始 MPS', url=entry['url'], status='optimal',
                            objective=value, mps_rehashed=True)
        elif name in legacy or case.get('mps'):
            error(f'无法读取已冻结的原始 MPS {path}')
            return None
        return dict(source='MIPLIB 2017 v37', status=entry['solution_status'],
                    objective=entry.get('solution_value'), mps_rehashed=path.exists())

    source = ROOT / case['source_path']
    companion = source.with_suffix('.opt')
    if companion.exists():
        payload = companion.read_bytes()
        value = float(payload)
        if not finite(value):
            error(f'{companion} 的动态规划最优值无效')
            return None
        return dict(source='动态规划 .opt', status='optimal', objective=value, sha256=sha256(payload).hexdigest())
    if name in ('solve-1000', 'solve-8000', 'solve-100000'):
        return dict(source='容量 LP 构造的已知最优值', status='optimal', objective=0.0)
    return None


def load_case(base, case, public, error):
    payload = read_bytes(base / case['model'])
    if sha256(payload).hexdigest() != case['model_sha256']:
        error('导出模型的 SHA-256 与冻结值不同')
    model = json.loads(payload)
    original = read_bytes(base / case['input'])
    if sha256(original).hexdigest() != case['milp_sha256']:
        error('原始 .milp 的 SHA-256 与冻结值不同')
    compare_model(model, parse_milp(original), error)
    return model, case_truth(base, case, public, error)


def valid_time(value):
    return finite(value) and value >= 0


def read_events(model, row, error):
    first, goals, rejected = None, [], []
    previous = 0.0
    elapsed = row['elapsed_secs']
    for index, event in enumerate(row.get('incumbent_events', []), 1):
        seen = event.get('elapsed_secs')
        if not valid_time(elapsed) or not valid_time(seen) or not previous <= seen <= elapsed:
            error(f'第 {index} 个回调的时间无效或顺序错误')
            continue
        previous = seen
        issues = []
        claims = dict(status='limit_reached', bound=None, auto_degraded=False,
                      objective=event.get('objective'), values=event.get('values'))
        checked = verify_solution(model, claims, issues.append)
        if checked is None or issues:
            rejected.append(dict(event=index, elapsed_secs=seen, issues=issues[:5] or ['缺少解或目标'],
                                 violations=len(issues)))
            continue
        goals.append(checked)
        if first is None:
            first = dict(event=index, elapsed_secs=seen, objective=checked)
    return first, goals, rejected


def better(model, a, b):
    return exceeds(b, a) if model['sense'] == 'min' else exceeds(a, b)


def audit_truth(model, truth, row, witnesses, error):
    if not truth:
        return
    expected, optimum = truth['status'], truth['objective']
    status = row['status']
    if expected == 'infeasible' and (status in ('optimal', 'unbounded') or witnesses):
        error('终态或可行解与官方不可行结论矛盾')
    if expected == 'unbounded' and status in ('infeasible', 'optimal'):
        error('终态与官方无界结论矛盾')
    if expected == 'best_known' and status == 'infeasible':
        error('不可行终态与官方已知可行目标矛盾')
    if expected in ('optimal', 'best_known'):
        verify_bound(model, row, optimum, '外部参考目标', error)
    if expected == 'optimal':
        if status in ('infeasible', 'unbounded', 'infeasible_or_unbounded'):
            error(f'终态与参考已证最优值 {optimum} 矛盾')
        if status == 'optimal' and (not finite(row['objective']) or not close(row['objective'], optimum)):
            error(f'最优目标 {row["objective"]} 与参考已证最优值 {optimum} 不符')
        if any(better(model, value, optimum) for value in witnesses):
            error(f'已复核的可行目标优于参考已证最优值 {optimum}')


def audit_run(model, kind, truth, freeze, case, raw, issues):
    error = issues.append
    missing = [key for key in ROW_KEYS if key not in raw]
    if missing:
        error(f'运行记录缺少字段 {", ".join(missing)}')
    row = dict.fromkeys(ROW_KEYS, None) | raw
    for key in ('case_id', 'instance', 'suite', 'model', 'input', *HASH_KEYS):
        if row.get(key) != case.get(key):
            error(f'运行记录与冻结清单的 {key} 不符')
    for key in RUN_KEYS:
        if row[key] != freeze.get(key):
            error(f'运行记录与冻结参数 {key} 不符')
    if row['schema'] != 1 or row['rep'] != 0 or row['order'] not in (0, 1):
        error('运行的版本、遍数或求解器顺序无效')
    solver, status, elapsed = row['solver'], row['status'], row['elapsed_secs']
    if status not in STATUSES:
        error(f'未知终态 {status!r}')
    if solver == 'gurobi' and row['version'] != freeze.get('gurobi_version'):
        error('运行记录中的 Gurobi 版本与冻结版本不符')
    if solver == 'mic' and row['version'] not in (freeze.get('source_commit'), freeze.get('mic_binary_sha256')):
        error('运行记录中的 MIC 版本与冻结提交和二进制均不符')
    if status != 'error' and not valid_time(elapsed):
        error('求解器正常返回却没有有效的求解时间')
    if 'proof_elapsed_secs' in row and row['proof_elapsed_secs'] != (elapsed if status == 'optimal' else None):
        error('原始证优时间与已记录的求解终态或耗时不同')
    if solver == 'gurobi' and status != 'error':
        if row.get('model_roundtrip_passed') is not True:
            error('缺少 Gurobi 模型逐项读回通过记录')
        for key, value in {'Threads': freeze['threads'], 'TimeLimit': freeze['time_limit_secs'],
                           'NodeLimit': freeze['node_limit'], 'MIPGap': freeze['mip_gap'],
                           'MIPGapAbs': freeze['mip_gap_abs']}.items():
            if (row.get('parameters') or {}).get(key) != value:
                error(f'Gurobi 生效参数 {key} 与冻结值不符')

    # error 终态有独立计数，仍复核该行的输入、参数及任何附带的数据。
    checked = verify_solution(model, row | {'status': 'limit_reached'} if status == 'error' else row, error)
    if status == 'error' and (row['values'] is not None or row['objective'] is not None):
        error('错误终态附带最终解')
    objective = checked if not issues else None
    event, goals, rejected = read_events(model, row, error)
    if status == 'infeasible' and goals:
        error('回调有已复核的可行解，但返回不可行')
    kind_of_first, first_secs = 'no_solution', None
    if kind == 'LP':
        kind_of_first = 'lp_return_upper' if objective is not None else 'lp_not_applicable'
        first_secs = elapsed if objective is not None else None
    elif event is not None:
        kind_of_first, first_secs = 'callback', event['elapsed_secs']
    elif objective is not None:
        kind_of_first, first_secs = 'return_upper', elapsed
    witnesses = goals + ([objective] if objective is not None else [])
    audit_truth(model, truth, row, witnesses, error)
    return dict(solver=solver, status=status, raw_status=row['raw_status'], objective=row['objective'],
                bound=row['bound'], objective_recomputed=objective, return_secs=elapsed,
                first=dict(source=kind_of_first, elapsed_secs=first_secs), first_checked_callback=event,
                events_accepted=len(goals), events_rejected=rejected, callback_errors=row.get('callback_errors') or [],
                reported_error=row.get('error'), order=row['order'], rep=row['rep'], issues=issues), witnesses


def compare_runs(model, members, error):
    witnesses = [goal for _, goals in members for goal in goals]
    if not witnesses:
        statuses = {run['status'] for run, _ in members}
        if {'infeasible', 'unbounded'} <= statuses:
            error('双方的不可行和无界结论冲突')
        return
    best = (min if model['sense'] == 'min' else max)(witnesses)
    for run, _ in members:
        issue = run['issues'].append
        if run['status'] == 'infeasible':
            issue('不可行结论与原模型已复核可行解矛盾')
        verify_bound(model, run, best, '双方已复核的最好可行目标', issue)
        if run['status'] == 'optimal' and finite(run['objective']) and better(model, best, run['objective']):
            issue('已报告最优值被原模型可行解改进')
        if run['status'] == 'unbounded' and any(other['status'] == 'optimal' for other, _ in members):
            error('双方的无界和证优结论冲突')
    optima = [run['objective'] for run, _ in members if run['status'] == 'optimal' and finite(run['objective'])]
    if len(optima) == 2 and not close(*optima):
        error(f'双方的已报告最优目标不一致 {optima}')


def time_stats(pairs):
    faster = Counter('mic' if m < g else 'gurobi' if m > g else 'tie' for m, g in pairs)
    ratios = [m / g for m, g in pairs if m > 0 and g > 0]
    return dict(pairs=len(pairs), mic_faster=faster['mic'], gurobi_faster=faster['gurobi'],
                equal=faster['tie'], ratio_count=len(ratios), median_ratio=median(ratios) if ratios else None,
                geometric_mean_ratio=geometric_mean(ratios) if ratios else None)


def summarize(cases):
    by_solver = {solver: Counter() for solver in SOLVERS}
    proofs, firsts = [], []
    for case in cases:
        runs = case['runs']
        for solver, run in runs.items():
            counts = by_solver[solver]
            counts[run['status']] += 1
            if run['proof_return_upper_secs'] is not None:
                counts['verified_optimal'] += 1
                counts['verified_optimal_within_budget'] += run['proof_return_upper_secs'] <= case['time_limit_secs']
            if run['objective_recomputed'] is not None:
                counts['checked_final_vectors'] += 1
            counts[run['first']['source']] += 1
            counts['rejected_callbacks'] += len(run['events_rejected'])
        if any(solver not in runs for solver in SOLVERS):
            continue
        a, b = (runs[solver] for solver in SOLVERS)
        optimal = a['proof_return_upper_secs'], b['proof_return_upper_secs']
        if all(value is not None for value in optimal):
            proofs.append(optimal)
        if not case['issues'] and not a['issues'] and not b['issues'] and all(
            run['first']['source'] == 'callback' for run in (a, b)
        ):
            firsts.append((a['first']['elapsed_secs'], b['first']['elapsed_secs']))
    return dict(cases=len(cases), lp=sum(case['kind'] == 'LP' for case in cases),
                milp=sum(case['kind'] == 'MILP' for case in cases), by_solver=by_solver,
                proof=time_stats(proofs), proof_at_least_one_ms=time_stats([p for p in proofs if min(p) >= MICROSECONDS]),
                first_callback=time_stats(firsts))


def verify(base):
    manifest = json.loads((base / 'manifest.json').read_bytes())
    freeze = json.loads((base / 'freeze.json').read_bytes())
    public = public_truths()
    issues, raw_by_case, cases = [], defaultdict(dict), []
    expected = {case['case_id']: case for case in manifest['cases']}
    if not expected or len(expected) != len(manifest['cases']):
        issues.append('冻结清单没有算例或含重复编号')
    raw_path = base / 'timings.jsonl.gz'
    with gzip.open(raw_path, 'rt', encoding='utf-8') as stream:
        for line, text in enumerate(stream, 1):
            try:
                row = json.loads(text)
                case_id, solver = row['case_id'], row['solver']
                if case_id not in expected or solver not in SOLVERS:
                    issues.append(f'原始记录第 {line} 行出现未冻结的算例或求解器 {case_id}/{solver}')
                elif solver in raw_by_case[case_id]:
                    issues.append(f'原始记录第 {line} 行重复运行 {case_id}/{solver}')
                else:
                    raw_by_case[case_id][solver] = (line, row)
            except (ValueError, KeyError, TypeError) as err:
                issues.append(f'原始记录第 {line} 行结构无效 {err}')

    for case_id, case in expected.items():
        case_issues = []
        item = dict(case_id=case_id, instance=case['instance'], suite=case['suite'], kind=None,
                    time_limit_secs=freeze['time_limit_secs'], truth=None, issues=case_issues, runs={})
        cases.append(item)
        try:
            model, truth = load_case(base, case, public, case_issues.append)
            item['kind'] = 'MILP' if any(v['kind'] != 'continuous' for v in model['variables']) else 'LP'
            item['truth'] = truth
        except (OSError, EOFError, KeyError, IndexError, TypeError, ValueError) as err:
            case_issues.append(f'无法复核原模型或参考值 {err}')
            continue
        members = []
        for solver in SOLVERS:
            raw = raw_by_case[case_id].get(solver)
            if raw is None:
                case_issues.append(f'缺少 {SOLVER_NAMES[solver]} 运行记录')
                continue
            line, row = raw
            try:
                verified, goals = audit_run(model, item['kind'], truth, freeze, case, row, [])
                verified['raw_line'] = line
                item['runs'][solver] = verified
                members.append((verified, goals))
            except (ArithmeticError, AttributeError, KeyError, TypeError, ValueError) as err:
                case_issues.append(f'无法复核原始记录第 {line} 行 {err}')
        compare_runs(model, members, case_issues.append)
        if len(members) == 2 and members[0][0]['order'] == members[1][0]['order']:
            case_issues.append('双方运行次序重复')
        for run, _ in members:
            trusted = not case_issues and not run['issues'] and valid_time(run['return_secs'])
            status, elapsed = run['status'], run['return_secs']
            run['proof_return_upper_secs'] = elapsed if trusted and status == 'optimal' else None
            run['conclusion_return_upper_secs'] = elapsed if trusted and status in CONCLUSIONS else None
            run['verified'] = trusted and status != 'error'

    issues.extend(f'{case["instance"]} {issue}' for case in cases for issue in case['issues'])
    issues.extend(f'{case["instance"]}/{solver} {issue}' for case in cases for solver, run in case['runs'].items()
                  for issue in run['issues'])
    counts = summarize(cases)
    return dict(passed=not issues, issues=issues, raw_records=sum(map(len, raw_by_case.values())),
                expected_records=len(expected) * len(SOLVERS), stats=counts, cases=cases,
                raw_sha256=sha256(raw_path.read_bytes()).hexdigest(),
                manifest_sha256=sha256((base / 'manifest.json').read_bytes()).hexdigest(),
                freeze_sha256=sha256((base / 'freeze.json').read_bytes()).hexdigest())
