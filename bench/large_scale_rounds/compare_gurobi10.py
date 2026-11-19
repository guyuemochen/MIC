#!/usr/bin/env python3
"""在冻结的 20 例上串行测量本机安装的 Gurobi 10 C API。"""
import argparse
import ctypes as ct
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

from check_results import ABS_TOL
from compare import invoke, preflight_rss, require, result, sha
from run_probe import AUTO_CASES, METADATA
from run_proof_diagnostics import MODERN, baseline_rows
from run_root_work import CASES, NODE_LIMIT


SECONDS = 30
GRB_INFINITY = 1e100
DATA_NOT_AVAILABLE = 10005
NO_LICENSE = 10009
SIZE_LIMIT = 10010
SENSES = {'min': 1, 'max': -1}
VAR_TYPES = {'continuous': b'C', 'integer': b'I', 'binary': b'B'}
ROW_TYPES = {'<=': b'<', '>=': b'>', '=': b'='}
STATUSES = {
    2: ('optimal', 'OPTIMAL'), 3: ('infeasible', 'INFEASIBLE'),
    4: ('infeasible_or_unbounded', 'INF_OR_UNBD'), 5: ('unbounded', 'UNBOUNDED'),
    6: ('limit_reached', 'CUTOFF'), 7: ('limit_reached', 'ITERATION_LIMIT'),
    8: ('limit_reached', 'NODE_LIMIT'), 9: ('limit_reached', 'TIME_LIMIT'),
    10: ('limit_reached', 'SOLUTION_LIMIT'), 11: ('limit_reached', 'INTERRUPTED'),
    12: ('limit_reached', 'NUMERIC'), 13: ('limit_reached', 'SUBOPTIMAL'),
    15: ('limit_reached', 'USER_OBJ_LIMIT'), 16: ('limit_reached', 'WORK_LIMIT'),
    17: ('limit_reached', 'MEM_LIMIT'),
}
RELIABLE_BOUND_STATUSES = {2, 7, 8, 9, 10, 11, 15, 16, 17}


class GurobiError(RuntimeError):
    def __init__(self, operation, code, reason):
        self.operation = operation
        self.code = code
        self.reason = reason
        super().__init__(f'{operation}: Gurobi 错误 {code}，{reason}')


def load_library(path):
    library = ct.CDLL(str(path))
    pointer, integer, real, string = ct.c_void_p, ct.c_int, ct.c_double, ct.c_char_p
    ip, dp, pp = map(ct.POINTER, (integer, real, pointer))
    names = ct.POINTER(string)
    types = ct.POINTER(ct.c_char)
    prototypes = {
        'GRBversion': (None, ip, ip, ip),
        'GRBemptyenv': (integer, pp), 'GRBstartenv': (integer, pointer),
        'GRBfreeenv': (None, pointer), 'GRBfreemodel': (integer, pointer),
        'GRBresetparams': (integer, pointer),
        'GRBgeterrormsg': (string, pointer), 'GRBgetenv': (pointer, pointer),
        'GRBsetintparam': (integer, pointer, string, integer),
        'GRBsetdblparam': (integer, pointer, string, real),
        'GRBgetintparam': (integer, pointer, string, ip),
        'GRBgetdblparam': (integer, pointer, string, dp),
        'GRBloadmodel': (integer, pointer, pp, string, integer, integer, integer, real,
                         dp, types, dp, ip, ip, ip, dp, dp, dp, types, names, names),
        'GRBgetintattr': (integer, pointer, string, ip),
        'GRBgetdblattr': (integer, pointer, string, dp),
        'GRBgetdblattrarray': (integer, pointer, string, integer, integer, dp),
        'GRBgetcharattrarray': (integer, pointer, string, integer, integer, types),
        'GRBgetvars': (integer, pointer, ip, ip, ip, dp, integer, integer),
        'GRBoptimize': (integer, pointer),
    }
    for name, (restype, *arguments) in prototypes.items():
        function = getattr(library, name)
        function.restype = restype
        function.argtypes = arguments
    version = (ct.c_int(), ct.c_int(), ct.c_int())
    library.GRBversion(*(ct.byref(item) for item in version))
    require(version[0].value == 10, f'须使用 Gurobi 10，实际主版本为 {version[0].value}')
    return library, '.'.join(str(item.value) for item in version)


def check(library, env, code, operation):
    if code == 0:
        return
    if code == NO_LICENSE:
        reason = '无法启动经授权的 Gurobi 环境'
    elif code == SIZE_LIMIT:
        reason = '本机许可证限制该模型规模'
    elif operation in ('建立环境', '启动环境'):
        reason = 'Gurobi 环境初始化失败'
    else:
        raw = library.GRBgeterrormsg(env) if env else None
        reason = raw.decode('utf-8', errors='replace') if raw else '无错误详情'
    raise GurobiError(operation, code, reason)


def configure(library, env, seconds, nodes):
    expected_ints = {'OutputFlag': 0, 'Threads': 1}
    expected_reals = {'TimeLimit': float(seconds), 'NodeLimit': float(nodes),
                      'MIPGap': 0.0, 'MIPGapAbs': 0.0, 'IntFeasTol': ABS_TOL}
    for setter, values in ((library.GRBsetintparam, expected_ints), (library.GRBsetdblparam, expected_reals)):
        for name, value in values.items():
            check(library, env, setter(env, name.encode(), value), f'设置 {name}')
    return expected_ints, expected_reals


def scalar(library, env, function, handle, name, kind, optional=False):
    value = kind()
    code = function(handle, name.encode(), ct.byref(value))
    if optional and code == DATA_NOT_AVAILABLE:
        return None
    check(library, env, code, f'读取 {name}')
    return value.value


def effective_params(library, env, expected_ints, expected_reals):
    actual = {}
    for getter, kind, expected in ((library.GRBgetintparam, ct.c_int, expected_ints),
                                   (library.GRBgetdblparam, ct.c_double, expected_reals)):
        for name, value in expected.items():
            actual[name] = scalar(library, env, getter, env, name, kind)
            require(actual[name] == value, f'模型环境参数 {name} 与冻结值不同')
    actual['DualReductions'] = scalar(library, env, library.GRBgetintparam, env, 'DualReductions', ct.c_int)
    for name in ('FeasibilityTol', 'OptimalityTol'):
        actual[name] = scalar(library, env, library.GRBgetdblparam, env, name, ct.c_double)
    return actual


def as_array(kind, values):
    return (kind * len(values))(*values)


def model_arrays(data):
    variables, rows = data['variables'], data['rows']
    columns = [[] for _ in variables]
    for row_index, row in enumerate(rows):
        for column, value in row['terms']:
            columns[column].append((row_index, value))
    starts, lengths, row_indices, coefficients = [], [], [], []
    for entries in columns:
        starts.append(len(coefficients))
        lengths.append(len(entries))
        for index, coefficient in entries:
            row_indices.append(index)
            coefficients.append(coefficient)
    objective = dict(data['objective'])
    lower = [var['lb'] for var in variables]
    upper = [GRB_INFINITY if var['ub'] is None else var['ub'] for var in variables]
    return (
        as_array(ct.c_double, [objective.get(j, 0.0) for j in range(len(variables))]),
        as_array(ct.c_char, [ROW_TYPES[row['op']] for row in rows]),
        as_array(ct.c_double, [row['rhs'] for row in rows]),
        as_array(ct.c_int, starts), as_array(ct.c_int, lengths), as_array(ct.c_int, row_indices),
        as_array(ct.c_double, coefficients), as_array(ct.c_double, lower), as_array(ct.c_double, upper),
        as_array(ct.c_char, [VAR_TYPES[var['kind']] for var in variables]),
    )


def model_array(library, env, function, model, name, kind, count):
    values = (kind * count)()
    if count:
        check(library, env, function(model, name.encode(), 0, count, values), f'读取 {name}')
    return values[:]


def verify_model(library, env, model, data, arrays):
    def int_attr(name):
        return scalar(library, env, library.GRBgetintattr, model, name, ct.c_int)

    objective, senses, rhs, starts, lengths, row_indices, terms, lower, upper, types = arrays
    require((int_attr('NumVars'), int_attr('NumConstrs'), int_attr('NumNZs'), int_attr('ModelSense'))
            == (len(objective), len(rhs), len(terms), SENSES[data['sense']]),
            'Gurobi 读回模型的维度、非零计数或目标方向改变')
    constant = scalar(library, env, library.GRBgetdblattr, model, 'ObjCon', ct.c_double)
    require(constant == data['offset'], 'Gurobi 读回的目标常数改变')
    for name, expected in (('Obj', objective), ('LB', lower), ('UB', upper), ('RHS', rhs)):
        observed = model_array(library, env, library.GRBgetdblattrarray, model, name, ct.c_double, len(expected))
        require(observed == expected[:], f'Gurobi 读回的 {name} 与原模型不同')
    for name, expected in (('Sense', senses), ('VType', types)):
        observed = model_array(library, env, library.GRBgetcharattrarray, model, name, ct.c_char, len(expected))
        require(observed == expected[:], f'Gurobi 读回的 {name} 与原模型不同')
    count = ct.c_int()
    actual_starts, actual_rows, actual_terms = ((ct.c_int * len(starts))(),
                                                (ct.c_int * len(terms))(), (ct.c_double * len(terms))())
    if len(starts):
        check(library, env, library.GRBgetvars(model, ct.byref(count), actual_starts, actual_rows,
                                               actual_terms, 0, len(starts)), '读回原模型矩阵')
    require(count.value == len(terms) and actual_starts[:] == starts[:], 'Gurobi 列内非零计数改变')
    for start, length in zip(starts, lengths, strict=True):
        end = start + length
        actual = sorted(zip(actual_rows[start:end], actual_terms[start:end], strict=True))
        expected = list(zip(row_indices[start:end], terms[start:end], strict=True))
        require(actual == expected, 'Gurobi 读回的原模型约束系数不同')


def maybe_bound(value):
    return value if value is not None and math.isfinite(value) and abs(value) < GRB_INFINITY else None


def effort(value):
    if value is None:
        return -1
    require(math.isfinite(value) and value >= 0 and value.is_integer(), '求解器返回未知的计数')
    return int(value)


def optimize(library, env, model, version, data):
    def int_attr(name):
        return scalar(library, env, library.GRBgetintattr, model, name, ct.c_int)

    def dbl_attr(name, optional=False):
        return scalar(library, env, library.GRBgetdblattr, model, name, ct.c_double, optional)

    start = time.perf_counter()
    code = library.GRBoptimize(model)
    elapsed = time.perf_counter() - start
    check(library, env, code, '求解')
    status_code, count = int_attr('Status'), int_attr('SolCount')
    status, raw_status = STATUSES.get(status_code, ('error', f'UNKNOWN_STATUS_{status_code}'))
    has_solution = count > 0 and status in ('optimal', 'limit_reached')
    if status == 'optimal':
        require(has_solution, '最优状态没有解向量')
    values = (model_array(library, env, library.GRBgetdblattrarray, model, 'X', ct.c_double,
                          len(data['variables'])) if has_solution else None)
    objective = dbl_attr('ObjVal') if has_solution else None
    bound = maybe_bound(dbl_attr('ObjBound', optional=True)) if status_code in RELIABLE_BOUND_STATUSES else None
    bound_source = 'ObjBound' if bound is not None else None
    if status == 'optimal' and bound is None and not int_attr('IsMIP'):
        bound, bound_source = objective, 'certified_optimal_objective'
    row = result(status, raw_status, elapsed, objective, bound, values, effort(dbl_attr('NodeCount', True)),
                 effort(dbl_attr('IterCount', True)), version)
    barrier_iterations = scalar(library, env, library.GRBgetintattr, model, 'BarIterCount', ct.c_int, True)
    row.update(raw_status_code=status_code, solution_count=count, bound_source=bound_source,
               optimizer_runtime_secs=dbl_attr('Runtime'), barrier_iterations=barrier_iterations,
               model_roundtrip_passed=True)
    return row


def solve(path, data, seconds, nodes):
    library, version = load_library(path)
    env, model = ct.c_void_p(), ct.c_void_p()
    parameters = None
    try:
        check(library, env, library.GRBemptyenv(ct.byref(env)), '建立环境')
        check(library, env, library.GRBsetintparam(env, b'OutputFlag', 0), '关闭环境日志')
        check(library, env, library.GRBstartenv(env), '启动环境')
        check(library, env, library.GRBresetparams(env), '恢复 Gurobi 默认参数')
        ints, reals = configure(library, env, seconds, nodes)
        arrays = model_arrays(data)
        check(library, env, library.GRBloadmodel(env, ct.byref(model), b'frozen-mic-model',
                                                 len(data['variables']), len(data['rows']), SENSES[data['sense']],
                                                 data['offset'], *arrays, None, None), '载入模型')
        model_env = library.GRBgetenv(model)
        require(model_env, 'Gurobi 未返回模型环境')
        parameters = effective_params(library, model_env, ints, reals)
        verify_model(library, model_env, model, data, arrays)
        row = optimize(library, model_env, model, version, data)
    except GurobiError as error:
        row = result('error', f'GUROBI_ERROR_{error.code}', None, None, None, None, -1, -1, version)
        row.update(error=str(error), error_code=error.code, error_stage=error.operation)
    finally:
        if model:
            library.GRBfreemodel(model)
        if env:
            library.GRBfreeenv(env)
    row['parameters'] = parameters
    return row


def tiny_model(size):
    variables = [{'name': f'x{j}', 'kind': 'binary', 'lb': 0.0, 'ub': 1.0} for j in range(size)]
    rows = [{'op': '<=', 'rhs': 1.0, 'terms': [[j, 1.0]]} for j in range(size)]
    return {'variables': variables, 'rows': rows, 'sense': 'max', 'offset': 2.25,
            'objective': [[size - 1, 1.5]]}


def preflight(library):
    results = {}
    for size in (2, 2001):
        row = solve(library, tiny_model(size), 1, 1)
        valid = (row['status'] == 'optimal' and row['objective'] == row['bound'] == 3.75
                 and len(row['values']) == size and row['values'][-1] == 1.0)
        require(valid or row.get('error_code') in (NO_LICENSE, SIZE_LIMIT),
                f'{size} 行列微型探针返回无效结果：{row["raw_status"]}')
        results[str(size)] = {key: row.get(key) for key in ('version', 'status', 'raw_status', 'error_code',
                                                           'objective', 'bound', 'parameters')}
        if size == 2 and not valid:
            return results, True
    return results, False


def stage(legacy, modern):
    staged = {}
    for baseline, cases in ((legacy, AUTO_CASES), (modern, MODERN)):
        for name, row in baseline_rows(baseline, cases).items():
            model = (baseline.parent / row['model']).resolve()
            require(row['schema'] == 1 and sha(model) == row['model_sha256'], f'{name} 原模型摘要不同')
            staged[name] = model, {key: row[key] for key in METADATA}
    return staged


def run(library, legacy, modern, output):
    require(not output.exists(), f'已有 Gurobi 结果 {output}')
    preflight_rss()
    _, version = load_library(library)
    library_digest, host = sha(library), platform.platform()
    models = stage(legacy, modern)
    script = str(Path(__file__).resolve())
    output.parent.mkdir(parents=True, exist_ok=True)
    failed = False
    with output.open('x', encoding='utf-8') as stream:
        for index, name in enumerate(CASES):
            model, reference = models[name]
            command = [sys.executable, script, '_worker', str(library), str(model)]
            try:
                row = invoke(command, SECONDS + 120, rss=True)
                require(row['version'] == version, '子进程 Gurobi 版本不一致')
            except (ValueError, RuntimeError, OSError, subprocess.TimeoutExpired) as error:
                row = result('error', type(error).__name__, None, None, None, None, -1, -1, version)
                row.update(error=str(error), peak_rss_bytes=None)
            failed |= row['status'] == 'error'
            row = dict(reference, **row)
            row.update(solver='gurobi', variant='gurobi', rep=0, order=0, case_order=index,
                       threads=1, time_limit_secs=SECONDS, node_limit=NODE_LIMIT, host=host,
                       model=os.path.relpath(model, output.parent), command=command,
                       library_sha256=library_digest, model_source='archive_json')
            stream.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
            stream.flush()
            print(name, row['status'], row['raw_status'], row['elapsed_secs'], flush=True)
    return failed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest='action', required=True)
    probe = actions.add_parser('probe')
    probe.add_argument('--library', type=Path, required=True)
    matrix = actions.add_parser('run')
    for field in ('library', 'legacy-baseline', 'modern-baseline', 'output'):
        matrix.add_argument(f'--{field}', type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.action == 'probe':
            report, failed = preflight(args.library.resolve())
            print(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False))
            return failed
        paths = (args.library, args.legacy_baseline, args.modern_baseline, args.output)
        return run(*(path.resolve() for path in paths))
    except (KeyError, TypeError, ValueError, RuntimeError, OSError, subprocess.SubprocessError) as error:
        print(f'Gurobi 10 对照失败：{error}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    if len(sys.argv) == 4 and sys.argv[1] == '_worker':
        library, model = map(Path, sys.argv[2:])
        row = solve(library, json.loads(model.read_text(encoding='utf-8')), SECONDS, NODE_LIMIT)
        print(json.dumps(row, ensure_ascii=False, allow_nan=False))
    else:
        sys.exit(main())
