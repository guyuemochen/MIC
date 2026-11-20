#!/usr/bin/env python3
"""串行对比 MIC、HiGHS、SCIP，保留完整解、原模型和进程峰值内存。"""
import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import signal
import subprocess
import sys
import tempfile
import time


def row_bounds(row):
    rhs = row['rhs']
    return (-math.inf if row['op'] == '<=' else rhs,
            math.inf if row['op'] == '>=' else rhs)


def finite(value):
    return value if value is not None and math.isfinite(value) else None


def result(status, raw_status, elapsed, objective, bound, values, nodes, iterations, version):
    objective, bound = finite(objective), finite(bound)
    absolute = abs(objective - bound) if objective is not None and bound is not None else None
    return {
        'status': status, 'raw_status': raw_status, 'elapsed_secs': elapsed,
        'objective': objective, 'bound': bound, 'gap_absolute': finite(absolute),
        'gap_percent': (0.0 if absolute == 0 else
                        finite(100 * absolute / abs(objective)) if absolute is not None and objective != 0 else None),
        'values': values, 'nodes': nodes if nodes >= 0 else None,
        'lp_iterations': iterations if iterations >= 0 else None, 'version': version,
        'cuts': None, 'rc_fixed': None, 'rc_reapplies': None, 'rc_refixed': None,
        'presolve_reductions': None, 'auto_degraded': None, 'gap_terminated': None,
    }


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_original(model, columns, costs, bounds, integers, rows, sense, offset):
    """校对独立 MPS 读取结果，返回 MIC 原列顺序在外部模型中的位置。"""
    by_name = {name: j for j, name in enumerate(columns)}
    expected_names = model['mps_columns']
    require(len(columns) == len(by_name) and set(columns) == set(expected_names),
            'MPS 的列名与转换文件不同')
    order = [by_name[name] for name in expected_names]
    mapping = {j: i for i, j in enumerate(order)}
    objective = dict(model['objective'])
    require((sense, offset) == (model['sense'], model['offset']), 'MPS 目标方向或常数不同')
    for i, (j, var) in enumerate(zip(order, model['variables'], strict=True)):
        upper = math.inf if var['ub'] is None else var['ub']
        require((bounds[j], integers[j], costs[j]) ==
                ((var['lb'], upper), var['kind'] != 'continuous', objective.get(i, 0)),
                f'MPS 列 {columns[j]} 的变量界、类型或目标系数不同')

    def normalized(lo, hi, terms):
        terms = tuple(sorted(terms))
        if lo == hi:
            yield '=', lo, terms
        else:
            if math.isfinite(lo):
                yield '>=', lo, terms
            if math.isfinite(hi):
                yield '<=', hi, terms

    actual_rows = Counter(
        item for lo, hi, terms in rows
        for item in normalized(lo, hi, ((mapping[j], a) for j, a in terms if a != 0))
    )
    expected_rows = Counter(
        (row['op'], row['rhs'], tuple(map(tuple, row['terms']))) for row in model['rows']
    )
    require(actual_rows == expected_rows, 'MPS 约束与转换文件不同')
    return order


def highs_read(highspy, path):
    solver = highspy.Highs()
    require(solver.setOptionValue('output_flag', False) == highspy.HighsStatus.kOk,
            'HiGHS 设置静默失败')
    if path:
        require(solver.readModel(str(path)) == highspy.HighsStatus.kOk, 'HiGHS 读取 MPS 失败')
    return solver


def highs_model(highspy, solver, model, source):
    if source:
        lp = solver.getLp()
        require(lp.a_matrix_.format_ == highspy.MatrixFormat.kColwise, 'MPS 需要按列存储')
        matrix = lp.a_matrix_
        starts, indices, values = matrix.start_, matrix.index_, matrix.value_
        terms = [[] for _ in range(lp.num_row_)]
        for j in range(lp.num_col_):
            for k in range(starts[j], starts[j + 1]):
                terms[indices[k]].append((j, values[k]))
        integers = [v == highspy.HighsVarType.kInteger for v in lp.integrality_]
        integers = integers or [False] * lp.num_col_
        sense = 'min' if lp.sense_ == highspy.ObjSense.kMinimize else 'max'
        return verify_original(model, lp.col_names_, lp.col_cost_,
                               list(zip(lp.col_lower_, lp.col_upper_)), integers,
                               zip(lp.row_lower_, lp.row_upper_, terms), sense, lp.offset_)

    lp = highspy.HighsLp()
    lp.num_col_, lp.num_row_ = len(model['variables']), len(model['rows'])
    costs = dict(model['objective'])
    lp.col_cost_ = [costs.get(j, 0) for j in range(lp.num_col_)]
    lp.col_lower_ = [v['lb'] for v in model['variables']]
    lp.col_upper_ = [math.inf if v['ub'] is None else v['ub'] for v in model['variables']]
    lp.sense_ = highspy.ObjSense.kMinimize if model['sense'] == 'min' else highspy.ObjSense.kMaximize
    lp.offset_ = model['offset']
    lp.integrality_ = [highspy.HighsVarType.kContinuous if v['kind'] == 'continuous'
                       else highspy.HighsVarType.kInteger for v in model['variables']]
    lower, upper, starts, indices, values = [], [], [0], [], []
    for row in model['rows']:
        lo, hi = row_bounds(row)
        lower.append(lo)
        upper.append(hi)
        for j, value in row['terms']:
            indices.append(j)
            values.append(value)
        starts.append(len(indices))
    lp.row_lower_, lp.row_upper_ = lower, upper
    matrix = highspy.HighsSparseMatrix()
    matrix.format_ = highspy.MatrixFormat.kRowwise
    matrix.num_col_, matrix.num_row_ = lp.num_col_, lp.num_row_
    matrix.start_, matrix.index_, matrix.value_ = starts, indices, values
    lp.a_matrix_ = matrix
    require(solver.passModel(lp) == highspy.HighsStatus.kOk, 'HiGHS 载入模型失败')
    return range(lp.num_col_)


def run_highs(model, source, seconds, node_limit):
    import highspy

    solver = highs_read(highspy, source)
    for key, value in {'threads': 1, 'parallel': 'off', 'time_limit': float(seconds),
                       'mip_max_nodes': node_limit, 'mip_rel_gap': 0.0, 'mip_abs_gap': 0.0}.items():
        require(solver.setOptionValue(key, value) == highspy.HighsStatus.kOk,
                f'HiGHS 设置 {key} 失败')
    order = highs_model(highspy, solver, model, source)
    start = time.perf_counter()
    outcome = solver.run()
    elapsed = time.perf_counter() - start
    require(outcome != highspy.HighsStatus.kError, 'HiGHS 求解出错')
    raw = solver.getModelStatus()
    status = {
        'kOptimal': 'optimal', 'kInfeasible': 'infeasible', 'kUnbounded': 'unbounded',
        'kUnboundedOrInfeasible': 'infeasible_or_unbounded',
        'kTimeLimit': 'limit_reached', 'kIterationLimit': 'limit_reached',
        'kSolutionLimit': 'limit_reached', 'kObjectiveBound': 'limit_reached',
        'kObjectiveTarget': 'limit_reached', 'kMemoryLimit': 'limit_reached',
        'kInterrupt': 'limit_reached',
    }.get(raw.name, 'error')
    info = solver.getInfo()
    solution = solver.getSolution()
    has_solution = (solution.value_valid and
                    info.primal_solution_status == int(highspy.kSolutionStatusFeasible))
    columns = solution.col_value
    values = [columns[j] for j in order] if has_solution else None
    objective = info.objective_function_value if has_solution else None
    bound = objective if status == 'optimal' else None
    is_mip = any(v['kind'] != 'continuous' for v in model['variables'])
    if status == 'limit_reached' and is_mip:
        bound = info.mip_dual_bound
    return result(status, raw.name, elapsed, objective, bound, values,
                  info.mip_node_count, info.simplex_iteration_count, solver.version())


def scip_model(solver, model, source):
    if source:
        from ortools.linear_solver.python import model_builder
        imported = model_builder.Model()
        require(imported.import_from_mps_file(str(source)), 'OR-Tools 读取 MPS 失败')
        proto = imported.export_to_proto()
        order = verify_original(
            model, [v.name for v in proto.variable], [v.objective_coefficient for v in proto.variable],
            [(v.lower_bound, v.upper_bound) for v in proto.variable], [v.is_integer for v in proto.variable],
            ((c.lower_bound, c.upper_bound, zip(c.var_index, c.coefficient)) for c in proto.constraint),
            'max' if proto.maximize else 'min', proto.objective_offset,
        )
        error = solver.LoadModelFromProtoKeepNames(proto)
        require(not error, f'SCIP 载入原 MPS 失败：{error}')
        return order

    variables = []
    for v in model['variables']:
        upper = solver.infinity() if v['ub'] is None else v['ub']
        add = solver.NumVar if v['kind'] == 'continuous' else solver.IntVar
        variables.append(add(v['lb'], upper, v['name']))
    objective = solver.Objective()
    objective.SetOffset(model['offset'])
    objective.SetOptimizationDirection(model['sense'] == 'max')
    for j, value in model['objective']:
        objective.SetCoefficient(variables[j], value)
    for row in model['rows']:
        cons = solver.RowConstraint(*row_bounds(row), '')
        for j, value in row['terms']:
            cons.SetCoefficient(variables[j], value)
    return range(len(variables))


def run_scip(model, source, seconds, node_limit):
    from ortools.linear_solver import pywraplp

    solver = pywraplp.Solver.CreateSolver('SCIP')
    require(solver is not None, 'OR-Tools 未提供 SCIP')
    solver.SuppressOutput()
    require(solver.SetNumThreads(1), 'SCIP 设置线程失败')
    solver.SetTimeLimit(seconds * 1000)
    parameters = (f'limits/totalnodes = {node_limit}\nlimits/gap = 0\n'
                  'limits/absgap = 0\ntiming/clocktype = 2')
    require(solver.SetSolverSpecificParametersAsString(parameters), 'SCIP 设置节点上限或 gap 失败')
    order = scip_model(solver, model, source)
    start = time.perf_counter()
    raw = solver.Solve()
    elapsed = time.perf_counter() - start
    statuses = {
        solver.OPTIMAL: ('optimal', 'OPTIMAL'),
        solver.FEASIBLE: ('limit_reached', 'FEASIBLE'),
        solver.INFEASIBLE: ('infeasible', 'INFEASIBLE'),
        solver.UNBOUNDED: ('unbounded', 'UNBOUNDED'),
        solver.NOT_SOLVED: ('limit_reached', 'NOT_SOLVED'),
        solver.ABNORMAL: ('error', 'ABNORMAL'),
        solver.MODEL_INVALID: ('error', 'MODEL_INVALID'),
    }
    status, raw_status = statuses.get(raw, ('error', str(raw)))
    has_solution = raw in (solver.OPTIMAL, solver.FEASIBLE)
    variables = solver.variables()
    values = [variables[j].solution_value() for j in order] if has_solution else None
    objective = solver.Objective().Value() if has_solution else None
    bound = solver.Objective().BestBound() if has_solution else None
    if bound is not None and abs(bound) >= 1e20:
        bound = None
    return result(status, raw_status, elapsed, objective, bound, values, solver.nodes(),
                  solver.iterations(), solver.SolverVersion())


def child():
    backend, data, source, seconds, nodes = sys.argv[2:]
    if backend == 'columns':
        import highspy
        solver = highs_read(highspy, Path(source))
        print(json.dumps(solver.getLp().col_names_, allow_nan=False))
        return
    model = json.loads(Path(data).read_text(encoding='utf-8'))
    run = {'highs': run_highs, 'scip': run_scip}[backend]
    row = run(model, Path(source) if source else None, int(seconds), int(nodes))
    print(json.dumps(row, allow_nan=False))


def invoke(command, seconds, rss=False):
    prefix, pattern, unit = [], '', 1
    if rss and sys.platform == 'darwin':
        prefix, pattern = ['/usr/bin/time', '-l'], r'(\d+)\s+maximum resident set size'
    elif rss and sys.platform.startswith('linux'):
        prefix, pattern, unit = ['/usr/bin/time', '-v'], r'Maximum resident set size \(kbytes\):\s+(\d+)', 1024
    env = dict(os.environ, LC_ALL='C', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    with subprocess.Popen(prefix + command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True, env=env, start_new_session=True) as proc:
        try:
            stdout, stderr = proc.communicate(timeout=seconds)
        except subprocess.TimeoutExpired:
            # 计时器和求解进程一起停止，避免残留进程干扰后续串行测量。
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.communicate()
            raise
    if proc.returncode:
        raise RuntimeError(f'进程退出码 {proc.returncode}\n{stderr[-4000:]}\n{stdout[-1000:]}')
    lines = stdout.strip().splitlines()
    require(lines, '子进程未返回 JSON')
    row = json.loads(lines[-1])
    if rss:
        match = re.search(pattern, stderr) if pattern else None
        row['peak_rss_bytes'] = int(match[1]) * unit if match else None
    return row


class RssPreflightError(RuntimeError):
    """当前环境无法可靠测量独立进程的峰值 RSS。"""


def preflight_rss():
    try:
        sample = invoke([sys.executable, '-c', "print('{}')"], 10, rss=True)
        rss = sample.get('peak_rss_bytes')
        require(type(rss) in (int, float) and math.isfinite(rss) and rss > 0,
                f'未取得有限的正值峰值 RSS：{rss!r}')
    except (OSError, RuntimeError, ValueError, TypeError, subprocess.TimeoutExpired) as error:
        raise RssPreflightError(f'RSS 计时环境预检失败：{error}') from error


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def store_model(directory, model):
    payload = (json.dumps(model, ensure_ascii=False, allow_nan=False,
                          sort_keys=True, separators=(',', ':')) + '\n').encode('utf-8')
    digest = hashlib.sha256(payload).hexdigest()
    path = directory / f'{digest}.json'
    if not path.exists():
        # 完整写入后再用硬链接发布，已存在的模型不会被覆盖。
        with tempfile.TemporaryDirectory(dir=directory) as staging:
            candidate = Path(staging) / path.name
            candidate.write_bytes(payload)
            try:
                os.link(candidate, path)
            except FileExistsError:
                pass
    require(path.read_bytes() == payload, f'已保存的原模型与哈希不一致：{path}')
    return path, digest


def positive_integer(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError('须为正整数')
    return number


def cases(args):
    for path in args.models:
        yield path.stem, 'synthetic', path.resolve(), None
    if args.public:
        manifest = json.loads(args.public.read_text(encoding='utf-8'))
        for item in manifest['instances']:
            if args.public_group and item['group'] != args.public_group:
                continue
            if not item['included']:
                print(item['instance'], 'excluded', item['reason'], flush=True)
                continue
            milp = (args.public.parent / item['milp']).resolve()
            mps = (args.public.parent / item['mps']).resolve()
            require(sha(milp) == item['milp_sha256'] and sha(mps) == item['mps_sha256'],
                    f'公共算例 {item["instance"]} 文件与清单哈希不同')
            yield item['instance'], item['group'], milp, mps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('models', type=Path, nargs='*')
    parser.add_argument('--mic', required=True, type=Path)
    parser.add_argument('--mic-label', default='mic')
    parser.add_argument('--mic-option', action='append', default=[], metavar='--OPTION=VALUE')
    parser.add_argument('--solvers', nargs='+', choices=('mic', 'highs', 'scip'),
                        default=['mic', 'highs', 'scip'])
    parser.add_argument('--public', type=Path, help='fetch_public.py 生成的 sources.json')
    parser.add_argument('--public-group', choices=('regression', 'holdout'))
    parser.add_argument('--seconds', required=True, type=positive_integer)
    parser.add_argument('--node-limit', type=positive_integer, default=2_000_000)
    parser.add_argument('--reps', type=positive_integer, default=3)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.node_limit > 2_147_483_647:
        parser.error('节点上限须不超过 2147483647')
    if len(args.solvers) != len(set(args.solvers)):
        parser.error('求解器不可重复')
    mic_options = []
    for option in args.mic_option:
        key, sep, value = option.partition('=')
        if not key.startswith('--') or not sep or not value:
            parser.error('MIC 选项须为 --选项=值')
        if key in {'--node-limit', '--time-limit-secs', '--threads'}:
            parser.error(f'--mic-option 不接受 {key}，统一使用外层预算和单线程')
        mic_options.extend((key, value))
    roster = list(cases(args))
    if not roster or len({name for name, *_ in roster}) != len(roster):
        parser.error('必须给出算例，且每个算例名只能出现一次')
    output_path = args.output.resolve()
    model_dir = output_path.parent / 'models'
    if output_path.exists():
        parser.error('已有结果文件，请换新的输出路径')
    model_dir.mkdir(parents=True, exist_ok=True)
    mic = str(args.mic.resolve())
    mic_hash = sha(Path(mic))
    external = [sys.executable, str(Path(__file__).resolve()), '_child']
    staged = []
    for name, group, milp, mps in roster:
        exported = invoke([mic, 'export', str(milp)], 180)
        if mps:
            columns = invoke([*external, 'columns', '', str(mps), '0', '0'], 60)
            require(len(columns) == len(exported['variables']) and
                    all(v['name'] == f'x{j}' for j, v in enumerate(exported['variables'])),
                    'MPS 转换列无法映射回原文件')
            exported['mps_columns'] = columns
        path, digest = store_model(model_dir, exported)
        staged.append((name, group, milp, mps, path, digest, sha(milp), sha(mps) if mps else None))

    failures = 0
    with output_path.open('x', encoding='utf-8') as output:
        for rep in range(args.reps):
            for index, (name, group, milp, mps, data, model_hash, milp_hash, mps_hash) in enumerate(staged):
                offset = (rep + index) % len(args.solvers)
                order = args.solvers[offset:] + args.solvers[:offset]
                for position, backend in enumerate(order):
                    if backend == 'mic':
                        command = [mic, 'solve', str(milp), *mic_options, '--threads', '1',
                                   '--node-limit', str(args.node_limit), '--time-limit-secs',
                                   str(args.seconds)]
                    else:
                        command = [*external, backend, str(data), str(mps or ''),
                                   str(args.seconds), str(args.node_limit)]
                    try:
                        row = invoke(command, args.seconds + 120, rss=True)
                    except (ValueError, RuntimeError, OSError, subprocess.TimeoutExpired) as error:
                        row = result('error', type(error).__name__, None, None, None, None, -1, -1, None)
                        row.update(error=str(error), peak_rss_bytes=None)
                    if row['status'] == 'error':
                        failures += 1
                    if backend == 'mic':
                        row['version'] = mic_hash
                    row.update(schema=1, instance=name, group=group, rep=rep, solver=backend,
                               variant=args.mic_label if backend == 'mic' else backend,
                               order=position, threads=1, time_limit_secs=args.seconds,
                               node_limit=args.node_limit, model=str(data.relative_to(output_path.parent)),
                               model_sha256=model_hash, milp_sha256=milp_hash, mps_sha256=mps_hash,
                               mic_options=mic_options if backend == 'mic' else None,
                               host=platform.platform(), command=command)
                    output.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
                    output.flush()
                    print(name, backend, rep, row['status'], row['elapsed_secs'], flush=True)
    return bool(failures)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '_child':
        child()
    else:
        sys.exit(main())
