#!/usr/bin/env python3
"""串行对照两例运输模型的首个可行解与证优耗时。"""

import argparse
from collections import defaultdict
import ctypes as ct
from datetime import datetime, UTC
import json
import os
from pathlib import Path
import platform
from statistics import median
import subprocess
import sys
import time

from check_results import close, finite, inspect_results, verify_solution
from compare import invoke, require, result, sha
from compare_gurobi10 import check, load_library, read_outcome, solve, tiny_model
from run_probe import METADATA
from run_proof_diagnostics import baseline_rows


CASES = ('drayage-100-23', 'drayage-25-23')
RESULTS = Path(__file__).resolve().parent / 'results'
BASELINE = RESULTS / 'miplib2017-baseline.jsonl'
TRUTH = RESULTS / 'miplib2017-sources.json'
VERSION = '10.0.3'
NODES = 2_000_000
SECONDS = 30
REPETITIONS = 3
MIPSOL = 4
MIPSOL_SOL = 4001
MIPSOL_OBJ = 4002
MIPSOL_SOLCNT = 4006
RUNTIME = 6002
CALLBACK = ct.CFUNCTYPE(ct.c_int, ct.c_void_p, ct.c_void_p, ct.c_int, ct.c_void_p)


def observe_gurobi(library, env, model, version, data):
    require(version == VERSION, f'冻结的动态库版本为 {VERSION}，实际版本为 {version}')
    library.GRBsetcallbackfunc.argtypes = (ct.c_void_p, CALLBACK, ct.c_void_p)
    library.GRBsetcallbackfunc.restype = ct.c_int
    library.GRBcbget.argtypes = (ct.c_void_p, ct.c_int, ct.c_int, ct.c_void_p)
    library.GRBcbget.restype = ct.c_int
    events, callback_errors = [], []
    vector = (ct.c_double * len(data['variables']))()

    @CALLBACK
    def callback(_model, cbdata, where, _userdata):
        if where != MIPSOL:
            return 0
        elapsed = time.perf_counter() - start
        try:
            objective, runtime, count = ct.c_double(), ct.c_double(), ct.c_int()
            queries = ((MIPSOL_OBJ, ct.byref(objective)), (MIPSOL_SOL, vector),
                       (RUNTIME, ct.byref(runtime)), (MIPSOL_SOLCNT, ct.byref(count)))
            for key, output in queries:
                code = library.GRBcbget(cbdata, where, key, output)
                if code:
                    callback_errors.append(dict(elapsed_secs=elapsed, query=key, error_code=code))
                    return 0
            events.append(dict(elapsed_secs=elapsed, optimizer_runtime_secs=runtime.value,
                               solver_solution_count=count.value, objective=objective.value,
                               values=vector[:]))
        except Exception as error:
            callback_errors.append(dict(elapsed_secs=elapsed, error=str(error)))
        return 0

    check(library, env, library.GRBsetcallbackfunc(model, callback, None), '安装 MIPSOL 回调')
    start = time.perf_counter()
    code = library.GRBoptimize(model)
    elapsed = time.perf_counter() - start
    check(library, env, code, '求解')
    row = read_outcome(library, env, model, version, data, elapsed)
    # 保留异常候选给离线复核，单个坏向量不应中断余下回调。
    for event in events:
        if not finite(event['objective']):
            event['objective'] = None
        event['values'] = [value if finite(value) else None for value in event['values']]
    row.update(incumbent_events=events, callback_errors=callback_errors)
    return row


def probe(library):
    results = {}
    for size in (2, 2001):
        model = tiny_model(size)
        row = solve(library, model, 3, NODES, observe_gurobi)
        checked = check_events(model, row)
        require(row['status'] == 'optimal' and row['objective'] == row['bound'] == 3.75
                and row['values'][-1] == 1 and checked['final_vector_event'] is not None
                and checked['first'] is not None and not checked['errors'],
                f'{size} 行列的授权或回调预检失败：{row["raw_status"]}')
        results[str(size)] = dict(version=row['version'], final_objective=row['objective'],
                                  parameters=row['parameters'], verified=checked)
    return results


def skip_columns(model):
    """从原时间差分行的正二元项确定各任务的自环列。"""
    kind = [var['kind'] for var in model['variables']]
    clocks = kind.count('continuous')
    skipped = defaultdict(set)
    for row in model['rows']:
        if row['op'] == '=' or len(row['terms']) != 4:
            continue
        sign = 1 if row['op'] == '>=' else -1
        continuous = [(j, sign * value) for j, value in row['terms'] if kind[j] == 'continuous']
        binary = [(j, sign * value) for j, value in row['terms'] if kind[j] == 'binary']
        if sorted(value for _, value in continuous) != [-1, 1] or len(binary) != 2:
            continue
        (earlier,) = (j for j, value in continuous if value == -1)
        (skip,) = (j for j, value in binary if value > 0)
        require(any(value < 0 for _, value in binary), '时间差分行须含行程激活列')
        skipped[earlier].add(skip)
    require(clocks == 65 and len(skipped) == clocks and all(len(columns) == 1 for columns in skipped.values()),
            '每个运输时钟须对应唯一的任务自环')
    columns = tuple(sorted(next(iter(columns)) for columns in skipped.values()))
    require(len(set(columns)) == clocks, '多个运输任务共用同一自环')
    return columns


def check_events(model, row, skipped=None):
    first = first_route = final_event = final_tasks = None
    accepted, rejected, errors = 0, [], []
    elapsed = row['elapsed_secs']
    if not finite(elapsed) or elapsed < 0:
        return dict(first=first, first_route=first_route, final_vector_event=final_event,
                    final_tasks_served=final_tasks, accepted=accepted, rejected=rejected,
                    errors=['没有有效求解返回时间'])
    previous = 0.0
    for index, event in enumerate(row.get('incumbent_events', []), 1):
        seen = event['elapsed_secs']
        if not finite(seen) or not previous <= seen <= elapsed:
            errors.append(f'第 {index} 次回调的时间顺序有误')
            continue
        previous = seen
        candidate_errors = []
        claims = dict(status='limit_reached', bound=None, auto_degraded=False,
                      objective=event['objective'], values=event['values'])
        goal = verify_solution(model, claims, candidate_errors.append)
        if goal is None or candidate_errors:
            if not candidate_errors:
                candidate_errors.append('缺少可复核的解与目标')
            rejected.append(dict(event=index, elapsed_secs=seen, errors=candidate_errors[:8],
                                 violations=len(candidate_errors)))
            continue
        accepted += 1
        tasks = sum(event['values'][j] < 0.5 for j in skipped) if skipped is not None else None
        metrics = dict(event=index, elapsed_secs=seen, objective=goal, tasks_served=tasks)
        if first is None:
            first = metrics
        if first_route is None and tasks:
            first_route = metrics
        if final_event is None and row['values'] is not None and finite(row['objective']) \
                and close(goal, row['objective']):
            if len(event['values']) == len(row['values']) and all(
                    close(a, b) for a, b in zip(event['values'], row['values'])):
                final_event = index
    if row.get('callback_errors'):
        errors.append(f'回调读取失败 {row["callback_errors"]}')
    if row['values'] is not None and first is None:
        errors.append('返回了可行解，但没有通过原模型核验的回调向量')
    if skipped is not None and row['values'] is not None and all(j < len(row['values']) for j in skipped):
        final_tasks = sum(row['values'][j] < 0.5 for j in skipped)
    return dict(first=first, first_route=first_route, final_vector_event=final_event,
                final_tasks_served=final_tasks, accepted=accepted, rejected=rejected, errors=errors)


def mic_command(probe_binary, case, seconds, mode):
    return [str(probe_binary), mode, str(case), '--time-limit-secs', str(seconds),
            '--node-limit', str(NODES), '--threads', '1', '--early-stop', 'off']


def run(args):
    require(args.seconds > 0 and args.repetitions > 0, '时限和重复次数须为正')
    _, version = load_library(args.library)
    require(version == VERSION, f'冻结的动态库版本为 {VERSION}，实际版本为 {version}')
    source = Path(__file__).resolve().parents[2]
    source_commit = subprocess.check_output(['git', '-C', source, 'rev-parse', 'HEAD'], text=True).strip()
    source_status = subprocess.check_output(
        ['git', '-C', source, 'status', '--porcelain', '--untracked-files=all'], text=True)
    require(not source_status, '冻结前须提交并核对全部源码')
    staged = {}
    for name, record in baseline_rows(BASELINE, CASES).items():
        model = (BASELINE.parent / record['model']).resolve()
        instance = args.instances / f'{name}.milp'
        require(sha(model) == record['model_sha256'] and sha(instance) == record['milp_sha256'],
                f'{name} 原模型或 MIC 输入的摘要不同')
        parsed = invoke([str(args.mic), 'export', str(instance)], 10)
        data = json.loads(model.read_text(encoding='utf-8'))
        require(all(parsed[key] == data[key] for key in parsed), f'{name} 完整 MIC 输入与归档模型不同')
        staged[name] = model, instance, {key: record[key] for key in METADATA}

    args.output.mkdir(parents=True)
    host = platform.platform()
    script = Path(__file__).resolve()
    binary_sha, library_sha = sha(args.mic), sha(args.library)
    metadata = dict(source_commit=source_commit, host=host, python=sys.version,
                    freeze_command=[sys.executable, *sys.argv],
                    rustc=subprocess.check_output(['rustc', '-Vv'], text=True).strip(),
                    ref_manifest_sha256=sha(TRUTH), baseline_sha256=sha(BASELINE),
                    frozen_at_utc=datetime.now(UTC).isoformat(), library=str(args.library),
                    library_sha256=library_sha, mic_binary=str(args.mic), mic_binary_sha256=binary_sha,
                    gurobi_version=version, benchmark_sha256=sha(script),
                    instances_directory=str(args.instances), time_limit_secs=args.seconds, node_limit=NODES,
                    threads=1, mip_gap=0, mip_gap_abs=0, reps=args.repetitions,
                    input_hashes={name: {key: data[key] for key in
                                         ('model_sha256', 'milp_sha256', 'mps_sha256')}
                                  for name, (_, _, data) in staged.items()})
    (args.output / 'freeze.json').write_text(render(metadata), encoding='utf-8')
    failed = False
    with (args.output / 'timings.jsonl').open('x', encoding='utf-8') as stream:
        for rep in range(args.repetitions):
            cases = CASES if rep % 2 == 0 else tuple(reversed(CASES))
            for case_order, name in enumerate(cases):
                model, instance, reference = staged[name]
                gurobi = [sys.executable, str(script), '_worker', str(args.library), str(model), str(args.seconds)]
                mic = mic_command(args.mic, instance, args.seconds, 'timed-solve')
                commands = [('gurobi', gurobi), ('mic', mic)]
                if (rep + CASES.index(name)) % 2:
                    commands.reverse()
                for order, (solver, command) in enumerate(commands):
                    try:
                        row = invoke(command, args.seconds + 30)
                    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as error:
                        row = result('error', type(error).__name__, None, None, None, None,
                                     -1, -1, version if solver == 'gurobi' else source_commit)
                        row['error'] = str(error)
                    merged = dict(reference, **row)
                    merged.update(solver=solver, variant='drayage-timing', rep=rep, order=order,
                                  case_order=case_order, time_limit_secs=args.seconds, node_limit=NODES,
                                  host=host, model=os.path.relpath(model, args.output), command=command,
                                  proof_elapsed_secs=row['elapsed_secs'] if row['status'] == 'optimal' else None)
                    if solver == 'mic':
                        merged['version'] = source_commit
                    failed |= row['status'] == 'error'
                    stream.write(json.dumps(merged, ensure_ascii=False, allow_nan=False) + '\n')
                    stream.flush()
                    print(name, solver, row['raw_status'], row['elapsed_secs'], flush=True)
    return failed


def render(obj):
    return json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n'


def verify(path, truth):
    final = inspect_results([path], truth)
    rows, models = [], {}
    with path.open(encoding='utf-8') as stream:
        for line in stream:
            row = json.loads(line)
            model_path = path.parent / row['model']
            if model_path not in models:
                model = json.loads(model_path.read_text(encoding='utf-8'))
                models[model_path] = model, skip_columns(model)
            model, skipped = models[model_path]
            observed = check_events(model, row, skipped)
            first = observed['first']
            elapsed = row['elapsed_secs']
            cap = row['time_limit_secs']
            if row['status'] != 'optimal' and row['proof_elapsed_secs'] is not None:
                observed['errors'].append('非最优终态记录了证优时间')
            if row['status'] == 'optimal' and row['proof_elapsed_secs'] != elapsed:
                observed['errors'].append('证优的正常返回时间记录不一致')
            rows.append(dict(
                instance=row['instance'], solver=row['solver'], rep=row['rep'],
                status=row['status'], raw_status=row['raw_status'], return_secs=elapsed,
                time_limit_secs=cap, final_objective=row['objective'], final_bound=row['bound'],
                proof_return_upper_secs=row['proof_elapsed_secs'],
                proof_right_censored_at_secs=min(elapsed, cap) if row['status'] == 'limit_reached' and finite(elapsed) else None,
                first_right_censored_at_secs=(min(elapsed, cap) if first is None and finite(elapsed) and
                                              row['status'] == 'limit_reached' else None),
                **observed,
            ))
    groups = defaultdict(list)
    for row in rows:
        groups[row['instance'], row['solver'], row['time_limit_secs']].append(row)
    summary = []
    for (instance, solver, cap), runs in sorted(groups.items()):
        first = [row['first']['elapsed_secs'] for row in runs
                 if row['first'] is not None and row['first']['elapsed_secs'] <= cap]
        route = [row['first_route']['elapsed_secs'] for row in runs
                 if row['first_route'] is not None and row['first_route']['elapsed_secs'] <= cap]
        optimal = [row['proof_return_upper_secs'] for row in runs
                   if row['proof_return_upper_secs'] is not None and row['proof_return_upper_secs'] <= cap]
        summary.append(dict(instance=instance, solver=solver, time_limit_secs=cap, n=len(runs),
                            first_within_limit=len(first), route_within_limit=len(route),
                            proof_within_limit=len(optimal),
                            first_median_secs=median(first) if len(first) == len(runs) else None,
                            first_route_median_secs=median(route) if len(route) == len(runs) else None,
                            proof_return_median_secs=median(optimal) if len(optimal) == len(runs) else None))
    return dict(passed=final['passed'] and not any(row['errors'] for row in rows),
                final_checks=final, rows=rows, summary=summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest='action', required=True)
    license_probe = subcommands.add_parser('probe')
    license_probe.add_argument('--library', type=Path, required=True)
    matrix = subcommands.add_parser('run')
    for parameter in ('library', 'mic', 'instances', 'output'):
        matrix.add_argument(f'--{parameter}', type=lambda arg: Path(arg).resolve(), required=True)
    matrix.add_argument('--seconds', type=int, default=SECONDS)
    matrix.add_argument('--repetitions', type=int, default=REPETITIONS)
    check = subcommands.add_parser('verify')
    check.add_argument('timings', type=Path)
    check.add_argument('--truth', type=Path, default=TRUTH)
    args = parser.parse_args()
    if args.action == 'probe':
        print(render(probe(args.library.resolve())), end='')
        return False
    if args.action == 'run':
        return run(args)
    checked = verify(args.timings, args.truth)
    print(render(checked), end='')
    return not checked['passed']


if __name__ == '__main__':
    if len(sys.argv) == 5 and sys.argv[1] == '_worker':
        library, model = map(Path, sys.argv[2:4])
        data = json.loads(model.read_text(encoding='utf-8'))
        row = solve(library, data, int(sys.argv[4]), NODES, observe_gurobi)
        print(json.dumps(row, ensure_ascii=False, allow_nan=False))
    else:
        sys.exit(main())
