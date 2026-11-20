#!/usr/bin/env python3
"""串行测量 MIC 与 Gurobi 10.0.3，逐次归档完整向量和原始输出。"""

import argparse
from datetime import UTC, datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time

from compare import invoke, require, result, sha
from compare_drayage_timing import NODES, SECONDS, VERSION, mic_command, observe_gurobi, probe
from compare_gurobi10 import load_library, solve


SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[2]
RESULTS = SCRIPT.parent / 'results'
BUILD_ARGS = ['build', '--manifest-path', 'milp/Cargo.toml', '--release', '--example', 'round_probe', '--locked']
SOURCE_GROUPS = (
    ('repository', ROOT / 'milp/instances'),
    ('benchmark', ROOT / 'bench/instances'),
    ('public', Path('/tmp/mic-rounds-public')),
    ('miplib2017', Path('/tmp/mic-rounds-miplib2017')),
    ('diagnostic', RESULTS / 'phase1-progress-validation/two-equalities.milp'),
)
ENV_OVERRIDES = {'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'}


def now():
    return datetime.now(UTC).isoformat()


def json_bytes(data, **kwargs):
    return (json.dumps(data, allow_nan=False, **kwargs) + '\n').encode('utf-8')


def write_json(path, data):
    path.write_bytes(json_bytes(data, indent=2, ensure_ascii=False))


def append_gzip(path, row):
    # 每条记录都是完整的 gzip 成员，外部复核可读取已完成的行。
    with path.open('ab') as stream:
        stream.write(gzip.compress(json_bytes(row), compresslevel=3, mtime=0))


def read_model(path):
    payload = gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes()
    return json.loads(payload)


def archived_models():
    models = {}
    for name in ('external-baseline.jsonl', 'miplib2017-baseline.jsonl'):
        with (RESULTS / name).open(encoding='utf-8') as stream:
            for line in stream:
                row = json.loads(line)
                models.setdefault((row['instance'], row['milp_sha256']), row)
    return models


def provenance():
    dependencies = ('compare.py', 'compare_drayage_timing.py', 'compare_gurobi10.py', 'check_results.py',
                    'run_probe.py', 'run_proof_diagnostics.py', 'run_root_work.py')
    return {
        str(path.relative_to(ROOT)): sha(path)
        for path in sorted([*ROOT.glob('milp/src/**/*.rs'), ROOT / 'milp/Cargo.toml',
                            ROOT / 'milp/Cargo.lock', ROOT / 'milp/examples/round_probe.rs', SCRIPT,
                            *(SCRIPT.parent / name for name in dependencies)])
    }


def stage(source_groups, mic, output, requested=None):
    archived = archived_models()
    (output / 'models').mkdir()
    (output / 'inputs').mkdir()
    (output / 'mps').mkdir()
    cases = []
    names = set()
    for suite, source in source_groups:
        paths = sorted(source.glob('*.milp')) if source.is_dir() else [source]
        if requested:
            paths = [path for path in paths if path.stem in requested]
            if not paths:
                continue
        require(paths and all(path.is_file() for path in paths), f'{suite} 缺少算例：{source}')
        for path in paths:
            source_path = path.absolute()
            name = path.stem
            require(name not in names, f'算例名重复：{name}')
            names.add(name)
            payload = path.read_bytes()
            digest = hashlib.sha256(payload).hexdigest()
            command = [str(mic), 'export', str(source_path)]
            parsed = invoke(command, 180)
            baseline = archived.get((name, digest))
            mps = source_path.with_suffix('.mps')
            mps_digest = sha(mps) if mps.exists() else None
            mps_archive = None
            if mps_digest is not None:
                compressed = source_path.with_suffix('.mps.gz').read_bytes()
                require(hashlib.sha256(gzip.decompress(compressed)).hexdigest() == mps_digest,
                        f'{name} 原始 MPS 压缩档摘要不同')
                mps_archive = output / 'mps' / f'{name}.mps.gz'
                mps_archive.write_bytes(compressed)
            if baseline:
                model = RESULTS / baseline['model']
                model_digest = sha(model)
                original = read_model(model)
                require(model_digest == baseline['model_sha256']
                        and all(parsed[field] == original[field] for field in parsed),
                        f'{name} 当前导出与归档模型不一致')
                require(mps_digest == baseline['mps_sha256'], f'{name} 的原始 MPS 摘要不符')
            else:
                model_payload = json_bytes(parsed, ensure_ascii=False, separators=(',', ':'), sort_keys=True)
                model_digest = hashlib.sha256(model_payload).hexdigest()
                model = output / 'models' / f'{name}.json.gz'
                model.write_bytes(gzip.compress(model_payload, compresslevel=3, mtime=0))
            if source_path.is_relative_to(ROOT):
                original_path = source_path.relative_to(ROOT)
                durable_input = source_path
            else:
                original_path = source_path
                durable_input = output / 'inputs' / f'{name}.milp.gz'
                durable_input.write_bytes(gzip.compress(payload, compresslevel=3, mtime=0))
            cases.append(dict(case_id=f'{suite}/{name}', instance=name, suite=suite,
                              source_path=str(original_path), milp_sha256=digest,
                              input=os.path.relpath(durable_input, output),
                              mps=os.path.relpath(mps_archive, output) if mps_archive else None,
                              mps_sha256=mps_digest,
                              model=os.path.relpath(model, output), model_sha256=model_digest,
                              variables=len(parsed['variables']), rows=len(parsed['rows']),
                              integer_variables=sum(var['kind'] != 'continuous' for var in parsed['variables']),
                              export_command=command))
    if requested:
        missing = sorted(requested - names)
        require(not missing, f'找不到指定算例：{", ".join(missing)}')
    return cases


def capture(command, timeout):
    started, begin = now(), time.perf_counter()
    killed, code, output, errors = False, None, b'', b''
    spawn_error = None
    try:
        with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              env=dict(os.environ, **ENV_OVERRIDES), start_new_session=True) as proc:
            try:
                output, errors = proc.communicate(timeout=timeout)
            except subprocess.TimeoutExpired:
                killed = True
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                output, errors = proc.communicate()
            code = proc.returncode
    except OSError as error:
        spawn_error = str(error)
    return dict(command=command, started_at_utc=started, finished_at_utc=now(),
                process_elapsed_secs=time.perf_counter() - begin, returncode=code,
                timeout=killed, spawn_error=spawn_error,
                stdout=output.decode('utf-8', 'surrogateescape'),
                stderr=errors.decode('utf-8', 'surrogateescape'))


def parse_outcome(raw, version):
    try:
        require(not raw['timeout'], '外层进程超时并已终止')
        require(raw['spawn_error'] is None, raw['spawn_error'])
        require(raw['returncode'] == 0, f'进程退出码 {raw["returncode"]}')
        lines = raw['stdout'].strip().splitlines()
        require(lines, '子进程未输出 JSON')
        row = json.loads(lines[-1])
        require({'status', 'raw_status', 'elapsed_secs', 'objective', 'bound', 'values', 'nodes'} <= row.keys(),
                '求解器 JSON 缺少必要字段')
        require(not row.get('proof_profile'), 'MIC 编译了诊断计时特性')
    except (TypeError, ValueError, AttributeError) as error:
        row = result('error', type(error).__name__, None, None, None, None, -1, -1, version)
        row['error'] = str(error)
    row.setdefault('incumbent_events', [])
    row.setdefault('callback_errors', [])
    return row


def run(args):
    env = dict(os.environ, **ENV_OVERRIDES)
    os.environ.update(ENV_OVERRIDES)
    _, version = load_library(args.library)
    require(version == VERSION, f'冻结的动态库版本为 {VERSION}，实际版本为 {version}')
    args.output.mkdir(parents=True)
    checkout = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    source_tree = provenance()
    mic = args.target_dir / 'release/examples/round_probe'
    build = ['cargo', *BUILD_ARGS, '--target-dir', str(args.target_dir)]
    preflight = capture(build, 180)
    write_json(args.output / 'build.json', preflight)
    require(preflight['returncode'] == 0, '生产 release 构建失败，见 build.json')
    require(source_tree == provenance(), '构建期间源码改变')
    frozen_at = now()
    freeze = dict(schema=1, source_commit=checkout, source_sha256=source_tree, host=platform.platform(),
                  python=sys.version, rustc=subprocess.check_output(['rustc', '-Vv'], env=env, text=True).strip(),
                  cargo=subprocess.check_output(['cargo', '-V'], env=env, text=True).strip(),
                  freeze_command=[sys.executable, *sys.argv], build_command=build,
                  child_environment=ENV_OVERRIDES, frozen_at_utc=frozen_at,
                  benchmark_sha256=sha(SCRIPT), library=str(args.library), library_sha256=sha(args.library),
                  mic_binary=str(mic), mic_binary_sha256=sha(mic), gurobi_version=version,
                  time_limit_secs=args.seconds, node_limit=NODES, threads=1, mip_gap=0, mip_gap_abs=0,
                  mic_early_stop='off', repetitions=1)
    write_json(args.output / 'freeze.json', freeze)
    write_json(args.output / 'gurobi-preflight.json', probe(args.library))
    cases = stage([*SOURCE_GROUPS, *args.extra_source], mic, args.output, set(args.case))
    pair_order = ['gurobi', 'mic']
    write_json(args.output / 'manifest.json', dict(schema=1, cases=cases, first_pair_order=pair_order,
                                                  model_sha256_scope='uncompressed_bytes',
                                                  expected_runs=2 * len(cases)))
    print(f'[{now()}] 已完成预检，开始串行测量 {len(cases)} 例、{len(cases) * 2} 次', flush=True)
    errors = []
    for index, case in enumerate(cases):
        commands = dict(
            gurobi=[sys.executable, str(SCRIPT), '_worker', str(args.library),
                    str(args.output / case['model']), str(args.seconds)],
            mic=mic_command(mic, ROOT / case['source_path'], args.seconds, 'timed-solve'),
        )
        for order, solver in enumerate(pair_order):
            raw = capture(commands[solver], args.seconds + 60)
            raw.update(case_id=case['case_id'], instance=case['instance'], solver=solver,
                       case_order=index, order=order, rep=0)
            append_gzip(args.output / 'raw.jsonl.gz', raw)
            solver_version = checkout if solver == 'mic' else version
            row = parse_outcome(raw, solver_version)
            row.update(case, solver=solver, version=solver_version, schema=1, rep=0, order=order,
                       case_order=index, time_limit_secs=args.seconds, node_limit=NODES, threads=1,
                       variant='all-gurobi10-timing', host=freeze['host'], command=commands[solver],
                       started_at_utc=raw['started_at_utc'], finished_at_utc=raw['finished_at_utc'],
                       raw_returncode=raw['returncode'],
                       proof_elapsed_secs=row['elapsed_secs'] if row['status'] == 'optimal' else None)
            append_gzip(args.output / 'timings.jsonl.gz', row)
            if row['status'] == 'error':
                errors.append(dict(case_id=case['case_id'], solver=solver, raw_status=row['raw_status'],
                                   error=row.get('error')))
            print(f'[{now()}] {index + 1}/{len(cases)} {case["instance"]} {solver} '
                  f'{row["raw_status"]} {row["elapsed_secs"]}', flush=True)
        pair_order.reverse()
    write_json(args.output / 'completed.json', dict(finished_at_utc=now(), runs=len(cases) * 2, errors=errors,
                                                   timings_sha256=sha(args.output / 'timings.jsonl.gz'),
                                                   raw_sha256=sha(args.output / 'raw.jsonl.gz')))
    return bool(errors)


def extra_source(value):
    suite, separator, path = value.partition('=')
    if not separator or not suite or not path:
        raise argparse.ArgumentTypeError('请使用 SUITE=PATH')
    return suite, Path(path).absolute()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('library', 'output', 'target-dir'):
        parser.add_argument(f'--{option}', type=lambda value: Path(value).absolute(), required=True)
    parser.add_argument('--extra-source', action='append', type=extra_source, default=[])
    parser.add_argument('--case', action='append', default=[], help='只测指定名称，可重复传入')
    parser.add_argument('--seconds', type=int, default=SECONDS)
    args = parser.parse_args()
    if args.seconds <= 0:
        parser.error('--seconds 必须大于 0')
    return run(args)


if __name__ == '__main__':
    if len(sys.argv) == 5 and sys.argv[1] == '_worker':
        library_path, model_path = map(Path, sys.argv[2:4])
        observed = solve(library_path, read_model(model_path), int(sys.argv[4]), NODES, observe_gurobi)
        print(json.dumps(observed, ensure_ascii=False, allow_nan=False))
    else:
        sys.exit(main())
