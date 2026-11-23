#!/usr/bin/env python3
"""在同一台机器轮换两个 MIC 版本和 Gurobi，逐条复核可行解、回调及最优界。"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import platform
import statistics
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROUNDS = HERE.parents[1]
ROOT = ROUNDS.parents[1]
sys.path.insert(0, str(ROUNDS))
from check_results import verify_solution, verify_bound
from compare_all_timing import read_model

MATRIX = ROUNDS / 'results/main-gurobi10-100-60s'
CASES = {'setcover-400x800': 113, 'setcover-300x600': 144}
ORDER = [('baseline', 'current', 'gurobi'), ('current', 'gurobi', 'baseline'),
         ('gurobi', 'baseline', 'current')]
FLAGS = ['--time-limit-secs', '60', '--node-limit', '2000000', '--threads', '1',
         '--early-stop', 'off']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(model, result, goal):
    errors = []
    verified = verify_solution(model, result, errors.append)
    for event in result['incumbent_events']:
        verify_solution(model, dict(status='limit_reached', objective=event['objective'],
                                    values=event['values'], bound=None, auto_degraded=False),
                        errors.append)
    verify_bound(model, result, goal, '基线证优值', errors.append)
    if (errors or result.get('callback_errors') or result['status'] != 'optimal'
            or abs(verified - goal) > 1e-6 or abs(result['bound'] - goal) > 1e-6):
        raise ValueError(f'最终解、界或回调核验失败 {errors}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ('baseline', 'current', 'library', 'output'):
        parser.add_argument('--' + field, type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    binaries = {'baseline': args.baseline.resolve(), 'current': args.current.resolve()}
    library = args.library.resolve()
    manifest = json.loads((MATRIX / 'manifest.json').read_text())['cases']
    inputs = {c['instance']: c for c in manifest if c['instance'] in CASES}
    results = []
    summary = {}
    for n, goal in CASES.items():
        entry = inputs[n]
        source = MATRIX / entry['input']
        model_path = MATRIX / entry['model']
        model = read_model(model_path)
        if sha(source) != entry['milp_sha256'] or model['sense'] != 'min':
            raise ValueError(f'{n} 输入与存档不符')
        for binary in binaries.values():
            exported = json.loads(subprocess.check_output([binary, 'export', source]))
            if any(value != model[key] for key, value in exported.items()):
                raise ValueError(f'{n} 两种输入表示不同')
        for rep, order in enumerate(ORDER, 1):
            for solver in order:
                cmd = ([sys.executable, ROUNDS / 'compare_all_timing.py', '_worker', library,
                        model_path, '60'] if solver == 'gurobi' else
                       [binaries[solver], 'timed-solve', source, *FLAGS])
                data = subprocess.check_output(cmd, text=True, timeout=100)
                r = json.loads(data.strip().splitlines()[-1])
                check(model, r, goal)
                results.append(dict(instance=n, solver=solver, rep=rep, **r))
                print(n, solver, rep, r['elapsed_secs'], r['nodes'], flush=True)
        summary[n] = {}
        for solver in ORDER[0]:
            samples = [r for r in results if r['solver'] == solver and r['instance'] == n]
            seconds = [r['elapsed_secs'] for r in samples]
            summary[n][solver] = dict(
                status='optimal', objective=goal, bound=goal, seconds=seconds,
                median_seconds=statistics.median(seconds),
                nodes=[r['nodes'] for r in samples],
                lp_iterations=[r['lp_iterations'] for r in samples],
                first_feasible_seconds=[r['incumbent_events'][0]['elapsed_secs'] for r in samples],
                optimal_incumbent_seconds=[next(e['elapsed_secs'] for e in r['incumbent_events']
                                                if abs(e['objective'] - goal) < 1e-6) for r in samples])
    raw = ''.join(json.dumps(r, ensure_ascii=False) + '\n' for r in results)
    (args.output / 'runs.jsonl.gz').write_bytes(gzip.compress(raw.encode(), mtime=0))
    grb = next(r for r in results if r['solver'] == 'gurobi')
    metadata = dict(
        platform=platform.platform(), baseline_commit='7e2c76671f3cb3647458bde72ddf6e43bfc629e1',
        binary_sha256={k: sha(v) for k, v in binaries.items()},
        source_sha256={str(p.relative_to(ROOT)): sha(p) for p in sorted((ROOT / 'milp/src').rglob('*.rs'))},
        input_sha256={n: inputs[n]['milp_sha256'] for n in CASES},
        gurobi_library_sha256=sha(library), gurobi_version=grb['version'],
        gurobi_parameters=grb['parameters'], mic_flags=FLAGS, order=ORDER, results=summary)
    (args.output / 'timing.json').write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + '\n')


if __name__ == '__main__':
    main()
