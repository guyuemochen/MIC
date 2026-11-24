#!/usr/bin/env python3
"""在冻结模型上交错运行两个正式入口，核验所有原解、回调和已知目标。"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROUNDS = Path(__file__).resolve().parents[2]
FROZEN = ROUNDS / 'results/main-gurobi10-100-60s'
sys.path.insert(0, str(ROUNDS))
from check_results import verify_bound, verify_solution
from compare_drayage_timing import NODES
from summarize_all_timing import audit_truth, load_case, public_truths, read_bytes, read_events


def measure(command, name, model, truth, known, seconds):
    row = json.loads(subprocess.check_output(command, timeout=seconds + 60))
    errors = []
    final = verify_solution(model, row, errors.append)
    _, goals, rejected = read_events(model, row, errors.append)
    if known['objective'] is not None:
        verify_bound(model, row, known['objective'], '归档可行解', errors.append)
    witnesses = goals + ([final] if final is not None else [])
    audit_truth(model, truth, row, witnesses, errors.append)
    if known['status'] in ('optimal', 'infeasible', 'unbounded'):
        audit_truth(model, known, row, witnesses, errors.append)
    assert not (rejected or errors or row.get('callback_errors')), (name, errors, rejected)
    fields = ('status', 'objective', 'bound', 'elapsed_secs', 'nodes', 'lp_iterations',
              'cuts', 'presolve_reductions', 'version')
    return dict(**{k: row[k] for k in fields if k in row}, checked_callbacks=len(goals))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('old', type=Path)
    p.add_argument('new', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--seconds', type=int, default=2)
    p.add_argument('--reps', type=int, default=1)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--nodes', type=int, default=NODES)
    p.add_argument('--gurobi-library', type=Path)
    p.add_argument('--case', nargs='+', dest='cases')
    args = p.parse_args()
    if args.gurobi_library and (args.nodes != NODES or args.threads != 1):
        p.error('Gurobi 对照使用 1 个线程、200 万节点')
    binaries = {k: getattr(args, k).resolve() for k in ('old', 'new')}
    results = dict(sha256={k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in binaries.items()},
                   seconds=args.seconds, threads=args.threads, node_limit=args.nodes, cases={})
    if args.gurobi_library:
        results['sha256']['gurobi'] = hashlib.sha256(args.gurobi_library.read_bytes()).hexdigest()
    with gzip.open(ROUNDS/'results/remaining-after-weighted-setcover/mic-runs.jsonl.gz', 'rt') as f:
        known = {r['instance']: r for r in map(json.loads, f)}
    public = public_truths()
    manifest = json.loads((FROZEN/'manifest.json').read_bytes())['cases']
    selected = [c for c in manifest if not args.cases or c['instance'] in args.cases]
    assert not args.cases or {c['instance'] for c in selected} == set(args.cases)
    with tempfile.TemporaryDirectory() as scratch:
        for i, c in enumerate(selected):
            name = c['instance']
            errors = []
            model, truth = load_case(FROZEN, c, public, errors.append)
            assert not errors, errors
            prior = known[name]
            path = Path(scratch)/(name+'.milp')
            path.write_bytes(read_bytes(FROZEN/c['input']))
            mode = 'timed-solve' if args.threads == 1 else 'solve'
            opts = ['--time-limit-secs', str(args.seconds), '--node-limit', str(args.nodes),
                    '--threads', str(args.threads), '--early-stop', 'off']
            commands = {k: [v, mode, path, *opts] for k, v in binaries.items()}
            if args.gurobi_library:
                commands['gurobi'] = [sys.executable, ROUNDS/'compare_all_timing.py', '_worker',
                                      args.gurobi_library, FROZEN/c['model'], str(args.seconds)]
            pairs = results['cases'][name] = []
            for rep in range(args.reps):
                pair = {}
                for solver in (('old', 'new') if (rep+i)%2 == 0 else ('new', 'old')):
                    pair[solver] = measure(commands[solver], name, model, truth, prior, args.seconds)
                if args.gurobi_library:
                    pair['gurobi'] = measure(commands['gurobi'], name, model, truth, prior, args.seconds)
                for solver, result in pair.items():
                    for other, witness in pair.items():
                        if solver != other and witness['objective'] is not None:
                            verify_bound(model, result, witness['objective'], other, errors.append)
                assert not errors, (name, errors)
                pairs.append(pair)
            args.output.write_text(json.dumps(results, indent=2, ensure_ascii=False)+'\n')


if __name__ == '__main__':
    main()
