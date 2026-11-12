#!/usr/bin/env python3
"""按事前矩阵串行测量缩放，交替新旧顺序，并检查位摘要和容量门槛。"""
import argparse
import hashlib
import json
import re
import statistics
import subprocess
import sys
from pathlib import Path


def invoke(binary, arguments):
    command = [str(binary), *arguments]
    if sys.platform == 'darwin':
        command = ['/usr/bin/time', '-l', *command]
    proc = subprocess.run(command, capture_output=True, text=True, timeout=180)
    if proc.returncode:
        raise RuntimeError(f'{command}: {proc.returncode}\n{proc.stdout}\n{proc.stderr}')
    row = json.loads(proc.stdout.strip().splitlines()[-1])
    rss = re.search(r'(\d+)\s+maximum resident set size', proc.stderr)
    row['peak_rss_bytes'] = int(rss.group(1)) if rss else None
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    cases = []
    for n in (1000, 3000, 8000):
        for kind in ('continuous', 'integer', 'mixed'):
            cases.append((f'scale-{n}-{kind}', ['scale', str(n), '8', kind], True))
    for name in ('lp-wide-1000x2000', 'lp-sparse-4000x800', 'milp-mixed-250i40c', 'setcover-300x600'):
        cases.append((name, ['file', str(root / 'bench' / 'instances' / f'{name}.milp')], True))
    for n in (1000, 8000):
        cases.append((f'solve-{n}', ['solve', str(n)], True))
    for n in (25000, 50000, 100000):
        cases.append((f'capacity-{n}', ['scale', str(n), '8', 'mixed'], False))
    cases.append(('solve-100000', ['solve', '100000'], False))
    binaries = {'baseline': args.baseline.resolve(), 'candidate': args.candidate.resolve()}
    hashes = {k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in binaries.items()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    with args.output.open('x') as output:
        for rep in range(3):
            for case, arguments, paired in cases:
                order = ['baseline', 'candidate'] if paired else ['candidate']
                if rep % 2:
                    order.reverse()
                for variant in order:
                    row = invoke(binaries[variant], arguments)
                    row.update(case=case, variant=variant, rep=rep, binary_sha256=hashes[variant], arguments=arguments)
                    output.write(json.dumps(row, ensure_ascii=False) + '\n')
                    output.flush()
                    rows.append(row)
                    print(case, variant, rep, round(row['operation_secs'], 6), row['peak_rss_bytes'], flush=True)
    summary = []
    for case, _, paired in cases:
        grouped = {v: [r for r in rows if r['case'] == case and r['variant'] == v] for v in binaries}
        candidate = grouped['candidate']
        assert len({r['digest'] for r in candidate}) == 1, case
        item = {'case': case, 'candidate_secs': statistics.median(r['operation_secs'] for r in candidate),
                'candidate_rss': statistics.median(r['peak_rss_bytes'] for r in candidate)}
        if candidate[0]['mode'] == 'solve':
            assert all(r['status'] == 'Optimal' and r['objective'] == 0 for r in candidate), case
        if paired:
            baseline = grouped['baseline']
            fields = ('mode', 'rows', 'cols', 'nnz', 'digest', 'status', 'objective', 'nodes', 'lp_iterations')
            assert all(tuple(r.get(k) for k in fields) == tuple(candidate[0].get(k) for k in fields)
                       for r in baseline + candidate), case
            item['baseline_secs'] = statistics.median(r['operation_secs'] for r in baseline)
            item['baseline_rss'] = statistics.median(r['peak_rss_bytes'] for r in baseline)
            item['time_ratio'] = item['candidate_secs'] / item['baseline_secs']
            if candidate[0]['mode'] == 'scale' and item['baseline_secs'] >= 0.01:
                assert item['time_ratio'] <= 1.25, item
        summary.append(item)
    by_case = {r['case']: r for r in summary}
    target = by_case['scale-8000-mixed']
    assert target['time_ratio'] <= 0.5, target
    assert target['candidate_rss'] <= target['baseline_rss'] * 0.5, target
    for n in (25000, 50000):
        assert by_case[f'capacity-{2*n}']['candidate_rss'] <= 2.5 * by_case[f'capacity-{n}']['candidate_rss']
    assert by_case['capacity-100000']['candidate_rss'] < 512 * 1024**2
    result = {'passed': True, 'runs': len(rows), 'summary': summary, 'binaries': hashes}
    args.output.with_suffix('.summary.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
