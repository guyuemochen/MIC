#!/usr/bin/env python3
"""新旧父域实现固定节点预算串行对照；每轮反转先后顺序。"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--baseline', type=Path, required=True)
p.add_argument('--candidate', type=Path, required=True)
p.add_argument('--directory', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--plan', choices=['regression', 'memory'], required=True)
a = p.parse_args()
if a.plan == 'memory':
    names = ['packing_correlated_s12017', 'packing_correlated_s12043', 'sparse_packing']
    modes, repeats, budget = ['hybrid'], 3, 200000
else:
    names = ['packing_many_rows_s12017', 'packing_many_rows_s12043', 'packing_correlated_s12017', 'packing_correlated_s12043', 'sparse_packing', 'covering_120x200_s12109', 'negative_coefficient', 'boundary_n200']
    modes, repeats, budget = ['dfs', 'hybrid', 'hybrid_rc'], 1, 50000
with a.output.open('w') as out:
    for rep in range(repeats):
        for name in names:
            for mode in modes:
                variants = [('baseline', a.baseline), ('candidate', a.candidate)]
                if rep % 2: variants.reverse()
                for variant, binary in variants:
                    cmd = [str(binary), str(a.directory / (name + '.milp')), mode, str(budget)]
                    if sys.platform == 'darwin': cmd = ['/usr/bin/time', '-l'] + cmd
                    try:
                        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180, check=True)
                        row = json.loads(proc.stdout.strip().splitlines()[-1])
                        match = re.search(r'(\d+)\s+maximum resident set size', proc.stderr)
                        row['peak_rss_bytes'] = int(match.group(1)) if match else None
                    except Exception as exc:
                        row = dict(status='driver_error', error=str(exc), stdout=str(getattr(exc, 'stdout', '') or ''), stderr=str(getattr(exc, 'stderr', '') or ''))
                    row.update(instance=name, mode=mode, variant=variant, rep=rep, node_limit=budget, threads=1, time_limit_secs=0)
                    out.write(json.dumps(row) + '\n'); out.flush()
                    print(name, mode, variant, rep, row.get('status'), row.get('wall_secs'), flush=True)
