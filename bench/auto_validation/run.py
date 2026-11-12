#!/usr/bin/env python3
"""串行自动配置对照，交替执行顺序；保留限额、失败与峰值RSS。"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--directory', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--seconds', type=int, default=5)
p.add_argument('--repeats', type=int, default=1)
p.add_argument('--names', nargs='*')
p.add_argument('--modes', nargs='+', choices=['existing','auto','dfs_controls','hybrid'], default=['existing','auto'])
a = p.parse_args()
files = [a.directory / (n + '.milp') for n in a.names] if a.names else sorted(a.directory.glob('*.milp'))
with a.output.open('w') as out:
    for rep in range(a.repeats):
        for path in files:
            for mode in (a.modes if rep % 2 == 0 else a.modes[::-1]):
                cmd = [str(a.binary), str(path), mode, str(a.seconds)]
                if sys.platform == 'darwin':
                    cmd = ['/usr/bin/time', '-l'] + cmd
                try:
                    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=a.seconds + 30, check=True)
                    row = json.loads(proc.stdout.strip().splitlines()[-1])
                    match = re.search(r'(\d+)\s+maximum resident set size', proc.stderr)
                    row['peak_rss_bytes'] = int(match.group(1)) if match else None
                except Exception as exc:
                    row = {'status': 'driver_error', 'error': str(exc),
                           'stdout': str(getattr(exc, 'stdout', '') or ''),
                           'stderr': str(getattr(exc, 'stderr', '') or '')}
                row.update(instance=path.stem, mode=mode, rep=rep, time_limit_secs=a.seconds, threads=1)
                out.write(json.dumps(row, ensure_ascii=False) + '\n'); out.flush()
                print(path.stem, mode, rep, row.get('status'), row.get('objective'), row.get('wall_secs'), flush=True)
