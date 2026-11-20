#!/usr/bin/env python3
"""沿用旧验收门槛汇总复查；限额退出不按耗时算加速。"""
import argparse
from collections import defaultdict
import json
import math
import statistics
from pathlib import Path
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--screen', type=Path, required=True)
p.add_argument('--paired', type=Path, required=True)
a = p.parse_args()
summary = []
proof_losses = []
ratios = []
for path, phase in [(a.screen, 'screen'), (a.paired, 'paired')]:
    groups = defaultdict(dict)
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    for r in rows:
        assert r['status'] != 'driver_error', r
        groups[r['instance']].setdefault(r['mode'], []).append(r)
    for name, modes in sorted(groups.items()):
        assert set(modes) == {'existing', 'auto'}, name
        before, after = modes['existing'], modes['auto']
        by_rep = {r['rep']: r for r in after}
        for r in before:
            if r['status'] == 'Optimal' and by_rep[r['rep']]['status'] != 'Optimal':
                proof_losses.append(dict(phase=phase, instance=name, rep=r['rep']))
        if phase == 'screen': continue
        item = dict(instance=name)
        for mode, rs in modes.items():
            item[mode] = dict(runs=len(rs), median_secs=statistics.median(r['wall_secs'] for r in rs), min_secs=min(r['wall_secs'] for r in rs), max_secs=max(r['wall_secs'] for r in rs), peak_rss_bytes=max(r['peak_rss_bytes'] for r in rs), statuses=sorted(set(r['status'] for r in rs)), objectives=sorted(set(r['objective'] for r in rs if r['objective'] is not None)))
        comparable = all(r['status'] == 'Optimal' for r in before + after) and item['existing']['median_secs'] >= .1
        item['time_ratio'] = item['auto']['median_secs']/item['existing']['median_secs'] if comparable else None
        if comparable: ratios.append(item['time_ratio'])
        item['rss_limit_bytes'] = min(256*2**20, max(64*2**20, 2*item['existing']['peak_rss_bytes']))
        item['rss_passed'] = item['auto']['peak_rss_bytes'] <= item['rss_limit_bytes']
        summary.append(item)
geo = math.exp(statistics.mean(math.log(r) for r in ratios)) if ratios else None
print(json.dumps(dict(proof_losses=proof_losses, comparable_instances=len(ratios), geo_time_ratio=geo,
    full_acceptance_established=False, tested_gates_passed=bool(ratios) and geo <= .9 and max(ratios) <= 1.25 and not proof_losses and all(r['rss_passed'] for r in summary),
    summary=summary), indent=2))
