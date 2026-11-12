#!/usr/bin/env python3
"""核查不同父域表示下的确定性结果，汇总固定工作量时间与逐进程RSS。"""
import argparse
from collections import defaultdict
import json
import math
import statistics
from pathlib import Path
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('results', type=Path, nargs='+')
a = p.parse_args()
fields = ['status', 'objective', 'bound', 'gap_percent', 'nodes', 'lp_iterations', 'cuts', 'rc_fixed', 'rc_reapplies', 'rc_refixed', 'auto_degraded', 'values']
groups = defaultdict(list)
count = 0
for path in a.results:
    for line in path.read_text().splitlines():
        r = json.loads(line)
        assert r['status'] != 'driver_error', r
        assert not r['auto_degraded'], r
        groups[r['instance'], r['mode'], r['node_limit']].append(r)
        count += 1
summary = []
for (name, mode, budget), rows in sorted(groups.items()):
    signature = [rows[0][f] for f in fields]
    assert all([r[f] for f in fields] == signature for r in rows), (name, mode, budget, '结果或工作量改变')
    by_variant = {v: [r for r in rows if r['variant'] == v] for v in ('baseline', 'candidate')}
    assert all(by_variant.values()), (name, '缺少新旧配对')
    row = dict(instance=name, mode=mode, node_limit=budget, status=rows[0]['status'], nodes=rows[0]['nodes'])
    for v, rs in by_variant.items():
        row[v] = dict(runs=len(rs), median_secs=statistics.median(r['wall_secs'] for r in rs), min_secs=min(r['wall_secs'] for r in rs), max_secs=max(r['wall_secs'] for r in rs), median_rss_bytes=statistics.median(r['peak_rss_bytes'] for r in rs))
    row['time_ratio'] = row['candidate']['median_secs'] / row['baseline']['median_secs']
    row['rss_ratio'] = row['candidate']['median_rss_bytes'] / row['baseline']['median_rss_bytes']
    summary.append(row)
primary = [r for r in summary if r['node_limit'] == 200000]
geo = math.exp(statistics.mean(math.log(r['time_ratio']) for r in primary)) if primary else None
passed = bool(primary) and geo <= 1.10 and all(r['time_ratio'] <= 1.25 and r['rss_ratio'] <= .5 for r in primary)
print(json.dumps(dict(checked=count, deterministic_equal=True, acceptance_passed=passed, primary_geo_time_ratio=geo, summary=summary), indent=2))
