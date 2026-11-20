#!/usr/bin/env python3
"""独立检查本轮生成模型的解、目标、界和自动选择；不调用MIC求解代码。"""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from python_bench import parse_milp

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--directory', type=Path, required=True)
p.add_argument('--results', type=Path, nargs='+', required=True)
p.add_argument('--reference', type=Path)
a = p.parse_args()
references = {} if not a.reference else {r['instance']: r for r in map(json.loads, a.reference.read_text().splitlines())}
manifest = json.loads((Path(__file__).parent / 'manifest.json').read_text())
hashes = {r['file']: r['sha256'] for r in manifest['instances']}
models = {}
count = 0
optimal_objectives = {}
completed = {}
for result in a.results:
    for line in result.read_text().splitlines():
        r = json.loads(line); name = r['instance']
        assert r['status'] != 'driver_error', r
        if name not in models:
            path = a.directory / (name + '.milp')
            assert hashlib.sha256(path.read_bytes()).hexdigest() == hashes[path.name], name
            models[name] = parse_milp(path)
        m = models[name]; idx = {v: j for j, (v, _) in enumerate(m.vars)}
        def value(terms, xs): return sum(c * xs[idx[v]] for c, v in terms)
        selected = (200 <= len(m.vars) <= 512 and 2 <= len(m.cons) <= 16 and m.sense == 'max'
                    and all(t == 'B' for _, t in m.vars) and len(m.obj) == len(m.vars)
                    and all(c > 0 for c, _ in m.obj)
                    and all(rel == '<=' and rhs > 0 and terms and all(c >= 0 for c, _ in terms) for _, terms, rel, rhs in m.cons))
        if r['mode'] == 'auto':
            assert r['profile'] == ('PackingHybrid' if selected else 'Existing'), r
        assert r.get('values_error') is None, r
        obj, xs = r.get('objective'), r.get('values')
        if obj is not None:
            assert math.isfinite(obj) and xs is not None and len(xs) == len(m.vars), r
            for (_, kind), x in zip(m.vars, xs):
                assert math.isfinite(x) and x >= -1e-6, r
                if kind == 'B': assert x <= 1 + 1e-6, r
                if kind != 'C': assert abs(x - round(x)) <= 1e-6, r
            for _, terms, rel, rhs in m.cons:
                lhs = value(terms, xs); tol = 1e-7 * max(1, abs(rhs))
                assert (lhs <= rhs + tol if rel == '<=' else lhs >= rhs - tol if rel == '>=' else abs(lhs-rhs) <= tol), (name, lhs, rhs)
            assert abs(value(m.obj, xs)-obj) <= 1e-7 * max(1, abs(obj)), name
        if r['status'] == 'Optimal':
            assert obj is not None and not r['auto_degraded'], r
            if name in optimal_objectives:
                assert abs(obj-optimal_objectives[name]) <= 1e-7 * max(1, abs(obj)), r
            optimal_objectives[name] = obj
            # 无墙钟截断时，重复求解的工作量及解向量应完全一致。
            key = (name, r['mode'])
            signature = (obj, xs, r['nodes'], r['lp_iterations'])
            if key in completed: assert signature == completed[key], (name, r['mode'], '非确定性')
            completed[key] = signature
        if r['auto_degraded']: assert r['status'] == 'LimitReached', r
        ref = references.get(name)
        if ref and ref['status'] == 'optimal':
            optimum = ref['objective']; tol = 1e-7 * max(1, abs(optimum))
            assert r['status'] not in ('Infeasible', 'Unbounded'), r
            if r['status'] == 'Optimal': assert abs(obj-optimum) <= tol, r
            if r.get('bound') is not None:
                assert (r['bound'] >= optimum-tol if m.sense == 'max' else r['bound'] <= optimum+tol), r
        if obj is not None and abs(obj) > 1e-10 and r.get('gap_percent') is not None:
            diff = r['bound']-obj if m.sense == 'max' else obj-r['bound']
            assert abs(r['gap_percent']-max(0, diff/abs(obj)*100)) <= 1e-6, r
        count += 1
print(json.dumps({'checked': count, 'models': len(models), 'reference_models': len(references), 'passed': True}))
