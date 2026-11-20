#!/usr/bin/env python3
"""只为默认 Gurobi 运行的模糊状态补测关闭 DualReductions 的诊断。"""

import gzip
import json
from pathlib import Path
import sys

SCRIPT = Path(__file__).resolve()
OUT = SCRIPT.parent
sys.path.insert(0, str(SCRIPT.parents[2]))

from compare import require, sha
from compare_all_timing import capture, now, parse_outcome, read_model, write_json
from compare_drayage_timing import observe_gurobi
import compare_gurobi10 as grb


def diagnostic_parameters(library, env, seconds, nodes):
    ints, reals = original_parameters(library, env, seconds, nodes)
    grb.check(library, env, library.GRBsetintparam(env, b'DualReductions', 0), '设置 DualReductions')
    return dict(ints, DualReductions=0), reals


original_parameters = grb.configure


def main():
    freeze = json.loads((OUT / 'freeze.json').read_text(encoding='utf-8'))
    completion = json.loads((OUT / 'completed.json').read_text(encoding='utf-8'))
    timing = OUT / 'timings.jsonl.gz'
    require(sha(timing) == completion['timings_sha256'], '主矩阵原始档摘要不同')
    library = Path(freeze['library'])
    require(sha(library) == freeze['library_sha256'], '动态库摘要与主矩阵不同')
    with gzip.open(timing, 'rt', encoding='utf-8') as stream:
        source = [row for line in stream if (row := json.loads(line))['solver'] == 'gurobi'
                  and row['raw_status'] == 'INF_OR_UNBD']
    results = []
    for row in source:
        command = [sys.executable, str(SCRIPT), '_worker', str(library), str(OUT / row['model']),
                   str(freeze['time_limit_secs']), str(freeze['node_limit'])]
        raw = capture(command, freeze['time_limit_secs'] + 60)
        outcome = parse_outcome(raw, freeze['gurobi_version'])
        results.append(dict(case_id=row['case_id'], instance=row['instance'], model=row['model'],
                            model_sha256=row['model_sha256'], solver='gurobi',
                            primary_record={key: row[key] for key in ('status', 'raw_status', 'elapsed_secs',
                                                                      'parameters', 'version')},
                            diagnostic=outcome, raw=raw))
    write_json(OUT / 'gurobi-ambiguity-diagnostics.json', dict(
        schema=1, recorded_at_utc=now(), command=[sys.executable, *sys.argv], script_sha256=sha(SCRIPT),
        timings_sha256=completion['timings_sha256'], library_sha256=freeze['library_sha256'],
        changed_parameter={'DualReductions': 0}, count=len(results), results=results))
    return any(row['diagnostic']['status'] == 'error' for row in results)


if __name__ == '__main__':
    if len(sys.argv) == 6 and sys.argv[1] == '_worker':
        grb.configure = diagnostic_parameters
        library_path, model_path = map(Path, sys.argv[2:4])
        observed = grb.solve(library_path, read_model(model_path), int(sys.argv[4]), int(sys.argv[5]), observe_gurobi)
        print(json.dumps(observed, ensure_ascii=False, allow_nan=False))
    else:
        sys.exit(main())
