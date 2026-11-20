#!/usr/bin/env python3
"""公开算例的有界串行对比；输入目录包含同名原始 MPS 与经审核转换的 .milp。"""
import argparse
import json
import math
import subprocess
import sys
from pathlib import Path


def highs(path, seconds):
    import highspy
    h = highspy.Highs()
    for k, v in [('output_flag', False), ('threads', 1), ('time_limit', seconds), ('mip_rel_gap', 0.0), ('mip_abs_gap', 0.0)]:
        assert h.setOptionValue(k, v) == highspy.HighsStatus.kOk
    assert h.readModel(str(path)) == highspy.HighsStatus.kOk
    h.run()
    info = h.getInfo()
    return dict(status=h.modelStatusToString(h.getModelStatus()), objective=info.objective_function_value if h.getSolution().value_valid else None,
                bound=info.mip_dual_bound, gap=info.mip_gap, nodes=info.mip_node_count,
                wall_secs=h.getRunTime(), version=h.version())


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--mic', type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--seconds', type=float, default=10)
    p.add_argument('--worker')
    a = p.parse_args()
    if a.worker:
        print(json.dumps(highs(a.directory / (a.worker + '.mps'), a.seconds)))
        sys.exit(0)
    with a.output.open('w') as out:
        for rep in range(3):
            for name in ['p0033', 'p0201', 'enigma', 'stein27']:
                for solver in ['mic', 'highs']:
                    cmd = ([str(a.mic), str(a.directory / (name + '.milp')), '--threads', '1', '--time-limit-secs', str(int(a.seconds))]
                           if solver == 'mic' else [sys.executable, __file__, '--directory', str(a.directory), '--seconds', str(a.seconds), '--worker', name])
                    try:
                        proc = subprocess.run(cmd, text=True, capture_output=True, timeout=a.seconds + 60, check=True)
                        row = json.loads(proc.stdout.strip().splitlines()[-1])
                    except Exception as e:
                        row = {'status': 'driver_error', 'error': str(e)}
                    row = {k: (None if isinstance(v, float) and not math.isfinite(v) else v) for k, v in row.items()}
                    row.update(instance=name, solver=solver, rep=rep)
                    out.write(json.dumps(row) + '\n'); out.flush()
                    print(row, flush=True)
