#!/usr/bin/env python3
"""用外部SciPy-HiGHS核对生成集；依赖仅用于独立实验，不进入MIC运行依赖。"""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from python_bench import parse_milp, run_scipy
import scipy

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--directory', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--seconds', type=int, default=15)
a = p.parse_args()
with a.output.open('w') as out:
    for path in sorted(a.directory.glob('*.milp')):
        build, elapsed, status, obj = run_scipy(parse_milp(path), a.seconds)
        row = dict(instance=path.stem, solver='SciPy-HiGHS', scipy_version=scipy.__version__,
                   threads=1, time_limit=a.seconds, status=status, objective=obj,
                   build_secs=build, solve_secs=elapsed)
        out.write(json.dumps(row) + '\n'); out.flush()
        print(path.stem, status, obj, flush=True)
