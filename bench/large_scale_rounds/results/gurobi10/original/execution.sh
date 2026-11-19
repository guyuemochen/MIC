#!/bin/bash
set -euo pipefail
set -o noclobber

repo=$1
source=$2
commit=$3
work=$4
library=$5
python=${6:-python3.13}
refs="$repo/bench/large_scale_rounds/results"
script="$source/bench/large_scale_rounds/compare_gurobi10.py"
legacy="$refs/external-baseline.jsonl"
modern="$refs/miplib2017-baseline.jsonl"
truth="$refs/miplib2017-sources.json"

test "$(git -C "$source" rev-parse HEAD)" = "$commit"
test -z "$(git -C "$source" status --porcelain --untracked-files=all)"
test -s "$work/frozen-permission-probe.json"
run=("$python" "$script" run --library "$library" --legacy-baseline "$legacy" --modern-baseline "$modern" --output "$work/gurobi10.jsonl")
"$python" - "$source" "$library" "$legacy" "$modern" "$truth" "${run[@]}" > "$work/freeze.json" <<'PY'
import datetime
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

source, library, *refs = map(Path, sys.argv[1:6])

def git(*args):
    return subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

rounds = source / 'bench/large_scale_rounds'
value = {
    'source_commit': git('rev-parse', 'HEAD'), 'bench_tree': git('rev-parse', 'HEAD:bench/large_scale_rounds'),
    'library': str(library), 'library_sha256': digest(library), 'host': platform.platform(),
    'python': sys.version, 'command': sys.argv[6:], 'frozen_before_measurement_utc': datetime.datetime.now(datetime.UTC).isoformat(),
    'inputs_sha256': {path.name: digest(path) for path in refs},
    'tools_sha256': {name: digest(rounds / name) for name in ('compare_gurobi10.py', 'check_results.py', 'compare.py')},
}
print(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False))
PY
rc=0
{
    printf '%q ' "${run[@]}"
    printf '\n'
    "${run[@]}"
} > "$work/run.log" 2>&1 || rc=$?
printf '%s\n' "$rc" > "$work/run-exit-code.txt"
exit "$rc"
