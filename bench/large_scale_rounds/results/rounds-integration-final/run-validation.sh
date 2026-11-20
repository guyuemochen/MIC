#!/bin/bash
set -euo pipefail
set -o noclobber

repo=$1
root=$2
commit=$3
python=${4:-python3.13}
source="$root/source"
evidence="$root/evidence"
features=binary_mutex_row_fixing,sparse_phase1_recovery

mkdir "$evidence"
cp "$0" "$evidence/run-validation.sh"
git -C "$repo" worktree add --detach "$source" "$commit" > "$evidence/worktree.log" 2>&1
test "$(git -C "$source" rev-parse HEAD)" = "$commit"
test -z "$(git -C "$source" status --porcelain --untracked-files=all)"
cd "$source"

runlog() {
    local name="$1" result=0
    shift
    {
        printf '命令'
        printf ' %q' "$@"
        printf '\n'
        "$@"
    } > "$evidence/$name.log" 2>&1 || result=$?
    printf '%s\t%s\n' "$name" "$result" >> "$evidence/status.tsv"
    printf '%s %s\n' "$name" "$result"
    awk '/test result:|^Ran [0-9]+ tests|^error:|^FAILED|^OK$/' "$evidence/$name.log"
}

runlog 01-milp-default cargo test --manifest-path milp/Cargo.toml
runlog 02-milp-experimental-combination cargo test --manifest-path milp/Cargo.toml --lib --features "$features"
runlog 03-acceptance-default cargo test --manifest-path milp/Cargo.toml --release --test acceptance -- --include-ignored
runlog 04-acceptance-experimental-combination cargo test --manifest-path milp/Cargo.toml --release --test acceptance --features "$features" -- --include-ignored
runlog 05-planner cargo test --manifest-path planner/Cargo.toml
runlog 06-ffi-default cargo test --manifest-path milp-ffi/Cargo.toml --target-dir "$source/milp-ffi/target"
runlog 07-ffi-build cargo build --manifest-path milp-ffi/Cargo.toml --release --target-dir "$source/milp-ffi/target"
if awk -F '\t' '$1=="07-ffi-build" && $2==0 {found=1} END {exit !found}' "$evidence/status.tsv"; then
    runlog 08-python-default env -u MIC_FFI_LIB "$python" -m unittest discover -s micpython/tests -t . -v
    runlog 09-python-knapsack env -u MIC_FFI_LIB "$python" micpython/examples/knapsack.py
else
    printf '%s\t%s\n' 08-python-default skipped 09-python-knapsack skipped >> "$evidence/status.tsv"
fi

rounds=bench/large_scale_rounds
refs="$rounds/results"
r12="$refs/root-cut-recheck/paired/root_cut_recheck-paired.jsonl"
legacy="$refs/external-baseline.jsonl"
modern="$refs/miplib2017-baseline.jsonl"
truth="$refs/miplib2017-sources.json"
runlog 10-bench-unit "$python" -m unittest discover -s "$rounds" -p 'test_*.py' -q
runlog 11-r12-summary "$python" "$rounds/summarize_feature_trial.py" "$r12" \
    --legacy-baseline "$legacy" --modern-baseline "$modern" --public-truth "$truth" --output "$evidence/r12-replay-summary.json"
runlog 12-r12-gate "$python" "$rounds/gate_feature_trial.py" "$evidence/r12-replay-summary.json" --output "$evidence/r12-replay-gate.json"
runlog 13-public-independent-check "$python" "$rounds/check_results.py" "$legacy" "$modern" "$r12" \
    "$refs/gurobi10/gurobi10.jsonl" --public-truth "$truth" --output "$evidence/public-independent-check.json"

test -z "$(git -C "$source" status --porcelain --untracked-files=all)"
"$python" - "$evidence" "$features" <<'PY'
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
import tomllib

archive = Path(sys.argv[1])
experimental = sys.argv[2].split(',')


def command(*args):
    return subprocess.check_output(args, text=True).strip()


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


features = tomllib.loads(Path('milp/Cargo.toml').read_text())['features']
defaults = features['default']
assert defaults == ['sparse_row_products', 'binary_coefficient_strengthening', 'root_cut_recheck']
assert all(features[name] == [] and name not in defaults for name in experimental)
status = [line.split('\t') for line in (archive / 'status.tsv').read_text().splitlines()]
assert len(status) == 13 and len(dict(status)) == 13
results = {}
for name, code in status:
    log = archive / f'{name}.log'
    if not log.exists():
        results[name] = dict(status=code)
        continue
    payload = log.read_text()
    rust_counts = re.findall(r'test result: \w+\. (\d+) passed; (\d+) failed; (\d+) ignored', payload)
    python_counts = re.findall(r'^Ran (\d+) tests?', payload, re.M)
    results[name] = {
        'exit_code': int(code), 'sha256': digest(log),
        'rust_passed': sum(int(passed) for passed, _, _ in rust_counts) if rust_counts else None,
        'rust_failed': sum(int(failed) for _, failed, _ in rust_counts) if rust_counts else None,
        'rust_ignored': sum(int(ignored) for _, _, ignored in rust_counts) if rust_counts else None,
        'python_total': int(python_counts[0]) if python_counts else None,
    }
manifest = {
    'source_commit': command('git', 'rev-parse', 'HEAD'),
    'milp_tree': command('git', 'rev-parse', 'HEAD:milp'),
    'default_features': defaults, 'combination_explicit_features': experimental,
    'cargo': command('cargo', '-vV'), 'rustc': command('rustc', '-vV'),
    'python': sys.version, 'host': platform.platform(), 'results': results,
    'support_sha256': {p.name: digest(p) for p in sorted(archive.iterdir()) if p.name not in {f'{name}.log' for name in results}},
    'passed': all(code == '0' for _, code in status),
}
with (archive / 'manifest.json').open('x', encoding='utf-8') as stream:
    json.dump(manifest, stream, ensure_ascii=False, allow_nan=False, indent=2)
    stream.write('\n')
print(json.dumps({'passed': manifest['passed'], 'result_counts': results}, ensure_ascii=False))
sys.exit(not manifest['passed'])
PY
