#!/bin/bash
set -euo pipefail

repo=/Users/guyuemochen/Desktop/projects/MIC
source_root=/tmp/mic-dense-basis-refactor
screen=/tmp/mic-r11-screen
archive=/tmp/mic-r11-paired
binaries=/tmp/mic-r11-binaries
python=/opt/homebrew/bin/python3.13
rounds="$repo/bench/large_scale_rounds"
legacy="$rounds/results/external-baseline.jsonl"
modern="$rounds/results/miplib2017-baseline.jsonl"
truth="$rounds/results/miplib2017-sources.json"
feature=certified_row_fallback
screen_builds="$screen/$feature-builds.json"
screen_manifest="$screen/$feature-binaries.json"
screen_gate="$screen/$feature-screen-gates.json"

test ! -e "$archive"
mkdir "$archive"
cp "$0" "$archive/execution.sh"
cp "$screen_builds" "$screen_manifest" "$screen_gate" "$archive/"
manifest="$archive/$feature-binaries.json"

runlog() {
    local log="$1"
    shift
    {
        printf '%q ' "$@"
        printf '\n'
        "$@"
    } > "$archive/$log" 2>&1
}

runlog freeze-check.log "$python" - "$source_root" "$archive" "$binaries" "$rounds" "$feature" "$legacy" "$modern" "$truth" <<'PY'
import hashlib
import json
from pathlib import Path
import subprocess
import sys

source, archive, binaries, tools = map(Path, sys.argv[1:5])
feature = sys.argv[5]
legacy, modern, truth = map(Path, sys.argv[6:])

def git(*args):
    return subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

assert not git('status', '--porcelain', '--untracked-files=all'), '冻结工作树不干净'
commit, tree = git('rev-parse', 'HEAD'), git('rev-parse', 'HEAD:milp')
assert commit == '3d136d6905516292b5885121c5c6cf7b29d72eea'
assert tree == '23d628b71db35e3e5b846dcc66457842e5926402'
manifest = archive / f'{feature}-binaries.json'
builds = json.loads((archive / f'{feature}-builds.json').read_text())
gate = json.loads((archive / f'{feature}-screen-gates.json').read_text())
entries = json.loads(manifest.read_text())
assert sha(manifest) == builds['manifest_sha256'] == gate['source']['binary_manifest_sha256']
assert commit == builds['source_commit'] == gate['source']['source_commit']
assert tree == builds['source_tree'] == gate['source']['source_tree']
assert gate['feature'] == feature and gate['phase'] == 'screen'
assert gate['time_limit_secs'] == 10 and gate['repetitions'] == 1
assert gate['overall_passed'] is True and gate['eligible_for_paired'] is True
assert set(entries) == {'base', 'candidate'}
for variant, features in (('base', []), ('candidate', [feature])):
    expected = binaries / f'round_probe_{variant}'
    assert entries[variant] == {'binary': str(expected), 'features': features}
    assert sha(expected) == builds['builds'][variant]['binary_sha256'] == gate['variants'][variant]['binary_sha256']
print(json.dumps({
    'source_commit': commit, 'source_tree': tree, 'repository_tree': git('rev-parse', 'HEAD^{tree}'),
    'binary_manifest_sha256': sha(manifest),
    'binaries_sha256': {variant: sha(Path(entry['binary'])) for variant, entry in entries.items()},
    'inputs_sha256': {str(path): sha(path) for path in (legacy, modern, truth)},
    'tools_sha256': {name: sha(tools / name) for name in ('run_feature_trial.py', 'check_results.py',
                                                        'summarize_feature_trial.py', 'gate_feature_trial.py')},
    'time_limit_secs': 30, 'repetitions': 3, 'expected_rows': 120,
}, indent=2))
PY

runlog rss-preflight.log env PYTHONPATH="$rounds" "$python" -c 'from compare import preflight_rss; preflight_rss(); print("RSS 预检通过")'
stem="$archive/$feature-paired"
refs=(--legacy-baseline "$legacy" --modern-baseline "$modern")
run_rc=0
runlog paired-run.log "$python" "$rounds/run_feature_trial.py" --feature "$feature" \
    --binaries "$manifest" --source "$source_root" "${refs[@]}" --seconds 30 --reps 3 \
    --output "$stem.jsonl" || run_rc=$?
test -s "$stem.jsonl"

check_rc=0
runlog paired-check.log "$python" "$rounds/check_results.py" "$legacy" "$modern" "$stem.jsonl" \
    --public-truth "$truth" --output "$stem-check.json" || check_rc=$?
summary_rc=0
runlog paired-summary.log "$python" "$rounds/summarize_feature_trial.py" "$stem.jsonl" "${refs[@]}" \
    --public-truth "$truth" --output "$stem-summary.json" || summary_rc=$?
gate_rc=0
runlog paired-gate.log "$python" "$rounds/gate_feature_trial.py" "$stem-summary.json" \
    --output "$stem-gates.json" || gate_rc=$?
printf '执行 %s 独立复核 %s 汇总 %s 门禁 %s\n' "$run_rc" "$check_rc" "$summary_rc" "$gate_rc" \
    > "$archive/paired-exit-codes.txt"
cat "$archive/paired-exit-codes.txt"
test "$run_rc" = 0 && test "$check_rc" = 0 && test "$summary_rc" = 0 && test "$gate_rc" = 0
