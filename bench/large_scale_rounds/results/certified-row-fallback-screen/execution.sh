#!/bin/bash
set -euo pipefail

repo=/Users/guyuemochen/Desktop/projects/MIC
source_root=/tmp/mic-dense-basis-refactor
archive=/tmp/mic-r11-screen
binaries=/tmp/mic-r11-binaries
target=/tmp/mic-r11-target
python=/opt/homebrew/bin/python3.13
rounds="$repo/bench/large_scale_rounds"
legacy="$rounds/results/external-baseline.jsonl"
modern="$rounds/results/miplib2017-baseline.jsonl"
truth="$rounds/results/miplib2017-sources.json"
feature=certified_row_fallback

test "$(git -C "$source_root" rev-parse HEAD)" = 3d136d6905516292b5885121c5c6cf7b29d72eea
test "$(git -C "$source_root" rev-parse HEAD:milp)" = 23d628b71db35e3e5b846dcc66457842e5926402
test -z "$(git -C "$source_root" status --porcelain --untracked-files=all)"
test ! -e "$archive" && test ! -e "$binaries"
mkdir "$archive" "$binaries"
cp "$0" "$archive/execution.sh"
cd "$source_root"

runlog() {
    local log="$1"
    shift
    {
        printf '%q ' "$@"
        printf '\n'
        "$@"
    } > "$archive/$log" 2>&1
}

manifest=milp/Cargo.toml
runlog candidate-debug.log cargo test --manifest-path "$manifest" --features "$feature" --target-dir "$target"
runlog candidate-acceptance.log cargo test --manifest-path "$manifest" --features "$feature" --target-dir "$target" --release --test acceptance -- --include-ignored

runlog base-build.log cargo build --manifest-path "$manifest" --target-dir "$target" --release --example round_probe
cp "$target/release/examples/round_probe" "$binaries/round_probe_base"
runlog candidate-build.log cargo build --manifest-path "$manifest" --features "$feature" --target-dir "$target" --release --example round_probe
cp "$target/release/examples/round_probe" "$binaries/round_probe_candidate"

"$python" - "$source_root" "$archive" "$binaries" "$feature" <<'PY'
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import tomllib

source, archive, binaries = map(Path, sys.argv[1:4])
feature = sys.argv[4]
with (source / 'milp/Cargo.toml').open('rb') as stream:
    declared = tomllib.load(stream)['features']

def output(*command):
    return subprocess.check_output(command, text=True).strip()

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def closure(names):
    enabled = set(names)
    pending = list(names)
    while pending:
        for implied in declared[pending.pop()]:
            if implied not in enabled:
                enabled.add(implied)
                pending.append(implied)
    return enabled

defaults = declared['default']
assert feature not in closure(defaults)
assert 'binary_coefficient_strengthening' in closure(defaults)
assert declared[feature] == ['binary_coefficient_strengthening']
entries = {}
builds = {}
for variant, explicit in (('base', []), ('candidate', [feature])):
    binary = binaries / f'round_probe_{variant}'
    enabled = closure([*defaults, *explicit])
    assert enabled == set(defaults) | set(explicit)
    entries[variant] = {'binary': str(binary), 'features': explicit}
    builds[variant] = dict(
        entries[variant], binary_sha256=sha(binary), cargo_default_features=defaults,
        cargo_dependencies_of_explicit_features=sorted(closure(explicit) - set(explicit)),
        effective_features=sorted(enabled), build_log=f'{variant}-build.log',
        build_log_sha256=sha(archive / f'{variant}-build.log'),
    )
manifest = archive / 'certified_row_fallback-binaries.json'
with manifest.open('x', encoding='utf-8') as stream:
    json.dump(entries, stream, indent=2)
    stream.write('\n')
with (archive / 'certified_row_fallback-builds.json').open('x', encoding='utf-8') as stream:
    json.dump({
        'source': str(source), 'source_commit': output('git', '-C', str(source), 'rev-parse', 'HEAD'),
        'source_tree': output('git', '-C', str(source), 'rev-parse', 'HEAD:milp'),
        'repository_tree': output('git', '-C', str(source), 'rev-parse', 'HEAD^{tree}'),
        'manifest_sha256': sha(manifest), 'builds': builds, 'host': platform.platform(),
        'cargo': output('cargo', '-vV'), 'rustc': output('rustc', '-vV'), 'python': sys.version,
    }, stream, ensure_ascii=False, indent=2)
    stream.write('\n')
PY

runlog rss-preflight.log env PYTHONPATH="$rounds" "$python" -c 'from compare import preflight_rss; preflight_rss(); print("RSS 预检通过")'
stem="$archive/certified_row_fallback-screen"
refs=(--legacy-baseline "$legacy" --modern-baseline "$modern")
run_rc=0
runlog screen-run.log "$python" "$rounds/run_feature_trial.py" --feature "$feature" \
    --binaries "$archive/certified_row_fallback-binaries.json" --source "$source_root" \
    "${refs[@]}" --seconds 10 --reps 1 --output "$stem.jsonl" || run_rc=$?
test -s "$stem.jsonl"

check_rc=0
runlog screen-check.log "$python" "$rounds/check_results.py" \
    "$legacy" "$modern" "$stem.jsonl" --public-truth "$truth" --output "$stem-check.json" || check_rc=$?
summary_rc=0
runlog screen-summary.log "$python" "$rounds/summarize_feature_trial.py" "$stem.jsonl" \
    "${refs[@]}" --public-truth "$truth" --output "$stem-summary.json" || summary_rc=$?
gate_rc=0
runlog screen-gate.log "$python" "$rounds/gate_feature_trial.py" "$stem-summary.json" \
    --output "$stem-gates.json" || gate_rc=$?
printf '执行 %s 独立复核 %s 汇总 %s 门禁 %s\n' "$run_rc" "$check_rc" "$summary_rc" "$gate_rc" \
    > "$archive/screen-exit-codes.txt"
cat "$archive/screen-exit-codes.txt"
test "$run_rc" = 0 && test "$check_rc" = 0 && test "$summary_rc" = 0 && test "$gate_rc" = 0
