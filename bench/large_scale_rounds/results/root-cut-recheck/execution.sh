#!/bin/bash
set -euo pipefail
set -o noclobber

phase=$1
source_root=$2
commit=$3
rounds=$4
archive=$5
python=${6:-python3.13}
refs="$rounds/results"
legacy="$refs/external-baseline.jsonl"
modern="$refs/miplib2017-baseline.jsonl"
truth="$refs/miplib2017-sources.json"
feature=root_cut_recheck
binaries="$archive/binaries"
builds="$archive/builds"
target="$archive/target"
manifest="$builds/$feature-binaries.json"

runlog() {
    local log="$1"
    shift
    {
        printf '%q ' "$@"
        printf '\n'
        "$@"
    } > "$log" 2>&1
}

test "$(git -C "$source_root" rev-parse HEAD)" = "$commit"
test -z "$(git -C "$source_root" status --porcelain --untracked-files=all)"

build_probes() {
    local kind="$1"
    local candidate_features="$feature"
    local header=(cargo build --manifest-path "$source_root/milp/Cargo.toml" --target-dir "$target" --release --example round_probe)
    test ! -e "$binaries/round_probe_${kind}_base"
    test ! -e "$binaries/round_probe_${kind}_candidate"
    if test "$kind" = diagnostics; then
        runlog "$builds/$kind-base.log" "${header[@]}" --features proof_profile
        candidate_features="$feature,proof_profile"
    else
        runlog "$builds/$kind-base.log" "${header[@]}"
    fi
    cp "$target/release/examples/round_probe" "$binaries/round_probe_${kind}_base"
    runlog "$builds/$kind-candidate.log" "${header[@]}" --features "$candidate_features"
    cp "$target/release/examples/round_probe" "$binaries/round_probe_${kind}_candidate"
    "$python" "$archive/record-build.py" "$source_root" "$rounds" "$archive" "$kind"
}

case "$phase" in
    build)
        mkdir "$binaries" "$builds"
        runlog "$builds/rss-preflight.log" env PYTHONPATH="$rounds" "$python" -c \
            'import sys; from compare import invoke, preflight_rss; preflight_rss(); print(invoke([sys.executable, "-c", "print(\"{}\")"], 10, rss=True))'
        build_probes timing
        build_probes diagnostics
        exit
        ;;
    screen)
        seconds=10
        reps=1
        ;;
    paired)
        test -f "$archive/screen/exit-codes.json"
        "$python" - "$archive/screen/exit-codes.json" <<'PY'
import json
import sys
with open(sys.argv[1]) as stream:
    codes = json.load(stream)
if codes != dict(execute=0, check=0, summary=0, gate=0):
    sys.exit('筛查未全部通过，停止晋级')
PY
        seconds=30
        reps=3
        ;;
    diagnostics)
        run="$archive/diagnostics"
        mkdir "$run"
        stem="$run/$feature-diagnostics"
        run_rc=0 summary_rc=0
        runlog "$run/run.log" env PYTHONPATH="$rounds" "$python" "$archive/root-cut-diagnostics.py" run             --binaries "$builds/diagnostics-$feature-binaries.json" --source "$source_root"             --baseline "$legacy" --output "$stem.jsonl" || run_rc=$?
        runlog "$run/summary.log" env PYTHONPATH="$rounds" "$python" "$archive/root-cut-diagnostics.py" summarize             --results "$stem.jsonl" --baseline "$legacy" --output "$stem-summary.json" || summary_rc=$?
        printf '{"execute": %s, "summary": %s}\n' "$run_rc" "$summary_rc" > "$run/exit-codes.json"
        cat "$run/exit-codes.json"
        test "$run_rc" = 0 && test "$summary_rc" = 0
        exit
        ;;
    *) printf '未知阶段 %s\n' "$phase" >&2; exit 2 ;;
esac

"$python" "$archive/record-build.py" "$source_root" "$rounds" "$archive" verify
run="$archive/$phase"
mkdir "$run"
stem="$run/$feature-$phase"
run_rc=0 check_rc=0 summary_rc=0 gate_rc=0
runlog "$run/run.log" "$python" "$rounds/run_feature_trial.py" --feature "$feature" \
    --binaries "$manifest" --source "$source_root" --legacy-baseline "$legacy" \
    --modern-baseline "$modern" --seconds "$seconds" --reps "$reps" \
    --output "$stem.jsonl" || run_rc=$?
runlog "$run/check.log" "$python" "$rounds/check_results.py" "$legacy" "$modern" \
    "$stem.jsonl" --public-truth "$truth" --output "$stem-check.json" || check_rc=$?
runlog "$run/summary.log" "$python" "$rounds/summarize_feature_trial.py" "$stem.jsonl" \
    --legacy-baseline "$legacy" --modern-baseline "$modern" --public-truth "$truth" \
    --output "$stem-summary.json" || summary_rc=$?
runlog "$run/gate.log" "$python" "$rounds/gate_feature_trial.py" "$stem-summary.json" \
    --output "$stem-gates.json" || gate_rc=$?
printf '{"execute": %s, "check": %s, "summary": %s, "gate": %s}\n' \
    "$run_rc" "$check_rc" "$summary_rc" "$gate_rc" > "$run/exit-codes.json"
cat "$run/exit-codes.json"
test "$run_rc" = 0 && test "$check_rc" = 0 && test "$summary_rc" = 0 && test "$gate_rc" = 0
