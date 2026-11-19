#!/bin/sh
set -u
logs=/tmp/mic-r12-verify
run() {
    name=$1
    shift
    printf '开始 %s\n' "$name"
    "$@" > "$logs/$name.log" 2>&1
    rc=$?
    printf '%s\t%s\n' "$name" "$rc" | tee -a "$logs/status.tsv"
    rg 'test result:|error:|^warning:' "$logs/$name.log" | tail -14
}
run 03-default-unit cargo test --manifest-path milp/Cargo.toml
run 04-r12-unit cargo test --manifest-path milp/Cargo.toml --features root_cut_recheck
run 05-planner cargo test --manifest-path planner/Cargo.toml
run 06-default-release-acceptance cargo test --manifest-path milp/Cargo.toml --release --test acceptance -- --include-ignored
run 07-r12-release-acceptance cargo test --manifest-path milp/Cargo.toml --release --features root_cut_recheck --test acceptance -- --include-ignored
run 08-r12-profile-unit cargo test --manifest-path milp/Cargo.toml --lib --features root_cut_recheck,proof_profile
run 09-r12-snapshot-legacy cargo test --manifest-path milp/Cargo.toml --lib --features root_cut_recheck,root_snapshot_restore local_cuts_fire_deterministically_and_toggle_restores
