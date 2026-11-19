#!/bin/sh
set -eu
archive=${1:?需要归档目录}
source=${2:?需要冻结源码树}
expected=b4b2d1000ec34ae68b4e91b53c28442892e2c7f1
rounds="$source/bench/large_scale_rounds"
test "$(git -C "$source" rev-parse HEAD)" = "$expected"
test -z "$(git -C "$source" status --porcelain --untracked-files=all)"
PYTHONDONTWRITEBYTECODE=1 python3.13 "$rounds/run_binary_mutex_point.py" --source "$source" \
  --binaries "$archive/binaries.json" \
  --reference "$rounds/results/miplib2017-baseline.jsonl" \
  --public-truth "$rounds/results/miplib2017-sources.json" \
  --output "$archive/point.jsonl" --gate "$archive/point-gate.json"
