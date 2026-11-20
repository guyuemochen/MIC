# 最终整合验收

干净独立工作树 `bf0d21d44c3d61d65366bc1aea1ade6f97bcb6c8` 的 13 项命令全部退出 0。
MILP 子树为 `c6871fb661e882125ab35f4087533aeeba35ad9c`。默认包含稀疏行乘积、二元系数加强
与根割复查。组合组显式开启 `binary_mutex_row_fixing,sparse_phase1_recovery`。

| 验收项目 | 通过项 | 忽略项 |
| --- | ---: | ---: |
| 默认 MILP 全量 | 539 | 6 |
| 两项候选组合完整库 | 508 | 4 |
| 默认 release acceptance | 16 | 0 |
| 两项候选组合 release acceptance | 16 | 0 |
| planner | 127 | 2 |
| 默认 FFI | 27 | 0 |
| 默认 Python | 65 | 0 |
| 基准门禁 | 9 | 0 |

默认 FFI release 动态库用于 Python 装载、ABI 版本和背包示例检查。R12 正式门禁与 Gurobi
交叉核验也复现，后者覆盖 272 条记录、233 个可行解、126 条有外部最优值可对照的 MIC 界。
`manifest.json` 登记源码、工具链、解压后的原日志和输出摘要，`status.tsv` 保留逐项退出码。
`log-archive.json` 登记 gzip 文件及其解压后的原件摘要。

重新运行整合矩阵需要从仓库根目录创建干净工作树和新日志目录。

```sh
refs=bench/large_scale_rounds/results
repo=$(git rev-parse --show-toplevel)
out=$(mktemp -d)
/bin/bash "$refs/rounds-integration-final/run-validation.sh" "$repo" "$out" \
  bf0d21d44c3d61d65366bc1aea1ade6f97bcb6c8 "$(command -v python3.13)"
```

以下命令从已有结果重算 R12 门禁、外部对照和两项候选的定点证据，并核对归档摘要。

```sh
rounds=bench/large_scale_rounds
refs="$rounds/results"
legacy="$refs/external-baseline.jsonl"
modern="$refs/miplib2017-baseline.jsonl"
truth="$refs/miplib2017-sources.json"
r12="$refs/root-cut-recheck/paired/root_cut_recheck-paired.jsonl"
out=$(mktemp -d)
python3.13 "$rounds/summarize_feature_trial.py" "$r12" --legacy-baseline "$legacy" \
  --modern-baseline "$modern" --public-truth "$truth" --output "$out/r12-summary.json"
python3.13 "$rounds/gate_feature_trial.py" "$out/r12-summary.json" --output "$out/r12-gate.json"
python3.13 "$rounds/check_results.py" "$legacy" "$modern" "$r12" "$refs/gurobi10/gurobi10.jsonl" \
  --public-truth "$truth" --output "$out/public-check.json"
(cd "$refs/bnatt-phase1-recovery" && shasum -a 256 -c SHA256SUMS)
python3.13 "$refs/bnatt-phase1-recovery/verify.py"
(cd "$refs/binary-mutex-row-fixing" && shasum -a 256 -c SHA256SUMS)
PYTHONPATH="$rounds" python3.13 - "$refs" <<'PY_REPLAY'
import json
from pathlib import Path
import sys
from check_results import inspect_results
from run_binary_mutex_point import CASES, VARIANTS, evaluate

refs = Path(sys.argv[1])
archive = refs / 'binary-mutex-row-fixing'
point = archive / 'point.jsonl'
records = list(map(json.loads, point.read_text().splitlines()))
rows = {(row['case'], row['variant']): row for row in records}
assert len(records) == 4 and set(rows) == {(case, variant) for case in CASES for variant in VARIANTS}
check = inspect_results([point, refs / 'miplib2017-baseline.jsonl'], refs / 'miplib2017-sources.json')
actual = evaluate(rows, check)
expected = json.loads((archive / 'point-gate.json').read_text())
for report in (actual, expected):
    report['correctness']['public_truth'].pop('manifest')
assert actual == expected and actual['passed']
print(json.dumps(actual, ensure_ascii=False, indent=2))
PY_REPLAY
```

`bnatt-phase1-recovery/README.md` 还记录了重采现场后核对原始字节的命令。
`binary-mutex-row-fixing.md` 记录了重新构建探针并经精确有理活动度核对固定集的命令。
