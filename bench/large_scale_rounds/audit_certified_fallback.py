#!/usr/bin/env python3
"""对冻结20例只做预处理/认证审计，与既有R8丢失认证记录逐例对账，不求解LP。"""

import argparse
import json
from pathlib import Path

from audit_coefficient_strengthening import FROZEN_LEGACY, HERE, inputs, run
from compare import require


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--legacy', type=Path, default=FROZEN_LEGACY)
    parser.add_argument('--legacy-public-dir', type=Path, required=True)
    parser.add_argument('--miplib-manifest', type=Path, required=True)
    parser.add_argument('--expected-r8-audit', type=Path,
                        default=HERE / 'results' / 'coefficient-strengthening-audit.jsonl')
    args = parser.parse_args()
    paths = inputs(args.legacy, args.legacy_public_dir, args.miplib_manifest)
    expected = {row['instance']: row for row in map(json.loads, args.expected_r8_audit.read_text().splitlines())}
    require(set(expected) == set(paths), 'R8审计须覆盖同一冻结20例')
    raw = run(args.probe, *(path for path, _ in paths.values())).decode().splitlines()
    records = {}
    for line in raw:
        fields = line.split('\t')
        require(len(fields) == 9, '探针须输出9列')
        path, variables, rows, before, skipped_before, after, skipped_after, lp_same, map_same = fields
        require(Path(path) not in records, '探针模型重复')
        require(lp_same == map_same == 'true', 'LP与Postsolve映射须逐位相同')
        records[Path(path)] = tuple(map(int, (variables, rows, before, skipped_before, after, skipped_after)))
    require(set(records) == {path for path, _ in paths.values()}, '探针须完整覆盖冻结输入')
    total_rescued = 0
    for name, (path, digest) in paths.items():
        variables, rows, before, skipped_before, after, skipped_after = records[path]
        reference = expected[name]
        require(digest == reference['sha256'], f'{name} 的冻结哈希变化')
        require(rows == reference['rows'], f'{name} 的预处理行数变化')
        require((before, skipped_before) == (reference['propagation_accepted_candidate'],
                reference['propagation_skipped_candidate']), f'{name} 的无回退认证计数变化')
        require((after, skipped_after) == (reference['propagation_accepted_base'],
                reference['propagation_skipped_base']), f'{name} 未恢复既有R8丢失认证行')
        rescued = after - before
        total_rescued += rescued
        print(json.dumps({'instance': name, 'sha256': digest, 'variables': variables,
                          'rows': rows, 'lp_bitwise_identical': True,
                          'postsolve_bitwise_identical': True, 'accepted_without_fallback': before,
                          'skipped_without_fallback': skipped_before, 'accepted_with_fallback': after,
                          'skipped_with_fallback': skipped_after, 'rescued_rows': rescued},
                         allow_nan=False), flush=True)
    require(total_rescued == 597, '冻结20例须实际救回既有597个认证行')


if __name__ == '__main__':
    main()
