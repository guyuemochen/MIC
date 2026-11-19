#!/usr/bin/env python3
"""核对固定二十例的生产预处理覆盖。仅运行缩放及预处理。"""

import argparse
import json
from pathlib import Path
import re
import subprocess

from compare import require, sha
from run_root_work import AUTO_CASES, MODERN, baseline_rows


SUMMARY = re.compile(r'(\d+)列/(\d+)行 → (\d+)列/(\d+)行 固定(\d+) 聚合(\d+)')


def measure(probe, instance):
    observed = subprocess.run([probe, instance], check=True, capture_output=True, text=True, timeout=20)
    match = SUMMARY.search(observed.stderr)
    require(match is not None, f'缺少预处理摘要：{instance}')
    orig_cols, orig_rows, cols, rows, fixed, aggregates = map(int, match.groups())
    assignments = {}
    for line in observed.stdout.splitlines():
        name, j, value = line.split('\t')
        index = int(j)
        require(index not in assignments and 0 <= index < orig_cols, '固定的原列下标重复或越界')
        assignments[index] = (name, float(value))
    require(len(assignments) == fixed, f'固定摘要与列表不一致：{instance}')
    return dict(orig_cols=orig_cols, orig_rows=orig_rows, cols=cols, rows=rows,
                fixed=fixed, aggregates=aggregates), assignments


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['base', 'candidate', 'legacy', 'modern', 'output']:
        parser.add_argument(f'--{name}', type=Path, required=True)
    args = parser.parse_args()
    references = baseline_rows(args.legacy, AUTO_CASES) | baseline_rows(args.modern, MODERN)
    evidence = dict(base_sha256=sha(args.base), candidate_sha256=sha(args.candidate), cases={})
    for case in AUTO_CASES + MODERN:
        reference = references[case]
        instance = Path(reference['command'][2])
        checksum = sha(instance)
        require(checksum == reference['milp_sha256'], f'{case} 与冻结原模型哈希不符')
        base, old = measure(args.base, instance)
        candidate, new = measure(args.candidate, instance)
        require((base['orig_cols'], base['orig_rows']) == (candidate['orig_cols'], candidate['orig_rows']),
                f'{case} 两组原模型规模不同')
        added = new.keys() - old.keys()
        lost = old.keys() - new.keys()
        replaced = [j for j in old.keys() & new.keys() if old[j] != new[j]]
        evidence['cases'][case] = dict(milp_sha256=checksum, base=base, candidate=candidate,
                                       added=len(added), lost=len(lost), replaced=len(replaced),
                                       net_fixed=len(added) - len(lost),
                                       added_zero=sum(new[j][1] == 0.0 for j in added))
    with args.output.open('x', encoding='utf-8') as output:
        json.dump(evidence, output, ensure_ascii=False, indent=2)
        output.write('\n')


if __name__ == '__main__':
    main()
