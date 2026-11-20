#!/usr/bin/env python3
r"""审计两套缩放后预处理导出，精确验证变更系数并统计生产传播认证。

在仓库根目录分别构建基线和候选，各用自己的 target 目录。候选的 Cargo
命令追加 --features binary_coefficient_strengthening，基线使用默认配置。

    cargo build --manifest-path milp/Cargo.toml --release --lib --example round_probe --target-dir "$target"
    rustc --edition=2021 bench/large_scale_rounds/coefficient_certification_probe.rs \
        --extern milp="$target/release/libmilp.rlib" -L "dependency=$target/release/deps" \
        -O -o "$target/release/coefficient_certification_probe"

将两个 release/examples/round_probe 和两个认证探针分别传给对应参数。
历史公开模型来自 fetch_public.py，新公开模型的清单来自 fetch_miplib2017.py。
"""

import argparse
import json
import math
import struct
import subprocess
from fractions import Fraction
from pathlib import Path

from compare import require, sha


HERE = Path(__file__).resolve().parent
FROZEN_LEGACY = HERE / 'results' / 'external-baseline.jsonl'


def inputs(legacy, legacy_public, modern):
    paths = {}
    with legacy.open(encoding='utf-8') as file:
        for line in file:
            row = json.loads(line)
            if row['solver'] != 'mic' or row['rep'] != 0:
                continue
            source = Path(row['command'][2])
            directory = HERE.parent / 'instances' if source.parent.name == 'instances' else legacy_public
            paths[row['instance']] = (directory / source.name, row['milp_sha256'])
    require(len(paths) == 12, '历史基线须含 12 个模型')
    for row in json.loads(modern.read_text(encoding='utf-8'))['instances']:
        if row['included']:
            paths[row['instance']] = (modern.parent / row['milp'], row['milp_sha256'])
    require(len(paths) == 20, '完整矩阵须含 20 个模型')
    for name, (path, expected_digest) in paths.items():
        require(sha(path) == expected_digest, f'{name} .milp 与冻结哈希不同')
    return paths


def run(binary, *args):
    return subprocess.check_output([binary.resolve(), *args], timeout=180)


def certificates(binary, paths):
    lines = run(binary, *(path for path, _ in paths.values())).decode().splitlines()
    counts = {}
    for line in lines:
        path, accepted, skipped = line.split('\t')
        counts[Path(path)] = int(accepted), int(skipped)
    require(len(lines) == len(paths) and set(counts) == {path for path, _ in paths.values()},
            f'{binary} 的认证结果不完整或含重复项')
    return counts


def same(left, right):
    if type(left) is not type(right):
        return False
    if isinstance(left, float):
        return struct.pack('!d', left) == struct.pack('!d', right)
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(same(value, right[key]) for key, value in left.items())
    if isinstance(left, list):
        return len(left) == len(right) and all(same(a, b) for a, b in zip(left, right))
    return left == right


def compare(name, baseline, candidate):
    base_rows = baseline.pop('rows')
    cand_rows = candidate.pop('rows')
    require(same(baseline, candidate), f'{name} 的非约束数据变化')
    require(len(base_rows) == len(cand_rows), f'{name} 的约束行数变化')
    variables = baseline['variables']
    modified_rows = modified_coefficients = mixed_rows = 0
    for index, (old, new) in enumerate(zip(base_rows, cand_rows)):
        if same(old, new):
            continue
        label = f'{name} 第 {index} 行'
        require(old['op'] == new['op'] and same(old['rhs'], new['rhs']), f'{label} 的关系或右端变化')
        require(old['op'] in ('<=', '>='), f'{label} 为等式')
        sign = 1 if old['op'] == '<=' else -1
        require(len(old['terms']) == len(new['terms']), f'{label} 的项数变化')
        modified_rows += 1
        mixed_rows += any(variables[j]['kind'] != 'binary' for j, _ in old['terms'])
        rhs = sign * Fraction(old['rhs'])
        maximum = Fraction()
        for j, coeff in old['terms']:
            normalized = sign * Fraction(coeff)
            bound = variables[j]['ub' if normalized >= 0 else 'lb']
            require(bound is not None and math.isfinite(bound), f'{label} 第 {j} 列无有限活动界')
            maximum += normalized * Fraction(bound)
        require(maximum > rhs, f'{label} 的原始活动上界不大于右端')
        for (j, old_coeff), (k, new_coeff) in zip(old['terms'], new['terms']):
            require(j == k, f'{label} 的列下标变化')
            if same(old_coeff, new_coeff):
                continue
            var = variables[j]
            require((var['kind'], var['lb'], var['ub']) == ('binary', 0.0, 1.0),
                    f'{label} 第 {j} 列并非自由二元列')
            old_weight = sign * Fraction(old_coeff)
            new_weight = sign * Fraction(new_coeff)
            # 这列取 1 时，其他列各自取旧行最大项，新行仍须成立。
            require(old_weight < new_weight < 0 and maximum + new_weight <= rhs,
                    f'{label} 第 {j} 列的新权重未通过精确有理数核验')
            modified_coefficients += 1
    return modified_rows, modified_coefficients, mixed_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--baseline-certification', type=Path, required=True)
    parser.add_argument('--candidate-certification', type=Path, required=True)
    parser.add_argument('--legacy', type=Path, default=FROZEN_LEGACY)
    parser.add_argument('--legacy-public-dir', type=Path, required=True)
    parser.add_argument('--miplib-manifest', type=Path, required=True)
    args = parser.parse_args()

    paths = inputs(args.legacy, args.legacy_public_dir, args.miplib_manifest)
    base_certified = certificates(args.baseline_certification, paths)
    cand_certified = certificates(args.candidate_certification, paths)
    for name, (path, digest) in paths.items():
        raw_base = run(args.baseline, 'export-presolved', path)
        raw_cand = run(args.candidate, 'export-presolved', path)
        require(raw_cand == run(args.candidate, 'export-presolved', path), f'{name} 候选重复导出不同')
        baseline = json.loads(raw_base)
        candidate = json.loads(raw_cand)
        total_rows = len(baseline['rows'])
        changed, coefficients, mixed = compare(name, baseline, candidate)
        print(json.dumps({
            'instance': name, 'sha256': digest, 'rows': total_rows,
            'changed_rows': changed, 'changed_coefficients': coefficients, 'changed_mixed_rows': mixed,
            'propagation_accepted_base': base_certified[path][0],
            'propagation_accepted_candidate': cand_certified[path][0],
            'propagation_skipped_base': base_certified[path][1],
            'propagation_skipped_candidate': cand_certified[path][1],
        }, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
