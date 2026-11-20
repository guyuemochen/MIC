#!/usr/bin/env python3
"""复核根阶段对照，并仅比较共同完成最优性证明的耗时。"""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
from statistics import median
import sys

from check_results import close, finite, inspect_results
from compare import require
from run_root_work import CASES, FEATURES, NODE_LIMIT, schedule, settings, validate_variants


PROOFS = {'optimal', 'infeasible', 'unbounded'}


def load_runs(path):
    with path.open(encoding='utf-8') as stream:
        rows = [json.loads(line) for line in stream]
    require(rows, f'{path} 没有求解记录')
    variants = rows[0]['requested_variants']
    reps, seconds = rows[0]['repetitions'], rows[0]['time_limit_secs']
    validate_variants(variants)
    require(type(reps) is int and reps > 0 and type(seconds) is int and seconds > 0,
            '重复次数或时间预算非法')
    expected = list(schedule(variants, reps))
    require([tuple(row[key] for key in ('instance', 'rep', 'variant', 'order')) for row in rows] == expected,
            f'{path} 记录缺失、重复或顺序与冻结矩阵不同')
    source = (rows[0]['source_commit'], rows[0]['source_tree'], rows[0]['binary_manifest_sha256'],
              rows[0]['host'])
    builds = {}
    for row in rows:
        name, variant = row['instance'], row['variant']
        command = row['command']
        options = dict(zip(command[3::2], command[4::2], strict=True))
        require(row['solver'] == 'mic' and row['requested_variants'] == variants and row['repetitions'] == reps
                and (row['threads'], row['time_limit_secs'], row['node_limit']) == (1, seconds, NODE_LIMIT)
                and command[1] == 'solve' and len(options) * 2 == len(command) - 3
                and options == settings(seconds) and row['features'] == list(FEATURES[variant])
                and 'proof_profile' not in row, f'{name}/{variant} 的设置与冻结矩阵不同')
        require((row['source_commit'], row['source_tree'], row['binary_manifest_sha256'], row['host']) == source,
                f'{name}/{variant} 的源码、二进制清单或运行环境与同组结果不同')
        build = (command[0], row['binary_sha256'])
        require(row['version'] == row['binary_sha256'] and builds.setdefault(variant, build) == build,
                f'{name}/{variant} 的二进制哈希与同组结果不同')
        if row.get('error_stage') != 'export':
            require(row['export_model_sha256'] == row['model_sha256'],
                    f'{name}/{variant} 未核对新导出模型')
    return rows, variants, reps, seconds, source, builds


def common_time(base, candidate):
    runs = base + candidate
    if any(row['status'] != 'optimal' or not finite(row['objective']) for row in runs):
        return None
    if not all(close(row['objective'], base[0]['objective']) for row in runs):
        return None
    if not all(finite(row['elapsed_secs']) and row['elapsed_secs'] > 0 for row in runs):
        return None
    base_secs = median(row['elapsed_secs'] for row in base)
    candidate_secs = median(row['elapsed_secs'] for row in candidate)
    return {'base_median_secs': base_secs, 'candidate_median_secs': candidate_secs,
            'candidate_over_base': candidate_secs / base_secs}


def summarize(path, legacy=None, modern=None, public_truth=None):
    rows, variants, reps, seconds, source, builds = load_runs(path)
    inputs = [item for item in (legacy, modern, path) if item is not None]
    checked = inspect_results(inputs, public_truth=public_truth)
    errors = list(checked['errors'])
    grouped = defaultdict(list)
    proof_statuses = defaultdict(set)
    for row in rows:
        name, variant = row['instance'], row['variant']
        grouped[name, variant].append(row)
        if row['status'] in PROOFS:
            proof_statuses[name].add(row['status'])
        if row['auto_degraded'] or row['gap_terminated']:
            errors.append(f'{name}/{variant}/{row["rep"]} 关闭早停后仍记录提前退出')
    for name, statuses in proof_statuses.items():
        if len(statuses) > 1:
            errors.append(f'{name} 的证明终态矛盾：{sorted(statuses)}')

    counts, comparisons = {}, {}
    for variant in variants:
        members = [row for name in CASES for row in grouped[name, variant]]
        counts[variant] = {
            'features': list(FEATURES[variant]), 'binary': builds[variant][0],
            'binary_sha256': builds[variant][1],
            'status_counts': dict(Counter(row['status'] for row in members)),
            'proof_runs': sum(row['status'] in PROOFS for row in members),
            'proved_models_in_all_reps': sum(all(row['status'] in PROOFS for row in grouped[name, variant])
                                             for name in CASES),
            'missing_rss_runs': sum(row['peak_rss_bytes'] is None for row in members),
        }
        if variant == 'base':
            continue
        times, excluded = {}, {}
        for name in CASES:
            base, candidate = grouped[name, 'base'], grouped[name, variant]
            timing = common_time(base, candidate)
            if timing is None:
                excluded[name] = {'base_statuses': [row['status'] for row in base],
                                  'candidate_statuses': [row['status'] for row in candidate]}
            else:
                times[name] = timing
        comparisons[variant] = {'common_optimal_models': len(times), 'times': times,
                                'excluded_from_timing': excluded}
    return {
        'correctness_passed': checked['passed'] and not errors, 'correctness_errors': errors,
        'independent_check': checked, 'models': len(CASES), 'repetitions': reps, 'time_limit_secs': seconds,
        'source_commit': source[0], 'source_tree': source[1], 'binary_manifest_sha256': source[2],
        'host': source[3], 'variants': counts, 'comparisons_to_base': comparisons,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--legacy-baseline', type=Path)
    parser.add_argument('--modern-baseline', type=Path)
    parser.add_argument('--public-truth', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f'已有汇总文件 {args.output}')
    try:
        report = summarize(args.results, args.legacy_baseline, args.modern_baseline, args.public_truth)
    except (OSError, KeyError, TypeError, ValueError) as error:
        report = {'correctness_passed': False, 'input_error': str(error)}
        print(f'根阶段汇总输入无效：{error}', file=sys.stderr)
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(report, target, ensure_ascii=False, allow_nan=False, indent=2)
        target.write('\n')
    print(json.dumps({key: report[key] for key in ('correctness_passed', 'variants') if key in report},
                     ensure_ascii=False, indent=2))
    return not report['correctness_passed']


if __name__ == '__main__':
    sys.exit(main())
