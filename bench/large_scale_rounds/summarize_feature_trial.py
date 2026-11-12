#!/usr/bin/env python3
"""复核单特性试验的证明、有效界、耗时与内存。"""
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import sys

from check_results import close, finite, inspect_results
from compare import require
from run_feature_trial import FEATURES, VARIANTS, WINDOWS, expected_features
from run_root_work import CASES, NODE_LIMIT, schedule, settings
from summarize_root_work import PROOFS, common_time


MIN_BASE_SECS = 0.05
MAX_CASE_RATIO = 1.25
MIN_RSS_LIMIT_BYTES = 64 * 1024**2


def load_runs(path):
    with path.open(encoding='utf-8') as stream:
        rows = [json.loads(line) for line in stream]
    require(rows, f'{path} 没有求解记录')
    first = rows[0]
    feature, seconds, reps = first['feature'], first['time_limit_secs'], first['repetitions']
    require(feature in FEATURES and (seconds, reps) in WINDOWS, '特性或试验窗口无效')
    phase = WINDOWS[seconds, reps]
    keys = ('instance', 'rep', 'variant', 'order')
    require([tuple(row[key] for key in keys) for row in rows] == list(schedule(VARIANTS, reps)),
            f'{path} 记录缺失、重复或顺序与冻结矩阵不同')
    context = ('source_commit', 'source_tree', 'binary_manifest_sha256', 'host')
    source = {key: first[key] for key in context}
    builds = {}
    for row in rows:
        label = f'{row["instance"]}/{row["variant"]}'
        command = row['command']
        options = dict(zip(command[3::2], command[4::2], strict=True))
        require(row['solver'] == 'mic' and row['feature'] == feature and row['phase'] == phase
                and row['requested_variants'] == VARIANTS and row['repetitions'] == reps
                and (row['threads'], row['time_limit_secs'], row['node_limit']) == (1, seconds, NODE_LIMIT)
                and command[1] == 'solve' and len(options) * 2 == len(command) - 3 and options == settings(seconds)
                and row['features'] == expected_features(feature, row['variant'])
                and 'proof_profile' not in row, f'{label} 的设置与冻结矩阵不同')
        require(all(row[key] == source[key] for key in context), f'{label} 的运行来源或主机不同')
        binary = {'binary': command[0], 'binary_sha256': row['binary_sha256'], 'features': row['features']}
        require(row['version'] == row['binary_sha256'] and builds.setdefault(row['variant'], binary) == binary,
                f'{label} 的二进制与同组结果不同')
        if row.get('error_stage') != 'export':
            require(row['export_model_sha256'] == row['model_sha256'], f'{label} 未核对新导出模型')
    return rows, feature, phase, seconds, reps, source, builds


def paired_evidence(base, candidate):
    old, new = base['status'], candidate['status']
    proof_lost = old in PROOFS and (
        old != new or (old == 'optimal' and (
            not finite(base['objective']) or not finite(candidate['objective'])
            or not close(base['objective'], candidate['objective']))))
    base_rss = base['peak_rss_bytes']
    candidate_rss = candidate['peak_rss_bytes']
    limit = max(MIN_RSS_LIMIT_BYTES, 2 * base_rss) if finite(base_rss) and base_rss > 0 else None
    rss_ok = limit is not None and finite(candidate_rss) and 0 < candidate_rss <= limit
    fields = ('status', 'objective', 'bound', 'nodes', 'lp_iterations', 'peak_rss_bytes')
    off, on = (dict({key: row[key] for key in fields}, has_incumbent=row['values'] is not None)
               for row in (base, candidate))
    return {'rep': base['rep'], 'base': off, 'candidate': on,
            'proof_lost': proof_lost, 'rss_limit_bytes': limit, 'rss_passed': rss_ok}


def summarize(path, legacy=None, modern=None, public_truth=None):
    rows, feature, phase, seconds, reps, source, builds = load_runs(path)
    checked = inspect_results([p for p in (legacy, modern, path) if p is not None],
                              public_truth=public_truth)
    errors = list(checked['errors'])
    by_case = defaultdict(list)
    for row in rows:
        by_case[row['instance'], row['variant']].append(row)
        if row['auto_degraded'] or row['gap_terminated']:
            errors.append(f'{row["instance"]}/{row["variant"]}/{row["rep"]} 关闭早停后仍记录提前退出')

    cases = {}
    eligible = []
    slow = []
    lost = rss_failures = 0
    for name in CASES:
        base, candidate = (by_case[name, variant] for variant in VARIANTS)
        terminal_proofs = {row['status'] for row in base + candidate if row['status'] in PROOFS}
        if len(terminal_proofs) > 1:
            errors.append(f'{name} 的证明终态矛盾：{sorted(terminal_proofs)}')
        pairs = [paired_evidence(off, on) for off, on in zip(base, candidate, strict=True)]
        timing = common_time(base, candidate)
        timed_case_ok = None
        eligible_case = False
        if timing is not None:
            ratio = timing['candidate_over_base']
            timed_case_ok = ratio <= MAX_CASE_RATIO
            eligible_case = timing['base_median_secs'] >= MIN_BASE_SECS
            if eligible_case:
                eligible.append(ratio)
            if not timed_case_ok:
                slow.append(name)
        lost += sum(pair['proof_lost'] for pair in pairs)
        rss_failures += sum(not pair['rss_passed'] for pair in pairs)
        cases[name] = {'pairs': pairs, 'common_optimal_timing': timing, 'included_in_geomean': eligible_case,
                       'within_case_time_ratio': timed_case_ok}

    counts = {}
    for variant in VARIANTS:
        members = [row for name in CASES for row in by_case[name, variant]]
        counts[variant] = dict(builds[variant], status_counts=dict(Counter(row['status'] for row in members)),
                               proof_runs=sum(row['status'] in PROOFS for row in members),
                               proved_models_in_all_reps=sum(all(row['status'] in PROOFS
                                                                 for row in by_case[name, variant])
                                                             for name in CASES))
    geometric = math.exp(math.fsum(map(math.log, eligible)) / len(eligible)) if eligible else None
    return {
        'correctness_passed': checked['passed'] and not errors, 'correctness_errors': errors,
        'independent_check': checked, 'feature': feature, 'phase': phase, 'time_limit_secs': seconds,
        'repetitions': reps, 'source': source, 'variants': counts, 'cases': cases,
        'timing': {'common_optimal_models': sum(case['common_optimal_timing'] is not None for case in cases.values()),
                   'eligible_models': len(eligible), 'min_base_median_secs': MIN_BASE_SECS,
                   'geometric_candidate_over_base': geometric,
                   'max_case_ratio': MAX_CASE_RATIO, 'models_over_case_ratio': slow},
        'paired_proof_losses': lost, 'rss_failures': rss_failures,
        'rss_floor_bytes': MIN_RSS_LIMIT_BYTES, 'rss_base_multiplier': 2,
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
        print(f'单特性汇总输入无效：{error}', file=sys.stderr)
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(report, target, ensure_ascii=False, allow_nan=False, indent=2)
        target.write('\n')
    print(json.dumps({key: report[key] for key in
                      ('correctness_passed', 'variants', 'timing', 'paired_proof_losses', 'rss_failures')
                      if key in report}, ensure_ascii=False, indent=2))
    return not report['correctness_passed']


if __name__ == '__main__':
    sys.exit(main())
