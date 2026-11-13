#!/usr/bin/env python3
"""根据单特性试验汇总判定能否进入正式配对或推广默认。"""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import sys

from check_results import close, finite
from compare import require
from run_feature_trial import FEATURES, VARIANTS, WINDOWS, expected_features
from run_proof_diagnostics import STATUSES
from run_root_work import CASES
from summarize_feature_trial import MAX_CASE_RATIO, MIN_BASE_SECS, MIN_RSS_LIMIT_BYTES
from summarize_root_work import PROOFS


RSS_BASE_MULTIPLIER = 2
MIN_ELIGIBLE = 3
MAX_SPEEDUP_RATIO = 0.90
MIN_NEW_PROOFS = 2
MAX_NEW_PROOF_RATIO = 1.05
SOURCE_FIELDS = ('source_commit', 'source_tree', 'binary_manifest_sha256', 'host')


def finite_float(value):
    number = float(value)
    require(math.isfinite(number), f'不允许非有限数 {value}')
    return number


def unique_object(items):
    result = {}
    for key, value in items:
        require(key not in result, f'JSON 字段重复 {key}')
        result[key] = value
    return result


def check_count(value, expected, label):
    require(type(value) is int and value == expected, f'{label} 与逐次记录不一致')


def check_pair(pair, rep, name):
    label = f'{name}/{rep}'
    check_count(pair['rep'], rep, f'{label} 的重复编号')
    base, candidate = (pair[variant] for variant in VARIANTS)
    for row in (base, candidate):
        require(row['status'] in STATUSES | {'error'}, f'{label} 终态无效')
        for key in ('objective', 'bound', 'peak_rss_bytes'):
            require(row[key] is None or finite(row[key]), f'{label} 的 {key} 无效')
        for key in ('nodes', 'lp_iterations'):
            require(row[key] is None or (type(row[key]) is int and row[key] >= 0),
                    f'{label} 的 {key} 无效')
        require(type(row['has_incumbent']) is bool, f'{label} 缺少有效的可行解标记')

    old, new = base['status'], candidate['status']
    proof_lost = old in PROOFS and (old != new or (old == 'optimal' and (
        base['objective'] is None or candidate['objective'] is None
        or not close(base['objective'], candidate['objective']))))
    base_rss, candidate_rss = base['peak_rss_bytes'], candidate['peak_rss_bytes']
    limit = (max(MIN_RSS_LIMIT_BYTES, RSS_BASE_MULTIPLIER * base_rss)
             if base_rss is not None and base_rss > 0 else None)
    rss_passed = limit is not None and candidate_rss is not None and 0 < candidate_rss <= limit
    require(pair['proof_lost'] is proof_lost and pair['rss_passed'] is rss_passed
            and (pair['rss_limit_bytes'] is None if limit is None else
                 finite(pair['rss_limit_bytes']) and pair['rss_limit_bytes'] == limit),
            f'{label} 的证明丢失或 RSS 判定与记录不一致')
    return proof_lost, rss_passed


def check_timing(case, name, correctness):
    runs = [pair[variant] for pair in case['pairs'] for variant in VARIANTS]
    objective = runs[0]['objective']
    common_optimal = objective is not None and all(
        row['status'] == 'optimal' and row['objective'] is not None
        and close(row['objective'], objective) for row in runs)
    timing = case['common_optimal_timing']
    eligible = False
    within_case_ratio = None
    ratio = None
    if timing is not None:
        base, candidate, reported = (timing[key] for key in
                                     ('base_median_secs', 'candidate_median_secs', 'candidate_over_base'))
        require(common_optimal and all(finite(value) and value > 0 for value in (base, candidate, reported)),
                f'{name} 的共同证优耗时无效')
        ratio = candidate / base
        require(ratio == reported, f'{name} 的耗时比与中位数不一致')
        eligible = base >= MIN_BASE_SECS
        within_case_ratio = ratio <= MAX_CASE_RATIO
    else:
        require(not (correctness and common_optimal), f'{name} 缺少共同证优耗时')
    require(case['included_in_geomean'] is eligible and case['within_case_time_ratio'] is within_case_ratio,
            f'{name} 的计时门槛标记与耗时不一致')
    return ratio, eligible


def check_correctness(summary, reps):
    independent = summary['independent_check']
    for key, minimum in (('runs', len(CASES) * len(VARIANTS) * reps), ('instances', len(CASES))):
        require(type(independent[key]) is int and independent[key] >= minimum,
                f'独立复核的 {key} 未覆盖完整矩阵')
    checked, errors = independent['passed'], independent['errors']
    all_errors = summary['correctness_errors']
    require(type(checked) is bool and isinstance(errors, list) and isinstance(all_errors, list)
            and all(isinstance(error, str) for error in errors + all_errors),
            '正确性复核字段无效')
    require(checked == (not errors) and all_errors[:len(errors)] == errors,
            '独立复核与汇总错误记录不一致')
    correctness = checked and not all_errors
    require(summary['correctness_passed'] is correctness, '正确性标记与复核结果不一致')
    return correctness


def decide(summary):
    feature, phase = summary['feature'], summary['phase']
    seconds, reps = summary['time_limit_secs'], summary['repetitions']
    require(feature in FEATURES and type(seconds) is int and type(reps) is int
            and (seconds, reps) in WINDOWS and WINDOWS[seconds, reps] == phase,
            '特性、阶段或试验窗口无效')
    source = summary['source']
    require(all(isinstance(source[key], str) and source[key].strip() for key in SOURCE_FIELDS),
            '运行来源字段无效')
    variants, cases = summary['variants'], summary['cases']
    require(isinstance(variants, dict) and set(variants) == set(VARIANTS), '试验分组与冻结矩阵不同')
    require(isinstance(cases, dict) and set(cases) == set(CASES), '算例集合与冻结矩阵不同')
    require(summary['rss_floor_bytes'] == MIN_RSS_LIMIT_BYTES
            and summary['rss_base_multiplier'] == RSS_BASE_MULTIPLIER, 'RSS 门槛与冻结规则不同')
    correctness = check_correctness(summary, reps)

    proof_losses = rss_failures = common_optimal = 0
    eligible_ratios, slow, new_proofs = [], [], []
    status_counts = {variant: Counter() for variant in VARIANTS}
    proved_models = Counter()
    for name in CASES:
        case = cases[name]
        pairs = case['pairs']
        require(isinstance(pairs, list) and len(pairs) == reps, f'{name} 缺少完整配对')
        for rep, pair in enumerate(pairs):
            lost, rss_passed = check_pair(pair, rep, name)
            proof_losses += lost
            rss_failures += not rss_passed

        proofs = {}
        terminal = set()
        for variant in VARIANTS:
            statuses = [pair[variant]['status'] for pair in pairs]
            status_counts[variant].update(statuses)
            terminal.update(set(statuses) & PROOFS)
            require(not (correctness and 'error' in statuses), f'{name} 的错误终态与正确性标记矛盾')
            proofs[variant] = sum(status in PROOFS for status in statuses)
            proved_models[variant] += proofs[variant] == reps
        require(not correctness or len(terminal) <= 1, f'{name} 的证明终态与正确性标记矛盾')
        if correctness:
            for pair in pairs:
                for variant in VARIANTS:
                    row = pair[variant]
                    if row['status'] == 'optimal':
                        require(row['has_incumbent'] and row['objective'] is not None and row['bound'] is not None
                                and close(row['objective'], row['bound']), f'{name} 的最优记录与正确性标记矛盾')
        if phase == 'paired' and proofs['base'] == 0 and proofs['candidate'] == reps:
            new_proofs.append(name)

        ratio, eligible = check_timing(case, name, correctness)
        if ratio is not None:
            common_optimal += 1
            if ratio > MAX_CASE_RATIO:
                slow.append(name)
            if eligible:
                eligible_ratios.append(ratio)

    check_count(summary['paired_proof_losses'], proof_losses, '证明丢失计数')
    check_count(summary['rss_failures'], rss_failures, 'RSS 失败计数')
    for variant in VARIANTS:
        group, counts = variants[variant], status_counts[variant]
        require(all(isinstance(group[key], str) and group[key].strip() for key in ('binary', 'binary_sha256'))
                and group['features'] == expected_features(feature, variant), f'{variant} 的构建记录无效')
        recorded = group['status_counts']
        require(isinstance(recorded, dict) and all(type(n) is int for n in recorded.values())
                and recorded == counts, f'{variant} 的终态计数与逐次记录不一致')
        check_count(group['proof_runs'], sum(counts[status] for status in PROOFS), f'{variant} 的证明次数')
        check_count(group['proved_models_in_all_reps'], proved_models[variant], f'{variant} 的稳定证明数')

    eligible_count = len(eligible_ratios)
    geo = math.exp(math.fsum(map(math.log, eligible_ratios)) / eligible_count) if eligible_count else None
    timing = summary['timing']
    check_count(timing['common_optimal_models'], common_optimal, '共同证优数')
    check_count(timing['eligible_models'], eligible_count, '几何平均有效例数')
    reported = timing['geometric_candidate_over_base']
    require(reported is None if geo is None else finite(reported) and math.isclose(reported, geo, rel_tol=1e-12),
            '几何平均与有效例耗时比不一致')
    require(timing['min_base_median_secs'] == MIN_BASE_SECS and timing['max_case_ratio'] == MAX_CASE_RATIO
            and timing['models_over_case_ratio'] == slow, '单例耗时汇总与冻结门槛或逐例数据不一致')

    thresholds = {'max_proof_losses': 0}
    gates = {'correctness': correctness, 'no_proof_loss': proof_losses == 0}
    if phase == 'paired' or feature == 'root_heuristic_work_limit':
        thresholds.update(rss_floor_bytes=MIN_RSS_LIMIT_BYTES, rss_base_multiplier=RSS_BASE_MULTIPLIER,
                          max_case_ratio=MAX_CASE_RATIO)
        gates.update(rss=rss_failures == 0, case_median=not slow)
    if phase == 'paired':
        thresholds.update(min_eligible_models=MIN_ELIGIBLE, min_base_median_secs=MIN_BASE_SECS,
                          max_geomean_ratio=MAX_SPEEDUP_RATIO, min_new_stable_proof_models=MIN_NEW_PROOFS,
                          max_geomean_ratio_with_new_proofs=MAX_NEW_PROOF_RATIO)
        gates.update(enough_eligible=eligible_count >= MIN_ELIGIBLE,
                     performance=geo is not None and (geo <= MAX_SPEEDUP_RATIO or (
                         len(new_proofs) >= MIN_NEW_PROOFS and geo <= MAX_NEW_PROOF_RATIO)))

    passed = all(gates.values())
    return {
        'feature': feature, 'phase': phase, 'time_limit_secs': seconds, 'repetitions': reps,
        'source': source, 'variants': variants, 'thresholds': thresholds, 'gates': gates,
        'overall_passed': passed, 'eligible_for_paired' if phase == 'screen' else 'eligible_for_default': passed,
        'correctness_errors': summary['correctness_errors'],
        'paired_proof_losses': proof_losses, 'rss_failures': rss_failures,
        'timing': dict(timing, geometric_candidate_over_base=geo),
        'new_stable_proof_models': new_proofs,
        'new_proof_evidence': {name: cases[name] for name in new_proofs},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary', type=Path, help='summarize_feature_trial.py 生成的 JSON')
    parser.add_argument('--output', type=Path, help='另存门禁 JSON，已有文件不覆盖')
    args = parser.parse_args()
    try:
        summary = json.loads(args.summary.read_text(encoding='utf-8'),
                             parse_float=finite_float, parse_constant=finite_float,
                             object_pairs_hook=unique_object)
        report = decide(summary)
        report['input_summary'] = str(args.summary)
        payload = json.dumps(report, ensure_ascii=False, allow_nan=False, indent=2)
        if args.output is not None:
            with args.output.open('x', encoding='utf-8') as target:
                target.write(payload + '\n')
    except (KeyError, TypeError, ValueError, OverflowError, OSError) as error:
        print(f'单特性门禁输入或输出无效：{error}', file=sys.stderr)
        return 2
    print(payload)
    return 0 if report['overall_passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
