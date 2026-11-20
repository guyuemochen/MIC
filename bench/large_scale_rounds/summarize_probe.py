#!/usr/bin/env python3
"""核对强分支复用结果并按冻结门槛裁决是否默认推广。"""
import argparse
import json
import math
from pathlib import Path
from statistics import median
import sys

from check_results import close, finite, inspect_results
from compare import require
from run_probe import PHASES, case_name


PROOF = {'optimal', 'infeasible', 'unbounded'}
MIN_SECS = 0.05
MIN_CASES = 3
MAX_GEOMETRIC_RATIO = 0.90
MAX_CASE_RATIO = 1.25
RSS_FLOOR = 64 * 1024**2


def load_pairs(path, phase):
    reps, cases = PHASES[phase]
    rows = {}
    with path.open(encoding='utf-8') as stream:
        for line in stream:
            row = json.loads(line)
            key = (row['case'], row['rep'], row['variant'])
            require(key not in rows and row['phase'] == phase and row['solver'] == 'mic',
                    f'{path} 包含重复记录或阶段错误 {key}')
            rows[key] = row
    result = []
    for instance, branch in cases:
        name = case_name(instance, branch)
        pairs = []
        for rep in range(reps):
            keys = [(name, rep, f'{branch}-{reuse}') for reuse in ('off', 'on')]
            require(all(key in rows for key in keys), f'{path} 缺少对照 {name}/{rep}')
            off, on = (rows.pop(key) for key in keys)
            require(off['instance'] == on['instance'] == instance, f'{name} 算例与冻结矩阵不同')
            pairs.append((off, on))
        result.append((name, pairs))
    require(not rows, f'{path} 含冻结矩阵外的记录 {list(rows)}')
    return result


def difference(off, on, key):
    a, b = off[key], on[key]
    return b - a if finite(a) and finite(b) else None


def pair_summary(name, off, on, errors):
    a, b = off['status'], on['status']
    off_proof, on_proof = a in PROOF, b in PROOF
    both_optimal = a == b == 'optimal'
    common_optimum = (both_optimal and finite(off['objective']) and finite(on['objective'])
                      and close(off['objective'], on['objective']))
    if (off_proof and on_proof and a != b) or (both_optimal and not common_optimum):
        errors.append(f'{name}/{off["rep"]} off/on 的证明结果互相矛盾')
    if off['version'] != on['version']:
        errors.append(f'{name}/{off["rep"]} off/on 使用了不同 MIC 二进制')
    if off['auto_degraded'] or on['auto_degraded']:
        errors.append(f'{name}/{off["rep"]} 关闭早停的记录出现自动提前停止')
    off_rss, on_rss = off['peak_rss_bytes'], on['peak_rss_bytes']
    rss_limit = max(RSS_FLOOR, 2 * off_rss) if finite(off_rss) and off_rss > 0 else None
    rss_passed = (rss_limit is not None and finite(on_rss) and 0 < on_rss <= rss_limit)
    return {
        'rep': off['rep'], 'binary_sha256': off['version'],
        'off_status': a, 'on_status': b, 'common_optimum': common_optimum,
        'proof_lost': off_proof and (a != b or (both_optimal and not common_optimum)),
        'proof_gained': on_proof and not off_proof,
        'off_objective': off['objective'], 'on_objective': on['objective'],
        'objective_delta': difference(off, on, 'objective'),
        'off_bound': off['bound'], 'on_bound': on['bound'], 'bound_delta': difference(off, on, 'bound'),
        'off_secs': off['elapsed_secs'], 'on_secs': on['elapsed_secs'],
        'off_peak_rss_bytes': off_rss, 'on_peak_rss_bytes': on_rss,
        'rss_limit_bytes': rss_limit, 'rss_passed': rss_passed,
    }


def phase_summary(path, phase, errors):
    cases = []
    for name, matched in load_pairs(path, phase):
        pairs = [pair_summary(name, off, on, errors) for off, on in matched]
        common = all(pair['common_optimum'] for pair in pairs)
        off_secs = on_secs = ratio = None
        if common:
            off_times = [pair['off_secs'] for pair in pairs]
            on_times = [pair['on_secs'] for pair in pairs]
            if all(finite(secs) and secs > 0 for secs in off_times + on_times):
                off_secs, on_secs = median(off_times), median(on_times)
                if phase == 'paired':
                    ratio = on_secs / off_secs
            else:
                errors.append(f'{name} 共同证优记录缺少有效时间')
        cases.append({
            'case': name, 'pairs': pairs, 'all_common_optimum': common,
            'off_median_secs': off_secs, 'on_median_secs': on_secs, 'time_ratio': ratio,
            'counts_for_geometric_mean': ratio is not None and off_secs >= MIN_SECS,
            'proof_lost': any(pair['proof_lost'] for pair in pairs),
            'rss_passed': all(pair['rss_passed'] for pair in pairs),
        })
    return cases


def summarize(screen, paired=None, baseline=None):
    inputs = [path for path in (screen, paired, baseline) if path is not None]
    checks = inspect_results(inputs)
    errors = list(checks['errors'])
    screen_cases = phase_summary(screen, 'screen', errors)
    paired_cases = phase_summary(paired, 'paired', errors) if paired else []
    versions = {pair['binary_sha256'] for case in screen_cases + paired_cases for pair in case['pairs']}
    if len(versions) != 1:
        errors.append('强分支筛查和重复测量记录使用了不同 MIC 二进制')
    result = {
        'correctness_passed': checks['passed'] and not errors, 'correctness_errors': errors,
        'independent_check': checks, 'screen': screen_cases, 'paired': paired_cases,
    }
    result['eligible_for_paired'] = result['correctness_passed']
    if paired:
        eligible = [case['time_ratio'] for case in paired_cases if case['counts_for_geometric_mean']]
        geo = math.exp(math.fsum(map(math.log, eligible)) / len(eligible)) if eligible else None
        gates = {
            'correctness': result['correctness_passed'],
            'no_proof_loss': not any(case['proof_lost'] for case in screen_cases + paired_cases),
            'enough_eligible_cases': len(eligible) >= MIN_CASES,
            'geometric_ratio': geo is not None and geo <= MAX_GEOMETRIC_RATIO,
            'per_case_ratio': all(case['time_ratio'] is None or case['time_ratio'] <= MAX_CASE_RATIO
                                  for case in paired_cases),
            'peak_rss': all(case['rss_passed'] for case in screen_cases + paired_cases),
        }
        result['promotion'] = {'passed': all(gates.values()), 'gates': gates,
                               'eligible_cases': len(eligible), 'geometric_time_ratio': geo}
    else:
        result['promotion'] = None
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--screen', type=Path, required=True)
    parser.add_argument('--paired', type=Path)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), f'已有汇总文件 {args.output}')
    result = summarize(args.screen, args.paired, args.baseline)
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, ensure_ascii=False, allow_nan=False, indent=2)
        target.write('\n')
    print(json.dumps({key: result[key] for key in ('correctness_passed', 'eligible_for_paired', 'promotion')},
                     ensure_ascii=False, indent=2))
    return not result['correctness_passed']


if __name__ == '__main__':
    sys.exit(main())
