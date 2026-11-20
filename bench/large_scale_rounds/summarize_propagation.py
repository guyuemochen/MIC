#!/usr/bin/env python3
"""核对第五轮证明和可行解，并按冻结的两条性能路径裁决。"""
import argparse
import json
import math
from pathlib import Path
from statistics import median
import sys

from check_results import close, finite, inspect_results
from compare import require
from run_probe import AUTO_CASES, BASE_OPTIONS
from run_propagation import PHASES, SWITCH


PROOF = {'optimal', 'infeasible', 'unbounded'}
COUNTERS = (
    'prop_accepted_rows', 'prop_skipped_rows', 'prop_calls', 'prop_row_visits',
    'prop_term_visits', 'prop_fixes', 'prop_conflicts',
)
MIN_SECS = 0.05
MIN_CASES = 3
MAX_CASE_RATIO = 1.25
RSS_FLOOR = 64 * 1024**2
PATH_A_RATIO = 0.90
PATH_B_RATIO = 1.05
PATH_B_MIN_GAIN = 2


def load_phase(path, phase):
    with path.open(encoding='utf-8') as stream:
        rows = [json.loads(line) for line in stream]
    expected = []
    for rep in range(PHASES[phase]):
        for index, name in enumerate(AUTO_CASES):
            order = ('off', 'on') if (rep + index) % 2 == 0 else ('on', 'off')
            expected.extend((name, rep, variant, position) for position, variant in enumerate(order))
    keys = ('instance', 'rep', 'variant', 'order')
    require([tuple(row[key] for key in keys) for row in rows] == expected,
            f'{path} 记录缺失、重复或与冻结顺序不同')
    for row in rows:
        command = row['command']
        options = dict(zip(command[3::2], command[4::2], strict=True))
        require(row['phase'] == phase and row['solver'] == 'mic' and command[1] == 'solve'
                and len(options) * 2 == len(command) - 3
                and (row['threads'], row['time_limit_secs'], row['node_limit']) == (1, 10, 2_000_000)
                and all(options.get(key) == value for key, value in BASE_OPTIONS.items())
                and options.get('--branch-rule', 'auto') == 'auto'
                and options.get(SWITCH) == row['variant'],
                f'{path} 记录的设置与冻结矩阵不同：{row["instance"]}/{row["rep"]}/{row["variant"]}')
    return rows


def brief(row, errors):
    details = {key: row[key] for key in ('status', 'objective', 'bound', 'nodes', 'lp_iterations',
                                        'elapsed_secs', 'peak_rss_bytes')}
    counters = {key: row.get(key) for key in COUNTERS}
    if (row['status'] not in PROOF | {'limit_reached'} or row['auto_degraded']
            or row['gap_terminated']
            or any(type(value) is not int or value < 0 for value in counters.values())):
        errors.append(f'{row["instance"]}/{row["phase"]}/{row["rep"]}/{row["variant"]} 状态、提前停止或传播统计异常')
    return dict(details, propagation=counters)


def difference(off, on, key):
    return on[key] - off[key] if finite(off[key]) and finite(on[key]) else None


def describe_pair(off_row, on_row, errors):
    off, on = brief(off_row, errors), brief(on_row, errors)
    a, b = off['status'], on['status']
    both_optimal = a == b == 'optimal'
    common_optimum = (both_optimal and finite(off['objective']) and finite(on['objective'])
                      and close(off['objective'], on['objective']))
    if (a in PROOF and b in PROOF and a != b) or (both_optimal and not common_optimum):
        errors.append(f'{off_row["instance"]}/{off_row["phase"]}/{off_row["rep"]} 开关两组的证明结果矛盾')
    base_rss, on_rss = off['peak_rss_bytes'], on['peak_rss_bytes']
    rss_limit = max(RSS_FLOOR, 2 * base_rss) if finite(base_rss) and base_rss > 0 else None
    return {
        'rep': off_row['rep'], 'off': off, 'on': on, 'common_optimum': common_optimum,
        'proof_lost': a in PROOF and (a != b or (both_optimal and not common_optimum)),
        'proof_gained': a == 'limit_reached' and b in PROOF,
        'optimal_proof_gained': a == 'limit_reached' and b == 'optimal',
        'objective_delta': difference(off, on, 'objective'), 'bound_delta': difference(off, on, 'bound'),
        'nodes_delta': difference(off, on, 'nodes'), 'lp_iterations_delta': difference(off, on, 'lp_iterations'),
        'rss_limit_bytes': rss_limit,
        'rss_passed': rss_limit is not None and finite(on_rss) and 0 < on_rss <= rss_limit,
    }


def describe_phase(rows, phase, errors):
    indexed = {(row['instance'], row['rep'], row['variant']): row for row in rows}
    cases = []
    for name in AUTO_CASES:
        pairs = [describe_pair(indexed[(name, rep, 'off')], indexed[(name, rep, 'on')], errors)
                 for rep in range(PHASES[phase])]
        common = all(pair['common_optimum'] for pair in pairs)
        off_median = on_median = ratio = None
        if common:
            off_times = [pair['off']['elapsed_secs'] for pair in pairs]
            on_times = [pair['on']['elapsed_secs'] for pair in pairs]
            if all(finite(t) and t > 0 for t in off_times + on_times):
                off_median, on_median = median(off_times), median(on_times)
                if phase == 'paired':
                    ratio = on_median / off_median
            else:
                errors.append(f'{name}/{phase} 共同证优记录缺少有效耗时')
        cases.append({
            'instance': name, 'pairs': pairs, 'all_common_optimum': common,
            'off_median_secs': off_median, 'on_median_secs': on_median, 'time_ratio': ratio,
            'counts_for_geometric_mean': ratio is not None and off_median >= MIN_SECS,
            'proof_lost': any(pair['proof_lost'] for pair in pairs),
            'new_proofs': sum(pair['proof_gained'] for pair in pairs),
            'new_optimal_proofs': sum(pair['optimal_proof_gained'] for pair in pairs),
            'new_optimum_in_all_three': phase == 'paired' and all(pair['optimal_proof_gained'] for pair in pairs),
            'rss_passed': all(pair['rss_passed'] for pair in pairs),
        })
    return cases


def promotion(screen, paired, correct):
    eligible = [case['time_ratio'] for case in paired if case['counts_for_geometric_mean']]
    geometric = math.exp(math.fsum(map(math.log, eligible)) / len(eligible)) if eligible else None
    new_models = [case['instance'] for case in paired if case['new_optimum_in_all_three']]
    path_a = geometric is not None and geometric <= PATH_A_RATIO
    path_b = geometric is not None and geometric <= PATH_B_RATIO and len(new_models) >= PATH_B_MIN_GAIN
    gates = {
        'correctness': correct,
        'no_proof_loss': not any(case['proof_lost'] for case in screen + paired),
        'enough_eligible_cases': len(eligible) >= MIN_CASES,
        'per_case_ratio': all(case['time_ratio'] is None or case['time_ratio'] <= MAX_CASE_RATIO
                              for case in paired),
        'peak_rss': all(case['rss_passed'] for case in screen + paired),
        'performance': path_a or path_b,
    }
    return {'passed': all(gates.values()), 'gates': gates, 'eligible_cases': len(eligible),
            'geometric_time_ratio': geometric, 'performance_paths': {'a': path_a, 'b': path_b},
            'newly_proven_models': new_models,
            'new_proofs': sum(case['new_proofs'] for case in paired),
            'new_optimal_proofs': sum(case['new_optimal_proofs'] for case in paired)}


def summarize(baseline, screen, paired=None):
    screen_rows = load_phase(screen, 'screen')
    paired_rows = load_phase(paired, 'paired') if paired else []
    checks = inspect_results([baseline, screen, *([paired] if paired else [])])
    errors = list(checks['errors'])
    builds = {(row['version'], row['source_commit'], row['source_tree'], row['host'])
              for row in screen_rows + paired_rows}
    if len(builds) != 1:
        errors.append('第五轮记录使用了不同的二进制、源码或主机')
    screen_cases = describe_phase(screen_rows, 'screen', errors)
    paired_cases = describe_phase(paired_rows, 'paired', errors) if paired else []
    correct = checks['passed'] and not errors
    return {
        'correctness_passed': correct, 'eligible_for_paired': correct, 'correctness_errors': errors,
        'independent_check': checks, 'screen': screen_cases, 'paired': paired_cases,
        'promotion': promotion(screen_cases, paired_cases, correct) if paired else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--screen', type=Path, required=True)
    parser.add_argument('--paired', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), f'已有汇总文件 {args.output}')
    report = summarize(args.baseline, args.screen, args.paired)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, ensure_ascii=False, allow_nan=False, indent=2)
        stream.write('\n')
    print(json.dumps({key: report[key] for key in ('correctness_passed', 'eligible_for_paired', 'promotion')},
                     ensure_ascii=False, indent=2))
    return not report['correctness_passed']


if __name__ == '__main__':
    sys.exit(main())
