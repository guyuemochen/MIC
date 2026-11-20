#!/usr/bin/env python3
"""检查第四轮配对结果并按冻结门槛计算是否保留候选。"""
import argparse
import json
import math
from pathlib import Path
from statistics import median
import sys

from compare import require
from run_row_products import CASES, COUNTERS, REPS


TARGET = 'lp-sparse-4000x800'
TARGET_LIMIT = 0.90
GEOMETRIC_LIMIT = 0.95
CASE_LIMIT = 1.10
RSS_RATIO = 1.10
RSS_BYTES = 8 * 1024**2


def load_runs(directory):
    manifest = json.loads((directory / 'manifest.json').read_text(encoding='utf-8'))
    with (directory / 'formal.jsonl').open(encoding='utf-8') as stream:
        formal = [json.loads(line) for line in stream]
    with (directory / 'diagnostic.jsonl').open(encoding='utf-8') as stream:
        diagnostic = [json.loads(line) for line in stream]
    expected = []
    for rep in range(REPS):
        for index, name in enumerate(CASES):
            variants = ('off', 'on') if (rep + index) % 2 == 0 else ('on', 'off')
            expected.extend((name, variant, rep, order) for order, variant in enumerate(variants))
    fields = ('case', 'variant', 'rep', 'order')
    require([tuple(row[key] for key in fields) for row in formal] == expected,
            '正式记录缺失、重复或顺序与冻结矩阵不同')
    require([tuple(row[key] for key in fields) for row in diagnostic] ==
            [(name, 'diagnostic', 0, 0) for name in CASES], '诊断记录缺失、重复或顺序错误')

    for row in formal + diagnostic:
        name, variant = row['case'], row['variant']
        model, binary = manifest['inputs'][name], manifest['binaries'][variant]
        require((row['input_sha256'], row['binary_sha256'], row['source_commit'], row['source_tree'])
                == (model['sha256'], binary['sha256'], manifest['source']['commit'],
                    manifest['source']['tree']), f'{name}/{variant} 的输入或来源哈希与清单不同')
        require(row['command'] == ['/usr/bin/time', '-l', binary['path'], model['path']],
                f'{name}/{variant} 的实际命令与清单不同')
        require(math.isfinite(row['wall_secs']) and row['wall_secs'] > 0 and row['peak_rss_bytes'] > 0,
                f'{name}/{variant} 缺少有效耗时或峰值 RSS')
    require(all('counters' not in row for row in formal), '正式记录混入诊断')
    for row in diagnostic:
        counts = row['counters']
        require(set(counts) == set(COUNTERS) and 0 <= counts['row_product_hits'] <= counts['row_product_calls'],
                f'{row["case"]} 的诊断计数错误')
    return manifest, formal, diagnostic


def signature(row):
    return row['status'], row['objective_text'], row['iterations']


def check_solution(name, row, reference, errors):
    run = f'{name}/{row["variant"]}/{row["rep"]}'
    objective = row['objective']
    if row['status'] != 'Optimal' or objective is None or not math.isfinite(objective):
        errors.append(f'{run} 未返回带有限目标的最优解，状态 {row["status"]}，目标 {row["objective_text"]}')
    if signature(row) != reference:
        errors.append(f'{run} 状态、目标或迭代数不同 {signature(row)}，基准 {reference}')


def rss_passed(off, on):
    return on <= RSS_RATIO * off or on - off <= RSS_BYTES


def summarize(directory):
    manifest, formal, diagnostic = load_runs(directory)
    errors = []
    cases = []
    by_diagnostic = {row['case']: row for row in diagnostic}
    for name in CASES:
        runs = [row for row in formal if row['case'] == name]
        ref = signature(runs[0])
        for row in runs + [by_diagnostic[name]]:
            check_solution(name, row, ref, errors)
        by_run = {(row['rep'], row['variant']): row for row in runs}
        off = [by_run[(rep, 'off')] for rep in range(REPS)]
        on = [by_run[(rep, 'on')] for rep in range(REPS)]
        off_secs = median(row['wall_secs'] for row in off)
        on_secs = median(row['wall_secs'] for row in on)
        rss_pairs = []
        for a, b in zip(off, on):
            base_rss, new_rss = a['peak_rss_bytes'], b['peak_rss_bytes']
            rss_pairs.append({'rep': a['rep'], 'off_bytes': base_rss, 'on_bytes': new_rss,
                              'ratio': new_rss / base_rss, 'delta_bytes': new_rss - base_rss,
                              'passed': rss_passed(base_rss, new_rss)})
        counts = by_diagnostic[name]['counters']
        calls, hits = counts['row_product_calls'], counts['row_product_hits']
        cases.append({
            'case': name, 'status': ref[0], 'objective_text': ref[1], 'iterations': ref[2],
            'off_median_secs': off_secs, 'on_median_secs': on_secs, 'time_ratio': on_secs / off_secs,
            'off_median_rss_bytes': median(row['peak_rss_bytes'] for row in off),
            'on_median_rss_bytes': median(row['peak_rss_bytes'] for row in on),
            'rss_pairs': rss_pairs, 'rss_passed': all(pair['passed'] for pair in rss_pairs),
            'diagnostic': {'counters': counts, 'hit_rate': hits / calls if calls else None},
        })
    by_name = {case['case']: case for case in cases}
    target_ratio = by_name[TARGET]['time_ratio']
    geometric = math.exp(math.fsum(math.log(case['time_ratio']) for case in cases) / len(cases))
    gates = {
        'correctness': not errors,
        'target_time_ratio': target_ratio <= TARGET_LIMIT,
        'geometric_time_ratio': geometric <= GEOMETRIC_LIMIT,
        'per_case_time_ratio': all(case['time_ratio'] <= CASE_LIMIT for case in cases),
        'peak_rss': all(case['rss_passed'] for case in cases),
    }
    return {
        'passed': all(gates.values()), 'gates': gates, 'correctness_errors': errors,
        'target_time_ratio': target_ratio, 'geometric_time_ratio': geometric,
        'thresholds': {'target_time_ratio': TARGET_LIMIT, 'geometric_time_ratio': GEOMETRIC_LIMIT,
                       'per_case_time_ratio': CASE_LIMIT, 'peak_rss_ratio': RSS_RATIO,
                       'peak_rss_delta_bytes': RSS_BYTES},
        'source': manifest['source'], 'binaries': manifest['binaries'], 'cases': cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path, help='run_row_products.py 生成的结果目录')
    args = parser.parse_args()
    output = args.directory / 'summary.json'
    require(not output.exists(), f'汇总文件已存在：{output}')
    report = summarize(args.directory)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, ensure_ascii=False, allow_nan=False, indent=2)
        stream.write('\n')
    print(json.dumps({key: report[key] for key in
                      ('passed', 'gates', 'target_time_ratio', 'geometric_time_ratio')},
                     ensure_ascii=False, indent=2))
    return not report['gates']['correctness']


if __name__ == '__main__':
    sys.exit(main())
