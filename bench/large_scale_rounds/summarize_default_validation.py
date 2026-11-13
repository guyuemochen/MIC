#!/usr/bin/env python3
"""复核二十例默认回归，并逐模型对照第八轮正式配对的证明。"""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys

from check_results import close, finite, inspect_results
from compare import require, sha
from run_default_validation import SECONDS, VARIANT
from run_proof_diagnostics import STATUSES
from run_root_work import CASES, NODE_LIMIT, settings
from summarize_root_work import PROOFS


ARCHIVE = Path(__file__).resolve().parent / 'results'
R8_REFERENCE = ARCHIVE / 'binary_coefficient_strengthening-paired.jsonl'
PUBLIC_TRUTH = ARCHIVE / 'miplib2017-sources.json'
REFERENCE_REPS = 3
REFERENCE_STABLE_OPTIMAL = 10
SOURCE_FIELDS = ('source_commit', 'source_tree', 'binary_sha256', 'default_features', 'host')


def read_rows(path):
    with path.open(encoding='utf-8') as stream:
        return [json.loads(line) for line in stream]


def load_validation(path):
    rows = read_rows(path)
    require(tuple(row['instance'] for row in rows) == CASES, '回归记录与冻结的二十例顺序不同')
    source = {key: rows[0][key] for key in SOURCE_FIELDS}
    options = [token for pair in settings(SECONDS).items() for token in pair]
    for row in rows:
        name = row['instance']
        require(row['solver'] == 'mic' and row['variant'] == VARIANT and row['rep'] == 0
                and (row['threads'], row['time_limit_secs'], row['node_limit']) == (1, SECONDS, NODE_LIMIT)
                and row['command'][1] == 'solve' and row['command'][3:] == options
                and row['status'] in STATUSES | {'error'} and 'proof_profile' not in row,
                f'{name} 的试验设置或终态无效')
        require(all(row[key] == source[key] for key in SOURCE_FIELDS)
                and row['version'] == source['binary_sha256'] and row['command'][0] == rows[0]['command'][0],
                f'{name} 的源码或二进制来源与同批记录不同')
    return rows, source


def load_reference(path):
    candidates = defaultdict(list)
    for row in read_rows(path):
        if row['variant'] == 'candidate':
            candidates[row['instance']].append(row)
    require(set(candidates) == set(CASES), '第八轮候选组模型集合不符')
    source = {key: candidates[CASES[0]][0][key] for key in ('source_commit', 'source_tree', 'binary_sha256')}
    for name, runs in candidates.items():
        require([row['rep'] for row in runs] == list(range(REFERENCE_REPS))
                and all(row['solver'] == 'mic' and row['feature'] == 'binary_coefficient_strengthening'
                        and row['phase'] == 'paired' and row['time_limit_secs'] == SECONDS
                        and row['repetitions'] == REFERENCE_REPS and row['status'] in STATUSES
                        and all(row[key] == value for key, value in source.items()) for row in runs),
                f'{name} 的第八轮候选组缺少完整的 30 秒配对记录')
    return candidates, source


def summarize(results, legacy, modern, public_truth=PUBLIC_TRUTH, reference=R8_REFERENCE):
    rows, source = load_validation(results)
    historical, historical_source = load_reference(reference)
    independent = inspect_results([legacy, modern, results], public_truth=public_truth)
    errors = list(independent['errors'])
    cases = {}
    stable, not_reproduced = [], []
    for row in rows:
        name, status = row['instance'], row['status']
        prior = historical[name]
        require(all(item['model_sha256'] == row['model_sha256'] for item in prior),
                f'{name} 的原模型哈希与第八轮参考不同')
        prior_statuses = [item['status'] for item in prior]
        stable_optimal = all(item == 'optimal' for item in prior_statuses)
        optimum = None
        reproduced = None
        if stable_optimal:
            stable.append(name)
            optimum = prior[0]['objective']
            require(finite(optimum) and all(finite(item['objective']) and close(item['objective'], optimum)
                                            for item in prior),
                    f'{name} 的第八轮参考最优值无效或不一致')
            reproduced = status == 'optimal' and finite(row['objective']) and close(row['objective'], optimum)
            if not reproduced:
                not_reproduced.append(name)
                if status == 'optimal':
                    errors.append(f'{name} 的已证最优值与第八轮已证最优值 {optimum} 不符')
                elif status in PROOFS:
                    errors.append(f'{name} 的终态 {status} 与第八轮已证最优矛盾')
        if row['auto_degraded'] or row['gap_terminated']:
            errors.append(f'{name} 在关闭提前停止后仍记录提前退出')
        if row.get('error_stage') != 'export' and row['export_model_sha256'] != row['model_sha256']:
            errors.append(f'{name} 未通过重新导出原模型核验')
        if status != 'error' and (type(row['peak_rss_bytes']) is not int or row['peak_rss_bytes'] <= 0):
            errors.append(f'{name} 缺少有效峰值 RSS')
        cases[name] = {
            'status': status, 'objective': row['objective'], 'bound': row['bound'],
            'proof_completed': status in PROOFS, 'r8_candidate_statuses': prior_statuses,
            'r8_stable_optimal': stable_optimal, 'r8_optimum': optimum,
            'r8_stable_proof_reproduced': reproduced,
        }
    require(len(stable) == REFERENCE_STABLE_OPTIMAL, '第八轮候选组稳定证优模型数与冻结参考不同')
    return {
        'source': source, 'results': str(results), 'results_sha256': sha(results),
        'time_limit_secs': SECONDS, 'repetitions': 1, 'cases': cases,
        'status_counts': dict(Counter(row['status'] for row in rows)),
        'proofs_completed': sum(case['proof_completed'] for case in cases.values()),
        'r8_reference': {
            'path': str(reference), 'sha256': sha(reference), 'source': historical_source,
            'time_limit_secs': SECONDS, 'repetitions': REFERENCE_REPS,
            'stable_optimal_models': stable, 'stable_proofs_not_reproduced': not_reproduced,
        },
        'proof_regression_passed': all(cases[name]['r8_stable_proof_reproduced'] for name in stable),
        'correctness_passed': independent['passed'] and not errors,
        'correctness_errors': errors, 'independent_check': independent,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--legacy-baseline', type=Path, required=True)
    parser.add_argument('--modern-baseline', type=Path, required=True)
    parser.add_argument('--public-truth', type=Path, default=PUBLIC_TRUTH)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        require(not args.output.exists(), f'已有汇总文件 {args.output}')
        report = summarize(args.results, args.legacy_baseline, args.modern_baseline, args.public_truth)
        with args.output.open('x', encoding='utf-8') as target:
            json.dump(report, target, ensure_ascii=False, allow_nan=False, indent=2)
            target.write('\n')
    except (KeyError, TypeError, ValueError, OSError) as error:
        print(f'默认回归汇总输入或输出无效：{error}', file=sys.stderr)
        return 2
    print(json.dumps({key: report[key] for key in ('status_counts', 'proofs_completed', 'correctness_passed',
                                                  'correctness_errors', 'proof_regression_passed')},
                     ensure_ascii=False, indent=2))
    return not (report['correctness_passed'] and report['proof_regression_passed'])


if __name__ == '__main__':
    sys.exit(main())
