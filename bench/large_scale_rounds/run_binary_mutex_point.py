#!/usr/bin/env python3
"""按冻结顺序执行二元互斥固定的两例定点对照。"""

import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from check_results import exceeds, inspect_results
from compare import RssPreflightError, invoke, preflight_rss, require, sha
from run_feature_trial import load_binaries
from run_probe import METADATA
from run_proof_diagnostics import baseline_rows
from run_propagation import canonical_hash
from run_root_work import NODE_LIMIT, execute, schedule, settings
from run_row_products import source_info


FEATURE = 'binary_mutex_row_fixing'
CASES = ('drayage-100-23', 'drayage-25-23')
VARIANTS = ('base', 'candidate')
SECONDS = 30


def evaluate(rows, correctness):
    failures = []
    improvements = {}
    for case in CASES:
        base, candidate = (rows[case, variant] for variant in VARIANTS)
        base_objective, candidate_objective = base['objective'], candidate['objective']
        base_bound, candidate_bound = base['bound'], candidate['bound']
        base_status, candidate_status = base['status'], candidate['status']
        if base_status != 'limit_reached' and candidate_status != base_status:
            failures.append(f'{case} 丢失基础组终态证明')
        if base_objective is not None and (candidate_objective is None
                or exceeds(candidate_objective, base_objective)):
            failures.append(f'{case} 候选可行目标退化')
        if base_bound is not None and (candidate_bound is None or exceeds(base_bound, candidate_bound)):
            failures.append(f'{case} 候选对偶界退化')
        base_rss, candidate_rss = base['peak_rss_bytes'], candidate['peak_rss_bytes']
        if base_rss is None or candidate_rss is None or candidate_rss > max(64 << 20, 2 * base_rss):
            failures.append(f'{case} 超过 RSS 门槛或证据缺失')
        if base_status == candidate_status == 'optimal' and candidate['elapsed_secs'] > 1.25 * base['elapsed_secs']:
            failures.append(f'{case} 两组证优但超过耗时门槛')
        improvements[case] = {
            'new_optimal': base_status == 'limit_reached' and candidate_status == 'optimal',
            'new_incumbent': base_objective is None and candidate_objective is not None,
            'stronger_bound': candidate_bound is not None and (
                base_bound is None or exceeds(candidate_bound, base_bound)),
        }
    proof_or_incumbent = any(item['new_optimal'] or item['new_incumbent'] for item in improvements.values())
    both_bounds = all(item['stronger_bound'] for item in improvements.values())
    return dict(passed=correctness['passed'] and not failures and (proof_or_incumbent or both_bounds),
                correctness=correctness, failures=failures, improvements=improvements)


def run(args):
    require(not args.output.exists() and not args.gate.exists(), '结果或判定文件已存在')
    preflight_rss()
    source = source_info(args.source)
    binaries = load_binaries(args.binaries, FEATURE)
    baseline = baseline_rows(args.reference, CASES)
    models = {}
    for case in CASES:
        row = baseline[case]
        instance = Path(row['command'][2])
        model_path = args.reference.parent / row['model']
        require(sha(instance) == row['milp_sha256'] and sha(model_path) == row['model_sha256'],
                f'{case} 原模型哈希不符')
        archived = json.loads(model_path.read_text())
        require(archived['sense'] == 'min', f'{case} 不符合冻结的最小化口径')
        for binary, _ in binaries.values():
            exported = invoke([binary, 'export', str(instance)], 180)
            if 'mps_columns' in archived:
                exported['mps_columns'] = archived['mps_columns']
            require(canonical_hash(exported) == row['model_sha256'], f'{case} 导出模型哈希不符')
        metadata = {key: row[key] for key in METADATA}
        metadata['model'] = os.path.relpath(model_path, args.output.parent)
        models[case] = metadata, str(instance)

    options = [token for pair in settings(SECONDS).items() for token in pair]
    results = {}
    manifest_digest = sha(args.binaries)
    with args.output.open('x', encoding='utf-8') as out:
        for case, rep, variant, order in schedule(list(VARIANTS), 1):
            if case not in CASES:
                continue
            metadata, instance = models[case]
            binary, digest = binaries[variant]
            command = [binary, 'solve', instance, *options]
            row = execute(command, SECONDS)
            row.update(metadata, case=case, rep=rep, variant=variant, order=order, threads=1,
                       node_limit=NODE_LIMIT, time_limit_secs=SECONDS, feature=FEATURE,
                       source_commit=source['commit'], source_tree=source['tree'],
                       binary_sha256=digest, binary_manifest_sha256=manifest_digest,
                       host=platform.platform(), version=digest, command=command, options=options,
                       features=[] if variant == 'base' else [FEATURE])
            out.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
            out.flush()
            print(case, variant, row['status'], row['elapsed_secs'], flush=True)
            results[case, variant] = row

    correctness = inspect_results([args.output, args.reference], public_truth=args.public_truth)
    gate = evaluate(results, correctness)
    with args.gate.open('x', encoding='utf-8') as out:
        json.dump(gate, out, ensure_ascii=False, indent=2)
        out.write('\n')
    return 0 if gate['passed'] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['source', 'binaries', 'reference', 'public-truth', 'output', 'gate']:
        parser.add_argument(f'--{name}', type=Path, required=True)
    args = parser.parse_args()
    try:
        return run(args)
    except (KeyError, ValueError, TypeError, OSError, RssPreflightError, subprocess.SubprocessError) as error:
        print(f'二元互斥定点试验失败：{error}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
