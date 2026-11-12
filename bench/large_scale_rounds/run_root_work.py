#!/usr/bin/env python3
"""在固定二十例上串行比较根阶段工作量的 Cargo 特性组合。"""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from compare import invoke, positive_integer, require, result, sha
from run_probe import AUTO_CASES, METADATA
from run_proof_diagnostics import MODERN, STATUSES, baseline_rows
from run_propagation import canonical_hash
from run_row_products import source_info


FEATURES = {
    'base': (),
    'snapshot': ('root_snapshot_restore',),
    'budget': ('root_heuristic_work_budget',),
    'both': ('root_snapshot_restore', 'root_heuristic_work_budget'),
}
CASES = AUTO_CASES + MODERN
NODE_LIMIT = 2_000_000
RUN_ERRORS = (ValueError, RuntimeError, TypeError, OSError, subprocess.TimeoutExpired)


def validate_variants(variants):
    require(variants == list(FEATURES) or (
        len(variants) == 2 and variants[0] == 'base'
        and variants[1] in FEATURES and variants[1] != 'base'
    ), '使用固定四组，或显式指定 base 加一组候选')


def schedule(variants, reps):
    for rep in range(reps):
        for index, name in enumerate(CASES):
            shift = (index + rep) % len(variants)
            order = variants[shift:] + variants[:shift]
            for position, variant in enumerate(order):
                yield name, rep, variant, position


def settings(seconds):
    return {
        '--early-stop': 'off', '--probe-reuse': 'off', '--binary-propagation': 'auto',
        '--threads': '1', '--time-limit-secs': str(seconds), '--node-limit': str(NODE_LIMIT),
    }


def load_binaries(manifest, variants):
    entries = json.loads(manifest.read_text(encoding='utf-8'))
    require(isinstance(entries, dict) and not set(entries) - set(FEATURES)
            and all(variant in entries for variant in variants), '二进制清单缺少所需组合或含未知组合')
    binaries = {}
    for variant, entry in entries.items():
        features = entry['features']
        require(isinstance(features, list) and len(features) == len(FEATURES[variant])
                and set(features) == set(FEATURES[variant]), f'{variant} 的 Cargo 特性与冻结组合不同')
        if variant in variants:
            path = (manifest.parent / entry['binary']).resolve()
            binaries[variant] = (str(path), sha(path))
    require(len({path for path, _ in binaries.values()}) == len(variants), '不同组合须使用各自的二进制')
    return binaries


def stage(legacy, modern, binaries, output):
    models = {}
    export_errors = {}
    for baseline, names in ((legacy, AUTO_CASES), (modern, MODERN)):
        rows = baseline_rows(baseline, names)
        for name in names:
            row = rows[name]
            command = row['command']
            require(row['schema'] == 1 and isinstance(command, list) and len(command) >= 3
                    and command[1] == 'solve', f'{name} 基线缺少有效的模型记录')
            milp = Path(command[2])
            model = baseline.parent / row['model']
            require(sha(milp) == row['milp_sha256'] and sha(model) == row['model_sha256'],
                    f'{name} 模型文件与基线哈希不同')
            archived = json.loads(model.read_text(encoding='utf-8'))
            for variant, (binary, _) in binaries.items():
                try:
                    exported = invoke([binary, 'export', str(milp)], 180)
                    require(isinstance(exported, dict), '探针导出未返回 JSON 对象')
                    if 'mps_columns' in archived:
                        exported['mps_columns'] = archived['mps_columns']
                    require(canonical_hash(exported) == row['model_sha256'],
                            '新探针导出与基线原模型哈希不同')
                except RUN_ERRORS as error:
                    export_errors[name, variant] = (type(error).__name__, str(error))
            metadata = {key: row[key] for key in METADATA}
            metadata['model'] = os.path.relpath(model, output.parent)
            models[name] = metadata, str(milp)
    return models, export_errors


def error_row(stage, error_type, message):
    row = result('error', error_type, None, None, None, None, -1, -1, None)
    row.update(error=message, error_stage=stage, peak_rss_bytes=None)
    return row


def execute(command, seconds):
    try:
        row = invoke(command, seconds + 120, rss=True)
        require(isinstance(row, dict) and row.get('status') in STATUSES,
                '探针返回无效的 JSON 结果或终态')
        require('proof_profile' not in row, '正式计时二进制启用了 proof_profile')
        fields = ('objective', 'bound', 'values', 'nodes', 'lp_iterations', 'auto_degraded',
                  'gap_terminated', 'elapsed_secs', 'peak_rss_bytes')
        require(all(key in row for key in fields), '探针缺少统一结果字段')
        json.dumps(row, allow_nan=False)
        return row
    except RUN_ERRORS as error:
        return error_row('solve', type(error).__name__, str(error))


def run(manifest, source, legacy, modern, seconds, reps, variants, output):
    require(not output.exists(), f'已有结果文件 {output}')
    validate_variants(variants)
    src = source_info(source)
    manifest_hash = sha(manifest)
    binaries = load_binaries(manifest, variants)
    models, export_errors = stage(legacy, modern, binaries, output)
    host = platform.platform()
    options = [token for pair in settings(seconds).items() for token in pair]
    output.parent.mkdir(parents=True, exist_ok=True)
    failed = False
    with output.open('x', encoding='utf-8') as target:
        for name, rep, variant, position in schedule(variants, reps):
            metadata, instance = models[name]
            binary, digest = binaries[variant]
            command = [binary, 'solve', instance, *options]
            failure = export_errors.get((name, variant))
            row = error_row('export', *failure) if failure else execute(command, seconds)
            failed |= row['status'] == 'error'
            row.update(metadata, case=name, rep=rep, variant=variant, order=position,
                       threads=1, time_limit_secs=seconds, node_limit=NODE_LIMIT, command=command,
                       mic_options=options, version=digest, host=host, binary_sha256=digest,
                       binary_manifest_sha256=manifest_hash, features=list(FEATURES[variant]),
                       source_commit=src['commit'], source_tree=src['tree'],
                       requested_variants=variants, repetitions=reps,
                       export_model_sha256=None if failure else metadata['model_sha256'])
            target.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
            target.flush()
            print(name, rep, variant, row['status'], row['elapsed_secs'], flush=True)
    return failed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binaries', type=Path, required=True, help='二进制路径及 Cargo 特性清单')
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--legacy-baseline', type=Path, required=True)
    parser.add_argument('--modern-baseline', type=Path, required=True)
    parser.add_argument('--seconds', type=positive_integer, default=10)
    parser.add_argument('--reps', type=positive_integer, default=1)
    parser.add_argument('--variants', choices=FEATURES, nargs='+', default=list(FEATURES))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        return run(args.binaries.resolve(), args.source.resolve(), args.legacy_baseline.resolve(),
                   args.modern_baseline.resolve(), args.seconds, args.reps, args.variants,
                   args.output.resolve())
    except (KeyError, TypeError, ValueError, OSError, subprocess.SubprocessError) as error:
        print(f'根阶段对照输入无效：{error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
