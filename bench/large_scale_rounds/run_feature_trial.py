#!/usr/bin/env python3
"""在固定二十例上串行测量单个可选特性。"""
import argparse
import json
from pathlib import Path
import platform
import subprocess
import sys

from compare import RssPreflightError, positive_integer, preflight_rss, require, sha
from run_root_work import NODE_LIMIT, error_row, execute, schedule, settings, stage
from run_row_products import source_info


FEATURES = (
    'binary_coefficient_strengthening', 'risk_pivot_refresh', 'root_heuristic_work_limit',
)
VARIANTS = ['base', 'candidate']
WINDOWS = {(10, 1): 'screen', (30, 3): 'paired'}


def expected_features(feature, variant):
    return [feature] if variant == 'candidate' else []


def load_binaries(manifest, feature):
    entries = json.loads(manifest.read_text(encoding='utf-8'))
    require(isinstance(entries, dict) and set(entries) == set(VARIANTS),
            '二进制清单须包含 base 和 candidate')
    binaries = {}
    for variant in VARIANTS:
        entry = entries[variant]
        require(entry['features'] == expected_features(feature, variant),
                f'{variant} 的 Cargo 特性与试验组合不同')
        path = (manifest.parent / entry['binary']).resolve()
        binaries[variant] = (str(path), sha(path))
    require(binaries['base'][0] != binaries['candidate'][0], '两组须使用各自的二进制')
    return binaries


def run(feature, manifest, source, legacy, modern, seconds, reps, output):
    require(feature in FEATURES, f'未知试验特性：{feature}')
    require((seconds, reps) in WINDOWS, '筛查使用 10 秒一次，配对使用 30 秒三次')
    require(not output.exists(), f'已有结果文件 {output}')
    preflight_rss()
    source_data = source_info(source)
    manifest_hash = sha(manifest)
    binaries = load_binaries(manifest, feature)
    models, export_errors = stage(legacy, modern, binaries, output)
    host = platform.platform()
    options = [token for pair in settings(seconds).items() for token in pair]
    phase = WINDOWS[seconds, reps]
    output.parent.mkdir(parents=True, exist_ok=True)
    failed = False
    with output.open('x', encoding='utf-8') as target:
        for name, rep, variant, position in schedule(VARIANTS, reps):
            metadata, instance = models[name]
            binary, digest = binaries[variant]
            command = [binary, 'solve', instance, *options]
            failure = export_errors.get((name, variant))
            row = error_row('export', *failure) if failure else execute(command, seconds)
            failed |= row['status'] == 'error'
            row.update(metadata, case=name, feature=feature, phase=phase, rep=rep, variant=variant,
                       order=position, threads=1, time_limit_secs=seconds, node_limit=NODE_LIMIT,
                       command=command, mic_options=options, version=digest, host=host,
                       binary_sha256=digest, binary_manifest_sha256=manifest_hash,
                       features=expected_features(feature, variant), source_commit=source_data['commit'],
                       source_tree=source_data['tree'], requested_variants=VARIANTS, repetitions=reps,
                       export_model_sha256=None if failure else metadata['model_sha256'])
            target.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
            target.flush()
            print(name, rep, variant, row['status'], row['elapsed_secs'], flush=True)
    return failed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--feature', choices=FEATURES, required=True)
    parser.add_argument('--binaries', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--legacy-baseline', type=Path, required=True)
    parser.add_argument('--modern-baseline', type=Path, required=True)
    parser.add_argument('--seconds', type=positive_integer, default=10)
    parser.add_argument('--reps', type=positive_integer, default=1)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        return run(args.feature, args.binaries.resolve(), args.source.resolve(),
                   args.legacy_baseline.resolve(), args.modern_baseline.resolve(), args.seconds,
                   args.reps, args.output.resolve())
    except RssPreflightError as error:
        print(error, file=sys.stderr)
        return 1
    except (KeyError, TypeError, ValueError, OSError, subprocess.SubprocessError) as error:
        print(f'单特性对照输入无效：{error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
