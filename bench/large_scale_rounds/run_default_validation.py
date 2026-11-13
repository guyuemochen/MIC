#!/usr/bin/env python3
"""在冻结的二十例上复核默认求解器，每例串行运行一次。"""
import argparse
import json
from pathlib import Path
import platform
import subprocess
import sys
import tomllib

from compare import RssPreflightError, preflight_rss, require, sha
from run_root_work import CASES, NODE_LIMIT, error_row, execute, settings, stage
from run_row_products import source_info


SECONDS = 30
VARIANT = 'default'


def default_features(source):
    with (source / 'milp/Cargo.toml').open('rb') as stream:
        features = tomllib.load(stream)['features']['default']
    require(isinstance(features, list) and all(isinstance(name, str) for name in features),
            'Cargo 默认特性声明无效')
    return features


def run(binary, source, legacy, modern, output):
    require(not output.exists(), f'已有结果文件 {output}')
    preflight_rss()
    src = source_info(source)
    features = default_features(source)
    digest = sha(binary)
    models, export_errors = stage(legacy, modern, {VARIANT: (str(binary), digest)}, output)
    host = platform.platform()
    options = [token for pair in settings(SECONDS).items() for token in pair]
    output.parent.mkdir(parents=True, exist_ok=True)
    failed = False
    with output.open('x', encoding='utf-8') as target:
        for name in CASES:
            metadata, instance = models[name]
            command = [str(binary), 'solve', instance, *options]
            export_error = export_errors.get((name, VARIANT))
            row = error_row('export', *export_error) if export_error else execute(command, SECONDS)
            if row['status'] != 'error':
                rss = row['peak_rss_bytes']
                if type(rss) is not int or rss <= 0:
                    row['rss_error'] = '未取得有效的峰值 RSS'
            failed |= row['status'] == 'error' or 'rss_error' in row
            row.update(metadata, case=name, rep=0, variant=VARIANT, threads=1,
                       time_limit_secs=SECONDS, node_limit=NODE_LIMIT, command=command,
                       version=digest, host=host, binary_sha256=digest, default_features=features,
                       source_commit=src['commit'], source_tree=src['tree'],
                       export_model_sha256=None if export_error else metadata['model_sha256'])
            target.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
            target.flush()
            print(name, row['status'], row['elapsed_secs'], flush=True)
    return failed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True, help='以 Cargo 默认特性构建的 release round_probe')
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--legacy-baseline', type=Path, required=True)
    parser.add_argument('--modern-baseline', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        return run(args.binary.resolve(), args.source.resolve(), args.legacy_baseline.resolve(),
                   args.modern_baseline.resolve(), args.output.resolve())
    except RssPreflightError as error:
        print(error, file=sys.stderr)
        return 1
    except (KeyError, TypeError, ValueError, OSError, subprocess.SubprocessError) as error:
        print(f'默认回归输入或输出无效：{error}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
