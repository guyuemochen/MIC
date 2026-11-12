#!/usr/bin/env python3
"""按固定的十二例矩阵测量二元静态行传播。"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from compare import invoke, require, result, sha
from run_probe import AUTO_CASES, BASE_OPTIONS, METADATA
from run_row_products import source_info


PHASES = {'screen': 1, 'paired': 3}
SWITCH = '--binary-propagation'


def canonical_hash(model):
    payload = (json.dumps(model, ensure_ascii=False, allow_nan=False,
                          sort_keys=True, separators=(',', ':')) + '\n').encode('utf-8')
    return hashlib.sha256(payload).hexdigest()


def stage(baseline, output, mic):
    with baseline.open(encoding='utf-8') as stream:
        rows = [row for line in stream if (row := json.loads(line))['solver'] == 'mic' and row['rep'] == 0]
    by_name = {row['instance']: row for row in rows}
    require(len(rows) == len(by_name) and set(by_name) == set(AUTO_CASES),
            '外部对照 MIC 首轮与冻结的十二例名单不同')

    staged = {}
    for name in AUTO_CASES:
        row = by_name[name]
        require((row['schema'], row['threads'], row['time_limit_secs'], row['node_limit'])
                == (1, 1, 10, 2_000_000), f'{name} 基线预算与冻结设置不一致')
        original = row['command']
        options = dict(zip(original[3::2], original[4::2], strict=True))
        require(original[1] == 'solve' and len(options) * 2 == len(original) - 3
                and all(options.get(key) == value for key, value in BASE_OPTIONS.items())
                and options.get('--branch-rule', 'auto') == 'auto' and SWITCH not in options,
                f'{name} 基线命令与冻结设置不一致')
        path = Path(original[2])
        require(sha(path) == row['milp_sha256'], f'{name} .milp 与基线哈希不同')
        model_path = baseline.parent / row['model']
        require(sha(model_path) == row['model_sha256'], f'{name} 原模型与基线哈希不同')
        archived = json.loads(model_path.read_text(encoding='utf-8'))
        exported = invoke([str(mic), 'export', str(path)], 180)
        if 'mps_columns' in archived:
            exported['mps_columns'] = archived['mps_columns']
        require(canonical_hash(exported) == row['model_sha256'],
                f'{name} 用新二进制导出的原模型与基线不同')
        metadata = {key: row[key] for key in METADATA}
        metadata['model'] = os.path.relpath(model_path, output.parent)
        staged[name] = metadata, [str(mic), *original[1:]], options
    return staged


def run(phase, baseline, mic, source, output):
    require(not output.exists(), f'已有结果文件 {output}')
    src = source_info(source)
    digest = sha(mic)
    host = platform.platform()
    staged = stage(baseline, output, mic)
    output.parent.mkdir(parents=True, exist_ok=True)
    failed = False
    with output.open('x', encoding='utf-8') as target:
        for rep in range(PHASES[phase]):
            for index, name in enumerate(AUTO_CASES):
                metadata, start, options = staged[name]
                order = ('off', 'on') if (rep + index) % 2 == 0 else ('on', 'off')
                for position, variant in enumerate(order):
                    command = [*start, SWITCH, variant]
                    try:
                        row = invoke(command, 130, rss=True)
                    except (ValueError, RuntimeError, OSError, subprocess.TimeoutExpired) as error:
                        row = result('error', type(error).__name__, None, None, None, None, -1, -1, None)
                        row.update(error=str(error), peak_rss_bytes=None)
                    failed |= row['status'] == 'error'
                    mic_options = [token for pair in options.items() if pair[0] not in
                                   ('--threads', '--node-limit', '--time-limit-secs') for token in pair]
                    row.update(metadata, phase=phase, case=name, rep=rep, variant=variant,
                               order=position, command=command, mic_options=[*mic_options, SWITCH, variant],
                               version=digest, host=host, source_commit=src['commit'], source_tree=src['tree'])
                    target.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
                    target.flush()
                    print(name, rep, variant, row['status'], row['elapsed_secs'], flush=True)
    return failed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase', choices=PHASES, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--mic', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    return run(args.phase, args.baseline.resolve(), args.mic.resolve(), args.source.resolve(),
               args.output.resolve())


if __name__ == '__main__':
    sys.exit(main())
