#!/usr/bin/env python3
"""使用外部对照的 MIC 二进制、预算和原模型，串行测量强分支复用。"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from compare import invoke, require, result, sha


AUTO_CASES = (
    'knapsack-5000', 'mknap-5x300', 'setcover-300x600', 'lp-sparse-4000x800',
    'p0033', 'p0201', 'enigma', 'stein27', 'p0282', 'p0548', 'lseu', 'air03',
)
PHASES = {
    'screen': (1, tuple((name, 'auto') for name in AUTO_CASES) +
               tuple((name, 'pseudo') for name in ('p0033', 'p0201', 'lseu'))),
    'paired': (3, (('mknap-5x300', 'auto'), ('knapsack-5000', 'auto'),
                   ('p0033', 'pseudo'), ('lseu', 'pseudo'))),
}
METADATA = (
    'schema', 'instance', 'group', 'solver', 'threads', 'time_limit_secs', 'node_limit',
    'model_sha256', 'milp_sha256', 'mps_sha256', 'host', 'version',
)
BASE_OPTIONS = {
    '--early-stop': 'off', '--probe-reuse': 'off', '--threads': '1',
    '--node-limit': '2000000', '--time-limit-secs': '10',
}


def case_name(instance, branch):
    return f'{instance}-{branch}'


def stage(baseline, output, cases):
    with baseline.open(encoding='utf-8') as stream:
        rows = [row for line in stream if (row := json.loads(line))['solver'] == 'mic' and row['rep'] == 0]
    by_name = {row['instance']: row for row in rows}
    require(len(by_name) == len(rows), '外部对照 MIC 首轮有重复算例')

    staged = {}
    binaries = set()
    for name, _ in cases:
        if name in staged:
            continue
        require(name in by_name, f'外部对照缺少 {name}')
        row = by_name[name]
        require((row['schema'], row['threads'], row['time_limit_secs'], row['node_limit'])
                == (1, 1, 10, 2_000_000), f'{name} 基线预算与冻结设置不一致')
        command = row['command']
        options = dict(zip(command[3::2], command[4::2], strict=True))
        require(command[1] == 'solve' and len(options) * 2 == len(command) - 3
                and all(options.get(key) == value for key, value in BASE_OPTIONS.items())
                and options.get('--branch-rule', 'auto') == 'auto', f'{name} 基线命令与冻结设置不一致')
        model = baseline.parent / row['model']
        require(sha(model) == row['model_sha256'], f'{name} 原模型与哈希不一致')
        require(sha(Path(command[2])) == row['milp_sha256'], f'{name} .milp 与哈希不一致')
        binaries.add((Path(command[0]).resolve(), row['version']))
        metadata = {key: row[key] for key in METADATA}
        metadata['model'] = os.path.relpath(model, output.parent)
        staged[name] = metadata, command[:3], options
    require(len(binaries) == 1, '外部对照包含多个 MIC 二进制')
    binary, digest = next(iter(binaries))
    require(sha(binary) == digest, 'MIC 二进制与外部对照的哈希不一致')
    return staged


def run(phase, baseline, output):
    reps, cases = PHASES[phase]
    require(not output.exists(), f'已有结果文件 {output}')
    staged = stage(baseline, output, cases)
    output.parent.mkdir(parents=True, exist_ok=True)
    failed = False
    with output.open('x', encoding='utf-8') as target:
        for rep in range(reps):
            for index, (name, branch) in enumerate(cases):
                metadata, start, options = staged[name]
                order = ('off', 'on') if (rep + index) % 2 == 0 else ('on', 'off')
                for position, reuse in enumerate(order):
                    actual = options.copy()
                    actual['--probe-reuse'] = reuse
                    if branch == 'pseudo':
                        actual['--branch-rule'] = 'pseudo'
                    command = [*start, *(token for pair in actual.items() for token in pair)]
                    try:
                        row = invoke(command, 130, rss=True)
                    except (ValueError, RuntimeError, OSError, subprocess.TimeoutExpired) as error:
                        row = result('error', type(error).__name__, None, None, None, None, -1, -1, None)
                        row.update(error=str(error), peak_rss_bytes=None)
                    failed |= row['status'] == 'error'
                    row.update(metadata, phase=phase, case=case_name(name, branch), rep=rep,
                               variant=f'{branch}-{reuse}', order=position, command=command,
                               mic_options=[token for pair in actual.items() if pair[0] not in
                                            ('--threads', '--node-limit', '--time-limit-secs')
                                            for token in pair])
                    target.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
                    target.flush()
                    print(name, branch, rep, reuse, row['status'], row['elapsed_secs'], flush=True)
    return failed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase', choices=PHASES, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    return run(args.phase, args.baseline.resolve(), args.output.resolve())


if __name__ == '__main__':
    sys.exit(main())
