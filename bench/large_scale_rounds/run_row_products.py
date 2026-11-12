#!/usr/bin/env python3
"""按第四轮冻结矩阵串行测量稀疏 LP，把诊断单独记录。"""
import argparse
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys

from compare import require, sha


HERE = Path(__file__).resolve().parent
ORIGINALS = ('lp-sparse-4000x800', 'lp-wide-1000x2000', 'setcover-300x600')
HOLDOUTS = ('scan2-sparse-wide', 'scan2-sparse-tall', 'scan2-medium', 'scan2-dense-wide')
CASES = ORIGINALS + HOLDOUTS
REPS = 3
SOLVE = re.compile(r'^status=(\w+) objective=(\S+) iterations=(\d+) wall_secs=(\S+)$')
RSS = re.compile(r'(\d+)\s+maximum resident set size')
COUNTERS = (
    'row_product_calls', 'row_product_hits', 'column_entries', 'touched_entries',
    'estimated_row_work', 'hit_column_entries', 'hit_touched_entries',
)


def source_info(path):
    def git(*args):
        return subprocess.check_output(['git', '-C', str(path), *args], text=True).strip()

    require(not git('status', '--porcelain', '--untracked-files=all', '--', 'milp'),
            f'MILP 源码树含未提交文件：{path}')
    return {'path': str(path), 'commit': git('rev-parse', 'HEAD'), 'tree': git('rev-parse', 'HEAD:milp')}


def input_info(holdout_dir, holdout_manifest):
    inputs = {}
    original_dir = HERE.parent / 'instances'
    for name in ORIGINALS:
        path = original_dir / f'{name}.milp'
        inputs[name] = {'path': str(path), 'sha256': sha(path)}
    frozen = json.loads(holdout_manifest.read_text(encoding='utf-8'))
    require(tuple(row['name'] for row in frozen) == HOLDOUTS, '留出清单与冻结矩阵不同')
    for entry in frozen:
        name = entry['name']
        path = holdout_dir / f'{name}.milp'
        require(sha(path) == entry['sha256'], f'{name} 与冻结留出哈希不同')
        inputs[name] = dict(entry, path=str(path))
    return inputs


def parse_output(stdout, stderr, diagnostic):
    lines = stdout.splitlines()
    solve = SOLVE.fullmatch(lines[0]) if lines else None
    require(solve is not None, f'缺少稀疏 LP 结果行：{stdout[-1000:]}')
    status, token, iterations, elapsed = solve.groups()
    objective = float(token)
    elapsed = float(elapsed)
    require(math.isfinite(elapsed) and elapsed >= 0, '稀疏 LP 耗时非法')
    memory = RSS.search(stderr)
    require(memory is not None and int(memory[1]) > 0, '缺少正值峰值 RSS')
    row = {
        'status': status, 'objective': objective if math.isfinite(objective) else None,
        'objective_text': token, 'iterations': int(iterations), 'wall_secs': elapsed,
        'peak_rss_bytes': int(memory[1]),
    }
    if diagnostic:
        counters = [line for line in lines[1:] if line.startswith('row_product_calls=')]
        require(len(counters) == 1, '诊断二进制缺少唯一的行乘积统计')
        fields = dict(item.split('=', 1) for item in counters[0].split())
        require(set(fields) == set(COUNTERS), '行乘积统计字段与诊断探针不一致')
        row['counters'] = {key: int(fields[key]) for key in COUNTERS}
        require(all(value >= 0 for value in row['counters'].values()), '行乘积统计含负值')
        row['extra_diagnostics'] = [line for line in lines[1:] if line != counters[0]]
    else:
        require(len(lines) == 1, '正式计时二进制启用了额外诊断输出')
    return row


def invoke(command, diagnostic=False):
    env = dict(os.environ, LC_ALL='C')
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                          env=env, start_new_session=True) as proc:
        try:
            stdout, stderr = proc.communicate(timeout=180)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.communicate()
            raise
    if proc.returncode:
        raise RuntimeError(f'进程退出码 {proc.returncode}\n{stderr[-3000:]}\n{stdout[-1000:]}')
    return parse_output(stdout, stderr, diagnostic)


def record(target, manifest, name, variant, rep, order):
    model = manifest['inputs'][name]
    binary = manifest['binaries'][variant]
    command = ['/usr/bin/time', '-l', binary['path'], model['path']]
    row = invoke(command, diagnostic=variant == 'diagnostic')
    row.update(case=name, variant=variant, rep=rep, order=order, command=command,
               input_sha256=model['sha256'], binary_sha256=binary['sha256'],
               source_commit=manifest['source']['commit'], source_tree=manifest['source']['tree'])
    target.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
    target.flush()
    print(name, variant, rep, row['status'], row['wall_secs'], row['peak_rss_bytes'], flush=True)


def run(output, binaries, source, holdout_dir, holdout_manifest):
    require(not output.exists(), f'结果目录已存在：{output}')
    manifest = {
        'source': source_info(source),
        'inputs': input_info(holdout_dir, holdout_manifest),
        'holdout_manifest_sha256': sha(holdout_manifest),
        'binaries': {name: {'path': str(path), 'sha256': sha(path)} for name, path in binaries.items()},
    }
    require(len({entry['sha256'] for entry in manifest['binaries'].values()}) == 3,
            '正式关闭、正式开启和诊断须为三个不同二进制')
    output.mkdir(parents=True)
    with (output / 'manifest.json').open('x', encoding='utf-8') as stream:
        json.dump(manifest, stream, ensure_ascii=False, allow_nan=False, indent=2)
        stream.write('\n')
    with (output / 'formal.jsonl').open('x', encoding='utf-8') as target:
        for rep in range(REPS):
            for index, name in enumerate(CASES):
                variants = ('off', 'on') if (rep + index) % 2 == 0 else ('on', 'off')
                for order, variant in enumerate(variants):
                    record(target, manifest, name, variant, rep, order)
    with (output / 'diagnostic.jsonl').open('x', encoding='utf-8') as target:
        for name in CASES:
            record(target, manifest, name, 'diagnostic', 0, 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', '--base', dest='baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--diagnostic', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--holdout-dir', type=Path, required=True)
    parser.add_argument('--holdout-manifest', type=Path, default=HERE / 'results' / 'scan-holdouts.json')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    run(args.output_dir.resolve(), {
        'off': args.baseline.resolve(),
        'on': args.candidate.resolve(),
        'diagnostic': args.diagnostic.resolve(),
    }, args.source.resolve(), args.holdout_dir.resolve(), args.holdout_manifest.resolve())


if __name__ == '__main__':
    sys.exit(main())
