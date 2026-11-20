#!/usr/bin/env python3
"""用同一个诊断探针串行执行冻结的根节点与树搜索矩阵。"""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from compare import invoke, require, result, sha
from run_probe import METADATA
from run_propagation import canonical_hash
from run_row_products import source_info


LEGACY = ('setcover-300x600', 'p0201', 'p0282', 'p0548', 'lseu')
MODERN = (
    'bnatt500', 'seymour', 'neos17', 'neos-4954672-berkel',
    'gen-ip054', 'icir97_tension', 'drayage-100-23', 'drayage-25-23',
)
ROOTS = {
    'root-cuts-off': ('dfs', 'mostfrac', '0', 'off'),
    'root-cuts-on': ('dfs', 'mostfrac', '0', 'on'),
}
TREE = {
    'A': ('dfs', 'mostfrac', '0', 'on'),
    'B': ('hybrid', 'mostfrac', '0', 'on'),
    'C': ('dfs', 'pseudo', '0', 'on'),
    'D': ('dfs', 'pseudo', '8', 'on'),
    'E': ('hybrid', 'pseudo', '8', 'on'),
}
STATUSES = {'optimal', 'infeasible', 'unbounded', 'limit_reached'}
PROFILE_FIELDS = ('requested_threads', 'chosen_omitted', 'failures_omitted',
                  'tree_lp_calls', 'tree_lp_iterations', 'unresolved_lp_nodes', 'unresolved_est',
                  'processed_min_depth', 'processed_max_depth')


def settings(phase, variant):
    config = ROOTS if phase == 'roots' else TREE
    node_rule, branch_rule, probes, root_cuts = config[variant]
    return {
        '--threads': '1', '--time-limit-secs': '10',
        '--node-limit': '1' if phase == 'roots' else '2000000',
        '--early-stop': 'off', '--binary-propagation': 'on', '--probe-reuse': 'off',
        '--rc-fixing': 'on', '--local-cuts': 'off', '--root-cuts': root_cuts,
        '--node-rule': node_rule, '--branch-rule': branch_rule, '--pc-probes': probes,
    }


def schedule(phase):
    names = LEGACY if phase == 'roots' else LEGACY + MODERN
    variants = list(ROOTS if phase == 'roots' else TREE)
    rep = 0
    for index, name in enumerate(names):
        shift = (index + rep) % len(variants)
        order = variants[shift:] + variants[:shift]
        for position, variant in enumerate(order):
            yield name, rep, variant, position


def baseline_rows(path, names):
    wanted = set(names)
    rows = {}
    with path.open(encoding='utf-8') as stream:
        for number, line in enumerate(stream, 1):
            try:
                row = json.loads(line)
                name = row['instance']
                if row['solver'] == 'mic' and row['rep'] == 0 and name in wanted:
                    require(name not in rows, f'首轮 MIC 记录重复：{name}')
                    rows[name] = row
            except (TypeError, KeyError, ValueError) as error:
                raise ValueError(f'{path}:{number} {error}') from error
    require(set(rows) == wanted, f'{path} 缺少 MIC 首轮记录：{sorted(wanted - set(rows))}')
    return rows


def stage(phase, legacy, modern, mic, output):
    groups = [(legacy, LEGACY)]
    if phase == 'tree':
        groups.append((modern, MODERN))
    staged = {}
    for baseline, names in groups:
        rows = baseline_rows(baseline, names)
        for name in names:
            row = rows[name]
            require(row['schema'] == 1, f'{name} 基线模型记录版本不支持')
            original = row['command']
            require(isinstance(original, list) and len(original) >= 3 and original[1] == 'solve',
                    f'{name} 基线缺少 MIC 求解命令')
            instance = Path(original[2])
            require(sha(instance) == row['milp_sha256'], f'{name} .milp 与基线哈希不同')
            model = baseline.parent / row['model']
            require(sha(model) == row['model_sha256'], f'{name} 已导出的模型与基线哈希不同')
            archived = json.loads(model.read_text(encoding='utf-8'))
            exported = invoke([str(mic), 'export', str(instance)], 180)
            require(isinstance(exported, dict), f'{name} 新探针模型导出格式无效')
            if 'mps_columns' in archived:
                exported['mps_columns'] = archived['mps_columns']
            require(canonical_hash(exported) == row['model_sha256'],
                    f'{name} 新探针导出的模型哈希与外部基线不同')
            metadata = {key: row[key] for key in METADATA}
            metadata['model'] = os.path.relpath(model, output.parent)
            staged[name] = metadata, str(instance)
    return staged


def execute(command):
    try:
        row = invoke(command, 130, rss=True)
        require(isinstance(row, dict), '求解探针没有返回 JSON 对象')
        require(row.get('status') in STATUSES, f'求解探针返回未知状态：{row.get("status")!r}')
        required = ('objective', 'bound', 'values', 'nodes', 'lp_iterations', 'cuts',
                    'auto_degraded', 'gap_terminated', 'elapsed_secs', 'peak_rss_bytes')
        require(all(key in row for key in required), '求解探针缺少统一结果字段')
        profile = row.get('proof_profile')
        if not isinstance(profile, dict):
            row['proof_profile'] = None
            row['profile_missing_reason'] = '探针未返回 proof_profile，需用 proof_profile 特性构建'
        elif profile.get('available') is not True or profile.get('requested_threads') != 1:
            row['profile_missing_reason'] = '证明诊断在本次单线程运行中不可用'
        else:
            missing = [key for key in (*PROFILE_FIELDS, 'events', 'probes') if key not in profile]
            if missing:
                row['profile_missing_reason'] = f'证明诊断字段不完整：{", ".join(missing)}'
        json.dumps(row, allow_nan=False)
        return row
    except (ValueError, RuntimeError, TypeError, OSError, subprocess.TimeoutExpired) as error:
        row = result('error', type(error).__name__, None, None, None, None, -1, -1, None)
        row.update(error=str(error), peak_rss_bytes=None, proof_profile=None,
                   profile_missing_reason='求解探针失败，缺少诊断快照')
        return row


def run(phase, legacy, modern, mic, source, output):
    require(not output.exists(), f'已有结果文件 {output}')
    src = source_info(source)
    digest = sha(mic)
    host = platform.platform()
    staged = stage(phase, legacy, modern, mic, output)
    output.parent.mkdir(parents=True, exist_ok=True)
    failed = False
    with output.open('x', encoding='utf-8') as target:
        for name, rep, variant, position in schedule(phase):
            metadata, instance = staged[name]
            options = settings(phase, variant)
            args = [token for pair in options.items() for token in pair]
            command = [str(mic), 'solve', instance, *args]
            row = execute(command)
            failed |= row['status'] == 'error' or row.get('profile_missing_reason') is not None
            row.update(metadata, phase=phase, case=name, rep=rep, variant=variant, order=position,
                       node_limit=int(options['--node-limit']), time_limit_secs=10, threads=1,
                       command=command, mic_options=args, version=digest, host=host,
                       source_commit=src['commit'], source_tree=src['tree'])
            target.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
            target.flush()
            print(name, variant, row['status'], row['nodes'], row['lp_iterations'], flush=True)
    return failed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase', choices=('roots', 'tree'), required=True)
    parser.add_argument('--legacy-baseline', type=Path, required=True)
    parser.add_argument('--modern-baseline', type=Path, required=True)
    parser.add_argument('--mic', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        return run(args.phase, args.legacy_baseline.resolve(), args.modern_baseline.resolve(),
                   args.mic.resolve(), args.source.resolve(), args.output.resolve())
    except (OSError, TypeError, KeyError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f'诊断输入或运行失败：{error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
