#!/usr/bin/env python3
"""执行或汇总第十轮两例 drayage 的根阶段诊断。"""
import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from compare import invoke, preflight_rss, require, sha
from run_probe import METADATA
from run_proof_diagnostics import baseline_rows, execute
from run_propagation import canonical_hash
from run_root_work import settings
from run_row_products import source_info
from summarize_proof_diagnostics import EVENT_FIELDS, snapshot


CASES = ('drayage-100-23', 'drayage-25-23')
FEATURES = {
    'old': ('root_heuristic_work_budget', 'proof_profile'),
    'new': ('root_heuristic_work_limit', 'proof_profile'),
}
BUDGET_FIELD = {'old': 'remaining_iterations', 'new': 'remaining_work_steps'}
ROOT_STAGES = ('root_lp', 'root_restart', 'cut_entry', 'cut_chunk', 'cut_rollback',
               'cut_outer_rollback', 'root_after_cuts', 'root_ready')


def schedule():
    for index, name in enumerate(CASES):
        order = tuple(FEATURES) if index == 0 else tuple(reversed(FEATURES))
        for position, variant in enumerate(order):
            yield name, variant, position


def options():
    opts = settings(10)
    opts['--node-limit'] = '1'
    return opts


def load_builds(path):
    entries = json.loads(path.read_text(encoding='utf-8'))
    require(set(entries) == set(FEATURES), '二进制清单须包含 old 和 new')
    builds = {}
    for variant, features in FEATURES.items():
        entry = entries[variant]
        require(isinstance(entry['features'], list) and len(entry['features']) == len(features)
                and set(entry['features']) == set(features), f'{variant} 的诊断特性与冻结配置不同')
        binary = (path.parent / entry['binary']).resolve()
        builds[variant] = (str(binary), sha(binary))
    require(builds['old'][0] != builds['new'][0], '诊断组须使用各自的二进制')
    return builds


def prepare(baseline, output, builds):
    models = {}
    for name, row in baseline_rows(baseline, CASES).items():
        command = row['command']
        require(row['schema'] == 1 and len(command) >= 3 and command[1] == 'solve',
                f'{name} 基线缺少原模型记录')
        milp = Path(command[2])
        model = baseline.parent / row['model']
        require(sha(milp) == row['milp_sha256'] and sha(model) == row['model_sha256'],
                f'{name} 原模型与基线哈希不同')
        archived = json.loads(model.read_text(encoding='utf-8'))
        for variant, (binary, _) in builds.items():
            exported = invoke([binary, 'export', str(milp)], 180)
            require(isinstance(exported, dict), f'{name}/{variant} 未返回模型 JSON')
            if 'mps_columns' in archived:
                exported['mps_columns'] = archived['mps_columns']
            require(canonical_hash(exported) == row['model_sha256'],
                    f'{name}/{variant} 导出的模型哈希与基线不同')
        metadata = {key: row[key] for key in METADATA}
        metadata['model'] = os.path.relpath(model, output.parent)
        models[name] = str(milp), metadata
    return models


def run(binaries, source, baseline, output):
    require(not output.exists(), f'已有诊断结果 {output}')
    preflight_rss()
    src = source_info(source)
    builds = load_builds(binaries)
    models = prepare(baseline, output, builds)
    host = platform.platform()
    args = [item for pair in options().items() for item in pair]
    output.parent.mkdir(parents=True, exist_ok=True)
    failed = False
    with output.open('x', encoding='utf-8') as stream:
        for name, variant, order in schedule():
            binary, digest = builds[variant]
            milp, metadata = models[name]
            command = [binary, 'solve', milp, *args]
            row = execute(command)
            failed |= row['status'] == 'error' or row.get('profile_missing_reason') is not None
            row.update(metadata, variant=variant, order=order, rep=0, threads=1, node_limit=1,
                       time_limit_secs=10, version=digest, binary_sha256=digest, host=host,
                       features=list(FEATURES[variant]), command=command, mic_options=args,
                       source_commit=src['commit'], source_tree=src['tree'],
                       export_model_sha256=metadata['model_sha256'])
            stream.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
            stream.flush()
            print(name, variant, row['status'], row['nodes'], row['lp_iterations'], flush=True)
    return failed


def observed(snapshots):
    return {'observed': bool(snapshots), 'snapshots': snapshots,
            'missing_reason': None if snapshots else '探针未记录该阶段'}


def describe(row, errors):
    groups = defaultdict(list)
    profile = row.get('proof_profile')
    label = f'{row["instance"]}/{row["variant"]}'
    if row['status'] == 'error' or row.get('profile_missing_reason'):
        errors.append(f'{label} {row.get("error") or row.get("profile_missing_reason")}')
    if not isinstance(profile, dict) or profile.get('available') is not True:
        errors.append(f'{label} 缺少可用的证明诊断')
        profile = {}
    budget_key = BUDGET_FIELD[row['variant']]
    other_units = set(BUDGET_FIELD.values()) - {budget_key}
    after_cuts = False
    for index, event in enumerate(profile.get('events', ())):
        kind = event['kind']
        if kind == 'root_after_cuts':
            after_cuts = True
        if kind == 'root_heuristic_budget':
            stage = 'post_cut_heuristic_budget' if after_cuts else 'pre_cut_heuristic_budget'
            keys = (budget_key, 'remaining_lp_calls')
            if other_units.intersection(event) or any(
                    type(event.get(key)) is not int or event[key] < 0 for key in keys):
                errors.append(f'{label} {stage} 的预算计数缺失或无效')
        elif kind in (*ROOT_STAGES, 'lp_failure'):
            stage, keys = kind, EVENT_FIELDS[kind]
        else:
            continue
        groups[stage].append(dict(snapshot(event, keys), event_index=index))
    fields = ('status', 'objective', 'bound', 'nodes', 'lp_iterations', 'cuts')
    stages = (*ROOT_STAGES, 'pre_cut_heuristic_budget', 'post_cut_heuristic_budget')
    stage_reports = {key: observed(groups[key]) for key in stages}
    restore = observed([event for event in groups['lp_failure']
                        if event.get('stage') == 'root_heuristic_restore'])
    if not restore['observed']:
        restore['missing_reason'] = '探针只记录割后根启发式恢复的失败事件'
    stage_reports['root_heuristic_restore_failure'] = restore
    return {
        **{key: row[key] for key in fields}, 'budget_unit_field': budget_key,
        'stages': stage_reports,
        'lp_failures': groups['lp_failure'], 'omitted_lp_failures': profile.get('failures_omitted'),
        'tree_lp_calls': profile.get('tree_lp_calls'), 'tree_lp_iterations': profile.get('tree_lp_iterations'),
    }


def summarize(results, output):
    require(not output.exists(), f'已有诊断汇总 {output}')
    with results.open(encoding='utf-8') as stream:
        rows = [json.loads(line) for line in stream]
    require([tuple(row[key] for key in ('instance', 'variant', 'order')) for row in rows]
            == list(schedule()), '诊断结果与冻结的四次运行不一致')
    source = {key: rows[0][key] for key in ('source_commit', 'source_tree', 'host')}
    builds = {}
    reports = defaultdict(dict)
    errors = []
    for row in rows:
        name, variant = row['instance'], row['variant']
        command = row['command']
        actual = dict(zip(command[3::2], command[4::2]))
        require(row['features'] == list(FEATURES[variant]) and command[1] == 'solve'
                and actual == options() and len(actual) * 2 == len(command) - 3
                and all(row[key] == value for key, value in source.items())
                and row['export_model_sha256'] == row['model_sha256'],
                f'{name}/{variant} 的诊断设置、源码或模型校验不一致')
        build = {'binary': command[0], 'binary_sha256': row['binary_sha256'], 'features': row['features']}
        require(builds.setdefault(variant, build) == build, f'{variant} 诊断二进制不一致')
        reports[name][variant] = describe(row, errors)
    report = {'source': source, 'builds': builds, 'cases': reports, 'errors': errors}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, ensure_ascii=False, allow_nan=False, indent=2)
        stream.write('\n')
    return bool(errors)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    phases = parser.add_subparsers(dest='phase', required=True)
    launch = phases.add_parser('run')
    for option in ('--binaries', '--source', '--modern-baseline', '--output'):
        launch.add_argument(option, type=Path, required=True)
    report = phases.add_parser('summarize')
    report.add_argument('results', type=Path)
    report.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.phase == 'run':
            return run(args.binaries.resolve(), args.source.resolve(), args.modern_baseline.resolve(),
                       args.output.resolve())
        return summarize(args.results, args.output)
    except (OSError, TypeError, KeyError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f'根启发式诊断失败：{error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
