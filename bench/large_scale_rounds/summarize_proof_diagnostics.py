#!/usr/bin/env python3
"""复核诊断求解结果，展示真实根节点和搜索前沿的观测。"""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

from check_results import inspect_results
from compare import require
from run_proof_diagnostics import PROFILE_FIELDS, STATUSES, schedule, settings


ROOT_EVENTS = ('root_lp', 'root_after_cuts', 'root_ready', 'root_restart')
CUT_EVENTS = ('cut_entry', 'cut_round', 'cut_chunk', 'cut_rollback',
              'cut_outer_rollback_begin', 'cut_outer_rollback')
SAMPLE_EVENTS = ('frontier', 'node', 'chosen')
LP_FIELDS = ('status', 'bound', 'iterations', 'iterations_scope', 'seconds', 'retained_cuts')
FRONTIER_EXTENSIONS = ('unresolved_count', 'unresolved_bound', 'combined_bound')
EVENT_FIELDS = {
    'root_lp': LP_FIELDS, 'root_after_cuts': LP_FIELDS, 'root_restart': LP_FIELDS,
    'root_ready': ('inherited_lp_bound', 'incumbent', 'fractional_integer_count',
                   'root_rc_fixed', 'resolved_branch'),
    'context': ('work_variables', 'work_rows', 'effective_engine', 'resolved_branch',
                'root_cuts_enabled'),
    'cut_entry': LP_FIELDS, 'cut_outer_rollback': LP_FIELDS,
    'cut_outer_rollback_begin': ('discarded', 'status'),
    'cut_round': ('round', 'entry_bound', 'generated', 'sanitize_survived', 'truncated', 'planned'),
    'cut_chunk': ('round', 'chunk_start', 'attempted', 'applied_before', 'status', 'bound',
                  'iterations', 'committed_optimal', 'rolled_back'),
    'cut_rollback': ('round', 'discarded', 'retained_cuts', 'status', 'bound', 'iterations'),
    'frontier': ('nodes', 'exit', 'incumbent', 'stack_bound', 'size', 'unknown_est_count',
                 'weakest_ties', 'weakest_min_depth', 'depth1_pending', 'min_depth', 'max_depth',
                 *FRONTIER_EXTENSIONS),
    'node': ('depth', 'var', 'up', 'inherited_est', 'parent_lp_bound', 'best_bound_selection',
             'outcome', 'status', 'actual_bound', 'iterations'),
    'chosen': ('depth', 'var', 'observations', 'probed_this_node'),
    'lp_failure': ('stage', 'status', 'numerical_error', 'phase1_pending'),
}
PROBE_FIELDS = ('depth', 'calls', 'iterations', 'seconds', 'positive_gain', 'zero_gain',
                'negative_gain', 'nonfinite', 'infeasible', 'iteration_limit', 'time_limit',
                'unbounded')
SOLVE_FIELDS = ('status', 'objective', 'bound', 'gap_absolute', 'gap_percent', 'elapsed_secs',
                'nodes', 'lp_iterations', 'cuts', 'peak_rss_bytes')


def load_phase(path, phase):
    rows = []
    with path.open(encoding='utf-8') as stream:
        for number, line in enumerate(stream, 1):
            try:
                row = json.loads(line)
                require(isinstance(row, dict), '记录不是 JSON 对象')
                rows.append(row)
            except (TypeError, ValueError) as error:
                raise ValueError(f'{path}:{number} {error}') from error
    keys = ('instance', 'rep', 'variant', 'order')
    require([tuple(row[key] for key in keys) for row in rows] == list(schedule(phase)),
            f'{path} 记录缺失、重复或顺序与冻结矩阵不同')
    for row in rows:
        command = row['command']
        require(isinstance(command, list) and len(command) >= 3 and command[1] == 'solve',
                f'{path} 记录缺少实际求解命令')
        options = dict(zip(command[3::2], command[4::2], strict=True))
        budget = 1 if phase == 'roots' else 2_000_000
        require(row['phase'] == phase and row['solver'] == 'mic'
                and len(options) * 2 == len(command) - 3 and options == settings(phase, row['variant'])
                and (row['threads'], row['time_limit_secs'], row['node_limit']) == (1, 10, budget),
                f'{path} 的预算或开关与冻结矩阵不同：{row["instance"]}/{row["variant"]}')
    return rows


def snapshot(event, fields):
    output = event.copy()
    reasons = {}
    for field, value in event.items():
        if value is not None:
            continue
        status = event.get('status')
        if field == 'stack_bound' and event.get('unknown_est_count', 0):
            reasons[field] = '待处理栈有未知估计，未提供有限栈界'
        elif field == 'stack_bound' and event.get('size') == 0:
            reasons[field] = '待处理栈为空'
        elif field == 'unresolved_bound' and event.get('unresolved_count') == 0:
            reasons[field] = '没有未完成的树节点 LP'
        elif field == 'unresolved_bound':
            reasons[field] = '未完成的树节点 LP 没有提供有限继承界'
        elif (field == 'combined_bound' and event.get('size') == 0
              and event.get('unresolved_count') == 0 and event.get('incumbent') is None):
            reasons[field] = '待处理栈为空，且没有未完成节点或可行目标'
        elif 'bound' in field and status is not None and str(status).lower() != 'optimal':
            reasons[field] = f'该阶段 LP 状态为 {status}，无最优界'
        elif field == 'incumbent':
            reasons[field] = '快照没有记录可行目标'
        elif field == 'unresolved_est' and event.get('unresolved_lp_nodes') == 0:
            reasons[field] = '没有记录求解未完成的树节点 LP'
        elif field == 'numerical_error':
            reasons[field] = 'LP 引擎未提供数值错误文本'
        elif field == 'phase1_pending':
            reasons[field] = 'LP 引擎未提供第一阶段的待完成状态'
        elif field in ('processed_min_depth', 'processed_max_depth') and event.get('tree_lp_calls') == 0:
            reasons[field] = '没有执行树节点 LP'
        else:
            reasons[field] = '探针在此快照未提供有限观测'
    if reasons:
        output['null_reasons'] = reasons
    missing = {}
    unobserved = {}
    for field in fields:
        if field in event:
            continue
        outcome = event.get('outcome')
        if (event.get('kind') == 'node' and field in ('status', 'actual_bound', 'iterations')
                and outcome and outcome != 'lp'):
            unobserved[field] = f'本次节点处理结果为 {outcome}，未记录子 LP 样本'
        else:
            missing[field] = '探针事件未包含该字段'
    if missing:
        output['missing_fields'] = missing
    if unobserved:
        output['unobserved_fields'] = unobserved
    return output


def report_profile(row, error):
    profile = row.get('proof_profile')
    reason = row.get('profile_missing_reason')
    if not isinstance(profile, dict):
        error(reason or '缺少 proof_profile 对象')
        return {'available': False, 'complete': False, 'reason': reason or '求解记录中没有证明诊断'}
    if profile.get('available') is not True or profile.get('requested_threads') != 1:
        error(reason or '单线程证明诊断不可用')
        return {'available': False, 'complete': False, 'reason': reason or '求解记录未启用单线程证明诊断'}

    events = profile.get('events')
    probes = profile.get('probes')
    valid_events = isinstance(events, list)
    if not valid_events:
        error('证明诊断缺少 events 数组')
        events = []
    valid_probes = (isinstance(probes, list) and all(isinstance(bucket, dict) for bucket in probes)
                    and [bucket.get('depth') for bucket in probes] == list(range(9)))
    if not valid_probes:
        error('证明诊断缺少深度 0 至 8 的探测桶')
    valid_counters = all(type(profile.get(key)) is int and profile[key] >= 0
                         for key in ('chosen_omitted', 'failures_omitted'))
    if not valid_counters:
        error('证明诊断缺少有效的省略选择或 LP 失败计数')
    grouped = defaultdict(list)
    for index, event in enumerate(events):
        if not isinstance(event, dict) or not isinstance(event.get('kind'), str):
            error('证明诊断含无种类的事件')
            valid_events = False
            continue
        kind = event['kind']
        grouped[kind].append(dict(snapshot(event, EVENT_FIELDS.get(kind, ())), event_index=index))

    roots = {
        kind: {'snapshots': grouped[kind],
               'missing_reason': None if grouped[kind] else '探针未观测到该阶段'}
        for kind in ROOT_EVENTS
    }
    samples = {}
    for kind in SAMPLE_EVENTS:
        items = grouped[kind]
        samples[kind] = {
            'sample_count': len(items), 'snapshots': items,
            'missing_reason': None if items else '探针未观测到该类事件',
        }
    probe_snapshots = [snapshot(bucket, PROBE_FIELDS) for bucket in probes] if valid_probes else None
    metadata = snapshot({key: value for key, value in profile.items() if key not in ('events', 'probes')},
                        PROFILE_FIELDS)
    snapshots = [metadata, *(probe_snapshots or ()),
                 *(item for items in grouped.values() for item in items)]
    complete = (valid_events and valid_probes and valid_counters
                and all('missing_fields' not in item for item in snapshots))
    recognized = set(ROOT_EVENTS + CUT_EVENTS + SAMPLE_EVENTS + ('lp_failure',))
    return {
        'available': True, 'complete': complete, 'roots': roots, 'samples': samples,
        'cut_events': {kind: grouped[kind] for kind in CUT_EVENTS if grouped[kind]},
        'lp_failures': {'snapshots': grouped['lp_failure'], 'omitted': profile.get('failures_omitted')},
        'event_counts': {kind: len(items) for kind, items in grouped.items() if items},
        'other_events': {kind: items for kind, items in grouped.items() if kind not in recognized},
        'probes': probe_snapshots,
        'probes_missing_reason': None if valid_probes else '探针未提供完整的探测桶',
        'profile_metadata': metadata,
    }


def describe(row, errors):
    def error(message):
        errors.append(f'{row["phase"]}/{row["instance"]}/{row["variant"]} {message}')

    status = row['status']
    if status not in STATUSES:
        error(f'求解失败或终态无效：{row.get("error", status)}')
    if row['auto_degraded'] or row['gap_terminated']:
        error('关闭早停和 gap 目标后仍记录提前退出')
    fields = {key: row.get(key) for key in SOLVE_FIELDS}
    null_reasons = {'objective': '求解记录没有可行目标', 'bound': '求解记录没有有限终态对偶界',
                    'gap_absolute': '求解记录没有有限绝对差', 'gap_percent': '求解记录没有有限相对差'}
    reasons = {key: null_reasons.get(key, '探针或计时器未提供该值')
               for key, value in fields.items() if value is None}
    return dict(fields, variant=row['variant'], order=row['order'], null_reasons=reasons,
                profile=report_profile(row, error))


def summarize(legacy, modern, roots=None, tree=None, public_truth=None):
    inputs = [path for path in (roots, tree) if path is not None]
    require(inputs, '至少需要一个诊断阶段')
    phases = {phase: load_phase(path, phase) for phase, path in (('roots', roots), ('tree', tree)) if path}
    checked = inspect_results([legacy, modern, *inputs], public_truth=public_truth)
    errors = list(checked['errors'])
    builds = {(row['version'], row['source_commit'], row['source_tree'], row['host'])
              for rows in phases.values() for row in rows}
    if len(builds) != 1:
        errors.append('诊断记录使用了不同的二进制、源码或主机')
    cases = {}
    for phase, rows in phases.items():
        by_name = defaultdict(dict)
        for row in rows:
            by_name[row['instance']][row['variant']] = describe(row, errors)
        cases[phase] = dict(by_name)
    return {'correctness_passed': checked['passed'] and not errors,
            'diagnostics_complete': all(run['profile']['complete'] for phase in cases.values()
                                        for group in phase.values() for run in group.values()),
            'correctness_errors': errors, 'independent_check': checked,
            'diagnostic_builds': [dict(zip(('binary_sha256', 'source_commit', 'source_tree', 'host'), build,
                                           strict=True)) for build in sorted(builds)],
            'cases': cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--legacy-baseline', type=Path, required=True)
    parser.add_argument('--modern-baseline', type=Path, required=True)
    parser.add_argument('--roots', type=Path)
    parser.add_argument('--tree', type=Path)
    parser.add_argument('--public-truth', type=Path, help='冻结的公开模型来源与官方证明清单')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f'已有汇总文件 {args.output}')
    try:
        report = summarize(args.legacy_baseline, args.modern_baseline, args.roots, args.tree,
                           public_truth=args.public_truth)
    except (OSError, ValueError, KeyError, TypeError) as error:
        report = {'correctness_passed': False, 'input_error': str(error)}
        print(f'诊断汇总输入无效：{error}', file=sys.stderr)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, ensure_ascii=False, allow_nan=False, indent=2)
        stream.write('\n')
    print(json.dumps({'correctness_passed': report['correctness_passed'],
                      'diagnostics_complete': report.get('diagnostics_complete', False),
                      'independent_check': report.get('independent_check')}, ensure_ascii=False, indent=2))
    return not report['correctness_passed']


if __name__ == '__main__':
    sys.exit(main())
