"""执行并复核 R12 两例根节点诊断，记录每轮候选池、实际施加顺序和回退。"""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from check_results import inspect_results
from compare import invoke, preflight_rss, require, sha
from run_probe import METADATA
from run_proof_diagnostics import baseline_rows, execute
from run_propagation import canonical_hash
from run_root_work import settings
from run_row_products import source_info


CASES = ('p0201', 'setcover-300x600')
FEATURES = {'base': ['proof_profile'], 'candidate': ['root_cut_recheck', 'proof_profile']}
ORDER = (('p0201', 'base'), ('p0201', 'candidate'),
         ('setcover-300x600', 'candidate'), ('setcover-300x600', 'base'))
ROOT_SETTINGS = settings(30)
ROOT_SETTINGS['--node-limit'] = '1'
OPTIONS = [token for pair in ROOT_SETTINGS.items() for token in pair]
CHUNK = 10


def write_json(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, allow_nan=False, indent=2)
        stream.write('\n')


def build_inputs(manifest, source):
    entries = json.loads(manifest.read_text())
    require(set(entries) == set(FEATURES), '诊断二进制分组不符')
    builds = {}
    for variant, features in FEATURES.items():
        entry = entries[variant]
        require(entry['features'] == features, f'{variant} 诊断特性不符')
        binary = (manifest.parent / entry['binary']).resolve()
        builds[variant] = (str(binary), sha(binary))
    require(len({digest for _, digest in builds.values()}) == len(FEATURES), '诊断二进制重复')
    return source_info(source), builds


def stage_models(baseline, builds, output):
    models = {}
    for name, row in baseline_rows(baseline, CASES).items():
        model = baseline.parent / row['model']
        command = row['command']
        require(row['schema'] == 1 and command[1] == 'solve', f'{name} 缺少基线原模型记录')
        instance = Path(command[2])
        require(sha(instance) == row['milp_sha256'] and sha(model) == row['model_sha256'],
                f'{name} 原模型与基线哈希不同')
        archived = json.loads(model.read_text())
        for binary, _ in builds.values():
            exported = invoke([binary, 'export', str(instance)], 180)
            if 'mps_columns' in archived:
                exported['mps_columns'] = archived['mps_columns']
            require(canonical_hash(exported) == row['model_sha256'], f'{name} 诊断导出哈希不同')
        metadata = {key: row[key] for key in METADATA}
        metadata['model'] = os.path.relpath(model, output.parent)
        models[name] = (str(instance), metadata)
    return models


def run(manifest, source, baseline, output):
    require(not output.exists(), f'已有诊断记录 {output}')
    preflight_rss()
    src, builds = build_inputs(manifest, source)
    models = stage_models(baseline, builds, output)
    host, manifest_hash = platform.platform(), sha(manifest)
    failed = False
    with output.open('x', encoding='utf-8') as target:
        for order, (name, variant) in enumerate(ORDER):
            binary, digest = builds[variant]
            instance, metadata = models[name]
            command = [binary, 'solve', instance, *OPTIONS]
            row = execute(command)
            failed |= bool(row.get('profile_missing_reason')) or row['status'] == 'error'
            row.update(metadata, features=FEATURES[variant], variant=variant, order=order % 2,
                       command=command, rep=0, threads=1, time_limit_secs=30, node_limit=1,
                       source_commit=src['commit'], source_tree=src['tree'], binary_sha256=digest,
                       binary_manifest_sha256=manifest_hash, export_model_sha256=metadata['model_sha256'],
                       host=host, version=digest, purpose='根割候选及回退机制诊断')
            target.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
            target.flush()
            print(name, variant, row['status'], row['cuts'], flush=True)
    return failed


def verify_profile(row):
    require(not row.get('profile_missing_reason'), row.get('profile_missing_reason') or '诊断不可用')
    profile = row['proof_profile']
    require(profile['available'] and profile['requested_threads'] == 1, '单线程根诊断不可用')
    rounds, pool, attempted, reviewed = [], [], [], []
    start = applied = 0
    expected = None
    stopped = False
    outer_discarded = 0
    root_events = []
    for index, event in enumerate(profile['events']):
        kind = event['kind']
        if kind in ('cut_entry', 'cut_outer_rollback', 'root_after_cuts', 'root_ready'):
            root_events.append(dict(event, event_index=index))
        if kind == 'cut_outer_rollback_begin':
            require(not outer_discarded and event['discarded'] == applied,
                    '整轮回退累计割数有误')
            root_events.append(dict(event, event_index=index))
            outer_discarded = applied
            applied = 0
        if kind == 'cut_round':
            require(not stopped and expected is None and not attempted and start == len(pool),
                    '上一轮候选未完成即开始新一轮')
            require(event['round'] == len(rounds), '根割轮次不连续')
            require(0 <= event['planned'] <= event['sanitize_survived'] <= event['generated']
                    and event['planned'] + event['truncated'] == event['sanitize_survived'],
                    '根割卫生检查及截断计数不一致')
            pool = list(range(event['planned']))
            start = 0
            rounds.append(dict(event, event_index=index, chunks=[]))
        elif kind == 'cut_chunk_candidate':
            require(not stopped and expected is None and event['round'] == len(rounds) - 1,
                    '根割添加事件发生在待复查或回退阶段')
            position = start + len(attempted)
            require(position < len(pool) and event['pool_position'] == position
                    and event['ordinal'] == pool[position], '实际添加顺序与候选池不一致')
            attempted.append(event['ordinal'])
        elif kind == 'cut_chunk':
            require(not stopped and expected is None and event['round'] == len(rounds) - 1
                    and event['chunk_start'] == start and event['applied_before'] == applied,
                    '实际小批位置或累计前缀不一致')
            count = min(CHUNK, len(pool) - start)
            require(event['attempted'] == len(attempted) == count and count > 0, '实际小批数量错误')
            status = event['status']
            require(event['committed_optimal'] is (status == 'optimal')
                    and event['rolled_back'] is (status == 'iteration_limit'), '小批回退标志与 LP 终态不符')
            rounds[-1]['chunks'].append(dict(event, event_index=index, ordinals=attempted))
            attempted = []
            start += count
            applied += count
            if status == 'iteration_limit':
                expected = 'rollback'
            elif status != 'optimal':
                stopped = True
            elif row['variant'] == 'candidate' and start < len(pool):
                expected = 'recheck'
        elif kind == 'cut_recheck_candidate':
            require(expected == 'recheck' and event['round'] == len(rounds) - 1,
                    '仅允许在候选组小批证优后复查')
            position = start + len(reviewed)
            require(position < len(pool) and position == event['pool_position']
                    and pool[position] == event['ordinal'] and type(event['retained']) is bool,
                    '复查顺序或保留标记不一致')
            reviewed.append(event)
        elif kind == 'cut_recheck':
            require(expected == 'recheck' and event['round'] == len(rounds) - 1,
                    '复查汇总没有相应的已证优小批')
            chunk = rounds[-1]['chunks'][-1]
            kept = [item['ordinal'] for item in reviewed if item['retained']]
            skipped = [item['ordinal'] for item in reviewed if not item['retained']]
            require(event['chunk_start'] == chunk['chunk_start']
                    and event['pending_before'] == len(reviewed) == len(pool) - start
                    and event['pending_after'] == len(kept) and event['skipped'] == len(skipped),
                    '复查计数与逐条决定不一致')
            chunk['recheck'] = dict(event, event_index=index, kept_ordinals=kept, skipped_ordinals=skipped)
            pool[start:] = kept
            reviewed = []
            expected = None
        elif kind == 'cut_rollback':
            require(expected == 'rollback' and event['round'] == len(rounds) - 1,
                    '小批回退没有相应的迭代失败')
            chunk = rounds[-1]['chunks'][-1]
            require(event['discarded'] == chunk['attempted']
                    and event['retained_cuts'] == chunk['applied_before'] == applied - event['discarded'],
                    '回退未正确保留先前提交的割前缀')
            applied -= event['discarded']
            chunk['rollback'] = dict(event, event_index=index)
            expected = None
            stopped = True
    require(expected is None and not attempted and not reviewed and (stopped or start == len(pool)),
            '根割诊断记录不完整')
    require(row['cuts'] == applied, '根结束时保留割数与逐批账不符')
    for event in root_events:
        if event['kind'] in ('cut_outer_rollback', 'root_after_cuts'):
            require(event['retained_cuts'] == applied, '根阶段快照与最终保留割数不符')
    fields = ('status', 'objective', 'bound', 'nodes', 'lp_iterations', 'cuts')
    return {**{key: row[key] for key in fields}, 'rounds': rounds, 'root_events': root_events,
            'outer_discarded': outer_discarded}


def summarize(results, baseline, output):
    rows = [json.loads(line) for line in results.read_text().splitlines()]
    require([(row['instance'], row['variant']) for row in rows] == list(ORDER), '诊断矩阵或顺序不符')
    source = {key: rows[0][key] for key in ('source_commit', 'source_tree', 'host', 'binary_manifest_sha256')}
    check = inspect_results([baseline, results])
    report = {'source': source, 'independent_check': check, 'cases': {}, 'errors': []}
    for index, row in enumerate(rows):
        name, variant = row['instance'], row['variant']
        try:
            require(all(row[key] == value for key, value in source.items())
                    and row['features'] == FEATURES[variant] and row['command'][1] == 'solve'
                    and row['command'][3:] == OPTIONS and row['order'] == index % 2
                    and (row['threads'], row['node_limit'], row['time_limit_secs']) == (1, 1, 30)
                    and row['export_model_sha256'] == row['model_sha256'], '冻结配置、来源或导出哈希不符')
            report['cases'].setdefault(name, {})[variant] = verify_profile(row)
        except (ValueError, KeyError, TypeError) as error:
            report['errors'].append(f'{name}/{variant} {error}')
    write_json(output, report)
    print(json.dumps({'independent_check_passed': check['passed'], 'errors': report['errors']}, ensure_ascii=False))
    return not check['passed'] or bool(report['errors'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    run_parser = commands.add_parser('run')
    for flag in ('binaries', 'source', 'baseline', 'output'):
        run_parser.add_argument(f'--{flag}', type=Path, required=True)
    summary_parser = commands.add_parser('summarize')
    for flag in ('results', 'baseline', 'output'):
        summary_parser.add_argument(f'--{flag}', type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.action == 'run':
            return run(args.binaries.resolve(), args.source.resolve(), args.baseline.resolve(), args.output.resolve())
        return summarize(args.results.resolve(), args.baseline.resolve(), args.output.resolve())
    except (KeyError, TypeError, ValueError, OSError, RuntimeError, subprocess.SubprocessError) as error:
        print(f'R12 根割诊断失败：{error}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
