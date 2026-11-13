#!/usr/bin/env python3
"""归档 p0201 两次根割诊断。旧源码有已知安全比值缺陷，只用于历史定位。"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from check_results import finite, inspect_results
from compare import invoke, preflight_rss, require, sha, store_model
from run_default_validation import default_features
from run_proof_diagnostics import execute
from run_propagation import canonical_hash
from run_root_work import settings
from run_row_products import source_info


RESULTS = Path(__file__).resolve().parent / 'results'
BASELINE = RESULTS / 'safe-primal-ratio-default.jsonl'
EXTERNAL = RESULTS / 'external-baseline.jsonl'
EVENTS = {'cut_round', 'cut_chunk', 'cut_rollback', 'cut_outer_rollback_begin',
          'cut_outer_rollback', 'root_after_cuts', 'root_ready'}


def matching_rows(path, solvers, variant=None):
    with path.open(encoding='utf-8') as stream:
        rows = [row for line in stream if (row := json.loads(line))['instance'] == 'p0201'
                and row['rep'] == 0 and row['solver'] in solvers
                and (variant is None or row['variant'] == variant)]
    require(len(rows) == len(solvers) and {row['solver'] for row in rows} == solvers,
            f'{path} 缺少唯一的 p0201 首次记录')
    return rows


def build_info(binary, source, features):
    info = source_info(source)
    return {
        'source_commit': info['commit'], 'source_tree': info['tree'], 'source_path': info['path'],
        'binary': str(binary), 'binary_sha256': sha(binary), 'default_features': default_features(source),
        'explicit_features': features,
        'expected_build_flags': ['--manifest-path', 'milp/Cargo.toml', '--release', '--example', 'round_probe',
                                 '--features', ','.join(features)],
    }


def dump(path, rows, *, lines=False):
    with path.open('x', encoding='utf-8') as stream:
        if lines:
            for row in rows:
                stream.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
        else:
            json.dump(rows, stream, indent=2, ensure_ascii=False, allow_nan=False)
            stream.write('\n')


def root_events(profile):
    result = []
    before = None
    for event in profile['events']:
        kind = event['kind']
        if kind not in EVENTS:
            continue
        saved = dict(event)
        if kind == 'cut_round':
            before = event['entry_bound']
        if kind == 'cut_chunk':
            bound = event['bound']
            saved['bound_before_chunk'] = before
            saved['committed_bound_gain'] = None
            if event['committed_optimal'] and finite(bound):
                if finite(before):
                    saved['committed_bound_gain'] = bound - before
                before = bound
        if kind == 'cut_rollback':
            before = event['bound']
        result.append(saved)
    return result


def run(old_binary, old_source, new_binary, new_source, baseline, output):
    require(not output.exists(), f'已有归档目录 {output}')
    preflight_rss()
    base, = matching_rows(baseline, {'mic'}, 'default')
    references = matching_rows(EXTERNAL, {'highs', 'scip'})
    options = [item for pair in settings(30).items() for item in pair]
    command = base['command']
    require(command[1] == 'solve' and command[3:] == options, '默认回归的 p0201 命令与冻结设置不同')
    old = build_info(old_binary, old_source, ['binary_coefficient_strengthening', 'proof_profile'])
    new = build_info(new_binary, new_source, ['proof_profile'])
    require(old['source_commit'].startswith('07f4d23'), '旧版必须来自第八轮冻结源码 07f4d23')
    require('binary_coefficient_strengthening' in new['default_features']
            and 'proof_profile' not in new['default_features'], '新版必须使用当前生产默认并显式启用诊断')
    require(new['source_tree'] == base['source_tree'] and new['default_features'] == base['default_features'],
            '新版 MILP 源树和默认特性须与安全比值回归记录一致')
    require(old['binary_sha256'] != new['binary_sha256'], '两组需要各自的诊断二进制')
    instance = Path(command[2])
    model_path = baseline.parent / base['model']
    model_digest = base['model_sha256']
    require(sha(instance) == base['milp_sha256'] and sha(model_path) == model_digest
            and base['export_model_sha256'] == model_digest, '默认回归原模型哈希不符')
    model = json.loads(model_path.read_text(encoding='utf-8'))
    require(model['sense'] == 'min', 'p0201 应为最小化模型')
    for ref in references:
        require(sha(EXTERNAL.parent / ref['model']) == model_digest
                and all(ref[key] == base[key] for key in ('model_sha256', 'milp_sha256', 'mps_sha256')),
                '外部参考原模型哈希不符')

    output.mkdir(parents=True)
    archived, digest = store_model(output, model)
    require(digest == model_digest, '归档模型哈希不符')
    for ref in references:
        ref['model'] = archived.name
    refs_path = output / 'p0201-external.jsonl'
    dump(refs_path, references, lines=True)

    report = {
        'instance': 'p0201', 'purpose': '根阶段差异定位，旧版存在已知安全比值缺陷，仅作历史参照',
        'performance_promotion_evidence': False, 'baseline': str(baseline), 'baseline_sha256': sha(baseline),
        'external_sha256': sha(EXTERNAL), 'runs': {},
    }
    for variant, binary, build in (('r8_coefficient', old_binary, old), ('safe_default', new_binary, new)):
        export_command = [str(binary), 'export', str(instance)]
        exported = invoke(export_command, 180)
        require(isinstance(exported, dict), '探针导出未返回 JSON 对象')
        if 'mps_columns' in model:
            exported['mps_columns'] = model['mps_columns']
        require(canonical_hash(exported) == model_digest, f'{variant} 新导出原模型与基线不同')
        solve = [str(binary), *command[1:]]
        solve[solve.index('--node-limit', 3) + 1] = '1'
        row = execute(solve)
        row.update({key: base[key] for key in ('schema', 'instance', 'group', 'solver', 'threads',
                                              'time_limit_secs', 'model_sha256', 'milp_sha256', 'mps_sha256')})
        row.update(build, variant=variant, rep=0, node_limit=1, model=archived.name, command=solve,
                   baseline_command=command, export_command=export_command, export_model_sha256=model_digest,
                   version=build['binary_sha256'], performance_promotion_evidence=False)
        path = output / f'{variant}.jsonl'
        dump(path, [row], lines=True)
        check = inspect_results([refs_path, path])
        events = root_events(row['proof_profile']) if row['proof_profile'] else []
        rss = row['peak_rss_bytes']
        profile_passed = (not row.get('profile_missing_reason') and type(rss) is int and rss > 0
                          and row['nodes'] == 1 and not row['auto_degraded'] and not row['gap_terminated']
                          and all(any(event['kind'] == kind for event in events)
                                  for kind in ('cut_round', 'root_after_cuts', 'root_ready')))
        report['runs'][variant] = {
            **build, 'profile_path': path.name, 'status': row['status'], 'retained_cuts': row['cuts'],
            'profile_passed': profile_passed, 'independent_check': check, 'root_events': events,
        }
    dump(output / 'summary.json', report)
    print(json.dumps({label: {key: run[key] for key in ('status', 'retained_cuts', 'profile_passed')}
                      for label, run in report['runs'].items()}, ensure_ascii=False, indent=2))
    return not all(run['profile_passed'] and run['independent_check']['passed'] for run in report['runs'].values())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('old-binary', 'old-source', 'new-binary', 'new-source', 'output-dir'):
        parser.add_argument('--' + flag, type=Path, required=True)
    parser.add_argument('--baseline', type=Path, default=BASELINE)
    args = parser.parse_args()
    try:
        return run(*(path.resolve() for path in (args.old_binary, args.old_source, args.new_binary,
                                                 args.new_source, args.baseline, args.output_dir)))
    except (KeyError, TypeError, ValueError, OSError, RuntimeError, subprocess.SubprocessError) as error:
        print(f'p0201 根阶段诊断未完成：{error}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
