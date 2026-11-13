"""归档安全比值修复后的 lseu 单例诊断。"""
import argparse
from collections import Counter
import json
from pathlib import Path
import platform
import struct
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from compare import invoke, require, sha, store_model
from run_row_products import source_info


BASE_SOURCE = 'ba8f06f23dfd570b2d197a9fe4907e1a06ce1da4'
TIMING_FIELDS = {'elapsed_secs', 'peak_rss_bytes'}


def exact(left, right):
    if type(left) is not type(right):
        return False
    if isinstance(left, float):
        return struct.pack('!d', left) == struct.pack('!d', right)
    if isinstance(left, list):
        return len(left) == len(right) and all(exact(a, b) for a, b in zip(left, right))
    return left == right


def capture(binary, source, baseline, output):
    for name in ('lseu-profile.jsonl', 'lseu-parity.json'):
        require(not (output / name).exists(), f'已有诊断记录 {output / name}')
    with baseline.open(encoding='utf-8') as stream:
        rows = [row for line in stream if (row := json.loads(line))['instance'] == 'lseu']
    require(len(rows) == 1, '默认回归须包含唯一 lseu 记录')
    base = rows[0]
    require(base['source_commit'] == BASE_SOURCE and base['variant'] == 'default'
            and (base['threads'], base['time_limit_secs'], base['node_limit']) == (1, 30, 2_000_000),
            'lseu 的基线来源或窗口与冻结配置不同')
    archived = baseline.parent / base['model']
    original = Path(base['command'][2])
    require(sha(archived) == base['model_sha256'] and sha(original) == base['milp_sha256'],
            '基线模型哈希不符')
    diagnostic = source_info(source)
    export_command = [str(binary), 'export', str(original)]
    model = invoke(export_command, 180)
    original_model = json.loads(archived.read_text(encoding='utf-8'))
    if 'mps_columns' in original_model:
        model['mps_columns'] = original_model['mps_columns']
    model_path, model_digest = store_model(output, model)
    require(model_digest == base['model_sha256'], '诊断二进制重新导出模型的哈希与基线不同')

    command = [str(binary), *base['command'][1:]]
    row = invoke(command, 150, rss=True)
    profile = row['proof_profile']
    require(profile['available'] is True and profile['requested_threads'] == 1, '单线程证明轨迹缺失')
    fields = {key: exact(base[key], row[key]) for key in row.keys() - TIMING_FIELDS - {'proof_profile'}}
    binary_digest = sha(binary)
    row.update({key: base[key] for key in ('schema', 'instance', 'group', 'solver', 'threads',
                                           'time_limit_secs', 'node_limit', 'model_sha256',
                                           'milp_sha256', 'mps_sha256', 'rep', 'default_features')})
    row.update(host=platform.platform(), command=command, features=['proof_profile'], variant='profile',
               source_commit=diagnostic['commit'], source_tree=diagnostic['tree'], binary_sha256=binary_digest,
               version=binary_digest, model=model_path.name, export_model_sha256=model_digest)
    with (output / 'lseu-profile.jsonl').open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')

    failures = [dict(profile_event_index=index, **event) for index, event in enumerate(profile['events'])
                if event['kind'] == 'lp_failure']
    report = {
        'baseline': {
            'path': str(baseline), 'sha256': sha(baseline), 'command': base['command'],
            'source_commit': base['source_commit'], 'source_tree': base['source_tree'],
            'binary_sha256': base['binary_sha256'], 'elapsed_secs': base['elapsed_secs'],
        },
        'diagnostic': {
            **diagnostic, 'binary_sha256': binary_digest, 'source_diff_sha256': sha(output / 'source.diff'),
            'build_log_sha256': sha(output / 'build.log'), 'export_command': export_command,
            'export_model_sha256': model_digest, 'command': command, 'elapsed_secs': row['elapsed_secs'],
        },
        'fields_equal': dict(sorted(fields.items())),
        'passed': all(fields.values()),
        'changed': {key: {'baseline': base[key], 'diagnostic': row[key]} for key in sorted(fields) if not fields[key]},
        'failure_samples': len(failures),
        'failure_statuses': dict(Counter(event['status'] for event in failures)),
        'failure_events': failures,
        **{key: profile[key] for key in ('failures_omitted', 'chosen_omitted', 'tree_lp_calls',
                                        'tree_lp_iterations', 'unresolved_lp_nodes', 'unresolved_est')},
    }
    payload = json.dumps(report, ensure_ascii=False, allow_nan=False, indent=2)
    with (output / 'lseu-parity.json').open('x', encoding='utf-8') as stream:
        stream.write(payload + '\n')
    print(payload)
    return not report['passed']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    args = parser.parse_args()
    return capture(args.binary.resolve(), args.source.resolve(), args.baseline.resolve(), Path(__file__).resolve().parent)


if __name__ == '__main__':
    sys.exit(main())
