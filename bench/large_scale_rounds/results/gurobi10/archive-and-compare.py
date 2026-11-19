"""归档 Gurobi 10 运行原件，并与 R12 同一原模型的结果交叉核验。"""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import shutil
import sys

from check_results import inspect_results
from compare import require, sha
from compare_gurobi10 import SECONDS, SIZE_LIMIT, stage
from run_root_work import CASES, NODE_LIMIT


FIELDS = ('status', 'raw_status', 'objective', 'bound', 'elapsed_secs', 'peak_rss_bytes', 'nodes', 'lp_iterations')


def read_rows(path):
    with path.open(encoding='utf-8') as source:
        return [json.loads(line) for line in source]


def dump(path, data):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def relocate(rows, raw, output, refs):
    with output.open('x', encoding='utf-8') as target:
        for row in rows:
            original = (raw.parent / row['model']).resolve()
            archive = refs / 'models' / f'{row["model_sha256"]}.json'
            require(sha(original) == sha(archive) == row['model_sha256'], '重定位的原模型摘要不同')
            target.write(json.dumps(dict(row, model=os.path.relpath(archive, output.parent)),
                                    ensure_ascii=False, allow_nan=False) + '\n')


def compare(gurobi_rows, r12, legacy, modern):
    expected = stage(legacy, modern)
    mic_rows = read_rows(r12)
    require([row['instance'] for row in gurobi_rows] == list(CASES), 'Gurobi 20 例的顺序或完整性不同')
    require(len(mic_rows) == 120, 'R12 缺少完整的三遍配对')
    results = {}
    for index, row in enumerate(gurobi_rows):
        name = row['instance']
        _, model_ref = expected[name]
        for key in ('model_sha256', 'milp_sha256', 'mps_sha256'):
            require(row[key] == model_ref[key], f'{name} Gurobi 与冻结模型摘要不同')
        require(row['solver'] == row['variant'] == 'gurobi' and row['case_order'] == index
                and (row['threads'], row['time_limit_secs'], row['node_limit'], row['rep']) == (1, SECONDS, NODE_LIMIT, 0),
                f'{name} Gurobi 预算或分组不同')
        if row['status'] != 'error':
            params = row['parameters']
            require(row['version'] == '10.0.3' and row['model_roundtrip_passed'] is True
                    and type(row['peak_rss_bytes']) is int and row['peak_rss_bytes'] > 0
                    and (params['Threads'], params['TimeLimit'], params['NodeLimit'], params['MIPGap'],
                         params['MIPGapAbs'], params['DualReductions'], params['IntFeasTol'])
                        == (1, 30, 2_000_000, 0.0, 0.0, 1, 1e-6), f'{name} Gurobi 版本、参数或读回证据不同')
        matching = [trial for trial in mic_rows if trial['instance'] == name]
        require(all(all(trial[key] == row[key] for key in ('model_sha256', 'milp_sha256', 'mps_sha256',
                                                         'threads', 'time_limit_secs', 'node_limit', 'host'))
                    for trial in matching), f'{name} Gurobi 与 R12 的模型、试验主机或预算不同')
        require(Counter((trial['variant'], trial['rep']) for trial in matching)
                == Counter((variant, rep) for variant in ('base', 'candidate') for rep in range(3)),
                f'{name} R12 逐次记录不完整')
        by_variant = {}
        for variant in ('base', 'candidate'):
            by_variant[variant] = [{**{key: trial[key] for key in FIELDS},
                                    'cuts': trial['cuts'], 'rep': trial['rep']}
                                   for trial in matching if trial['variant'] == variant]
        results[name] = {'model_sha256': row['model_sha256'],
                         'gurobi': {key: row[key] for key in FIELDS},
                         'r12': by_variant}
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('run', 'r12', 'refs', 'output'):
        parser.add_argument(f'--{flag}', type=Path, required=True)
    args = parser.parse_args()
    work, r12, refs, output = (item.resolve() for item in (args.run, args.r12, args.refs, args.output))
    legacy = refs / 'external-baseline.jsonl'
    modern = refs / 'miplib2017-baseline.jsonl'
    truth = refs / 'miplib2017-sources.json'
    raw = work / 'gurobi10.jsonl'
    raw_rows = read_rows(raw)
    pairs = compare(raw_rows, r12, legacy, modern)
    output.mkdir()
    originals = output / 'original'
    originals.mkdir()
    sources = ('gurobi10.jsonl', 'run.log', 'run-exit-code.txt', 'frozen-permission-probe.json', 'freeze.json',
               'execution.sh')
    for name in sources:
        shutil.copyfile(work / name, originals / name)
    shutil.copyfile(Path(__file__), output / 'archive-and-compare.py')
    archive = output / raw.name
    relocate(raw_rows, raw, archive, refs)
    check = inspect_results([legacy, modern, r12, archive], truth)
    products = {'independent_check': check, 'status_counts': Counter(row['status'] for row in raw_rows),
                'license_size_limits': [row['instance'] for row in raw_rows if row.get('error_code') == SIZE_LIMIT],
                'r12_sha256': sha(r12), 'instances': pairs}
    dump(output / 'independent-check.json', check)
    dump(output / 'comparison.json', products)
    dump(output / 'archive-relocation.json', {
        'changed_row_field': 'model', 'originals_sha256': {name: sha(originals / name) for name in sources},
        'archived_record_sha256': sha(archive), 'states_match_20_models': len(pairs) == 20,
    })
    print(json.dumps({'check_passed': check['passed'], 'status_counts': products['status_counts'],
                      'license_size_limits': products['license_size_limits']}, ensure_ascii=False))
    return not check['passed']


if __name__ == '__main__':
    sys.exit(main())
