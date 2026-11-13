#!/usr/bin/env python3
"""按默认回归参数运行三个定点证明探针，并独立复核原模型。"""
import argparse
import json
from pathlib import Path
import platform
import subprocess
import sys

from check_results import inspect_results
from compare import invoke, preflight_rss, require, sha, store_model
from run_default_validation import SECONDS, default_features
from run_root_work import NODE_LIMIT, RUN_ERRORS, error_row, execute, settings
from run_propagation import canonical_hash
from run_row_products import source_info
from summarize_root_work import PROOFS


CASES = ('p0033', 'p0201', 'lseu')
EXTERNAL = Path(__file__).resolve().parent / 'results/external-baseline.jsonl'
METADATA = ('schema', 'instance', 'group', 'threads', 'time_limit_secs', 'node_limit',
            'model_sha256', 'milp_sha256', 'mps_sha256')


def selected_rows(path, accept, solvers):
    selected = {}
    with path.open(encoding='utf-8') as stream:
        for line in stream:
            row = json.loads(line)
            if row['instance'] not in CASES or not accept(row):
                continue
            key = row['instance'], row['solver']
            require(key not in selected, f'{path} 重复记录 {key}')
            selected[key] = row
    require(set(selected) == {(name, solver) for name in CASES for solver in solvers},
            f'{path} 缺少定点模型记录')
    return selected


def write_json(path, value, lines=False):
    with path.open('x', encoding='utf-8') as stream:
        if lines:
            for row in value:
                stream.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
        else:
            json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
            stream.write('\n')


def git(source, *args):
    return subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()


def stage_case(name, base, references, baseline, output):
    options = [item for pair in settings(SECONDS).items() for item in pair]
    require((base['threads'], base['time_limit_secs'], base['node_limit']) == (1, SECONDS, NODE_LIMIT)
            and base['command'][1] == 'solve' and base['command'][3:] == options,
            f'{name} 的基线命令与冻结窗口不同')
    instance = Path(base['command'][2])
    archived = baseline.parent / base['model']
    require(sha(instance) == base['milp_sha256'] and sha(archived) == base['model_sha256']
            and base['export_model_sha256'] == base['model_sha256'], f'{name} 的基线模型哈希不符')
    model = json.loads(archived.read_text(encoding='utf-8'))
    path, digest = store_model(output, model)
    require(digest == base['model_sha256'], f'{name} 的归档原模型哈希不符')
    checked_refs = []
    for ref in references:
        require(all(ref[key] == base[key] for key in ('schema', 'model_sha256', 'milp_sha256', 'mps_sha256'))
                and sha(EXTERNAL.parent / ref['model']) == digest, f'{name} 的外部参考模型哈希不符')
        checked_refs.append(dict(ref, model=path.name))
    return instance, model, path.name, checked_refs


def probe(binary, command, instance, original, expected_hash):
    export_command = [str(binary), 'export', str(instance)]
    try:
        exported = invoke(export_command, 180)
        require(isinstance(exported, dict), '探针导出未返回 JSON 对象')
        if 'mps_columns' in original:
            exported['mps_columns'] = original['mps_columns']
        digest = canonical_hash(exported)
        require(digest == expected_hash, '探针重新导出原模型的哈希与基线不同')
    except RUN_ERRORS as error:
        return error_row('export', type(error).__name__, str(error)), export_command, None
    return execute(command, SECONDS), export_command, digest


def run(binary, source, output, baseline):
    require(not output.exists(), f'已有归档目录 {output}')
    preflight_rss()
    info = source_info(source)
    candidate = {
        'commit': info['commit'], 'milp_tree': info['tree'], 'tree': git(source, 'rev-parse', 'HEAD^{tree}'),
        'binary_sha256': sha(binary), 'default_features': default_features(source),
    }
    originals = selected_rows(baseline, lambda row: row['rep'] == 0 and row['variant'] == 'default'
                              and row['solver'] == 'mic', ('mic',))
    external = selected_rows(EXTERNAL, lambda row: row['rep'] == 0 and row['solver'] in ('highs', 'scip'),
                             ('highs', 'scip'))
    first = originals[CASES[0], 'mic']
    baseline_commit = first['source_commit']
    require(all(all(row[key] == first[key] for key in ('source_commit', 'source_tree', 'binary_sha256',
                                                       'default_features')) for row in originals.values()),
            '基线记录来源不一致')
    require(git(source, 'rev-parse', f'{baseline_commit}:milp') == first['source_tree'],
            '基线提交的 MILP 源树哈希不符')
    base_source = {
        'commit': baseline_commit, 'milp_tree': first['source_tree'],
        'tree': git(source, 'rev-parse', f'{baseline_commit}^{{tree}}'),
        'binary_sha256': first['binary_sha256'], 'default_features': first['default_features'],
    }
    require(candidate['default_features'] == base_source['default_features'],
            '候选与基线的 Cargo 默认特性不同')
    changed = git(source, 'diff', '--name-status', baseline_commit, candidate['commit'], '--', 'milp').splitlines()
    difference = {'baseline': base_source, 'candidate': candidate, 'milp_changes': changed}
    host = platform.platform()

    output.mkdir(parents=True)
    cases = {}
    for name in CASES:
        base = originals[name, 'mic']
        refs = [external[name, solver] for solver in ('highs', 'scip')]
        instance, model, model_name, archived_refs = stage_case(name, base, refs, baseline, output)
        command = [str(binary), *base['command'][1:]]
        row, export_command, export_digest = probe(binary, command, instance, model, base['model_sha256'])
        row.update({key: base[key] for key in METADATA})
        row.update(solver='mic', rep=0, variant='recovery', host=host, model=model_name,
                   baseline_command=base['command'], command=command, export_command=export_command,
                   export_model_sha256=export_digest, source_difference=difference,
                   source_commit=candidate['commit'], source_tree=candidate['milp_tree'],
                   binary_sha256=candidate['binary_sha256'], version=candidate['binary_sha256'],
                   default_features=candidate['default_features'], performance_promotion_evidence=False)
        result_path, external_path = output / f'{name}.jsonl', output / f'{name}-external.jsonl'
        write_json(result_path, [row], lines=True)
        write_json(external_path, archived_refs, lines=True)
        check = inspect_results([external_path, result_path])
        errors = list(check['errors'])
        if row['auto_degraded'] or row['gap_terminated']:
            errors.append('关闭提前停止后仍记录提前退出')
        rss = row['peak_rss_bytes']
        rss_passed = type(rss) is int and rss > 0
        cases[name] = {
            'baseline_status': base['status'], 'status': row['status'], 'objective': row['objective'],
            'bound': row['bound'], 'nodes': row['nodes'], 'lp_iterations': row['lp_iterations'],
            'correctness_passed': check['passed'] and not errors, 'correctness_errors': errors,
            'proof_completed': row['status'] in PROOFS, 'rss_passed': rss_passed, 'independent_check': check,
        }
        print(name, row['status'], '正确性', cases[name]['correctness_passed'], flush=True)

    summary = {
        'performance_promotion_evidence': False, 'source_difference': difference, 'case_order': list(CASES),
        'baseline': str(baseline), 'baseline_sha256': sha(baseline),
        'external': str(EXTERNAL), 'external_sha256': sha(EXTERNAL), 'cases': cases,
        'correctness_passed': all(case['correctness_passed'] for case in cases.values()),
        'proofs_completed': sum(case['proof_completed'] for case in cases.values()),
        'rss_passed': all(case['rss_passed'] for case in cases.values()),
    }
    write_json(output / 'summary.json', summary)
    print(f"完成证明 {summary['proofs_completed']}/{len(CASES)}，"
          f"正确性 {summary['correctness_passed']}，RSS 记录 {summary['rss_passed']}", flush=True)
    return not (summary['correctness_passed'] and summary['rss_passed'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True, help='候选 Cargo 默认 release round_probe')
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    args = parser.parse_args()
    try:
        return run(args.binary.resolve(), args.source.resolve(), args.output_dir.resolve(), args.baseline.resolve())
    except (KeyError, TypeError, ValueError, OSError, RuntimeError, subprocess.SubprocessError) as error:
        print(f'定点回归无法完成：{error}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
