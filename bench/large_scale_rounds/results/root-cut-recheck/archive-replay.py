"""保存 R12 运行原件，按仓库相对模型路径重算独立复核与门禁。"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


FEATURE = 'root_cut_recheck'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def run(log, command, rounds):
    with log.open('x', encoding='utf-8') as stream:
        env = dict(os.environ, PYTHONPATH=str(rounds))
        return subprocess.run([sys.executable, *map(str, command)], stdout=stream,
                              stderr=subprocess.STDOUT, env=env).returncode


def relocate_models(raw, destination, models):
    with raw.open(encoding='utf-8') as source, destination.open('x', encoding='utf-8') as target:
        count = 0
        for line in source:
            row = json.loads(line)
            model = (raw.parent / row['model']).resolve()
            if model.parent != models or model.stem != row['model_sha256'] or digest(model) != row['model_sha256']:
                raise ValueError(f'{row["instance"]} 原模型路径或摘要不符')
            row['model'] = os.path.relpath(model, destination.parent)
            target.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')
            count += 1
    return count


def replay(work, phase, rounds, archive):
    refs = rounds / 'results'
    legacy = refs / 'external-baseline.jsonl'
    modern = refs / 'miplib2017-baseline.jsonl'
    truth = refs / 'miplib2017-sources.json'
    original = work / phase
    original_archive = archive / 'original' / phase
    shutil.copytree(original, original_archive)
    replay_dir = archive / phase
    replay_dir.mkdir()
    stem = f'{FEATURE}-{phase}'
    raw = original / f'{stem}.jsonl'
    relocated = replay_dir / raw.name
    count = relocate_models(raw, relocated, (refs / 'models').resolve())
    summary = replay_dir / f'{stem}-summary.json'
    if phase == 'diagnostics':
        calls = {'summary': [archive / 'root-cut-diagnostics.py', 'summarize', '--results', relocated,
                             '--baseline', legacy, '--output', summary]}
    else:
        check = replay_dir / f'{stem}-check.json'
        gate = replay_dir / f'{stem}-gates.json'
        calls = {
            'check': [rounds / 'check_results.py', legacy, modern, relocated, '--public-truth', truth,
                      '--output', check],
            'summary': [rounds / 'summarize_feature_trial.py', relocated, '--legacy-baseline', legacy,
                        '--modern-baseline', modern, '--public-truth', truth, '--output', summary],
            'gate': [rounds / 'gate_feature_trial.py', summary, '--output', gate],
        }
    codes = {kind: run(replay_dir / f'{kind}.log', command, rounds) for kind, command in calls.items()}
    same_results = {}
    for kind, command in calls.items():
        output = command[-1]
        first = original / output.name
        if first.exists() and output.exists():
            original_data = json.loads(first.read_text().replace(str(raw), str(relocated)).replace(
                str(original / summary.name), str(summary)))
            same_results[kind] = original_data == read(output)
        else:
            same_results[kind] = not first.exists() and not output.exists()
    original_codes = read(original / 'exit-codes.json')
    match = all(codes[kind] == original_codes[kind] for kind in calls) and all(same_results.values())
    return {
        'rows': count, 'changed_row_field': 'model',
        'original_sha256': {path.name: digest(path) for path in sorted(original_archive.iterdir())},
        'archive_sha256': {path.name: digest(path) for path in sorted(replay_dir.iterdir())},
        'original_exit_codes': original_codes, 'replay_exit_codes': codes,
        'replay_matches_original': same_results, 'matched': match,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('work', 'rounds', 'validation', 'output'):
        parser.add_argument('--' + flag, type=Path, required=True)
    args = parser.parse_args()
    work, rounds, archive = (path.resolve() for path in (args.work, args.rounds, args.output))
    archive.mkdir()
    (archive / 'original').mkdir()
    for name in ('execution.sh', 'record-build.py', 'root-cut-diagnostics.py', 'input-recovery.json',
                 'external-reference-check.json', 'archive-replay.py'):
        shutil.copyfile(work / name, archive / name)
    shutil.copytree(work / 'builds', archive / 'original' / 'builds')
    validation = args.validation.resolve()
    validation_archive = archive / 'original' / 'validation'
    validation_archive.mkdir()
    for path in (*sorted(validation.glob('0*.log')), validation / 'run.sh', validation / 'status.tsv'):
        shutil.copyfile(path, validation_archive / path.name)
    reports = {phase: replay(work, phase, rounds, archive) for phase in ('screen', 'paired', 'diagnostics')
               if (work / phase / 'exit-codes.json').exists()}
    with (archive / 'archive-relocation.json').open('x', encoding='utf-8') as stream:
        support = {str(path.relative_to(archive)): digest(path)
                   for folder in (archive / 'original' / 'builds', validation_archive)
                   for path in sorted(folder.iterdir())}
        json.dump({'phases': reports, 'support_files_sha256': support}, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    print(json.dumps({phase: data['matched'] for phase, data in reports.items()}, ensure_ascii=False))
    return not reports or not all(entry['matched'] for entry in reports.values())


if __name__ == '__main__':
    sys.exit(main())
