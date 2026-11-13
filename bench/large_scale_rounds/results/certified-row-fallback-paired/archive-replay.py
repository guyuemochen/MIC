"""保留 R11 原始证据并重放仓库相对路径副本。"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

repo = Path('/Users/guyuemochen/Desktop/projects/MIC')
rounds = repo / 'bench/large_scale_rounds'
results = rounds / 'results'
raw = Path('/tmp/mic-r11-paired')
destination = results / 'certified-row-fallback-paired'
archive = Path('/tmp/mic-r11-paired-repo-archive')
stem = 'certified_row_fallback-paired'
raw_jsonl = raw / f'{stem}.jsonl'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):
    return json.loads(path.read_text())

def run(log, command):
    with log.open('x', encoding='utf-8') as stream:
        return subprocess.run([sys.executable, *map(str, command)], stdout=stream,
                              stderr=subprocess.STDOUT).returncode

def verify_and_relocate():
    original_rows = [json.loads(line) for line in raw_jsonl.open(encoding='utf-8')]
    assert len(original_rows) == 120
    archive.mkdir()
    (archive / 'original').mkdir()
    shutil.copyfile(Path(__file__), archive / 'archive-replay.py')
    original_hashes = {}
    for item in sorted(raw.iterdir()):
        assert item.is_file()
        original_hashes[item.name] = sha(item)
        shutil.copyfile(item, archive / 'original' / item.name)
    replay_workspace = Path('/tmp/mic-r11-paired-replay')
    replay_workspace.mkdir()
    (replay_workspace / 'models').symlink_to(results / 'models', target_is_directory=True)
    replay = replay_workspace / destination.name
    replay.mkdir()
    replay_data = replay / raw_jsonl.name
    with replay_data.open('x', encoding='utf-8') as target:
        for row in original_rows:
            model = (raw / row['model']).resolve()
            assert model.parent == (results / 'models').resolve()
            assert model.stem == row['model_sha256'] and sha(model) == row['model_sha256']
            target.write(json.dumps(dict(row, model=f'../models/{model.name}'),
                                    ensure_ascii=False, allow_nan=False) + '\n')
    shutil.copyfile(replay_data, archive / replay_data.name)
    legacy = results / 'external-baseline.jsonl'
    modern = results / 'miplib2017-baseline.jsonl'
    truth = results / 'miplib2017-sources.json'
    refs = ['--legacy-baseline', legacy, '--modern-baseline', modern, '--public-truth', truth]
    rc = {}
    outputs = {kind: archive / f'{stem}-archive-replay-{kind}.json' for kind in ('check', 'summary', 'gates')}
    rc['check'] = run(archive / 'archive-replay-check.log', [
        rounds / 'check_results.py', legacy, modern, replay_data, '--public-truth', truth,
        '--output', outputs['check'],
    ])
    rc['summary'] = run(archive / 'archive-replay-summary.log', [
        rounds / 'summarize_feature_trial.py', replay_data, *refs, '--output', outputs['summary'],
    ])
    rc['gates'] = run(archive / 'archive-replay-gate.log', [
        rounds / 'gate_feature_trial.py', outputs['summary'], '--output', outputs['gates'],
    ])
    gate_original = read(raw / f'{stem}-gates.json')
    gate_replay = read(outputs['gates'])
    gate_original.pop('input_summary')
    gate_replay.pop('input_summary')
    check_original = read(raw / f'{stem}-check.json')
    check_replay = read(outputs['check'])
    summary_original = read(raw / f'{stem}-summary.json')
    summary_replay = read(outputs['summary'])
    comparison = {
        'check_exact': check_original == check_replay,
        'summary_exact': summary_original == summary_replay,
        'gates_exact_except_input_summary': gate_original == gate_replay,
    }
    ledger = {
        'destination': str(destination.relative_to(repo)),
        'changed_field': 'model', 'raw_sha256': original_hashes[raw_jsonl.name],
        'repo_relative_sha256': sha(archive / raw_jsonl.name),
        'original_sha256': original_hashes,
        'archive_replay_sha256': {kind: sha(path) for kind, path in outputs.items()},
        'archive_replay_exit_codes': rc, 'archive_replay_comparison': comparison,
    }
    with (archive / 'archive-relocation.json').open('x', encoding='utf-8') as target:
        json.dump(ledger, target, ensure_ascii=False, indent=2)
        target.write('\n')
    print(json.dumps(ledger, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    verify_and_relocate()
