"""记录或验证 R12 构建、源码及输入摘要。"""
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import tomllib


FEATURE = 'root_cut_recheck'
TIMING_TOOLS = ('run_feature_trial.py', 'compare.py', 'run_root_work.py', 'run_row_products.py',
                'run_proof_diagnostics.py', 'check_results.py', 'summarize_feature_trial.py',
                'gate_feature_trial.py')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command(*parts):
    return subprocess.check_output(parts, text=True).strip()


def require(condition, message):
    if not condition:
        sys.exit(message)


def dump(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2)
        stream.write('\n')


def build_info(source, rounds, archive, kind):
    declared = tomllib.loads((source / 'milp/Cargo.toml').read_text())['features']

    def closure(names):
        enabled = set(names)
        pending = list(names)
        while pending:
            for item in declared[pending.pop()]:
                if item not in enabled:
                    enabled.add(item)
                    pending.append(item)
        return enabled

    defaults = declared['default']
    effective_defaults = closure(defaults)
    require(not {FEATURE, 'proof_profile'} & effective_defaults, '试验特性或证明诊断已在生产默认中')
    entries, builds = {}, {}
    for variant in ('base', 'candidate'):
        explicit = [FEATURE] if variant == 'candidate' else []
        if kind == 'diagnostics':
            explicit.append('proof_profile')
        enabled = closure(defaults + explicit)
        require(enabled == effective_defaults | set(explicit), '构建出现额外 Cargo 特性')
        binary = archive / 'binaries' / f'round_probe_{kind}_{variant}'
        log = archive / 'builds' / f'{kind}-{variant}.log'
        entries[variant] = {'binary': str(binary), 'features': explicit}
        builds[variant] = dict(entries[variant], binary_sha256=digest(binary),
                               effective_features=sorted(enabled), build_log=log.name,
                               build_log_sha256=digest(log))
    require(builds['base']['binary_sha256'] != builds['candidate']['binary_sha256'],
            '两组二进制摘要不应相同')
    repository = ['git', '-C', str(source)]
    refs = rounds / 'results'
    return entries, {
        'source': str(source), 'source_commit': command(*repository, 'rev-parse', 'HEAD'),
        'source_tree': command(*repository, 'rev-parse', 'HEAD:milp'),
        'repository_tree': command(*repository, 'rev-parse', 'HEAD^{tree}'),
        'default_features': defaults, 'builds': builds, 'host': platform.platform(),
        'cargo': command('cargo', '-vV'), 'rustc': command('rustc', '-vV'), 'python': sys.version,
        'refs_sha256': {name: digest(refs / name) for name in (
            'external-baseline.jsonl', 'miplib2017-baseline.jsonl', 'miplib2017-sources.json',
            'public-sources.json')},
        'tools_sha256': {name: digest(rounds / name) for name in TIMING_TOOLS},
    }


def main():
    source, rounds, archive = (Path(value).resolve() for value in sys.argv[1:4])
    phase = sys.argv[4]
    kind = 'timing' if phase == 'verify' else phase
    entries, builds = build_info(source, rounds, archive, kind)
    prefix = 'diagnostics-' if kind == 'diagnostics' else ''
    manifest = archive / 'builds' / f'{prefix}{FEATURE}-binaries.json'
    record = archive / 'builds' / f'{prefix}{FEATURE}-builds.json'
    if phase != 'verify':
        dump(manifest, entries)
    builds['manifest_sha256'] = digest(manifest)
    if phase == 'verify':
        require(json.loads(manifest.read_text()) == entries, '两阶段的二进制清单不一致')
        require(json.loads(record.read_text()) == builds, '两阶段的源码、二进制或工具摘要不一致')
    else:
        dump(record, builds)


if __name__ == '__main__':
    main()
