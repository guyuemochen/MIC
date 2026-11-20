#!/usr/bin/env python3
"""生成与调参样本分离的 LP 扫描复核集，不读写原基准或 .opt 真值。"""
import argparse
import hashlib
import json
import random
from pathlib import Path

CASES = [
    ('scan-sparse-newseed', 2400, 500, 0.05, 20261111),
    ('scan-sparse-tall', 1000, 1500, 0.08, 20261112),
    ('scan-medium', 1600, 800, 0.30, 20261113),
    ('scan-dense-wide', 2000, 500, 1.0, 20261114),
]


def generate(directory):
    directory.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name, n, m, density, seed in CASES:
        rng = random.Random(seed)
        lines = [f'name {name}', 'sense max', 'vars']
        lines += [f'x{j} continuous' for j in range(n)]
        lines += ['obj', ' '.join(f'{rng.randint(3, 19)} x{j}' for j in range(n)), 'cons']
        covered = set()
        for i in range(m):
            columns = sorted(rng.sample(range(n), max(1, round(n * density))))
            covered.update(columns)
            row = ' '.join(f'{rng.randint(1, 9)} x{j}' for j in columns)
            lines.append(f'r{i}: {row} <= {rng.randint(50, 400)}')
        assert len(covered) == n, '拒绝有游离无界列的样本'
        lines.append('end')
        payload = ('\n'.join(lines) + '\n').encode()
        (directory / (name + '.milp')).write_bytes(payload)
        manifest.append(dict(name=name, variables=n, rows=m, density=density, seed=seed, sha256=hashlib.sha256(payload).hexdigest()))
    (directory / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', type=Path)
    generate(p.parse_args().directory)
