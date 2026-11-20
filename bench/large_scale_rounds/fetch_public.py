#!/usr/bin/env python3
"""下载固定的 MIPLIB3 清单并转换，记录原文件、转换文件和排除理由。"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import urlopen

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'experiments'))
from import_mps import convert

ARCHIVE = 'https://miplib2010.zib.de/miplib3/miplib3'
CASES = {
    'regression': ('p0033', 'p0201', 'enigma', 'stein27'),
    'holdout': ('p0282', 'p0548', 'lseu', 'air03'),
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def fetch(directory, group, name):
    url = f'{ARCHIVE}/{name}.mps.gz'
    item = {'instance': name, 'group': group, 'url': url}
    try:
        with urlopen(url, timeout=30) as response:
            compressed = response.read()
        (directory / f'{name}.mps.gz').write_bytes(compressed)
        item['gzip_sha256'] = sha(compressed)
        original = gzip.decompress(compressed)
        mps = directory / f'{name}.mps'
        mps.write_bytes(original)
        item.update(mps=mps.name, mps_sha256=sha(original))
    except (OSError, ValueError, EOFError) as error:
        return dict(item, included=False, stage='download', reason=str(error))
    milp = directory / f'{name}.milp'
    try:
        convert(mps, milp)
    except (OSError, ValueError) as error:
        return dict(item, included=False, stage='conversion', reason=str(error))
    return dict(item, included=True, milp=milp.name, milp_sha256=sha(milp.read_bytes()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    args.directory.mkdir(parents=True)
    manifest = {'source': ARCHIVE, 'instances': []}
    for group, names in CASES.items():
        for name in names:
            item = fetch(args.directory, group, name)
            manifest['instances'].append(item)
            print(name, 'included' if item['included'] else item['reason'], flush=True)
    (args.directory / 'sources.json').write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
