#!/usr/bin/env python3
"""冻结 MIPLIB 2017 留出名单，按冻结清单下载与转换，无求解步骤。"""
import argparse
import gzip
import hashlib
import html
import json
import math
from pathlib import Path
import re
import sys
from urllib.parse import urljoin, urlparse
from urllib.request import urlopen


BASE = 'https://miplib.zib.de/'
PREFIX = 'MIC-MIPLIB2017-R6-20261112:'
ASSETS = {
    'benchmark_list': 'downloads/benchmark-v2.test',
    'benchmark_metadata': 'tag_benchmark.html',
    'solution_file': 'downloads/miplib2017-v37.solu',
}
TIERS = ('small_binary', 'small_mixed_binary', 'small_general_integer', 'large')
SIZE_KEYS = ('variables', 'binaries', 'integers', 'continuous', 'rows', 'nonzeros')
SIZE_LABELS = dict(zip(('Variables', 'Binaries', 'Integers', 'Continuous', 'Constraints', 'Nonzeroes'),
                       SIZE_KEYS, strict=True))


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write_new(path, data):
    with path.open('xb') as stream:
        stream.write(data)


def download(url, path):
    if not path.exists():
        for attempt in range(3):
            try:
                with urlopen(url, timeout=45) as response:
                    data = response.read()
                break
            except OSError:
                if attempt == 2:
                    raise
        write_new(path, data)
    return path.read_bytes()


def save_json(path, value, copy):
    payload = (json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n').encode('utf-8')
    write_new(path, payload)
    if copy is not None:
        write_new(copy, payload)


def table_rows(document):
    for row in re.findall(r'<tr\b[^>]*>(.*?)</tr>', document, re.S):
        yield re.findall(r'<td\b[^>]*>(.*?)</td>', row, re.S)


def plain(cell):
    return html.unescape(re.sub(r'<[^>]+>', '', cell)).strip()


def list_names(data):
    lines = data.decode().splitlines()
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+\.mps\.gz', name) for name in lines):
        raise ValueError('官方 v2 名单存在非模型文件名')
    names = {name.removesuffix('.mps.gz') for name in lines}
    if len(lines) != 240 or len(names) != 240:
        raise ValueError('官方 v2 名单须有 240 个唯一实例')
    return names


def benchmark_metadata(data):
    metadata = {}
    for cells in table_rows(data.decode()):
        if not cells or 'instance_details_' not in cells[0]:
            continue
        name = plain(cells[0])
        detail = f'instance_details_{name}.html'
        if name in metadata or len(cells) != 12 or detail not in cells[0]:
            raise ValueError(f'官方元数据行结构异常或名字重复：{name}')
        sizes = dict(zip(SIZE_KEYS, (int(plain(cell)) for cell in cells[2:8]), strict=True))
        if any(value < 0 for value in sizes.values()) or (
                sizes['variables'] != sum(sizes[key] for key in ('binaries', 'integers', 'continuous'))):
            raise ValueError(f'官方元数据变量计数不一致：{name}')
        metadata[name] = dict(sizes, detail_url=urljoin(BASE, detail))
    return metadata


def tier_of(model):
    if model['binaries'] + model['integers'] == 0:
        return None
    n, m, nnz = (model[key] for key in ('variables', 'rows', 'nonzeros'))
    if n <= 10_000 and m <= 10_000 and nnz <= 200_000:
        if model['integers']:
            return 'small_general_integer'
        return 'small_mixed_binary' if model['continuous'] else 'small_binary'
    if (n > 10_000 or m > 10_000) and n <= 50_000 and m <= 50_000 and nnz <= 500_000:
        return 'large'
    return None


def selection_from_assets(assets, data):
    names = list_names(data['benchmark_list'])
    metadata = benchmark_metadata(data['benchmark_metadata'])
    if set(metadata) != names:
        raise ValueError(f'官方 v2 名单与 HTML 名字不同：{sorted(set(metadata) ^ names)}')
    eligible = {tier: [] for tier in TIERS}
    for name in names:
        model = metadata[name]
        tier = tier_of(model)
        if tier:
            rank = sha((PREFIX + name).encode('utf-8'))
            eligible[tier].append({
                'instance': name, 'tier': tier, 'selection_sha256': rank, 'official_metadata': model,
            })
    chosen = []
    for tier in TIERS:
        candidates = sorted(eligible[tier], key=lambda item: (item['selection_sha256'], item['instance']))
        if len(candidates) < 2:
            raise ValueError(f'分层 {tier} 少于两个候选')
        chosen.extend(candidates[:2])
    return {
        'source': BASE, 'benchmark_count': len(names), 'prefix': PREFIX, 'per_tier': 2,
        'small_caps': {'variables': 10_000, 'rows': 10_000, 'nonzeros': 200_000},
        'large_caps': {'variables': 50_000, 'rows': 50_000, 'nonzeros': 500_000},
        'assets': assets, 'eligible_counts': {tier: len(eligible[tier]) for tier in TIERS},
        'instances': chosen,
    }


def select(directory, copy):
    path = directory / 'selection.json'
    if path.exists() or (copy is not None and copy.exists()):
        raise FileExistsError('已有冻结名单')
    directory.mkdir(parents=True, exist_ok=True)
    assets, data = {}, {}
    for key, relative in ASSETS.items():
        url = urljoin(BASE, relative)
        name = Path(relative).name
        content = download(url, directory / name)
        assets[key] = {'url': url, 'file': name, 'sha256': sha(content)}
        data[key] = content
    assets['solution_file']['release_date'] = '2026-09-24'
    frozen = selection_from_assets(assets, data)
    save_json(path, frozen, copy)
    for item in frozen['instances']:
        print(item['tier'], item['instance'], item['selection_sha256'], flush=True)


def solution_info(data, names):
    solutions = {}
    statuses = {'=opt=': 'optimal', '=best=': 'best_known', '=inf=': 'infeasible',
                '=unbd=': 'unbounded', '=unkn=': 'unknown'}
    for line in data.decode().splitlines():
        status, name, *value = line.split()
        if name not in names:
            continue
        if name in solutions or status not in statuses or len(value) != (status in ('=opt=', '=best=')):
            raise ValueError(f'官方解文件行异常：{line}')
        number = float(value[0]) if value else None
        if number is not None and not math.isfinite(number):
            raise ValueError(f'官方解文件的目标非有限：{line}')
        solutions[name] = {
            'solution_status': statuses[status], 'solution_value': number,
            'solution_value_text': value[0] if value else None,
        }
    if set(solutions) != names:
        raise ValueError(f'官方解文件缺少冻结算例：{sorted(names - set(solutions))}')
    return solutions


def instance_url(details, name, metadata):
    counts = {}
    for cells in table_rows(details):
        if len(cells) == 3 and (label := plain(cells[0])) in SIZE_LABELS:
            key = SIZE_LABELS[label]
            counts[key] = int(plain(cells[1]))
    if counts != {key: metadata[key] for key in SIZE_KEYS}:
        raise ValueError(f'算例详情与基准元数据的原始计数不同：{counts}')
    filename = f'{name}.mps.gz'
    links = {urljoin(metadata['detail_url'], html.unescape(href))
             for href in re.findall(r'href=["\']([^"\']+)["\']', details)
             if urlparse(href).path.rsplit('/', 1)[-1] == filename}
    if len(links) != 1:
        raise ValueError(f'官方详情页缺少唯一的 {filename} 下载链接')
    url = links.pop()
    parsed = urlparse(url)
    if (parsed.scheme, parsed.netloc) != ('https', 'miplib.zib.de'):
        raise ValueError(f'MPS 下载链接不在官方 HTTPS 域名：{url}')
    return url


def unsupported_details(mps):
    import highspy

    reader = highspy.Highs()
    reader.setOptionValue('output_flag', False)
    if reader.readModel(str(mps)) != highspy.HighsStatus.kOk:
        return {}
    model = reader.getLp()
    return {
        'objective_offset': model.offset_,
        'negative_finite_lower_bounds': sum(math.isfinite(v) and v < 0 for v in model.col_lower_),
        'negative_infinite_lower_bounds': sum(v == -math.inf for v in model.col_lower_),
    }


def fetch_instance(directory, item, solution, convert):
    name = item['instance']
    metadata = item['official_metadata']
    result = dict(item, group='holdout', **solution)
    try:
        detail_file = directory / f'instance_details_{name}.html'
        details = download(metadata['detail_url'], detail_file)
        result.update(detail_file=detail_file.name, detail_sha256=sha(details))
        url = instance_url(details.decode(), name, metadata)
        result['url'] = url
        archive = directory / f'{name}.mps.gz'
        compressed = download(url, archive)
        result.update(gzip=archive.name, gzip_sha256=sha(compressed))
        raw = gzip.decompress(compressed)
        mps = directory / f'{name}.mps'
        if mps.exists():
            if mps.read_bytes() != raw:
                raise ValueError(f'现有 MPS 与官方压缩文件不同：{mps}')
        else:
            write_new(mps, raw)
        result.update(mps=mps.name, mps_sha256=sha(raw))
    except (OSError, ValueError, EOFError) as error:
        return dict(result, included=False, stage='download', reason=str(error))
    milp = directory / f'{name}.milp'
    if milp.exists():
        raise FileExistsError(f'已有转换文件：{milp}')
    try:
        convert(mps, milp)
    except (OSError, ValueError) as error:
        return dict(result, included=False, stage='conversion', reason=str(error),
                    conversion_details=unsupported_details(mps))
    return dict(result, included=True, milp=milp.name, milp_sha256=sha(milp.read_bytes()))


def convert_frozen(directory, copy):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'experiments'))
    from import_mps import convert

    path = directory / 'sources.json'
    if path.exists() or (copy is not None and copy.exists()):
        raise FileExistsError('已有下载结果清单')
    selection_path = directory / 'selection.json'
    selection_data = selection_path.read_bytes()
    frozen = json.loads(selection_data)
    data = {}
    for key, asset in frozen['assets'].items():
        content = (directory / asset['file']).read_bytes()
        if sha(content) != asset['sha256']:
            raise ValueError(f'冻结数据的哈希不同：{asset["file"]}')
        data[key] = content
    if frozen != selection_from_assets(frozen['assets'], data):
        raise ValueError('当前冻结名单与按原始数据重算的名单不同')
    chosen = frozen['instances']
    solutions = solution_info(data['solution_file'], {item['instance'] for item in chosen})
    manifest = {
        'source': BASE, 'selection': selection_path.name, 'selection_sha256': sha(selection_data),
        'assets': frozen['assets'], 'instances': [],
    }
    for item in chosen:
        result = fetch_instance(directory, item, solutions[item['instance']], convert)
        manifest['instances'].append(result)
        print(result['instance'], 'included' if result['included'] else result['reason'], flush=True)
    save_json(path, manifest, copy)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('select', 'download'))
    parser.add_argument('directory', type=Path)
    parser.add_argument('--copy', type=Path, help='额外保存同内容的清单，路径须未占用')
    args = parser.parse_args()
    action = select if args.phase == 'select' else convert_frozen
    action(args.directory, args.copy)


if __name__ == '__main__':
    main()
