#!/usr/bin/env python3
"""核对六批 Gurobi 对照的保留文件，按需恢复原路径上的重复副本。"""

import argparse
import csv
from hashlib import sha256
from pathlib import Path


HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "gurobi-duplicate-inputs.tsv"


def checked_path(root: Path, name: str) -> Path:
    relative = Path(name)
    path = root / relative
    if relative.is_absolute() or ".." in relative.parts or path.resolve() != path:
        raise ValueError(f"归档清单中的路径不安全：{name}")
    return path


def checked_bytes(path: Path, size: int, digest: str) -> bytes:
    data = path.read_bytes()
    if len(data) != size or sha256(data).hexdigest() != digest:
        raise ValueError(f"文件大小或 SHA-256 与清单不符：{path}")
    return data


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results",
                        help="历史结果目录，默认使用本仓库的结果目录")
    parser.add_argument("--restore", action="store_true", help="恢复缺失的副本，保留现有文件")
    args = parser.parse_args()
    root = args.results_dir.resolve(strict=True)
    retained = {}
    missing = []
    present = 0
    with MANIFEST.open(encoding="utf-8", newline="") as index:
        for row in csv.DictReader(index, delimiter="\t"):
            source = checked_path(root, row["retained"])
            target = checked_path(root, row["original"])
            expected = (int(row["bytes"]), row["sha256"])
            if source not in retained:
                checked_bytes(source, *expected)
                retained[source] = expected
            if retained[source] != expected:
                raise ValueError(f"同一保留文件在清单中有不同摘要：{source}")
            if target.exists():
                checked_bytes(target, *expected)
                present += 1
            else:
                missing.append((source, target))

    if args.restore:
        for source, target in missing:
            data = checked_bytes(source, *retained[source])
            target.parent.mkdir(parents=True, exist_ok=True)
            output = target.open("xb")
            try:
                with output:
                    output.write(data)
            except BaseException:
                target.unlink()
                raise

    print(f"已核对 {len(retained)} 份保留文件和 {present} 份现有副本。")
    print(f"{'已恢复' if args.restore else '可恢复'} {len(missing)} 份副本。")


if __name__ == "__main__":
    main()
