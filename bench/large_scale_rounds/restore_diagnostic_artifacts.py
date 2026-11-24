#!/usr/bin/env python3
"""核对 bnatt500 和 p0201 压缩诊断证据，按需恢复原件。"""

import argparse
import csv
import gzip
from hashlib import sha256
from pathlib import Path


HERE = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results",
                        help="历史结果目录，默认使用本仓库的结果目录")
    parser.add_argument("--restore", action="store_true", help="恢复缺失的原件，保留现有文件")
    args = parser.parse_args()
    root = args.results_dir.resolve(strict=True)
    missing = []
    present = 0
    with (HERE / "diagnostic-artifacts.tsv").open(encoding="utf-8", newline="") as index:
        for row in csv.DictReader(index, delimiter="\t"):
            relative = Path(row["original"])
            target = root / relative
            source = target.with_name(target.name + ".gz")
            if relative.is_absolute() or ".." in relative.parts or target.resolve() != target or source.resolve() != source:
                raise ValueError(f"归档清单中的路径不安全：{relative}")
            data = gzip.decompress(source.read_bytes())
            if len(data) != int(row["bytes"]) or sha256(data).hexdigest() != row["sha256"]:
                raise ValueError(f"解压后大小或 SHA-256 与归档清单不符：{source}")
            if target.exists():
                if target.read_bytes() != data:
                    raise ValueError(f"现有原件与归档内容不同：{target}")
                present += 1
            else:
                missing.append((target, data))

    if args.restore:
        for target, data in missing:
            output = target.open("xb")
            try:
                with output:
                    output.write(data)
            except BaseException:
                target.unlink()
                raise

    print(f"已核对 {present + len(missing)} 份压缩诊断证据，{'已恢复' if args.restore else '可恢复'} {len(missing)} 份原件。")


if __name__ == "__main__":
    main()
