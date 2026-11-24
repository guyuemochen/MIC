"""在临时副本恢复上一轮 MIC 源码，逐文件核对后编译旧版基准程序。"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory


ROUND = Path(__file__).resolve().parent
BENCH = ROUND.parents[1]
ROOT = BENCH.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path, help="保存旧版可执行文件的新路径")
    parser.add_argument("--round", choices=("paired-scoring", "local-search"),
                        default="paired-scoring", help="还原指定轮次开始前的版本")
    args = parser.parse_args()
    dest = args.out
    if dest.exists():
        parser.error("输出文件已存在")
    stages = [(ROUND.parent / "local-search/before.patch", ROUND / "results/environment.json",
               "workspace_source_sha256")]
    if args.round == "paired-scoring":
        stages.append((ROUND / "before.patch", BENCH / "results/2026-11-24/environment.json",
                       "source_sha256"))
    with TemporaryDirectory(prefix="mic-rebuild-before-") as tmp:
        root = Path(tmp)
        crate = root / "planner"
        shutil.copytree(ROOT / "planner/src", crate / "src")
        (crate / "examples").mkdir()
        for name in ("Cargo.lock", "Cargo.toml", "examples/compare_worker.rs"):
            shutil.copy2(ROOT / "planner" / name, crate / name)
        for patch, manifest, key in stages:
            frozen = {name: expected for name, expected in json.loads(manifest.read_text())[key].items()
                      if name.startswith("planner/")}
            def matches(name):
                path = root / name
                return path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() == frozen[name]
            if not all(matches(name) for name in frozen):
                subprocess.run(["git", "apply", "--no-index", str(patch)], cwd=root, check=True)
            for name in frozen:
                if not matches(name):
                    raise ValueError("还原源码指纹不符：" + name)
            print(f"{patch.parent.name} 旧版 {len(frozen)} 个 Rust 文件均与原始归档一致。", flush=True)
        target = root / "target"
        subprocess.run(["cargo", "build", "--release", "--locked", "--offline", "--example",
                        "compare_worker", "--target-dir", str(target)], cwd=crate, check=True)
        name = "compare_worker.exe" if os.name == "nt" else "compare_worker"
        worker = target / "release/examples" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("xb") as out:
            out.write(worker.read_bytes())
        shutil.copymode(worker, dest)
    print("SHA-256", hashlib.sha256(dest.read_bytes()).hexdigest(), dest)


if __name__ == "__main__":
    main()
