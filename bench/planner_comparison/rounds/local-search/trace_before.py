"""在上一轮源码的临时副本中记录每次提交的工作分数；使用固定步数，可独立重放。"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH))
from mic_model import compile_model
from run import SEEDS, digest

# 只观察每个禁忌搜索步骤，沿用原来的采样、判分及提交路径。
ENTRY = """use planner::business::{optimize, wire};
fn main() {
    let input = std::fs::read_to_string(std::env::args().nth(1).unwrap()).unwrap();
    let (m, cfg) = wire::parse(&input).unwrap();
    let best = optimize(&m, &cfg).unwrap();
    println!("{} {} {}", best.steps, best.score.hard, best.score.soft);
}
"""
STEPS = 10000


def trace(crate):
    previous = json.loads((BENCH / "rounds/paired-scoring/results/environment.json").read_text())
    with TemporaryDirectory(prefix="mic-trace-") as name:
        tmp = Path(name)
        shutil.copytree(crate / "src", tmp / "src")
        for filename in ("Cargo.toml", "Cargo.lock"):
            shutil.copy2(crate / filename, tmp / filename)
        sources = {tmp / Path(key).relative_to("planner"): sha
                   for key, sha in previous["workspace_source_sha256"].items()
                   if key.startswith("planner/src/")}
        if any(digest(path) != sha for path, sha in sources.items()):
            patch = Path(__file__).resolve().with_name("before.patch")
            subprocess.run(["git", "apply", "--no-index", "-p2", str(patch)], cwd=tmp, check=True)
        for path, sha in sources.items():
            if digest(path) != sha:
                raise ValueError("源码与 109 版本不符：" + str(path))
        search = tmp / "src/business/search.rs"
        code = search.read_text()
        marker = "        if work.scorer.score() < best {"
        if code.count(marker) != 1:
            raise ValueError("没有找到唯一的提交点")
        search.write_text(code.replace(marker, '        eprintln!("{} {}", work.scorer.score().0, work.scorer.score().1);\n' + marker))
        (tmp / "examples").mkdir()
        (tmp / "examples/trace.rs").write_text(ENTRY)
        subprocess.run(["cargo", "build", "--offline", "--locked", "--release", "--example", "trace"],
                       cwd=tmp, check=True)
        wire = compile_model(BENCH / "instances/timetable.txt")
        model = tmp / "timetable.wire"
        model.write_text(wire)
        if digest(model) != previous["mic_wire_sha256"]["timetable"]:
            raise ValueError("输入与 109 版本不符")
        report = []
        for seed in SEEDS:
            model.write_text(wire + f"\nC {seed} {STEPS} -1 1 tabu bf\n")
            result = subprocess.run([tmp / "target/release/examples/trace", model], check=True,
                                    text=True, capture_output=True)
            scores = [tuple(map(int, line.split())) for line in result.stderr.splitlines()]
            steps, hard, soft = map(int, result.stdout.split())
            if not len(scores) == steps == STEPS:
                raise ValueError("跟踪步数不一致")
            report.append(dict(seed=seed, steps=steps, best=[hard, soft],
                               first_worse=next((i for i, rank in enumerate(scores, 1)
                                                 if rank > (0, 109)), None),
                               feasible_steps=sum(h == 0 for h, _ in scores),
                               steps_equal_to_initial=sum(rank == (0, 109) for rank in scores),
                               checkpoints={i: scores[i-1] for i in (1, 10, 100, 1000, STEPS)}))
        return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("crate", nargs="?", default=BENCH.parents[1] / "planner", type=Path,
                        help="本轮或上一轮的 planner 目录")
    print(json.dumps(trace(parser.parse_args().crate), indent=2, ensure_ascii=False))
