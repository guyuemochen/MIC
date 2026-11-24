"""依次运行两个常驻进程，核查期限内首个可行解和最终解，保留逐次原始记录。"""
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
from queue import Queue, Empty
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from threading import Thread

from cases import CASES, check, is_complete
from mic_model import compile_model

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OPTA = HERE / "opta"
SEEDS = [7, 19, 43, 73, 101]
JVM_FLAGS = ["-XX:ActiveProcessorCount=1", "-XX:+UseSerialGC", "-Xms512m", "-Xmx512m"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def output(command):
    return subprocess.check_output(command, cwd=ROOT, stderr=subprocess.STDOUT,
                                   text=True).strip()


class Worker:
    def __init__(self, cmd, log):
        self.process = subprocess.Popen(cmd, cwd=ROOT, text=True, bufsize=1,
                                        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=log)
        self.lines = Queue()
        def collect():
            for line in self.process.stdout:
                self.lines.put(line.rstrip("\n"))
            self.lines.put(None)
        Thread(target=collect, daemon=True).start()

    def solve(self, path, seed, budget):
        self.process.stdin.write("%s %d %d\n" % (path, seed, budget))
        self.process.stdin.flush()
        while True:
            try:
                line = self.lines.get(timeout=max(60, 3 * budget / 1000))
            except Empty as error:
                raise RuntimeError("求解进程超时") from error
            if line is None:
                raise RuntimeError("求解进程退出，退出码 %s；详见 stderr" % self.process.poll())
            if line.startswith("RESULT "):
                return json.loads(line[7:])
            print(line, file=sys.stderr, flush=True)

    def close(self):
        try:
            self.process.stdin.close()
            self.process.wait(timeout=3)
        except (BrokenPipeError, subprocess.TimeoutExpired):
            self.process.kill()
            self.process.wait()
        finally:
            self.process.stdout.close()


def validate(path, r, budget):
    complete = is_complete(path, r["values"], r["lists"])
    if complete != r["complete"]:
        raise ValueError("完整性标志不一致")
    if complete:
        try:
            independent = check(path, r["values"], r["lists"])
        except ValueError as e:
            # Java 工序搜索允许访问成环候选，报告时将完整的循环等待记为硬违反。
            if path.stem != "factory" or "循环等待" not in str(e):
                raise
            independent = (1, 0)
        if independent != (r["hard"], r["soft"]):
            raise ValueError("最终分数不符：独立复算 %s，求解器 %s" % (independent, r))
    first = r["first_ms"]
    if (first is not None) != (complete and r["hard"] == 0):
        raise ValueError("首次可行时间与历史最好解不一致")
    if not math.isfinite(r["elapsed_ms"]) or r["elapsed_ms"] < 0:
        raise ValueError("返回耗时无效")
    if first is not None:
        if not math.isfinite(first) or not 0 <= first <= min(budget, r["elapsed_ms"]) + 0.001:
            raise ValueError("首次可行时刻越过同一截止时间")
        hard, soft = check(path, r["first_values"], r["first_lists"])
        if hard != 0 or r["soft"] > soft:
            raise ValueError("首个可行快照无法复核，或最终解变差")
        r["first_soft_checked"] = soft
    r["checked_feasible"] = complete and r["hard"] == 0


def build():
    subprocess.run(["cargo", "build", "--release", "--manifest-path", ROOT / "planner/Cargo.toml",
                    "--example", "compare_worker", "--target-dir", ROOT / "planner/target"], check=True)
    subprocess.run(["mvn", "-q", "-f", OPTA / "pom.xml", "compile", "dependency:build-classpath",
                    "-Dmdep.outputFile=target/classpath"], check=True)
    classpath = (OPTA / "target/classpath").read_text().strip()
    optajar = next(Path(p) for p in classpath.split(os.pathsep) if "optaplanner-core-impl-" in p)
    mic = ROOT / "planner/target/release/examples/compare_worker"
    java = shutil.which("java")
    return {"MIC": [str(mic)],
            "OptaPlanner": [java, *JVM_FLAGS, "-cp",
                            str(OPTA / "target/classes") + os.pathsep + classpath, "benchmark.Worker"]
            }, {"mic_binary": digest(mic), "optaplanner_jar": digest(optajar)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--budget-ms", type=int, default=10000)
    parser.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    parser.add_argument("--warmup-ms", type=int, default=2000)
    args = parser.parse_args()
    if args.budget_ms <= 0 or args.warmup_ms <= 0 or len(set(args.seeds)) != len(args.seeds):
        parser.error("预算必须为正，种子不能重复")
    cmd, artifact_hashes = build()
    dest = args.out.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        raw = stack.enter_context((dest / "runs.jsonl").open("x", encoding="utf8"))
        workers = {}
        for engine, command in cmd.items():
            log = stack.enter_context((dest / (engine + ".stderr.log")).open("w", encoding="utf8"))
            workers[engine] = Worker(command, log)
            stack.callback(workers[engine].close)
        tmp = Path(stack.enter_context(TemporaryDirectory(prefix="mic-planner-")))
        inputs = {k: HERE / "instances" / (k + ".txt") for k in CASES}
        wires = {k: tmp / (k + ".wire") for k in CASES}
        for name, case in inputs.items():
            wires[name].write_text(compile_model(case), encoding="ascii")
        paths = [*list((ROOT / "planner/src").rglob("*.rs")), ROOT / "planner/examples/compare_worker.rs",
                 *list((OPTA / "src").rglob("*.java")), OPTA / "pom.xml",
                 *(HERE / n for n in ("cases.py", "mic_model.py", "run.py"))]
        meta = dict(time_utc=datetime.now(timezone.utc).isoformat(), platform=platform.platform(),
                    cpu=output(["sysctl", "-n", "machdep.cpu.brand_string"]) if sys.platform == "darwin" else platform.processor(),
                    git=output(["git", "rev-parse", "HEAD"]), rust=output(["rustc", "-Vv"]),
                    java=output([cmd["OptaPlanner"][0], "-version"]), jvm_flags=JVM_FLAGS,
                    solver_threads=1, optaplanner="8.44.0.Final", budget_ms=args.budget_ms,
                    seeds=args.seeds, warmup_ms=args.warmup_ms, warmup_seeds=[991, 997],
                    artifact_sha256=artifact_hashes, cases_sha256={k:digest(v) for k,v in inputs.items()},
                    mic_wire_sha256={k:digest(v) for k,v in wires.items()},
                    source_sha256={str(p.relative_to(ROOT)):digest(p) for p in sorted(paths)})
        (dest / "environment.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2)+"\n")
        def paired(stage, round_index, seed, budget):
            # 轮转算例，配对交换先后次序，任何时刻只给一个引擎发送任务。
            for j in range(len(CASES)):
                name = CASES[(j + round_index) % len(CASES)]
                order = list(workers)
                if (CASES.index(name) + round_index) % 2:
                    order.reverse()
                for engine in order:
                    request = wires[name] if engine == "MIC" else inputs[name]
                    result = workers[engine].solve(request, seed, budget)
                    validate(inputs[name], result, budget)
                    result.update(stage=stage, case=name, seed=seed, budget_ms=budget,
                                  engine=engine, position=order.index(engine)+1)
                    raw.write(json.dumps(result, separators=(",", ":"), allow_nan=False)+"\n")
                    raw.flush()
                    print(stage, name, seed, engine, result["first_ms"], result["hard"],
                          result["soft"], "returned_ms", result["elapsed_ms"], flush=True)
        for i, seed in enumerate(meta["warmup_seeds"]):
            paired("warmup", i, seed, args.warmup_ms)
        for i, seed in enumerate(args.seeds):
            paired("measured", i, seed, args.budget_ms)


if __name__ == "__main__":
    main()
