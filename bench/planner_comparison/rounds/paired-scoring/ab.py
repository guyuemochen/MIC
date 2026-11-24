"""串行轮换两版 MIC 与可选的 OptaPlanner，独立复算冻结输入上的首解和终解。"""
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import sys
from tempfile import TemporaryDirectory

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH))
from cases import CASES
from mic_model import compile_model
from run import JVM_FLAGS, OPTA, ROOT, SEEDS, Worker, build, digest, output, validate


BASELINE = BENCH / "results/2026-11-24/environment.json"
VERSIONS = ("before", "after")
OPTA_NAME = "OptaPlanner"


def schedule(meta):
    """排课轮换出场顺序，两版 MIC 在相邻种子中互换先后，三个短测单列。"""
    phases = [("warmup", "timetable", seed, meta["warmup_ms"])
              for seed in meta["warmup_seeds"]]
    phases.extend(("measured", "timetable", seed, meta["budget_ms"])
                  for seed in meta["seeds"])
    phases.extend(("regression", case, meta["regression_seed"], meta["regression_ms"])
                  for case in CASES if case != "timetable")
    for i, (stage, case, seed, budget) in enumerate(phases):
        if meta.get("with_opta") and stage != "regression":
            # 每组三次轮流占据第一、二、三位，再反向轮换。
            order = (*VERSIONS, OPTA_NAME)
            j = i % 3
            if i // 3 % 2:
                order = order[::-1]
                j = -j
            order = order[j:] + order[:j]
        else:
            order = VERSIONS if i % 2 == 0 else reversed(VERSIONS)
        for position, engine in enumerate(order, 1):
            yield dict(stage=stage, case=case, seed=seed, budget_ms=budget,
                       position=position, engine=engine)


def case_path(name):
    return BENCH / "instances" / (name + ".txt")


def prepare(wire_dir, baseline, frozen=None):
    for name in ("cases.py", "mic_model.py", "run.py"):
        script = BENCH / name
        if digest(script) != baseline["source_sha256"][str(script.relative_to(ROOT))]:
            raise ValueError("算例编译器或验分程序与上轮归档不符：" + name)
    for case in CASES:
        path = case_path(case)
        wire = wire_dir / (case + ".wire")
        wire.write_bytes(compile_model(path).encode("ascii"))
        if (digest(path) != baseline["cases_sha256"][case] or
                digest(wire) != baseline["mic_wire_sha256"][case]):
            raise ValueError(f"{case} 的算例或模型与上轮归档不符")
        if frozen is not None and (frozen / wire.name).read_bytes() != wire.read_bytes():
            raise ValueError(f"{case} 的独立冻结副本与重建模型不同")


def run(args):
    baseline = json.loads(BASELINE.read_text())
    binary = {v: getattr(args, v).resolve(strict=True) for v in VERSIONS}
    commands = {v: [str(path)] for v, path in binary.items()}
    opta_settings = {}
    java_sources = []
    if args.with_opta:
        java_sources = [*sorted((OPTA / "src").rglob("*.java")), OPTA / "pom.xml"]
        for path in java_sources:
            if digest(path) != baseline["source_sha256"][str(path.relative_to(ROOT))]:
                raise ValueError("OptaPlanner 对照实现与上轮归档不符：" + path.name)
        runner, artifacts = build()
        if artifacts["optaplanner_jar"] != baseline["artifact_sha256"]["optaplanner_jar"]:
            raise ValueError("OptaPlanner 核心 JAR 与上轮归档不符")
        commands[OPTA_NAME] = runner[OPTA_NAME]
        opta_settings = dict(optaplanner=baseline["optaplanner"], jvm_flags=JVM_FLAGS,
                             java=output([runner[OPTA_NAME][0], "-version"]),
                             optaplanner_jar_sha256=artifacts["optaplanner_jar"])
    hashes = {v: digest(path) for v, path in binary.items()}
    expected = args.before_sha256 or baseline["artifact_sha256"]["mic_binary"]
    if hashes["before"] != expected:
        raise ValueError("旧版程序指纹与预期不符；在其他工具链重建旧版后可指定 --before-sha256")
    dest = args.out.resolve()
    with TemporaryDirectory(prefix="mic-pair-ab-") as tmp:
        wire_dir = Path(tmp)
        prepare(wire_dir, baseline, args.frozen_wires)
        dest.mkdir(parents=True, exist_ok=False)
        source = [*sorted((ROOT / "planner/src").rglob("*.rs")),
                  ROOT / "planner/examples/compare_worker.rs", Path(__file__).resolve(),
                  *(BENCH / n for n in ("cases.py", "mic_model.py", "run.py")), *java_sources]
        meta = dict(time_utc=datetime.now(timezone.utc).isoformat(), platform=platform.platform(),
                    cpu=output(["sysctl", "-n", "machdep.cpu.brand_string"])
                    if sys.platform == "darwin" else platform.processor(),
                    rust=output(["rustc", "-Vv"]), git=output(["git", "rev-parse", "HEAD"]),
                    seeds=args.seeds, budget_ms=args.budget_ms,
                    warmup_seeds=[991, 997], warmup_ms=args.warmup_ms,
                    regression_ms=args.regression_ms, regression_seed=7,
                    artifact_sha256=hashes, expected_before_sha256=expected,
                    cases_sha256=baseline["cases_sha256"], mic_wire_sha256=baseline["mic_wire_sha256"],
                    baseline_manifest_sha256=digest(BASELINE),
                    workspace_source_sha256={str(p.relative_to(ROOT)): digest(p) for p in source},
                    solver_threads=1, construction="BestFit", algorithm="TabuSearch", tabu_size=8,
                    with_opta=args.with_opta, **opta_settings)
        (dest / "environment.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n")
        with ExitStack() as stack:
            raw = stack.enter_context((dest / "runs.jsonl").open("x", encoding="utf8"))
            workers = {}
            for version, cmd in commands.items():
                log = stack.enter_context((dest / (version + ".stderr.log")).open("x", encoding="utf8"))
                workers[version] = Worker(cmd, log)
                stack.callback(workers[version].close)
            for request in schedule(meta):
                case, budget = request["case"], request["budget_ms"]
                path = case_path(case) if request["engine"] == OPTA_NAME else wire_dir / (case + ".wire")
                result = workers[request["engine"]].solve(path, request["seed"], budget)
                validate(case_path(case), result, budget)
                if request["stage"] == "regression" and not result["checked_feasible"]:
                    raise ValueError(f"{case} 短测未取得完整可行解：{request}")
                result.update(request)
                raw.write(json.dumps(result, allow_nan=False, separators=(",", ":")) + "\n")
                raw.flush()
                print(request["stage"], case, request["seed"], request["engine"],
                      result["first_ms"], result["hard"], result["soft"], flush=True)
    verify(dest)


def verify(dest):
    meta = json.loads((dest / "environment.json").read_text())
    baseline = json.loads(BASELINE.read_text())
    if (meta["baseline_manifest_sha256"] != digest(BASELINE) or
            meta["artifact_sha256"]["before"] != meta["expected_before_sha256"]):
        raise ValueError("基线清单或旧版程序记录不匹配")
    for key in ("cases_sha256", "mic_wire_sha256"):
        if meta[key] != baseline[key]:
            raise ValueError("归档指纹不匹配：" + key)
    if meta.get("with_opta") and (
            meta["optaplanner_jar_sha256"] != baseline["artifact_sha256"]["optaplanner_jar"]
            or meta["jvm_flags"] != baseline["jvm_flags"]
            or meta["optaplanner"] != baseline["optaplanner"]):
        raise ValueError("OptaPlanner 归档与上轮依赖或 JVM 设置不符")
    with TemporaryDirectory(prefix="mic-pair-check-") as tmp:
        prepare(Path(tmp), baseline)
    records = [json.loads(line) for line in (dest / "runs.jsonl").read_text().splitlines()]
    expected = list(schedule(meta))
    if len(records) != len(expected):
        raise ValueError("归档缺少配对记录，或含额外记录")
    for row, req in zip(records, expected):
        if any(row[k] != v for k, v in req.items()):
            raise ValueError("配对次序、引擎、种子或预算不符")
        checked = dict(row)
        validate(case_path(req["case"]), checked, req["budget_ms"])
        if checked != row or (req["stage"] == "regression" and not checked["checked_feasible"]):
            raise ValueError("已存验分与独立复算不符")
    print(f"已核对 {len(records)} 条试验记录。")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    trial = sub.add_parser("run", help="运行两版并归档")
    trial.add_argument("--before", required=True, type=Path)
    trial.add_argument("--after", required=True, type=Path)
    trial.add_argument("--with-opta", action="store_true", help="在排课预热和计时轮中加入 OptaPlanner")
    trial.add_argument("--before-sha256", help="重建旧版时可指定已验证的二进制指纹")
    trial.add_argument("--frozen-wires", type=Path, help="可选的独立冻结模型副本所在目录")
    trial.add_argument("--out", required=True, type=Path, help="尚不存在的归档目录")
    trial.add_argument("--seeds", type=int, nargs="+", default=SEEDS)
    trial.add_argument("--budget-ms", type=int, default=10000)
    trial.add_argument("--warmup-ms", type=int, default=2000)
    trial.add_argument("--regression-ms", type=int, default=500)
    sub.add_parser("verify", help="独立复查已有归档").add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "verify":
        verify(args.out)
    else:
        if min(args.budget_ms, args.warmup_ms, args.regression_ms) <= 0:
            parser.error("求解预算须为正")
        limit = 2**63 if args.with_opta else 2**64
        if len(set(args.seeds)) != len(args.seeds) or any(s < 0 or s >= limit for s in args.seeds):
            parser.error("种子超出求解器支持的非负整数范围，或存在重复")
        run(args)


if __name__ == "__main__":
    main()
