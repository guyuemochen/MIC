#!/usr/bin/env python3
"""用现有探针比较诊断关闭与实际采集时的完整求解结果。"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
FEATURE = "lu_failure_diagnostics"
SOLVE_CALL = "    let result = solve_controlled(problem, &opts, &ctrl, None)?;"


def run(args):
    return subprocess.run(args, cwd=ROOT, capture_output=True, check=True, text=True)


def probe_source(directory):
    source = (ROOT / "milp/examples/round_probe.rs").read_text(encoding="utf-8")
    if source.count(SOLVE_CALL) != 1:
        raise ValueError("探针求解入口已变，请核对采集位置")
    source = source.replace(
        SOLVE_CALL,
        "    #[cfg(feature = \"lu_failure_diagnostics\")]\n"
        "    milp::simplex_sparse::lu_failure_diagnostics::reset();\n"
        f"{SOLVE_CALL}\n"
        "    #[cfg(feature = \"lu_failure_diagnostics\")]\n"
        '    eprintln!("lu_failure_captured={}", '
        "milp::simplex_sparse::lu_failure_diagnostics::take().is_some());",
    )
    path = directory / "lu_parity_probe.rs"
    path.write_text(source, encoding="utf-8")
    return path


def build(source, enabled):
    cargo = [
        "cargo", "build", "--manifest-path", "milp/Cargo.toml", "--release",
        "--lib", "--message-format=json",
    ]
    if enabled:
        cargo.extend(["--features", FEATURE])
    messages = (json.loads(line) for line in run(cargo).stdout.splitlines())
    libraries = [
        Path(name)
        for message in messages
        if message["reason"] == "compiler-artifact" and message["target"]["name"] == "milp"
        for name in message["filenames"]
        if name.endswith(".rlib")
    ]
    library, = libraries
    output = source.with_name("lu_probe_on" if enabled else "lu_probe_off")
    rustc = [
        "rustc", "--edition=2021", "-O", "-L", f"dependency={library.parent / 'deps'}",
        "--extern", f"milp={library}", "-o", str(output), str(source),
    ]
    if enabled:
        rustc.extend(["--cfg", f'feature="{FEATURE}"'])
    run(rustc)
    return output


def solve(probe, model, options, enabled):
    result = run([str(probe), "solve", str(model), *options])
    report = json.loads(result.stdout)
    del report["elapsed_secs"]
    captures = [line for line in result.stderr.splitlines() if line.startswith("lu_failure_captured=")]
    if enabled and captures in (["lu_failure_captured=true"], ["lu_failure_captured=false"]):
        return report, captures[0].endswith("=true")
    if not enabled and not captures:
        return report, None
    raise ValueError(f"{model.name} 采集开关输出异常：{captures}")


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("models", nargs="+", type=Path)
    parser.add_argument("--node-limit", type=int, default=500)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    options = [
        "--threads", "1", "--lp-engine", "revised", "--time-limit-secs", "0",
        "--node-limit", str(args.node_limit), "--early-stop", "off",
    ]
    if args.node_limit <= 0:
        parser.error("--node-limit 必须为正")

    cases = []
    with tempfile.TemporaryDirectory(prefix="mic-lu-parity-") as temp:
        source = probe_source(Path(temp))
        on = build(source, enabled=True)
        off = build(source, enabled=False)
        for model in args.models:
            model = model.resolve(strict=True)
            baseline, _ = solve(off, model, options, enabled=False)
            traced, captured = solve(on, model, options, enabled=True)
            if canonical(baseline) != canonical(traced):
                fields = sorted(baseline.keys() | traced.keys())
                changed = [
                    field for field in fields
                    if field not in baseline or field not in traced
                    or canonical(baseline[field]) != canonical(traced[field])
                ]
                raise ValueError(f"{model.name} 的求解结果不一致：{changed}")
            values = baseline["values"]
            cases.append({
                "model": model.name,
                "input_sha256": hashlib.sha256(model.read_bytes()).hexdigest(),
                **{field: baseline[field] for field in (
                    "status", "objective", "bound", "nodes", "lp_iterations",
                )},
                "solution_length": len(values) if values is not None else None,
                "solution_sha256": hashlib.sha256(canonical(values).encode()).hexdigest(),
                "lu_failure_captured": captured,
            })

    summary = {
        "source": run(["git", "rev-parse", "HEAD"]).stdout.strip(),
        "rustc": run(["rustc", "--version"]).stdout.strip(),
        "feature": FEATURE,
        "options": options,
        "cases": cases,
    }
    args.output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        if error.stdout:
            print(error.stdout, end="")
        if error.stderr:
            print(error.stderr, end="")
        raise
