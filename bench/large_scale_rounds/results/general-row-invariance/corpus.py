#!/usr/bin/env python3
"""在冻结的百例清单上逐例验解，核对已有证明及可行目标，输出短预算的状态和工作量。"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROUNDS = Path(__file__).resolve().parents[2]
FROZEN = ROUNDS / "results/main-gurobi10-100-60s"
sys.path.insert(0, str(ROUNDS))
from compare_all_timing import read_model
from check_results import verify_bound, verify_solution
from summarize_all_timing import audit_truth, read_bytes, read_events


def run(binary, seconds, scratch, instances):
    previous = ROUNDS / "results/remaining-after-weighted-setcover/mic-runs.jsonl.gz"
    with gzip.open(previous, "rt") as source:
        reference = {row["instance"]: row for row in map(json.loads, source)}
    official = json.loads((FROZEN / "verification.json").read_bytes())["extra_official_truth"]
    results = {}
    for case in json.loads((FROZEN / "manifest.json").read_bytes())["cases"]:
        name = case["instance"]
        if instances and name not in instances:
            continue
        model = read_model(FROZEN / case["model"])
        instance = scratch / (name + ".milp")
        content = read_bytes(FROZEN / case["input"])
        assert hashlib.sha256(content).hexdigest() == case["milp_sha256"]
        instance.write_bytes(content)
        command = [binary, "timed-solve", instance, "--time-limit-secs", str(seconds),
                   "--threads", "1", "--node-limit", "2000000", "--early-stop", "off"]
        result = json.loads(subprocess.check_output(command, timeout=seconds + 60))
        errors = []
        final = verify_solution(model, result, errors.append)
        _, callbacks, rejected = read_events(model, result, errors.append)
        witnesses = callbacks + ([final] if final is not None else [])
        prior = reference[name]
        if prior["objective"] is not None:
            verify_bound(model, result, prior["objective"], "已核验解", errors.append)
        truth = official.get(name)
        if prior["status"] in ("optimal", "infeasible", "unbounded"):
            truth = {"status": prior["status"], "objective": prior["objective"]}
        audit_truth(model, truth, result, witnesses, errors.append)
        assert not (errors or rejected or result.get("callback_errors")), (name, errors, rejected)
        results[name] = {k: result[k] for k in ("status", "objective", "bound", "nodes",
                          "lp_iterations", "elapsed_secs", "presolve_reductions")}
        results[name]["checked_callbacks"] = len(callbacks)
    if instances:
        assert set(results) == set(instances)
    return {"probe_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
            "seconds": seconds, "threads": 1, "node_limit": 2000000, "results": results}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--seconds", type=int, default=2)
    parser.add_argument("--instances", nargs="+")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as folder:
        data = run(args.binary.resolve(), args.seconds, Path(folder), args.instances)
    args.output.write_text(json.dumps(data, indent=2) + "\n")
