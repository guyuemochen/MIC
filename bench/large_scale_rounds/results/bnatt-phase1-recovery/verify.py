#!/usr/bin/env python3
"""复核同源定点报告、原模型摘要和首错证据。"""

import argparse
import gzip
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / "bnatt-current-lu"


def load(path):
    return json.loads(path.read_bytes())


def sha(data):
    return hashlib.sha256(data).hexdigest()


def result(name):
    return load(HERE / f"{name}.json")


def report_fields(value):
    return {key: val for key, val in value.items() if key not in {"elapsed_secs", "proof_profile"}}


def verify(live):
    manifest = result("manifest")
    fields = len(report_fields(result("bnatt-default")))
    assert fields == 23
    assert report_fields(result("bnatt-default")) == report_fields(load(PREVIOUS / "baseline.json"))
    for variant in ("default", "candidate"):
        assert report_fields(result(f"bnatt-{variant}")) == report_fields(result(f"bnatt-{variant}-traced"))
    healthy = report_fields(result("assign-default"))
    assert healthy == report_fields(result("assign-candidate"))
    # 四个行最小成本 1、2、3、4 能同时取到，且对应不同列。
    assert healthy["status"] == "optimal" and healthy["objective"] == healthy["bound"] == 10
    assert healthy["values"] == [0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0]

    capture = gzip.decompress((PREVIOUS / "capture.json.gz").read_bytes())
    model_hash = manifest["canonical_original_model_sha256"]
    model_bytes = (HERE.parent / "models" / f"{model_hash}.json").read_bytes()
    assert sha(model_bytes) == model_hash
    model = json.loads(model_bytes)
    assert (model["sense"], model["offset"], model["objective"]) == ("min", 0, [[0, 1]])
    assert model["mps_columns"][0] == "x1" and model["variables"][0]["lb"] == 0
    assert result("bnatt-candidate")["bound"] == 0

    for variant, status, iterations, bound, stop in (
        ("default", "iteration_limit", 3904, None, ("root_lp", "iteration_limit")),
        ("candidate", "optimal", 7211, 0, ("root_heuristic_restore", "time_limit")),
    ):
        events = result(f"bnatt-{variant}-traced")["proof_profile"]["events"]
        root = next(event for event in events if event["kind"] == "root_lp")
        assert (root["status"], root["iterations"], root["bound"]) == (status, iterations, bound)
        last_failure = [event for event in events if event["kind"] == "lp_failure"][-1]
        assert (last_failure["stage"], last_failure["status"]) == stop
        assert manifest[f"live_{variant}_capture_raw_sha256"] == sha(capture)
        if live is not None:
            assert (live / f"bnatt-{variant}-capture.json").read_bytes() == capture
            export = load(live / f"{variant}-original.json")
            export["mps_columns"] = model["mps_columns"]
            encoded = json.dumps(export, sort_keys=True, ensure_ascii=False, allow_nan=False,
                                 separators=(",", ":")) + "\n"
            assert encoded.encode() == model_bytes

    print(json.dumps({
        "equal_nontime_fields_per_comparison": fields,
        "first_failure_raw_sha256": sha(capture),
        "original_model_sha256": model_hash,
        "original_bound_checked": 0,
        "verified_live_bytes": live is not None,
    }, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--live", type=Path, help="也核对重采原始文件与重新导出模型")
    verify(parser.parse_args().live)
