#!/usr/bin/env python3
"""bench/run_timings.py —— 全算例计时矩阵驱动（第 3 轮基线，2026-10-09）。

口径（与 profile_stats 工具一致，复刻 CLI 默认选项）：
- 每算例 × threads ∈ {1,4}；快算例 3 次重复、慢算例（SLOW 集，>20s 量级）2 次；
- 统一 --time-limit-secs 90（本轮协议：终态须 ≤90s 给出，最优或 gap 类终止
  皆可；limit_reached 记录为越界/缺陷证据）；
- 每次运行一行 JSON 追加写入 bench/raw_runs.jsonl（覆盖旧文件）；
- wall_secs 为 solve_control 周期内计时（不含算例加载，与上一轮口径一致）；
- 入口先 `cargo build --release --example profile_stats`（失败即退出）：
  杜绝源码-二进制漂移——旧版仅检查 exe 存在，陈旧二进制会让一切 bench
  数字失去合法来源。
"""
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).parent
EXE = ROOT.parent / "milp" / "target" / "release" / "examples" / "profile_stats.exe"
OUT = ROOT / "raw_runs.jsonl"
TIME_LIMIT_SECS = 90

ORDER = [
    "assign-16x16",
    "lp-wide-200x300",
    "milp-mixed-120i30c",
    "milp-mixed-250i40c",
    "knapsack-200",
    "knapsack-500",
    "knapsack-2000",
    "mknap-3x100",
    "mknap-5x300",
    "marketsplit-2x20",
    "marketsplit-4x30",
    "setcover-60x100",
    "setcover-100x150",
]
SLOW = {"knapsack-500", "mknap-5x300", "marketsplit-4x30"}


def run_one(name: str, threads: int) -> dict:
    path = ROOT / "instances" / f"{name}.milp"
    cmd = [str(EXE), str(path), "--threads", str(threads),
           "--time-limit-secs", str(TIME_LIMIT_SECS)]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    outer = time.perf_counter() - t0
    row = json.loads(proc.stdout.strip().splitlines()[-1])
    row["outer_secs"] = round(outer, 6)
    return row


def ensure_probe_built() -> None:
    """重建探针二进制（杜绝源码-二进制漂移），失败即退出。"""
    cmd = ["cargo", "build", "--release", "--example", "profile_stats"]
    print(f"$ {' '.join(cmd)} (cwd={ROOT.parent / 'milp'})", flush=True)
    proc = subprocess.run(cmd, cwd=str(ROOT.parent / "milp"))
    if proc.returncode != 0:
        sys.exit(f"cargo build 失败（exit {proc.returncode}）：探针源码编译不过，中止计时")
    if not EXE.exists():
        sys.exit(f"cargo build 成功但 {EXE} 不存在（产物路径异常），中止计时")


def main() -> None:
    ensure_probe_built()
    with OUT.open("w", encoding="utf-8") as f:
        pass  # 覆盖旧文件
    with OUT.open("a", encoding="utf-8") as f:
        for name in ORDER:
            for threads in (1, 4):
                reps = 2 if name in SLOW else 3
                for rep in range(reps):
                    row = run_one(name, threads)
                    f.write(json.dumps(row, ensure_ascii=False) + "\n")
                    f.flush()
                    print(f"{name:22s} t={threads} rep={rep + 1}/{reps} "
                          f"wall={row['wall_secs']:9.3f}s nodes={row['nodes']:>9d} "
                          f"status={row['status']}", flush=True)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
