#!/usr/bin/env python3
"""bench/run_timings.py —— 计时矩阵驱动（第 4 轮，2026-10-09）。

用法：
  python bench/run_timings.py --set gate --threads 1 --json bench/gate_runs.json
  python bench/run_timings.py --set full --json bench/full_runs.json

口径（测量工具 milp/examples/profile_stats.rs，复刻 CLI 默认选项）：
- gate 集（快速回归，单次总耗时控制在约 3 分钟内）：
  knapsack-200 / setcover-100x150 / knapsack-2000 + 三个大算例
  （lp-sparse-4000x800 / setcover-400x800 / mknap-5x300）；
  --time-limit-secs 30；threads 只跑 --threads 指定值（gate 必填）；
  重复次数：快算例 3 次、三个大算例各 1 次——实测总求解 113.9s
  （2026-10-09，t1）。注意：时限只在 B&B 节点边界检查
  （branch_bound.rs:711，simplex 内部无时限），单节点纯 LP 的根 LP
  会跑完——lp-sparse 在 30s 时限下仍以 ~34s optimal 终态结束，
  大算例最坏 ≈ 36s+17s+30s+小算例 <5s ≈ 1.6min < 3min；
- full 集：原 13 算例（第 3 轮顺序不变）+ 第 4 轮新 5 个，
  threads ∈ {1,4} 都跑（--threads 给出时提示并忽略）；--time-limit-secs
  90；快算例 3 次、SLOW（≥20s 量级，按 2026-10-09 t1 标定）2 次；
- --json PATH：写中位数行数组，每行字段 instance/threads/wall_secs/
  status/objective/nodes；wall_secs 取重复中位；status 取众数；
  objective/nodes 取墙钟最接近中位数的那次（与 analyze_baseline.py
  同口径；limit_reached 无 incumbent 时 objective 如实为 null）；
- 逐次原始记录：gate → bench/raw_runs_gate.jsonl、full →
  bench/raw_runs_full.jsonl（各自独立文件互不覆盖；每次调用覆盖本集
  自己的文件；第 3 轮的 bench/raw_runs.jsonl 保留不动，供
  analyze_baseline.py 继续消费）；
- wall_secs 为 solve_control 周期内计时（不含算例加载，与既往口径一致）；
- 入口先 `cargo build --release --example profile_stats`（失败即退出）：
  杜绝源码-二进制漂移——陈旧二进制会让一切 bench 数字失去合法来源。
"""
import argparse
import json
import statistics
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent
EXE = ROOT.parent / "milp" / "target" / "release" / "examples" / "profile_stats.exe"

GATE_ORDER = [
    "knapsack-200",
    "setcover-100x150",
    "knapsack-2000",
    "lp-sparse-4000x800",
    "setcover-400x800",
    "mknap-5x300",
]
GATE_TIME_LIMIT_SECS = 30
# gate 三个大算例只跑 1 次：实测（2026-10-09 t1）lp-sparse ~34s optimal
# （时限只在节点边界检查，单节点 LP 跑完可越过 30s 名义值）、
# setcover-400x800 ~13-17s limit（迭代上限提前终止）、mknap-5x300 ~15s
# auto_degraded（30s 预算下自适应降级提前停）；合计 + 小算例 ×3 ≈ 2min < 3min。
GATE_SINGLE_REP = {"lp-sparse-4000x800", "setcover-400x800", "mknap-5x300"}

FULL_ORDER = [
    # —— 原 13 算例（第 3 轮基准，顺序不变）——
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
    # —— 第 4 轮新增 5（2026-10-09 t1+90s 标定：lp-sparse 35.8s optimal /
    #    setcover-400x800 17.3s limit（迭代上限 24000）/ knapsack-5000
    #    2.2s optimal（DP 真值 88526 吻合）/ mknap-10x200 26.4s optimal /
    #    assign-40x40 0.12s optimal）——
    "lp-sparse-4000x800",
    "setcover-400x800",
    "knapsack-5000",
    "mknap-10x200",
    "assign-40x40",
]
FULL_TIME_LIMIT_SECS = 90
FULL_SLOW = {
    "knapsack-500",
    "mknap-5x300",
    "marketsplit-4x30",
    "lp-sparse-4000x800",  # t1 实测 ~36s（>20s 量级）
    "setcover-400x800",    # t1 实测 ~17s（贴近 20s，保守入慢集）
    "mknap-10x200",        # t1 实测 ~26s
}

RAW_BY_SET = {"gate": ROOT / "raw_runs_gate.jsonl", "full": ROOT / "raw_runs_full.jsonl"}
JSON_DEFAULT = {"gate": ROOT / "gate_runs.json", "full": ROOT / "full_runs.json"}


def run_one(name: str, threads: int, time_limit_secs: int) -> dict:
    path = ROOT / "instances" / f"{name}.milp"
    cmd = [str(EXE), str(path), "--threads", str(threads),
           "--time-limit-secs", str(time_limit_secs)]
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


def median_row(runs: list) -> dict:
    """中位数行：wall_secs 取重复中位；status 取众数；objective/nodes 取
    墙钟最接近中位数的那次（并列取低者）——与 analyze_baseline.py 同口径。
    limit_reached 无 incumbent ⇒ objective 为 null，如实保留。"""
    med_wall = statistics.median(r["wall_secs"] for r in runs)
    med_run = min(runs, key=lambda r: (abs(r["wall_secs"] - med_wall), r["wall_secs"]))
    status = Counter(r["status"] for r in runs).most_common(1)[0][0]
    return {
        "instance": runs[0]["instance"],
        "threads": runs[0]["threads"],
        "wall_secs": round(med_wall, 6),
        "status": status,
        "objective": med_run["objective"],
        "nodes": med_run["nodes"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="milp 计时矩阵（gate=快速回归 / full=全量）")
    parser.add_argument("--set", required=True, choices=("gate", "full"),
                        help="gate：6 算例×指定 threads×30s（约 3 分钟内）；full：18 算例×{1,4}×90s")
    parser.add_argument("--threads", type=int, default=None,
                        help="gate 必填：只跑该线程数；full 固定跑 {1,4}，给出则忽略")
    parser.add_argument("--json", type=Path, default=None,
                        help="中位数行数组输出路径（默认 gate_runs.json / full_runs.json）")
    args = parser.parse_args()

    if args.set == "gate":
        if args.threads is None:
            sys.exit("gate 集必须指定 --threads N（threads 只跑该值）")
        order, threads_list = GATE_ORDER, [args.threads]
        time_limit, single_rep, slow = GATE_TIME_LIMIT_SECS, GATE_SINGLE_REP, set()
    else:
        if args.threads is not None:
            print(f"note: full 集固定跑 threads {{1, 4}}，忽略 --threads {args.threads}")
        order, threads_list = FULL_ORDER, [1, 4]
        time_limit, single_rep, slow = FULL_TIME_LIMIT_SECS, set(), FULL_SLOW

    ensure_probe_built()
    raw_path = RAW_BY_SET[args.set]
    json_path = args.json if args.json is not None else JSON_DEFAULT[args.set]

    with raw_path.open("w", encoding="utf-8"):
        pass  # 覆盖本集旧文件（gate/full 独立文件，互不影响）
    med_rows = []
    t_start = time.perf_counter()
    with raw_path.open("a", encoding="utf-8") as f:
        for name in order:
            for threads in threads_list:
                reps = 1 if name in single_rep else (2 if name in slow else 3)
                runs = []
                for rep in range(reps):
                    row = run_one(name, threads, time_limit)
                    runs.append(row)
                    f.write(json.dumps(row, ensure_ascii=False) + "\n")
                    f.flush()
                    print(f"{name:22s} t={threads} rep={rep + 1}/{reps} "
                          f"wall={row['wall_secs']:9.3f}s nodes={row['nodes']:>9d} "
                          f"status={row['status']}", flush=True)
                med_rows.append(median_row(runs))
    total = time.perf_counter() - t_start
    json_path.write_text(
        json.dumps(med_rows, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"set={args.set} 总求解耗时 {total:.1f}s；wrote {raw_path}；"
          f"{len(med_rows)} 行中位数 -> {json_path}")


if __name__ == "__main__":
    main()
