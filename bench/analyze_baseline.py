#!/usr/bin/env python3
"""汇总 bench/raw_runs.jsonl（profile_stats 每行一次运行）为 bench/baseline.json。

口径（自述）：
- 第 3 轮（2026-10-09）：13 算例（旧 7 + 新 6）× threads ∈ {1,4}，快算例 3 次、
  慢算例（knapsack-500/mknap-5x300/marketsplit-4x30）2 次重复；
- 统一 --time-limit-secs 90（终态协议：≤90s 内给出最优或 gap 类终止）；
- wall_secs 取中位数；nodes / lp_iterations / cuts / objective 取中位样本；
- derived: node_throughput = nodes / wall_secs；iter_density = lp_iterations / nodes；
  per_node_us = wall_secs/nodes × 1e6；iter_per_sec = lp_iterations / wall_secs。
- status: 众数（t1 亦可能因自适应控制器含墙钟判定而跨次摆动，见 experiments；
  并行路径工作窃取顺序影响 incumbent 发现时机）。
"""
import json
import statistics
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent
RAW = ROOT / "raw_runs.jsonl"
OUT = ROOT / "baseline.json"

rows = [json.loads(l) for l in RAW.read_text(encoding="utf-8").splitlines() if l.strip()]

CLASS = {
    "knapsack-200": "knapsack(强相关单约束背包)",
    "knapsack-500": "knapsack(强相关单约束背包)",
    "knapsack-2000": "knapsack(强相关单约束背包·大规模)",
    "mknap-3x100": "mknap(多维背包)",
    "mknap-5x300": "mknap(多维背包·大规模)",
    "assign-16x16": "assign(指派)",
    "milp-mixed-120i30c": "mixed(一般MILP)",
    "milp-mixed-250i40c": "mixed(一般MILP)",
    "lp-wide-200x300": "lp-wide(纯LP)",
    "marketsplit-2x20": "marketsplit(等式市场分割)",
    "marketsplit-4x30": "marketsplit(等式市场分割)",
    "setcover-60x100": "setcover(集合覆盖)",
    "setcover-100x150": "setcover(集合覆盖)",
}

by_key = {}
for r in rows:
    by_key.setdefault((r["instance"], r["threads"]), []).append(r)

results = []
for (inst, threads), runs in sorted(by_key.items()):
    walls = sorted(r["wall_secs"] for r in runs)
    med_wall = statistics.median(walls)
    # 中位样本：wall 最接近中位数的那次（秩一致时取低者）
    med_run = min(runs, key=lambda r: (abs(r["wall_secs"] - med_wall), r["wall_secs"]))
    status_modes = Counter(r["status"] for r in runs)
    nodes = med_run["nodes"]
    obj = med_run["objective"]
    lpi = med_run["lp_iterations"]
    results.append({
        "instance": inst,
        "class": CLASS[inst],
        "threads": threads,
        "reps": len(runs),
        "wall_secs_all": walls,
        "wall_secs": round(med_wall, 6),
        "nodes": nodes,
        "objective": obj,
        "status": status_modes.most_common(1)[0][0],
        "status_all": dict(status_modes),
        "lp_iterations": lpi,
        "cuts": med_run["cuts"],
        "stats_time_ms": med_run["stats_time_ms"],
        "bound": med_run["bound"],
        "gap_percent_at_stop": med_run["gap_percent"],
        "node_throughput_per_sec": round(nodes / med_wall, 1) if med_wall > 0 and nodes else None,
        "iter_density_per_node": round(lpi / nodes, 3) if nodes else None,
        "per_node_us": round(med_wall / nodes * 1e6, 2) if nodes else None,
        "lp_iter_per_sec": round(lpi / med_wall, 0) if med_wall > 0 else None,
    })

# 背包类 DP 核对（k2000 亦有 .opt 真值）
dp_check = {}
for name in ("knapsack-200", "knapsack-500", "knapsack-2000"):
    opt = int((ROOT / "instances" / f"{name}.opt").read_text().strip())
    r1 = next(r for r in results if r["instance"] == name and r["threads"] == 1)
    obj = int(r1["objective"])
    dp_check[name] = {
        "dp_optimal": opt,
        "solver_objective_t1": obj,
        "matches": obj == opt,
        "proved_optimal": r1["status"] == "optimal",
    }

doc = {
    "_format": (
        "bench/baseline.json —— 第 3 轮基准画像（2026-10-09，Windows/14 逻辑核）。"
        "13 算例 = 旧 7 + 新 6（marketsplit-2x20/4x30、setcover-60x100/100x150、"
        "knapsack-2000、mknap-5x300，规模经 90s/threads=1 终态标定，种子确定）。"
        "timings[]: 每算例×threads{1,4}，快 3 次/慢 2 次重复，wall_secs 取中位；"
        "统一 time-limit 90s（本轮协议变更；旧基线的 knapsack-500 无上限跑满 "
        "200 万节点 125s——当前工作区代码已能在 776 节点证最优，旧数字全面过时）。"
        "测量工具：milp/examples/profile_stats.rs（临时，画像后删除）——库 API "
        "solve_controlled 复刻 CLI 默认选项（MostFractional/DepthFirst/"
        "node_limit=2M/root_cuts/heuristics=auto/dive 4096/presolve on）+ 90s 上限，"
        "输出 stats.lp_iterations；CLI 抽查逐字段一致（marketsplit-4x30 t1: "
        "395095 节点/obj 15；setcover-60x100 t1: 1 节点/obj 243）。"
        "原始逐次数据：bench/raw_runs.jsonl。experiments[]: 隔离实验（非计时基线）。"
    ),
    "machine": {"logical_cores": 14, "os": "win32", "toolchain": "cargo release (default profile)"},
    "timings": results,
    "knapsack_dp_check": dp_check,
    "experiments": {
        "branch_rule_pseudocost_isolation_t1": {
            "desc": "BranchRule::Pseudocost vs 默认 MostFractional（profile_stats "
                    "--branch-rule；CLI 默认硬编码 MostFractional 于 main.rs:298）。",
            "marketsplit-4x30": {
                "mostfractional_median": {"wall_secs": 8.281, "nodes": 395095, "lp_iterations": 3348919},
                "pseudocost": {"wall_secs": 4.587, "nodes": 299280, "lp_iterations": 2849031,
                               "status": "optimal", "objective": 15},
                "conclusion": "节点 −24%（395k→299k）、墙钟 −45%（8.3s→4.6s）、迭代 −15%。",
            },
            "mknap-5x300": {
                "mostfractional_90s": {"wall_secs": 82.0, "nodes": 336610,
                                        "status": "limit_reached/auto_degraded 摆动",
                                        "objective": 4752, "bound": 4773.624886, "gap_percent": 0.455},
                "pseudocost_90s": {"wall_secs": 74.577, "nodes": 371665, "status": "optimal",
                                    "objective": 4762, "bound": 4762},
                "conclusion": "Pseudocost 在同预算内证得最优 4762（MostFractional 90s 未达该 "
                              "incumbent 且 gap 0.46% 不闭合）——分支规则是当前树规模的第一杠杆。",
            },
        },
        "mknap-5x300_adaptive_dive_30s_t1": {
            "desc": "30s 预算、threads=1：adaptive=0（仅根启发式） vs adaptive=1（默认 auto）。",
            "adaptive0_root_only": {"wall_secs": 30.000492, "nodes": 156110, "lp_iterations": 1728380,
                                    "node_rate_per_sec": 5203.7, "iter_per_node": 11.07,
                                    "status": "limit_reached", "objective": 4737, "bound": 4773.624886},
            "adaptive1_auto": {"wall_secs": 19.841771, "nodes": 100422, "lp_iterations": 1093843,
                                "node_rate_per_sec": 5061.1, "iter_per_node": 10.89,
                                "status": "auto_degraded", "objective": 4732, "bound": 4773.624886,
                                "gap_percent": 0.879647},
            "conclusion": "节点吞吐持平（5.2k vs 5.1k/s），但 adaptive 的 dive 反而拿到更差 "
                          "incumbent（4732 < 4737）并提前降级——多维背包类上 dive 投入接近零收益；"
                          "所有运行的 bound 冻结在根割值 4773.62 ⇒ 缺口在树/割不在启发式。",
        },
        "setcover_root_iterationlimit_abort": {
            "desc": "setcover-100x150 t1：默认（root_cuts=1）vs --root-cuts 0。",
            "default_root_cuts": {"wall_secs": 0.525, "nodes": 1, "lp_iterations": 13288, "cuts": 37,
                                    "status": "limit_reached", "objective": None},
            "root_cuts_off": {"wall_secs": 0.133305, "nodes": 45, "lp_iterations": 3406, "cuts": 0,
                               "status": "optimal", "objective": 276, "bound": 276},
            "conclusion": "根 LP 本体可解（3406 迭代）；是割循环加行后的 resolve 打满单次调用"
                          "迭代上限 max(10k, 20×(m+n))（simplex.rs:1312）触发 IterationLimit，"
                          "cut_loop 将其作为终态上报（cuts.rs:563），branch_bound.rs:364 即整个"
                          "求解以 limit_reached 终结——1 个节点、无 incumbent、无 bound。"
                          "割关闭后 45 节点即最优 ⇒ 该类被『迭代上限传播』整体锁死，而非难度高。",
        },
        "mknap-5x300_t1_nondeterminism": {
            "desc": "同参数 t1 两次：rep1 90s/336610 节点 limit_reached，rep2 74s/392468 节点 "
                    "auto_degraded；另一次探针 76s/356460。原因：自适应控制器与自动降级的判定"
                    "含墙钟比较（search_ctrl.rs:255 elapsed>frac×budget、:508/:528 plateau 窗），"
                    "给了 time budget 后 t1 也随时序抖动改变 dive/降级时机 ⇒ 节点数与终态跨次"
                    "漂移。bench 以 wall_secs_all/status_all 记录离散。",
        },
        "marketsplit_size_ladder_t1": {
            "desc": "90s 上限探针（planted 子集保证可行；系数 [100,200]）。",
            "results": {
                "marketsplit-2x20": "0.012s optimal 2936 节点",
                "marketsplit-3x30": "1.44s optimal 122234 节点",
                "marketsplit-3x40": "1.74s optimal 126768 节点",
                "marketsplit-4x30": "5.5~9.3s optimal 395095 节点（入列）",
                "marketsplit-4x40": "49.1s 打满 200 万节点上限，零可行解（超规格）",
                "marketsplit-5x50": "90s 时间上限，1.56M 节点，零可行解（超规格）",
            },
            "conclusion": "m≥4 且 n≥40 时 B&B 在 200 万节点内找不到任何可行解（等式行下 "
                          "round/dive 均无能为力，见 heuristics.rs:29/:103）——可行解发现"
                          "（feasibility pump 类启发式）是该类的第一瓶颈。",
        },
    },
}
OUT.write_text(json.dumps(doc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"wrote {OUT}")
for r in results:
    print(f"{r['instance']:22s} t={r['threads']} wall={r['wall_secs']:9.3f}s nodes={r['nodes']:9d} "
          f"lpi={r['lp_iterations']:9d} thr={r['node_throughput_per_sec']} it/node={r['iter_density_per_node']} "
          f"status={r['status']}")
