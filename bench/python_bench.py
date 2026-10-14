#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bench/python_bench.py —— micpython 与开源求解器（CBC / HiGHS / SCIP）同机对比。

环境：conda env `mic-bench`（python 3.11 + numpy scipy pulp ortools）。
创建（一次性）：conda create -n mic-bench python=3.11 -y &&
                conda run -n mic-bench pip install numpy scipy pulp ortools
复现：conda run --no-capture-output -n mic-bench python -u bench/python_bench.py

架构：驱动进程把每个（算例 × 求解器）对作为子进程执行（--pair 模式），
硬看门狗 = 时限 + 120s 宽限——单个求解器挂死（实测 CBC 在 marketsplit 类
等式模型上会在不受 -sec 约束的阶段卡死）只损失自己的宽限窗，记
watchdog_kill，不拖垮全表。

口径：
- 算例来自 bench/instances/*.milp（统一解析后接入各 API）；
- 全部求解器 threads=1、统一时限；build/solve 分计（perf_counter）；
- 目标值交叉核对（DP .opt 真值一并核对），不一致显式 MISMATCH；
- 结果写 bench/python_bench_results.{json,md}。
"""
import json
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INST = ROOT / "bench" / "instances"
GRACE_SECS = 120.0  # 看门狗宽限：超过时限这么多仍不返回 ⇒ 强杀

INSTANCES = [
    ("knapsack-200", 60), ("knapsack-500", 60), ("knapsack-2000", 60),
    ("knapsack-5000", 60), ("mknap-3x100", 60), ("mknap-5x300", 120),
    ("setcover-100x150", 60), ("setcover-300x600", 120),
    ("milp-mixed-250i40c", 60), ("marketsplit-4x30", 60),
    ("marketsplit-4x40", 60), ("assign-40x40", 60),
    ("lp-wide-200x300", 60), ("lp-sparse-4000x800", 60),
]

VAR_TYPES = {"binary": "B", "integer": "I", "continuous": "C"}


class Problem:
    def __init__(self):
        self.name = ""
        self.sense = "max"
        self.vars = []
        self.obj = []
        self.cons = []


def parse_milp(path):
    p = Problem()
    section = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        low = line.lower()
        if low == "vars":
            section = "vars"; continue
        if low == "obj":
            section = "obj"; continue
        if low == "cons":
            section = "cons"; continue
        if low == "end":
            break
        if low.startswith("name "):
            p.name = line[5:].strip(); continue
        if low.startswith("sense "):
            p.sense = line[6:].strip().lower(); continue
        if section == "vars":
            toks = line.split()
            if toks[1].lower() not in VAR_TYPES:
                raise ValueError(f"未知变量类型 {toks[1]} @ {path.name}")
            p.vars.append((toks[0], VAR_TYPES[toks[1].lower()]))
        elif section == "obj":
            p.obj.extend(_terms(line))
        elif section == "cons":
            cname, _, rest = line.partition(":")
            m = re.match(r"^(.*?)\s*(<=|>=|=)\s*(-?[\d.eE+]+)\s*$", rest.strip())
            if not m:
                raise ValueError(f"约束无法解析 @ {path.name}: {line[:60]}")
            p.cons.append((cname.strip(), _terms(m.group(1)), m.group(2), float(m.group(3))))
    return p


def _terms(s):
    toks = s.split()
    if len(toks) % 2 != 0:
        raise ValueError(f"项数不是偶数: {s[:60]}")
    return [(float(toks[i]), toks[i + 1]) for i in range(0, len(toks), 2)]


# ---------------------------------------------------------------- micpython --
def run_micpython(p, limit):
    sys.path.insert(0, str(ROOT))
    import micpython as mp

    t0 = time.perf_counter()
    m = mp.Model(p.name, sense=mp.MAXIMIZE if p.sense == "max" else mp.MINIMIZE)
    vs = {nm: m.addVar(vtype=vt, name=nm) for nm, vt in p.vars}
    m.setObjective(mp.quicksum(c * vs[n] for c, n in p.obj),
                   mp.MAXIMIZE if p.sense == "max" else mp.MINIMIZE)
    for cn, terms, rel, rhs in p.cons:
        e = mp.quicksum(c * vs[n] for c, n in terms)
        m.addConstr(e <= rhs if rel == "<=" else e >= rhs if rel == ">=" else e == rhs, name=cn)
    m.setParam("TimeLimit", float(limit))
    m.setParam("Threads", 1)
    build = time.perf_counter() - t0

    t0 = time.perf_counter()
    rc = m.optimize()
    solve = time.perf_counter() - t0
    status = "optimal" if rc == mp.Status.OPTIMAL else f"status{rc}"
    obj = m.ObjVal if m.SolCount > 0 else None
    return build, solve, status, obj


# ------------------------------------------------------------------ PuLP-CBC --
def run_pulp(p, limit):
    import pulp

    t0 = time.perf_counter()
    cat = {"B": "Binary", "I": "Integer", "C": "Continuous"}
    vs = {nm: pulp.LpVariable(nm, cat=cat[vt]) for nm, vt in p.vars}
    prob = pulp.LpProblem(p.name, pulp.LpMaximize if p.sense == "max" else pulp.LpMinimize)
    prob += pulp.lpSum(c * vs[n] for c, n in p.obj)
    for cn, terms, rel, rhs in p.cons:
        e = pulp.lpSum(c * vs[n] for c, n in terms)
        prob += (e <= rhs if rel == "<=" else e >= rhs if rel == ">=" else e == rhs), cn
    build = time.perf_counter() - t0

    t0 = time.perf_counter()
    prob.solve(pulp.PULP_CBC_CMD(msg=False, threads=1, timeLimit=limit))
    solve = time.perf_counter() - t0
    status = "optimal" if pulp.LpStatus[prob.status] == "Optimal" else pulp.LpStatus[prob.status]
    obj = pulp.value(prob.objective) if prob.objective is not None else None
    return build, solve, status, obj


# ---------------------------------------------------------------- SciPy-HiGHS --
def run_scipy(p, limit):
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds

    t0 = time.perf_counter()
    idx = {nm: i for i, (nm, _) in enumerate(p.vars)}
    n = len(p.vars)
    c = np.zeros(n)
    for coef, nm in p.obj:
        c[idx[nm]] += coef
    if p.sense == "max":
        c = -c
    A = np.zeros((len(p.cons), n))
    lb = np.full(len(p.cons), -np.inf)
    ub = np.full(len(p.cons), np.inf)
    for r, (_cn, terms, rel, rhs) in enumerate(p.cons):
        for coef, nm in terms:
            A[r, idx[nm]] += coef
        if rel == "<=":
            ub[r] = rhs
        elif rel == ">=":
            lb[r] = rhs
        else:
            lb[r] = ub[r] = rhs
    integ = np.array([0 if vt == "C" else 1 for _, vt in p.vars])
    vlb = np.zeros(n)
    vub = np.array([1.0 if vt == "B" else np.inf for _, vt in p.vars])
    build = time.perf_counter() - t0

    t0 = time.perf_counter()
    res = milp(c=c, constraints=LinearConstraint(A, lb, ub), integrality=integ,
               bounds=Bounds(vlb, vub),
               options={"time_limit": limit, "mip_rel_gap": 0.0})
    solve = time.perf_counter() - t0
    status = {0: "optimal", 1: "limit", 4: "limit"}.get(res.status, f"status{res.status}")
    obj = None
    if res.fun is not None:
        obj = -res.fun if p.sense == "max" else res.fun
    return build, solve, status, obj


# -------------------------------------------------------------- OR-Tools-SCIP --
def run_ortools(p, limit):
    from ortools.linear_solver import pywraplp

    t0 = time.perf_counter()
    s = pywraplp.Solver.CreateSolver("SCIP")
    if s is None:
        raise RuntimeError("SCIP 不可用")
    mk = {"B": lambda nm: s.BoolVar(nm),
          "I": lambda nm: s.IntVar(0, s.infinity(), nm),
          "C": lambda nm: s.NumVar(0, s.infinity(), nm)}
    vs = {nm: mk[vt](nm) for nm, vt in p.vars}
    s.SetTimeLimit(int(limit * 1000))
    try:
        s.SetNumThreads(1)
    except AttributeError:
        pass  # SCIP 默认串行
    for cn, terms, rel, rhs in p.cons:
        lo, hi = (-s.infinity(), rhs) if rel == "<=" else \
                 (rhs, s.infinity()) if rel == ">=" else (rhs, rhs)
        ct = s.Constraint(lo, hi, cn)
        for coef, nm in terms:
            ct.SetCoefficient(vs[nm], coef)
    obj = s.Objective()
    for coef, nm in p.obj:
        obj.SetCoefficient(vs[nm], coef)
    obj.SetMaximization() if p.sense == "max" else obj.SetMinimization()
    build = time.perf_counter() - t0

    t0 = time.perf_counter()
    rc = s.Solve()
    solve = time.perf_counter() - t0
    status = "optimal" if rc == pywraplp.Solver.OPTIMAL else (
        "limit" if rc == pywraplp.Solver.FEASIBLE else f"status{rc}")
    val = obj.Value() if rc in (pywraplp.Solver.OPTIMAL, pywraplp.Solver.FEASIBLE) else None
    return build, solve, status, val


SOLVERS = [("micpython", run_micpython), ("cbc", run_pulp),
           ("highs", run_scipy), ("scip", run_ortools)]


def run_pair(name, sname):
    limit = dict(INSTANCES)[name]
    p = parse_milp(INST / f"{name}.milp")
    fn = dict(SOLVERS)[sname]
    build, solve, status, obj = fn(p, limit)
    return {"instance": name, "solver": sname, "build_secs": round(build, 4),
            "solve_secs": round(solve, 3), "status": status, "objective": obj,
            "limit": limit}


def write_results(rows):
    mismatch = []
    by_inst = {}
    for r in rows:
        by_inst.setdefault(r["instance"], []).append(r)
    for name, rs in by_inst.items():
        objs = [r["objective"] for r in rs if r["objective"] is not None]
        if len(objs) >= 2:
            lo, hi = min(objs), max(objs)
            if hi - lo > 1e-6 * (1 + abs(hi)):
                mismatch.append((name, [(r["solver"], r["objective"], r["status"]) for r in rs]))
        opt_file = INST / f"{name}.opt"
        if opt_file.exists():
            dp = float(opt_file.read_text().strip())
            for r in rs:
                if r["objective"] is not None and abs(r["objective"] - dp) > 1e-6:
                    mismatch.append((name + " [DP]", [(r["solver"], r["objective"], dp)]))

    (ROOT / "bench" / "python_bench_results.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=1), encoding="utf-8")

    lines = ["# micpython × 开源求解器（CBC / HiGHS / SCIP）同机对比", "",
             "环境：conda `mic-bench`（py3.11，pip 的 pulp/scipy/ortools 官方轮子），threads=1，"
             "时限见正文；机器 15W 笔记本（墙钟噪声 ±20-60%，见 docs/parallel-profiling.md §1.3）。"
             "solve 单位秒；`watchdog_kill` = 求解器时限+120s 宽限内未返回被强杀"
             "（实测 CBC 在 mknap/marketsplit 类会进入不受时限约束的阶段）。", "",
             "| 算例 | micpython | CBC | HiGHS | SCIP | 目标一致 |",
             "|---|---|---|---|---|---|"]
    for name, _limit in INSTANCES:
        cells, objs = [], []
        for sname, _fn in SOLVERS:
            r = next(x for x in rows if x["instance"] == name and x["solver"] == sname)
            cells.append(f"{r['solve_secs']!s}/{r['status']}" if r["solve_secs"] is not None
                         else r["status"])
            if r["objective"] is not None:
                objs.append(r["objective"])
        uniq = len(set(objs))
        if len(objs) < 2 or uniq == 1:
            agree = "✓" if uniq == 1 or len(objs) == 4 else "—"
        else:
            lo, hi = min(objs), max(objs)
            agree = "✓" if hi - lo <= 1e-6 * (1 + abs(hi)) else "△"
        lines.append(f"| {name} | " + " | ".join(cells) + f" | {agree} |")
    if mismatch:
        lines += ["", "## MISMATCH（目标不一致，需人工核查）", ""]
        for name, detail in mismatch:
            lines.append(f"- {name}: {detail}")
    (ROOT / "bench" / "python_bench_results.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")
    print(f"done; mismatch={len(mismatch)}", flush=True)


def main_driver():
    rows = []
    for name, limit in INSTANCES:
        for sname, _fn in SOLVERS:
            cap = limit + GRACE_SECS
            try:
                proc = subprocess.run(
                    [sys.executable, "-u", str(Path(__file__).resolve()), "--pair", name, sname],
                    capture_output=True, text=True, timeout=cap, cwd=str(ROOT))
                row = json.loads(proc.stdout.strip().splitlines()[-1])
            except subprocess.TimeoutExpired:
                row = {"instance": name, "solver": sname, "build_secs": None,
                       "solve_secs": round(cap, 1), "status": "watchdog_kill",
                       "objective": None, "limit": limit}
            except Exception as e:  # noqa: BLE001 —— 解析失败/非零退出等
                row = {"instance": name, "solver": sname, "build_secs": None,
                       "solve_secs": None, "status": f"error: {type(e).__name__}",
                       "objective": None, "limit": limit}
            rows.append(row)
            print(f"{name:20s} {sname:9s} solve={row['solve_secs']!s:>8} "
                  f"{row['status']:13s} obj={row['objective']}", flush=True)
    write_results(rows)


if __name__ == "__main__":
    if len(sys.argv) >= 4 and sys.argv[1] == "--pair":
        print(json.dumps(run_pair(sys.argv[2], sys.argv[3]), ensure_ascii=False))
    elif len(sys.argv) >= 3 and sys.argv[1] == "--merge":
        # --merge <jsonl>：以现有 json 为底，用 jsonl 中的行（--pair 输出）替换同键行后重写
        base = json.loads((ROOT / "bench" / "python_bench_results.json").read_text(encoding="utf-8"))
        key = lambda r: (r["instance"], r["solver"])  # noqa: E731
        repl = {}
        for line in Path(sys.argv[2]).read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                repl[key(r)] = r
        rows = [repl.get(key(r), r) for r in base]
        write_results(rows)
    else:
        main_driver()
