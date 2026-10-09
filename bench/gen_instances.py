#!/usr/bin/env python3
"""生成 milp 大规模基准算例（确定性种子，可复现）。

规模轴：变量数（LP 规模）与组合难度（B&B 节点数）双重覆盖。
背包类用整数重量/价值，便于用 DP 独立核对最优解。
"""
import random
from pathlib import Path

OUT = Path(__file__).parent / "instances"


def write(name: str, sense: str, var_decls, obj_terms, cons):
    lines = [f"name {name}", f"sense {sense}", "vars"]
    lines += [f"  {v}" for v in var_decls]
    lines.append("obj")
    for i in range(0, len(obj_terms), 8):
        lines.append("  " + " ".join(f"{c} {v}" for c, v in obj_terms[i:i + 8]))
    lines.append("cons")
    for cname, terms, rel, rhs in cons:
        t = " ".join(f"{c} {v}" for c, v in terms)
        lines.append(f"  {cname}: {t} {rel} {rhs}")
    lines.append("end")
    (OUT / f"{name}.milp").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"{name}: {len(var_decls)} vars, {len(cons)} cons")


def knapsack(n: int, seed: int):
    rng = random.Random(seed)
    vals, weights = [], []
    for _ in range(n):
        w = rng.randint(7, 47)
        weights.append(w)
        vals.append(w + rng.randint(1, 10))  # 强相关背包（难）
    cap = sum(weights) // 2
    write(f"knapsack-{n}", "max", [f"x{i+1} binary" for i in range(n)],
          [(vals[i], f"x{i+1}") for i in range(n)],
          [("cap", [(weights[i], f"x{i+1}") for i in range(n)], "<=", cap)])
    # DP 交叉验证最优值
    dp = [0] * (cap + 1)
    for i in range(n):
        w, v = weights[i], vals[i]
        for c in range(cap, w - 1, -1):
            cand = dp[c - w] + v
            if cand > dp[c]:
                dp[c] = cand
    (OUT / f"knapsack-{n}.opt").write_text(str(dp[cap]), encoding="utf-8")
    print(f"  DP 最优 = {dp[cap]}, cap = {cap}")


def mknap(n: int, m: int, seed: int):
    rng = random.Random(seed)
    vals = [rng.randint(10, 99) for _ in range(n)]
    ws = [[rng.randint(7, 47) for _ in range(n)] for _ in range(m)]
    cons = [(f"k{j+1}", [(ws[j][i], f"x{i+1}") for i in range(n)], "<=",
             sum(ws[j]) // (m + 1)) for j in range(m)]
    write(f"mknap-{m}x{n}", "max", [f"x{i+1} binary" for i in range(n)],
          [(vals[i], f"x{i+1}") for i in range(n)], cons)


def assign(k: int, seed: int):
    rng = random.Random(seed)
    names = lambda i, j: f"x_{i}_{j}"
    write(f"assign-{k}x{k}", "min",
          [f"{names(i, j)} binary" for i in range(1, k + 1) for j in range(1, k + 1)],
          [(rng.randint(1, 99), names(i, j)) for i in range(1, k + 1) for j in range(1, k + 1)],
          [(f"r{i}", [(1, names(i, j)) for j in range(1, k + 1)], "=", 1) for i in range(1, k + 1)]
          + [(f"c{j}", [(1, names(i, j)) for i in range(1, k + 1)], "=", 1) for j in range(1, k + 1)])


def mixed(n_int: int, n_cont: int, m: int, seed: int):
    """一般 MILP：整数变量与连续变量耦合（等式含连续 ⇒ LP 松弛分数，逼真分支）。"""
    rng = random.Random(seed)
    vs = [f"y{i+1} integer" for i in range(n_int)] + [f"t{i+1} continuous" for i in range(n_cont)]
    terms = [(rng.randint(3, 19), f"y{i+1}") for i in range(n_int)] + \
            [(rng.choice([5, 10, 25, 50]), f"t{i+1}") for i in range(n_cont)]
    cons = []
    for j in range(m):
        row = [(rng.randint(2, 15), f"y{i+1}") for i in range(n_int)]
        row += [(rng.choice([1, 2, 4]), f"t{i+1}") for i in range(n_cont)]
        rhs = rng.randint(20, 40) + m * 3
        rel = "=" if j % 2 == 0 else "<="
        cons.append((f"row{j+1}", row, rel, rhs))
    write(f"milp-mixed-{n_int}i{n_cont}c", "max", vs, terms, cons)


def lp_wide(n: int, m: int, seed: int):
    """纯 LP（无整数变量）：考验单纯形本体吞吐。"""
    rng = random.Random(seed)
    vs = [f"z{i+1} continuous" for i in range(n)]
    obj = [(rng.randint(3, 19), f"z{i+1}") for i in range(n)]
    cons = []
    for j in range(m):
        row = [(rng.randint(1, 9), f"z{i+1}") for i in range(n)]
        cons.append((f"r{j+1}", row, "<=", rng.randint(30, 60)))
    write(f"lp-wide-{n}x{m}", "max", vs, obj, cons)


def market_split(m: int, n: int, seed: int, lo: int = 100, hi: int = 200):
    """等式型市场分割（market-split，等式背包；Cornuejols/Aardal 经典硬类）。

    m 条等式行 sum_i a_ij x_i = b_j，x ∈ {0,1}，系数 a_ij ∈ [lo,hi] 整数；
    可行性由 planted 子集保证：随机选 |S| = n//2 列，b_j = sum_{i∈S} a_ij，
    目标 min sum x_i（求最小支撑的市场分割）。LP 松弛对等式行界极弱
    （0≤x≤1 内几乎总能配平），B&B 树大——Gurobi 时代经典Benchmark 类。
    """
    rng = random.Random(seed)
    a = [[rng.randint(lo, hi) for _ in range(n)] for _ in range(m)]
    planted = sorted(rng.sample(range(n), n // 2))
    cons = []
    for j in range(m):
        rhs = sum(a[j][i] for i in planted)
        cons.append((f"m{j+1}", [(a[j][i], f"x{i+1}") for i in range(n)], "=", rhs))
    write(f"marketsplit-{m}x{n}", "min", [f"x{i+1} binary" for i in range(n)],
          [(1, f"x{i+1}") for i in range(n)], cons)
    print(f"  planted |S| = {len(planted)}（可行解目标上界）")


def set_cover(m: int, n: int, density: float, seed: int):
    """集合覆盖（set-cover）：min c·x s.t. A x >= 1，A 为 0/1 随机矩阵。

    m 条覆盖行 × n 个二元列，密度 density（每行期望 density*n 个 1）；
    目标系数 c ∈ [1,100]（随机费用版 SCP；单位费用版 LP 界过强）。
    每行至少补 1 个 1 保证可行。
    """
    rng = random.Random(seed)
    cons = []
    for j in range(m):
        cols = [i for i in range(n) if rng.random() < density]
        if not cols:
            cols = [rng.randrange(n)]
        cons.append((f"c{j+1}", [(1, f"x{i+1}") for i in cols], ">=", 1))
    costs = [rng.randint(1, 100) for _ in range(n)]
    write(f"setcover-{m}x{n}", "min", [f"x{i+1} binary" for i in range(n)],
          [(costs[i], f"x{i+1}") for i in range(n)], cons)
    print(f"  density = {density}（每行期望 {density * n:.0f} 个 1）")


OUT.mkdir(parents=True, exist_ok=True)
# —— 既有 7 算例（第 1~2 轮优化基准；种子不变即字节不变）——
knapsack(200, 20261008)
knapsack(500, 20261009)
mknap(100, 3, 20261010)
assign(16, 20261011)
mixed(120, 30, 14, 20261012)
mixed(250, 40, 18, 20261014)
lp_wide(200, 300, 20261013)
# —— 第 3 轮新增（Gurobi 时代经典难类；规模经 90s/threads=1 终态标定，
#    2026-10-09 探针实测：ms-2x20 0.012s 最优 / ms-3x30 1.4s / ms-3x40 1.7s /
#    ms-4x30 5.5s 最优（39.5 万节点）/ ms-4x40 49s 打满 200 万节点上限且
#    零可行解（超规格，不入列）；setcover ≤60x100 LP 常直接 integral，
#    ≥80x120 根 LP 触单纯形迭代上限 max(10k, 20(m+n)) 整体终止——
#    setcover-100x150 作为该失败模式的存档样本保留）——
market_split(2, 20, 20261015)
market_split(4, 30, 20261030)
set_cover(60, 100, 0.10, 20261022)
set_cover(100, 150, 0.10, 20261017)
knapsack(2000, 20261019)
mknap(300, 5, 20261020)
