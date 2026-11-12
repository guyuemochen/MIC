#!/usr/bin/env python3
"""生成自动配置独立验证集；不求解、不改旧算例、不生成未经证明的最优值。"""
import argparse
import hashlib
import json
import random
from pathlib import Path

# 参数与种子在首次性能评测前冻结；难度标签仅为结构预期，不是测得结果。
PACKING = [
    ("packing_small", 80, 2, 1, 2, False),
    ("packing_medium", 160, 4, 1, 5, False),
    ("packing_many_rows", 240, 7, 1, 8, False),
    ("packing_correlated", 300, 5, 1, 6, True),
]
SEEDS = [12017, 12043]


def encode(name, sense, n, objective, constraints):
    lines = [f"name {name}", f"sense {sense}", "vars"]
    lines += [f"  x{i + 1} binary" for i in range(n)]
    lines += ["obj", "  " + " ".join(f"{v} x{i + 1}" for i, v in enumerate(objective)), "cons"]
    for j, (terms, op, rhs) in enumerate(constraints):
        body = " ".join(f"{v} x{i + 1}" for i, v in terms)
        lines.append(f"  c{j + 1}: {body} {op} {rhs}")
    return ("\n".join(lines + ["end"]) + "\n").encode()


def packing(name, n, m, numerator, denominator, correlated, seed, density=1.0, sense="max", negative_objective=False, negative_coefficient=False, greater=False):
    rng = random.Random(seed)
    weights = [[rng.randint(7, 47) for _ in range(n)] for _ in range(m)]
    values = ([sum(row[i] for row in weights) // m + rng.randint(1, 8) for i in range(n)]
              if correlated else [rng.randint(10, 99) for _ in range(n)])
    if negative_objective:
        values[0] = -values[0]
    rows = []
    for j, row in enumerate(weights):
        cols = list(range(n)) if density == 1.0 else sorted(rng.sample(range(n), max(1, round(n * density))))
        terms = [(i, -row[i] if negative_coefficient and j == 0 and i == cols[0] else row[i]) for i in cols]
        rows.append((terms, ">=" if greater else "<=", sum(row[i] for i in cols) * numerator // denominator))
    return encode(name, sense, n, values, rows)


def covering(name, n, m, seed):
    rng = random.Random(seed)
    values = [rng.randint(1, 20) for _ in range(n)]
    rows = [(sorted((i, 1) for i in rng.sample(range(n), max(2, n // 10))), ">=", 1)
            for _ in range(m)]
    return encode(name, "min", n, values, rows)


def payloads():
    for family, n, m, num, den, corr in PACKING:
        for seed in SEEDS:
            name = f"{family}_s{seed}"
            params = dict(family=family, n=n, m=m, capacity_fraction=[num, den], correlated=corr, seed=seed)
            yield name, params, packing(name, n, m, num, den, corr, seed)
    for family, n, m in [("boundary_n199", 199, 2), ("boundary_n200", 200, 2), ("boundary_n512_m16", 512, 16), ("boundary_n513_m16", 513, 16), ("boundary_m1", 240, 1), ("boundary_m17", 240, 17)]:
        seed = 12071
        params = dict(family=family, n=n, m=m, capacity_fraction=[1, 3], correlated=False, seed=seed)
        yield family, params, packing(family, n, m, 1, 3, False, seed)
    for family, extra in [("sparse_packing", {"density": 0.25}), ("minimization", {"sense": "min"}),
                          ("negative_objective", {"negative_objective": True}),
                          ("negative_coefficient", {"negative_coefficient": True}), ("ge_rows", {"greater": True})]:
        seed = 12203
        params = dict(family=family, n=220, m=4, capacity_fraction=[1, 5], correlated=False, seed=seed, **extra)
        yield family, params, packing(family, 220, 4, 1, 5, False, seed, **extra)
    for n, m, seed in [(100, 60, 12101), (200, 120, 12109)]:
        name = f"covering_{m}x{n}_s{seed}"
        yield name, dict(family="covering", n=n, m=m, row_nnz=n // 10, seed=seed), covering(name, n, m, seed)


def write_once(path, data):
    if path.exists():
        if path.read_bytes() != data:
            raise SystemExit(f"拒绝覆盖不同内容：{path}")
    else:
        with path.open("xb") as target:
            target.write(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    if any(p.name == "instances" and p.parent.name == "bench" for p in [output, *output.parents]):
        raise SystemExit("拒绝写入原 bench/instances")
    output.mkdir(parents=True, exist_ok=True)
    manifest = {"schema": 1, "generator": "bench/auto_validation/generate.py", "instances": []}
    for name, params, data in payloads():
        path = output / f"{name}.milp"
        write_once(path, data)
        manifest["instances"].append(dict(file=path.name, sha256=hashlib.sha256(data).hexdigest(),
                                           parameters=params, known_optimum=None))
    encoded = (json.dumps(manifest, ensure_ascii=False, indent=2) + "\n").encode()
    write_once(output / "manifest.json", encoded)
    print(f"已生成或逐字节核验 {len(manifest['instances'])} 个算例：{output}")


if __name__ == "__main__":
    main()
