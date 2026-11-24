"""重新验分并输出正式试验的中位数、范围和配对结果。"""
import argparse
from collections import Counter
import json
from pathlib import Path
from statistics import median

from cases import CASES, read
from run import HERE, digest, validate


NAMES = dict(roster="排班", timetable="排课", delivery="配送", factory="生产排程")
ENGINES = ("MIC", "OptaPlanner")


def describe(case, path):
    _, head, rows = read(path)
    if case == "roster":
        return f"{len(rows)} 班次／{head[0]} 人"
    if case == "timetable":
        return f"{len(rows)} 门课／{head[0]} 时段／{head[1]} 间教室"
    if case == "delivery":
        return f"{len(rows)} 客户／{head[0]} 辆车"
    return f"{len(rows)} 道工序／{head[0]} 台机器"


def summary(folder):
    meta = json.loads((folder / "environment.json").read_text())
    paths = {k: HERE / "instances" / (k + ".txt") for k in CASES}
    for case, path in paths.items():
        if digest(path) != meta["cases_sha256"][case]:
            raise ValueError("算例已变更：" + case)
    rows = {}
    warmups = set()
    elapsed = []
    for line in (folder / "runs.jsonl").read_text().splitlines():
        r = json.loads(line)
        key = r["case"], r["engine"], r["seed"]
        validate(paths[r["case"]], r, r["budget_ms"])
        if r["stage"] == "warmup":
            if key in warmups or r["budget_ms"] != meta["warmup_ms"]:
                raise ValueError("预热重复或超出约定预算")
            warmups.add(key)
        elif r["stage"] == "measured":
            if key in rows or r["budget_ms"] != meta["budget_ms"]:
                raise ValueError("正式记录重复或超出约定预算")
            rows[key] = r
            elapsed.append(r["elapsed_ms"])
        else:
            raise ValueError("未知阶段")
    seeds = meta["seeds"]
    def keys(seed_list):
        return {(case, engine, seed) for case in CASES for engine in ENGINES for seed in seed_list}
    if set(rows) != keys(seeds) or warmups != keys(meta["warmup_seeds"]):
        raise ValueError("配对次数或预热次数不足")
    print("# MIC 与 OptaPlanner 等时实测\n")
    print(f'{meta["cpu"]}，单求解线程。每类一个固定算例，每个引擎每次 '
          f'{meta["budget_ms"]/1000:g} 秒，同一算例使用 {len(seeds)} 个配对种子。')
    print(f'正式求解 {len(rows)} 次，其中 {sum(r["checked_feasible"] for r in rows.values())} 次'
          '在截止前找到完整可行解。每个首解和每个最终可行解均已独立验分。\n')
    print("## 首次可行和截止解质量\n")
    print("首次可行列为毫秒，中位数后面的括号给出最小值至最大值。截止代价取各自成功试验的中位数，"
          "越小越好。差值等于 MIC 减 OptaPlanner，百分比以 OptaPlanner 的中位代价为分母。\n")
    print("| 算例与规模 | 可行次数 MIC／Opta | MIC 首次可行 ms | Opta 首次可行 ms | MIC 截止代价 | Opta 截止代价 | 差值 |")
    print("|---|---:|---:|---:|---:|---:|---:|")
    for case in CASES:
        groups = [[rows[case, eng, s] for s in seeds if rows[case, eng, s]["checked_feasible"]]
                  for eng in ENGINES]
        times = [[r["first_ms"] for r in group] for group in groups]
        tm = ["%.2f（%.2f–%.2f）" % (median(v), min(v), max(v)) if v else "—" for v in times]
        costs = [median([r["soft"] for r in group]) if group else None for group in groups]
        gap = "—"
        if None not in costs:
            difference = costs[0] - costs[1]
            gap = "%+g" % difference if difference else "0"
            if costs[1] != 0:
                gap += "（%+.1f%%）" % (100 * difference / costs[1]) if difference else "（0%）"
        cols = [NAMES[case] + "，" + describe(case, paths[case]),
                "／".join(f"{len(g)}/{len(seeds)}" for g in groups), *tm,
                *["—" if c is None else "%g" % c for c in costs], gap]
        print("| " + " | ".join(cols) + " |")
    print("\n## 同一种子的最终代价\n")
    print("每格依次为 MIC／OptaPlanner。若截止前未找到可行解，记为“未可行”。\n")
    print("| 算例 | " + " | ".join(map(str, seeds)) + " | MIC 更低／持平／Opta 更低 |")
    print("|---|" + "---:|" * (len(seeds) + 1))
    for case in CASES:
        cols, wins = [], Counter()
        for seed in seeds:
            a, b = [rows[case, e, seed] for e in ENGINES]
            cols.append("／".join(str(r["soft"]) if r["checked_feasible"] else "未可行" for r in (a,b)))
            if a["checked_feasible"] and b["checked_feasible"]:
                wins[(a["soft"] > b["soft"]) - (a["soft"] < b["soft"])] += 1
        print("| " + " | ".join([NAMES[case], *cols, "／".join(str(wins[i]) for i in (-1,0,1))]) + " |")
    # 每台机器必须依次执行分给它的全部工序，其总负荷给出可独立核验的完工时间下界。
    _, _, operations = read(paths["factory"])
    load = Counter()
    for machine, duration, _ in operations:
        load[machine] += duration
    lower = max(load.values())
    hits = [sum(rows["factory", e, s]["checked_feasible"] and rows["factory", e, s]["soft"] == lower
                for s in seeds) for e in ENGINES]
    print(f"\n生产实例最忙的机器需要连续加工 {lower} 个时间单位，这给出完工时间下界。"
          f"MIC 有 {hits[0]} 次、OptaPlanner 有 {hits[1]} 次达到该下界。")
    print("\n## 从首个可行解到截止时的改进\n")
    print("首解和截止代价分别取各次成功试验的中位数。\n")
    print("| 算例 | MIC 首解 → 截止 | OptaPlanner 首解 → 截止 |")
    print("|---|---:|---:|")
    for case in CASES:
        cols = []
        for engine in ENGINES:
            ok = [rows[case, engine, s] for s in seeds if rows[case, engine, s]["checked_feasible"]]
            cols.append("%g → %g" % (median(r["first_soft_checked"] for r in ok),
                                      median(r["soft"] for r in ok)) if ok else "未可行")
        print("| " + " | ".join([NAMES[case], *cols]) + " |")
    print("\n## 每次首次可行时间\n")
    print("单位为毫秒，每格依次为 MIC／OptaPlanner。\n")
    print("| 算例 | " + " | ".join(map(str, seeds)) + " |")
    print("|---|" + "---:|" * len(seeds))
    for case in CASES:
        cols = []
        for seed in seeds:
            times = [rows[case, e, seed]["first_ms"] for e in ENGINES]
            cols.append("／".join("未可行" if t is None else "%.2f" % t for t in times))
        print("| " + " | ".join([NAMES[case], *cols]) + " |")
    print(f'\n正式调用从求解入口到返回实际用时 {min(elapsed)/1000:.6f} 至 {max(elapsed)/1000:.6f} 秒。'
          '仅使用在预算内观测到的解。每类仅一个实例，范围反映重复搜索及墙钟波动。\n')
    print("[方法与复现](../../README.md) · [逐次原始解及时间](runs.jsonl) · [环境与输入指纹](environment.json)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    summary(parser.parse_args().folder)
