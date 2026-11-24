"""从归档首解独立枚举单课改位置及双课互换，并用整表评分复核每一步。"""
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
sys.path.insert(0, str(BENCH))
from cases import check, read


def fingerprint(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def neighborhood(rows, shape, positions):
    """从可行解出发，目标时段的教师、班级和教室均需空闲或由参与移动的课程占用。"""
    used = {}
    def keys(i, time, room):
        teacher, cohort, _ = rows[i]
        return ((time, "teacher", teacher), (time, "cohort", cohort), (time, "room", room))

    for i, (time, room) in enumerate(positions):
        for key in keys(i, time, room):
            assert key not in used
            used[key] = i

    def allowed(i, target, moving):
        return all(k not in used or used[k] in moving for k in keys(i, *target))

    slots, rooms = shape
    for i, (_, _, pref) in enumerate(rows):
        old_time, old_room = positions[i]
        for t, r in product(range(slots), range(rooms)):
            if (t, r) != (old_time, old_room):
                delta = abs(t - pref) - abs(old_time - pref)
                yield "relocate", ((i, t, r),), delta, allowed(i, (t, r), (i,))

    for i, j in combinations(range(len(rows)), 2):
        a, b = positions[i], positions[j]
        p, q = rows[i][2], rows[j][2]
        delta = abs(b[0] - p) + abs(a[0] - q) - abs(a[0] - p) - abs(b[0] - q)
        # 互换时，目标位置彼此不同。若目标时段相同，两个参与课程原本已共处该时段。
        ok = allowed(i, b, (i, j)) and allowed(j, a, (i, j))
        yield "swap", ((i, *b), (j, *a)), delta, ok


def candidates(rows, shape, positions):
    stats = {k: {"total": 0, "feasible": 0, "improving": 0, "equal": 0, "worse": 0}
             for k in ("relocate", "swap")}
    gains = []
    for kind, change, delta, ok in neighborhood(rows, shape, positions):
        counts = stats[kind]
        counts["total"] += 1
        if ok:
            counts["feasible"] += 1
            if delta < 0:
                counts["improving"] += 1
                gains.append((delta, change))
            elif delta == 0:
                counts["equal"] += 1
            else:
                counts["worse"] += 1
    return stats, sorted(gains)


def scalar_swaps(instance, rows, shape, positions):
    """原求解器还可交换两个标量变量，包括时段和教室交叉交换。"""
    def keys(i, t, r):
        return ((0, t, rows[i][0]), (1, t, rows[i][1]), (2, t, r))
    occupied = {k: i for i, (t, r) in enumerate(positions) for k in keys(i, t, r)}
    counts = {"total_valid_nontrivial": 0, "feasible": 0, "improving": 0, "equal": 0, "worse": 0}
    gains = []
    start = sum(abs(t - row[2]) for (t, _), row in zip(positions, rows))
    for a, b in combinations(range(len(positions) * 2), 2):
        i, x = divmod(a, 2)
        j, y = divmod(b, 2)
        va, vb = positions[i][x], positions[j][y]
        if va == vb or vb >= shape[x] or va >= shape[y]:
            continue
        counts["total_valid_nontrivial"] += 1
        new = {k: list(positions[k]) for k in (i, j)}
        new[i][x], new[j][y] = vb, va
        targets = [k for idx, (t, r) in new.items() for k in keys(idx, t, r)]
        if (len(set(targets)) != len(targets) or
                any(k in occupied and occupied[k] not in new for k in targets)):
            continue
        counts["feasible"] += 1
        delta = sum(abs(t - rows[k][2]) - abs(positions[k][0] - rows[k][2])
                    for k, (t, _) in new.items())
        if delta < 0:
            counts["improving"] += 1
            gains.append({"variables": (a, b), "delta": delta})
            apply_checked(instance, positions, [(k, *v) for k, v in new.items()], start + delta)
        elif delta == 0:
            counts["equal"] += 1
        else:
            counts["worse"] += 1
    return counts, gains


def apply_checked(instance, positions, change, expected):
    result = positions.copy()
    for i, time, room in change:
        result[i] = (time, room)
    assert check(instance, [x for pair in result for x in pair], []) == (0, expected)
    return result


def main():
    instance = BENCH / "instances/timetable.txt"
    archive = HERE.parent / "paired-scoring/results/runs.jsonl"
    records = [r for r in map(json.loads, archive.read_text().splitlines())
               if (r["stage"], r["engine"], r["case"]) == ("measured", "after", "timetable")]
    assert [r["seed"] for r in records] == [7, 19, 43, 73, 101]
    start = records[0]["first_values"]
    for r in records:
        assert r["first_values"] == r["values"] == start
        assert check(instance, r["first_values"], r["first_lists"]) == (0, 109)
    _, shape, rows = read(instance)
    positions = list(zip(start[::2], start[1::2]))
    output = {"seeds": [r["seed"] for r in records], "input_sha256": fingerprint(instance),
              "runs_sha256": fingerprint(archive), "same_start_and_final_for_all_seeds": True,
              "initial_soft": 109, "strict_descent": []}
    counts, wins = scalar_swaps(instance, rows, shape, positions)
    output["initial_scalar_swaps"] = counts
    output["initial_improving_scalar_swaps"] = wins
    score = 109
    while True:
        stats, gains = candidates(rows, shape, positions)
        if score == 109:
            output["initial_neighborhood"] = stats
            output["initial_improving_moves"] = [{"delta": d, "course_time_room": c}
                                                 for d, c in gains]
            for d, c in gains:
                apply_checked(instance, positions, c, score + d)
        if not gains:
            output["last_neighborhood"] = stats
            break
        delta, change = gains[0]
        positions = apply_checked(instance, positions, change, score + delta)
        score += delta
        output["strict_descent"].append({"course_time_room": change, "soft": score})
    output["final_soft"] = score
    # 旧版只在抽中成对移动并选择交换时才可能抽中这六项中的一项。
    assert output["initial_scalar_swaps"]["improving"] == 0
    assert output["initial_neighborhood"]["relocate"]["improving"] == 0
    swaps = output["initial_neighborhood"]["swap"]
    probability = 0.25 * 0.75 * swaps["improving"] / swaps["total"]
    output["old_sampler_at_initial_state"] = {
        "single_draw_improves_probability": probability,
        "at_least_one_in_32_independent_draws_probability": 1 - (1 - probability) ** 32,
    }
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
