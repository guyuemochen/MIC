"""冻结算例的生成器和独立验分器。见证只用于自检，不写入求解器输入。"""
import argparse
from collections import Counter, deque
import hashlib
from pathlib import Path
import random


SEED = 20261124
HERE = Path(__file__).resolve().parent
CASES = ("roster", "timetable", "delivery", "factory")


def _roster(rng):
    workers, limit, slots = 18, 9, 24
    shifts = []
    skills = [{w % 5, (w + 1 + w // 5) % 5} for w in range(workers)]
    favored = [rng.randrange(4) for _ in range(workers)]
    for block in range(slots // 3):
        assigned = list(range(workers))
        rng.shuffle(assigned)
        for i, witness in enumerate(assigned):
            slot = 3 * block + i // 6
            skill = rng.choice(sorted(skills[witness]))
            costs = []
            for w in range(workers):
                available = skill in skills[w] and (w == witness or rng.random() < 0.84)
                costs.append(2 * abs(slot % 4 - favored[w]) + rng.randrange(6)
                             if available else -1)
            shifts.append(([slot] + costs, witness))
    rng.shuffle(shifts)
    return [workers, limit], [row for row, _ in shifts], [w for _, w in shifts], []


def _timetable(rng):
    slots, rooms, teachers, cohorts = 20, 8, 16, 14
    lessons = []
    for slot in range(slots):
        teacher_ids = rng.sample(range(teachers), 7)
        cohort_ids = rng.sample(range(cohorts), 7)
        room_ids = rng.sample(range(rooms), 7)
        for teacher, cohort, room in zip(teacher_ids, cohort_ids, room_ids):
            preferred = (rng.choice((4, 5, 10, 11)) if rng.random() < 0.30
                         else max(0, min(slots - 1, slot + rng.randint(-4, 4))))
            lessons.append(([teacher, cohort, preferred], [slot, room]))
    rng.shuffle(lessons)
    return [slots, rooms], [row for row, _ in lessons], [v for _, pair in lessons for v in pair], []


def _delivery(rng):
    vehicles, cap = 8, 48
    points = []
    coords = {(50, 50)}
    centers = [(16, 15), (26, 72), (75, 29), (79, 82), (49, 40), (12, 48)]
    for route in range(vehicles):
        loads = [2, 3, 4, 4, 5, 6, 6, 7, 8]
        rng.shuffle(loads)
        for d in loads:
            while True:
                cx, cy = rng.choice(centers)
                pt = (cx + rng.randint(-10, 10), cy + rng.randint(-10, 10))
                if pt not in coords:
                    break
            coords.add(pt)
            points.append(([pt[0], pt[1], d], route))
    rng.shuffle(points)
    routes = [[] for _ in range(vehicles)]
    for i, (_, owner) in enumerate(points):
        routes[owner].append(i)
    for route in routes:
        rng.shuffle(route)
    return [vehicles, cap, 50, 50], [row for row, _ in points], [], routes


def _job_order(rng, n_jobs, n_ops):
    # 每次只从已经就绪的作业中选一道工序，所得全序一定满足工艺前置关系。
    positions = [0] * n_jobs
    active = list(range(n_jobs))
    out = []
    while active:
        j = rng.choice(active)
        out.append((j, positions[j]))
        positions[j] += 1
        if positions[j] == n_ops:
            active.remove(j)
    return out


def _factory(rng):
    n_jobs, machines = 12, 8
    jobs = [[(m, rng.randint(1, 20)) for m in rng.sample(range(machines), machines)]
            for _ in range(n_jobs)]
    numbering = _job_order(rng, n_jobs, machines)
    ids = {op: i for i, op in enumerate(numbering)}
    rows = [list(jobs[j][k]) + [ids[j, k - 1] if k else -1] for j, k in numbering]
    routes = [[] for _ in range(machines)]
    for j, k in _job_order(rng, n_jobs, machines):
        routes[jobs[j][k][0]].append(ids[j, k])
    return [machines], rows, [], routes


def generate(kind):
    """返回文件字节和单独的可行见证。每类实例使用独立随机流。"""
    generators = (_roster, _timetable, _delivery, _factory)
    tags = ("S", "T", "D", "J")
    i = CASES.index(kind)
    args, rows, values, lists = generators[i](random.Random(SEED + i))
    lines = [[kind.upper()] + args] + [[tags[i]] + r for r in rows]
    body = "".join(" ".join(map(str, row)) + "\n" for row in lines).encode("ascii")
    return body, values, lists


def read(path):
    """读取生成器的数字表，返回小写类型名、表头数字、所有项目行的数字。"""
    lines = [s.split() for s in Path(path).read_text(encoding="ascii").splitlines() if s.strip()]
    kind = lines[0][0].lower()
    tag = dict(zip(CASES, ("S", "T", "D", "J")))[kind]
    if any(row[0] != tag for row in lines[1:]):
        raise ValueError("算例项目类型不匹配")
    return kind, list(map(int, lines[0][1:])), [list(map(int, row[1:])) for row in lines[1:]]


def _complete(kind, args, rows, values, lists):
    if kind in ("roster", "timetable"):
        n = len(rows) * (2 if kind == "timetable" else 1)
        return len(values) == n and not lists and all(type(x) is int for x in values)
    if values or len(lists) != args[0] or any(not isinstance(r, (tuple, list)) for r in lists):
        return False
    flat = [x for route in lists for x in route]
    return (all(type(x) is int for x in flat) and
            sorted(flat) == list(range(len(rows))))


def is_complete(path, values, lists):
    """检查完整性。合法值域、资格、硬违反量及循环等待另外由 check 核对。"""
    return _complete(*read(path), values, lists)


def _pairs(counts):
    return sum(n * (n - 1) // 2 for n in counts.values())


def _factory_score(rows, lists):
    # 工艺依赖 + 同机器相邻依赖构成有向图，拓扑顺序最长路径即最早完工时间。
    n = len(rows)
    preds = [set() for _ in rows]
    for i, (_, _, prev) in enumerate(rows):
        if prev != -1:
            preds[i].add(prev)
    for machine, route in enumerate(lists):
        for pos, i in enumerate(route):
            if rows[i][0] != machine:
                raise ValueError("工序分配到错误机器")
            if pos:
                preds[i].add(route[pos - 1])
    succ = [[] for _ in rows]
    degree = [len(p) for p in preds]
    for i, before in enumerate(preds):
        for k in before:
            succ[k].append(i)
    ready = deque(i for i in range(n) if degree[i] == 0)
    ends = [0] * n
    done = 0
    while ready:
        i = ready.popleft()
        ends[i] = rows[i][1] + max((ends[k] for k in preds[i]), default=0)
        done += 1
        for k in succ[i]:
            degree[k] -= 1
            if not degree[k]:
                ready.append(k)
    if done != n:
        raise ValueError("工艺和机器顺序导致循环等待")
    return 0, max(ends, default=0)


def check(path, values, lists):
    """独立重算完整方案的 (硬罚分, 软代价)，两项均越小越好。"""
    kind, args, rows = read(path)
    if not _complete(kind, args, rows, values, lists):
        raise ValueError("缺少变量、列表项目重复/缺失或解的维度不符")
    if kind == "roster":
        n, cap = args
        loads, concurrent, soft = Counter(), Counter(), 0
        for worker, row in zip(values, rows):
            if not 0 <= worker < n or row[worker + 1] < 0:
                raise ValueError("员工无资格承担该班次")
            loads[worker] += 1
            concurrent[row[0], worker] += 1
            soft += row[worker + 1]
        return _pairs(concurrent) + sum(max(0, v - cap) for v in loads.values()), soft
    if kind == "timetable":
        slots, rooms = args
        teachers, cohorts, classrooms = Counter(), Counter(), Counter()
        soft = 0
        for i, (teacher, cohort, pref) in enumerate(rows):
            time, room = values[2 * i:2 * i + 2]
            if not (0 <= time < slots and 0 <= room < rooms):
                raise ValueError("教室或时段越界")
            teachers[time, teacher] += 1
            cohorts[time, cohort] += 1
            classrooms[time, room] += 1
            soft += abs(time - pref)
        return sum(map(_pairs, (teachers, cohorts, classrooms))), soft
    if kind == "delivery":
        _, cap, dx, dy = args
        hard = soft = 0
        for route in lists:
            x, y = dx, dy
            load = 0
            for i in route:
                nx, ny, demand = rows[i]
                load += demand
                soft += abs(nx - x) + abs(ny - y)
                x, y = nx, ny
            soft += abs(x - dx) + abs(y - dy)
            hard += max(0, load - cap)
        return hard, soft
    return _factory_score(rows, lists)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="重新生成并覆盖四个算例")
    write = parser.parse_args().write
    for kind in CASES:
        data, values, lists = generate(kind)
        path = HERE / "instances" / (kind + ".txt")
        if write:
            path.write_bytes(data)
        if data != path.read_bytes():
            raise ValueError("冻结算例与生成器不符：" + str(path))
        hard, soft = check(path, values, lists)
        if hard:
            raise ValueError("生成器的见证违反硬约束")
        print(kind, "sha256", hashlib.sha256(data).hexdigest(), "见证分数", hard, soft)


if __name__ == "__main__":
    main()
