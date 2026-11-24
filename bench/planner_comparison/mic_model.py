"""从共同算例构造 MIC 的建模协议，文件解析和建模在计时区外完成。"""
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from micpython.planner import Model


def roster(m, header, rows):
    employees, limit = header
    variables = []
    costs = []
    slots = defaultdict(list)
    for slot, *cost in rows:
        eligible = [e for e in range(employees) if cost[e] >= 0]
        x = m.addVar(eligible)
        slots[slot].append(x)
        variables.append(x)
        costs.append(m.lookup(x, cost))
    zero = m.constant(0)
    for e in range(employees):
        total = m.sum(x == e for x in variables)
        m.addPenalty(m.max([zero, total - limit]), "workload%d" % e, hard=True)
    collisions = [x == y for group in slots.values() for x, y in combinations(group, 2)]
    m.addPenalty(m.sum(collisions), "overlap", hard=True)
    m.setObjective(m.sum(costs))


def timetable(m, header, rows):
    slots, rooms = header
    times, room_times, costs = [], [], []
    for _, _, preferred in rows:
        time = m.addVar(range(slots))
        room = m.addVar(range(rooms))
        times.append(time)
        room_times.append(time * rooms + room)
        costs.append(m.lookup(time, [abs(t - preferred) for t in range(slots)]))
    room_clashes, teacher_clashes, class_clashes = [], [], []
    for i, j in combinations(range(len(rows)), 2):
        room_clashes.append(room_times[i] == room_times[j])
        if rows[i][0] == rows[j][0]:
            teacher_clashes.append(times[i] == times[j])
        if rows[i][1] == rows[j][1]:
            class_clashes.append(times[i] == times[j])
    for name, clashes in [("room", room_clashes), ("teacher", teacher_clashes),
                          ("cohort", class_clashes)]:
        m.addPenalty(m.sum(clashes), name, hard=True)
    m.setObjective(m.sum(costs))


def delivery(m, header, rows):
    vehicles, capacity, depot_x, depot_y = header
    routes = m.addListVar(vehicles, [range(vehicles) for _ in rows])
    points = [(x, y) for x, y, _ in rows] + [(depot_x, depot_y)]
    matrix = [[abs(x - u) + abs(y - v) for u, v in points] for x, y in points]
    costs = []
    for vehicle in range(vehicles):
        load = routes.load(vehicle, [d for _, _, d in rows])
        m.addPenalty(m.max([0, load - capacity]), "load%d" % vehicle, hard=True)
        costs.append(routes.distance(vehicle, matrix))
    m.setObjective(m.sum(costs))


def factory(m, header, rows):
    machines, = header
    schedule = m.addListVar(machines, [[machine] for machine, _, _ in rows],
                            [duration for _, duration, _ in rows])
    for item, (_, _, prev) in enumerate(rows):
        if prev >= 0:
            schedule.precede(prev, item)
    m.setObjective(m.max(schedule.finish(i) for i in range(len(rows))))


def compile_model(path):
    with path.open(encoding="utf8") as f:
        kind, *header = f.readline().split()
        rows = [list(map(int, line.split()[1:])) for line in f if line.strip()]
    with Model(path.stem) as m:
        {"ROSTER": roster, "TIMETABLE": timetable, "DELIVERY": delivery,
         "FACTORY": factory}[kind](m, list(map(int, header)), rows)
        # 交给原生计时入口，避免把 SDK 序列化和跨语言往返算入求解时间。
        return "\n".join(m._lines + ["O %d" % m._objective]) + "\n"
