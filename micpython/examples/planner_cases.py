"""四种业务问题的建模示例。运行 python micpython/examples/planner_cases.py。"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from micpython.planner import Model


def solve_and_print(model, case, variables=(), group=None):
    result = model.optimize(seed=7, steps=120, threads=1, algo="tabu", construct="bf")
    values = ",".join(str(result.value(v)) for v in variables)
    routes = "" if group is None else "/".join(
        ",".join(map(str, result.sequence(group, owner))) for owner in range(group.owners)
    )
    penalties = ",".join(map(str, result.penalties.values()))
    print(case, result.hard, result.soft, result.steps, result.status, result.reason,
          int(result.complete), values, routes, penalties, sep="\t")


def roster():
    with Model("排班") as model:
        shifts = [model.addVar([0, 1]), model.addVar([0, 1]), model.addVar([1])]
        for employee in range(2):
            count = model.count(shift == employee for shift in shifts)
            model.addConstr(count <= 2, "capacity%d" % employee)
        for i, shift in enumerate(shifts):
            model.addConstr(shift == 1, "pref%d" % i, weight=2, hard=False)
        solve_and_print(model, "roster", shifts)


def timetable():
    with Model("排课") as model:
        lessons = [(model.addVar([0, 1]), model.addVar([0, 1])) for _ in range(3)]
        times, rooms = zip(*lessons)
        model.addConstr(times[0] != times[1], "teacher")
        for i in range(len(lessons)):
            for j in range(i + 1, len(lessons)):
                conflict = (times[i] == times[j]) * (rooms[i] == rooms[j])
                model.addPenalty(conflict, "room%d%d" % (i, j), hard=True)
        model.setObjective(model.sum(times))
        solve_and_print(model, "timetable", [v for lesson in lessons for v in lesson])


def delivery():
    with Model("配送") as model:
        trucks = model.addListVar(2, [[0, 1], [0, 1], [0, 1]])
        # 三位客户在直线上的位置分别为 1、3、8，最后一个矩阵下标对应仓库。
        locations = [1, 3, 8, 0]
        matrix = [[abs(a - b) for b in locations] for a in locations]
        for vehicle in range(2):
            model.addConstr(trucks.load(vehicle, [1, 1, 1]) <= 2, "capacity%d" % vehicle)
        model.setObjective(model.sum(trucks.distance(v, matrix) for v in range(2)))
        solve_and_print(model, "delivery", group=trucks)


def factory():
    with Model("生产排程") as model:
        machines = model.addListVar(2, [[0], [1], [1], [0]], durations=[2, 1, 2, 1])
        machines.precede(0, 1)
        machines.precede(2, 3)
        model.setObjective(model.max(machines.finish(i) for i in range(4)))
        solve_and_print(model, "factory", group=machines)


if __name__ == "__main__":
    roster()
    timetable()
    delivery()
    factory()
