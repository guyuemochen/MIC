"""使用真实动态库验证建模、重复求解及错误处理。"""
import unittest

from micpython.planner import Model


class PlannerTests(unittest.TestCase):
    def test_solve_and_reuse(self):
        with Model() as model:
            a, b = model.addVar([0, 1]), model.addVar([0, 1])
            model.addConstr(a != b, "同一时段不可重合")
            model.addPenalty(model.count([a == 0, b == 0]), "零班", weight=2)
            model.setObjective(model.lookup(a, [5, 0]))
            first = model.optimize(steps=20, construct="bf")
            self.assertEqual((first.value(a), first.value(b), first.hard, first.soft), (1, 0, 0, 2))
            self.assertEqual(first.penalties, {"零班": 2, "同一时段不可重合": 0})
            self.assertEqual((first.reason, first.status), ("Steps", "LimitReached"))
            model.setObjective(a * 5)
            second = model.optimize(steps=20, construct="bf")
            self.assertEqual((second.value(a), second.value(b)), (0, 1))
        self.assertTrue(first.feasible)
        self.assertEqual(first.value(a), 1)
        model.close()
        with self.assertRaises(RuntimeError):
            model.optimize()

    def test_list_and_timing(self):
        with Model() as m:
            group = m.addListVar(2, [[0], [1]], [2, 3])
            group.precede(1, 0)
            m.addConstr(group.load(1, [1, 3]) == 3, "载重")
            m.setObjective(m.max([group.finish(0), group.finish(1)]) +
                           group.distance(0, [[0, 2, 1], [2, 0, 8], [3, 9, 0]]))
            result = m.optimize(steps=10)
            self.assertTrue(result.feasible)
            self.assertEqual((result.soft, result.sequence(group, 0), result.sequence(group, 1)),
                             (9, (0,), (1,)))

    def test_group_aggregation(self):
        with Model() as m:
            keys = [m.addVar([0]), m.addVar([1]), m.addVar([0])]
            m.addConstr(m.groupCount(keys, 0) == 2, "人数")
            m.setObjective(m.groupSum(keys, [7, 11, 13], 0))
            result = m.optimize(steps=0)
            self.assertEqual((result.hard, result.soft), (0, 20))

    def test_directed_adjacency(self):
        with Model() as m:
            x = m.addVar([2])
            route = m.addListVar(2, [[0], [0]])
            route.precede(0, 1)
            forward = route.adjacent(0, 0, 1)
            reverse = route.adjacent(0, 1, 0)
            empty = route.adjacent(1, 0, 1)
            m.addConstr(forward, "顺序")
            m.addConstr(reverse, "逆序", weight=3, hard=False)
            m.setObjective(m.choose(empty, x * 4, x) + m.choose(forward, x, 0))
            with self.assertRaises(ValueError):
                route.adjacent(0, 0, 2)
            result = m.optimize(steps=10)
            self.assertEqual((result.sequence(route, 0), result.sequence(route, 1)), ((0, 1), ()))
            self.assertEqual((result.hard, result.soft, result.penalties),
                             (0, 7, {"顺序": 0, "逆序": 3}))

    def test_native_error_is_sticky(self):
        with Model() as m:
            m.addVar([])
            with self.assertRaises(ValueError) as first:
                m.optimize(steps=0)
            with self.assertRaises(ValueError) as second:
                m.optimize(steps=0)
            self.assertEqual(str(first.exception), str(second.exception))

    def test_invalid_builders(self):
        with Model() as a, Model() as b:
            x, y = a.addVar([0]), b.addVar([1])
            with self.assertRaises(ValueError):
                a.addConstr(x == y, "混合")
            with self.assertRaises(TypeError):
                a.addConstr(x, "整数")
            with self.assertRaises(TypeError):
                bool(x <= 3)
            a.addConstr(x < 1, "重复规则")
            with self.assertRaises(ValueError):
                a.addConstr(x < 2, "重复规则")


if __name__ == "__main__":
    unittest.main()
