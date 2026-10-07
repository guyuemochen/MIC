# -*- coding: utf-8 -*-
"""求解终态与结果读取测试（公共 API；数值容差 1e-6，与 acceptance.rs 同款）。

重建向量均取自 milp/instances/ 对应 .milp 文件：
- knapsack-10 → 17（x2/x3/x5=1，即 0 基下标 1/2/4）
- assign-4x4 → 10；milp-wolsey-2 → 40；milp-mixed-5 → 29；lp-prod-3 → 13
- LP+offset（solver.rs:279 用例）→ 111.0
- dive 模型（solver.rs:481 的 dive-gap 问题，FFI 默认选项下根割闭合）→ 21
"""

import unittest

from micpython import (
    BINARY,
    INF,
    INTEGER,
    MAXIMIZE,
    MINIMIZE,
    MicError,
    Model,
    NoSolutionError,
    Status,
    quicksum,
)

TOL = 1e-6


def build_knapsack_10():
    """等价 milp/instances/knapsack-10.milp。"""
    m = Model('knapsack-10')
    x = m.addVars(5, vtype=BINARY, name='x')
    values = [7, 5, 9, 6, 3]
    weights = [4, 3, 5, 6, 2]
    m.setObjective(quicksum(values[i] * x[i] for i in range(5)), MAXIMIZE)
    m.addConstr(quicksum(weights[i] * x[i] for i in range(5)) <= 10)
    return m, x, values


def build_assign_4x4():
    """等价 milp/instances/assign-4x4.milp。"""
    cost = [[9, 8, 1, 9],
            [2, 8, 8, 9],
            [9, 8, 9, 3],
            [8, 4, 8, 9]]
    m = Model('assign-4x4')
    x = {}
    for i in range(1, 5):
        for j in range(1, 5):
            x[i, j] = m.addVar(vtype=BINARY, name='x%d%d' % (i, j))
    m.setObjective(quicksum(cost[i - 1][j - 1] * x[i, j]
                            for i in range(1, 5) for j in range(1, 5)),
                   MINIMIZE)
    for i in range(1, 5):
        m.addConstr(quicksum(x[i, j] for j in range(1, 5)) == 1)
    for j in range(1, 5):
        m.addConstr(quicksum(x[i, j] for i in range(1, 5)) == 1)
    return m, x


def build_lp_offset():
    """solver.rs:279 用例：max 3x+2y+100，x+y<=4，x+3y<=6，[0,3]^2 → 111。"""
    m = Model('lp-offset')
    x = m.addVar(lb=0.0, ub=3.0, name='x')
    y = m.addVar(lb=0.0, ub=3.0, name='y')
    m.setObjective(3 * x + 2 * y + 100, MAXIMIZE)
    m.addConstr(x + y <= 4)
    m.addConstr(x + 3 * y <= 6)
    return m, x, y


def build_dive_model():
    """solver.rs:481 的 dive-gap 问题：max 10x1+x2，20x1+3x2<=43，整数 → 21。"""
    m = Model('dive')
    x1 = m.addVar(vtype=INTEGER, ub=10.0, name='x1')
    x2 = m.addVar(vtype=INTEGER, ub=10.0, name='x2')
    m.setObjective(10 * x1 + x2, MAXIMIZE)
    m.addConstr(20 * x1 + 3 * x2 <= 43)
    return m, x1, x2


def build_wolsey_2():
    """等价 milp/instances/milp-wolsey-2.milp。"""
    m = Model('milp-wolsey-2')
    x1 = m.addVar(vtype=INTEGER, name='x1')
    x2 = m.addVar(vtype=INTEGER, name='x2')
    m.setObjective(5 * x1 + 8 * x2, MAXIMIZE)
    m.addConstr(x1 + x2 <= 6)
    m.addConstr(5 * x1 + 9 * x2 <= 45)
    return m


def build_mixed_5():
    """等价 milp/instances/milp-mixed-5.milp。"""
    m = Model('milp-mixed-5')
    y1 = m.addVar(vtype=BINARY, name='y1')
    y2 = m.addVar(vtype=BINARY, name='y2')
    z = m.addVar(vtype=INTEGER, ub=10.0, name='z')
    c = m.addVar(vtype='C', ub=10.0, name='c')
    m.setObjective(10 * y1 + 8 * y2 + 4 * z + 1 * c, MAXIMIZE)
    m.addConstr(6 * y1 + 5 * y2 + 3 * z + c <= 20)
    m.addConstr(y1 + y2 + z <= 4)
    m.addConstr(c - 6 * y1 <= 0)
    m.addConstr(z - 2 * y2 >= 0)
    return m


def build_lp_prod_3():
    """等价 milp/instances/lp-prod-3.milp。"""
    m = Model('lp-prod-3')
    x = m.addVars(3, name='x')  # 连续变量，默认 [0, +inf)
    m.setObjective(5 * x[0] + 4 * x[1] + 3 * x[2], MAXIMIZE)
    m.addConstr(2 * x[0] + 3 * x[1] + x[2] <= 5)
    m.addConstr(4 * x[0] + x[1] + 2 * x[2] <= 11)
    m.addConstr(3 * x[0] + 4 * x[1] + 2 * x[2] <= 8)
    return m, x


class SolveStatusTests(unittest.TestCase):

    def test_knapsack_10_optimal(self):
        m, x, _ = build_knapsack_10()
        self.addCleanup(m.dispose)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(m.ObjVal, 17.0, delta=TOL)
        for i in (1, 2, 4):
            self.assertAlmostEqual(x[i].X, 1.0, delta=TOL,
                                   msg='x%d 应入选' % i)
        for i in (0, 3):
            self.assertAlmostEqual(x[i].X, 0.0, delta=TOL,
                                   msg='x%d 不应入选' % i)
        self.assertGreaterEqual(m.NodeCount, 1, '根节点计 1')
        self.assertGreaterEqual(m.Runtime, 0.0)
        self.assertEqual(m.SolCount, 1)

    def test_assign_4x4_optimal(self):
        m, x = build_assign_4x4()
        self.addCleanup(m.dispose)
        self.assertEqual(m.NumVars, 16)
        self.assertEqual(m.NumConstrs, 8)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(m.ObjVal, 10.0, delta=TOL)
        self.assertAlmostEqual(sum(v.X for v in x.values()), 4.0, delta=TOL,
                               msg='指派矩阵每行每列恰一个 1')

    def test_lp_with_offset_111(self):
        m, x, y = build_lp_offset()
        self.addCleanup(m.dispose)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(m.ObjVal, 111.0, delta=TOL)
        self.assertAlmostEqual(x.X, 3.0, delta=TOL)
        self.assertAlmostEqual(y.X, 1.0, delta=TOL)

    def test_infeasible_status_and_x_raises(self):
        m = Model('infeasible')
        self.addCleanup(m.dispose)
        x = m.addVar(name='x')  # 连续，[0, +inf)
        m.addConstr(x >= 3)
        m.addConstr(x <= 1)
        self.assertEqual(m.optimize(), Status.INFEASIBLE)
        self.assertEqual(m.Status, Status.INFEASIBLE)
        with self.assertRaises(NoSolutionError) as ctx:
            x.X
        self.assertIn("Unable to retrieve attribute 'X'",
                      str(ctx.exception))
        with self.assertRaises(NoSolutionError):
            m.ObjVal
        self.assertEqual(m.SolCount, 0)

    def test_unbounded_status(self):
        m = Model('unbounded', MAXIMIZE)
        self.addCleanup(m.dispose)
        x = m.addVar(name='x')  # 连续无上界
        m.setObjective(x)
        self.assertEqual(m.optimize(), Status.UNBOUNDED)
        with self.assertRaises(NoSolutionError):
            m.ObjVal
        self.assertEqual(m.SolCount, 0)

    def test_dive_model_default_params_optimal_21(self):
        # 该模型在 FFI 默认选项（presolve+根割开）下根节点闭合整数 hull，
        # gap 永不触发（默认 MIPGap=0 即关闭）——正向佐证 objective=21。
        m, x1, x2 = build_dive_model()
        self.addCleanup(m.dispose)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(m.ObjVal, 21.0, delta=TOL)
        self.assertFalse(m.GapTerminated)

    def test_wolsey_2_optimal_40(self):
        m = build_wolsey_2()
        self.addCleanup(m.dispose)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(m.ObjVal, 40.0, delta=TOL)

    def test_mixed_5_optimal_29(self):
        m = build_mixed_5()
        self.addCleanup(m.dispose)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(m.ObjVal, 29.0, delta=TOL)

    def test_lp_prod_3_optimal_13(self):
        m, _ = build_lp_prod_3()
        self.addCleanup(m.dispose)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(m.ObjVal, 13.0, delta=TOL)

    def test_objval_matches_quicksum_recompute(self):
        m, x, values = build_knapsack_10()
        self.addCleanup(m.dispose)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        recomputed = quicksum(values[i] * x[i] for i in range(5)).getValue()
        self.assertAlmostEqual(recomputed, m.ObjVal, delta=TOL,
                               msg='ObjVal 与按解重算的目标一致')


class LifecycleTests(unittest.TestCase):

    def test_dispose_guard_and_idempotence(self):
        m = Model('dispose-guard')
        v = m.addVar(vtype=BINARY, name='b')
        m.setObjective(v, MAXIMIZE)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        m.dispose()
        # dispose 后一切经 _call 的使用抛 MicError（Python 侧守卫，非段错误）。
        with self.assertRaises(MicError):
            m.optimize()
        with self.assertRaises(MicError):
            m.addVar(name='c')
        with self.assertRaises(MicError):
            m.ObjVal
        with self.assertRaises(MicError):
            v.X
        with self.assertRaises(MicError):
            m.setParam('Threads', 2)
        with self.assertRaises(MicError):
            m.NumVars
        # 二次 dispose：先于 _call 短路，静默返回不抛。
        m.dispose()
        # getParam 读取本地参数镜像：不触碰 FFI（setParam 已被守卫拒绝，
        # 镜像保持默认值），因此 dispose 后仍可安全读取。
        self.assertEqual(m.getParam('Threads'), 1)
        self.assertEqual(m.getParam('HeuristicMode'), 'on')
        # del 已 dispose 的模型无异常（__del__ 幂等路径）。
        m2 = Model('dispose-del')
        m2.addVar(name='z')
        m2.dispose()
        del m2
        # 上下文管理器：退出即 dispose，此后使用同样被守卫拦截。
        with Model('dispose-ctx') as m3:
            m3.addVar(name='w')
        with self.assertRaises(MicError):
            m3.NumVars

    def test_x_reflects_latest_optimize_and_dimension_anchor(self):
        m = Model('reopt')
        self.addCleanup(m.dispose)
        x = m.addVar(lb=0.0, ub=1.0, name='x')
        m.setObjective(x, MAXIMIZE)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(x.X, 1.0, delta=TOL)
        # 结构性修改后可重复 solve；X/ObjVal 反映最近一次 optimize。
        m.addConstr(x <= 0.5)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(x.X, 0.5, delta=TOL)
        # solve 后新增变量：解维数锚定 n_solved，新变量取值 → NoSolutionError。
        y = m.addVar(lb=0.0, ub=1.0, name='y')
        with self.assertRaises(NoSolutionError):
            y.X
        self.assertAlmostEqual(x.X, 0.5, delta=TOL, msg='旧变量仍可读')
        # 再次 optimize 后新变量进入解向量。
        m.setObjective(x + y, MAXIMIZE)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(y.X, 1.0, delta=TOL)
        self.assertAlmostEqual(m.ObjVal, 1.5, delta=TOL)


if __name__ == '__main__':
    unittest.main()
