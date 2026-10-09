# -*- coding: utf-8 -*-
"""参数语义测试：镜像三铁律、关断编码、gap 向量、非法值（公共 API）。

gap 统一向量 hard-knapsack-12（前轮实跑钉死：CLI ``--mipgap 20`` →
``status=gap_reached, nodes=1, gap=3.11%, best=511.000000``；``--mipgap 0``
→ ``status=optimal, objective=523.000000``。2026-10 割族增强后根 LP 终
顶点漂移：旧值 512/2.83% 出自旧 GMI 终顶点的取整复核；新割终顶点中面
化，根 incumbent 改由割前顶点 dive 供应 = 511）。
"""

import unittest

from micpython import (
    BINARY,
    ArgumentError,
    INF,
    MAXIMIZE,
    Model,
    NoSolutionError,
    Status,
    UnknownParameterError,
    quicksum,
)

TOL = 1e-6
# 时间相关断言留足松弛（求解器内部 time_ms 受 5s 上限约束，另有调度开销）。
TIME_SLACK_SECS = 10.0


def build_knapsack_10():
    """等价 milp/instances/knapsack-10.milp（最优 17）。"""
    m = Model('knapsack-10')
    x = m.addVars(5, vtype=BINARY, name='x')
    values = [7, 5, 9, 6, 3]
    weights = [4, 3, 5, 6, 2]
    m.setObjective(quicksum(values[i] * x[i] for i in range(5)), MAXIMIZE)
    m.addConstr(quicksum(weights[i] * x[i] for i in range(5)) <= 10)
    return m


def build_hard_knapsack_12():
    """等价 milp/instances/hard-knapsack-12.milp（最优 523）。"""
    m = Model('hard-knapsack-12')
    x = m.addVars(12, vtype=BINARY, name='x')
    values = [11, 20, 32, 44, 55, 70, 88, 99, 111, 122, 133, 144]
    weights = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    m.setObjective(quicksum(values[i] * x[i] for i in range(12)), MAXIMIZE)
    m.addConstr(quicksum(weights[i] * x[i] for i in range(12)) <= 172)
    return m


class ParamMirrorTests(unittest.TestCase):
    """镜像三铁律：恒调 FFI、仅成功后更新、getParam 返回原值。"""

    def test_initial_mirror_defaults(self):
        m = Model('param-defaults')
        self.addCleanup(m.dispose)
        self.assertEqual(m.getParam('TimeLimit'), INF)
        self.assertEqual(m.getParam('Threads'), 1)
        self.assertEqual(m.getParam('MIPGap'), 0.0)
        self.assertEqual(m.getParam('HeuristicMode'), 'on')
        self.assertEqual(m.getParam('LogPeriodSecs'), 0.0)

    def test_mipgap_exact_readback(self):
        m = Model('param-readback')
        self.addCleanup(m.dispose)
        m.setParam('MIPGap', 0.07)
        # 精确相等：镜像保存原始 Python 值，不做 ÷100 逆运算
        # （0.07×100=7.000000000000001 的浮点噪声不进入回读）。
        self.assertEqual(m.getParam('MIPGap'), 0.07)

    def test_rejected_value_does_not_pollute_mirror(self):
        m = Model('param-reject')
        self.addCleanup(m.dispose)
        with self.assertRaises(ArgumentError):
            m.setParam('Threads', 0)
        self.assertEqual(m.getParam('Threads'), 1, '默认值未被被拒值污染')
        m.setParam('Threads', 4)
        with self.assertRaises(ArgumentError):
            m.setParam('Threads', 0)
        self.assertEqual(m.getParam('Threads'), 4, '前次成功值保持')

    def test_shutdown_encoding_time_limit_inf(self):
        # TimeLimit 5.0 → INF：FFI 收到 0.0（不限），真正解除限速，
        # 无『FFI 残留 Some(5.0) 照旧限速』的失真。
        for builder, expected in ((build_knapsack_10, 17.0),
                                  (build_hard_knapsack_12, 523.0)):
            m = builder()
            self.addCleanup(m.dispose)
            m.setParam('TimeLimit', 5.0)
            self.assertEqual(m.getParam('TimeLimit'), 5.0)
            m.setParam('TimeLimit', INF)
            self.assertEqual(m.getParam('TimeLimit'), INF)
            self.assertEqual(m.optimize(), Status.OPTIMAL,
                             '解除限速后应完整最优（%s）' % m.ModelName)
            self.assertAlmostEqual(m.ObjVal, expected, delta=TOL)

    def test_shutdown_encoding_log_period_secs(self):
        # 2.5 → 0.0：FFI 收到 0.0 真正关断进度钩子（无 stdout 断言依赖）。
        m = build_knapsack_10()
        self.addCleanup(m.dispose)
        m.setParam('LogPeriodSecs', 2.5)
        self.assertEqual(m.getParam('LogPeriodSecs'), 2.5)
        m.setParam('LogPeriodSecs', 0.0)
        self.assertEqual(m.getParam('LogPeriodSecs'), 0.0)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(m.ObjVal, 17.0, delta=TOL)


class ParamBehaviorTests(unittest.TestCase):
    """参数经 FFI 透传后的求解行为。"""

    def test_time_limit_knapsack_still_optimal(self):
        m = build_knapsack_10()
        self.addCleanup(m.dispose)
        m.setParam('TimeLimit', 5.0)
        self.assertEqual(m.getParam('TimeLimit'), 5.0, '精确回读')
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(m.ObjVal, 17.0, delta=TOL)
        # stat time_ms ≤ 5000 + 松弛（Runtime = time_ms/1000）。
        self.assertGreaterEqual(m.Runtime, 0.0)
        self.assertLessEqual(m.Runtime, 5.0 + TIME_SLACK_SECS)

    def test_threads_parallel_hard_knapsack(self):
        m = build_hard_knapsack_12()
        self.addCleanup(m.dispose)
        m.setParam('Threads', 4)
        self.assertEqual(m.getParam('Threads'), 4)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        # 并行统计量可变（纪律同 acceptance.rs）：仅断言终态与目标。
        self.assertAlmostEqual(m.ObjVal, 523.0, delta=TOL)

    def test_mipgap_triggers_limit_reached(self):
        m = build_hard_knapsack_12()
        self.addCleanup(m.dispose)
        m.setParam('MIPGap', 0.20)  # 小数 → FFI 20.0 百分比
        self.assertEqual(m.getParam('MIPGap'), 0.20)
        self.assertEqual(m.optimize(), Status.LIMIT_REACHED,
                         '达阈值提前终止的诚实终态（非 OPTIMAL）')
        self.assertTrue(m.GapTerminated)
        self.assertAlmostEqual(m.ObjVal, 511.0, delta=TOL,
                               msg='LIMIT_REACHED 携带 incumbent')
        gap = m.MIPGapPercent
        self.assertGreater(gap, 2.0)
        self.assertLess(gap, 4.0, '实跑 CLI 佐证 3.11%')
        self.assertEqual(m.SolCount, 1)

    def test_mipgap_zero_still_optimal(self):
        m = build_hard_knapsack_12()
        self.addCleanup(m.dispose)
        m.setParam('MIPGap', 0.0)
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(m.ObjVal, 523.0, delta=TOL)
        self.assertFalse(m.GapTerminated)
        with self.assertRaises(NoSolutionError):
            m.MIPGapPercent  # 未因 gap 终止 → 无 gap 证据可读

    def test_heuristic_mode_off_knapsack_still_optimal(self):
        m = build_knapsack_10()
        self.addCleanup(m.dispose)
        m.setParam('HeuristicMode', 'off')
        self.assertEqual(m.getParam('HeuristicMode'), 'off')
        self.assertEqual(m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(m.ObjVal, 17.0, delta=TOL)


class ParamRejectionTests(unittest.TestCase):
    """非法值逐一断言异常类型（FFI 码 9 → ArgumentError；未知名 → 10）。"""

    def setUp(self):
        self.m = Model('param-invalid')
        self.addCleanup(self.m.dispose)

    def test_threads_zero_raises_argument_error(self):
        with self.assertRaises(ArgumentError):
            self.m.setParam('Threads', 0)

    def test_threads_above_1024_raises_argument_error(self):
        with self.assertRaises(ArgumentError):
            self.m.setParam('Threads', 1025)

    def test_mipgap_negative_raises_argument_error(self):
        with self.assertRaises(ArgumentError):
            self.m.setParam('MIPGap', -1)

    def test_time_limit_negative_raises_argument_error(self):
        with self.assertRaises(ArgumentError):
            self.m.setParam('TimeLimit', -2)

    def test_heuristic_mode_invalid_string_raises_argument_error(self):
        with self.assertRaises(ArgumentError):
            self.m.setParam('HeuristicMode', 'fast')

    def test_unknown_parameter_name_raises_unknown_parameter_error(self):
        with self.assertRaises(UnknownParameterError):
            self.m.setParam('Foo', 1)


if __name__ == '__main__':
    unittest.main()
