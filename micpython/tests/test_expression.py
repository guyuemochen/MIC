# -*- coding: utf-8 -*-
"""LinExpr / quicksum / TempConstr 构建语义测试（全部经 micpython 公共 API）。

数值容差 1e-6（与 milp/tests/acceptance.rs 同款）。
"""

import unittest

from micpython import (
    ArgumentError,
    LinExpr,
    MAXIMIZE,
    Model,
    Status,
    quicksum,
)

TOL = 1e-6


class ExpressionTests(unittest.TestCase):

    def setUp(self):
        self.m = Model('expr-tests')

    def tearDown(self):
        self.m.dispose()

    def test_temp_constr_terms_rhs_direction(self):
        x = self.m.addVar(name='x')
        y = self.m.addVar(name='y')
        c = x + 2 * y <= 5
        self.assertEqual(c._terms, {x._index: 1.0, y._index: 2.0})
        self.assertEqual(c._rhs, 5.0)
        self.assertEqual(c._op, 1, '方向编码：1=LE')
        self.assertEqual((x + 2 * y >= 5)._op, 2, '方向编码：2=GE')
        self.assertEqual((x + 2 * y == 5)._op, 3, '方向编码：3=EQ')
        # 双端表达式比较：移项合并常数到右端（x <= y + 2 ⇒ x - y <= 2）。
        moved = x <= y + 2
        self.assertEqual(moved._terms, {x._index: 1.0, y._index: -1.0})
        self.assertEqual(moved._rhs, 2.0)
        self.assertEqual(moved._op, 1)

    def test_quicksum_empty_sequence(self):
        empty = quicksum([])
        self.assertIsInstance(empty, LinExpr)
        self.assertEqual(len(empty), 0)
        self.assertEqual(empty.getConstant(), 0.0)
        self.assertEqual(empty.getValue(), 0.0, '空表达式在任意解上取值 0')

    def test_quicksum_mixed_terms_and_constant(self):
        x = self.m.addVar(name='x')
        y = self.m.addVar(name='y')
        e = quicksum([1, x, 2 * y, 3])
        self.assertEqual(e._terms, {x._index: 1.0, y._index: 2.0})
        self.assertEqual(e.getConstant(), 4.0, '数值项 1 与 3 累入常数')
        self.assertEqual(len(e), 2)

    def test_quicksum_negative_coefficient_merges(self):
        x = self.m.addVar(name='x')
        merged = quicksum([x, x, -x])
        self.assertEqual(merged._terms, {x._index: 1.0}, 'x + x - x = x')
        vanished = quicksum([2 * x, -2 * x])
        self.assertEqual(len(vanished), 0, '归并为零系数的项被丢弃')
        self.assertEqual(vanished.getConstant(), 0.0)

    def test_scalar_multiplication_commutative(self):
        x = self.m.addVar(name='x')
        left = 2 * x
        right = x * 2
        self.assertIsInstance(left, LinExpr)
        self.assertEqual(left._terms, right._terms)
        self.assertEqual(left.getConstant(), right.getConstant())
        self.assertEqual(left._terms, {x._index: 2.0})
        zero = 0 * x
        self.assertEqual(len(zero), 0, '零系数项丢弃（milp LinearExpr 不变量）')

    def test_expression_constant_enters_objective_offset(self):
        # 移植 milp/src/solver.rs:279 用例：
        # max 3x+2y(+100), x+y<=4, x+3y<=6, [0,3]^2 → (3,1)，目标 11+100。
        x = self.m.addVar(lb=0.0, ub=3.0, name='x')
        y = self.m.addVar(lb=0.0, ub=3.0, name='y')
        self.m.setObjective(3 * x + 2 * y + 100, MAXIMIZE)
        self.m.addConstr(x + y <= 4)
        self.m.addConstr(x + 3 * y <= 6)
        self.assertEqual(self.m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(self.m.ObjVal, 111.0, delta=TOL,
                               msg='常数 100 应作为目标偏移透传')
        self.assertAlmostEqual(x.X, 3.0, delta=TOL)
        self.assertAlmostEqual(y.X, 1.0, delta=TOL)
        self.assertAlmostEqual((3 * x + 2 * y + 100).getValue(), 111.0,
                               delta=TOL, msg='getValue 按最近一次解求值')

    def test_cross_model_mixing_raises_argument_error(self):
        other = Model('expr-other')
        self.addCleanup(other.dispose)
        a = self.m.addVar(name='a')
        b = other.addVar(name='a')
        with self.assertRaises(ArgumentError):
            a + b
        with self.assertRaises(ArgumentError):
            a - b
        with self.assertRaises(ArgumentError):
            a <= b
        with self.assertRaises(ArgumentError):
            quicksum([a, b])
        with self.assertRaises(ArgumentError):
            self.m.addConstr(b <= 1, '引用其他 Model 的变量')
        with self.assertRaises(ArgumentError):
            self.m.setObjective(b)


if __name__ == '__main__':
    unittest.main()
