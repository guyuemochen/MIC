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

    def test_quicksum_inplace_matches_operator_chain(self):
        """S1 回归：就地累加与逐项 ``total + term`` 算子链逐位等价。

        旧实现（O(n²)）经 LinExpr.__add__ 逐步新建整拷；新实现共享 dict
        就地累加。同序浮点加法 + 零系数终态一致 ⇒ terms/constant/model
        三者必须完全相同。
        """
        xs = self.m.addVars(4, name='q')
        terms = [xs[0], 3, xs[1] * 2, xs[0] * -1.5, xs[2], xs[3] * 0.5, -7,
                 xs[2] * -1, xs[3] * 0.5]
        inplace = quicksum(terms)
        chained = LinExpr()
        for term in terms:
            chained = chained + term
        self.assertEqual(inplace._terms, chained._terms)
        self.assertEqual(inplace.getConstant(), chained.getConstant())
        self.assertIs(inplace._model, chained._model)
        # 中途归零再复活：x - x + x 的终态系数与逐步删零的旧实现一致。
        x = xs[0]
        revived = quicksum([x, -x, x])
        self.assertEqual(revived._terms, {x._index: 1.0})

    def test_quicksum_non_numeric_term_raises_type_error(self):
        # 与 total + term 返回 NotImplemented 后解释器报 TypeError 同款。
        x = self.m.addVar(name='x')
        with self.assertRaises(TypeError):
            quicksum([x, 'not-a-term'])

    def test_quicksum_large_n_linear(self):
        """S1 验收挂钩：5000 项 quicksum 的项数与终值正确（O(n) 化后

        python_bench knapsack-5000 的 build_secs 大头即两条 5000 项
        quicksum；这里只钉语义，计时归微基准/基准脚本。
        """
        xs = self.m.addVars(5000, name='z')
        e = quicksum(1.0 * v for v in xs)
        self.assertEqual(len(e), 5000)
        self.assertEqual(e._terms[0], 1.0)
        self.assertEqual(e._terms[4999], 1.0)

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
