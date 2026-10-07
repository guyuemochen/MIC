# -*- coding: utf-8 -*-
"""Model 建模语义测试：计数精确性、命名规则、界校验、目标方向（公共 API）。

- 重名 / 撞保留字 / lb>ub / Binary 界越界均在 FFI→milp 层拒绝
  （model.rs:203-247）→ micpython 映射 ModelError；
- lb 非有限（-∞）在 FFI 层拒绝（码 9）→ ArgumentError。
"""

import unittest

from micpython import (
    BINARY,
    ArgumentError,
    INF,
    MAXIMIZE,
    MINIMIZE,
    Model,
    ModelError,
    Status,
)

# 与 milp/src/model.rs:62-75 的 RESERVED_WORDS 一致（.milp 语法保留字）。
RESERVED_WORDS = (
    'name', 'sense', 'vars', 'obj', 'cons', 'bounds', 'end',
    'continuous', 'integer', 'binary', 'max', 'min',
)


class ModelBuildingTests(unittest.TestCase):

    def setUp(self):
        self.m = Model('model-tests')

    def tearDown(self):
        self.m.dispose()

    def test_num_vars_and_num_constrs_exact(self):
        m = self.m
        self.assertEqual(m.NumVars, 0)
        self.assertEqual(m.NumConstrs, 0)
        v = m.addVar(name='a')
        self.assertEqual(m.NumVars, 1)
        m.addVars(5, vtype=BINARY, name='x')
        self.assertEqual(m.NumVars, 6)
        m.addConstr(v <= 3)
        self.assertEqual(m.NumConstrs, 1)
        m.addConstr(v >= 0)
        self.assertEqual(m.NumConstrs, 2)
        with self.assertRaises(ModelError):
            m.addVar(name='a')  # 重名被拒
        self.assertEqual(m.NumVars, 6, '被拒绝的 addVar 不留残余列')

    def test_auto_names_unique_and_reserved_free(self):
        names = [self.m.addVar().VarName for _ in range(20)]
        self.assertEqual(len(set(names)), 20, '自动命名必须唯一')
        for name in names:
            self.assertNotIn(name, RESERVED_WORDS,
                             '自动命名必须避开 .milp 保留字')
        # 显式占用首个自动名后，自动命名跳过该名字（不自作聪明去重）。
        m2 = Model('auto-skip')
        self.addCleanup(m2.dispose)
        m2.addVar(name='x0')
        self.assertEqual(m2.addVar().VarName, 'x1')

    def test_duplicate_explicit_name_raises_model_error(self):
        self.m.addVar(name='dup')
        with self.assertRaises(ModelError) as ctx:
            self.m.addVar(name='dup')
        self.assertEqual(ctx.exception.code, 2, 'FFI 码 2=ERR_INVALID_MODEL')
        self.assertIn('重复', str(ctx.exception))

    def test_reserved_word_name_raises_model_error(self):
        for word in ('name', 'sense', 'end', 'binary', 'max', 'min'):
            with self.assertRaises(ModelError):
                self.m.addVar(name=word)

    def test_lb_greater_than_ub_raises_model_error(self):
        with self.assertRaises(ModelError):
            self.m.addVar(lb=5.0, ub=1.0, name='bad')
        # lb 非有限（-∞）：能力边界，FFI 层码 9 → ArgumentError。
        with self.assertRaises(ArgumentError):
            self.m.addVar(lb=-INF, name='neg-inf')

    def test_binary_bounds_outside_unit_interval_raises_model_error(self):
        with self.assertRaises(ModelError):
            self.m.addVar(vtype=BINARY, ub=2.0, name='b1')
        with self.assertRaises(ModelError):
            self.m.addVar(vtype=BINARY, lb=-0.5, name='b2')

    def test_get_var_by_name_hit_and_miss(self):
        v = self.m.addVar(name='hit')
        self.assertIs(self.m.getVarByName('hit'), v)
        self.assertIsNone(self.m.getVarByName('miss'))

    def test_set_objective_both_senses(self):
        x = self.m.addVar(lb=1.0, ub=3.0, name='x')
        self.m.setObjective(x, MINIMIZE)
        self.assertEqual(self.m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(self.m.ObjVal, 1.0, delta=1e-6)
        self.m.setObjective(x, MAXIMIZE)
        self.assertEqual(self.m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(self.m.ObjVal, 3.0, delta=1e-6)
        # sense=None 保持现方向（上一次显式给出的 MAXIMIZE）。
        self.m.setObjective(2 * x)
        self.assertEqual(self.m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(self.m.ObjVal, 6.0, delta=1e-6)

    def test_addvars_same_prefix_twice_second_raises_model_error(self):
        first = self.m.addVars(3, name='x')
        self.assertEqual([v.VarName for v in first], ['x0', 'x1', 'x2'])
        with self.assertRaises(ModelError):
            self.m.addVars(2, name='x')
        self.assertEqual(self.m.NumVars, 3, '第二次 addVars 整体失败不留列')


if __name__ == '__main__':
    unittest.main()
