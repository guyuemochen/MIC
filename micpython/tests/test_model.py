# -*- coding: utf-8 -*-
"""Model 建模语义测试：计数精确性、命名规则、界校验、目标方向（公共 API）。

- 重名 / 撞保留字 / lb>ub / Binary 界越界均在 FFI→milp 层拒绝
  （model.rs:203-247）→ micpython 映射 ModelError；
- lb 非有限（-∞）在 FFI 层拒绝（码 9）→ ArgumentError。
"""

import ctypes
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
from micpython import _ffi

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


class AddVarsBatchContractTests(unittest.TestCase):
    """S1 批量契约：addVars 经 milp_ffi_add_vars_batch 单次跨界，语义与
    逐个 addVar 等价；任一元素非法 → ArgumentError/ModelError 且零突变
    （FFI 侧 validate-then-mutate，Python 侧入口预检）。"""

    def setUp(self):
        self.m = Model('addvars-batch')

    def tearDown(self):
        self.m.dispose()

    def test_batch_matches_per_var_addvar_semantics(self):
        # 批量路径与逐个 addVar 双模型对照：名字/界/类型/obj 镜像一致。
        ref = Model('addvars-ref')
        self.addCleanup(ref.dispose)
        got = self.m.addVars(3, lb=0.0, ub=None, vtype=BINARY,
                             obj=[7.0, 5.0, 0.0], name='b')
        want = [ref.addVar(vtype=BINARY, obj=o, name='b%d' % i)
                for i, o in enumerate([7.0, 5.0, 0.0])]
        self.assertEqual([v.VarName for v in got], [v.VarName for v in want])
        self.assertEqual([v.lb for v in got], [0.0, 0.0, 0.0])
        self.assertEqual([v.ub for v in got], [1.0, 1.0, 1.0],
                         'BINARY 的 ub=None → 1 规则不变')
        self.assertEqual([v.VType for v in got], ['B'] * 3)
        self.assertEqual([v.Obj for v in got], [7.0, 5.0, 0.0])
        self.assertEqual(self.m.NumVars, 3)
        self.assertIs(self.m.getVarByName('b2'), got[2])
        # 批量建出的列可解：obj 镜像并入规则不变（表达式为纯常数 ⇒ 目标
        # 恰为镜像 7*b0+5*b1；b0+b1<=1 ⇒ 最优 7）。
        self.m.setObjective(0, MAXIMIZE)
        self.m.addConstr(got[0] + got[1] <= 1)
        self.assertEqual(self.m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(self.m.ObjVal, 7.0, delta=1e-6)
        # 连续变量 ub=None → +INF 默认；obj 标量路径。
        cs = self.m.addVars(2, obj=2.5, name='c')
        self.assertEqual([v.ub for v in cs], [INF, INF])
        self.assertEqual([v.Obj for v in cs], [2.5, 2.5])
        # 下标连续：3 个 b + 2 个 c ⇒ c0._index == 3。
        self.assertEqual([v._index for v in cs], [3, 4])

    def test_batch_invalid_element_zero_mutation(self):
        # vtype 非法（FFI 码 9）。
        with self.assertRaises(ArgumentError):
            self.m.addVars(2, vtype='X', name='x')
        # lb 非有限（FFI 码 9）。
        with self.assertRaises(ArgumentError):
            self.m.addVars(2, lb=-INF, name='x')
        # obj 标量非数值：与逐个 addVar 同款 ArgumentError。
        with self.assertRaises(ArgumentError):
            self.m.addVars(2, obj=None, name='x')
        # obj 序列含非数值元素：入口预检全量拒绝（零突变时点更早）。
        with self.assertRaises(ArgumentError):
            self.m.addVars(2, obj=[1.0, 'bad'], name='x')
        # obj 序列长度不符。
        with self.assertRaises(ArgumentError):
            self.m.addVars(2, obj=[1.0], name='x')
        # 模型语义（lb>ub → ModelError，FFI 入模阶段回滚）。
        with self.assertRaises(ModelError):
            self.m.addVars(2, lb=5.0, ub=1.0, name='x')
        self.assertEqual(self.m.NumVars, 0, '全部非法路径零突变')

    def test_batch_count_zero_returns_empty(self):
        self.assertEqual(self.m.addVars(0), [])
        self.assertEqual(self.m.NumVars, 0)
        # 收紧语义（docstring 已注明）：count=0 也校验标量值域。
        with self.assertRaises(ArgumentError):
            self.m.addVars(0, vtype='X')

    def test_batch_large_count_single_crossing(self):
        # 微基准验收的同款规模：5000 变量一次跨界建列 + 名字/镜像完整。
        vs = self.m.addVars(5000, vtype=BINARY, name='k')
        self.assertEqual(len(vs), 5000)
        self.assertEqual(self.m.NumVars, 5000)
        self.assertEqual(vs[0].VarName, 'k0')
        self.assertEqual(vs[4999].VarName, 'k4999')
        self.assertEqual(vs[4999]._index, 4999)
        self.assertIs(self.m.getVarByName('k2500'), vs[2500])


class SetObjectiveAtomicityTests(unittest.TestCase):
    """B2 回归：setObjective 被 FFI 拒绝时目标方向零突变（FFI 与镜像同步）。"""

    def setUp(self):
        self.m = Model('objective-atomicity', MINIMIZE)

    def tearDown(self):
        self.m.dispose()

    def test_ffi_set_objective_out_of_range_keeps_sense_minimize(self):
        # ctypes 直调 FFI（绕过 core.py）：MINIMIZE 模型上越界 set_objective
        # 携 sense=2(MAXIMIZE) → 码 2；旧实现在校验前已翻 problem.sense。
        lib = _ffi._lib()
        x = self.m.addVar(lb=0.0, ub=1.0, name='x')
        bad_idx = (ctypes.c_uint32 * 1)(7)
        coef = (ctypes.c_double * 1)(1.0)
        rc = lib.milp_ffi_set_objective(self.m._handle, bad_idx, coef, 1, 2, 0.0)
        self.assertEqual(rc, 2, '越界目标项 → FFI 码 2（ERR_INVALID_MODEL）')
        # 合法重设也用 ctypes 直调且 sense=0（保持现状）——不走公共 API
        # （公共 API 会重发显式 sense，掩盖被翻的方向）。
        ok_idx = (ctypes.c_uint32 * 1)(x._index)
        rc2 = lib.milp_ffi_set_objective(self.m._handle, ok_idx, coef, 1, 0, 0.0)
        self.assertEqual(rc2, 0)
        self.assertEqual(self.m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(self.m.ObjVal, 0.0, delta=1e-9,
                               msg='sense 须保持 MINIMIZE（ObjVal=1.0 即已被翻为 MAXIMIZE）')

    def test_rejected_set_objective_keeps_python_sense_mirror(self):
        # 公共 API 路径：inf 系数被 FFI 拒（码 9 → ArgumentError）时，
        # Python 镜像 _sense 也不得先翻（与 FFI problem.sense 同步零突变）。
        x = self.m.addVar(lb=0.0, ub=1.0, name='x')
        with self.assertRaises(ArgumentError):
            self.m.setObjective(INF * x, MAXIMIZE)
        # 镜像仍 MINIMIZE：后续 sense=None 的 setObjective 不得翻转方向。
        self.m.setObjective(x)
        self.assertEqual(self.m.optimize(), Status.OPTIMAL)
        self.assertAlmostEqual(self.m.ObjVal, 0.0, delta=1e-9,
                               msg='镜像 _sense 须保持 MINIMIZE（ObjVal=1.0 即镜像已被翻）')


class NulNameRejectionTests(unittest.TestCase):
    """B3 回归：含内嵌 NUL 的字符串入口统一拒绝（ArgumentError 且为 ValueError）。

    不拦截则 FFI 的 CStr 在首个 NUL 截断而 Python 镜像保留全名，两侧静默
    分叉（实证：addVar(name='a\\x00b') 后 FFI 实存截断名 'a'，
    getVarByName('a') 返回 None 但 addVar(name='a') 报重名）。
    """

    def setUp(self):
        self.m = Model('nul-tests')

    def tearDown(self):
        self.m.dispose()

    def assert_nul_rejected(self, callable_, *args, **kwargs):
        with self.assertRaises(ArgumentError) as ctx:
            callable_(*args, **kwargs)
        self.assertIsInstance(ctx.exception, ValueError,
                              '含 NUL 的参数错误须可按 ValueError 捕获')
        self.assertIn('NUL', str(ctx.exception), '消息须点名 NUL 原因')

    def test_model_name_with_nul_rejected(self):
        self.assert_nul_rejected(Model, 'a\x00b')

    def test_addvar_name_with_nul_rejected_no_residual(self):
        self.assert_nul_rejected(self.m.addVar, name='a\x00b')
        self.assertEqual(self.m.NumVars, 0, '被拒 addVar 不留残余列')
        # 反转原 bug 场景：合法名 'a' 不被截断残影顶撞（不再误报重名）。
        self.m.addVar(name='a')

    def test_addvars_prefix_with_nul_rejected(self):
        self.assert_nul_rejected(self.m.addVars, 3, name='x\x00')
        self.assertEqual(self.m.NumVars, 0, '被拒 addVars 不留残余列')

    def test_addconstr_name_with_nul_rejected(self):
        v = self.m.addVar(name='v')
        self.assert_nul_rejected(self.m.addConstr, v <= 1, 'c\x00d')
        self.assertEqual(self.m.NumConstrs, 0, '被拒 addConstr 不留残余行')

    def test_setparam_str_value_with_nul_rejected(self):
        self.assert_nul_rejected(self.m.setParam, 'HeuristicMode', 'on\x00off')
        self.assertEqual(self.m.getParam('HeuristicMode'), 'on',
                         '被拒值不污染镜像（默认 on 保持）')


if __name__ == '__main__':
    unittest.main()
