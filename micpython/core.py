# -*- coding: utf-8 -*-
"""micpython 公开 API：gurobipy 风格建模 + milp 求解（Py3.8+，零第三方依赖）。

能力边界与语义纪律
==================

- **X/ObjVal 反映最近一次 optimize**（FFI 侧每次 solve 整体覆盖存档）；
  未求解或 Infeasible / Unbounded / LIMIT_REACHED-无 incumbent 时取
  X/ObjVal/Runtime/NodeCount/Bound/MIPGapPercent → ``NoSolutionError``；
  LIMIT_REACHED 有 incumbent 时可用（Status 已诚实标明 9=LIMIT_REACHED）；
- **solve 后结构性修改**：再 addVar 的新变量取值 → ``NoSolutionError``
  （解维数锚定 n_solved），旧变量仍取最近一次 solve 的值；
- **混用不同 Model 的变量** → ``ArgumentError``（表达式/约束/目标构造时
  即拒绝）；
- **句柄守卫**：一切方法（addVar/addConstr/setParam/optimize/属性读取）
  仅经 ``Model._call`` 触碰 FFI；``dispose()`` 后使用 → ``MicError``。
  句柄非线程安全：禁止多线程并发使用同一 Model。

参数镜像三铁律（setParam/getParam）
==================================

1. **恒调 FFI，无跳过分支**：对任何合法值都调用 FFI setter。TimeLimit 的
   ``INF``（Python 语义『不限』）翻译为 0.0 传入（FFI 0=不限；FFI 拒绝
   非有限值，故需该翻译）；MIPGap 的小数 ×100 转百分比传入；LogPeriodSecs
   关断直接传 0.0（FFI 0=关）。先 ``TimeLimit=5.0`` 再 ``TimeLimit=INF``
   ⇒ FFI 收到 0.0 ⇒ 真正解除限速，无『残留 5 秒照旧限速』的失真；
2. **镜像仅在 FFI 返回 0 后更新**：被 FFI 拒绝的值（如 Threads=0）不污染
   ``_params`` 镜像，getParam 不会返回已被拒绝的值；
3. **getParam 返回镜像原值，不做逆运算**：MIPGap 0.07 必须精确回读 0.07
   （不做 ÷100 反算，避免 7.000000000000001 的浮点噪声）。

与 gurobipy / milp CLI 的已知差异（docstring 即文档，README 由另行落盘）：

- ``setObjective`` 的**合并约定**：目标项 = 表达式项 ∪ {未被表达式覆盖的
  ``addVar(obj=...)`` 系数项}（gurobipy 为完全替换语义）；
- ``HeuristicMode`` 默认 ``'on'``（而非 milp CLI 的 auto）：对齐
  SearchControl::default 库语义，避免 auto+TimeLimit 意外自动降级为
  LIMIT_REACHED；需要自适应行为请显式 setParam('HeuristicMode', 'auto')；
- 两次同名 ``addVars`` 第二次抛 ``ModelError``（FFI/milp 拒绝重名变量，
  不自作聪明去重）；
- 名字/字符串参数不得含内嵌 NUL（``'\\x00'``）——模型名、变量名（含
  addVars 前缀）、约束名、str 参数值均抛 ``ArgumentError``（亦为
  ``ValueError`` 子类）：FFI 的 C 字符串会在首个 NUL 截断而 Python 镜像
  保留全名，两侧将静默分叉，入口统一拦截（B3）。
"""

import ctypes

from . import _ffi
from ._ffi import (
    ArgumentError,
    FfiLoadError,
    InternalError,
    IoError,
    MicError,
    ModelError,
    NoSolutionError,
    NumericalError,
    PanicError,
    ParseError,
    UnknownParameterError,
    check,
)

# ---------------------------------------------------------------------------
# 常量（gurobipy 同值）
# ---------------------------------------------------------------------------

BINARY = 'B'
INTEGER = 'I'
CONTINUOUS = 'C'

MAXIMIZE = -1
MINIMIZE = 1

INF = float('inf')

# FFI 枚举编码（milp-ffi/src/api.rs）
_VTYPE_CODES = {CONTINUOUS: 1, INTEGER: 2, BINARY: 3}
_SENSE_CODES = {MINIMIZE: 1, MAXIMIZE: 2}
_RELOP_LE, _RELOP_GE, _RELOP_EQ = 1, 2, 3

# .milp 语法保留字（与 milp/src/model.rs:62-75 的 RESERVED_WORDS 一致）；
# 仅用于**自动命名**时避开——显式撞保留字交 FFI 报错，保持单一事实源。
_RESERVED_WORDS = frozenset((
    'name', 'sense', 'vars', 'obj', 'cons', 'bounds', 'end',
    'continuous', 'integer', 'binary', 'max', 'min',
))

# get_stat 的 stat_id（milp-ffi/src/api.rs:881-883）
_STAT_NODES = 1
_STAT_LP_ITERATIONS = 2
_STAT_CUTS = 3
_STAT_TIME_MS = 4
_STAT_PRESOLVE_REDUCTIONS = 5
_STAT_GAP_TERMINATED = 6
_STAT_AUTO_DEGRADED = 7


class Status(object):
    """求解状态（gurobipy 数值风格；FFI 紧凑码 0..4 在本层映射）。"""

    NOT_SOLVED = 0
    OPTIMAL = 2
    INFEASIBLE = 3
    UNBOUNDED = 5
    LIMIT_REACHED = 9
    TIME_LIMIT = 9  # gurobipy 别名（与 LIMIT_REACHED 同值）


_STATUS_FROM_FFI = {
    0: Status.NOT_SOLVED,
    1: Status.OPTIMAL,
    2: Status.INFEASIBLE,
    3: Status.UNBOUNDED,
    4: Status.LIMIT_REACHED,
}


class GRB(object):
    """gurobipy 兼容别名命名空间（GRB.BINARY / GRB.MAXIMIZE / GRB.OPTIMAL …）。"""

    BINARY = BINARY
    INTEGER = INTEGER
    CONTINUOUS = CONTINUOUS
    MAXIMIZE = MAXIMIZE
    MINIMIZE = MINIMIZE
    INF = INF
    INFINITY = INF
    NOT_SOLVED = Status.NOT_SOLVED
    OPTIMAL = Status.OPTIMAL
    INFEASIBLE = Status.INFEASIBLE
    UNBOUNDED = Status.UNBOUNDED
    LIMIT_REACHED = Status.LIMIT_REACHED
    TIME_LIMIT = Status.TIME_LIMIT


# ---------------------------------------------------------------------------
# 内部小工具
# ---------------------------------------------------------------------------

def _merged_model(m1, m2):
    """合并两个 Model 引用；冲突（混用不同 Model）抛 ArgumentError。"""
    if m1 is None:
        return m2
    if m2 is None or m1 is m2:
        return m1
    raise ArgumentError('不能在同一个表达式/约束中混用来自不同 Model 的变量')


def _as_float(value, what):
    """数值类型强制（int/float→float）；失败抛 ArgumentError。"""
    if isinstance(value, (int, float)):
        return float(value)
    raise ArgumentError('%s 需要数值（int/float），收到 %r' % (what, value))


def _as_int(value, what):
    """整数类型强制（整数或整值 float）；失败抛 ArgumentError。"""
    if isinstance(value, bool):
        raise ArgumentError('%s 需要整数，收到 %r' % (what, value))
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    raise ArgumentError('%s 需要整数，收到 %r' % (what, value))


def _drop_zeros(terms):
    """就地删除恰好为 0.0 的系数项（对齐 milp LinearExpr『无零系数』不变量）。"""
    for idx in [i for i, c in terms.items() if c == 0.0]:
        del terms[idx]


def _pack_terms(terms):
    """{index: coeff} → (升序 (c_uint32*n), (c_double*n), n)；空返回 (None, None, 0)。

    count=0 时两个数组指针传 NULL（milp-ffi 契约：count=0 允许 NULL，空行/
    空目标合法）。
    """
    n = len(terms)
    if n == 0:
        return None, None, 0
    items = sorted(terms.items())
    idx_arr = (ctypes.c_uint32 * n)(*[int(i) for i, _ in items])
    coef_arr = (ctypes.c_double * n)(*[float(c) for _, c in items])
    return idx_arr, coef_arr, n


def _check_no_nul(value, what):
    """含内嵌 NUL（``'\\x00'``）的字符串入口统一拒绝（B3）。

    FFI 的 C 字符串（CStr）在首个 NUL 处截断，而 Python 侧镜像保留全名
    ——不拦截则两侧静默分叉（实证：``addVar(name='a\\x00b')`` 后 FFI 实存
    截断名 ``'a'``，``getVarByName('a')`` 却返回 None）。抛 ``ArgumentError``
    （亦为 ``ValueError`` 子类），消息点名入参与原因。
    """
    if '\x00' in value:
        raise ArgumentError(
            '%s 不得包含 NUL 字符（\\x00）：FFI 的 C 字符串会在首个 NUL 处'
            '截断，Python 镜像与求解器将各持一半名字（静默分叉）。收到 %r'
            % (what, value))


# ---------------------------------------------------------------------------
# 表达式与约束
# ---------------------------------------------------------------------------

class LinExpr(object):
    """线性表达式：内部 ``{var_index: coeff}`` 稀疏项 + ``constant`` + 所属
    Model 引用。

    - 由 Var / LinExpr / 常量经 ``+``、``-``、``*``（双侧数乘）、``-x``
      组合；``Var*Var`` 等二次表达式不在 milp 能力边界内（TypeError）；
    - ``==``/``<=``/``>=``（对常量或另一表达式，双表达式移项合并常数）
      产生 :class:`TempConstr`；
    - ``len()`` 为变量项数；``getValue()`` 需模型已求解且有解；
    - 混用不同 Model 的变量 → ``ArgumentError``。
    """

    def __init__(self, expr=None, constant=0.0, model=None):
        self.constant = float(constant)
        self._terms = {}
        self._model = model
        if expr is None:
            return
        if isinstance(expr, Var):
            self._model = _merged_model(self._model, expr._model)
            self._terms[expr._index] = 1.0
        elif isinstance(expr, LinExpr):
            self._model = _merged_model(self._model, expr._model)
            self._terms.update(expr._terms)
            self.constant += expr.constant
        elif isinstance(expr, (int, float)):
            self.constant += float(expr)
        else:
            raise ArgumentError('无法由 %r 构造 LinExpr' % (expr,))

    @staticmethod
    def _coerce(other):
        """Var/LinExpr/数值 → (terms, constant, model)；其余返回 None。"""
        if isinstance(other, Var):
            return {other._index: 1.0}, 0.0, other._model
        if isinstance(other, LinExpr):
            return other._terms, other.constant, other._model
        if isinstance(other, (int, float)):
            return {}, float(other), None
        return None

    def _combine(self, other, sign):
        """self + sign*other → 新 LinExpr；不可组合返回 None（给 NotImplemented）。"""
        coerced = self._coerce(other)
        if coerced is None:
            return None
        terms, constant, model = coerced
        new = LinExpr()
        new.constant = self.constant + sign * constant
        new._model = _merged_model(self._model, model)
        for idx, coeff in self._terms.items():
            new._terms[idx] = new._terms.get(idx, 0.0) + coeff
        for idx, coeff in terms.items():
            new._terms[idx] = new._terms.get(idx, 0.0) + sign * coeff
        _drop_zeros(new._terms)
        return new

    def __add__(self, other):
        result = self._combine(other, 1.0)
        return result if result is not None else NotImplemented

    __radd__ = __add__

    def __sub__(self, other):
        result = self._combine(other, -1.0)
        return result if result is not None else NotImplemented

    def __rsub__(self, other):
        coerced = self._coerce(other)
        if coerced is None:
            return NotImplemented
        terms, constant, model = coerced  # other - self
        new = LinExpr()
        new.constant = constant - self.constant
        new._model = _merged_model(self._model, model)
        for idx, coeff in terms.items():
            new._terms[idx] = new._terms.get(idx, 0.0) + coeff
        for idx, coeff in self._terms.items():
            new._terms[idx] = new._terms.get(idx, 0.0) - coeff
        _drop_zeros(new._terms)
        return new

    def __neg__(self):
        new = LinExpr()
        new.constant = -self.constant
        new._model = self._model
        for idx, coeff in self._terms.items():
            new._terms[idx] = -coeff
        return new

    def __mul__(self, other):
        if not isinstance(other, (int, float)):
            return NotImplemented  # Var*Var / LinExpr*LinExpr → TypeError
        factor = float(other)
        new = LinExpr()
        new.constant = self.constant * factor
        new._model = self._model
        for idx, coeff in self._terms.items():
            new._terms[idx] = coeff * factor
        _drop_zeros(new._terms)
        return new

    __rmul__ = __mul__

    def __len__(self):
        return len(self._terms)

    def _constr(self, other, op):
        """self op other → TempConstr（双表达式移项合并常数：termsA-termsB
        op cB-cA）。不可比较返回 None（给 NotImplemented）。"""
        coerced = self._coerce(other)
        if coerced is None:
            return NotImplemented
        terms, constant, model = coerced
        merged = dict(self._terms)
        for idx, coeff in terms.items():
            merged[idx] = merged.get(idx, 0.0) - coeff
        _drop_zeros(merged)
        rhs = constant - self.constant
        return TempConstr(merged, rhs, op, _merged_model(self._model, model))

    def __le__(self, other):
        return self._constr(other, _RELOP_LE)

    def __ge__(self, other):
        return self._constr(other, _RELOP_GE)

    def __eq__(self, other):
        return self._constr(other, _RELOP_EQ)

    def getConstant(self):
        """gurobipy 兼容：返回常数项。"""
        return self.constant

    def getValue(self):
        """在最近一次 optimize 的解上求值（未求解/无解 → NoSolutionError）。"""
        if self._model is None:
            return self.constant
        total = self.constant
        for idx, coeff in self._terms.items():
            total += coeff * self._model._vars[idx].X
        return total

    def __repr__(self):
        names = {}
        if self._model is not None:
            names = {v._index: v._name for v in self._model._vars}
        parts = []
        for idx in sorted(self._terms):
            parts.append('%g*%s' % (self._terms[idx], names.get(idx, 'v%d' % idx)))
        if self.constant != 0.0 or not parts:
            parts.append('%g' % self.constant)
        return 'LinExpr(%s)' % ' + '.join(parts).replace('+ -', '- ')


class TempConstr(object):
    """临时约束（==/<=/>= 比较的产物），交 :meth:`Model.addConstr` 落模。

    内部已移项归一：``_terms op _rhs``（双表达式比较时常数已合并到右端）。
    """

    def __init__(self, terms, rhs, op, model):
        self._terms = terms  # {var_index: coeff}（升序无关，已去零）
        self._rhs = float(rhs)
        self._op = op        # 1=LE, 2=GE, 3=EQ
        self._model = model

    def __repr__(self):
        symbol = {1: '<=', 2: '>=', 3: '=='}[self._op]
        items = sorted(self._terms.items())
        lhs = ' + '.join('%g*v%d' % (c, i) for i, c in items) or '0'
        return '<TempConstr %s %s %g>' % (lhs, symbol, self._rhs)


def quicksum(terms):
    """顺序累加 iterable 中的 Var/LinExpr/数值；空 iterable → ``LinExpr(0)``。

    混用不同 Model 的变量 → ``ArgumentError``；不可累加元素 → ``TypeError``
    （与 ``LinExpr + x`` 返回 NotImplemented 后解释器报错同款语义）。

    实现注（S1）：**就地累加**共享同一 ``_terms`` dict，每项 O(1) 摊还、
    整体 O(n)。旧实现 ``total = total + term`` 经 ``LinExpr.__add__`` 每步
    新建 LinExpr 并整拷增长中的 dict（O(n²)：5000 项 ≈1.27s，是
    python_bench knapsack-5000 build 2.35s 的绝对大头——逐变量 ctypes
    调用仅 ≈0.03s）。数学语义与旧实现逐项等价：浮点加法同序（系数终值
    bit 级一致），仅零系数清理从『逐步删』改为『末尾一次删』（对运行和，
    中途归零项被删后再遇同下标时 ``get(idx, 0.0)`` 重取同一累加值，终态
    dict 内容不变）。
    """
    total = LinExpr()
    accumulate = total._terms
    model = None
    for term in terms:
        coerced = LinExpr._coerce(term)
        if coerced is None:
            raise TypeError("unsupported operand type(s) for +: "
                            "'LinExpr' and '%s'" % type(term).__name__)
        t_terms, t_const, t_model = coerced
        model = _merged_model(model, t_model)
        for idx, coeff in t_terms.items():
            accumulate[idx] = accumulate.get(idx, 0.0) + coeff
        total.constant += t_const
    _drop_zeros(accumulate)
    total._model = model
    return total


# ---------------------------------------------------------------------------
# Var / Constr
# ---------------------------------------------------------------------------

class Var(object):
    """决策变量：持有所属 Model 的强引用（保住 FFI 句柄存活）与列下标。

    属性：``VarName``/``lb``/``ub``/``VType``/``Obj`` 为建模镜像（只读）；
    ``X`` 为最近一次 optimize 的解值——未求解或无解 → ``NoSolutionError``
    （消息 "Unable to retrieve attribute 'X'"）；模型已 dispose →
    ``MicError``（经 Model._call 同守卫）。支持全部 LinExpr 运算符。
    """

    def __init__(self, model, index, name, lb, ub, vtype, obj):
        self._model = model
        self._index = index
        self._name = name
        self._lb = lb
        self._ub = ub
        self._vtype = vtype
        self._obj = obj

    # ---- 属性 ----
    @property
    def VarName(self):
        return self._name

    @property
    def lb(self):
        return self._lb

    @property
    def ub(self):
        return self._ub

    @property
    def VType(self):
        return self._vtype

    @property
    def Obj(self):
        return self._obj

    @property
    def X(self):
        return self._model._var_value(self._index)

    # ---- LinExpr 运算符（比较运算产生 TempConstr，gurobipy 语义）----
    def __add__(self, other):
        return LinExpr(self) + other

    __radd__ = __add__

    def __sub__(self, other):
        return LinExpr(self) - other

    def __rsub__(self, other):
        return other - LinExpr(self)

    def __mul__(self, other):
        return LinExpr(self) * other

    __rmul__ = __mul__

    def __neg__(self):
        return -LinExpr(self)

    def __le__(self, other):
        return LinExpr(self) <= other

    def __ge__(self, other):
        return LinExpr(self) >= other

    def __eq__(self, other):
        return LinExpr(self) == other

    # __eq__ 返回 TempConstr（gurobipy 同款语义），显式保留基于身份的哈希。
    __hash__ = object.__hash__

    def __repr__(self):
        return '<Var %s (%s)>' % (self._name, self._vtype)


class Constr(object):
    """已落模约束（``.ConstrName``；下标为 FFI 返回的行号）。"""

    def __init__(self, model, index, name):
        self._model = model
        self._index = index
        self._name = name

    @property
    def ConstrName(self):
        return self._name

    def __repr__(self):
        return '<Constr %s>' % (self._name or '#%d' % self._index)


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

# 参数路由表（类型强制）与初始镜像（与 FFI Params::default 语义一致：
# TimeLimit 0=不限 ↔ INF；HeuristicMode 默认 'on' 而非 CLI 的 auto）。
_PARAM_TYPES = {
    'TimeLimit': float,
    'Threads': int,
    'MIPGap': float,
    'HeuristicMode': str,
    'LogPeriodSecs': float,
}

_DEFAULT_PARAMS = {
    'TimeLimit': INF,
    'Threads': 1,
    'MIPGap': 0.0,
    'HeuristicMode': 'on',
    'LogPeriodSecs': 0.0,
}


class Model(object):
    """MILP 模型（gurobipy 语义，milp 能力边界内）。

    生命周期：构造时装载 FFI 库并创建句柄（``_handle`` 初 None → create 后
    非 None）；``dispose()`` 显式释放——第一行短路 ``if self._handle is
    None: return``（二次 dispose 静默返回，不抛、不经 _call），否则经
    ``_call`` 调 destroy 后置 ``_handle=None``；``__del__`` 同模式（仅
    ``_handle`` 非 None 时调 destroy 再置 None，与 dispose 互斥，永不多调；
    即使误调 destroy(NULL) 也安全）。支持 with 语句。

    一切方法仅经私有 ``_call(fn, *args)`` 触碰 FFI：``_handle`` 为 None →
    抛 ``MicError('Model 已 dispose，禁止再使用')`` 后才调用 + check。
    """

    def __init__(self, name='Model', sense=MINIMIZE):
        self._handle = None  # 第一行：保证 __init__ 任何失败路径下 __del__ 安全
        if sense not in (MINIMIZE, MAXIMIZE):
            raise ArgumentError('sense 必须为 MINIMIZE(1) 或 MAXIMIZE(-1)，收到 %r'
                                % (sense,))
        if not isinstance(name, str):
            raise ArgumentError('name 需为 str，收到 %r' % (name,))
        _check_no_nul(name, '模型名 name')
        self._model_name = name
        self._sense = sense
        self._vars = []
        self._vars_by_name = {}
        self._obj_coeff = {}
        self._params = dict(_DEFAULT_PARAMS)
        self._status = Status.NOT_SOLVED
        self._auto_counter = 0
        self._constr_count = 0
        lib = _ffi._lib()
        handle = lib.milp_ffi_model_create(name.encode('utf-8'), _SENSE_CODES[sense])
        if not handle:  # NULL：sense 非法（Python 已拦）或 UTF-8/分配失败
            message = _ffi.last_error_message()
            raise MicError('milp_ffi_model_create 失败：%s'
                           % (message or '返回 NULL 且无错误消息'), 0)
        self._handle = handle

    # ---- 句柄守卫与 dispose ----

    def _call(self, fn, *args, **kwargs):
        """统一 FFI 入口：dispose 守卫 → 调用 → check(rc)。

        返回 void 的符号（仅 ``milp_ffi_model_destroy``）按成功处理。
        私有关键字 ``_tolerate`` 列出的错误码原样返回而不抛（仅 SolCount
        的无解探测使用，不入公共语义）。
        """
        if self._handle is None:
            raise MicError('Model 已 dispose，禁止再使用')
        rc = fn(self._handle, *args)
        if rc is None:
            rc = 0  # milp_ffi_model_destroy 返回 void
        if rc != 0 and rc not in kwargs.get('_tolerate', ()):
            check(rc)
        return rc

    def dispose(self):
        """释放 FFI 句柄；幂等（二次调用静默返回），此后任何使用抛 MicError。"""
        if self._handle is None:
            return
        lib = _ffi._lib()
        self._call(lib.milp_ffi_model_destroy)
        self._handle = None

    def __del__(self):
        # 与 dispose 同模式：仅 _handle 非 None 时 destroy 再置 None（互斥、
        # 永不多调）；不触发新的 DLL 装载（_reset_cache 后的遗留句柄放弃
        # 回收，属测试专用路径）；吞掉一切异常（解释器关闭阶段）。
        handle = self._handle
        if handle is None:
            return
        try:
            lib = _ffi._lib_cache
            if lib is not None:
                lib.milp_ffi_model_destroy(handle)
        except Exception:
            pass
        finally:
            self._handle = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.dispose()
        return False

    # ---- 建模 ----

    def addVar(self, lb=0.0, ub=None, obj=0.0, vtype=CONTINUOUS, name=''):
        """添加变量，返回 :class:`Var`。

        - ``lb`` 默认 0.0；``ub=None`` 按类型默认（C/I→+INF，B→1.0）；
        - ``vtype``：'C'/'I'/'B'（= CONTINUOUS/INTEGER/BINARY = GRB.*）；
        - ``name`` 为空则自动生成唯一名（自增 + 查重，避开 .milp 保留字表
          ——仅用于自动命名；显式撞保留字交 FFI 报错，保持单一事实源）；
        - ``obj``≠0 记入目标系数镜像（并入规则见 :meth:`setObjective`）；
        - lb 非有限 / lb>ub / Binary 界越界 / 重名 / 撞保留字由 FFI 校验
          （码 9/2 → ArgumentError/ModelError）；
        - name 含内嵌 NUL → ``ArgumentError``（B3：FFI C 字符串会在首个
          NUL 截断，Python 镜像与求解器将各持一半名字，入口统一拒绝）。
        """
        if vtype not in _VTYPE_CODES:
            raise ArgumentError("vtype 必须为 'C'/'I'/'B'，收到 %r" % (vtype,))
        lb_f = _as_float(lb, 'lb')
        obj_f = _as_float(obj, 'obj')
        if ub is None:
            ub_f = 1.0 if vtype == BINARY else INF
        else:
            ub_f = _as_float(ub, 'ub')
        if not isinstance(name, str):
            raise ArgumentError('name 需为 str，收到 %r' % (name,))
        _check_no_nul(name, '变量名 name')
        if not name:
            name = self._next_auto_name()
        lib = _ffi._lib()
        index_box = ctypes.c_uint32()
        self._call(lib.milp_ffi_add_var, name.encode('utf-8'),
                   _VTYPE_CODES[vtype], lb_f, ub_f, ctypes.byref(index_box))
        index = index_box.value
        var = Var(self, index, name, lb_f, ub_f, vtype, obj_f)
        self._vars.append(var)
        self._vars_by_name[name] = var
        if obj_f != 0.0:
            self._obj_coeff[index] = obj_f
        return var

    def _next_auto_name(self):
        """自动生成唯一变量名（自增 + 查重 + 避开保留字）。"""
        while True:
            candidate = 'x%d' % self._auto_counter
            self._auto_counter += 1
            if candidate not in self._vars_by_name and candidate not in _RESERVED_WORDS:
                return candidate

    def addVars(self, count, lb=0.0, ub=None, obj=0.0, vtype=CONTINUOUS, name='x'):
        """批量添加变量，返回 ``list[Var]``（S1 起经 FFI 批量符号
        ``milp_ffi_add_vars_batch`` 单次跨界建 count 个变量）。

        - 第 i 个变量名为 ``f"{name}{i}"``——**唯一性仅保证单次调用内**；
        - **两次同名 addVars 第二次抛 ModelError**：批内重名/撞现有名由
          FFI 预检阶段整批拒绝（validate-then-mutate：任一元素非法零
          突变），micpython 不自作聪明去重——批量命名冲突请更换 name
          前缀；
        - ``obj`` 可为标量或长度等于 count 的序列（lb/ub/vtype 仅标量）；
          标量 obj 为非数值（含 None）时与旧逐个 addVar 路径同款抛
          ``ArgumentError``；
        - ``name`` 前缀含内嵌 NUL → ``ArgumentError``（B3，同 addVar）；
        - 语义收紧（原逐个路径的空缺）：count=0 时 vtype/lb/ub 的值域
          现在也会被校验（旧实现循环零次、静默放过非法标量）；
        - 界默认、BINARY 的 ub=1 规则、obj 镜像（``setObjective`` 并入
          规则）与逐个 addVar 完全一致（批量契约测试见 test_model.py）。
        """
        if isinstance(count, bool) or not isinstance(count, int):
            raise ArgumentError('count 需为非负 int，收到 %r' % (count,))
        if count < 0:
            raise ArgumentError('count 需为非负 int，收到 %d' % count)
        if vtype not in _VTYPE_CODES:
            raise ArgumentError("vtype 必须为 'C'/'I'/'B'，收到 %r" % (vtype,))
        lb_f = _as_float(lb, 'lb')
        if ub is None:
            ub_f = 1.0 if vtype == BINARY else INF
        else:
            ub_f = _as_float(ub, 'ub')
        if isinstance(name, str):
            _check_no_nul(name, 'addVars 的 name 前缀')
        if obj is None or isinstance(obj, (int, float)):
            # 标量路径与旧逐个 addVar 同款：None 亦经 _as_float 抛
            # ArgumentError（obj 需要数值）。
            obj_values = None
            obj_scalar = _as_float(obj, 'obj')
        elif isinstance(obj, (list, tuple)):
            if len(obj) != count:
                raise ArgumentError('obj 序列长度 %d 与 count %d 不符'
                                    % (len(obj), count))
            # 预检整个序列（旧路径在首个非法元素处的 addVar 抛，本侧在
            # 任何 FFI 调用前全量拒绝——零突变时点更早，消息同款）。
            obj_values = [_as_float(v, 'obj') for v in obj]
        else:
            raise ArgumentError('obj 需为标量或序列，收到 %r' % (obj,))
        if count == 0:
            return []
        names = ['%s%d' % (name, i) for i in range(count)]
        for nm in names:
            # 防御：非 str 的 name 经 %s 格式化仍可能引入 NUL（旧路径由
            # 逐个 addVar 的 _check_no_nul 兜住，批量路径在此等价拦截）。
            _check_no_nul(nm, '变量名 name')
        vtype_code = _VTYPE_CODES[vtype]
        # ctypes 数组类型严格声明（禁裸指针）：字节串由 c_char_p 数组持有
        # 引用，存活至调用返回。
        names_arr = (ctypes.c_char_p * count)(
            *[nm.encode('utf-8') for nm in names])
        vtypes_arr = (ctypes.c_int32 * count)(*([vtype_code] * count))
        lbs_arr = (ctypes.c_double * count)(*([lb_f] * count))
        ubs_arr = (ctypes.c_double * count)(*([ub_f] * count))
        lib = _ffi._lib()
        first_box = ctypes.c_uint32()
        self._call(lib.milp_ffi_add_vars_batch, count, names_arr, vtypes_arr,
                   lbs_arr, ubs_arr, ctypes.byref(first_box))
        first = first_box.value
        result = []
        for i, nm in enumerate(names):
            obj_i = obj_scalar if obj_values is None else obj_values[i]
            var = Var(self, first + i, nm, lb_f, ub_f, vtype, obj_i)
            result.append(var)
            self._vars.append(var)
            self._vars_by_name[nm] = var
            if obj_i != 0.0:
                self._obj_coeff[first + i] = obj_i
        return result

    def addConstr(self, constr, name=''):
        """落模一条 :class:`TempConstr`（由 <=/>=/== 比较产生），返回
        :class:`Constr`。

        约束引用其他 Model 的变量 → ``ArgumentError``；rhs/系数非法由 FFI
        校验。``name`` 为空则不命名（FFI 传 NULL）。
        """
        if not isinstance(constr, TempConstr):
            raise ArgumentError('addConstr 需要 TempConstr（由 <=/>=/== 比较产生），'
                                '收到 %s' % type(constr).__name__)
        if constr._model is not None and constr._model is not self:
            raise ArgumentError('约束引用了其他 Model 的变量')
        if not isinstance(name, str):
            raise ArgumentError('name 需为 str，收到 %r' % (name,))
        _check_no_nul(name, '约束名 name')
        indices, coeffs, n = _pack_terms(constr._terms)
        name_b = name.encode('utf-8') if name else None
        lib = _ffi._lib()
        row_box = ctypes.c_uint32()
        self._call(lib.milp_ffi_add_constraint, name_b, indices, coeffs, n,
                   constr._op, float(constr._rhs), ctypes.byref(row_box))
        self._constr_count += 1
        return Constr(self, row_box.value, name)

    def setObjective(self, expr, sense=None):
        """设置目标函数。

        - ``expr`` 可为 Var / LinExpr / 数值（数值即纯常数目标）；
        - ``sense`` 为 MAXIMIZE/MINIMIZE 时更新方向（透传 FFI，成功后才写
          problem.sense 与 Python 镜像——FFI 拒绝目标时方向零突变），
          None 保持现方向；
        - **合并约定**：最终目标项 = 表达式项 ∪ {未被表达式覆盖的
          ``addVar(obj=...)`` 系数项}——addVar 记录的 obj 系数在每次
          setObjective 时自动并入未被显式覆盖的列（与 gurobipy 的『完全
          替换』语义不同，为本绑定的明确约定）；``expr`` 的 constant 作为
          目标偏移 obj_offset 透传；
        - 表达式引用其他 Model 的变量 → ``ArgumentError``。
        """
        if isinstance(expr, Var):
            lexpr = LinExpr(expr)
        elif isinstance(expr, LinExpr):
            lexpr = expr
        elif isinstance(expr, (int, float)):
            lexpr = LinExpr(constant=float(expr))
        else:
            raise ArgumentError('目标需为 Var/LinExpr/数值，收到 %s'
                                % type(expr).__name__)
        if lexpr._model is not None and lexpr._model is not self:
            raise ArgumentError('目标表达式引用了其他 Model 的变量')
        if sense is not None and sense not in (MINIMIZE, MAXIMIZE):
            raise ArgumentError('sense 必须为 MINIMIZE(1) 或 MAXIMIZE(-1)，'
                                '收到 %r' % (sense,))
        # 目标方向本次应写入 FFI 的值：显式 sense 优先，None 保持现方向。
        target_sense = self._sense if sense is None else sense
        terms = dict(lexpr._terms)
        for idx, coeff in self._obj_coeff.items():
            if idx not in terms:
                terms[idx] = coeff
        indices, coeffs, n = _pack_terms(terms)
        lib = _ffi._lib()
        self._call(lib.milp_ffi_set_objective, indices, coeffs, n,
                   _SENSE_CODES[target_sense], float(lexpr.constant))
        # 镜像仅在 FFI 返回 0 后更新（参数三铁律同款纪律）：FFI 侧
        # set_objective 失败（如系数非有限→9）时 problem.sense 零突变
        # （B2），Python 镜像必须同步零突变，否则两侧方向分叉。
        if sense is not None:
            self._sense = sense

    # ---- 参数（三铁律见模块 docstring）----

    def setParam(self, name, value):
        """设置求解参数：路由表类型强制 → 翻译 → **恒调 FFI setter（无任何
        『跳过』分支）** → rc==0 才写镜像 ``_params[name]=原始 Python 值``。

        翻译规则：TimeLimit 的 INF→0.0（FFI 0=不限）；MIPGap 小数×100→
        百分比；LogPeriodSecs 关断传 0.0（FFI 0=关）。值域/类型非法由 FFI
        拒绝（码 9 → ArgumentError，被拒值不污染镜像）。未知名 →
        ``UnknownParameterError``。
        """
        if name not in _PARAM_TYPES:
            raise UnknownParameterError(
                '未知参数：%r（可用参数：%s）'
                % (name, ', '.join(sorted(_PARAM_TYPES))), 10)
        lib = _ffi._lib()
        name_b = name.encode('utf-8')
        if name == 'TimeLimit':
            v = _as_float(value, 'TimeLimit')
            ffi_value = 0.0 if v == INF else v  # INF=Python『不限』→ FFI 0=不限
            self._call(lib.milp_ffi_set_param_real, name_b, ffi_value)
        elif name == 'MIPGap':
            v = _as_float(value, 'MIPGap')
            self._call(lib.milp_ffi_set_param_real, name_b, v * 100.0)
        elif name == 'LogPeriodSecs':
            v = _as_float(value, 'LogPeriodSecs')
            self._call(lib.milp_ffi_set_param_real, name_b, v)  # 关断即传 0.0
        elif name == 'Threads':
            v = _as_int(value, 'Threads')
            self._call(lib.milp_ffi_set_param_int, name_b, v)
        else:  # HeuristicMode（'auto'/'on'/'off'，非法值交 FFI 拒绝→9）
            if not isinstance(value, str):
                raise ArgumentError("HeuristicMode 需为 str（'auto'/'on'/'off'），"
                                    '收到 %r' % (value,))
            _check_no_nul(value, "参数 HeuristicMode 的值 value")
            self._call(lib.milp_ffi_set_param_str, name_b, value.encode('utf-8'))
        # 仅 FFI 返回 0（_call 未抛）后到达：镜像写入原始 Python 值。
        self._params[name] = value

    def getParam(self, name):
        """读取参数镜像：返回 setParam 时的**原始 Python 值**（不做逆运算——
        MIPGap 0.07 精确回读 0.07）；未 setParam 过则返回初始镜像。未知名 →
        ``UnknownParameterError``。"""
        if name not in _PARAM_TYPES:
            raise UnknownParameterError(
                '未知参数：%r（可用参数：%s）'
                % (name, ', '.join(sorted(_PARAM_TYPES))), 10)
        return self._params[name]

    # ---- 求解 ----

    def optimize(self):
        """求解；返回 gurobipy 数值风格 Status（FFI 紧凑码映射：
        OPTIMAL=2/INFEASIBLE=3/UNBOUNDED=5/LIMIT_REACHED=9）。

        求解错误（FFI 码 2-7）抛对应异常（旧结果不覆盖：X/ObjVal 仍反映上
        一次成功 optimize）。
        """
        lib = _ffi._lib()
        self._call(lib.milp_ffi_solve)
        status_box = ctypes.c_int32()
        self._call(lib.milp_ffi_get_status, ctypes.byref(status_box))
        code = status_box.value
        try:
            self._status = _STATUS_FROM_FFI[code]
        except KeyError:
            raise InternalError('FFI 返回未知求解状态码 %d（契约值域 0..=4）'
                                % code, 6)
        return self._status

    # ---- 属性读取（全部经 _call，共享 dispose 守卫）----

    @property
    def Status(self):
        """最近一次 optimize 的终态；未求解 = NOT_SOLVED(0)。"""
        return self._status

    @property
    def ObjVal(self):
        """最近一次求解的目标值；无解 → NoSolutionError。"""
        return self._read_double('milp_ffi_get_objective', 'ObjVal')

    @property
    def SolCount(self):
        """当前解个数（0/1）：经 get_objective 探测（容忍码 11=无解）。"""
        lib = _ffi._lib()
        rc = self._call(lib.milp_ffi_get_objective,
                        ctypes.byref(ctypes.c_double()), _tolerate=(11,))
        return 1 if rc == 0 else 0

    @property
    def Runtime(self):
        """求解耗时秒 = stats.time_ms / 1000.0；未求解/无解 → NoSolutionError。"""
        return self._read_stat(_STAT_TIME_MS, 'Runtime') / 1000.0

    @property
    def NodeCount(self):
        """分支定界节点数；未求解/无解 → NoSolutionError。"""
        return self._read_stat(_STAT_NODES, 'NodeCount')

    @property
    def GapTerminated(self):
        """是否因达到 MIPGap 阈值提前终止（stat 6）。"""
        return bool(self._read_stat(_STAT_GAP_TERMINATED, 'GapTerminated'))

    @property
    def AutoDegraded(self):
        """是否因策略级自动降级提前终止（stat 7；仅 HeuristicMode='auto' 可触发）。"""
        return bool(self._read_stat(_STAT_AUTO_DEGRADED, 'AutoDegraded'))

    @property
    def Bound(self):
        """对偶界证据；无证据 → NoSolutionError。"""
        return self._read_double('milp_ffi_get_bound', 'Bound')

    @property
    def MIPGapPercent(self):
        """gap 触发时刻的百分比；未因 gap 终止 → NoSolutionError。"""
        return self._read_double('milp_ffi_get_gap_percent', 'MIPGapPercent')

    @property
    def NumVars(self):
        """变量数（当前模型，只读）。"""
        lib = _ffi._lib()
        out = ctypes.c_uint32()
        self._call(lib.milp_ffi_model_num_vars, ctypes.byref(out))
        return out.value

    @property
    def NumConstrs(self):
        """约束数（当前模型，只读）。"""
        lib = _ffi._lib()
        out = ctypes.c_uint32()
        self._call(lib.milp_ffi_model_num_constraints, ctypes.byref(out))
        return out.value

    @property
    def ModelName(self):
        """模型名（构造时指定；FFI 无改名符号，只读）。"""
        return self._model_name

    def getVarByName(self, name):
        """按名取 Var；不存在返回 None（Python 侧字典，不触碰 FFI）。"""
        return self._vars_by_name.get(name)

    # ---- 内部读取辅助（NoSolutionError 统一改写为 gurobipy 风格消息）----

    def _read_double(self, fn_name, attr):
        lib = _ffi._lib()
        out = ctypes.c_double()
        try:
            self._call(getattr(lib, fn_name), ctypes.byref(out))
        except NoSolutionError:
            raise NoSolutionError("Unable to retrieve attribute '%s'" % attr, 11)
        return out.value

    def _read_stat(self, stat_id, attr):
        lib = _ffi._lib()
        out = ctypes.c_int64()
        try:
            self._call(lib.milp_ffi_get_stat, stat_id, ctypes.byref(out))
        except NoSolutionError:
            raise NoSolutionError("Unable to retrieve attribute '%s'" % attr, 11)
        return out.value

    def _var_value(self, index):
        """Var.X 的实现：经 _call 守卫取单变量值。"""
        lib = _ffi._lib()
        out = ctypes.c_double()
        try:
            self._call(lib.milp_ffi_get_var_value, index, ctypes.byref(out))
        except NoSolutionError:
            raise NoSolutionError("Unable to retrieve attribute 'X'", 11)
        return out.value

    def __repr__(self):
        return '<Model %s (%s): %d vars, %d constraints>' % (
            self._model_name,
            'minimize' if self._sense == MINIMIZE else 'maximize',
            len(self._vars), self._constr_count)
