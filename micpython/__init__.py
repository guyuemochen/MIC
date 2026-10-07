# -*- coding: utf-8 -*-
"""micpython —— milp 求解器的 gurobipy 风格 Python 绑定。

零第三方依赖（ctypes + 标准库），Python 3.8+。目标求解器为本地 FFI 库
``milp-ffi``（见 ``micpython/_ffi.py`` 的装载规则与 ABI 契约）。

- ``import micpython`` **无副作用**（惰性装载 DLL）：首次构造
  :class:`Model` 或显式调用 :func:`load` 时才装载 milp-ffi；
- DLL 未构建/不可用时上述操作抛 :class:`FfiLoadError`（消息含全部候选
  路径与构建提示 ``cargo build --manifest-path milp-ffi/Cargo.toml
  --release``）；
- 常量与 Status 采用 gurobipy 数值风格（MAXIMIZE=-1/MINIMIZE=1；
  OPTIMAL=2/INFEASIBLE=3/UNBOUNDED=5/LIMIT_REACHED=9）；
- 参数镜像三铁律、setObjective 合并约定、dispose 守卫等语义详见
  ``micpython/core.py`` 模块 docstring。

用法示例::

    import micpython as gp

    m = gp.Model('demo')
    x = m.addVar(vtype=gp.BINARY, ub=1.0, name='x')
    m.setObjective(1.0 * x, gp.MAXIMIZE)
    m.addConstr(x <= 1)
    m.optimize()        # -> 2 (OPTIMAL)
    print(m.ObjVal)     # -> 1.0
    print(x.X)          # -> 1.0
"""

from .core import (
    BINARY,
    CONTINUOUS,
    GRB,
    INF,
    INTEGER,
    MAXIMIZE,
    MINIMIZE,
    Status,
    Constr,
    LinExpr,
    Model,
    TempConstr,
    Var,
    quicksum,
    MicError,
    FfiLoadError,
    ArgumentError,
    ModelError,
    NumericalError,
    ParseError,
    IoError,
    InternalError,
    PanicError,
    UnknownParameterError,
    NoSolutionError,
)


def load():
    """显式预热：立即装载 milp-ffi 动态库并校验 ABI（失败抛 FfiLoadError）。"""
    from . import _ffi
    return _ffi._lib()


__version__ = '1.0.0'

__all__ = [
    'Model', 'Var', 'LinExpr', 'TempConstr', 'Constr', 'quicksum',
    'BINARY', 'INTEGER', 'CONTINUOUS', 'MAXIMIZE', 'MINIMIZE', 'INF',
    'Status', 'GRB',
    'MicError', 'FfiLoadError', 'ArgumentError', 'ModelError',
    'NumericalError', 'ParseError', 'IoError', 'InternalError',
    'PanicError', 'UnknownParameterError', 'NoSolutionError',
    'load',
]
