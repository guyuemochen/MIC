#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""knapsack.py —— gurobipy 风格的 micpython 示例（0-1 背包）。

模型等价 ``milp/instances/knapsack-10.milp``（acceptance.rs CASES 同向量）::

    max  7x1 + 5x2 + 9x3 + 6x4 + 3x5
    s.t. 4x1 + 3x2 + 5x3 + 6x4 + 2x5 <= 10,  x ∈ {0,1}

最优解 x2 = x3 = x5 = 1（0 基下标 1/2/4），目标 5 + 9 + 3 = 17。

运行（工作区根）::

    python micpython/examples/knapsack.py

成功时逐行打印各变量取值，最后一行恰为 ``objective = 17.000000``（数值
格式与 milp CLI 的 objective 行一致），退出码 0；非 OPTIMAL 打印原因并
退出码 1。
"""

import sys
from pathlib import Path

# 脚本位于 micpython/examples/，仓库根为 parents[2]——包导入路径由脚本自理
# （micpython 为路径导入包，无需安装）。
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from micpython import BINARY, MAXIMIZE, Model, Status, quicksum


def _format_number_fallback(value):
    """纯 Python 兜底格式化：与 milp::search_ctrl::format_number 同规则。

    ``0 < |v| < 1e-6`` 或 ``|v| >= 1e15`` 用科学记数（尾数去零、指数无前导
    零，对齐 Rust ``{:e}`` 输出，如 ``1e-7`` / ``-2.5e18``），否则定点 6 位
    小数（``{:.6f}``）。非有限值按 Rust Display 输出（NaN/inf/-inf）。
    """
    v = float(value)
    if v != v:  # NaN
        return 'NaN'
    if v == float('inf'):
        return 'inf'
    if v == float('-inf'):
        return '-inf'
    magnitude = abs(v)
    if 0.0 < magnitude < 1e-6 or magnitude >= 1e15:
        text = '{:e}'.format(v)
        mantissa, _, exponent = text.partition('e')
        if '.' in mantissa:
            mantissa = mantissa.rstrip('0').rstrip('.')
        negative = exponent.startswith('-')
        digits = exponent.lstrip('+-').lstrip('0') or '0'
        return '%se%s%s' % (mantissa, '-' if negative else '', digits)
    return '{:.6f}'.format(v)


def format_number(value):
    """优先经 ABI ``milp_ffi_format_number``（与 milp CLI 逐字节一致，如
    ``17.000000``）；FFI 不可用时退回纯 Python 等价实现。"""
    try:
        from micpython import _ffi
        return _ffi.format_number(float(value))
    except Exception:
        return _format_number_fallback(value)


def main():
    m = Model("knapsack")

    x = m.addVars(5, vtype=BINARY, name="x")
    values = [7, 5, 9, 6, 3]
    weights = [4, 3, 5, 6, 2]
    capacity = 10

    m.setObjective(quicksum(values[i] * x[i] for i in range(5)), MAXIMIZE)
    m.addConstr(quicksum(weights[i] * x[i] for i in range(5)) <= capacity)

    status = m.optimize()
    if status != Status.OPTIMAL:
        print("求解未达最优：status = %d（期望 OPTIMAL = %d）"
              % (status, Status.OPTIMAL))
        return 1

    for v in x:
        print(f"{v.VarName} = {v.X:.0f}")
    print(f"objective = {format_number(m.ObjVal)}")

    m.dispose()
    return 0


if __name__ == '__main__':
    sys.exit(main())
