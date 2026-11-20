#!/usr/bin/env python3
"""用基准虚拟环境的 highspy 将可直接表达的 MPS 转为 MIC 格式；不改求解器依赖。"""
import argparse
import math
from pathlib import Path
import highspy


def convert(source, target):
    h = highspy.Highs()
    h.setOptionValue('output_flag', False)
    if h.readModel(str(source)) != highspy.HighsStatus.kOk:
        raise ValueError('MPS 读取失败')
    if h.getModel().hessian_.dim_ != 0:
        raise ValueError("当前转换器不支持二次目标")
    lp = h.getLp()
    lower, upper, costs, kinds = lp.col_lower_, lp.col_upper_, lp.col_cost_, lp.integrality_
    row_lower, row_upper = lp.row_lower_, lp.row_upper_
    if lp.offset_ != 0 or any(v < 0 or not math.isfinite(v) for v in lower):
        raise ValueError('当前转换器拒绝目标常数、负下界或自由变量，不能静默改变模型')
    matrix = lp.a_matrix_
    if matrix.format_ != highspy.MatrixFormat.kColwise:
        raise ValueError('需要按列存储的矩阵')
    # 绑定属性返回列表副本；缓存后逐项遍历只扫描一次稀疏矩阵。
    starts, indices, values = matrix.start_, matrix.index_, matrix.value_
    rows = [[] for _ in range(lp.num_row_)]
    for j in range(lp.num_col_):
        for k in range(starts[j], starts[j + 1]):
            rows[indices[k]].append((j, values[k]))
    lines = [f'name {source.stem}', 'sense min' if lp.sense_ == highspy.ObjSense.kMinimize else 'sense max', 'vars']
    for j in range(lp.num_col_):
        typ = kinds[j] if kinds else highspy.HighsVarType.kContinuous
        if typ not in (highspy.HighsVarType.kContinuous, highspy.HighsVarType.kInteger):
            raise ValueError('不支持半连续/半整数变量')
        kind = 'continuous' if typ == highspy.HighsVarType.kContinuous else 'integer'
        if kind == 'integer' and lower[j] == 0 and upper[j] == 1:
            kind = 'binary'
        lines.append(f'x{j} {kind}')
    lines += ['obj', ' '.join(f'{v:.17g} x{j}' for j, v in enumerate(costs) if v != 0), 'cons']
    for i, terms in enumerate(rows):
        body = ' '.join(f'{v:.17g} x{j}' for j, v in terms)
        lo, hi = row_lower[i], row_upper[i]
        if lo == hi:
            lines.append(f'r{i}: {body} = {lo:.17g}')
        else:
            if math.isfinite(lo): lines.append(f'r{i}l: {body} >= {lo:.17g}')
            if math.isfinite(hi): lines.append(f'r{i}u: {body} <= {hi:.17g}')
    lines.append('bounds')
    for j in range(lp.num_col_):
        lines.append(f'x{j} >= {lower[j]:.17g}')
        if math.isfinite(upper[j]): lines.append(f'x{j} <= {upper[j]:.17g}')
    lines.append('end')
    target.write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('target', type=Path)
    args = parser.parse_args()
    convert(args.source, args.target)
