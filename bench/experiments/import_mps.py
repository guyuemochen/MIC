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
    if lp.offset_ != 0 or any(v < 0 or not math.isfinite(v) for v in lp.col_lower_):
        raise ValueError('当前转换器拒绝目标常数、负下界或自由变量，不能静默改变模型')
    if lp.a_matrix_.format_ != highspy.MatrixFormat.kColwise:
        raise ValueError('需要按列存储的矩阵')
    rows = [[] for _ in range(lp.num_row_)]
    for j in range(lp.num_col_):
        for k in range(lp.a_matrix_.start_[j], lp.a_matrix_.start_[j + 1]):
            rows[lp.a_matrix_.index_[k]].append((j, lp.a_matrix_.value_[k]))
    lines = [f'name {source.stem}', 'sense min' if lp.sense_ == highspy.ObjSense.kMinimize else 'sense max', 'vars']
    for j in range(lp.num_col_):
        typ = lp.integrality_[j] if lp.integrality_ else highspy.HighsVarType.kContinuous
        if typ not in (highspy.HighsVarType.kContinuous, highspy.HighsVarType.kInteger):
            raise ValueError('不支持半连续/半整数变量')
        kind = 'continuous' if typ == highspy.HighsVarType.kContinuous else 'integer'
        if kind == 'integer' and lp.col_lower_[j] == 0 and lp.col_upper_[j] == 1:
            kind = 'binary'
        lines.append(f'x{j} {kind}')
    lines += ['obj', ' '.join(f'{v:.17g} x{j}' for j, v in enumerate(lp.col_cost_) if v != 0), 'cons']
    for i, terms in enumerate(rows):
        body = ' '.join(f'{v:.17g} x{j}' for j, v in terms)
        lo, hi = lp.row_lower_[i], lp.row_upper_[i]
        if lo == hi:
            lines.append(f'r{i}: {body} = {lo:.17g}')
        else:
            if math.isfinite(lo): lines.append(f'r{i}l: {body} >= {lo:.17g}')
            if math.isfinite(hi): lines.append(f'r{i}u: {body} <= {hi:.17g}')
    lines.append('bounds')
    for j in range(lp.num_col_):
        lines.append(f'x{j} >= {lp.col_lower_[j]:.17g}')
        if math.isfinite(lp.col_upper_[j]): lines.append(f'x{j} <= {lp.col_upper_[j]:.17g}')
    lines.append('end')
    target.write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('target', type=Path)
    args = parser.parse_args()
    convert(args.source, args.target)
