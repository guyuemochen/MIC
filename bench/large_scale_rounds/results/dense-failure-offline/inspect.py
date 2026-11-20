"""离线核对失败时的基方程；不修改或调用求解器。"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
import os
from pathlib import Path
import struct

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
import numpy as np

def f(bits):
    return struct.unpack('!d', struct.pack('!Q', bits))[0]

def vector(bits):
    return [f(b) for b in bits]

def largest(values):
    index = max(range(len(values)), key=lambda i: abs(values[i]))
    return {'index': index, 'value': float(values[index]), 'absolute': float(abs(values[index]))}

def inspect(s):
    assert s['snapshot_complete'] and s['current_B_reconstructable']
    n, ncs = s['num_structural'], s['num_static']
    basis = s['basis']; m = len(basis)
    assert len(set(basis)) == m and all(0 <= j < ncs for j in basis)
    lb, ub = vector(s['lb_bits']), vector(s['ub_bits'])
    rhs, b = vector(s['rhs_bits']), vector(s['row_rhs_bits'])
    pre = vector(s['pre_refresh_rhs_bits']) if s['pre_refresh_rhs_bits'] is not None else None
    inv = [vector(row) for row in s['binv_bits']]
    A = [[0.0] * ncs for _ in range(m)]
    for i, row in enumerate(s['row_terms_bits']):
        for j, bits in row:
            A[i][j] = f(bits)
        slack = s['row_slack'][i]
        if slack >= 0:
            A[i][slack] = f(s['slack_sign_bits'][slack - n])
    B = [[row[j] for j in basis] for row in A]
    assert len(rhs) == m and len(b) == m and len(inv) == m
    assert all(math.isfinite(v) for rows in (A, B, inv) for row in rows for v in row)
    assert all(math.isfinite(v) for v in rhs+b)
    assert all(s['in_basis'][j] == (j in basis) for j in range(len(s['in_basis'])))
    nonbasic = [0.0] * ncs
    for j in range(ncs):
        if j in basis:
            continue
        if j < n:
            nonbasic[j] = ub[j] if s['at_upper'][j] else lb[j]
        else:
            assert not s['at_upper'][j], '非基松弛不能贴无穷上界'
    assert all(math.isfinite(v) for v in nonbasic)
    # 完全重放生产刷新算序，并另算保存浮点输入的精确二进制有理右端。
    d, de = [], []
    for i, terms in enumerate(s['row_terms_bits']):
        acc = 0.0
        for j, bits in terms:
            if nonbasic[j] != 0.0:
                acc += f(bits) * nonbasic[j]
        d.append(b[i] - acc)
        de.append(F(b[i]) - sum((F(f(bits)) * F(nonbasic[j]) for j,bits in terms), F(0)))
    replay = []
    for row in inv:
        acc = 0.0
        for k, v in enumerate(row):
            if v != 0.0:
                acc += v * d[k]
        replay.append(acc)
    precise_stored_d = [sum((F(v)*F(d[k]) for k,v in enumerate(row)), F(0)) for row in inv]
    precise_exact_d = [sum((F(v)*de[k] for k,v in enumerate(row)), F(0)) for row in inv]
    fresh = np.linalg.solve(np.array(B), np.array(d)).tolist()
    fresh_exact_d = np.linalg.solve(np.array(B), np.array([float(v) for v in de])).tolist()
    basis_bounds = [(lb[j], ub[j]) if j < n else (0.0, float('inf')) for j in basis]
    def report(values):
        vf = [float(v) for v in values]
        x = list(nonbasic[:n])
        for i,j in enumerate(basis):
            if j < n: x[j] = vf[i]
        bounds = [max(l-v, v-u, 0.0) for v,(l,u) in zip(vf,basis_bounds)]
        residual = [sum((F(B[i][k])*F(values[k]) for k in range(m)), F(0))-de[i] for i in range(m)]
        row_errors = []
        for i in range(m):
            residual_original = sum((F(A[i][j])*F(x[j]) for j in range(n)),F(0))-F(b[i])
            op=s['row_op'][i]
            row_errors.append(max(residual_original,F(0)) if op==0 else max(-residual_original,F(0)) if op==1 else abs(residual_original))
        return {'basis_bound_violation':largest(bounds),'basis_residual_exact_inputs':largest(residual),
                'original_normalized_row_violation':largest(row_errors),'violated_values':{str(i):vf[i] for i in s['violated_rows']}}
    inverse_error = np.array(inv) @ np.array(B) - np.eye(m)
    bad_i,bad_j = np.unravel_index(np.abs(inverse_error).argmax(),inverse_error.shape)
    exact_inverse_entry = sum((F(inv[bad_i][k])*F(B[k][bad_j]) for k in range(m)),F(0))-F(int(bad_i==bad_j))
    return {'event':s['event'],'call_id':s['call_id'],'attempt':s['attempt'],'site':s['site'],
            'iter_call':s['iter_call'],'basis_rows':m,'structural_columns':n,'feasibility':f(s['feasibility_bits']),
            'refresh_replay_bit_equal':all(struct.pack('!d',a)==struct.pack('!d',v) for a,v in zip(replay,rhs)),
            'd_rounding_error':largest([F(v)-e for v,e in zip(d,de)]),
            'stored_inverse_identity_error':{'row':int(bad_i),'column':int(bad_j),'float_scan_max':float(np.max(np.abs(inverse_error))),'exact_entry':float(exact_inverse_entry)},
            'basis_condition_estimate_inf':float(np.linalg.cond(np.array(B),np.inf)),
            'saved':report(rhs),'before_refresh':report(pre) if pre else None,
            'same_inverse_exact_arithmetic_stored_d':report(precise_stored_d),
            'same_inverse_exact_arithmetic_exact_d':report(precise_exact_d),
            'fresh_basis_stored_d':report(fresh),'fresh_basis_rounded_exact_d':report(fresh_exact_d),
            'fresh_difference_saved':largest([a-v for a,v in zip(fresh,rhs)]),
            'basis_columns_at_violations':{str(i):basis[i] for i in s['violated_rows']}}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path); parser.add_argument('output',type=Path)
    args=parser.parse_args()
    assert not args.output.exists()
    data=args.input.read_bytes()
    rows=[json.loads(line) for line in data.decode().splitlines()]
    result={'input':str(args.input),'input_sha256':hashlib.sha256(data).hexdigest(),'numpy':np.__version__,
            'scope':'保存的浮点输入的精确残差与独立双精度基解；不是精确最优或不可行证书',
            'cases':[inspect(s) for s in rows]}
    args.output.write_text(json.dumps(result,ensure_ascii=False,allow_nan=False,indent=2)+'\n')
    for c in result['cases']:
        print(c['event'],c['call_id'],c['attempt'],c['site'],c['iter_call'],
              '旧',c['saved']['basis_bound_violation'],'重解',c['fresh_basis_stored_d']['basis_bound_violation'])

if __name__=='__main__':main()
