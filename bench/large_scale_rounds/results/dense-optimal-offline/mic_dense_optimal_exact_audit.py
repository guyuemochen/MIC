import argparse
import gzip
import hashlib
import json
import struct
from fractions import Fraction as F
from pathlib import Path

parser = argparse.ArgumentParser(description='从归档快照重算 p0201 稠密单纯形离线核验')
parser.add_argument('--input', type=Path,
                    default=Path(__file__).resolve().parents[1] / 'dense-optimal-profile/p0201-root.stdout.gz',
                    help='原始快照或 gzip 归档')
parser.add_argument('--output', type=Path, help='新报告路径，省略时将 JSON 打印到标准输出')
args = parser.parse_args()
def fl(b): return struct.unpack('>d', struct.pack('>Q', b))[0]
def fs(bits): return [fl(b) for b in bits]
def qs(bits): return [F(fl(b)) for b in bits]
def num(q): return {'float': float(q), 'exact': str(q)}
def extreme(items):
    i,v=max(items, key=lambda p:abs(p[1]))
    return {'index': i, **num(v)}

payload = args.input.read_bytes()
if args.input.suffix == '.gz':
    payload = gzip.decompress(payload)
data = json.loads(payload)
results=[]
for s in data['dense_optimal_profile']['samples']:
    n=s['num_structural']; nc=s['num_static']; basis=s['basis']
    rows=[[(j,F(fl(a))) for j,a in row] for row in s['row_terms_bits']]
    rhs=qs(s['row_rhs_bits']); c=qs(s['c_internal_bits']); y=qs(s['y_internal_all_rows_bits'])
    lb=qs(s['lb_bits']); ub=qs(s['ub_bits'])
    signs=qs(s['slack_sign_bits']); vals=qs(s['raw_static_and_art_values_bits'])
    tol=F(fl(s['tolerance_bits'][0])); opt=F(fl(s['tolerance_bits'][1]))
    z=F(fl(s['reported_internal_objective_bits']))
    off=F(fl(s['obj_offset_bits'])) * (-1 if s['sense_max'] else 1)
    def assess(x):
        delta=[sum((a*x[j] for j,a in row),F(0))-b for row,b in zip(rows,rhs)]
        violations=[max(d,0) if op==0 else max(-d,0) if op==1 else abs(d) for d,op in zip(delta,s['row_op'])]
        bound=[max(lo-v,v-hi,0) for v,lo,hi in zip(x,lb,ub)]
        obj=sum((a*v for a,v in zip(c,x)),off)
        return delta,violations,bound,obj
    r=c+[F(0)]*(nc-n)
    for i,row in enumerate(rows):
        for j,a in row:r[j]-=a*y[i]
        if s['row_slack'][i] is not None:r[s['row_slack'][i]]-=signs[s['row_slack'][i]-n]*y[i]
    nonbasics=[j for j in range(nc) if not s['in_basis'][j]]
    side=[(j,max(r[j] if s['at_upper'][j] else -r[j],0)) for j in nonbasics if j>=n or lb[j]!=ub[j]]
    bad_inf=[(j,r[j]) for j in range(n,nc) if r[j]<0]
    stats={'ordinal':s['optimal_ordinal'], 'returned_internal_objective':num(z), 'rows':len(rows),
        'negative_slack_rc':[{ 'column':j, **num(v)} for j,v in bad_inf],
        'basis_stationarity':extreme((j,r[j]) for j in basis if j<nc),
        'side_optimality':extreme(side),
        'first_three_state_equal':None}
    for name,key in [('raw','raw_x_bits'),('reported','reported_x_bits')]:
        x=qs(s[key]);delta,viol,bnd,obj=assess(x)
        stats[name]={'rows':extreme(enumerate(viol)), 'bounds':extreme(enumerate(bnd)),
            'row_violations_positive':sum(v>0 for v in viol), 'rows_over_feas':sum(v>tol for v in viol),
            'bounds_over_feas':sum(v>tol for v in bnd), 'objective_exact':num(obj),
            'returned_minus_exact_objective':num(z-obj)}
        if name=='raw':
            aug=[d+(signs[j-n]*vals[j] if j is not None else 0) for d,j in zip(delta,s['row_slack'])]
            k=max(range(len(aug)),key=lambda i:abs(aug[i]))
            stats['augmented']={**extreme(enumerate(aug)), 'row_original_residual':num(delta[k]),
                'row_relation':s['row_op'][k], 'stored_slack':num(vals[s['row_slack'][k]]),
                'term_absolute_sum':num(sum((abs(a*x[j]) for j,a in rows[k]),F(0))), 'rhs':num(rhs[k])}
    attempts=[]
    x0=qs(s['reported_x_bits'])
    for max_den in [1,2,3,4,5,6,8,10,16,32,64,100,1000,10000,1000000]:
        x=[v.limit_denominator(max_den) for v in x0]
        delta,viol,bnd,obj=assess(x)
        maxv=max(viol); maxb=max(bnd)
        attempt={'max_denominator':max_den, 'max_row_violation':float(maxv), 'max_bound_violation':float(maxb),
            'exact_feasible':maxv==0 and maxb==0, 'objective':num(obj), 'return_minus_candidate':num(z-obj)}
        attempts.append(attempt)
        if maxv==0 and maxb==0:
            attempt['witness_x']=[str(v) for v in x]
            break
    stats['rationalized_candidates']=attempts
    results.append(stats)
out={'input':str(args.input), 'input_sha256':hashlib.sha256(payload).hexdigest(),
    'scope':'只核验保存的规范化当前LP binary64输入；无矩阵求解或优化器', 'samples':results}
report = json.dumps(out, ensure_ascii=False, indent=2)
if args.output is None:
    print(report)
else:
    with args.output.open('x', encoding='utf-8') as stream:
        stream.write(report + '\n')
