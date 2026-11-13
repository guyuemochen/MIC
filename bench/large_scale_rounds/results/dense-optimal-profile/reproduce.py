import json,pathlib,subprocess,hashlib,sys
root=pathlib.Path('/tmp/mic-dense-failure-profile');out=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else root/'bench/large_scale_rounds/results/dense-optimal-profile'
out.mkdir(parents=True,exist_ok=True)
assert not (out/'p0201-root.stdout').exists(), '复现必须指定新目录，不能覆盖原始证据'
sys.path.insert(0,str(root/'bench/large_scale_rounds'))
from check_results import verify_solution,inspect_results
from run_propagation import canonical_hash
main=pathlib.Path('/Users/guyuemochen/Desktop/projects/MIC')
binaries={side:pathlib.Path('/tmp/mic-dense-optimal-binaries')/('round_probe_'+side) for side in ('off','on')}
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def write(name,value): (out/name).write_text(json.dumps(value,ensure_ascii=False,allow_nan=False,indent=2)+'\n')
def invoke(name,cmd,timeout=90):
    r=subprocess.run(cmd,capture_output=True,text=True,timeout=timeout)
    (out/(name+'.stdout')).write_text(r.stdout);(out/(name+'.stderr')).write_text(r.stderr)
    assert r.returncode==0,(name,r.returncode,r.stderr)
    return json.loads(r.stdout)
source={'commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'milp_tree':subprocess.check_output(['git','rev-parse','HEAD:milp'],cwd=root,text=True).strip(),'clean_before_and_after_build':True,'binaries':{side:{'path':str(path),'sha256':sha(path),'features':['sparse_row_products','binary_coefficient_strengthening']+(['proof_profile'] if side=='on' else [])} for side,path in binaries.items()}}
small=root/'milp/instances/knapsack-10.milp'
args=['solve',str(small),'--early-stop','off','--probe-reuse','off','--binary-propagation','auto','--threads','1','--time-limit-secs','0','--node-limit','100']
models={side:invoke('small-export-'+side,[str(binary),'export',str(small)]) for side,binary in binaries.items()}
assert canonical_hash(models['off'])==canonical_hash(models['on'])
rows={side:invoke('small-'+side,[str(binary),*args]) for side,binary in binaries.items()}
keys=sorted(set(rows['off'])-{'elapsed_secs','proof_profile','dense_optimal_profile'})
diff={key:[rows['off'][key],rows['on'][key]] for key in keys if rows['off'][key]!=rows['on'][key]}
errors=[]
for side in rows:verify_solution(models[side],rows[side],errors.append)
write('02-off-on-parity.json',{'commands':{side:[str(binary),*args] for side,binary in binaries.items()},'time_limit_disabled':True,'node_limit':100,'non_time_field_count':len(keys),'fields':keys,'differences':diff,'model_sha256':canonical_hash(models['on']),'original_model_errors':errors,'natural_optimal_calls':rows['on']['dense_optimal_profile']['natural_optimal_calls']})
assert not diff and not errors
basefile=main/'bench/large_scale_rounds/results/dense-recovery-probes/p0201.jsonl';base=json.loads(basefile.read_text());original=json.loads((basefile.parent/base['model']).read_text())
command=[str(binaries['on']),*base['command'][1:]];command[command.index('--node-limit')+1]='1'
assert command[command.index('--time-limit-secs')+1]=='30'
assert sha(command[2])==base['milp_sha256']
export_command=[str(binaries['on']),'export',command[2]]
model=invoke('p0201-root-export',export_command)
if 'mps_columns' in original:model['mps_columns']=original['mps_columns']
digest=canonical_hash(model);assert digest==base['model_sha256']
(out/(digest+'.json')).write_text(json.dumps(model,ensure_ascii=False,allow_nan=False,sort_keys=True,separators=(',',':'))+'\n')
result=invoke('p0201-root',command)
profile=result['dense_optimal_profile']
summary={k:v for k,v in profile.items() if k!='samples'}
summary['sample_metadata']=[{k:v for k,v in sample.items() if k in ('optimal_ordinal','anomaly_categories_mask','iter_call','iterations_total','num_structural','num_static','num_art','snapshot_complete','estimated_bytes','augmented_reconstruction_available','internal_box_lower_bound_kind','basis_artificial_count')} for sample in profile['samples']]
write('p0201-root-summary.json',summary)
ledger={k:v for k,v in result.items() if k not in ('proof_profile','dense_optimal_profile')}
for key in ('schema','instance','group','threads','model_sha256','milp_sha256','mps_sha256'):ledger[key]=base[key]
ledger.update(solver='mic',rep=0,variant='dense-optimal-root-diagnostic',time_limit_secs=30,node_limit=1,model=digest+'.json',command=command,source_commit=source['commit'],source_tree=source['milp_tree'],binary_sha256=source['binaries']['on']['sha256'],performance_promotion_evidence=False,diagnostic_artifact='p0201-root.stdout')
(out/'p0201-root-result.jsonl').write_text(json.dumps(ledger,ensure_ascii=False,allow_nan=False)+'\n')
refs=[]
for line in (main/'bench/large_scale_rounds/results/external-baseline.jsonl').read_text().splitlines():
    r=json.loads(line)
    if r['instance']=='p0201' and r['rep']==0 and r['solver'] in ('highs','scip'):
        assert r['model_sha256']==digest;r['model']=digest+'.json';refs.append(r)
(out/'p0201-external.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=False,allow_nan=False)+'\n' for r in refs))
checked=inspect_results([out/'p0201-root-result.jsonl',out/'p0201-external.jsonl'])
write('p0201-root-independent-check.json',checked)
source.update(command=command,baseline_command=base['command'],export_command=export_command,export_model_sha256=digest,model_sha256=base['model_sha256'],milp_sha256=sha(command[2]),performance_promotion_evidence=False,diagnostic_artifact_sha256=sha(out/'p0201-root.stdout'))
write('run-manifest.json',source)
print(json.dumps({'parity_differences':diff,'parity_non_time_fields':len(keys),'status':result['status'],'elapsed_secs':result['elapsed_secs'],'nodes':result['nodes'],'lp_iterations':result['lp_iterations'],'cuts':result['cuts'],'natural_optimal_calls':profile['natural_optimal_calls'],'anomalous_calls':profile['anomalous_calls'],'samples':len(profile['samples']),'independent_check_passed':checked['passed']},ensure_ascii=False),flush=True)
