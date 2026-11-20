import json, pathlib, subprocess, os, hashlib, sys
root=pathlib.Path('/tmp/mic-dense-failure-profile')
out=root/'bench/large_scale_rounds/results/dense-failure-profile'
binary=pathlib.Path('/tmp/mic-dense-failure-binaries/round_probe_profile')
basefile=pathlib.Path('/Users/guyuemochen/Desktop/projects/MIC/bench/large_scale_rounds/results/dense-recovery-probes/lseu.jsonl')
base=json.loads(basefile.read_text())
original=json.loads((basefile.parent/base['model']).read_text())
def sha(p): return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def write(name,value): (out/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
command=[str(binary),*base['command'][1:]]
export_command=[str(binary),'export',base['command'][2]]
assert sha(base['command'][2])==base['milp_sha256']
export=subprocess.run(export_command,capture_output=True,text=True,check=True)
(out/'export.stderr').write_text(export.stderr)
model=json.loads(export.stdout)
if 'mps_columns' in original: model['mps_columns']=original['mps_columns']
model_bytes=(json.dumps(model,ensure_ascii=False,allow_nan=False,sort_keys=True,separators=(',',':'))+'\n').encode()
digest=hashlib.sha256(model_bytes).hexdigest()
assert digest==base['model_sha256']
(out/(digest+'.json')).write_bytes(model_bytes)
snapshot=out/'lseu-failures.jsonl'
assert not snapshot.exists()
env=os.environ.copy();env['MIC_DENSE_FAILURE_PROFILE']=str(snapshot)
run=subprocess.run(command,capture_output=True,text=True,env=env,timeout=90)
(out/'lseu.stdout').write_text(run.stdout);(out/'lseu.stderr').write_text(run.stderr)
assert run.returncode==0,run.returncode
result=json.loads(run.stdout)
keys=sorted(set(result)-{'elapsed_secs','proof_profile'})
differences={key:{'baseline':base.get(key),'diagnostic':result[key]} for key in keys if base.get(key)!=result[key]}
info={'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'source_milp_tree':subprocess.check_output(['git','rev-parse','HEAD:milp'],cwd=root,text=True).strip(),'binary_sha256':sha(binary),'features':['sparse_row_products','binary_coefficient_strengthening','proof_profile'],'command':command,'baseline_command':base['command'],'export_command':export_command,'export_model_sha256':digest,'milp_sha256':sha(base['command'][2]),'non_time_fields':keys,'non_time_field_count':len(keys),'differences':differences,'performance_promotion_evidence':False,'returncode':run.returncode}
write('run-manifest.json',info);write('lseu-result.json',result)
print(json.dumps({'elapsed_secs':result['elapsed_secs'],'non_time_field_count':len(keys),'differences':differences,'failures':sum(1 for _ in snapshot.open())},ensure_ascii=False),flush=True)
