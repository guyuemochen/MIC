from pathlib import Path
import sys,json,subprocess,hashlib,shutil,math
root=Path('/tmp/mic-safe-primal-ratio')
original=Path('/Users/guyuemochen/Desktop/projects/MIC')
sys.path.insert(0,str(root/'bench/large_scale_rounds'))
from compare import invoke,sha,store_model
from check_results import inspect_results
out=root/'results/safe-primal-ratio-validation'
base_path=original/'bench/large_scale_rounds/results/root_heuristic_work_limit-screen.jsonl'
bases=[json.loads(line) for line in base_path.read_text().splitlines()]
base=next(r for r in bases if r['instance']=='p0033' and r['variant']=='base')
def git(*args):return subprocess.check_output(['git',*args],cwd=root,text=True).strip()
assert git('status','--porcelain')=='', '构建与探针须使用干净源提交'
commit,tree=git('rev-parse','HEAD'),git('rev-parse','HEAD^{tree}')
binary=root/'milp/target/release/examples/round_probe'
binsha=sha(binary)
command=[str(binary)]+base['command'][1:]
assert command[1:]==base['command'][1:]
instance=Path(command[2])
assert sha(instance)==base['milp_sha256']
archived=json.loads((base_path.parent/base['model']).read_text())
fresh=invoke([str(binary),'export',str(instance)],180)
if 'mps_columns' in archived:fresh['mps_columns']=archived['mps_columns']
(out/'models').mkdir(exist_ok=True)
model_path,modelsha=store_model(out/'models',fresh)
assert modelsha==base['model_sha256']
row=invoke(command,130,rss=True)
row.update({k:base[k] for k in ['schema','instance','group','solver','threads','time_limit_secs','node_limit','model_sha256','milp_sha256','mps_sha256','host']})
row.update(model=str(model_path.relative_to(out)),variant='safe-primal-ratio-current-source',rep=0,phase='correctness-probe',command=command,mic_options=command[3:],features=[],default_features=['sparse_row_products','binary_coefficient_strengthening'],source_commit=commit,source_tree=tree,source_clean_before_probe=True,binary_sha256=binsha,version=binsha,export_model_sha256=modelsha,baseline_file=str(base_path),baseline_file_sha256=sha(base_path),performance_promotion_evidence=False)
result_path=out/'16-p0033-current-source.jsonl'
with result_path.open('x') as f:f.write(json.dumps(row,ensure_ascii=False)+'\n')
external_path=original/'bench/large_scale_rounds/results/external-baseline.jsonl'
external=[json.loads(line) for line in external_path.read_text().splitlines()]
external=[r for r in external if r['instance']=='p0033' and r['solver'] in ['highs','scip'] and r['rep']==0]
assert len(external)==2
for r in external:
 assert r['model_sha256']==modelsha
 r['model']=row['model']
witness=out/'17-p0033-external-witnesses.jsonl'
with witness.open('x') as f:
 for r in external:f.write(json.dumps(r,ensure_ascii=False)+'\n')
check=inspect_results([result_path,witness])
check.update(external_source=str(external_path),external_source_sha256=sha(external_path))
(out/'18-p0033-original-check.json').write_text(json.dumps(check,ensure_ascii=False,indent=2)+'\n')
manifest=dict(source_commit=commit,source_tree=tree,source_clean_before_build_and_probe=True,build_command=['cargo','build','--manifest-path','milp/Cargo.toml','--release','--example','round_probe'],features=[],default_features=row['default_features'],binary=str(binary),binary_sha256=binsha,command=command,original_command=base['command'],only_binary_path_replaced=command[1:]==base['command'][1:],model_sha256=modelsha,milp_sha256=sha(instance),baseline_status=base['status'],baseline_objective=base['objective'],baseline_bound=base['bound'],baseline_nodes=base['nodes'],baseline_lp_iterations=base['lp_iterations'],performance_promotion_evidence=False)
(out/'19-p0033-manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
shutil.copyfile('/tmp/mic-safe-ratio-p0033-build.log',out/'20-p0033-build.log')
print(json.dumps({'result':row,'check':check},ensure_ascii=False,indent=2))
assert check['passed']
