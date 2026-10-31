"""Track bounded complete evidence, freeze two commits and verify every immutable blob."""
import datetime,gzip,hashlib,importlib.util,json,os,subprocess,sys
from pathlib import Path
O=Path(__file__).resolve().parent;W=O.parents[2];E=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/next-qcd-migration';D=Path(str(E)+'-durability')
CODE=json.loads((O/'code.json').read_text())['code'];BASE='83276dc047815d1dcc7f003aeb512c204b26758e'
ENV={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')};ENV.update(GIT_OPTIONAL_LOCKS='0',GIT_CONFIG_NOSYSTEM='1',GIT_CONFIG_GLOBAL='/dev/null')
commands=[]
def sha(b):return hashlib.sha256(b).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2)+'\n')
def run(args):
 r=subprocess.run(args,cwd=W,env=ENV,capture_output=True);commands.append({'argv':args,'exit':r.returncode,'stdout':r.stdout.decode(),'stderr':r.stderr.decode()});write(O/'freeze-commands.json',commands);assert r.returncode==0,(args,r.stderr.decode());return r.stdout.decode().strip()
def jj(*a):return run(['jj','--config','user.name="Lucien Huber"','--config','user.email="im@lcnbr.ch"',*a])
def git(*a):return run(['git','--git-dir='+str(W.parent/'main/.git'),*a])
def track(files):
 assert all(p.stat().st_size<=1048576 for p in files)
 for i in range(0,len(files),150):jj('file','track','--include-ignored',*[str(p.relative_to(W)) for p in files[i:i+150]])
assert jj('log','--no-graph','-r','@-','-T','commit_id')==CODE
assert jj('log','--no-graph','-r','generation-integration','-T','commit_id')==BASE
cap=jj('config','get','snapshot.max-new-file-size');assert cap in ['1048576','1MiB','"1MiB"'],cap
files=[W/line.split('  ',1)[1] for line in (E/'SHA256SUMS').read_text().splitlines()]+[E/'SHA256SUMS',E/'SEAL.json'];assert len(files)==len(set(files));track(files)
jj('metaedit','--author','Lucien Huber <im@lcnbr.ch>','-m','freeze postapplication QCD migration qualification')
evidence=jj('log','--no-graph','-r','@','-T','commit_id');jj('bookmark','create','generation-next-qcd-migration-evidence-20261031','-r',evidence)
git('diff','--check',CODE,evidence,'--','docs/developers.toml','docs/architecture/generation-next-qcd-migration.typ')
spec=importlib.util.spec_from_file_location('committed',O/'verify_committed.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
audit=v.audit(evidence);write(O/'evidence-immutable-audit.json',audit)
assert all(not p[k] for p in audit['packets'] for k in ['missing','mismatched','unexpected'])
jj('new','-m','record immutable QCD evidence closure for independent review')
D.mkdir(exist_ok=False)
b=(O/'evidence-immutable-audit.json').read_bytes();z=gzip.compress(b,compresslevel=9,mtime=0);assert len(z)<=1048576;(D/'evidence-immutable-audit.json.gz').write_bytes(z)
# Bind preserved external inputs again after runtime. Mutable unrelated source/config is never repaired.
original=json.loads((O/'input-verification.json').read_text())['unique_verified_files'];moved=[]
for p,r in original.items():
 f=Path(p)
 if not f.is_file() or sha(f.read_bytes())!=r['sha256']:moved.append({'path':p,'expected':r,'actual_sha256':sha(f.read_bytes()) if f.is_file() else None})
assert not moved,moved
start=json.loads((O/'start.json').read_text());main=W.parent/'main';index=main/'.git/index'
end={'utc':datetime.datetime.now(datetime.UTC).isoformat(),'JJ_refs':jj('bookmark','list','--all'),'Git_refs':git('show-ref'),'main_head':git('rev-parse','HEAD'),'main_index_sha256':sha(index.read_bytes()),'main_status':run(['git','-C',str(main),'status','--porcelain=v1','--untracked-files=no']),'main_index_unchanged':sha(index.read_bytes())==start['main_index_sha256'],'main_status_unchanged':run(['git','-C',str(main),'status','--porcelain=v1','--untracked-files=no'])==start['main_status'].strip(),'all_original_input_bindings_unchanged':len(original),'global_integration':jj('log','--no-graph','-r','generation-integration','-T','commit_id'),'other_contributor_movements':'Compare exact before/after refs; no restoration or overwrite. Existing core union untouched.'}
write(O/'before-candidate-end-state.json',end);write(D/'before-candidate-end-state.json',end)
write(D/'tracking-and-freeze-commands.json',commands)
write(D/'receipt.json',{'schema':'next-QCD-immutable-evidence-closure/1','status':'NEW_INDEPENDENT_INTEGRATION_REVIEW_REQUIRED','code':CODE,'evidence_commit':evidence,'candidate_bookmark':'generation-next-qcd-migration-candidate-20261031','evidence_manifest_sha256':sha((E/'SHA256SUMS').read_bytes()),'evidence_seal_sha256':sha((E/'SEAL.json').read_bytes()),'immutable_blob_audit':{'compressed':'evidence-immutable-audit.json.gz','gzip_sha256':sha(z),'reconstructed_sha256':sha(b),'reconstructed_bytes':len(b),'matched':audit['packets'][0]['matched'],'missing':[],'mismatched':[],'unexpected':[],'lossless_reconstructions':audit['lossless_reconstructed_records']},'closure_protocol':'Primary manifest verified at exact evidence commit. This subordinate packet adds only immutable audit/record bytes. Final candidate is rechecked across BOTH exact manifests; final audit/handoff remains external to avoid circular self-commit claims. The committed verifier is replayable.','ignored_logs':'Explicit bounded jj file track --include-ignored on every primary manifest path; cap remains1MiB','global_integration':BASE,'global_advance':False,'self_approval':False,'external_dependencies_required':True})
subfiles=sorted(p for p in D.rglob('*') if p.is_file());(D/'SHA256SUMS').write_text(''.join(sha(p.read_bytes())+'  '+str(p.relative_to(W))+'\n' for p in subfiles));write(D/'SEAL.json',{'schema':'next-QCD-durability-seal/1','manifest_sha256':sha((D/'SHA256SUMS').read_bytes()),'files':len(subfiles),'evidence_commit':evidence,'primary_manifest_sha256':sha((E/'SHA256SUMS').read_bytes()),'primary_seal_sha256':sha((E/'SEAL.json').read_bytes()),'pending_independent_review':True})
track(sorted(p for p in D.rglob('*') if p.is_file()));jj('metaedit','--author','Lucien Huber <im@lcnbr.ch>','-m','record immutable QCD evidence closure for independent review')
candidate=jj('log','--no-graph','-r','@','-T','commit_id');jj('bookmark','create','generation-next-qcd-migration-candidate-20261031','-r',candidate)
final=v.audit(candidate);write(O/'candidate-immutable-audit.json',final);assert final['packets'][0]['vector_sha256']==audit['packets'][0]['vector_sha256'];assert len(final['packets'])==2
for rev in [CODE,evidence,candidate]:assert jj('log','--no-graph','-r',rev,'-T','author.name() ++ " <" ++ author.email() ++ ">"')=='Lucien Huber <im@lcnbr.ch>'
assert jj('log','--no-graph','-r','generation-integration','-T','commit_id')==BASE
assert not jj('resolve','--list','-r',candidate)
jj('new','-m','await new independent QCD migration integration review')
h={'schema':'next-QCD-postapplication-handoff/1','status':'NEW_INDEPENDENT_INTEGRATION_REVIEW_REQUIRED','code':CODE,'evidence':evidence,'candidate':candidate,'bookmarks':{'generation-next-qcd-migration-code-20261031':CODE,'generation-next-qcd-migration-evidence-20261031':evidence,'generation-next-qcd-migration-candidate-20261031':candidate},'receipt':str(E/'receipt.json'),'receipt_sha256':sha((E/'receipt.json').read_bytes()),'manifest':str(E/'SHA256SUMS'),'manifest_sha256':sha((E/'SHA256SUMS').read_bytes()),'seal':str(E/'SEAL.json'),'seal_sha256':sha((E/'SEAL.json').read_bytes()),'durability_receipt':str(D/'receipt.json'),'durability_seal_sha256':sha((D/'SEAL.json').read_bytes()),'closure':{'audit':str(O/'candidate-immutable-audit.json'),'audit_sha256':sha((O/'candidate-immutable-audit.json').read_bytes()),'packets':[{k:p[k] for k in ['root','manifest_entries','matched','missing','mismatched','unexpected','vector_sha256']} for p in final['packets']],'lossless_reconstructed_records':final['lossless_reconstructed_records']},'global_integration':BASE,'global_advance':False,'publication':False,'self_approval':False,'workspace':str(W),'local_artifacts':str(O),'main_and_inputs':end,'remaining_gates':json.loads((E/'receipt.json').read_text())['remaining_gates']}
write(O/'handoff.json',h);print(json.dumps(h,indent=2))
