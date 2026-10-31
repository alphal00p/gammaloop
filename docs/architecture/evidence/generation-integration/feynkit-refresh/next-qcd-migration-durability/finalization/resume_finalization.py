"""Preserve the frozen candidate and close the benign JJ conflict-list diagnostic."""
import datetime,gzip,hashlib,importlib.util,json,os,subprocess
from pathlib import Path
O=Path(__file__).resolve().parent;W=O.parents[2];E=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/next-qcd-migration';D=Path(str(E)+'-durability')
CODE=json.loads((O/'code.json').read_text())['code'];BASE='83276dc047815d1dcc7f003aeb512c204b26758e';OLD='d24fc097a25fd63dcb025371b5a2fc9a7587ac31'
ENV={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')};ENV.update(GIT_OPTIONAL_LOCKS='0',GIT_CONFIG_NOSYSTEM='1',GIT_CONFIG_GLOBAL='/dev/null');commands=[]
def sha(b):return hashlib.sha256(b).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,indent=2)+'\n')
def run(args):
 r=subprocess.run(args,cwd=W,env=ENV,capture_output=True);commands.append({'argv':args,'exit':r.returncode,'stdout':r.stdout.decode(),'stderr':r.stderr.decode()});save(O/'resume-finalization-commands.json',commands);assert r.returncode==0,(args,r.stderr.decode());return r.stdout.decode().strip()
def jj(*a):return run(['jj','--config','user.name="Lucien Huber"','--config','user.email="im@lcnbr.ch"',*a])
assert jj('log','--no-graph','-r','@','-T','commit_id')==OLD
assert jj('diff','--from',OLD,'--to','@','--name-only')==''
evidence=jj('log','--no-graph','-r','generation-next-qcd-migration-evidence-20261031','-T','commit_id')
assert not jj('log','--no-graph','-r','conflicts() & ('+CODE+' | '+evidence+' | '+OLD+')','-T','commit_id')
spec=importlib.util.spec_from_file_location('v',O/'verify_committed.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v);old=v.audit(OLD)
oldD={name:(D/name).read_bytes() for name in ['SHA256SUMS','SEAL.json','receipt.json']}
jj('new','-m','retain QCD finalization diagnostic and exact candidate closure')
F=D/'finalization';F.mkdir(exist_ok=False)
for name,b in oldD.items():(F/('original-candidate-'+name)).write_bytes(b)
raw=json.dumps(old,indent=2).encode()+b'\n';compressed=gzip.compress(raw,compresslevel=9,mtime=0);assert len(compressed)<=1048576;(F/'original-candidate-immutable-audit.json.gz').write_bytes(compressed)
for name in ['freeze_candidate.stdout','freeze_candidate.stderr','freeze-commands.json','resume_finalization.py','resume_finalization.before-diff-correction.py','resume-attempt1-resume_finalization.stdout','resume-attempt1-resume_finalization.stderr','resume-attempt1-resume-finalization-commands.json']:(F/name).write_bytes((O/name).read_bytes())
original_commands=json.loads((O/'freeze-commands.json').read_text());failed=original_commands[-1];assert failed['argv'][-4:]==['resolve','--list','-r',OLD]
assert failed['exit']!=0 and failed['stderr']=='Error: No conflicts found at this revision\n'
save(F/'resolution.json',{'schema':'QCD-finalization-only-diagnostic-resolution/1','original_candidate':OLD,'original_candidate_bookmark':'generation-next-qcd-migration-candidate-20261031','original_evidence':evidence,'original_probe':failed,'resume_setup_diagnostic':'First resume used jj diff default(parent-to-current) as a clean-tree check; it stopped before mutation. Replaced with explicit OLD-to-@ comparison. Original script/commands/error retained.','meaning':'JJ resolve --list returns a nonzero code for a conflict-free revision. Original immutable source/manifest audits had already succeeded. This is a finalization helper assumption error, not a source/build/native selected-test failure.','corrected_check':commands[3] if len(commands)>3 else commands,'all_code_evidence_candidate_conflicts_revset_empty':True,'original_primary_manifest_vector_sha256':old['packets'][0]['vector_sha256'],'original_durability_manifest_sha256':old['packets'][1]['manifest_sha256'],'original_audit_gzip_sha256':sha(compressed),'original_audit_reconstructed_sha256':sha(raw),'original_audit_reconstructed_bytes':len(raw),'new_builds':0,'new_native_selected_or_proof_runs':0,'ref_overwrites':0,'published':False})
r=json.loads((D/'receipt.json').read_text());r.update(candidate_bookmark='generation-next-qcd-migration-candidate-r2-20261031',preserved_original_candidate=OLD,finalization_record='finalization/resolution.json',finalization_scope='Retain diagnostic and use an exact empty conflict revset. No source/runtime/proof change; old candidate and old subordinate manifest remain immutable at their named revision.');save(D/'receipt.json',r)
save(F/'corrected-pre-freeze-commands.json',commands)
files=sorted(p for p in D.rglob('*') if p.is_file() and p not in [D/'SHA256SUMS',D/'SEAL.json']);assert all(p.stat().st_size<=1048576 for p in files)
(D/'SHA256SUMS').write_text(''.join(sha(p.read_bytes())+'  '+str(p.relative_to(W))+'\n' for p in files));s=json.loads((D/'SEAL.json').read_text());s.update(manifest_sha256=sha((D/'SHA256SUMS').read_bytes()),files=len(files),preserved_original_candidate=OLD,finalization_only=True);save(D/'SEAL.json',s)
files=sorted(p for p in D.rglob('*') if p.is_file())
for i in range(0,len(files),150):jj('file','track','--include-ignored',*[str(p.relative_to(W)) for p in files[i:i+150]])
jj('metaedit','--author','Lucien Huber <im@lcnbr.ch>','-m','retain QCD finalization diagnostic and exact candidate closure')
candidate=jj('log','--no-graph','-r','@','-T','commit_id');jj('bookmark','create','generation-next-qcd-migration-candidate-r2-20261031','-r',candidate);jj('bookmark','create','generation-next-qcd-migration-finalization-20261031','-r',candidate)
final=v.audit(candidate);save(O/'candidate-r2-immutable-audit.json',final);assert final['packets'][0]['vector_sha256']==old['packets'][0]['vector_sha256']
assert not jj('log','--no-graph','-r','conflicts() & ('+CODE+' | '+evidence+' | '+candidate+')','-T','commit_id')
for rev in [CODE,evidence,OLD,candidate]:assert jj('log','--no-graph','-r',rev,'-T','author.name() ++ " <" ++ author.email() ++ ">"')=='Lucien Huber <im@lcnbr.ch>'
assert jj('log','--no-graph','-r','generation-integration','-T','commit_id')==BASE
jj('new','-m','await new independent QCD migration integration review')
end=json.loads((O/'before-candidate-end-state.json').read_text())
h={'schema':'next-QCD-postapplication-handoff/1','status':'NEW_INDEPENDENT_INTEGRATION_REVIEW_REQUIRED','code':CODE,'evidence':evidence,'candidate':candidate,'bookmarks':{'generation-next-qcd-migration-code-20261031':CODE,'generation-next-qcd-migration-evidence-20261031':evidence,'generation-next-qcd-migration-candidate-r2-20261031':candidate,'generation-next-qcd-migration-finalization-20261031':candidate},'preserved_original_candidate':{'commit':OLD,'bookmark':'generation-next-qcd-migration-candidate-20261031','scope':'Immutable complete primary evidence and original durability packet; finalization helper exit2 retained, then closed without native reruns'},'receipt':str(E/'receipt.json'),'receipt_sha256':sha((E/'receipt.json').read_bytes()),'manifest':str(E/'SHA256SUMS'),'manifest_sha256':sha((E/'SHA256SUMS').read_bytes()),'seal':str(E/'SEAL.json'),'seal_sha256':sha((E/'SEAL.json').read_bytes()),'durability_receipt':str(D/'receipt.json'),'durability_seal':str(D/'SEAL.json'),'durability_seal_sha256':sha((D/'SEAL.json').read_bytes()),'closure':{'audit':str(O/'candidate-r2-immutable-audit.json'),'audit_sha256':sha((O/'candidate-r2-immutable-audit.json').read_bytes()),'packets':[{k:p[k] for k in ['root','manifest_entries','matched','missing','mismatched','unexpected','vector_sha256']} for p in final['packets']],'lossless_reconstructed_records':final['lossless_reconstructed_records']},'global_integration':BASE,'global_advance':False,'publication':False,'self_approval':False,'workspace':str(W),'local_artifacts':str(O),'main_and_inputs':end,'remaining_gates':json.loads((E/'receipt.json').read_text())['remaining_gates']}
save(O/'handoff.json',h);print(json.dumps(h,indent=2))
