from pathlib import Path
import json,hashlib,shutil,datetime
P=Path(__file__).resolve().parent;W=P.parents[2];D=Path(json.loads((P/'evidence-directory.json').read_text())['path'])
def sha(p):return hashlib.file_digest(Path(p).open('rb'),'sha256').hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
old=P/'evidence-seal-attempt1';old.mkdir()
for n in ('SHA256SUMS','SEAL.json'):shutil.copyfile(D/n,old/n)
large=[]
for name in ('initial-candidate-files.json','artifacts/ordinary/compiled-source-bindings.json'):
 f=D/'qualification'/name;local=P/name;assert f.stat().st_size>1048576 and sha(f)==sha(local)
 large.append({'original_evidence_path':str(f.relative_to(D)),'retained_path':str(local),'sha256':sha(local),'bytes':local.stat().st_size,'purpose':'Exact guard/compiled-source provenance retained for independent review; omitted from Git under unchanged 1MiB new-file limit.'});f.unlink()
save(D/'large-evidence-bindings.json',{'files':large,'partial_evidence_revision':'b6c3886feb702cbf120ec448356684ec5c9b7857','prior_seal_retained':str(old),'existing_tracking_limit_unchanged':True,'not_a_producer_or_runtime_failure':True})
ret=json.loads((D/'retained-artifacts.json').read_text());ret['large_evidence']=large;save(D/'retained-artifacts.json',ret)
r=json.loads((D/'receipt.json').read_text());r['evidence_packaging_correction']={'partial_revision_preserved':'b6c3886feb702cbf120ec448356684ec5c9b7857','initial_tracking_warning':'JJ refused two files above existing 1MiB new-file limit; initial candidate seal referenced those not-yet-tracked files.','correction':'Exact large files retained locally by path/size/hash; final Git packet references them in large-evidence-bindings.json. Existing limits unchanged.','producer_runtime_or_source_result_changed':False};save(D/'receipt.json',r)
f=json.loads((P/'retained-failures.json').read_text());f['retained_attempts'].append({'id':'evidence-track-attempt1','partial_revision':'b6c3886feb702cbf120ec448356684ec5c9b7857','error':'JJ refused two >1MiB files; nothing changed in those two file tracking entries.','retained_files':large,'prior_seal':str(old),'correction':'Add exact local hash bindings; no tracking limit/budget/source/build/runtime change.'});save(P/'retained-failures.json',f);shutil.copyfile(P/'retained-failures.json',D/'qualification/retained-failures.json');shutil.copyfile(Path(__file__),D/'qualification/fix_evidence_index.py')
n=W/'docs/architecture/generation-host-stub-rustred.typ';s=n.read_text();s=s.replace('retained locally for review, as listed in `retained-artifacts.json`.','retained locally for review, as listed in `retained-artifacts.json`. Two large\nguard/provenance JSON inputs are hash-bound by `large-evidence-bindings.json`;\nthe partial JJ tracking candidate and its prior seal remain preserved. Existing\ntracking limits are unchanged.');n.write_text(s)
c=json.loads((P/'canonical-structure-check.json').read_text());c['source_sha256']=sha(n);save(P/'canonical-structure-check.json',c);shutil.copyfile(P/'canonical-structure-check.json',D/'qualification/canonical-structure-check.json')
manifest={str(f.relative_to(D)):sha(f) for f in sorted(D.rglob('*')) if f.is_file() and f.name not in ('SHA256SUMS','SEAL.json')};assert all((D/k).stat().st_size<=1048576 for k in manifest)
(D/'SHA256SUMS').write_text(''.join(h+'  '+k+'\n' for k,h in manifest.items()));save(D/'SEAL.json',{'schema_version':1,'sealed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'code':'2cfb52b28d99f92305a9d571bbe3911cf04b8a19','files':len(manifest),'manifest_sha256':sha(D/'SHA256SUMS'),'receipt_sha256':sha(D/'receipt.json'),'scope':'Evidence integrity only; NEW independent review required, no self approval.','large_input_bindings':'large-evidence-bindings.json','prior_partial_candidate':'b6c3886feb702cbf120ec448356684ec5c9b7857','prior_seal':str(old),'owned_processes_remaining':0})
for f,h in manifest.items():assert sha(D/f)==h
save(P/'report.json',r);print(json.dumps({'files':len(manifest),'receipt_sha256':sha(D/'receipt.json'),'seal_sha256':sha(D/'SEAL.json'),'large_files':large}))
