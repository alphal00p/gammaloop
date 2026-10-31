import sys,hashlib,json
from pathlib import Path
O=Path(__file__).resolve().parent
sha=lambda b:hashlib.sha256(b).hexdigest()
def write(name,obj):(O/name).write_text(json.dumps(obj,indent=2)+'\n')
r=json.loads((O/'review.json').read_bytes());f=json.loads((O/'integrity-final.json').read_bytes());before=json.loads((O/'integrity-before.json').read_bytes());controls=json.loads((O/'fresh-controls.json').read_bytes())
assert r['decision']=='APPROVE' and len(r['per_row_verdicts'])==3 and all(x['decision']=='APPROVE' for x in r['per_row_verdicts'])
assert not r['blocking_findings'] and not f['failures'] and all(f['checks'][p]==v for p,v in before['checks'].items())
assert (controls['total'],controls['accepted'],controls['rejected'])==(26,4,22)
assert not list(O.rglob('__pycache__'))
files=sorted(p for p in O.rglob('*') if p.is_file() and p.name not in ('SHA256SUMS','SEAL.json','SEAL.sha256'))
assert all(not p.is_symlink() for p in files)
manifest=''.join(f'{sha(p.read_bytes())}  {p.relative_to(O)}\n' for p in files).encode();(O/'SHA256SUMS').write_bytes(manifest)
seal={'schema':'independent-next-three-raw-qcd-migration-review-seal/1','decision':'APPROVE_EXACT_THREE_TOKENS_ONLY','files':len(files),'manifest_sha256':sha(manifest),'review_sha256':sha((O/'review.json').read_bytes()),'findings_sha256':sha((O/'findings.json').read_bytes()),'fresh_controls_sha256':sha((O/'fresh-controls.json').read_bytes()),'candidate_seal_sha256':r['candidate']['SEAL_sha256'],'candidate_patch_sha256':r['candidate']['patch_sha256'],'root_scientific_acceptance_sha256':r['additive_scientific_join']['root_acceptance_sha256'],'scientific_review_sha256':r['additive_scientific_join']['review_sha256'],'target':r['candidate']['target'],'applied_to_tracked_source':False,'postapplication_runtime_approved':False,'broader_integration_approved':False,'external_dependencies_required':True,'owned_processes_remaining':0}
write('SEAL.json',seal);sealhash=sha((O/'SEAL.json').read_bytes());(O/'SEAL.sha256').write_text(sealhash+'  SEAL.json\n')
for line in manifest.decode().splitlines():
 h,n=line.split('  ',1);assert sha((O/n).read_bytes())==h
print(json.dumps({'review_directory':str(O),'review_sha256':seal['review_sha256'],'SEAL_sha256':sealhash,'SHA256SUMS_sha256':seal['manifest_sha256'],'files':len(files),'final_external_files_verified':f['unique_files'],'refs_equal':r['validation']['Git_refs_before_after_equal'],'owned_processes':0},indent=2))
