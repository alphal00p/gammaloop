import pathlib,json,hashlib,datetime
D=pathlib.Path(__file__).resolve().parent

def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
assert not(D/'SEAL.json').exists()and not(D/'SHA256SUMS').exists()
review=json.loads((D/'review.json').read_text());assert review['decision']=='APPROVE'and len(review['per_case_mode_verdicts'])==6
files=sorted(p for p in D.rglob('*')if p.is_file());assert not any(p.is_symlink()for p in files)
manifest=''.join(sha(p)+'  '+str(p.relative_to(D))+'\n'for p in files);(D/'SHA256SUMS').write_text(manifest)
seal={'schema':'independent-next-qcd-scientific-finite-interface-review-seal/1','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'decision':'APPROVE','accepted_scope':review['accepted_scope'],'manifest_sha256':sha(D/'SHA256SUMS'),'bound_files':len(files),'review_sha256':sha(D/'review.json'),'candidate_complete_seal_sha256':review['candidate']['pins']['COMPLETE-SEAL.json'],'original_input_integrity_final_sha256':sha(D/'integrity-final.json'),'no_new_native_execution':True,'no_migration_or_runtime_application_approval':True,'owned_background_processes_remaining':0}
(D/'SEAL.json').write_text(json.dumps(seal,indent=2)+'\n')
for line in (D/'SHA256SUMS').read_text().splitlines():
 h,p=line.split('  ',1);assert sha(D/p)==h
assert sha(D/'SHA256SUMS')==seal['manifest_sha256']
print(json.dumps({'review_directory':str(D),'decision':review['decision'],'review_sha256':seal['review_sha256'],'seal_sha256':sha(D/'SEAL.json'),'manifest_sha256':seal['manifest_sha256'],'bound_files':len(files)},indent=2))
