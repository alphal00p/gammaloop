from pathlib import Path
import hashlib,json,re,os,subprocess
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh');O=Path(__file__).resolve().parent;C=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/grouped-identical-migration-b1-r3'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,indent=2)+'\n')
for p in C.rglob('*.json'):json.loads(p.read_text())
assert re.findall(r'^= (.+)$',(C/'record.typ').read_text(),re.M)==['Grouped Identical evidence durability correction']
initial=json.loads((O/'initial.json').read_text());assert sha(W/'.gitignore')==initial['gitignore_sha256'] and sha(W/'docs/developers.toml')==initial['registry_sha256']
assert sha(C/'review/review.json')==initial['review']['sha256'] and sha(C/'review/seal.json')==initial['seal']['sha256']
# Archive the exact sealer source before sealing; its outcomes live in ignored scratch.
(C/'seal_correction.py').write_bytes(Path(__file__).read_bytes())
files=sorted(p for p in C.rglob('*') if p.is_file() and p not in [C/'SHA256SUMS',C/'seal.json']);assert all(p.stat().st_size<=1048576 for p in files)
(C/'SHA256SUMS').write_text(''.join(sha(p)+'  '+str(p.relative_to(W))+'\n' for p in files))
receipt=json.loads((C/'receipt.json').read_text());save(C/'seal.json',{'schema':'grouped-identical-B1-correction-seal/1','status':'NEW_INDEPENDENT_FOCUSED_REREVIEW_REQUIRED','base':'7c0e13a08786416a9f24b3bec6df8a50e4592346','log_repair_commit':'ed669e2beeb39f36b3e5753ab9c746770ba7653d','old_manifest_sha256':receipt['unchanged_old_manifest_sha256'],'scope_manifest_entries':425,'scope_manifest_vector_sha256':json.loads((C/'manifest-audit.json').read_text())['scoped_manifest_tree_sha256'],'manifest_sha256':sha(C/'SHA256SUMS'),'receipt_sha256':sha(C/'receipt.json'),'record_sha256':sha(C/'record.typ'),'review_sha256':initial['review']['sha256'],'review_lowercase_seal_sha256':initial['seal']['sha256'],'files':len(files),'bytes':sum(p.stat().st_size for p in files),'self_approval':False,'excluded':['SHA256SUMS bound by seal','seal.json bound by external final handoff'],'final_commit_audit':'Exact final immutable commit plus this packet checked after snapshot; full425 vector must equal manifest-audit.json. Actual final ID is external to avoid circular commit hashing.'})
env={k:os.environ[k] for k in ['PATH','HOME','USER','LANG'] if k in os.environ}
p=subprocess.run(['jj','--config','snapshot.auto-track=none()','file','track',str(C.relative_to(W))],cwd=W,env=env,capture_output=True);(O/'track-seal.stdout').write_bytes(p.stdout);(O/'track-seal.stderr').write_bytes(p.stderr);save(O/'track-seal.command.json',{'argv':['jj','--config','snapshot.auto-track=none()','file','track',str(C.relative_to(W))],'exit':p.returncode,'stdout_sha256':hashlib.sha256(p.stdout).hexdigest(),'stderr_sha256':hashlib.sha256(p.stderr).hexdigest()});assert p.returncode==0,(p.stdout,p.stderr)
# An explicit snapshot precedes the immutable read, unlike the retained first attempt.
p=subprocess.run(['jj','--config','snapshot.auto-track=none()','describe','-m','record grouped migration evidence durability correction'],cwd=W,env=env,capture_output=True);(O/'snapshot-packet.stdout').write_bytes(p.stdout);(O/'snapshot-packet.stderr').write_bytes(p.stderr);assert p.returncode==0
rev=subprocess.check_output(['jj','--ignore-working-copy','log','--no-graph','-r','@','-T','commit_id'],cwd=W,env=env,text=True).strip();(O/'candidate-revision.txt').write_text(rev+'\n')
print('Correction candidate',rev,'seal',sha(C/'seal.json'),'files',len(files))
