from pathlib import Path
import os,json,hashlib,shutil,subprocess,datetime,ast
P=Path(__file__).resolve().parent;W=P.parents[2];D=Path(json.loads((P/'evidence-directory.json').read_text())['path']);E={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')};E['GIT_OPTIONAL_LOCKS']='0'
def sha(p):return hashlib.file_digest(Path(p).open('rb'),'sha256').hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def copy(p,q):q=D/q;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
# Remove only our disposable empty mounts, home/config/state, install caches and temp browser profiles.
removed=[]
for lane in ('ordinary','private'):
 root=P/'runtime'/lane
 for parent in [root/'installed',*list((root/'work').iterdir())]:
  if not parent.is_dir():continue
  for name in ('cache','tmp','home','config','state','runtime'):
   f=parent/name
   if f.exists():shutil.rmtree(f);removed.append(str(f.relative_to(P)))
  if parent.parent.name=='work' and (parent/'venv').exists():
   assert not any((parent/'venv').iterdir());(parent/'venv').rmdir();removed.append(str((parent/'venv').relative_to(P)))
# Check known build/product processes without touching unrelated native work.
known=[]
for f in [P/'build-process.json',P/'build-r2-process.json',*list((P/'runtime').rglob('process.json'))]:
 d=json.loads(f.read_text());pid=d.get('pid');proc=Path('/proc')/str(pid)
 if proc.exists():
  args=(proc/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
  if 'host-stub-rustred' in args or 'bwrap' in args:known.append({'pid':pid,'argv':args})
assert not known,known
save(P/'cleanup.json',{'removed_own_disposable_directories':removed,'owned_build_or_runtime_processes_remaining':known,'retained':['immutable source/proposal/output store inputs','actual native extractions and two complete offline-installed private venvs','all command/process/stdout/stderr records, old and new failed setup/build attempts','exact prior review packets/seal bindings'],'reason':'Independent review and continuation need these exact inputs; no old packet, shared cache, contributor state, named ref or workspace removed.'})
# Record all changes to shared JJ bookmarks; do not revert unrelated Family Auto activity.
args=['jj','--ignore-working-copy','bookmark','list','--all','-T','name ++ " " ++ normal_target.commit_id() ++ "\\n"'];current=subprocess.check_output(args,cwd=W,env=E).decode();(P/'final-bookmarks-before-candidate.txt').write_text(current)
def parse(s):return dict(l.rsplit(' ',1) for l in s.splitlines() if l.strip())
a=parse((P/'initial-bookmarks.txt').read_text());b=parse(current);changes={k:{'old':a[k],'current':b.get(k)} for k in a if b.get(k)!=a[k]};save(P/'shared-ref-changes.json',{'initial_existing_bookmark_changes':changes,'new_bookmarks_before_final_evidence':{k:v for k,v in b.items() if k not in a},'not_reverted':True,'no_global_vcs_invariance_claim':True,'global_integration_preserved':json.loads((P/'invariants.json').read_text())['global_integration_unchanged']})
failures={'retained_attempts':[{'id':'guard-object-store-attempt1','exit_code':1,'error':'fatal: not a git repository (or any parent up to mount point /)\nStopping at filesystem boundary (GIT_DISCOVERY_ACROSS_FILESYSTEM not set).','cause':'Feature workspace has JJ metadata without a Git worktree; ls-tree object reads corrected to main Git store read-only.','script':'prepare_inputs.attempt1.py','source_or_runtime_started':False},{'id':'pure-graph-adapter-attempt1','error_log':'pure-stage.stderr','error':'KeyError for full /nix/store/...drv in Nix version4 derivations wrapper/basename-key mapping','script':'pure_stage.attempt1.py','correction':'finish_pure_stage.py adapter only; no producer input change.'},{'id':'private-locked-metadata-attempt1','exit_code':101,'logs':['build-completion.json','qualification/log-private/stdout'],'correction':'lock-reconciliation.json; no version/source/feature/semantic change; final locked gates and both builds pass.'}],'runtime_failures':0,'old_failures':'Exact prior packets and reviews retained unchanged; original ordinary QCD failure is not relabelled as a private-host result.'};save(P/'retained-failures.json',failures)
for f in ('invariants.json','canonical-structure-check.json','cleanup.json','shared-ref-changes.json','retained-failures.json','finalize-evidence.stdout','finalize-evidence.stderr','final-invariants.stdout','final-invariants.stderr','seal_packet.py'):
 copy(P/f,Path('qualification')/f)
# Exact old source seals remain verifiable after our work.
verified={}
for name in ('host-api-stub-correction-20261031-r1','host-browser-correction-20261031-r1','rustred-host-preflight-20261031-r1'):
 root=P.parent/name;count=0
 for l in (root/'SHA256SUMS').read_text().splitlines():
  h,r=l.split('  ',1);f=root/r;got=hashlib.sha256(os.readlink(f).encode()).hexdigest() if f.is_symlink() else sha(f);assert got==h,(name,r);count+=1
 verified[name]={'verified_entries':count,'manifest_sha256':sha(root/'SHA256SUMS')}
save(P/'prior-seals-final-verification.json',verified);copy(P/'prior-seals-final-verification.json',Path('qualification/prior-seals-final-verification.json'))
# Final syntax and readonly input/file hashes; no compiler/Nix operation here.
syntax={}
for f in P.glob('*.py'):ast.parse(f.read_text());syntax[f.name]=sha(f)
save(P/'syntax-final.json',{'python_AST':syntax,'no_runtime_or_compile_credit':True});copy(P/'syntax-final.json',Path('qualification/syntax-final.json'))
for f in D.rglob('*.json'):json.loads(f.read_text())
manifest={str(f.relative_to(D)):sha(f) for f in sorted(D.rglob('*')) if f.is_file() and f.name not in ('SHA256SUMS','SEAL.json')}
(D/'SHA256SUMS').write_text(''.join(h+'  '+p+'\n' for p,h in manifest.items()));save(D/'SEAL.json',{'schema_version':1,'sealed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'code':'2cfb52b28d99f92305a9d571bbe3911cf04b8a19','files':len(manifest),'manifest_sha256':sha(D/'SHA256SUMS'),'receipt_sha256':sha(D/'receipt.json'),'scope':'Candidate evidence integrity only; NEW independent review required, no self approval.','big_inputs':'retained-artifacts.json','all_original_packet_seals_verified':True,'licensed_runtime_finished':True,'owned_processes_remaining':0})
for rel,h in manifest.items():assert sha(D/rel)==h
print(json.dumps({'sealed':str(D),'files':len(manifest),'receipt_sha256':sha(D/'receipt.json'),'seal_sha256':sha(D/'SEAL.json'),'shared_ref_changes':changes,'removed_own_disposable_dirs':len(removed)}))
