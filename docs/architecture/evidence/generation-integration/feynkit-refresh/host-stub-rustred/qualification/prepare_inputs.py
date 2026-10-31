import os, pathlib, subprocess, hashlib, json, time
W=pathlib.Path.cwd(); P=W/'.local/feynkit-refresh/host-stub-rustred-integration-20261031-r1'
ENV={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')};ENV['GIT_OPTIONAL_LOCKS']='0'
def run(a): return subprocess.check_output(a,env=ENV,cwd=W)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(n,o): (P/n).write_text(json.dumps(o,indent=2)+'\n')
checks={'.local/feynkit-refresh/host-api-stub-correction-20261031-r1/candidate.patch':'0666eba15326649c0e1fb3dc035bfe0bc8cbc1238ec5fcec99622ecbf2d58e47','.local/review-host-browser-stub-astra-20261031-013106Z-5655cc61/review.json':'8bd31e92448a517479fb2f2ec5a12c9999c31752cea08502392f290c47e37627','.local/review-host-browser-stub-astra-20261031-013106Z-5655cc61/SEAL.json':'6e104cabc93cb839a498c8b79ecb048a17693707704238288619ff41930759d6','.local/review-raw-three-integration-astra-20261031-012811Z-323f4486/review.json':'535b792610c531b446d8beb1a9f300246e6d498e8126509b637189e8eb80a091','.local/feynkit-refresh/rustred-host-preflight-20261031-r1/report.json':'855ecafaced15a291a4427d1fea0f3f27ac936ca14b95c1edb4fb5151a68e735'}
for p,h in checks.items(): assert sha(W/p)==h,p
counts={}
for n in ('rustred-host-preflight-20261031-r1','host-api-stub-correction-20261031-r1','host-browser-correction-20261031-r1'):
 d=W/'.local/feynkit-refresh'/n; count=0
 for l in (d/'SHA256SUMS').read_text().splitlines():
  h,r=l.split('  ',1); f=d/r
  if f.is_symlink(): got=hashlib.sha256(os.readlink(f).encode()).hexdigest()
  else: got=sha(f)
  assert got==h,(n,r,got,h)
  count+=1
 counts[n]=count
write('verified-prior-packets.json',{'sha256':checks,'verified_entries':counts,'preflight_seal':json.loads((W/'.local/feynkit-refresh/rustred-host-preflight-20261031-r1/seal.json').read_text())})
refs=run(['jj','--ignore-working-copy','bookmark','list','--all','-T','name ++ " " ++ normal_target.commit_id() ++ "\\n"']).decode();(P/'initial-bookmarks.txt').write_text(refs)
(P/'initial-jj-op.txt').write_bytes(run(['jj','--ignore-working-copy','op','log','-n','1','--no-graph']))
commit=run(['jj','--ignore-working-copy','log','-r','@','--no-graph','-T','commit_id']).decode();assert commit=='d0fb4e906140a6a0fc8dd6401ca839cab0c2734d',commit
stub='examples/notebooks/symbolica-host/python/symbolica/community/hepkit/__init__.pyi'; assert sha(W/stub)=='bf9d4ce3f4b5f3871c39164d3c11db248123eecd8a1939a845ad7569da4de92f'
# Compare every existing candidate file against its committed Git blob without touching index.
def inventory(root,rev):
 entries={}
 raw=subprocess.check_output(['git','ls-tree','-rz',rev],cwd=W.parent/'main',env=ENV)
 for e in raw.split(b'\0'):
  if not e:continue
  header,name=e.split(b'\t',1); mode,kind,oid=header.decode().split(); r=os.fsdecode(name); f=root/r
  data=os.readlink(f).encode() if mode=='120000' else f.read_bytes()
  assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==oid,(str(root),r)
  entries[r]={'mode':mode,'oid':oid,'sha256':hashlib.sha256(data).hexdigest()}
 return entries
write('initial-candidate-files.json',inventory(W,commit));M=W.parent/'main';head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=M,env=ENV).decode().strip();write('initial-main-files.json',inventory(M,head))
gitdir=pathlib.Path(subprocess.check_output(['git','rev-parse','--absolute-git-dir'],cwd=M,env=ENV).decode().strip());write('initial-guard.json',{'source':commit,'main_head':head,'main_index_path':str(gitdir/'index'),'main_index_sha256':sha(gitdir/'index'),'global_integration':run(['jj','--ignore-working-copy','log','-r','generation-integration','--no-graph','-T','commit_id']).decode(),'author':'Lucien Huber <im@lcnbr.ch>','private_environment_removed_for_children':True,'started_epoch':time.time()})
print(json.dumps({'verified':counts,'source':commit,'main_head':head}))
