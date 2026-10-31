from pathlib import Path
import hashlib,json,subprocess,os,time
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh')
P=W/'.local/feynkit-refresh/combined-integration-preflight-20261031-r1'
def read(path):
 q=Path(path); q=q if q.is_absolute() else W/q
 try:
  b=q.read_bytes()
  record={'path':str(q),'sha256':hashlib.sha256(b).hexdigest(),'size':len(b)}
 except Exception as e:
  record={'path':str(q),'error':str(e)}; b=None
 a=json.loads((P/'read-inventory.json').read_text()); a.append(record); (P/'read-inventory.json').write_text(json.dumps(a,indent=2)+'\n')
 if b is None: raise OSError(record['error'])
 return b.decode()
def data(path): return json.loads(read(path))
def run(label,argv,cwd=W,timeout=60):
 env={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')}
 t=time.time()
 try:
  r=subprocess.run(argv,cwd=cwd,env=env,capture_output=True,text=True,timeout=timeout)
  v={'argv':argv,'cwd':str(cwd),'exit_code':r.returncode,'stdout':r.stdout,'stderr':r.stderr,'elapsed_seconds':time.time()-t}
 except subprocess.TimeoutExpired as e:
  v={'argv':argv,'cwd':str(cwd),'error':'timeout','stdout':str(e.stdout),'stderr':str(e.stderr)}
 v['license_variables_omitted']=['SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE']
 (P/(label+'.json')).write_text(json.dumps(v,indent=2)+'\n'); return v
