"""Verify complete exact packet and originals before explicit bounded tracking."""
from pathlib import Path
import gzip,hashlib,json,os,subprocess
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh');O=Path(__file__).resolve().parent;E=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/grouped-identical-migration'
KEYS=['SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE'];assert not any(k in os.environ for k in KEYS)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest={}
for line in (E/'SHA256SUMS').read_text().splitlines():
 digest,path=line.split('  ',1);assert sha(W/path)==digest,path;assert path not in manifest;manifest[path]=digest
seal=json.loads((E/'SEAL.json').read_text());assert sha(E/'SHA256SUMS')==seal['manifest_sha256'];assert sha(E/'receipt.json')==seal['receipt_sha256']
actual={str(p.relative_to(W)) for p in E.rglob('*') if p.is_file() and p not in [E/'SHA256SUMS',E/'SEAL.json']}
actual.update(['docs/developers.toml','docs/architecture/generation-grouped-identical-migration.typ']);assert set(manifest)==actual
for name,key in [('large-evidence-bindings.json','all_own_evidence_bindings'),('input-bindings.json','bindings')]:
 d=json.loads((E/name).read_text())
 for record in d[key]:
  durable=W/record['durable_copy'];assert sha(durable)==record['copy_sha256']
  kind=record['encoding']
  if kind=='gzip-chunks':
   pieces=[]
   for part in record['parts']:
    chunk_path=W/part['path'];assert sha(chunk_path)==part['sha256'];assert chunk_path.stat().st_size==part['bytes'];pieces.append(chunk_path.read_bytes())
   raw=gzip.decompress(b''.join(pieces))
  elif kind.startswith('gzip'):raw=gzip.decompress(durable.read_bytes())
  else:raw=durable.read_bytes()
  assert hashlib.sha256(raw).hexdigest()==record['sha256'] and len(raw)==record['bytes']
  original=Path(record.get('external_path',record.get('original')));assert sha(original)==record['sha256'],original
newfiles=[p for p in E.rglob('*') if p.is_file()]+[W/'docs/architecture/generation-grouped-identical-migration.typ'];assert all(p.stat().st_size<=1048576 for p in newfiles)
limit=subprocess.check_output(['jj','config','get','snapshot.max-new-file-size'],cwd=W,text=True).strip();assert limit=='1MiB'
r=subprocess.run(['git','--git-dir='+str(W.parent/'main/.git'),'--work-tree='+str(W),'diff','--check'],cwd=W,capture_output=True,text=True);assert r.returncode==0,(r.stdout,r.stderr)
args=['jj','--config','snapshot.auto-track=none()','file','track',str(E.relative_to(W)),'docs/architecture/generation-grouped-identical-migration.typ','docs/developers.toml']
r=subprocess.run(args,cwd=W,capture_output=True,text=True);assert r.returncode==0,(r.stdout,r.stderr)
(O/'tracking-record.json').write_text(json.dumps({'argv':args,'exit':r.returncode,'snapshot_auto_track':'none()','tracking_limit':limit,'max_new_file_bytes':max(p.stat().st_size for p in newfiles),'explicit_new_paths':len(newfiles),'manifest_complete':True,'all_original_and_durable_input_bytes_verified':True,'git_diff_check_exit':0},indent=2)+'\n')
print('PASS complete manifest, original and durable lossless bindings, limits and explicit tracking; no native/build rerun',flush=True)
