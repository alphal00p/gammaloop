"""Retain final bounded-check diagnostics and reseal exact durable evidence."""
from pathlib import Path
import gzip,hashlib,json
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh');O=Path(__file__).resolve().parent;E=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/grouped-identical-migration'
def load(p):return json.loads(Path(p).read_text())
def save(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
f=load(O/'setup-failures.json');f['retained'].append({'attempt':'verify-packet-attempt2.py','stderr':'verify-packet-attempt2.stderr','cause':'Broad Git diff --check used shared main index, exposing inherited whitespace outside this task. Corrected to immutable457->actualcode and actualcode->evidence authored canonical note/registry. No unrelated repair or native/build rerun; bounded checks pass.'});save(O/'setup-failures.json',f)
receipt=load(E/'receipt.json');receipt['retained_failures']=f;receipt['validation']['bounded_diff_checks']={'code':'457e361d7d39e0a3daae801651f8c1dfabbf47aa ->6f41807b064bcf2f579883134daef65c8ac279d3','authored_docs':'actual code -> evidence working tree, exact canonical note and registry','exit':0,'frozen_verbatim_evidence_not_reformatted':True,'broad_shared_index_attempt_retained':True};receipt['retention']['preflight_bookmark']='generation-grouped-identical-migration-evidence-preflight-20261031';save(E/'receipt.json',receipt)
bind=load(E/'large-evidence-bindings.json');rows=bind['all_own_evidence_bindings']
for name in ['setup-failures.json','verify_packet.py','verify-packet-attempt2.py','verify-packet-attempt2.stderr','verify-packet.stdout','verify-packet.stderr','tracking-record.json','postseal.py']:
 p=O/name;b=p.read_bytes();assert len(b)<=1048576;dest=E/'qualification'/name;dest.write_bytes(b)
 record={'external_path':str(p),'sha256':sha(p),'bytes':len(b),'durable_copy':str(dest.relative_to(W)),'copy_sha256':sha(dest),'copy_bytes':len(b),'encoding':'identity','parts':[],'durable_reference':'generation-grouped-identical-migration-evidence-r2-20261031'}
 rows[:]=[r for r in rows if r['external_path']!=str(p)];rows.append(record)
save(E/'large-evidence-bindings.json',bind)
for r in rows:
 p=W/r['durable_copy'];assert sha(p)==r['copy_sha256']
 if r['encoding']=='gzip-chunks':
  chunks=[]
  for part in r['parts']:
   path=W/part['path'];assert sha(path)==part['sha256'];chunks.append(path.read_bytes())
  raw=gzip.decompress(b''.join(chunks))
 elif r['encoding']=='gzip':raw=gzip.decompress(p.read_bytes())
 else:raw=p.read_bytes()
 assert hashlib.sha256(raw).hexdigest()==r['sha256'] and len(raw)==r['bytes'];assert sha(r['external_path'])==r['sha256']
files=sorted([p for p in E.rglob('*') if p.is_file() and p not in [E/'SHA256SUMS',E/'SEAL.json']]+[W/'docs/developers.toml',W/'docs/architecture/generation-grouped-identical-migration.typ']);assert all(p.stat().st_size<=1048576 for p in files)
(E/'SHA256SUMS').write_text(''.join(sha(p)+'  '+str(p.relative_to(W))+'\n' for p in files))
seal=load(E/'SEAL.json');seal.update(manifest_sha256=sha(E/'SHA256SUMS'),receipt_sha256=sha(E/'receipt.json'),files=len(files),bytes=sum(p.stat().st_size for p in files));save(E/'SEAL.json',seal)
for line in (E/'SHA256SUMS').read_text().splitlines():h,p=line.split('  ',1);assert sha(W/p)==h
for p in E.rglob('*.json'):json.loads(p.read_text())
print('PASS final complete packet',len(files),'manifest files;seal',sha(E/'SEAL.json'))
