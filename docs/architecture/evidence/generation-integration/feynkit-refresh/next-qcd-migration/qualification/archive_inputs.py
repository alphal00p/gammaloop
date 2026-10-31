"""Commit reviewable input bytes; large producer dependencies remain explicit external bindings."""
import gzip,hashlib,json
from pathlib import Path
O=Path(__file__).resolve().parent;W=O.parents[2];F=W.parent/'generation-oracle-numerator-families'
E=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/next-qcd-migration'
E.mkdir(parents=True,exist_ok=False);bindings=[]
def sha(b):return hashlib.sha256(b).hexdigest()
def cp(source,relative):
 source=Path(source);b=source.read_bytes();dest=E/relative;dest.parent.mkdir(parents=True,exist_ok=True);r={'original':str(source),'bytes':len(b),'sha256':sha(b)}
 if len(b)<=1048576:
  dest.write_bytes(b);r.update(encoding='identity',copies=[{'path':str(dest.relative_to(W)),'bytes':len(b),'sha256':sha(b)}])
 else:
  encoded=gzip.compress(b,compresslevel=9,mtime=0);need=1048576;parts=[];reassembled=[]
  for i in range(0,len(encoded),need):
   part=encoded[i:i+need];suffix='.gz' if len(encoded)<=need else f'.gz.part{i//need:03d}';p=Path(str(dest)+suffix);p.write_bytes(part);parts.append({'path':str(p.relative_to(W)),'bytes':len(part),'sha256':sha(part)});reassembled.append(part)
  assert gzip.decompress(b''.join(reassembled))==b
  r.update(encoding='gzip' if len(parts)==1 else 'concatenate-ordered-gzip-parts',compressed_bytes=len(encoded),compressed_sha256=sha(encoded),copies=parts)
 bindings.append(r)
def packet(root,manifest,dest):
 for line in (root/manifest).read_text().splitlines():
  h,name=line.split('  ',1);p=Path(name);assert not p.is_absolute() and '..' not in p.parts;assert sha((root/p).read_bytes())==h;cp(root/p,Path(dest)/p)
 cp(root/manifest,Path(dest)/manifest)
 for seal in ['SEAL.json','SEAL.sha256','COMPLETE-SEAL.json']:
  if (root/seal).exists() and not any(r['original']==str(root/seal) for r in bindings):cp(root/seal,Path(dest)/seal)
packet(F/'.local/review-next-qcd-migration-astra-20261031-efvezb1p','SHA256SUMS','reviews/migration')
packet(F/'.local/review-next-qcd-scientific-interface-astra-20261031-5a3KbP5x','SHA256SUMS','reviews/scientific-interface')
packet(F/'target/raw-next-qcd-migration-preparation-20261031-r1','SHA256SUMS','preparation')
packet(F/'target/raw-next-qcd-proof-native-execution-20261031-r1','COMPLETE-SHA256SUMS','qualified-proof')
N=F/'target/raw-next-qcd-proof-native-execution-20261031-r1';R3=F/'target/raw-next-qcd-proof-preparation-20261031-r3/execution'
for line in (N/'EXECUTION-SHA256SUMS').read_text().splitlines():
 h,p=line.split('  ',1);p=Path(p);assert p.is_relative_to(R3) and sha(p.read_bytes())==h;cp(p,Path('qualified-proof/actual-execution')/p.relative_to(R3))
for name in ['root-next-qcd-migration-acceptance-and-application-authorization-20261031.json','root-next-qcd-scientific-interface-acceptance-20261031.json','root-next-qcd-migration-preparation-receipt-20261031.json','root-next-qcd-proof-native-execution-assignment-20261031-r3.json']:
 cp(W/'.local/feynkit-refresh'/name,Path('roots')/name)
CAP=F/'target/raw-next-qcd-captures-20261031-r2'
for p in (CAP/'prepared/captured').rglob('*.json'):cp(p,Path('captures/retained-worlds')/p.relative_to(CAP/'prepared/captured'))
for name in ['SEAL.json','SHA256SUMS','handoff.json','handoff.typ']:cp(CAP/name,Path('captures/provenance')/name)
R=F/'.local/review-next-qcd-captures-astra-20261031-IYGt5Jrv'
for p in sorted(R.rglob('*')):
 if p.is_file():cp(p,Path('reviews/capture-provenance')/p.relative_to(R))
(O/'input-copy-bindings.json').write_text(json.dumps(bindings,indent=2)+'\n')
print('Preserved',len(bindings),'exact input records in',E)
