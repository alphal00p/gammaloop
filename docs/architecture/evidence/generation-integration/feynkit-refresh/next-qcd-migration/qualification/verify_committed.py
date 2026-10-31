"""Verify exact immutable Git blobs and reconstruct every lossless evidence binding."""
import gzip,hashlib,json,os,subprocess,sys
from pathlib import Path
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh');O=Path(__file__).resolve().parent
E='docs/architecture/evidence/generation-integration/feynkit-refresh/next-qcd-migration';D=E+'-durability'
ENV={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')};ENV.update(GIT_OPTIONAL_LOCKS='0',GIT_CONFIG_NOSYSTEM='1',GIT_CONFIG_GLOBAL='/dev/null')
def sha(b):return hashlib.sha256(b).hexdigest()
def blob(b):return hashlib.sha1(b'blob %d\0'%len(b)+b).hexdigest()
def git(*a,data=None):return subprocess.check_output(['git','--git-dir='+str(W.parent/'main/.git'),*a],cwd=W,env=ENV,input=data)
def tree(rev):
 out={}
 for r in git('ls-tree','-rz','--full-tree',rev).split(b'\0'):
  if r:
   h,p=r.split(b'\t',1);mode,kind,oid=h.decode().split();assert kind=='blob';out[os.fsdecode(p)]={'mode':mode,'blob':oid}
 return out

def audit(rev):
 t=tree(rev);cache={}
 def fetch(oids):
  order=list(dict.fromkeys(oids));raw=git('cat-file','--batch',data=''.join(o+'\n' for o in order).encode());pos=0
  for o in order:
   end=raw.index(b'\n',pos);oid,kind,n=raw[pos:end].decode().split();assert oid==o and kind=='blob';n=int(n);pos=end+1;b=raw[pos:pos+n];pos+=n;assert raw[pos:pos+1]==b'\n';pos+=1;assert blob(b)==oid;cache[oid]=b
  assert pos==len(raw)
 def read(p):
  assert p in t,('missing immutable path',p)
  o=t[p]['blob']
  if o not in cache:fetch([o])
  return cache[o]
 manifests=[E]
 if D+'/SHA256SUMS' in t:manifests.append(D)
 reports=[];allmembers={}
 for root in manifests:
  manifest=read(root+'/SHA256SUMS');seal=json.loads(read(root+'/SEAL.json'));assert sha(manifest)==seal['manifest_sha256']
  members={}
  for line in manifest.decode().splitlines():
   h,p=line.split('  ',1);assert p not in members and not Path(p).is_absolute() and '..' not in Path(p).parts;members[p]=h
  assert len(members)==seal['files']
  missing=sorted(set(members)-set(t));assert not missing,missing
  fetch([t[p]['blob'] for p in members]);rows=[];bad=[]
  for p,h in sorted(members.items()):
   b=read(p);r={'path':p,'mode':t[p]['mode'],'git_blob':t[p]['blob'],'sha256':sha(b),'bytes':len(b)};rows.append(r)
   if r['sha256']!=h:bad.append(r)
  scoped={p for p in t if p.startswith(root+'/')};unexpected=sorted(scoped-set(members)-{root+'/SHA256SUMS',root+'/SEAL.json'});assert not bad and not unexpected,(bad,unexpected)
  reports.append({'root':root,'manifest_sha256':sha(manifest),'seal_sha256':sha(read(root+'/SEAL.json')),'manifest_entries':len(members),'matched':len(rows),'missing':missing,'mismatched':bad,'unexpected':unexpected,'entries':rows,'vector_sha256':sha(json.dumps(rows,sort_keys=True,separators=(',',':')).encode())});allmembers.update(members)
 recon=[]
 for name in ['input-copy-bindings.json','qualification-copy-bindings.json']:
  records=json.loads(read(E+'/'+name))
  for r in records:
   parts=[]
   for p in r['copies']:
    b=read(p['path']);assert sha(b)==p['sha256'] and len(b)==p['bytes'];parts.append(b)
   b=b''.join(parts)
   if r['encoding']!='identity':
    assert sha(b)==r['compressed_sha256'];b=gzip.decompress(b)
   assert sha(b)==r['sha256'] and len(b)==r['bytes'];recon.append({'original':r['original'],'sha256':sha(b),'bytes':len(b),'copies':[p['path'] for p in r['copies']]})
 receipt=json.loads(read(E+'/receipt.json'));assert receipt['status']=='NEW_INDEPENDENT_INTEGRATION_REVIEW_REQUIRED';assert receipt['runtime']['actual_executions']==receipt['runtime']['passes']==48 and receipt['runtime']['failures']==0
 base=tree('83276dc047815d1dcc7f003aeb512c204b26758e');removed=sorted(set(base)-set(t));assert not removed
 changed=sorted(p for p in base.keys()&t.keys() if base[p]!=t[p]);assert changed==['crates/feynkit-generator/tests/generation_oracle/main.rs','docs/developers.toml'],changed
 added=sorted(set(t)-set(base));assert all(p.startswith(E+'/') or p.startswith(D+'/') or p=='docs/architecture/generation-next-qcd-migration.typ' for p in added)
 assert blob(read('crates/feynkit-generator/tests/generation_oracle/main.rs'))=='cc4505a6c4f6fc1687c075ed79e154a8e19f05a6'
 return {'schema':'next-qcd-immutable-Git-manifest-closure/1','immutable_commit':rev,'method':'git ls-tree and cat-file --batch exact immutable blobs; every Git object SHA1 and manifest SHA256; every lossless copy reconstructed from committed bytes','packets':reports,'lossless_reconstructed_records':len(recon),'reconstructions':recon,'missing':[],'mismatched':[],'unexpected':[],'base':'83276dc047815d1dcc7f003aeb512c204b26758e','changed_prior_paths':changed,'added_paths':added,'self_contained':False,'external_dependencies':'Original source/store/capture/proof paths remain explicit retained premises; committed closure does not claim universal provenance or independent integration acceptance'}
if __name__=='__main__':
 rev=sys.argv[1];r=audit(rev);output=Path(sys.argv[2]);output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'commit':rev,'packets':[{k:p[k] for k in ['root','manifest_entries','matched','missing','mismatched','unexpected']} for p in r['packets']],'reconstructed':r['lossless_reconstructed_records']}))
