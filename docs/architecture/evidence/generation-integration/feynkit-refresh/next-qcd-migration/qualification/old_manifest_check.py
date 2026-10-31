"""Revision-aware immutable old manifest validation; no current-file fallback."""
import hashlib,json,os,subprocess
from pathlib import Path
O=Path(__file__).resolve().parent;W=O.parents[2];REV='7f03ab8fd5834d7c78164429e2a44dab8d01ef74'
E='docs/architecture/evidence/generation-integration/feynkit-refresh/grouped-identical-migration'
env={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')};env['GIT_OPTIONAL_LOCKS']='0'
def git(*a,data=None):return subprocess.check_output(['git','--git-dir='+str(W.parent/'main/.git'),*a],input=data,env=env,cwd=W)
def sha(b):return hashlib.sha256(b).hexdigest()
manifest=git('show',REV+':'+E+'/SHA256SUMS');assert sha(manifest)=='f9a8dec1335e5d23f90561c919b6e59136da023f344ee5f120217afd7a1fffdf'
entries=[]
for line in manifest.decode().splitlines():
 h,p=line.split('  ',1);entries.append((h,p))
assert len(entries)==425
raw=git('cat-file','--batch',data=''.join(REV+':'+p+'\n' for h,p in entries).encode());pos=0;records=[]
for h,p in entries:
 end=raw.index(b'\n',pos);oid,kind,size=raw[pos:end].decode().split();assert kind=='blob';size=int(size);pos=end+1;b=raw[pos:pos+size];pos+=size;assert raw[pos:pos+1]==b'\n';pos+=1
 assert sha(b)==h and hashlib.sha1(b'blob %d\0'%len(b)+b).hexdigest()==oid;records.append({'path':p,'git_blob':oid,'sha256':h,'bytes':len(b)})
assert pos==len(raw)
reg=git('show','83276dc047815d1dcc7f003aeb512c204b26758e:docs/developers.toml');old=git('show',REV+':docs/developers.toml')
result={'schema':'revision-aware-prior-manifest-preservation/1','immutable_source_revision':REV,'manifest_sha256':sha(manifest),'matched':len(records),'missing':[],'mismatched':[],'entries':records,'no_working_tree_fallback':True,'old_false_closure':'7c0e13a0 remains false412/425; complete corrected7f and all prior refs preserved','registry_old_manifest_source':{'revision':REV,'sha256':sha(old),'blob':hashlib.sha1(b'blob %d\0'%len(old)+old).hexdigest()},'registry_before_new_note':{'revision':'83276dc047815d1dcc7f003aeb512c204b26758e','sha256':sha(reg),'blob':hashlib.sha1(b'blob %d\0'%len(reg)+reg).hexdigest()},'registry_history':'832 contains already accepted status attribution. New note will be additive. The old manifest is verified at its exact7f revision, never against rewritten current metadata.'}
(O/'prior-manifest-revision-aware.json').write_text(json.dumps(result,indent=2)+'\n');(O/'registry-before-832.toml').write_bytes(reg)
print('PASS all425 old manifest bytes from immutable7f, with exact832 pre-edit registry snapshot')
