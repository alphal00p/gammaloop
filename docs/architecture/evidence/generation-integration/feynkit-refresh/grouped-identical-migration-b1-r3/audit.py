"""Exact immutable Git object audit. Never uses working files as manifest evidence."""
from pathlib import Path
import gzip,hashlib,json,os,subprocess,sys
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh');O=Path(__file__).resolve().parent
OLD='7c0e13a08786416a9f24b3bec6df8a50e4592346';CODE='6f41807b064bcf2f579883134daef65c8ac279d3'
E='docs/architecture/evidence/generation-integration/feynkit-refresh/grouped-identical-migration'
CORRECTION='docs/architecture/evidence/generation-integration/feynkit-refresh/grouped-identical-migration-b1-r3'
ENV={k:os.environ[k] for k in ['PATH','HOME','USER','LOGNAME','LANG','LC_ALL'] if k in os.environ};ENV.update(GIT_OPTIONAL_LOCKS='0',GIT_CONFIG_NOSYSTEM='1',GIT_CONFIG_GLOBAL='/dev/null')
def digest(b):return hashlib.sha256(b).hexdigest()
def blob(b):return hashlib.sha1(b'blob %d\0'%len(b)+b).hexdigest()
def git(*args,input=None):
 p=subprocess.run(['git','--git-dir='+str(W.parent/'main/.git'),*args],cwd=W,env=ENV,input=input,capture_output=True);assert p.returncode==0,(args,p.returncode,p.stderr.decode());return p.stdout

def tree(rev):
 result={}
 for row in git('ls-tree','-rz','--full-tree',rev).split(b'\0'):
  if row:
   header,path=row.split(b'\t',1);mode,kind,oid=header.decode().split();assert kind=='blob';result[os.fsdecode(path)]={'mode':mode,'blob':oid}
 return result

def batch(oids):
 order=list(dict.fromkeys(oids));raw=git('cat-file','--batch',input=''.join(oid+'\n' for oid in order).encode());pos=0;out={}
 for expected in order:
  end=raw.index(b'\n',pos);oid,kind,size=raw[pos:end].decode().split();assert oid==expected and kind=='blob';size=int(size);pos=end+1;data=raw[pos:pos+size];pos+=size;assert raw[pos:pos+1]==b'\n';pos+=1;assert blob(data)==oid;out[oid]=data
 assert pos==len(raw);return out

def audit(rev,historical=False):
 initial=json.loads((O/'initial.json').read_text());oldtree=tree(OLD);newtree=tree(rev);missing13={r['path'] for r in initial['missing_paths']}
 def read(path):
  assert path in newtree,('not in immutable tree',rev,path)
  oid=newtree[path]['blob']
  if oid not in cache:cache.update(batch([oid]))
  return cache[oid]
 cache={};manifest=read(E+'/SHA256SUMS');assert digest(manifest)=='f9a8dec1335e5d23f90561c919b6e59136da023f344ee5f120217afd7a1fffdf'
 expected={}
 for i,line in enumerate(manifest.decode().splitlines(),1):
  h,p=line.split('  ',1);assert p not in expected and not Path(p).is_absolute() and '..' not in Path(p).parts;expected[p]=(h,i)
 assert len(expected)==425
 cache.update(batch([newtree[p]['blob'] for p in expected if p in newtree]))
 missing=[];mismatch=[];rows=[]
 for p,(h,line) in expected.items():
  if p not in newtree:missing.append(p);continue
  data=read(p);actual=digest(data);row={'path':p,'manifest_line':line,'mode':newtree[p]['mode'],'git_blob':newtree[p]['blob'],'sha256':actual,'bytes':len(data),'expected_sha256':h,'match':actual==h};rows.append(row)
  if actual!=h:mismatch.append(row)
 scoped={p for p in newtree if p.startswith(E+'/')}|{'docs/architecture/generation-grouped-identical-migration.typ','docs/developers.toml'}
 allowed_root={E+'/SHA256SUMS',E+'/SEAL.json'};unexpected=sorted(scoped-set(expected)-allowed_root)
 assert not mismatch and not unexpected,(mismatch,unexpected)
 if historical:assert set(missing)==missing13 and len(rows)==412
 else:assert not missing and len(rows)==425
 additions=sorted(set(newtree)-set(oldtree));removed=sorted(set(oldtree)-set(newtree));changed=sorted(p for p in set(oldtree)&set(newtree) if oldtree[p]!=newtree[p]);assert not removed and not changed,(removed,changed)
 if not historical:assert set(additions)==missing13|{p for p in newtree if p.startswith(CORRECTION+'/')},additions
 seal=json.loads(read(E+'/SEAL.json'));assert digest(manifest)==seal['manifest_sha256'] and digest(read(E+'/receipt.json'))==seal['receipt_sha256']
 result={'schema':'immutable-git-manifest-audit/1','immutable_commit':rev,'base':OLD,'manifest_sha256':digest(manifest),'manifest_entries':425,'matched':len(rows),'missing_paths':missing,'mismatched_paths':mismatch,'unexpected_scoped_paths':unexpected,'allowed_nonmanifest_packet_roots':sorted(allowed_root),'entries':rows,'all_preexisting_paths_mode_blob_invariant':True,'added_paths':additions,'removed_paths':removed,'changed_preexisting_paths':changed,'source_code':CODE,'new_runtime_or_builds':0,'method':'git ls-tree exact immutable revision; git cat-file --batch every manifest blob; Git object SHA1 and manifest SHA256 on actual blob bytes. Canonical note and registry included. No working-tree/git-ls-files inference.'}
 result['scoped_manifest_tree_sha256']=digest(json.dumps([{k:r[k] for k in ['path','mode','git_blob','sha256','bytes']} for r in sorted(rows,key=lambda r:r['path'])],sort_keys=True,separators=(',',':')).encode())
 if historical:return result
 # All bindings are reconstructed from immutable Git blobs, including chunk indexes.
 logical={};records=[]
 for name,key in [('input-bindings.json','bindings'),('large-evidence-bindings.json','all_own_evidence_bindings')]:
  for r in json.loads(read(E+'/'+name))[key]:
   path=r['durable_copy'];encoded=read(path);assert digest(encoded)==r['copy_sha256'];assert len(encoded)==r['copy_bytes'];kind=r['encoding']
   if kind=='gzip-chunks':
    index=json.loads(encoded);assert index['parts']==r['parts'];pieces=[]
    for part in r['parts']:
     data=read(part['path']);assert digest(data)==part['sha256'] and len(data)==part['bytes'];pieces.append(data)
    raw=gzip.decompress(b''.join(pieces));assert path.endswith('.gz.index.json');target=path[:-len('.gz.index.json')]
   elif kind.startswith('gzip'):
    raw=gzip.decompress(encoded);assert path.endswith('.gz');target=path[:-3]
   else:assert kind=='identity';raw=encoded;target=path
   assert digest(raw)==r['sha256'] and len(raw)==r['bytes'];assert target not in logical or logical[target]==raw;logical[target]=raw
   records.append({'binding_file':name,'durable_copy':path,'logical_path':target,'encoding':kind,'sha256':digest(raw),'bytes':len(raw),'committed_blob_reconstruction':True})
 def copied(path):return logical[path] if path in logical else read(path)
 migration_root=E+'/reviews/migration';spec=json.loads(copied(migration_root+'/seal.json'));migration=[]
 for name,r in spec['files'].items():
  data=copied(migration_root+'/'+name);assert digest(data)==r['sha256'] and len(data)==r['bytes'];migration.append({'path':migration_root+'/'+name,'sha256':digest(data),'bytes':len(data)})
 assert digest(copied(migration_root+'/review.json'))=='7d96505dcb700eeae42b5a4e87676551ba910cb660df586cafe34d6b90e944f1'
 prep=E+'/preparation';spec=json.loads(copied(prep+'/seal.json'));sums=copied(prep+'/SHA256SUMS');assert digest(sums)==spec['SHA256SUMS_sha256'];assert digest(copied(prep+'/manifest.json'))==spec['manifest_sha256'];preparation=[]
 for line in sums.decode().splitlines():
  h,p=line.split('  ',1);data=copied(prep+'/'+p);assert digest(data)==h;preparation.append({'path':prep+'/'+p,'sha256':h,'bytes':len(data)})
 result.update(lossless_bindings=records,lossless_bindings_count=len(records),complete_copied_migration_review_seal=migration,complete_copied_preparation_manifest=preparation,reviewed_runtime_reused={'review':'0f0851e068653202a7981e71b4d5ccf54316d1bf544caa0725f99cc1794066f4','review_passes':39,'public_captures':5,'owner_associations':2124,'fresh_executions_in_repair':0},no_selfcontained_external_artifact_claim=True)
 return result

if __name__=='__main__':
 rev=sys.argv[1];historical='--historical' in sys.argv;result=audit(rev,historical);path=O/('historical-7c-audit.json' if historical else 'immutable-'+rev+'-audit.json');path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['immutable_commit','manifest_entries','matched','missing_paths','mismatched_paths','unexpected_scoped_paths','scoped_manifest_tree_sha256']}))
