import sys,re,tomllib,collections
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import common as C

def tree(rev):
 return {p:tuple(meta.split()) for meta,p in (line.split('\t',1) for line in C.git('ls-tree','-r',rev).decode().splitlines())}

def main():
 assert C.git('rev-parse',C.REV+'^{tree}').decode().strip()=='c9de5c3fd962bd5d7c9375ab1f86f75699e8a69b'
 trees={r:tree(r) for r in (C.ACT,C.REV,C.CODE)}
 a,b,code=(trees[r] for r in (C.ACT,C.REV,C.CODE))
 domain=lambda t:{p:v for p,v in t.items() if p.startswith(('crates/','assets/')) or p in ('Cargo.toml','Cargo.lock')}
 da,db,dc=map(domain,(a,b,code));assert da.keys()==db.keys()==dc.keys()
 assert len(da)==1521
 assert [p for p in da if da[p]!=db[p]]==[C.MAIN]
 assert db==dc
 # Read every exact immutable source object; digest bytes independently, including config source objects.
 objs={};records=[]
 for p in db:
  for rev,t in ((C.ACT,da),(C.REV,db)):
   mode,typ,obj=t[p];assert typ=='blob'
   if obj not in objs:
    data=C.git('cat-file','blob',obj);assert C.blob(data)==obj
    objs[obj]={'sha256':C.sha(data),'bytes':len(data)}
   records.append({'revision':rev,'path':p,'mode':mode,'blob':obj,**objs[obj]})
 frozen=C.load(C.P/'frozen-source-bindings.json')
 for r in frozen['files']:
  data=C.git('show',r['revision']+':'+r['path']);assert data==C.read(C.P/r['packet_path'])
  assert (C.blob(data),C.sha(data),len(data))==(r['git_blob'],r['sha256'],r['bytes'])
  C.put(r['packet_path'],data)
 before=C.git('show',C.REV+':'+C.MAIN);assert before==C.read(C.P/'before.rs')
 after=C.read(C.P/'after.NOT-APPLIED.rs');actual=C.git('show',C.ACT+':'+C.MAIN)
 assert before==C.git('show',C.CODE+':'+C.MAIN)
 new_group=b'numerators: 0xa60df3f531d57ce3';old_group=b'numerators: 0x00fac82174b46b59'
 # Discover complete unrelated Row, not hard-coded spelling guessed from decimal.
 row=lambda x:next(line for line in x.splitlines(keepends=True) if b' xs_qed_eemm_identical_L2: Row ' in line)
 ar,br=row(actual),row(before)
 mask=lambda x:re.sub(rb'numerators: 0x[0-9a-f]{16}',b'numerators: TOKEN',x)
 assert mask(ar)==mask(br) and ar!=br and before.replace(br,ar,1)==actual
 C.save('source-objects-independent.json',{'source_objects_per_revision':len(da),'unique_objects_read':len(objs),'bindings':records,'actual_producer':C.ACT,'target':C.REV,'relevant_code':C.CODE,'single_unrelated_delta':{'case':'xs_qed_eemm_identical_L2','actual_row':ar.decode(),'target_row':br.decode()},'target_vs_relevant_code_domain_equal':True})
 full=[]
 for p in sorted(a.keys()|b.keys()):
  if a.get(p)!=b.get(p): full.append({'path':p,'actual':a.get(p),'target':b.get(p)})
 C.save('all-target-deltas-independent.json',full)
 print('Full-delta categories',dict(collections.Counter(p['path'].split('/')[0] for p in full)))
 print('All non-doc/non-main deltas',[p['path'] for p in full if not p['path'].startswith('docs/') and p['path']!=C.MAIN])
 registry=C.git('show',C.REV+':docs/developers.toml');parsed=tomllib.loads(registry.decode());bindings=[]
 for section in parsed['section']:
  for note in section.get('note',[]):
   for scope in note.get('scope',[]):
    original=C.git('show',C.REV+':'+scope['path']);candidate=after if scope['path']==C.MAIN else original
    def extract(data):
     text=data.decode();assert text.count(scope['anchor'])==1
     if 'end_anchor' in scope:
      assert text.count(scope['end_anchor'])==1
      lo=text.index(scope['anchor']);hi=text.index(scope['end_anchor']);assert hi>=lo+len(scope['anchor']);return text[lo:hi].encode()
     return data
    old,new=extract(original),extract(candidate);assert old==new and C.blob(new)==scope['digest']
    bindings.append({'position':len(bindings),'note':note['id'],**scope,'actual_digest':C.blob(old),'sha256':C.sha(old),'bytes':len(old)})
 declared=C.load(C.P/'ordered118-source-bindings.json');assert len(bindings)==118 and bindings==declared['ordered_scopes']
 C.save('ordered-scopes-independent.json',{'registry_blob':C.blob(registry),'registry_sha256':C.sha(registry),'count':118,'before_after_equal':True,'ordered_scopes':bindings})
 C.flush('source-consumed.json')
 print('PASS:1521 complete source objects,16 frozen source files,118 ordered scopes; actualffa preserved as actual producer.')
if __name__=='__main__':main()
