from pathlib import Path
import subprocess,json,tomllib,hashlib,tarfile,importlib.util
r=Path.cwd(); out=r/'target/wp3-source38'; base='100da20603a083d031c31bd4e8c6970fc44dda2f'
def git(*a):return subprocess.check_output(['git','--git-dir=../main/.git',*a])
def h(b):return hashlib.sha256(b).hexdigest()
report={}
for p in ['Cargo.lock','examples/notebooks/symbolica-host/Cargo.lock']:
 old=tomllib.loads(git('show',base+':'+p).decode()); new=tomllib.loads((r/p).read_text()); a=old['package'];b=new['package'];assert len(a)==len(b)
 delta=[]
 for x,y in zip(a,b):
  if x!=y:
   assert x['name']==y['name']=='kurbo' and x['version']==y['version']=='0.13.1'
   assert {k:v for k,v in x.items() if k not in ['source','checksum']}==y
   delta.append(x['name'])
 assert delta==['kurbo'];report[p]=delta
for p in ['Tydenso/Cargo.lock','uv.lock','flake.lock']:
 if (r/p).exists():assert git('show',base+':'+p)==(r/p).read_bytes()
archive=r.parent/'generation-drawing-ports/.local/mark-geometry-round2-producers/kurbo-0.13.1.crate';assert h(archive.read_bytes())=='4b60dfc32f652b926df6192e55525b16d186c69d47876c3ead4da5cc9f8450e2'
with tarfile.open(archive) as t:
 files=[x for x in t.getmembers() if x.isfile()];assert len(files)==53
 delta=[]
 for x in files:
  p='/'.join(x.name.split('/')[1:]); b=t.extractfile(x).read(); current=(r/'vendor/kurbo'/p).read_bytes()
  if b!=current:delta.append(p)
 assert delta==['src/stroke.rs'];report['published_vendor']={'files':53,'only_delta':delta,'archive_sha256':h(archive.read_bytes())}
spec=importlib.util.spec_from_file_location('canonical',r/'docs/architecture/evidence/generation-integration/wp3/oracle-registry-audit.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
reg=tomllib.loads((r/'docs/developers.toml').read_text());rows=[];titles=[]
for n,s in m.scopes(reg):
 result=m.check(s,(r/s['path']).read_text());assert result['status']=='valid',(n['id'],s,result);rows.append({'note':n['id'],**s,**result})
for sec in reg['section']:
 for n in sec.get('note',[]):
  p=r/n['source'];text=p.read_text();title=next(x[2:].strip() for x in text.splitlines() if x.startswith('= '));assert title==n['title'],(n,title);titles.append({'id':n['id'],'title':title})
report['registry']={'scopes':len(rows),'titles':len(titles),'rows':rows}
(out/'audit.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS lock/vendor/scopes/titles',len(rows),len(titles))
