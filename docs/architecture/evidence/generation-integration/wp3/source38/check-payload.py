from pathlib import Path
import json,tarfile,hashlib,importlib.util,copy,subprocess
r=Path.cwd();t=r/'target/wp3-source38';report={}
archive=next((t/'sdist').glob('*.tar.gz'))
with tarfile.open(archive) as a:
 names=a.getnames();prefix=names[0].split('/')[0]+'/'
 vendor=subprocess.check_output(['git','--git-dir=../main/.git','ls-tree','-r','--name-only','925f739799572d23f46a762c895ebb85e5bd88d6','vendor/kurbo'],text=True).splitlines();assert len(vendor)==53
 for p in vendor:assert a.extractfile(prefix+p).read()==(r/p).read_bytes(),p
 assert not any('tiptoe' in n.lower() for n in names)
report['sdist']={'path':str(archive.relative_to(r)),'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'vendor_files':len(vendor),'licenses':['vendor/kurbo/LICENSE-MIT','vendor/kurbo/LICENSE-APACHE'],'qualification':'Actual current source archive; not installation/wheel/Pyodide/RECORD/release certification.'}
spec=importlib.util.spec_from_file_location('c',r/'docs/architecture/evidence/generation-integration/wp3/oracle-registry-audit.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
s={'anchor':'START','digest':'bad'};assert m.check(s,'none')['status']=='missing-or-ambiguous-anchor';assert m.check(s,'START START')['status']=='missing-or-ambiguous-anchor';assert m.check(dict(s,end_anchor='END'),'END START')['status']=='invalid-range';assert m.check(dict(s,end_anchor='END'),'START END END')['status']=='missing-or-ambiguous-end';assert m.check(s,'START')['status']=='stale-digest';report['canonical_negative_controls']=5
j=json.loads((t/'public-adversarial.json').read_text())[0];tri=next(x for x in j if x['name']=='triangle-thinner');stealth=next(x for x in j if x['name']=='stealth-deeper');assert abs(tri['outline-bounds'][0]['right']-100)<1e-5;assert stealth['outline-bounds'][0]['left'] < -2.7e8
report['fresh_public_adversarial']={'skinny_triangle_right':tri['outline-bounds'][0]['right'],'deep_stealth_left':stealth['outline-bounds'][0]['left'],'qualification':'Geometric outline bounds; blank extreme skinny raster crop remains downstream limitation, not visual coverage.'}
(t/'payload.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
