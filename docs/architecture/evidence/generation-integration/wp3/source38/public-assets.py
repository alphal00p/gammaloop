from pathlib import Path
import json,subprocess,hashlib,shutil
T=Path('target/wp3-source38');R=Path('../generation-drawing-ports');rows=[]
def run(name,cmd):
 r=subprocess.run(cmd,capture_output=True);(T/(name+'.log')).write_bytes(r.stdout+r.stderr);assert r.returncode==0,(name,r.stderr.decode()[-2000:]);rows.append(dict(name=name,command=cmd,exit=r.returncode));return r.stdout
for name,src,label in [('adversarial',R/'.local/review-mark-geometry-round2-a-astra/adversarial-public.typ','boundary-results'),('original',R/'.local/review-mark-geometry-a/boundaries.typ','boundary-results'),('generic',R/'.local/review-mark-geometry-round3-a-astra-20261029/public-generic.typ','generic')]:
 dst=T/(name+'.typ');dst.write_bytes(src.read_bytes());out=run('public-'+name,['typst','query','--root','.',str(dst),'<'+label+'>','--field','value']);value=json.loads(out);(T/('public-'+name+'.json')).write_bytes(out)
run('public-fixture',['typst','compile','--creation-timestamp','0','--root','.','crates/kurvst/typst/tests/mark-geometry.typ',str(T/'public-fixture.pdf')])
for name,src in [('manual','crates/kurvst/typst/docs/manual.typ'),('soft-triangle','crates/linnest/typst/examples/soft-triangle.typ'),('thrown-out','crates/linnest/typst/examples/thrown-out.typ'),('xbox-sum','crates/linnest/typst/examples/xbox-sum.typ')]:
 dst=T/(name+'.pdf');run(name,['typst','compile','--features','html','--root','.','--creation-timestamp','0',src,str(dst)]);assert dst.read_bytes()==Path(src).with_suffix('.pdf').read_bytes(),name
for gallery in ['marks','boundaries']:
 for theme in ['light','dark']:
  for fmt in ['svg','png','pdf']:
   dst=T/(gallery+'-'+theme+'.'+fmt);run(dst.stem+'-'+fmt,['typst','compile','--root','.','--creation-timestamp','0','--input','theme='+theme,'--format',fmt,'crates/kurvst/typst/examples/'+('marks' if gallery=='marks' else 'mark-boundaries')+'.typ',str(dst)]);expected=Path('docs/architecture/evidence/generation-integration/mark-geometry/round2-correction')/dst.name;assert dst.read_bytes()==expected.read_bytes(),dst
run('website-assets',['bash','scripts/render-docs-svg-assets.sh',str(T/'reproduced')]);assets=list((T/'reproduced').rglob('*.svg'));assert len(assets)==34
for p in assets:assert p.read_bytes()==Path(p.relative_to(T/'reproduced')).read_bytes(),p
(T/'public-assets.json').write_text(json.dumps(dict(commands=rows,website_SVGs=34,manual_and_examples=4,gallery_outputs=12,all_reproduced=True,qualification='Actual current tracked WASMs and native public templates; numeric outline controls tested, but blank skinny triangle raster crop is retained downstream limitation, NOT visual coverage proof. No extreme-renderer guarantee.'),indent=2)+'\n');print('Public adversarial/original/generic controls;34SVG+manual/3PDF+12gallery outputs reproduce')
