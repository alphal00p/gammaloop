from pathlib import Path
import json,math,hashlib,tarfile
out=Path('target/mark-geometry-round3');old=json.loads((out/'old-public.json').read_text())[0];new=json.loads((out/'new-public.json').read_text())[0];rows=[]
for a,b in zip(old,new):
 assert a['name']==b['name'];bounds=lambda r:{'left':min(x['left'] for x in r['outline-bounds']),'right':max(x['right'] for x in r['outline-bounds']),'top':max(x['top'] for x in r['outline-bounds'])};aa=bounds(a);bb=bounds(b)
 if a['name'].startswith('triangle'):assert abs(bb['right']-100)<1e-6
 else:
  i=.99999999 if a['name']=='stealth-deeper' else .999999;h=4.;l=10.;radius=h*l*(1-i)/(math.hypot(l,h)+math.hypot(l*i,h));grow=1.01*.5/radius;assert abs(bb['left']-(100-grow*l))<1e-4*grow*l;assert abs(bb['top']-grow*h)<1e-4*grow*h
 rows.append({'name':a['name'],'old':aa,'corrected':bb,'old_wasm_sha256':'31f7fbb7da4520f7c200e5ae7ead5a462c20e1c1bea9c6c020458a6bad04e27d','corrected_wasm_sha256':hashlib.sha256(Path('crates/kurvst/typst/kurvst.wasm').read_bytes()).hexdigest()})
assert rows[0]['old']['right']<-1e7 and rows[1]['old']['left']>rows[1]['corrected']['left']/10
(out/'public-comparison.json').write_text(json.dumps(rows,indent=2)+'\n')
# Preserve an actual negative source-delivery control from before explicit include.
oldarchive=next((out/'sdist').glob('*.tar.gz'));newarchive=next((out/'sdist-final').glob('*.tar.gz'))
with tarfile.open(oldarchive) as t:assert not any('/vendor/kurbo/' in n for n in t.getnames())
with tarfile.open(newarchive) as t:assert len([n for n in t.getnames() if '/vendor/kurbo/' in n])==53
(out/'sdist-negative-control.json').write_text(json.dumps({'old_archive_sha256':hashlib.sha256(oldarchive.read_bytes()).hexdigest(),'old_kurbo_files':0,'new_archive_sha256':hashlib.sha256(newarchive.read_bytes()).hexdigest(),'new_kurbo_files':53,'scope':'Real source archives; no final wheels claimed.'},indent=2)+'\n')
print('Two old public failures and two controls now pass; real sdist omission negative control rejects old archive')
