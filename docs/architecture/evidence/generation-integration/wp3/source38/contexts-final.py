from pathlib import Path
import json,subprocess
out=Path('target/wp3-source38'); uri='git+file:///shared/lcnbr/dev/gammaloop/main?rev=0c65bf7c519764460a15a922340d958a9dfa7429&shallow=1&allRefs=1'
def cmd(a):return subprocess.check_output(a,text=True).strip()
p=cmd(['nix','eval','--raw',uri+'#packages.x86_64-linux.crate-gammaloop-api.src.outPath'])
d=cmd(['nix','eval','--raw',uri+'#packages.x86_64-linux.crate-test-binaries-gammaloop-api.drvPath'])
def drv(d):
 j=json.loads(cmd(['nix','derivation','show',d]));return next(iter(j.get('derivations',j).values()))
w=drv(d); candidates=[('/nix/store/'+k if not k.startswith('/') else k) for k in w.get('inputs',{}).get('drvs',w.get('inputDrvs',{})) if 'crate-test-binaries-gammaloop-api-' in k]
if 'src' in w.get('env',{}):t=w['env']['src']; actual=d
else:
 assert len(candidates)==1,candidates
 actual=candidates[0];t=drv(actual)['env']['src']
report={'git_input':uri,'production':p,'test_wrapper':d,'test_producer':actual,'test_source':t}
for key in ['production','test_source']:
 src=Path(report[key]); files=[x for x in (src/'vendor/kurbo').rglob('*') if x.is_file()]; assert len(files)==53,(key,len(files))
 for x in files:assert x.read_bytes()==Path(x.relative_to(src)).read_bytes()
 report[key+'_vendor_count']=len(files)
(out/'contexts-final.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
