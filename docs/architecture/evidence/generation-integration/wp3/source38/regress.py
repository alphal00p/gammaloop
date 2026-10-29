from pathlib import Path
import json,hashlib,subprocess
T=Path('target/wp3-source38');R=Path('../generation-drawing-ports');r=subprocess.run(['cargo','build','--locked','-p','kurvst','--profile','dev-optim','--message-format=json'],capture_output=True);(T/'regress-build.stdout').write_bytes(r.stdout);(T/'regress-build.stderr').write_bytes(r.stderr);assert r.returncode==0
for line in r.stdout.decode().splitlines():
 try:d=json.loads(line)
 except json.JSONDecodeError:continue
 if d.get('reason')=='compiler-artifact' and d['target']['name']=='kurvst':
  rlib=next(Path(p) for p in d['filenames'] if p.endswith('.rlib'));artifact=d
rows=[]
for name,source,old,oldsha in [('original',R/'.local/review-mark-geometry-a/regressions.rs',R/'.local/review-mark-geometry-a/regressions',None),('adversarial',R/'.local/review-mark-geometry-round2-a-astra/adversarial.rs',R/'.local/review-mark-geometry-round2-a-astra/adversarial','7ed4bbf644e79745aff224ced379cd4976944d689f96f463d7b040ec7a29842c')]:
 dst=T/(name+'.rs');dst.write_bytes(source.read_bytes());binary=T/('new-'+name);cmd=['rustc','--test','--edition=2024',str(dst),'-L','dependency='+str(Path(json.loads((T/'native-binaries.json').read_text())['rust-build-meta']['build-directory'])/'dev-optim/deps'),'--extern','kurvst='+str(rlib),'-o',str(binary)];c=subprocess.run(cmd,capture_output=True);(T/(name+'-compile.log')).write_bytes(c.stdout+c.stderr);assert c.returncode==0,c.stderr.decode()
 if oldsha:assert hashlib.sha256(old.read_bytes()).hexdigest()==oldsha
 for label,b,expect in [('old',old,101),('new',binary,0)]:
  result=subprocess.run([str(b),'--nocapture'],capture_output=True);(T/(name+'-'+label+'.log')).write_bytes(result.stdout+result.stderr);assert result.returncode==expect,(name,label,result.stdout.decode());rows.append(dict(name=name,version=label,source_sha256=hashlib.sha256(dst.read_bytes()).hexdigest(),binary=str(b),binary_sha256=hashlib.sha256(b.read_bytes()).hexdigest(),exit=result.returncode,qualification='Old failures intentionally reproduced; only current new pass gets pass credit.'))
(T/'regressions.json').write_text(json.dumps(dict(current_Cargo_artifact=artifact,current_rlib_sha256=hashlib.sha256(rlib.read_bytes()).hexdigest(),results=rows,source='Unchanged original4/five-test adversarial sources; actual new main-tree kurvst owner compiled/linked with normal source-aware Cargo caches, no seeding/overrides.'),indent=2)+'\n');print('Original4 and adversarial5 old-fail/new-pass on actual current owner')
