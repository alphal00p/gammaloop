from pathlib import Path
import json,hashlib,shutil,subprocess
out=Path('target/mark-geometry-round3');dest=Path('.local/mark-geometry-round3-producers');dest.mkdir(exist_ok=True);sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();rows=[]
for filename in ['native-test-binaries.json','clinnet-binaries.json']:
 for key,v in json.loads((out/filename).read_text())['rust-binaries'].items():
  p=Path(v['binary-path']);d=dest/p.name;shutil.copyfile(p,d);d.chmod(0o755);rows.append({'binary_id':key,'source':str(p),'retained':str(d.resolve()),'sha256':sha(d)})
artifacts=[json.loads(l) for l in (out/'native-build.json').read_text().splitlines() if l.startswith('{')];crates=[]
for name in ['kurvst','kurbo']:
 a=next(a for a in artifacts if a.get('reason')=='compiler-artifact' and a['target']['name']==name);p=Path(next(f for f in a['filenames'] if f.endswith('.rlib')));d=dest/('lib'+name+'.rlib');shutil.copyfile(p,d);crates.append({'package':a['package_id'],'retained':str(d.resolve()),'sha256':sha(d),'build_artifact':a})
for name in ['new-probe','new-adversarial','new-regressions']:
 p=out/name;d=dest/name;shutil.copyfile(p,d);d.chmod(0o755);rows.append({'binary_id':name,'retained':str(d.resolve()),'sha256':sha(d)})
link=dest/'wasm-bundle';subprocess.run(['nix-store','--add-root',str(link.absolute()),'--indirect','--realise',str((out/'final-result').resolve())],check=True,capture_output=True)
(out/'retained-producers.json').write_text(json.dumps({'binaries':rows,'rlibs':crates,'wasm_gc_root':str(link.absolute()),'wasm_output':str((out/'final-result').resolve()),'qualification':'Task-owned copies retain source-qualified native producers; shared Cargo outputs remain mutable and are not claimed immutable.'},indent=2)+'\n')
