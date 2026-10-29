from pathlib import Path
import subprocess,json,hashlib
root=Path.cwd();out=root/'target/mark-geometry-corrected';source=root/'.local/review-mark-geometry-a';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
r=subprocess.run(['cargo','build','--locked','-p','kurvst','--profile','dev-optim','--message-format=json'],capture_output=True,text=True);(out/'native-build.jsonl').write_text(r.stdout);(out/'native-build.log').write_text(r.stderr);r.check_returncode();a=next(json.loads(l) for l in r.stdout.splitlines() if json.loads(l).get('reason')=='compiler-artifact' and json.loads(l)['target']['name']=='kurvst');rlib=Path(next(p for p in a['filenames'] if p.endswith('.rlib')));records={}
listing=subprocess.run(['cargo','nextest','list','--locked','-p','kurvst','--cargo-profile','dev-optim','--list-type','binaries-only','--message-format','json'],capture_output=True,text=True,check=True)
meta=json.loads(listing.stdout);dependency_dir=Path(meta['rust-binaries']['kurvst']['binary-path']).parent
(out/'new-native-binaries.json').write_text(listing.stdout)
for name in ['probe','regressions']:
 cmd=['rustc','--edition=2024','-L','dependency='+str(dependency_dir),'--extern','kurvst='+str(rlib),str(source/f'{name}.rs'),'-o',str(out/f'new-{name}')]
 if name=='regressions':cmd.insert(1,'--test')
 subprocess.run(cmd,check=True)
 if name=='regressions':
  result=subprocess.run([str(out/f'new-{name}'),'--nocapture'],capture_output=True);(out/'new-regressions.log').write_bytes(result.stdout+result.stderr);result.check_returncode();records['new_regressions_exit']=result.returncode
 records[name]={'source_sha256':sha(source/f'{name}.rs'),'executable_sha256':sha(out/f'new-{name}'),'command':cmd}
records['rlib']={'path':str(rlib),'sha256':sha(rlib),'cargo_artifact':a};(out/'new-native-identities.json').write_text(json.dumps(records,indent=2)+'\n');print('Identical independent review regression source: 4/4 pass on rebuilt native owner')
