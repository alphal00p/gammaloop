import json,subprocess
from pathlib import Path
out=Path('target/mark-geometry-round3');a=[json.loads(l) for l in (out/'native-build.json').read_text().splitlines()];artifact=next(x for x in a if x.get('reason')=='compiler-artifact' and x['target']['name']=='kurvst');rlib=next(f for f in artifact['filenames'] if f.endswith('.rlib'));dep=str(Path(next(iter(json.loads((out/'binaries.json').read_text())['rust-binaries'].values()))['binary-path']).parent)
# Cargo's returned rlib can be in deps already.
if Path(rlib).parent.name=='deps':dep=str(Path(rlib).parent)
for source,name,test in [('.local/review-mark-geometry-round2-a-astra/adversarial.rs','new-adversarial',True),('target/mark-geometry-round3/probe.rs','new-probe',False),('.local/review-mark-geometry-a/regressions.rs','new-regressions',True)]:
 cmd=['rustc','--edition=2024','-L','dependency='+dep,'--extern','kurvst='+rlib,source,'-o',str(out/name)]
 if test:cmd.insert(1,'--test')
 r=subprocess.run(cmd,capture_output=True,text=True);print(r.stdout+r.stderr);r.check_returncode()
for name in ['new-adversarial','new-regressions']:
 r=subprocess.run([str(out/name),'--nocapture'],capture_output=True,text=True);(out/(name+'.log')).write_text(r.stdout+r.stderr);print(r.stdout+r.stderr);r.check_returncode()
