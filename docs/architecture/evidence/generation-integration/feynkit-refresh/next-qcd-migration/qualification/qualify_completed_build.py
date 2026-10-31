"""Qualify the already running single producer when its completed record is available."""
import json,os,subprocess,time
from pathlib import Path
O=Path(__file__).resolve().parent;P='/nix/store/4168xqs4l5jbfc7xyi9zrhdhszx0pnlb-python3-3.13.16/bin/python3'
env={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')}
for i in range(1800):
 m=json.loads((O/'build.command.json').read_text())
 if 'exit' in m:break
 time.sleep(0.5)
else:raise RuntimeError('Single producer remains incomplete; no duplicate build')
assert m['exit']==0
for script in ['qualify.py','dependency_invariance.py','finish_pure.py','runtime_prepare.py']:
 with (O/(script+'.stdout')).open('wb') as out,(O/(script+'.stderr')).open('wb') as err:
  rc=subprocess.call([P,'-B',str(O/script)],cwd=O.parents[2],env=env,stdout=out,stderr=err)
 assert rc==0,script+' failed; original diagnostics retained; no native execution'
 print(script,'PASS',flush=True)
