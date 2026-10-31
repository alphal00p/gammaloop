"""Wait for the already active evaluation, then invoke the single authorized realization."""
import json,os,subprocess,time
from pathlib import Path
O=Path(__file__).resolve().parent
for suffix in ['stdout','stderr']:(O/('realization.'+suffix)).rename(O/('realization-before-evaluation-finished.'+suffix))
(O/'setup-errors.json').write_text(json.dumps([{'step':'bare python3 lookup','effect':'shell python3 unavailable; reused already qualified pinned Python; no task script executed'},{'step':'initial producer invocation','error':'selected-drv.txt absent because evaluation was still active','artifact':'realization-before-evaluation-finished.stderr','effect':'No Nix build invocation or runtime. Existing evaluation reused; wait until it fully completes.'}],indent=2)+'\n')
for i in range(1200):
 if (O/'selected-drv.txt').exists() and (O/'producer-delta.json').exists():break
 if (O/'evaluation.stderr').stat().st_size:raise RuntimeError('Original evaluation failed; retained evaluation.stderr')
 time.sleep(0.5)
else:raise RuntimeError('Existing evaluation did not complete within600s; no build launched')
assert all(json.loads((O/'producer-delta.json').read_text())['invariant_env_fields'].values())
env={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')}
with (O/'realization.stdout').open('wb') as out,(O/'realization.stderr').open('wb') as err:
 rc=subprocess.call(['/nix/store/4168xqs4l5jbfc7xyi9zrhdhszx0pnlb-python3-3.13.16/bin/python3','-B',str(O/'pure.py'),'build'],cwd=O.parents[2],env=env,stdout=out,stderr=err)
assert rc==0, 'Single Nix build failed; retained artifacts'
print('One authorized producer build completed')
