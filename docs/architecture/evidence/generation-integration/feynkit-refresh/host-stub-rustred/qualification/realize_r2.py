from pathlib import Path
import os,subprocess,json,time,datetime
P=Path(__file__).resolve().parent
E={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')}
outputs=json.loads((P/'producer-outputs-r2.json').read_text());assert (P/'source-bindings-r2.json').exists()
argv=['/run/current-system/sw/bin/nix','build','--no-link','--keep-going','--json',*[outputs[k]['drv']+'^out' for k in ('ordinary','private')]]
t=time.monotonic();start={'argv':argv,'cwd':'/tmp','started':datetime.datetime.now(datetime.timezone.utc).isoformat(),'private_variables_absent':True,'normal_substituters':True,'resource_policy':'unchanged','concurrent_characterization':'timings are not isolated performance evidence'}
(P/'build-r2-start.json').write_text(json.dumps(start,indent=2)+'\n')
with (P/'build-r2.stdout').open('wb') as out,(P/'build-r2.stderr').open('wb') as err:
 p=subprocess.Popen(argv,cwd='/tmp',env=E,stdout=out,stderr=err)
 (P/'build-r2-process.json').write_text(json.dumps({'pid':p.pid,'parent':os.getpid()})+'\n');rc=p.wait()
(P/'build-r2-completion.json').write_text(json.dumps({**start,'exit_code':rc,'elapsed_s':time.monotonic()-t,'reaped':True,'finished':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2)+'\n')
print(json.dumps({'build_exit':rc,'elapsed_s':time.monotonic()-t}),flush=True)
raise SystemExit(rc)
