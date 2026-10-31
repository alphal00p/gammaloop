from pathlib import Path
import os,subprocess,json,hashlib,shutil,time,sys,ast,csv,io,base64
P=Path(__file__).resolve().parent;W=P.parents[2];R=W/'.local/review-host-browser-stub-astra-20261031-013106Z-5655cc61';B=P.parent/'host-browser-correction-20261031-r1'
PY='/nix/store/wj1f56qky4pmjc7l0r6b6l7whnnw7hp6-python3-3.13.15/bin/python3';BW='/nix/store/dk9qhjgg469lv6mriys7v4c59igarmvx-bubblewrap-0.11.1/bin/bwrap';UV='/nix/store/gpc5h0d71izbwdqxvvxcmwdfr7f2zflx-uv-0.12.11/bin/uv';PRIVATE=('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')
A=json.loads((P/'artifact-qualification.json').read_text());SOURCE=A['source_bindings']['full_source'];HOSTNS=os.readlink('/proc/self/ns/net');ROOT=P/'runtime';ROOT.mkdir(exist_ok=True)
def sha(p):return hashlib.file_digest(Path(p).open('rb'),'sha256').hexdigest()
def save(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n')
def env():
 d={'PATH':str(Path(PY).parent)+':/nix/store/xjl7p8dvyk2j53kqf7f43kdj4ypbxz7g-coreutils-9.11/bin','HOME':'/work/home','USER':'runtime','LOGNAME':'runtime','LANG':'C.UTF-8','LC_ALL':'C.UTF-8','TMPDIR':'/work/tmp','XDG_CACHE_HOME':'/work/cache','XDG_CONFIG_HOME':'/work/config','XDG_STATE_HOME':'/work/state','XDG_RUNTIME_DIR':'/work/runtime','PYTHONDONTWRITEBYTECODE':'1','PYTHONPATH':'/nix/store/k1kjicaywfi1bp79rpgggj2gcbyhvkn0-python3-3.13.15-env/lib/python3.13/site-packages','SYMBOLICA_DATA_DIRECTORY':'/nix/store/5p34pih2c8i4vr4yhj3xj4v7064lpxlb-gammaloop-render-data','LINNET_TEST_BROWSER_ENGINES':'chromium,webkit','LINNET_TEST_CHROMIUM':'/nix/store/ds6wq0gzyimhl7c1kmzdvql0w4d2rcv4-chromium-153.0.8010.36/bin/chromium','PLAYWRIGHT_BROWSERS_PATH':'/nix/store/h8903h8s8x5zr3la0pnbgh24sbl5pf62-playwright-browsers','PLAYWRIGHT_HOST_PLATFORM_OVERRIDE':'ubuntu-24.04','PLAYWRIGHT_SKIP_BROWSER_DOWNLOAD':'1','EXPECTED_HOST_NETNS':HOSTNS,'FONTCONFIG_FILE':'/proposal/fonts.conf'};d.update(json.loads((B/'egl-runtime.json').read_text())['env']);assert all(k not in d for k in PRIVATE);return d

def execute(lane,label,command,licensed=False,private=None,installed=True,work=None):
 root=ROOT/lane;work=work or root/'work'/label;ev=root/'evidence'/label;ev.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=True)
 for n in ('home/.config/symbolica','config/symbolica','cache','state','tmp','runtime','venv'):(work/n).mkdir(parents=True,exist_ok=True)
 argv=[BW,'--unshare-all','--die-with-parent','--new-session','--proc','/proc','--dev','/dev','--tmpfs','/tmp','--bind',str(work),'/work']
 for p in A['store_mounts']:argv+=['--ro-bind',p,p]
 argv+=['--ro-bind',SOURCE,'/build/source','--ro-bind',str(root/'helpers'),'/proposal']
 if installed:argv+=['--ro-bind',str(root/'installed/venv'),'/work/venv']
 argv+=['--ro-bind',str(work/'config/symbolica'),'/work/config/symbolica','--ro-bind',str(work/'home/.config/symbolica'),'/work/home/.config/symbolica']
 for src,dst in A['extra_browser_mounts']:argv+=['--ro-bind',src,dst]
 argv+=['--chdir','/work','--',*command];environment=env();record={'label':label,'lane':lane,'argv':argv,'environment_before_handoff':environment,'licensed':licensed,'source_revision':A['source_bindings']['new_code'],'source':SOURCE,'wheel':A['records'][lane]['wheel'],'native_sha256':A['records'][lane]['native']['sha256'],'evidence_not_mounted':str(ev),'writable':str(work),'controller_sha256':sha(__file__),'helpers':{p.name:sha(p) for p in (root/'helpers').iterdir() if p.is_file()},'timeout_retry_resource_policy':'unchanged; no new cap or retry','started_epoch':time.time()};save(ev/'start.json',record);print(json.dumps({'start':lane+'/'+label}),flush=True)
 t=time.monotonic();p=subprocess.Popen(argv,cwd='/tmp',env=environment,stdin=subprocess.PIPE if licensed else subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.PIPE);save(ev/'process.json',{'pid':p.pid,'parent':os.getpid()});out,err=p.communicate(private.encode() if licensed else None);redacted=False
 if private:
  secret=private.encode();redacted=secret in out or secret in err;out=out.replace(secret,b'[REDACTED]');err=err.replace(secret,b'[REDACTED]')
  for f in work.rglob('*'):
   if not f.is_file() or f.is_symlink():continue
   with f.open('rb') as stream:
    tail=b''
    for block in iter(lambda:stream.read(2**20),b''):
     data=tail+block;assert secret not in data,'Private material in child file';tail=data[-len(secret):]
 (ev/'stdout').write_bytes(out);(ev/'stderr').write_bytes(err);r={'exit_code':p.returncode,'elapsed_s':time.monotonic()-t,'reaped':True,'secret_redaction_required':redacted,'stdout_sha256':sha(ev/'stdout'),'stderr_sha256':sha(ev/'stderr')};save(ev/'completion.json',r);print(json.dumps({'done':lane+'/'+label,**r}),flush=True);return r,out,err

def prepare():
 assert not (ROOT/'licensed-phase-start.json').exists();assert (P/'build-r2-completion.json').exists();assert A['pure_artifact_qualification_complete']
 needed=sorted({e['needed'] for e in json.loads((P/'native-library-resolution.json').read_text())['resolved_NEEDED']})
 for lane,item in A['records'].items():
  root=ROOT/lane;helpers=root/'helpers';helpers.mkdir(parents=True,exist_ok=False)
  for f in ('handoff.py','preflight.py','product.py','installed_sector_signature.py'):
   s=(R/'helpers'/f).read_text().replace("'SYMBOLICA_LICENSE_SIGNED')","'SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')")
   if f=='product.py':s=s.replace('53c5efb6344b233dfba09612835c33940f4daf64b215462f4292335067c420a2',item['native']['sha256'])
   (helpers/f).write_text(s)
  shutil.copyfile(B/'helpers/fonts.conf',helpers/'fonts.conf');save(helpers/'mounts.json',A['store_mounts']);save(helpers/'expected.json',{'native_sha256':item['native']['sha256'],'stub_sha256':item['stub_sha256'],'private':lane=='private'});save(helpers/'cwd-library-policy.json',{'needed_names':needed,'source':SOURCE,'empty_runpath_exception':'/nix/store/wkiqjb44n5k0fwhmc523dvflspd47q3b-gcc-15.3.0-libgcc/lib/libgcc_s.so.1'})
  for cwd in (Path(SOURCE),Path(SOURCE)/'crates/feynkit-py'):assert all(not (cwd/n).exists() for n in needed)
  (helpers/'identity.py').write_text('''import hashlib,json,sys,os\nfrom pathlib import Path\nimport symbolica.core as core\nfrom symbolica.community import hepkit,feynkit_native\ne=json.loads(Path('/proposal/expected.json').read_text());p=Path(core.__file__);assert str(p).startswith('/work/venv/');assert hashlib.sha256(p.read_bytes()).hexdigest()==e['native_sha256'];assert hepkit.FeynmanDiagram is feynkit_native.FeynmanDiagram;stub=Path(hepkit.__file__).with_suffix('.pyi');assert hashlib.sha256(stub.read_bytes()).hexdigest()==e['stub_sha256'];assert hasattr(hepkit,'IBPFamily')==e['private'];assert __debug__\nPath('/work/identity.json').write_text(json.dumps({'core_path':str(p),'core_sha256':e['native_sha256'],'stub_path':str(stub),'stub_sha256':e['stub_sha256'],'IBPFamily':hasattr(hepkit,'IBPFamily'),'native_class_identity':True,'assertions_enabled':__debug__},indent=2)+'\\n');print('Exact installed wheel identity passed')\n''')
  q=(P.parent/'host-api-stub-correction-20261031-r1/helpers/qcd_boundary.py').read_text().replace("    runpy.run_path(script, run_name='__main__')","    scope = runpy.run_path(script, run_name='__main__')").replace("    receipt['complete_script_passed'] = True","    receipt.update({'targets':scope.get('targets'),'residuals':scope['solution'].residuals,'trees_count':len(scope.get('trees',{})),'loop_groups':len(scope.get('terms_by_diagram',{}))})\n    receipt['complete_script_passed'] = True")
  (helpers/'qcd.py').write_text(q)
  commands={'identity':['/work/venv/bin/python','-B','/proposal/identity.py'],'qcd':['/work/venv/bin/python','-B','/proposal/qcd.py'],'ufo':['/work/venv/bin/python','-B','/build/source/crates/feynkit-py/tests/installed_ufo_ghost_crossing.py'],'squared':['/work/venv/bin/python','-B','/build/source/crates/feynkit-py/tests/installed_squared_amplitude_display.py'],'sector':['/work/venv/bin/python','-B','/proposal/installed_sector_signature.py','/work/sector.json'],'browser':['/work/venv/bin/python','-B','/proposal/product.py'],'stub-check':[A['exporter']['path'],'feynkit-community','/build/source/docs/api/python/feynkit-community.pyi','--check']};save(helpers/'commands.json',commands)
  for f in helpers.glob('*.py'):ast.parse(f.read_text())
  installation=root/'installed';r,_,err=execute(lane,'create-venv',[UV,'venv','--python',PY,'--system-site-packages','/work/venv'],installed=False,work=installation);assert r['exit_code']==0,err.decode()
  r,_,err=execute(lane,'install-whole-wheel',[UV,'pip','install','--offline','--no-deps','--python','/work/venv/bin/python',item['wheel']['path']],installed=False,work=installation);assert r['exit_code']==0,err.decode()
  venv=installation/'venv';inv=[]
  for f in venv.rglob('*'):
   if f.is_file() and not f.is_symlink():inv.append({'path':str(f.relative_to(venv)),'sha256':sha(f),'bytes':f.stat().st_size})
  # Compare every installed wheel member, including RECORD/METADATA, with the whole selected wheel.
  site=venv/'lib/python3.13/site-packages'
  for f in item['members']:
   if not f['path'].endswith('/RECORD'):assert sha(site/f['path'])==f['sha256'],f['path']
  record_name,=[k for k in item['metadata'] if k.endswith('/RECORD')]
  original_rows=list(csv.reader(io.StringIO(item['metadata'][record_name])));installed_rows=list(csv.reader(io.StringIO((site/record_name).read_text())))
  installed_by_name={row[0]:row for row in installed_rows}
  for row in original_rows:assert installed_by_name[row[0]]==row,row
  for name,h,size in installed_rows:
   data=(site/name).read_bytes()
   if h:assert h=='sha256='+base64.urlsafe_b64encode(hashlib.sha256(data).digest()).decode().rstrip('=')
   if size:assert len(data)==int(size)
  record_observation={'wheel_RECORD_sha256':next(f['sha256'] for f in item['members'] if f['path']==record_name),'installed_RECORD_sha256':sha(site/record_name),'all_original_rows_preserved':True,'installed_rows_verified':len(installed_rows),'installer_added_members':[row[0] for row in installed_rows if row[0] not in {r[0] for r in original_rows}]}
  save(root/'installation.json',{'wheel':item['wheel'],'native':item['native']['sha256'],'files':inv,'all_wheel_payload_members_byte_identical':True,'RECORD':record_observation,'stub_sha256':sha(site/'symbolica/community/hepkit/__init__.pyi'),'native_imported':False})
  # Exact loader inferred from this native's resolved libc carrier, and executed with no activation.
  edges=json.loads((P/'native-library-resolution.json').read_text())['resolved_NEEDED'];libc=next(Path(e['resolved']) for e in edges if e['consumer']==item['native']['path'] and e['needed']=='libc.so.6');loader=libc.parent/'ld-linux-x86-64.so.2';assert loader.is_file()
  r,_,err=execute(lane,'loader-list',[str(loader),'--list','/work/venv/lib/python3.13/site-packages/symbolica/core.abi3.so']);assert r['exit_code']==0,err.decode()
  r,_,err=execute(lane,'containment-preflight',[PY,'-I','-S','-B','/proposal/preflight.py']);assert r['exit_code']==0,err.decode()
 save(ROOT/'preparation-complete.json',{'lanes':list(A['records']),'artifacts_sha256':sha(P/'artifact-qualification.json'),'all_pure_queries_builds_and_installations_complete':True,'licensed_execution_started':False,'helpers':{str(p.relative_to(ROOT)):sha(p) for p in ROOT.rglob('helpers/*') if p.is_file()}})

def licensed():
 assert (ROOT/'preparation-complete.json').exists();assert not (ROOT/'licensed-phase-start.json').exists();save(ROOT/'licensed-phase-start.json',{'epoch':time.time(),'no_more_compiler_or_nix':True,'resource_policy_unchanged':True,'artifact_qualification_sha256':sha(P/'artifact-qualification.json')})
 private=os.environ.get('SYMBOLICA_LICENSE')
 if not private:private=(Path(os.environ.get('XDG_CONFIG_HOME',str(Path.home()/'.config')))/'symbolica/license').read_text().strip()
 assert private
 results=[]
 for lane,item in A['records'].items():
  commands=json.loads((ROOT/lane/'helpers/commands.json').read_text());modes=['identity']
  if lane=='ordinary':modes+=['stub-check']
  else:modes+=['qcd','ufo','squared','sector','browser']
  for mode in modes:
   r,out,err=execute(lane,mode,[PY,'-I','-S','-B','/proposal/handoff.py',mode],licensed=True,private=private);passed=r['exit_code']==0 and not r['secret_redaction_required'];entry={'lane':lane,'mode':mode,**r,'passed':passed}
   if mode=='identity' and not passed:results.append(entry);break
   if mode in ('qcd','ufo','squared'):
    marker={'qcd':'Generated ghost and antighost UV: full tensors, crossing, two IBP targets and coupling consistency passed','ufo':'Passed 36 EW and 3 QCD ghost crossings, plus 15 Lorentz-linearity probes.','squared':'Squared-amplitude pager, display settings, navigation and static fallback passed'}[mode];entry['completion_marker_present']=marker.encode() in out
    entry['passed']=passed and entry['completion_marker_present']
    if mode=='qcd':receipt=ROOT/lane/'work/qcd/qcd-boundary.json';entry['boundary']=json.loads(receipt.read_text()) if receipt.exists() else None
   if mode=='browser':entry['product']=json.loads((ROOT/lane/'work/browser/product.json').read_text()) if (ROOT/lane/'work/browser/product.json').exists() else None;entry['passed']=passed and bool(entry['product'] and entry['product']['complete'])
   results.append(entry);save(ROOT/'results.json',results)
 del private
 save(ROOT/'results.json',results);print(json.dumps({'runtime_complete':True,'counts':{'checks':len(results),'passed':sum(x['passed'] for x in results),'failed':sum(not x['passed'] for x in results)}}),flush=True)
if __name__=='__main__':{'prepare':prepare,'licensed':licensed}[sys.argv[1]]()
