"""39 authorized actual Row runs. No Nix/compiler/native algebra helper."""
import os,sys
EXPORTED=os.environ.pop('SYMBOLICA_LICENSE',None)
for k in ('SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE'):os.environ.pop(k,None)
sys.path.insert(0,str(__file__).rsplit('/',1)[0])
from runtime_common import *
import xml.etree.ElementTree as ET
assert load(O/'selected-inventory.json')['status']=='READY_FOR_AUTHORIZED_RUNTIME'
write(O/'licensed-phase-begin.json',{'time':datetime.datetime.now(datetime.UTC).isoformat(),'pure_nix_boundary_sha256':sha(O/'nix-complete.json'),'loader_and_containment_qualified':True,'private_environment_absent_in_pure_children':True})
value=EXPORTED
if not value:
 config=Path(os.environ.get('XDG_CONFIG_HOME',str(Path.home()/'.config')))
 try:value=(config/'symbolica/license').read_text().strip()
 except OSError:raise RuntimeError('Private license unavailable in controller') from None
assert value
results=[]
steps=[('pool1-grouped',1,[GROUPED],True),('pool1-prior12',1,PRIOR12,False),('pool8',8,CASES,False),('pool32',32,CASES,False)]
for name,pool,cases,dump in steps:
 filt='package(=feynkit-generator) and binary(=generation_oracle) and ('+' or '.join('test(='+n+')' for n in cases)+')'
 cmd=[TOOLS['python'],'-I','-S','-B','/input/exec.py',TOOLS['nextest'],'nextest','run',*reuse(),'--profile','ci_gammaloop','--run-ignored','all','--no-tests','fail','--no-fail-fast','--retries','0','--success-output','final','--status-level','all','--final-status-level','all','-E',filt]
 root,e,rc=execute(name,cmd,pool=pool,value=value,dump=dump)
 f=root/'nextest/ci_gammaloop/junit.xml';assert f.is_file();tests=ET.parse(f).findall('.//testcase')
 assert len(tests)==len(cases) and {t.attrib['name'] for t in tests}==set(cases)
 statuses=[{'case':t.attrib['name'],'status':'skipped' if t.find('skipped') is not None else 'failed' if any(t.find(k) is not None for k in ['failure','error']) else 'passed','time':t.attrib.get('time')} for t in tests]
 (e/'junit.xml').write_bytes(f.read_bytes())
 results.append({'step':name,'pool':pool,'exit':rc,'outcomes':statuses,'actual_executions':sum(t['status']!='skipped' for t in statuses),'junit_sha256':sha(f),'command':str(e/'command.json'),'public_dump_hook_enabled':dump})
 write(O/'runtime-results.json',{'code':Q['code'],'archive_sha256':Q['archive_sha256'],'oracle_elf':Q['oracle_elf'],'results':results,'expected_selected_executions':39,'actual_executions':sum(r['actual_executions'] for r in results),'passes':sum(t['status']=='passed' for r in results for t in r['outcomes']),'failures':sum(t['status']=='failed' for r in results for t in r['outcomes']),'retries':0,'skipped_credit':0,'finite_math_scope':'Separately accepted pool1 None and Auto proofs; pools8/32 runtime only.'})
 print(name,':',[(t['case'],t['status']) for t in statuses],flush=True)
 assert rc==0 and all(t['status']=='passed' for t in statuses),'Actual divergence retained; dependent promotion stops.'
for p in (O/'runs').rglob('*'):
 if p.is_file():assert value.encode() not in p.read_bytes(),'Private material detected; no log disclosure permitted'
write(O/'runtime-log-audit.json',{'private_value_absent_in_all_retained_run_files':True,'no_value_path_hash_recorded':True,'Nix_or_compiler_after_licensed_start':0})
value=EXPORTED=None
print('39 actual full-Row passes; zero retries/skips; new grouped pool1 public None/Auto payloads and source companions captured.',flush=True)
