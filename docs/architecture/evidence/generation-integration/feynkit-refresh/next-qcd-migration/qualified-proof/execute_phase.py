"""Exact once/order controller; preserves exported private env only for supervisor.
No private reader here. Pure receipt children always remove all four names.
"""
import datetime
import json
import os
from pathlib import Path
import subprocess
import sys
P=Path(__file__).resolve().parent
W=P.parents[1]
O=W/'target/raw-next-qcd-proof-preparation-20261031-r3'
PYTHON='/nix/store/4168xqs4l5jbfc7xyi9zrhdhszx0pnlb-python3-3.13.16/bin/python3'
PRIVATE={'SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE'}

def load(p):return json.loads(Path(p).read_bytes())
def write(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2)

assert load(P/'pure-qualification.json')['status']=='PASS'
assert load(P/'elf-qualification-corrected.json')['status']=='PASS'
assert load(P/'runtime-qualification.json')['status']=='PASS'
context=load(P/'validated-context.json')
plan=load(P/'planned-invocations.json')
assert [r['id'] for r in plan]==[r['id'] for r in context['cards']]
assert not (O/'execution').exists()
invocations=P/'invocations';invocations.mkdir()
pure_environment={k:v for k,v in os.environ.items() if k not in PRIVATE}
actual=[]
for index,row in enumerate(plan):
    ident=row['id']
    assert index==row['index']
    if index>=5:
        assert all(r['validated'] and r['status']=='BASELINE_PASS' for r in actual[:5])
    item={'index':index,'id':ident,'argv':row['argv'],'cwd':str(W),'start':datetime.datetime.now(datetime.UTC).isoformat(),'retries':0}
    write(invocations/(ident+'.actual-command.json'),item)
    # No environment serialization. Existing export is inherited unchanged;
    # the reviewed supervisor alone may perform current-user XDG fallback.
    with (invocations/(ident+'.controller.stdout')).open('xb') as out,(invocations/(ident+'.controller.stderr')).open('xb') as err:
        child=subprocess.Popen(row['argv'],cwd=W,stdout=out,stderr=err)
        code=child.wait()
    item['controller_pid']=child.pid;item['controller_exit']=code
    item['end']=datetime.datetime.now(datetime.UTC).isoformat()
    outcome=O/'execution/outcomes'/(ident+'.json')
    result=load(outcome) if outcome.is_file() else {}
    item['status']=result.get('status','NO_PROTECTED_OUTCOME')
    item['validated']=False
    if code==0 and result.get('validated_predicate') is True:
        with (invocations/(ident+'.audit.stdout')).open('xb') as out,(invocations/(ident+'.audit.stderr')).open('xb') as err:
            audit=subprocess.run([PYTHON,'-I','-S','-B',str(P/'validate_outcome.py'),ident],cwd=W,env=pure_environment,stdout=out,stderr=err)
        item['audit_exit']=audit.returncode
        item['validated']=audit.returncode==0
    actual.append(item)
    write(invocations/(ident+'.completion.json'),item)
    print(json.dumps({'index':index+1,'id':ident,'controller_exit':code,'status':item['status'],'complete_reaudit':item['validated']}),flush=True)
    if not item['validated'] or item['status'] not in ['BASELINE_PASS','EXPECTED_SEMANTIC_REJECTION']:
        write(P/'actual-invocations.json',actual)
        write(P/'phase-completion.json',{'status':'STOPPED_AT_FIRST_INFRASTRUCTURE_OR_UNVALIDATED_OUTCOME','stopped_id':ident,'attempted':len(actual),'remaining_unattempted':[r['id'] for r in plan[index+1:]],'retry':False})
        raise SystemExit(3)
write(P/'actual-invocations.json',actual)
write(P/'phase-completion.json',{'status':'EXACT_21_EXECUTED_AND_REAUDITED','attempted':21,'remaining_unattempted':[],'retry':False})
