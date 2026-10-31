"""Exact readonly native mounts with separate NEW work; no private acquisition."""
import json
import os
from pathlib import Path
import sys
P=Path(__file__).resolve().parent
W=P.parents[1]
O=W/'target/raw-next-qcd-proof-preparation-20261031-r3'
sys.path.insert(0,str(P))
from qualification import load,save,sha,PYTHON,PRIVATE
sys.path.insert(0,str(O))
from lifecycle import bounded

assert not PRIVATE.intersection(os.environ)
assert load(P/'pure-qualification.json')['status']=='PASS'
elf=load(P/'elf-qualification-corrected.json')
assert elf['status']=='PASS'
cards=load(O/'native-commands.json')['commands']
data=load(O/'unchanged-r1-prepared-data.json')
summary=[]
for config in elf['configurations']:
    card=next(c for c in cards if c['id']==config['id'])
    root=P/'preflight'/card['id']
    work=root/'work'
    root.mkdir(parents=True)
    work.mkdir()
    for sub in ['home','tmp','cache']:(work/sub).mkdir()
    original=card['argv'];divider=original.index('--')
    prefix=[str(work) if s==card['writable_root'] else str(work/'home') if s==str(Path(card['writable_root'])/'home') else s for s in original[:divider+1]]
    assert prefix.count(str(work))==1
    interpreter=config['objects'][card['program']['path']]['interpreter']
    cfg={'id':card['id'],'readonly_mounts':config['readonly_mounts'],'input':card['native_argv'][1],
         'files':{'/input/program':card['program']['sha256']}}
    for rel,meta in data.items():
        cfg['files']['/data/'+str(Path(rel).relative_to('prepared'))]=meta['sha256']
    for path,info in config['objects'].items():
        if path!=card['program']['path']:cfg['files'][path]=info['sha256']
    cfg['files']['/input/exec.py']=sha(O.parent/'raw-next-qcd-proof-preparation-20261031-r1/prepared/runtime-input/exec.py')
    invocations=[('loader-verify',prefix+[interpreter,'--verify','/input/program']),
                 ('stdlib-containment',prefix+[PYTHON,'-I','-S','-B','-c',(P/'containment_probe.py').read_text(),json.dumps(cfg,separators=(',',':'))])]
    for label,argv in invocations:
        with (root/(label+'.command.json')).open('x') as stream:
            json.dump({'argv':argv,'cwd':str(O),'environment':card['environment_before_private_exec'],'scientific_helper_executed':False,'license_acquired':False},stream,indent=2)
        with (root/(label+'.stdout')).open('xb') as out,(root/(label+'.stderr')).open('xb') as err:
            result=bounded(argv,card['environment_before_private_exec'],O,out,err,seconds=30)
        with (root/(label+'.lifecycle.json')).open('x') as stream:json.dump(result,stream,indent=2)
        summary.append({'id':card['id'],'probe':label,'lifecycle':result,'root':str(root)})
        if result['exit']!=0 or result['timed_out'] or not result['process_reaped'] or 'exception_type' in result:
            save('runtime-qualification.json',{'status':'FAIL','probes':summary,'scientific_helper_executions':0,'private_acquisitions':0})
            raise SystemExit(3)
        print(card['id'],label,'PASS',flush=True)
    assert load(work/'probe.json')['status']=='PASS'
    assert not Path(card['writable_root']).exists()
    assert not Path(card['records_root']).exists()
save('runtime-qualification.json',{'status':'PASS','probes':summary,'scientific_helper_executions':0,'private_acquisitions':0,'private_pipe_runtime':'Pending first authorized native baseline; no dummy license used.'})
assert not (O/'execution').exists()
