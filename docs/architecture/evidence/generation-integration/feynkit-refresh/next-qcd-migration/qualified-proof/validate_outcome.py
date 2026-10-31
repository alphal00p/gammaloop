"""Pure full protected receipt and semantic re-audit; no runtime or private call."""
import json
import os
from pathlib import Path
import sys
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
from qualification import O,PRIVATE,load,sha
sys.path.insert(0,str(O))
import supervisor as S
import supervisor_guard as K

def validate(ident):
    assert not PRIVATE.intersection(os.environ)
    b=load(P/'validated-context.json')
    card=next(c for c in b['cards'] if c['id']==ident)
    b['card']=card;b['contract']=b['all_contracts'][ident]
    record=Path(card['records_root'])
    receipt=S.public_json(record/'result.json')
    central=S.public_json(O/'execution/outcomes'/(ident+'.json'))
    assert receipt==central
    command=S.public_json(record/'command.json')
    binding=S._binding(b,card,'NATIVE_ASSIGNED',b['assignment'])
    for key,value in binding.items():
        assert receipt[key]==value,(key,'receipt')
        assert command[key]==value,(key,'command')
    assert command['card']==card
    assert command['protected_inputs_before']==b['protected_inputs']
    assert receipt['inputs_after']==b['protected_inputs']
    assert {p:sha(p) for p in b['protected_inputs']}==b['protected_inputs']
    assert receipt['deadline_seconds']==240 and receipt['retries']==0
    for key,value in [('native_launched',True),('authenticated_execution',True),('validated_predicate',True),('process_reaped',True),('launch_attempted',True),('all_input_hashes_unchanged',True),('synthetic',False),('fixture_only',False),('fixture_python_launched',False)]:
        assert type(receipt[key]) is type(value) and receipt[key]==value,key
    expected_artifacts={'command.json','stdout','stderr'}
    if b['contract']['classification']!='math-negative':expected_artifacts.add('results.json')
    assert set(receipt['artifact_hashes'])==expected_artifacts
    for name,digest in receipt['artifact_hashes'].items():assert sha(record/name)==digest
    life=receipt['lifecycle']
    contract=S._launch_contract(card,False)
    assert command['launch_contract']==contract
    for key,value in contract.items():assert life[key]==value
    assert type(life['pid']) is int and life['pid']>0
    assert set(life)=={'exit','timed_out','pid','cleanup','seconds','process_reaped'}|set(contract)
    assert life['process_reaped'] is True and life['timed_out'] is False and life['cleanup']['remaining_pids']==[]
    try:os.kill(life['pid'],0)
    except ProcessLookupError:pass
    else:raise AssertionError('owned outer launch PID still present')
    status,audit=S.classify(b,life,record)
    assert receipt['status']==status
    assert receipt['audit']==audit
    assert receipt['semantic_negative_credit']==(status=='EXPECTED_SEMANTIC_REJECTION')
    if b['contract']['dependencies']:S._check_dependencies(b,False)
    return {'id':ident,'status':status,'validated_predicate':True,'authenticated_execution':True,
            'semantic_negative_credit':receipt['semantic_negative_credit'],'full_semantic_audit':audit,
            'receipt_sha256':sha(record/'result.json'),'outcome_sha256':sha(O/'execution/outcomes'/(ident+'.json')),
            'artifact_hashes':receipt['artifact_hashes'],'process_reaped':True,
            'runtime_private_pipe_observation':'Native process reached declared source predicate through unchanged reviewed anonymous-pipe bridge; source premise, no inner PID/secret attestation.'}

if __name__=='__main__':
    ident=sys.argv[1]
    value=validate(ident)
    with (P/'invocations'/(ident+'.audit.json')).open('x') as stream:json.dump(value,stream,indent=2)
    print(json.dumps({k:value[k] for k in ['id','status','authenticated_execution','semantic_negative_credit','process_reaped']}))
