"""Separate exact Auto inputs and physical correspondence; approval stays pending."""
import copy
import json
from pathlib import Path
import runpy
import traceback

O=Path(__file__).resolve().parent
U=runpy.run_path(str(O/'inspect.py'))
W,j,read,save=(U[k] for k in ['W','j','read','save'])
B=runpy.run_path(str(O/'bind_evidence.py'))
need,one,owners,project=(B[k] for k in ['need','one','owners','project'])

def main():
    aroot=W/'target/grouped-auto-fresh8eb-20261031-r1'
    ap=aroot/'observer-runtime/output/captures/ledger.jsonl'
    pp=W.parent/'generation-feynkit-refresh/.local/feynkit-refresh/qualification/runtime12/observer/xs_qed_eemm_L2/world.jsonl'
    ar=[json.loads(x) for x in read(ap).splitlines()]
    nr=[json.loads(x) for x in read(B['P']/'observer-runtime/captures/ledger.jsonl').splitlines()]
    pr=[json.loads(x) for x in read(pp).splitlines()]
    live=one(ar,'live-inputs')
    plive=next(x['value'] for x in pr if x['stage']=='batch-live-inputs' and x['value']['request']['classes']=='Auto')
    for k in ['model','options','Run','prepared_len','request']:
        need(live[k]==plive[k],'Auto-independent-raw30-input '+k)
    need(live['request']['classes']=='Auto','actual-Auto-request')
    pc=next(x['value']['census'] for x in pr if x['stage']=='batch-raw-census' and x['value']['mode']=='Auto')
    pg=next(x['value']['diagrams'] for x in pr if x['stage']=='batch-raw-diagrams' and x['value']['mode']=='Auto')
    need(one(ar,'raw-census')==pc,'Auto-independent-complete-census')
    old={x['id']:x for x in pg}
    actual=one(ar,'raw-diagrams')
    need(len(old)==len(actual)==30,'Auto-exact30')
    for g in actual:
        expected=copy.deepcopy(old[g['id']]);expected['key']=g['key']
        need(expected==g,'Auto-complete-540-owners-record')
        owners(g)
    # Full canonical key joins derive the actual permutation; no copied ordinal policy.
    canon=B['canon']
    nby={canon(c['key']):(i,c) for i,c in enumerate(one(nr,'raw-census')['diagrams'])}
    nc=one(nr,'raw-census')['diagrams'];ng=one(nr,'raw-diagrams')
    records=[]
    for ai,(g,c) in enumerate(zip(actual,one(ar,'raw-census')['diagrams'],strict=True)):
        ni,n=nby[canon(c['key'])]
        need(c['id']==n['id'] and c['payload']==n['payload'],'Auto-None-key-native-ID-physical-source')
        need(owners(g)==owners(ng[ni]),'Auto-None-complete-ordered-raw-owners')
        records.append({'Auto_raw_index':ai,'None_raw_index':ni,'exact_native_ID':c['id'],'full_key_sha256':B['digest'](c['key']),'payload_sha256':B['digest'](c['payload'])})
    mapping={x['Auto_raw_index']:x['None_raw_index'] for x in records}
    phases=[]
    np={p['raw_index']:p for p in [x['value'] for x in nr if x['stage']=='staged-pair-final']}
    for p in [x['value'] for x in ar if x['stage']=='staged-pair-final']:
        ni=mapping[p['raw_index']]
        for side in ['comparison','representative']:
            a,n=p['pair'][side],np[ni]['pair'][side]
            need(project(a)==project(n) and owners(a)==owners(n),'Auto-None-exact-phase-physical-owners')
            phases.append({'Auto_raw_index':p['raw_index'],'None_raw_index':ni,'side':side,'Auto_native_ID':a['id'],'None_native_ID':n['id'],'physical_payload_sha256':B['digest'](project(a))})
    save('Auto-separate-exact-rebind.json',{'status':'AUTO_AGGREGATION_REVIEW_PENDING','independent_Auto_raw30_path':str(pp),'independent_raw30_inputs_exact':True,'raw_sources':records,'phase_correspondence':phases,'raw_owners':540,'phase_owners':1080,'counts_reused_native_Auto_obligations':{'owner_covariance':1080,'aggregate':119,'grammar':325,'matrix_coefficients':1024,'faithfulness':76,'total':2624},'native_reuse':'Only existing exact program/input/result admissions, as bound by separate Auto ScopeA/ScopeB and native-binding. No native execution here. Historical and None native portions are not Auto approval.','rounding':'This new finite identity join uses exact native hexIDs/full keys. Existing validators compare rounded json.id redundantly; no universal rounded-ID independence.','None_approval_for_Auto':False})
    print(json.dumps({'Auto_rebind':'exact separate inputs','phase_graphs':len(phases),'Auto_aggregation':'PENDING'}))

if __name__=='__main__':
    try:main()
    except BaseException:
        with (O/'Auto-rebind-failed-attempts.log').open('a') as f:f.write(traceback.format_exc())
        raise
