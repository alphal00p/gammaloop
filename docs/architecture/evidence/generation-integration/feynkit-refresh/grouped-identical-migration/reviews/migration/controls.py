from pathlib import Path
import copy,json,hashlib,traceback
O=Path(__file__).resolve().parent;W=O.parent.parent
Q=W.parent/'generation-feynkit-refresh/.local/feynkit-refresh/grouped-none-qualification-20261030T191919'
rows=[json.loads(x)for x in (Q/'observer-runtime/captures/ledger.jsonl').read_bytes().splitlines()]
model=json.loads((W/'target/xs-qed-l2/runtime/current/captures/xs_qed_eemm_L2/classes-None-t1/model.json').read_bytes())
def one(r,s):
    v=[x['value']for x in r if x['stage']==s];assert len(v)==1;return v[0]
def need(x,s):
    if not x:raise AssertionError(s)
def validate(r):
    live=one(r,'live-inputs');need(live['model']==model,'independently-approved-model');need(live['request']['classes']=='None','separate-mode')
    raw=one(r,'raw-diagrams');census=one(r,'final-census');groups=one(r,'collapse-entry')['groups'];pairs=one(r,'group-entry')['pairs'];order=one(r,'staged-order')['order'];final=one(r,'final-diagrams')
    for g,cg,d in zip(groups,census['groups'],final,strict=True):
        for m,cm in zip(g['members'],cg['members'],strict=True):
            idx=m['source_index'];need(order[m['diagram']]==idx and cm['source_index']==idx and raw[idx]['id']==cm['id'],'complete-membership-binding')
            need(m['overall_factor']==cm['overall_factor']==raw[idx]['json']['overall_factor'],'source-member-factor-binding')
        p=pairs[g['master']];rep=p['representative']if p['kind']=='distinct'else p['comparison']
        ignored=['name','id','cuts','topology_threshold_candidates','overall_factor','loop_momentum_basis']
        need({k:v for k,v in d['json'].items()if k not in ignored}=={k:v for k,v in rep['json'].items()if k not in ignored},'physical-representative-binding')
        for o in d['owners']:
            name=o['owner'];t=o['tensor']
            if ':'in name:
                kind,idx=name.split(':');p=d['json']['vertices'][int(idx)]if kind=='vertex'else d['json']['edges'][int(idx)][1];ik='num_interface'
            else:p=d['json'];ik={'numerator':'num_interface','projector':'projector_interface','numerator_prefactor':'numerator_prefactor_interface'}[name]
            if ik in p:need(p[ik]==t['ports'],'declared-ordered-owner-binding')
validate(rows);results=[]
def neg(name,pred,change):
    r=copy.deepcopy(rows);change(r);need(one(r,'final-census')==one(rows,'final-census'),'census-preserved-control')
    try:validate(r)
    except AssertionError as e:
        need(str(e)==pred,'intended-predicate: '+str(e));results.append({'name':name,'predicate':pred,'rejected':True,'complete_final_census_unchanged':True})
    else:raise AssertionError('control not rejected: '+name)
def model_change(r):
    one(r,'live-inputs')['model']['name']='wrong-independent-model'
    for g in one(r,'raw-diagrams')+one(r,'final-diagrams'):g['json']['model']='wrong-independent-model'
neg('mirrored model','independently-approved-model',model_change)
neg('wrong mode','separate-mode',lambda r:one(r,'live-inputs')['request'].__setitem__('classes','Auto'))
def member_change(r):
    g=one(r,'collapse-entry')['groups'][0];g['members'][0]['source_index']=19;g['members'][0]['diagram']=one(r,'staged-order')['order'].index(19)
    next(x['value']for x in r if x['stage']=='collapse-plan')['group']=[copy.deepcopy(g)]
neg('mirrored membership','complete-membership-binding',member_change)
def weight_change(r):
    g=one(r,'collapse-entry')['groups'][0];g['members'][0]['overall_factor']='-1'
    next(x['value']for x in r if x['stage']=='collapse-plan')['group']=[copy.deepcopy(g)]
neg('mirrored member factor','source-member-factor-binding',weight_change)
neg('physical numerator owner','physical-representative-binding',lambda r:one(r,'final-diagrams')[0]['json']['vertices'][0].__setitem__('numerator','1'))
def order_change(r):
    for d in one(r,'final-diagrams'):
        for o in d['owners']:
            if o['owner'].startswith('edge:')and 'num_interface'in d['json']['edges'][int(o['owner'].split(':')[1])][1] and len(o['tensor']['ports'])>1:
                o['tensor']['ports'][0],o['tensor']['ports'][1]=o['tensor']['ports'][1],o['tensor']['ports'][0];return
    raise AssertionError('missing order fixture')
neg('order-only captured declaration','declared-ordered-owner-binding',order_change)
(O/'fresh-controls.json').write_text(json.dumps({'baseline':'PASS','six_focused_controls':results,'additional_output_axis_control':'finite-named-consumer.json','scope':'Pure finite proof-binding negatives; no native or producer execution'},indent=2)+'\n')
print(json.dumps({'baseline':'PASS','negative_controls':len(results),'all_final_censuses_unchanged':True}))
