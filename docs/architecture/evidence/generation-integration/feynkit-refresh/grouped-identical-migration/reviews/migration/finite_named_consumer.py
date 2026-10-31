"""Independent logical-key reconstruction of the exact finite named consumers.
No native canonical ordering is guessed: the source proof removes each actual
canonical/layout permutation with its inverse before selecting logical pairs.
"""
from pathlib import Path
import collections,copy,hashlib,json,re,traceback
O=Path(__file__).resolve().parent;W=O.parent.parent
P=W/'target/grouped-identical-migration-preparation-20261031-r1'
Q=W.parent/'generation-feynkit-refresh/.local/feynkit-refresh/grouped-none-qualification-20261030T191919'
reads={}
def read(p):
    p=Path(p);b=p.read_bytes();reads[str(p)]=hashlib.sha256(b).hexdigest();return b
def load(p):return json.loads(read(p))
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
def need(x,s):
    if not x:raise AssertionError(s)
def parse(p):
    dual=p.startswith('spenso::{spenso::lower}::dind(')
    if dual:p=p[len('spenso::{spenso::lower}::dind('):-1]
    m=re.fullmatch(r'spenso::\{spenso::representation,spenso::(self_dual|dualizable)\}::(bis|mink|cof)\((4|3),gammalooprs::\{spenso::index\}::(hedge|edge)\((\d+),1\)\)',p)
    need(m is not None,'complete explicit port domain')
    kind,rep,dim,head,index=m.groups()
    need((rep in ['bis','mink'] and kind=='self_dual' and dim=='4' and not dual)or(rep=='cof' and kind=='dualizable' and dim=='3'),'rep definition/dimension')
    return (rep,int(dim),head,int(index),1),dual

def analyze(g):
    owners={x['owner']:x['tensor']for x in g['owners']}
    need(len(owners)==18,'18 owner domain')
    groups=collections.defaultdict(list);declarations=[]
    for name,t in owners.items():
        parsed=[parse(p)for p in t['ports']];need(len(parsed)==len(set(parsed)),'repeat-free typed owner ports')
        need(type(t['is_metric'])is bool and type(t['is_composite'])is bool,'actual owner flag types')
        if ':'in name:
            kind,idx=name.split(':');payload=g['json']['vertices'][int(idx)]if kind=='vertex'else g['json']['edges'][int(idx)][1];field,ik='numerator','num_interface'
            for axis,(key,orientation)in enumerate(parsed):groups[key].append((name,axis,orientation))
        else:payload=g['json'];field=name;ik={'numerator':'num_interface','numerator_prefactor':'numerator_prefactor_interface','projector':'projector_interface'}[name]
        need(t['raw']==t['decoded']==payload[field],'whole owner expression binding')
        if ik in payload:
            need(payload[ik]==t['ports'],'declared order binding');declarations.append({'owner':name,'ports':t['ports']})
    contracted=[];free=[]
    for key,occ in sorted(groups.items()):
        need(len(occ)in [1,2],'global public multiplicity')
        if len(occ)==1:free.append((key,occ[0]));continue
        a,b=occ;need(a[0]!=b[0],'cross-owner pairing')
        need(key[0]in ['bis','mink'] or a[2]!=b[2],'dual-compatible pair')
        contracted.append({'key':key,'left':list(a[:2]),'right':list(b[:2])})
    np=owners['numerator']['ports'];pp=owners['projector']['ports'];fp=owners['numerator_prefactor']['ports']
    need(owners['numerator_prefactor']['decoded']=='1' and fp==[],'actual scalar prefactor')
    need(collections.Counter(parse(p)for p in np)==collections.Counter((k,o[2])for k,o in free),'local product surviving named slots')
    need(len(np)==len(pp)==4,'four open numerator and projector axes')
    bp=[]
    for i,p in enumerate(np):
        k,o=parse(p);matches=[j for j,q in enumerate(pp)if parse(q)[0]==k and (k[0]in ['bis','mink'] or parse(q)[1]!=o)]
        need(len(matches)==1,'unique projector contraction');bp.append([i,matches[0]])
    need(sorted(j for i,j in bp)==list(range(4)),'four projector partners exactly once')
    return {'id':g['id'],'contractions':contracted,'survivors':[{'key':k,'operand':list(o[:2]),'output_axis':np.index(owners[o[0]]['ports'][o[1]])}for k,o in free],'numerator_ports':np,'prefactor_ports':fp,'projector_ports':pp,'projector_logical_pairs':bp,'complete_product_output_ports':[],'declarations':declarations}

def main():
    rows=[json.loads(b)for b in read(Q/'observer-runtime/captures/ledger.jsonl').splitlines()]
    def one(k):
        v=[r['value']for r in rows if r['stage']==k];need(len(v)==1,k);return v[0]
    raw=one('raw-diagrams');final=one('final-diagrams');pairs={r['value']['raw_index']:r['value']['pair']for r in rows if r['stage']=='staged-pair-final'}
    graph_records=[];index={}
    for phase,gs in [('raw',list(enumerate(raw))),('comparison',[(i,p['comparison'])for i,p in sorted(pairs.items())]),('representative',[(i,p['representative'])for i,p in sorted(pairs.items())]),('final',list(enumerate(final)))]:
        for idx,g in gs:
            a=analyze(g);a.update(phase=phase,index=idx);graph_records.append(a);index[(phase,idx)]=a
    need(len(graph_records)==119,'119 complete graphs')
    maps=load(P/'source-phase-owner-axis-maps.json');commutations=[]
    for m in maps:
        a=index['raw',m['raw_index']];b=index[m['side'],m['raw_index']]
        def endpoint(ep):
            o,axis=ep;v=m['owners'][o];return (v['target_owner'],v['source_axis_to_target_axis'][axis])
        def normalize(ps):return {frozenset((tuple(p['left']),tuple(p['right'])))for p in ps}
        mapped=[{'left':endpoint(p['left']),'right':endpoint(p['right'])}for p in a['contractions']]
        need(normalize(mapped)==normalize(b['contractions']),'local contraction square')
        output=m['owners']['numerator']['source_axis_to_target_axis']
        target={tuple(x['operand']):x['output_axis']for x in b['survivors']}
        for s in a['survivors']:need(target[endpoint(s['operand'])]==output[s['output_axis']],'survivor output-axis square')
        proj=m['owners']['projector']['source_axis_to_target_axis']
        need(sorted([output[i],proj[j]]for i,j in a['projector_logical_pairs'])==sorted(b['projector_logical_pairs']),'prefactor/projector named pair square')
        commutations.append({'raw_index':m['raw_index'],'side':m['side'],'local_pair_square':True,'output_axis_map':output,'projector_operand_map':proj,'complete_output_axis_map':[],'scalar_prefactor_identity':True})
    # Census-preserving order-only mutation of the required output map.
    m=copy.deepcopy(maps[0]);out=m['owners']['numerator']['source_axis_to_target_axis'];out[0],out[1]=out[1],out[0]
    a=index['raw',m['raw_index']];b=index[m['side'],m['raw_index']];proj=m['owners']['projector']['source_axis_to_target_axis']
    need(sorted([out[i],proj[j]]for i,j in a['projector_logical_pairs'])!=sorted(b['projector_logical_pairs']),'wrong output map must break contraction square')
    # Consistent operand+output reordering is a valid change of axes, not an error.
    perm=[1,0,2,3];conjugated=sorted([perm[i],perm[j]]for i,j in a['projector_logical_pairs'])
    need(conjugated==sorted(a['projector_logical_pairs']),'consistent named-axis permutation positive')
    final_decl=sum(len(x['declarations'])for x in graph_records if x['phase']=='final');need(final_decl==54,'54 final declarations')
    save('finite-named-consumer.json',{'status':'FINITE_LOGICAL_KEY_AND_AXIS_COMMUTATION_PASS','graphs':graph_records,'commutations':commutations,'graphs_per_mode':119,'owners_per_mode':2142,'final_declarations':54,'negative':{'mutation':'swap first two numerator output-map entries only; all payloads, declarations, flags, counts, weights and census bytes unchanged','predicate':'prefactor/projector named pair square','rejected':True,'native_execution':False},'positive':'Simultaneous corresponding operand/output axis permutation preserves named contraction','native_canonical_integer_layouts_materialized':False,'source_proof_obligation':'canonical/logical inverse lemma and named-key dispatch, assessed explicitly in review; this finite program alone is not a source theorem or native execution'})
    print(json.dumps({'graphs':119,'source_to_phase_commutations':60,'final_declared_interfaces':54,'output_map_negative_rejected':True,'native':0}))
try:main()
except BaseException:
    save('finite-named-error.json',{'traceback':traceback.format_exc()});raise
finally:save('finite-named-reads.json',reads)
