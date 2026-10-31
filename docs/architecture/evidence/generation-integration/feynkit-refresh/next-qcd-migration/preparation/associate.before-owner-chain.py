"""Bind completed proof outputs to the exact candidate worlds; no proof execution."""
import collections
import copy
import json
import re
import context as C
import reconstruct as R
B=R.B
cache={}

def resolve(ref):
    p=ref['path']
    if p not in cache:
        b=C.read(p);cache[p]=(C.sha(b),len(b),json.loads(b,object_pairs_hook=C.unique))
    h,n,obj=cache[p]
    C.need(h==ref['sha256'] and ('bytes' not in ref or n==ref['bytes']),'proof-source-binding',p)
    for key in ref.get('json_pointer','').split('/')[1:]:
        key=key.replace('~1','/').replace('~0','~')
        obj=obj[int(key)] if isinstance(obj,list) else obj[key]
    return obj

def layout(record):
    p,cp,pi=(record[k] for k in ('ports','canonical_ports','canonical_to_logical'))
    C.need(all(type(i)is int for i in pi) and sorted(pi)==list(range(len(p))),'actual-output-axis-bijection')
    C.need(cp==[p[i] for i in pi],'actual-output-axis-map')
    return pi

def owner(rec,native):
    C.need(native['owner']==rec['owner'],'actual-owner')
    pi=rec['permutation'];C.need(sorted(pi)==list(range(len(pi))),'operand-bijection')
    C.need(native['transport']==pi,'actual-operand-map')
    C.need(native['transport_same_value'] is True and native['declared_admission'] is True and native['owner_roundtrip'] is True,'actual-owner-predicates')
    maps=[]
    for route in ('historical','current'):
        C.need(native[route]['ports']==rec[route],'actual-owner-ports')
        for flag in ('is_metric','is_composite'): C.need(type(native[route][flag])is bool and native[route][flag]==rec[flag],'actual-owner-flags')
        maps.append(layout(native[route]))
    inverse=[pi.index(i) for i in range(len(pi))]
    C.need(maps[0]==[inverse[i] for i in maps[1]],'operand-output-commuting-square')
    C.need(native['historical']['canonical_ports']==native['current']['canonical_ports'],'complete-canonical-fields')
    affected=pi!=list(range(len(pi)))
    C.need(native['affected'] is affected and native['same_value_without_transport'] is (not affected),'actual-positional-limit')
    if affected:
        for route in ('historical','current'):
            C.need(native[route+'_lowering']['dangling']==sorted(rec[route]),'actual-lowering-output')
    return affected

def slot(s):
    t=B.parse(s);dual=t[0]=='call' and t[1]=='spenso::{spenso::lower}::dind'
    if dual: t=t[2][0]
    C.need(t[0]=='call','named-slot')
    return (t[1],t[2]),dual

def multiply(left,right):
    pairs=[]
    for i,s in enumerate(left):
        k,d=slot(s);hits=[]
        for j,t in enumerate(right):
            kk,dd=slot(t)
            if k==kk and ('::cof' not in k[0] or d!=dd): hits.append(j)
        C.need(len(hits)<=1,'unique-named-contraction')
        if hits: pairs.append([i,hits[0]])
    C.need(len({j for i,j in pairs})==len(pairs),'injective-named-contraction')
    usedl={i for i,j in pairs};usedr={j for i,j in pairs}
    survivors=[s for i,s in enumerate(left) if i not in usedl]+[s for j,s in enumerate(right) if j not in usedr]
    return pairs,sorted(survivors)

def plan(record,left,right,pairs):
    if pairs:
        C.need(record['kind']=='Contract','actual-plan-kind')
        C.need(sorted(record['names'])==sorted([[left[i],right[j]] for i,j in pairs]),'actual-plan-endpoints')
    else:
        C.need(record['kind']==('Scalar' if not left or not right else 'Outer'),'actual-plan-kind')

def cell_identity(c,n):
    C.need(type(n['case'])is str and n['case']==c['case'],'proof-case')
    C.need(type(n['mode'])is str and n['mode']==c['mode'],'proof-mode')

def main():
    source=C.load(C.R1/'prepared/interface-inputs.json')['cells']
    result_path=C.R3/'execution/records/all-six-cells--typed8eb/results.json'
    native=C.load(result_path);nref=C.ref(result_path)
    C.need(native['schema']=='raw-three-native-interface/1' and len(source)==len(native['cells'])==6,'actual-six-cell-domain')
    consumers=C.load(C.R1/'prepared/consumer-obligations.json')
    owner_rows=[];graph_rows=[];outputs=[];counts=collections.Counter()
    for ci,(c,n) in enumerate(zip(source,native['cells'])):
        cell_identity(c,n)
        actual=C.load(C.CAP/f"prepared/captured/{c['case']}/current-{c['mode']}.json")
        C.need(c['model']==actual['live']['model'] and c['actual_request']==actual['live']['request'] and c['options']==actual['live']['options'] and c['Run']==actual['live']['Run'],'proof-physical-inputs')
        C.need(len(c['graphs'])==len(n['graphs'])==len(actual['payload']['diagrams']),'actual-graph-domain')
        for gi,(g,ng,pg) in enumerate(zip(c['graphs'],n['graphs'],actual['payload']['diagrams'])):
            C.need(g['id']==ng['id']==pg['id'] and g['name']==ng['name']==pg['name'] and g['full_key']==pg['key'] and g['current_payload']==pg['payload'],'exact-proof-graph-key')
            for f in ('whole_graph_import','whole_graph_roundtrip','fragment_numerator_multiply_many_same_value','public_numerator_of_same_value'): C.need(ng[f] is True,'actual-graph-predicate',f)
            C.need(len(g['owners'])==len(ng['owners']),'actual-owner-domain')
            by={r['owner']:r for r in g['owners']}
            for oi,(rec,nr) in enumerate(zip(g['owners'],ng['owners'])):
                affected=owner(rec,nr);counts['changed_owner_lowerings' if affected else 'unchanged_conditional_inputs']+=1
                counts['owners']+=1
                public=B.owners(pg['payload'])[rec['owner']]
                C.need(rec['expression']==public[0] and rec['explicit']==(public[1] is not None),'proof-public-owner')
                if public[1] is not None: C.need(rec['current']==public[1],'proof-public-interface')
                pointer=f'/cells/{ci}/graphs/{gi}/owners/{oi}'
                owner_rows.append({'case':c['case'],'mode':c['mode'],'id':g['id'],'owner':rec['owner'],'full_key':g['full_key'],'actual_native':{**nref,'json_pointer':pointer},'supplied':{'path':str(C.R1/'prepared/interface-inputs.json'),'json_pointer':pointer},'operand_map':nr['transport'],'historical_output_map':nr['historical']['canonical_to_logical'],'current_output_map':nr['current']['canonical_to_logical'],'changed':affected,'lowering_scope':'Actual retained native dangling/layout output' if affected else 'Identical supplied input and deterministic exact source; no separate native lowering output.'})
            exp=next(e for e in consumers if (e['case'],e['mode'],e['id'])==(c['case'],c['mode'],g['id']))
            C.need(exp['full_key']==g['full_key'] and len(ng['plans'])==len(exp['owner_pair_inputs']),'all-selected-owner-pairs')
            for p,e in zip(ng['plans'],exp['owner_pair_inputs']):
                C.need((p['left'],p['right'])==(e['left'],e['right']),'selected-pair-owner')
                left,right=by[p['left']]['current'],by[p['right']]['current'];pairs,survivors=multiply(left,right)
                C.need(pairs==e['current']['pairs'] and survivors==e['current']['neutral_output'],'prepared-current-pair')
                plan(p['physical_plan'],left,right,pairs);counts['selected_pairs']+=1
            previous=[];seen=[];fold=[]
            for si,s in enumerate(ng['actual_fold']):
                name=s['owner'];C.need(name in by and ':' in name and name not in seen,'actual-fold-owner');seen.append(name)
                pairs,survivors=multiply(previous,by[name]['current']);plan(s['physical_plan'],previous,by[name]['current'],pairs)
                C.need(sorted(s['output']['ports'])==survivors and s['same_value'] is True,'actual-fold-output')
                layout(s['output']);actualports=s['output']['ports'];outputmap=[actualports.index(x) for x in survivors]
                C.need(sorted(outputmap)==list(range(len(survivors))),'actual-fold-survivor-map')
                fold.append({'step':si,'owner':name,'contraction_pairs':pairs,'neutral_to_actual_output':outputmap,'canonical_to_actual_output':s['output']['canonical_to_logical'],'actual_output_ports':actualports})
                previous=actualports;counts['fold_steps']+=1
            C.need(set(seen)=={x for x in by if ':' in x},'complete-actual-fold')
            C.need(sorted(previous)==sorted(by['numerator']['current']),'rank8-numerator-output');C.need(len(previous)==8,'rank8')
            layout(ng['public_numerator_of']);C.need(ng['public_numerator_of']['ports']==previous,'actual-public-output')
            C.need(by['numerator_prefactor']['current']==[] and len(by['projector']['current'])==4,'actual-prefactor-projector')
            pp,full=multiply(previous,by['projector']['current']);C.need(len(full)==4 and all('::coad(8,' in s for s in full),'four-coad8-full-output')
            outputs.append({'case':c['case'],'mode':c['mode'],'id':g['id'],'native_graph':{**nref,'json_pointer':f'/cells/{ci}/graphs/{gi}'},'fold_outputs':fold,'projector_pairs':pp,'full_projected_neutral_output':full,'full_product_scope':'Logical finite named contraction from actual rank8 numerator and captured rank4 projector. No scalar rankzero or new component execution claim.'})
            graph_rows.append({'case':c['case'],'mode':c['mode'],'id':g['id'],'full_key':g['full_key'],'owners':len(by),'actual_output_map_required':True})
    C.need((counts['owners'],counts['changed_owner_lowerings'],counts['unchanged_conditional_inputs'],counts['selected_pairs'])==(600,210,390,1890),'exact-finite-native-association-counts',counts)
    math=C.load(C.N/'math-member-to-native-row.json');partitions=collections.defaultdict(set)
    for m in math:
        result=resolve(m['result']);ob=resolve(m['obligation']);resolve(m['definition_order'])
        C.need(m['native_predicates']==result and set(result)=={'parsed_equal','canonical_equal','sign_rejected','coefficient_rejected','power_rejected'} and all(v is True for v in result.values()),'exact-native-math-row')
        C.need(m['obligation_kind']==ob['obligation'],'exact-math-obligation-kind')
        if m['obligation_kind']=='unexpanded-aggregate':
            C.need(m['member'] is None and ob['case'] in C.CASES and ob['mode'] in C.MODES,'finite-aggregate-scope')
            C.need(any((g['case'],g['mode'],g['id'])==(ob['case'],ob['mode'],ob['id']) for g in graph_rows),'exact-aggregate-graph')
        else:
            C.need(m['member'] in ob['owners'],'exact-math-owner-association')
            C.need(m['member']['case'] in C.CASES and m['member']['mode'] in C.MODES,'finite-math-scope')
        partitions[m['batch'],m['kernel']].add(m['row'])
    C.need(sorted((k,len(v)) for k,v in partitions.items())==sorted([((batch,kernel),count) for batch,count in [('faithfulness-historical',166),('current-faithfulness-and-equivalence',118)] for kernel in ('c78','8eb')]),'166-118-partitions')
    for k,v in partitions.items(): C.need(v==set(range(len(v))),'complete-native-pair-row-partition',k)
    indexes=C.load(C.N/'outcome-index.json');C.need(len(indexes)==21,'exact21')
    C.need(sum(x['status']=='BASELINE_PASS' for x in indexes)==5 and sum(x['semantic_negative_credit'] is True for x in indexes)==16,'reported-native-results')
    for x in indexes:
        C.need(x['native_launched'] is True and x['authenticated_execution'] is True and x['validated_predicate'] is True and x['error'] is None,'completed-actual-outcome')
        receipt=resolve(x['protected_receipt']);resolve(x['outcome'])
    C.save('actual-typed-owner-associations.json',owner_rows);C.save('actual-output-axis-associations.json',outputs)
    C.save('proof-association.json',{'status':'SCIENTIFIC_INTERFACE_REVIEW_PENDING','approval':False,'new_native_or_math_credit':0,'actual_native_report':C.ref(C.N/'report.json'),'counts':dict(counts),'graphs':graph_rows,'math_member_associations':len(math),'unique_pair_partition_sizes':[{ 'batch':k[0],'kernel':k[1],'rows':len(v)} for k,v in partitions.items()],'full_math_members':C.ref(C.N/'math-member-to-native-row.json'),'actual_typed_members':C.ref(C.N/'typed-native-associations.json'),'preparation_inputs':C.ref(C.N/'prepared-data-associations.json'),'all21_actual_outcomes':C.ref(C.N/'outcome-index.json'),'whole_definitions_and_order':C.ref(C.R3/'ordered-input-bindings.json'),'scope':'Reassociation of completed retained native outputs only; no new proof execution, no scientific acceptance, no historical c78 typed consumer. Operand AND output-axis maps checked. All58? relation/count claims remain individually bound in source records.','native_report_unresolved_obligations':C.load(C.N/'report.json')['unresolved_obligations']})
    C.save('associate-consumed-inputs.json',C.used)
    print('PASS retained associations',dict(counts),'math member references',len(math),'; no new native execution or scientific approval')

if __name__=='__main__': main()
