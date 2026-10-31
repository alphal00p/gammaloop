import sys,re,copy,collections,functools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import common as C
from rows_review import require,Rejected,validate_world,owners,encoded
sys.path.insert(0,str(C.O/'pure'))
import pure_parser as B
import strict_owners as G

cache={}
def resolve(ref):
 p=ref['path']
 if p not in cache:
  b=C.read(p);cache[p]=(C.sha(b),len(b),C.json.loads(b,object_pairs_hook=C.unique) if str(p).endswith('.json') else None)
 h,n,obj=cache[p];require(ref['sha256']==h and ('bytes' not in ref or ref['bytes']==n),'exact-reference-hash')
 for k in ref.get('json_pointer','').split('/')[1:]:
  k=k.replace('~1','/').replace('~0','~');obj=obj[int(k)] if isinstance(obj,list) else obj[k]
 return obj

def mapped(x):
 if isinstance(x,str):return B.TOKEN.sub(lambda m:B.SPELLINGS.get(m[0],m[0]),x)
 if isinstance(x,list):return [mapped(y) for y in x]
 if isinstance(x,dict):return {k:mapped(v) for k,v in x.items()}
 return x

def normal(s,h=False):return B.render(B.parse(s,h))
@functools.cache
def port(s):
 pattern=r'(spenso::\{spenso::lower\}::dind\()?spenso::\{spenso::representation,spenso::(self_dual|dualizable)\}::(bis|mink|coad|cof)\((\d+),gammalooprs::\{spenso::index\}::(hedge|edge)\((\d+),(\d+)\)\)(\))?'
 m=re.fullmatch(pattern,s);require(m is not None,'full-port-definition')
 dual,tag,rep,dim,kind,idx,lane,end=m.groups();require(bool(dual)==bool(end),'dual-syntax')
 require(int(dim)=={'bis':4,'mink':4,'coad':8,'cof':3}[rep] and int(lane)==1,'port-dimension-lane')
 require(tag==('dualizable' if rep=='cof' else 'self_dual') and (not dual or rep=='cof'),'port-variance')
 return (rep,int(dim),kind,int(idx),int(lane)),bool(dual)
def product(l,r):
 pairs=[]
 for i,a in enumerate(l):
  ak,ad=port(a)
  for j,b in enumerate(r):
   bk,bd=port(b)
   if ak==bk and (ak[0]!='cof' or ad!=bd):pairs.append([i,j])
 require(len({i for i,j in pairs})==len(pairs)==len({j for i,j in pairs}),'unique-pair')
 survivors=[s for i,s in enumerate(l) if i not in {i for i,j in pairs}]+[s for j,s in enumerate(r) if j not in {j for i,j in pairs}]
 require(len(survivors)==len(set(survivors)),'unique-output');return pairs,sorted(survivors)
def layout(rec):
 p,cp,inv=rec['ports'],rec['canonical_ports'],rec['canonical_to_logical']
 require(all(type(i)is int for i in inv) and sorted(inv)==list(range(len(p))),'output-permutation')
 require(cp==[p[i] for i in inv],'output-axis-map');return inv

def plan(debug,physical,l,r):
 pairs,out=product(l,r);kind='Scalar' if not l or not r else 'Contract' if pairs else 'Outer'
 require(physical['kind']==kind,'actual-dispatch-kind')
 if kind!='Contract':require(debug==kind,'integer-debug-plan')
 else:
  expected='Contract(['+', '.join(f'PortPair {{ left: {i}, right: {j} }}' for i,j in pairs)+'])'
  require(debug==expected,'integer-debug-plan')
  require(physical['names']==sorted([[l[i],r[j]] for i,j in pairs]),'named-endpoint-map')
 return pairs,out

def cell_input(cell,worlds):
 case,mode=cell['case'],cell['mode'];require(case in C.CASES and mode in C.MODES,'scientific-cell-domain')
 for era in ('historical','current'):
  w=worlds[era];validate_world(w,case,mode,era)
  for key,actual in [('model','model'),('actual_request','request'),('options','options'),('Run','Run'),('grouping','grouping')]:require(cell[key]==w['live'][actual],'proof-'+actual)

def owner_input(rec,native,h,c,public):
 require(rec['owner']==native['owner'],'native-owner')
 require(rec['expression']==c['decoded']==public[0] and rec['historical_expression']==h['decoded'] and rec['historical_raw_expression']==h['raw'],'owner-expression')
 require(rec['current']==c['ports'] and rec['historical']==mapped(h['ports']),'owner-ports')
 require(rec['explicit']==(public[1] is not None),'interface-presence')
 if public[1] is not None:require(public[1]==c['ports'],'public-interface')
 hp,cp=rec['historical'],rec['current'];pi=rec['permutation']
 require(sorted(pi)==list(range(len(hp))) and hp==[cp[i] for i in pi] and len(set(hp))==len(hp),'operand-axis-map')
 require(native['transport']==pi,'native-operand-map')
 for era,tensor,ports in [('historical',h,hp),('current',c,cp)]:
  require(native[era]['ports']==ports,'native-owner-ports');layout(native[era])
  for flag in ('is_metric','is_composite'):require(type(native[era][flag])is bool and native[era][flag]==tensor[flag]==rec[flag],'actual-owner-flag')
 hi,ci=native['historical']['canonical_to_logical'],native['current']['canonical_to_logical']
 require([pi[i] for i in hi]==ci and native['historical']['canonical_ports']==native['current']['canonical_ports'],'both-axis-commuting-square')
 for k in ('transport_same_value','declared_admission','owner_roundtrip'):require(native[k] is True,'native-owner-'+k)
 changed=pi!=list(range(len(pi)));require(native['affected']is changed and native['same_value_without_transport']is (not changed),'positional-restriction')
 for era,ports in [('historical',hp),('current',cp)]:
  if changed:require(native[era+'_lowering']['dangling']==sorted(ports),'lowered-output-domain')
  else:require(era+'_lowering' not in native,'unchanged-lowering-source-condition')
 return changed

def whole_definitions(world,table,era):
 defs={k:{} for k in ('raw_ports_and_factors','decoded_ports_and_factors','mapped_decoded_ports_and_factors')};faith=[]
 for d in world['payload']['diagrams']:
  for label,t in table[d['id']].items():
   for kind,texts in [('raw',[t['raw']]),('decoded',[t['decoded']]),('ports',t['ports'])]:
    for s in texts:
     B.definitions(B.parse(s),defs['raw_ports_and_factors'])
     if kind!='raw':
      B.definitions(B.parse(s),defs['decoded_ports_and_factors']);B.definitions(B.parse(s,era=='historical'),defs['mapped_decoded_ports_and_factors'])
     faith.append(((s,normal(s)),{'case':world['case'],'mode':world['mode'],'route':era,'id':d['id'],'owner':label,'kind':kind}))
  for label,s in [('overall_factor',d['overall_factor']),('payload.overall_factor',d['payload']['overall_factor'])]:
   for k in defs:B.definitions(B.parse(s,k=='mapped_decoded_ports_and_factors' and era=='historical'),defs[k])
   faith.append(((s,normal(s)),{'case':world['case'],'mode':world['mode'],'route':era,'id':d['id'],'owner':label,'kind':'source_factor'}))
 defs['whole_definition_map']={s:B.SPELLINGS.get(s,s) if era=='historical' else s for s in defs['decoded_ports_and_factors'].values()}
 defs['callback_userdata_status']='NOT captured; retained source qualification premise, no absent-to-empty inference'
 require(len(set(defs['whole_definition_map'].values()))==len(defs['whole_definition_map']),'injective-full-definition-map')
 return defs,faith

def main():
 sr=C.load(C.S/'review.json');root=C.load(C.F/'root-next-qcd-scientific-interface-acceptance-20261031.json')
 require(sr['decision']=='APPROVE' and not sr['blocking_findings'] and sr['additional_native_requirement'] is None,'accepted-scientific-review')
 require(root['review']['sha256']==C.sha(C.read(C.S/'review.json')) and root['joined_migration_preparation']['seal_sha256']==C.sha(C.read(C.P/'SEAL.json')),'additive-scientific-join')
 for name,h in sr['candidate']['pins'].items():require(C.sha(C.read(C.N/name))==h,'review-native-packet')
 inputs=C.load(C.R1/'prepared/interface-inputs.json')['cells'];native=C.load(C.R3/'execution/records/all-six-cells--typed8eb/results.json')['cells']
 relations=C.load(C.R1/'prepared/all-local-relations.json');rmap={(r['case'],r['mode'],r['id'],r['owner']):r for r in relations}
 require(len(rmap)==len(relations),'unique-relation-domain')
 owner_assoc=C.load(C.P/'actual-typed-owner-associations.json');chain=C.load(C.P/'captured-owner-proof-chain.json');output_assoc=C.load(C.P/'actual-output-axis-associations.json')
 oa={(x['case'],x['mode'],x['id'],x['owner']):x for x in owner_assoc};ch={(x['case'],x['mode'],x['id'],x['owner']):x for x in chain};outa={(x['case'],x['mode'],x['id']):x for x in output_assoc}
 science_interface=C.load(C.S/'interface-independent.json');sf={(x['case'],x['mode'],x['id']):x for x in science_interface['actual_folds']}
 science_math=C.load(C.S/'math-physical-reconstruction.json');science_source=C.load(C.S/'source-runtime-independent.json')
 sl={(x['case'],x['mode'],x['id'],x['owner'],x['route']):x for x in science_source['lowering_layouts']}
 accepted={(x['case'],x['Classes']):x for x in sr['per_case_mode_verdicts']}
 require([(x['case'],x['mode']) for x in inputs]==[(x['case'],x['mode']) for x in native]==[(c,m) for c in C.CASES for m in C.MODES],'six-independent-cells')
 faith={e:{} for e in ('historical','current')};aggregate=[];chains=[];summaries=[];allfold=[];nons=0
 for cell,ncell in zip(inputs,native):
  case,mode=cell['case'],cell['mode'];approved=accepted[case,mode];require(approved['finite_physical_value_mathematics']==approved['finite_named_consumer_interface']=='APPROVE','scientific-cell-approved')
  worlds={e:C.load(C.CAP/f'prepared/captured/{case}/{e}-{mode}.json') for e in ('historical','current')};cell_input(cell,worlds)
  tables={e:G.owner_table(worlds[e],e) for e in worlds};G.physical(worlds,tables)
  for era in worlds:
   defs,refs=whole_definitions(worlds[era],tables[era],era);require(defs==cell['definitions'][era],'complete-whole-definitions')
   for key,ref in refs:faith[era].setdefault(key,[]).append(ref)
   for d in worlds[era]['payload']['diagrams']:
    p=d['payload'];factors=[p['numerator_prefactor'],*[v['numerator'] for v in p['vertices']],*[e[1]['numerator'] for e in p['edges']]]
    aggregate.append({'historical':normal('*'.join('('+s+')' for s in factors),era=='historical'),'current':normal(p['numerator'],era=='historical'),'label':f'{case}/{mode}/{era}/{d["id"]}/unexpanded-local-product','case':case,'mode':mode,'route':era,'id':d['id'],'factor_owners':['numerator_prefactor',*[f'vertex:{i}' for i in range(len(p['vertices']))],*[f'edge:{i}' for i in range(len(p['edges']))]],'projector':'separate captured owner; excluded by actual local numerator assembly','overall_factor':d['overall_factor'],'normalization':'No distribute/expand; each original sum remains a factor'})
  require(cell['definitions']['historical']['mapped_decoded_ports_and_factors']==cell['definitions']['current']['mapped_decoded_ports_and_factors'],'mapped-definition-equality')
  require(len(cell['graphs'])==len(ncell['graphs'])==len(worlds['current']['payload']['diagrams']),'complete-graph-domain')
  counts=collections.Counter()
  for g,ng,hd,cd in zip(cell['graphs'],ncell['graphs'],worlds['historical']['payload']['diagrams'],worlds['current']['payload']['diagrams']):
   gid=g['id'];gkey=(case,mode,gid)
   require(gid==ng['id']==hd['id']==cd['id'] and g['name']==ng['name']==cd['name'] and g['full_key']==hd['key']==cd['key'],'exact-ID-key-name')
   require(g['current_payload']==cd['payload'] and g['mapped_historical_payload']==mapped(hd['payload']),'exact-proof-payload')
   by={r['owner']:r for r in g['owners']};nb={r['owner']:r for r in ng['owners']};require(len(by)==len(g['owners']) and len(nb)==len(ng['owners']) and by.keys()==nb.keys()==tables['current'][gid].keys(),'exact-owner-domain')
   for label,rec in by.items():
    key=(case,mode,gid,label);h,c=tables['historical'][gid][label],tables['current'][gid][label];nr=nb[label];relation=rmap[key]
    require(relation['raw_decoded']=={'historical':h,'current':c} and relation['full_key']==g['full_key'],'relation-to-capture')
    for era in worlds:require(relation['capture_inputs'][era]==C.sha(C.read(C.CAP/f'prepared/captured/{case}/{era}-{mode}.json')),'relation-world-sha')
    require(relation['historical']==normal(h['decoded'],True) and relation['current']==normal(c['decoded']),'complete-expression-relation')
    lt,rt=B.parse(h['decoded'],True),B.parse(c['decoded']);require(B.definitions(lt)==B.definitions(rt),'owner-full-definition')
    require((lt==rt)==relation['literal'] and B.cyclic(lt)==B.cyclic(rt),'finite-accepted-relation');nons+=not relation['literal']
    changed=owner_input(rec,nr,h,c,owners(cd['payload'])[label]);counts['changed' if changed else 'unchanged']+=1;counts['owners']+=1
    ar=oa[key];require(resolve(ar['actual_native'])==nr,'exact-owner-result-pointer');require(ar['full_key']==g['full_key'] and ar['operand_map']==nr['transport'] and ar['historical_output_map']==nr['historical']['canonical_to_logical'] and ar['current_output_map']==nr['current']['canonical_to_logical'],'candidate-owner-maps')
    cr=ch[key];require(relations[cr['relation_index']]==relation and cr['full_key']==g['full_key'] and cr['relation_literal']==relation['literal'],'candidate-owner-chain')
    if changed:
     for era in worlds:
      lr=sl[case,mode,gid,label,era];require(lr['canonical_to_logical']==nr[era]['canonical_to_logical'] and lr['logical_to_canonical']==[nr[era]['canonical_to_logical'].index(i) for i in range(len(rec[era]))],'scientific-layout-applies-to-owner')
    chains.append({'case':case,'mode':mode,'id':gid,'owner':label,'full_key':g['full_key'],'current_serialized_expression_sha256':C.sha(encoded(c['decoded'])),'current_serialized_interface_sha256':C.sha(encoded(owners(cd['payload'])[label][1])) if rec['explicit'] else None,'relation_literal':relation['literal'],'native_result_reference':ar['actual_native'],'operand_map':nr['transport'],'historical_output_map':nr['historical']['canonical_to_logical'],'current_output_map':nr['current']['canonical_to_logical']})
   selected=[(r['owner'],label) for r in g['owners'] if r['explicit'] for label in by if ':' in label and r['owner']!=label]
   require(selected==[(x['left'],x['right']) for x in ng['plans']],'exact-selected-pair-domain')
   for p in ng['plans']:
    a,b=by[p['left']],by[p['right']]
    for era in worlds:plan(p[era+'_plan'],p['physical_plan'],a[era],b[era])
   counts['selected_pairs']+=len(selected)
   previous=[];seen=[];steps=[];scientific_steps=[]
   for si,step in enumerate(ng['actual_fold']):
    label=step['owner'];require(':' in label and label in by and label not in seen,'complete-fold-owner');seen.append(label)
    for era in worlds:plan(step[era+'_plan'],step['physical_plan'],previous,by[label][era])
    pairs,surv=product(previous,by[label]['current']);output=step['output'];layout(output);require(sorted(output['ports'])==surv and step['same_value'] is True,'full-fold-output')
    op=[output['ports'].index(s) for s in surv]
    steps.append({'step':si,'owner':label,'contraction_pairs':pairs,'neutral_to_actual_output':op,'canonical_to_actual_output':output['canonical_to_logical'],'actual_output_ports':output['ports']})
    scientific_steps.append({'owner':label,'kind':step['physical_plan']['kind'],'pairs':len(pairs),'neutral_to_actual_output':op,'canonical_to_logical':output['canonical_to_logical']})
    previous=output['ports'];counts['fold_endpoint_pairs']+=len(pairs)
   require(set(seen)=={k for k in by if ':' in k},'full-fold-factor-set');require(scientific_steps==sf[gkey]['steps'],'accepted-scientific-folds-exact')
   oa_graph=outa[gkey];require(steps==oa_graph['fold_outputs'] and resolve(oa_graph['native_graph'])==ng,'candidate-fold-output-maps')
   require(sorted(previous)==sorted(by['numerator']['current']) and len(previous)==8 and by['numerator_prefactor']['current']==[] and len(by['projector']['current'])==4,'rank8-rank4')
   pp,full=product(previous,by['projector']['current']);require(pp==oa_graph['projector_pairs'] and full==oa_graph['full_projected_neutral_output'] and len(full)==4 and all(port(s)[0][0]=='coad' for s in full),'four-coad-output')
   layout(ng['public_numerator_of']);require(ng['public_numerator_of']['ports']==previous,'public-output')
   for k in ('whole_graph_import','whole_graph_roundtrip','fragment_numerator_multiply_many_same_value','public_numerator_of_same_value'):require(ng[k] is True,'native-graph-'+k)
   counts['fold_steps']+=len(steps);allfold.append({'case':case,'mode':mode,'id':gid,'steps':steps,'projector_pairs':pp,'full_output':full})
  sciencecell=next(x for x in science_interface['per_cell'] if (x['case'],x['mode'])==(case,mode));require(counts['owners']==sciencecell['owners'] and counts['fold_steps']==sciencecell['fold_steps'],'scientific-domain-complete')
  factor2=[d['id'] for d in worlds['current']['payload']['diagrams'] if 'NumeratorIndependentSymmetryGrouping(2)' in d['overall_factor']]
  require(len(factor2)==(4 if case==C.CASES[0] else 0),'factor2-four-of-five')
  summaries.append({'case':case,'mode':mode,'graphs':len(cell['graphs']),**counts,'factor2_ids':factor2})
 require(len(chains)==600 and nons==22 and set(oa)==set(ch)==set(rmap),'entire-proof-owner-domain')
 require(aggregate==C.load(C.R1/'prepared/all-aggregate-relations.json'),'all96-unexpanded-products')
 # Independently reconstruct entire native partitions from their actual expression owners.
 hist=[{'historical':a,'current':b,'obligation':'parse/render-faithfulness','owners':refs} for (a,b),refs in faith['historical'].items() if a!=b]
 curr=[{'historical':a,'current':b,'obligation':'parse/render-faithfulness','owners':refs} for (a,b),refs in faith['current'].items() if a!=b]
 cross=[{'historical':r['historical'],'current':r['current'],'obligation':'cross-owner','owners':[{k:r[k] for k in ('case','mode','id','owner')}]} for r in relations if not r['literal']]+[{**a,'obligation':'unexpanded-aggregate'} for a in aggregate if a['historical']!=a['current']]
 math=[]
 for batch,assocs in [('faithfulness-historical',hist),('current-faithfulness-and-equivalence',curr+cross)]:
  dedup={}
  for r in assocs:dedup.setdefault((r['historical'],r['current']),[]).append(r)
  records=[{'historical':a,'current':b,'obligations':v} for (a,b),v in dedup.items()]
  require(records==C.load(C.R1/'prepared/math'/batch/'obligations.json'),'complete-native-partition-owners')
  portable=C.load(C.R1/'prepared/math'/batch/'portable-pairs.json');require(len(portable)==len(records),'portable-partition-domain')
  defs={};order=[];originaldefs={}
  for r in records:
   for side in ('historical','current'):B.definitions(B.parse(r[side]),originaldefs)
  for i,(r,pair) in enumerate(zip(records,portable)):
   z=f'next_qcd_proof::{{}}::{"fresh_witness_" if batch=="faithfulness-historical" else "merged_fresh_witness_"}{i}'
   require('next_qcd_proof::'+z.rsplit('::',1)[-1] not in originaldefs,'fresh-formal-witness')
   require(pair=={s:'('+r[s]+')+'+z for s in ('historical','current')},'exact-formal-witness-pair')
   for side in ('historical','current'):
    ds=B.definitions(B.parse(pair[side]));B.definitions(B.parse(pair[side]),defs)
    for token in ds.values():
     if token not in order:order.append(token)
  declaration=C.load(C.R1/'prepared/math'/batch/'definition-order.json');require(defs==declaration['definitions'] and order==declaration['first_use_order'] and declaration['pair_side_order']==['historical','current'],'full-native-definition-side-order')
  for kernel in ('c78','8eb'):
   data=C.load(C.R3/f'execution/records/{batch}--{kernel}/results.json');require(len(data)==len(records) and all(set(r)=={'parsed_equal','canonical_equal','sign_rejected','coefficient_rejected','power_rejected'} and all(v is True for v in r.values()) for r in data),'actual-native-partition-results')
  math.append({'batch':batch,'unique_pairs':len(records),'obligations':len(assocs),'definitions':len(defs)})
 inverse=C.load(C.N/'math-member-to-native-row.json')
 for r in inverse:
  result,ob=resolve(r['result']),resolve(r['obligation']);resolve(r['definition_order']);require(result==r['native_predicates'] and all(v is True for v in result.values()),'actual-math-row-binding')
  if r['member'] is None:require(ob['obligation']=='unexpanded-aggregate' and {k:ob[k] for k in ('case','mode','route','id','historical','current')} in [{k:a[k] for k in ('case','mode','route','id','historical','current')} for a in aggregate],'aggregate-member')
  else:require(r['member'] in ob['owners'],'exact-math-owner-member')
 # Full independent audits are dependencies of the accepted scientific result, never replaced by counts.
 for f in ['associations-independent.json','execution-independent.json','actual-negative-inputs-independent.json','parser-source-independent.json','supplemental-final-source-rehash.json','fresh-heldouts.json','failures.json']:
  C.load(C.S/f)
 C.save('proof-owner-join-independent.json',chains);C.save('fold-output-join-independent.json',allfold)
 C.save('proof-review-independent.json',{'scientific_review_sha256':C.sha(C.read(C.S/'review.json')),'root_additive_join_sha256':C.sha(C.read(C.F/'root-next-qcd-scientific-interface-acceptance-20261031.json')),'per_cell':summaries,'all_local_relations':600,'literal_relations':578,'nonliteral_even_cycles':nons,'unexpanded_products':len(aggregate),'native_math_partitions':math,'inverse_references':len(inverse),'scientific_lowering_layouts_bound_to_actual_owners':len(sl),'restriction':'Accepted finite physical mathematics and named interfaces attached to all exact pristine digest owners. No new native/library work, general validator or broader runtime theorem.'})
 C.flush('proof-consumed.json');print('PASS complete scientific join:600 exact owners,96 full factor products,166/118 native partitions and4268 references, both operand/output axes; actual integer plans/folds checked.')
if __name__=='__main__':main()
