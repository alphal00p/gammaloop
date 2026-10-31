import sys,json,pathlib,hashlib,copy,collections,re,functools,itertools
W=pathlib.Path('/shared/lcnbr/dev/gammaloop/generation-oracle-numerator-families');D=pathlib.Path(__file__).resolve().parent;R=W/'target/raw-next-qcd-proof-preparation-20261031-r3';P=W/'target/raw-next-qcd-proof-preparation-20261031-r1';E=W/'target/raw-next-qcd-proof-native-execution-20261031-r1'
sys.path.insert(0,str(R));import prepare as A;B=A.B

def j(p):return json.loads(pathlib.Path(p).read_bytes(),object_pairs_hook=A.unique)
def save(n,x):(D/n).write_text(json.dumps(x,indent=2)+'\n')
replayed={};A.save=lambda n,x:replayed.update({n:copy.deepcopy(x)})
faith={r:{}for r in A.ROUTES};allrelations=[];allproducts=[];cells=[];worldsall={};scopes=[]
for case in A.CASES:
 for mode in A.MODES:
  worlds={r:A.load(A.C/f'prepared/captured/{case}/{r}-{mode}.json')for r in A.ROUTES};worldsall[case,mode]=worlds
  cell,matrix,relations,products,ff,ll=A.analyze(case,mode,worlds)
  for route in A.ROUTES:
   for k,refs in ff[route].items():faith[route].setdefault(k,[]).extend(refs)
  cells.append(cell);allrelations+=relations;allproducts+=products
  scopes.append({k:matrix[k]for k in ['case','mode','graphs','owners','pure_literal','nonliteral','affected_port_orders','declared_interfaces','first_raw_boundary','first_decoded_boundary','first_public_interface_boundary','source_factor_trees']})
for name,data in replayed.items():assert data==j(P/name),name
assert allrelations==j(P/'prepared/all-local-relations.json') and allproducts==j(P/'prepared/all-aggregate-relations.json')
# Whole-tree finite congruence, replaying exact recorded paths without cyclic canonicalization.
traces=j(P/'prepared/gggg-even-cycle-witnesses.json');td={(t['case'],t['mode'],t['id'],t['owner']):t for t in traces};rewrites=0

def at(tree,path):
 if not path:return tree
 if path[0]=='args':return at(tree[2][path[1]],path[2:])
 return at(tree[path[0]],path[1:])
def replace(tree,path,replacement):
 if not path:return replacement
 q=list(tree)
 if path[0]=='args':
  args=list(tree[2]);args[path[1]]=replace(args[path[1]],path[2:],replacement);q[2]=tuple(args)
 else:q[path[0]]=replace(q[path[0]],path[1:],replacement)
 return tuple(q)
for r in allrelations:
 left=B.parse(r['historical']);right=B.parse(r['current']);assert B.definitions(left)==B.definitions(right)
 assert (left==right)==r['literal']
 if left==right:continue
 t=td[r['case'],r['mode'],r['id'],r['owner']]
 for step in t['rewrites']:
  x=at(left,step['path']);assert x[0]=='call' and x[1]==B.F==step['definition'] and len(x[2])==3
  s=step['left_rotation'];pi=list(range(3))[s:]+list(range(3))[:s];inv=sum(pi[i]>pi[k]for i in range(3)for k in range(i+1,3));assert s in (1,2) and inv%2==0
  assert step['parity']==step['coefficient']==1 and step['transpositions']==2
  assert step['source_arguments']==[B.render(a)for a in x[2]]
  y=('call',B.F,tuple(x[2][i]for i in pi));assert step['target_arguments']==[B.render(a)for a in y[2]]
  left=replace(left,step['path'],y);rewrites+=1
 assert left==right
assert len(td)==22 and rewrites==144
# Independently reconstruct dedup domains/inverse obligations from captured strings and relations.
hist=[{'historical':a,'current':b,'obligation':'parse/render-faithfulness','owners':refs}for (a,b),refs in faith['historical'].items()if a!=b]
current=[{'historical':a,'current':b,'obligation':'parse/render-faithfulness','owners':refs}for (a,b),refs in faith['current'].items()if a!=b]
cross=[{'historical':r['historical'],'current':r['current'],'obligation':'cross-owner','owners':[{k:r[k]for k in ('case','mode','id','owner')}]}for r in allrelations if not r['literal']]+[{**a,'obligation':'unexpanded-aggregate'}for a in allproducts if a['historical']!=a['current']]
math=[]
for name,associations in [('faithfulness-historical',hist),('current-faithfulness-and-equivalence',current+cross)]:
 dedup={}
 for a in associations:dedup.setdefault((a['historical'],a['current']),[]).append(a)
 records=[{'historical':a,'current':b,'obligations':v}for (a,b),v in dedup.items()];assert records==j(P/'prepared/math'/name/'obligations.json')
 portable=j(P/'prepared/math'/name/'portable-pairs.json');defs={};order=[];original_defs={};z_names=set()
 for record in records:
  for side in ('historical','current'):B.definitions(B.parse(record[side]),original_defs)
 for i,(record,pair) in enumerate(zip(records,portable)):
  z=f'next_qcd_proof::{{}}::{"fresh_witness_"if name=="faithfulness-historical"else"merged_fresh_witness_"}{i}'
  z_name='next_qcd_proof::'+z.rsplit('::',1)[1];assert z_name not in original_defs and z_name not in z_names;z_names.add(z_name)
  assert pair=={side:'('+record[side]+')+'+z for side in ('historical','current')}
  for side in ('historical','current'):
   B.definitions(B.parse(pair[side]),defs)
   for token in B.definitions(B.parse(pair[side])).values():
    if token not in order:order.append(token)
 recorded=j(P/'prepared/math'/name/'definition-order.json');assert defs==recorded['definitions'] and order==recorded['first_use_order'] and recorded['pair_side_order']==['historical','current']
 for kernel in ('c78','8eb'):
  out=j(R/'execution/records'/(name+'--'+kernel)/'results.json');assert len(out)==len(records) and all(set(x)=={'parsed_equal','canonical_equal','sign_rejected','coefficient_rejected','power_rejected'}and all(v is True for v in x.values()) for x in out)
 math.append({'batch':name,'unique_pairs':len(records),'obligation_groups':len(associations),'definitions':len(defs),'member_references_per_kernel':sum(len(a.get('owners',[None]))for a in associations)})
assert [m['unique_pairs']for m in math]==[166,118]
# Independent finite route dependency and mass ownership, checked against every stored witness.
local=[];calls=0

def nodes(t):
 yield t
 for c in (t[2]if t[0]=='call'else t[1:]if t[0]not in ('integer','symbol')else()):yield from nodes(c)
for case in A.CASES:
 for mode in A.MODES:
  for route in A.ROUTES:
   w=worldsall[case,mode][route];model=w['live']['model'];dkeys=[]
   for d in w['payload']['diagrams']:
    p=d['payload'];dkeys.append(json.dumps(d['key'],sort_keys=True));raw=next(g for g in w['raw']['diagrams']if g['id']==d['id'])
    for o in raw['owners']:
     name=o['owner']
     if ':'not in name:continue
     kind,idx=name.split(':');idx=int(idx);allowed={idx}if kind=='edge'else{i for i,(ends,e)in enumerate(p['edges'])if idx in ends.values()};tree=B.parse(o['tensor']['decoded'],route=='historical');actual=[]
     for n in nodes(tree):
      if n[0]=='call'and A.base(n[1]) in ('gammalooprs::Q','gammalooprs::K','gammalooprs::P'):
       assert A.base(n[1])=='gammalooprs::Q';assert n[2][0][0]=='integer';e=int(n[2][0][1]);assert e in allowed
       actual.append({'source_call':B.render(n),'edge':e,'signed_routing':d['basis']['edge_signatures'][str(e)]});calls+=1
     item={'case':case,'mode':mode,'route':route,'id':d['id'],'full_key':d['key'],'owner':name,'allowed_momentum_edges':sorted(allowed),'actual_momentum_calls':actual,'no_nonlocal_momentum_dependency':True}
     if kind=='edge':
      ed=p['edges'][idx][1];particle=model['particles'][ed['particle']];m=[n[1]for n in nodes(tree)if n[0]=='symbol'and n[1].startswith('UFO::')and not A.base(n[1]).split('::')[-1].startswith('GC_')];assert all(A.base(x).split('::')[-1]==particle['mass']for x in m)
      item.update(particle_id=ed['particle'],particle=particle,edge_record={k:v for k,v in ed.items()if k not in ('numerator','num_interface')},mass_symbols=m)
     local.append(item)
   assert len(dkeys)==len(set(dkeys))
assert local==j(P/'prepared/finite-locality-and-mass-inputs.json')['records'] and len(local)==912 and calls==1096
save('math-physical-reconstruction.json',{'scope':scopes,'recomputed_cell_files':len(replayed),'full_relations':len(allrelations),'literal':578,'even_cycles':22,'even_rewrites':rewrites,'unexpanded_products':len(allproducts),'math_batches':math,'locality_records':len(local),'momentum_calls':calls,'method':'Fresh reconstruction from twelve captured worlds; exact full pure grammar plus independent path/parity congruence and finite locality traversal. No expansion or numerical substitution. Native results read, not rerun.'})
print('math and physical reconstruction passed',flush=True)
