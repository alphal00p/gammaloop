import json,pathlib,re,functools,collections,itertools
W=pathlib.Path('/shared/lcnbr/dev/gammaloop/generation-oracle-numerator-families');D=pathlib.Path(__file__).resolve().parent;P=W/'target/raw-next-qcd-proof-preparation-20261031-r1';R=W/'target/raw-next-qcd-proof-preparation-20261031-r3'
def j(p):return json.loads(pathlib.Path(p).read_bytes())
def save(n,x):(D/n).write_text(json.dumps(x,indent=2)+'\n')
@functools.cache
def port(s):
 original=s;dual=False
 if s.startswith('spenso::{spenso::lower}::dind('):assert s.endswith(')');dual=True;s=s[len('spenso::{spenso::lower}::dind('):-1]
 m=re.fullmatch(r'spenso::\{spenso::representation,spenso::(self_dual|dualizable)\}::(bis|mink|coad|cof)\(([0-9]+),gammalooprs::\{spenso::index\}::(hedge|edge)\(([0-9]+),([0-9]+)\)\)',s);assert m,original
 tag,rep,dim,own,idx,lane=m.groups();assert int(dim)=={'bis':4,'mink':4,'coad':8,'cof':3}[rep] and int(lane)==1
 assert tag==('dualizable'if rep=='cof'else'self_dual') and (not dual or rep=='cof')
 return (rep,int(dim),own,int(idx),int(lane)),dual

def match(a,b):
 x,xd=port(a);y,yd=port(b);return x==y and (x[0]!='cof'or xd!=yd)
def product(l,r):
 pairs=[(i,k)for i,a in enumerate(l)for k,b in enumerate(r)if match(a,b)];assert len({i for i,k in pairs})==len(pairs)==len({k for i,k in pairs})
 survivors=[('left',i,a)for i,a in enumerate(l)if i not in {i for i,k in pairs}]+[('right',k,b)for k,b in enumerate(r)if k not in {k for i,k in pairs}];out=sorted(s for _,_,s in survivors);assert len(out)==len(set(out))
 return {'pairs':[list(p)for p in pairs],'neutral_output':out,'survivors':[{'operand':[side,i],'output_axis':out.index(s),'port':s}for side,i,s in survivors]}
def square(a,b,lp,rp,op):
 assert sorted([lp[i],rp[k]]for i,k in a['pairs'])==sorted(b['pairs']);assert sorted(op)==list(range(len(op)))
 q={tuple(s['operand']):s['output_axis']for s in b['survivors']}
 for s in a['survivors']:
  side,i=s['operand'];assert q[side,(lp if side=='left'else rp)[i]]==op[s['output_axis']]
def layout(n):
 l=n['ports'];c=n['canonical_ports'];inv=n['canonical_to_logical'];assert sorted(inv)==list(range(len(l))) and c==[l[i]for i in inv]
 forward=[inv.index(i)for i in range(len(inv))];assert l==[c[i]for i in forward]
 return forward,inv

def debug_plan(s,l,r):
 pred=product(l,r);kind='Scalar'if not l or not r else'Contract'if pred['pairs']else'Outer'
 if kind!='Contract':assert s==kind;return {'kind':kind}
 pattern=r'PortPair \{ left: ([0-9]+), right: ([0-9]+) \}'
 assert s.startswith('Contract([')and s.endswith('])');body=s[10:-2];parts=re.findall(pattern,body);assert body==', '.join(f'PortPair {{ left: {i}, right: {k} }}'for i,k in parts)
 pairs=[[int(i),int(k)]for i,k in parts];assert pairs==pred['pairs']
 return {'kind':'Contract','names':sorted([[l[i],r[k]]for i,k in pairs])}
inputs=j(P/'prepared/interface-inputs.json')['cells'];native=j(R/'execution/records/all-six-cells--typed8eb/results.json')['cells'];cons={(x['case'],x['mode'],x['id']):x for x in j(P/'prepared/consumer-obligations.json')};layouts={(x['id'],json.dumps(x['full_key'],sort_keys=True)):x for x in j(P/'prepared/all-layout-cells.json')}
count=collections.Counter();percell=[];prefixcert=[];builtins=[];foldcert=[]
for cell,ncell in zip(inputs,native):
 assert (cell['case'],cell['mode'])==(ncell['case'],ncell['mode']);localcount=collections.Counter()
 for g,ng in zip(cell['graphs'],ncell['graphs']):
  assert (g['id'],g['name'])==(ng['id'],ng['name']);by={r['owner']:r for r in g['owners']};nb={r['owner']:r for r in ng['owners']};con=cons[cell['case'],cell['mode'],g['id']];locals=[r['owner']for r in g['owners']if ':'in r['owner']]
  assert list(by)==list(nb);changed=0
  for name,r in by.items():
   h,c=r['historical'],r['current'];pi=r['permutation'];nr=nb[name];assert len(set(h))==len(h)and len(set(c))==len(c) and sorted(pi)==list(range(len(h)))and h==[c[i]for i in pi]
   for s in h+c:port(s)
   hp,hi=layout(nr['historical']);cp,ci=layout(nr['current'])
   assert nr['historical']['ports']==h and nr['current']['ports']==c and nr['transport']==pi
   assert nr['historical']['canonical_ports']==nr['current']['canonical_ports'];transport=[cp[pi[hi[k]]]for k in range(len(hi))];assert transport==list(range(len(hi)))
   for flag in ['is_metric','is_composite']:assert nr['historical'][flag]==nr['current'][flag]==r[flag]
   affected=pi!=list(range(len(pi)));assert nr['affected']is affected and nr['same_value_without_transport']is (not affected)
   for flag in ['transport_same_value','declared_admission','owner_roundtrip']:assert nr[flag]is True
   if affected:
    changed+=1
    for side in ['current','historical']:
     lr=nr[side+'_lowering'];assert lr['dangling']==sorted(r[side]) and lr['execution']=='lowering only, no component execution'
    controls=nr['controls'];assert controls['wrong_permutation_reject']is True
    inv=[pi.index(i)for i in range(len(pi))];assert controls['inverse_is_same_permutation']==(pi==inv)and controls['inverse_accepted']==(pi==inv)
    for k in ['omitted','duplicate','wrong_owner']:assert controls[k]['accepted']is False
    noop=c[0].replace('::bis(4,','::bis(5,',1)==c[0];assert controls['wrong_type']['accepted']is noop
    builtins.append({'case':cell['case'],'mode':cell['mode'],'id':g['id'],'owner':name,'wrong_type_is_noop':noop,'rejection_credit':False if noop else 'diagnostic only; external used-domain negatives separately authenticated'})
   else:assert 'current_lowering'not in nr and 'historical_lowering'not in nr
  localcount['owners']+=len(by);localcount['changed_lowerings']+=changed;localcount['unchanged_source_arguments']+=len(by)-changed
  # Complete selected ordered pair domain, both integer plans and full logical squares.
  selected=[(a['owner'],b)for a in g['owners']if a['explicit']for b in locals if a['owner']!=b]
  assert selected==[(x['left'],x['right'])for x in con['owner_pair_inputs']]==[(x['left'],x['right'])for x in ng['plans']]
  for (left,right),q,n in zip(selected,con['owner_pair_inputs'],ng['plans']):
   a,b=by[left],by[right];h=product(a['historical'],b['historical']);c=product(a['current'],b['current']);assert h==q['historical']and c==q['current'];square(h,c,a['permutation'],b['permutation'],q['output_map'])
   cp=debug_plan(n['current_plan'],a['current'],b['current']);hp=debug_plan(n['historical_plan'],a['historical'],b['historical']);assert cp==hp==n['physical_plan']
   if hp['kind']=='Contract':
    oldpairs=h['pairs'];stale=sorted([[a['current'][i],b['current'][k]]for i,k in oldpairs]);assert n['stale_positional_plan_changes_wiring']==(stale!=hp['names'])
   else:assert n['stale_positional_plan_changes_wiring']is False
  localcount['selected_pairs']+=len(selected)
  for q in con['prefactor_projector_maps']:
   a,b=by[q['left']],by[q['right']];h=product(a['historical'],b['historical']);c=product(a['current'],b['current']);assert h==q['historical']and c==q['current'];square(h,c,a['permutation'],b['permutation'],q['output_map'])
  assert len(by['numerator']['current'])==8 and len(by['projector']['current'])==4 and by['numerator_prefactor']['current']==[]
  projection=product(by['numerator']['current'],by['projector']['current']);assert len(projection['pairs'])==4 and len(projection['neutral_output'])==4 and all(port(s)[0][0]=='coad'for s in projection['neutral_output'])
  localcount['factor_products']+=2
  # All possible visited-local-factor sets: no native subset/component execution is inferred.
  groups=collections.defaultdict(list)
  for name in locals:
   for axis,s in enumerate(by[name]['current']):groups[port(s)[0]].append((name,axis,s))
  assert all(len(v)in(1,2)for v in groups.values());assert all(len(v)==1 or (v[0][0]!=v[1][0]and match(v[0][2],v[1][2]))for v in groups.values())
  frontier=con['fold_frontier'];assert frontier['factors']==locals
  gp=sorted(sorted([[v[0][0],v[0][1]],[v[1][0],v[1][1]]])for v in groups.values()if len(v)==2)
  assert gp==sorted(sorted([q['left'],q['right']])for q in frontier['global_pairs'])
  for mask in range(1<<len(locals)):
   visited={n for k,n in enumerate(locals)if mask>>k&1};expected=sorted([s for v in groups.values() for n,i,s in v if n in visited and sum(x[0]in visited for x in v)==1])
   for route in ('historical','current'):
    live=[]
    for n in locals:
     if n in visited:live=product(live,by[n][route])['neutral_output']
    assert live==expected
  localcount['visited_sets']+=1<<len(locals);prefixcert.append({'case':cell['case'],'mode':cell['mode'],'id':g['id'],'local_factors':len(locals),'visited_sets':1<<len(locals),'meaning':'all subset frontiers via sequential logical reduction and independent whole-incidence count, each route separately; not all native subgraphs or component arithmetic'})
  prev=[];visited=[];foldrows=[]
  for step in ng['actual_fold']:
   n=step['owner'];assert n in locals and n not in visited;visited.append(n)
   cp=debug_plan(step['current_plan'],prev,by[n]['current']);hp=debug_plan(step['historical_plan'],prev,by[n]['historical']);assert cp==hp==step['physical_plan']
   pred=product(prev,by[n]['current']);out=step['output'];layout(out);assert sorted(out['ports'])==pred['neutral_output']and step['same_value']is True
   mapout=[out['ports'].index(s)for s in pred['neutral_output']];assert sorted(mapout)==list(range(len(mapout)))
   for survivor in pred['survivors']:assert out['ports'][mapout[survivor['output_axis']]]==survivor['port']
   foldrows.append({'owner':n,'kind':cp['kind'],'pairs':len(pred['pairs']),'neutral_to_actual_output':mapout,'canonical_to_logical':out['canonical_to_logical']});prev=out['ports'];localcount['fold_endpoint_pairs']+=len(pred['pairs'])
  assert set(visited)==set(locals)and sorted(prev)==sorted(by['numerator']['current']);layout(ng['public_numerator_of']);assert ng['public_numerator_of']['ports']==prev
  for flag in ['whole_graph_import','whole_graph_roundtrip','fragment_numerator_multiply_many_same_value','public_numerator_of_same_value']:assert ng[flag]is True
  localcount['fold_steps']+=len(visited);localcount['graphs']+=1;foldcert.append({'case':cell['case'],'mode':cell['mode'],'id':g['id'],'steps':foldrows})
 count.update(localcount);percell.append({'case':cell['case'],'mode':cell['mode'],**localcount});print(cell['case'],cell['mode'],dict(localcount),flush=True)
assert count['graphs']==48 and count['owners']==600 and count['selected_pairs']==1890 and count['visited_sets']==41088 and count['fold_steps']==456 and count['fold_endpoint_pairs']==720
assert count['changed_lowerings']==210 and count['unchanged_source_arguments']==390 and sum(x['wrong_type_is_noop']for x in builtins)==26
save('interface-independent.json',{'counts':count,'per_cell':percell,'prefixes':prefixcert,'actual_folds':foldcert,'builtin_diagnostics':builtins,'canonical_argument':'p maps logical to canonical, native diagnostic is p inverse. For historical logical->current logical pi, canonical transport is p_current o pi o inverse(p_historical). All 600 canonical transports are identities on the matched canonical field arrays; logical pi is often nonidentity. Both operands and surviving outputs are explicitly checked.','prefix_qualification':'41088 includes empty/full visited-factor subsets across 48 graph occurrences; the preparer supplies a rule, this independent review enumerates its domain. 456 actual native fold steps only. No arbitrary component evaluation or new historical kernel consumer run.'})
print('interface independent reconstruction passed',dict(count),flush=True)
