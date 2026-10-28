from pathlib import Path
import json
O=Path('target/wp3-units19-21/frozen-row-diagnosis');env={};s=Path('target/wp3-drawing-transplant/frame-digest-boundary.py').read_text();exec(s[:s.index('expected=')],env);heads=env['spellings']
def modern(s):
 for a,b in heads.items():s=s.replace(b,a)
 return s
pairs={};records=[];bijections=[]
for name in ['ym_ggg_zero_color_L3','xs_qcdqed_addx_canonical_scalar_lr_L3']:
 old=json.loads((O/'baseline-r2/captures'/name/'raw.json').read_text());new=json.loads((O/'current-metadata-r3/captures'/name/'raw.json').read_text())
 key=lambda d:json.dumps(d['key'],sort_keys=True,separators=(',',':'))
 aa={key(d):(i,d) for i,d in enumerate(old['diagrams'])};bb={key(d):(i,d) for i,d in enumerate(new['diagrams'])};assert len(aa)==len(old['diagrams']) and len(bb)==len(new['diagrams']) and aa.keys()==bb.keys()
 for k in sorted(aa):
  oi,a=aa[k];ni,b=bb[k];assert a['id']==b['id'];p=a['payload'];q=b['payload'];assert a['basis']==b['basis']
  # Local owners are joined by their retained graph order only after every owner incidence field matches.
  assert len(p['vertices'])==len(q['vertices']) and len(p['edges'])==len(q['edges'])
  bijections.append({'case':name,'baseline_index':oi,'current_index':ni,'id':a['id'],'full_key':a['key']})
  values=[('projector',p['projector'],q['projector']),('prefactor',p['numerator_prefactor'],q['numerator_prefactor'])]
  for j,(x,y) in enumerate(zip(p['vertices'],q['vertices'])):
   assert {k:v for k,v in x.items() if k not in ['numerator','num_interface']}=={k:v for k,v in y.items() if k not in ['numerator','num_interface']}
   values.append((f'vertex:{j}',x['numerator'],y['numerator']))
  for j,(x,y) in enumerate(zip(p['edges'],q['edges'])):
   assert x[0]==y[0]
   assert {k:v for k,v in x[1].items() if k not in ['numerator','num_interface']}=={k:v for k,v in y[1].items() if k not in ['numerator','num_interface']}
   values.append((f'edge:{j}',x[1]['numerator'],y[1]['numerator']))
  for owner,x,y in values:
   x=modern(x)
   if x==y:proof='literal_equal_under_verified_head_map'
   else:
    pair=(x,y)
    if pair not in pairs:pairs[pair]=len(pairs)
    proof=pairs[pair]
   records.append({'case':name,'id':a['id'],'owner':owner,'proof':proof})
unique=[None]*len(pairs)
for (x,y),i in pairs.items():unique[i]={'baseline':x,'current':y}
(O/'local-pairs.json').write_text(json.dumps(unique)+'\n');(O/'local-factor-inventory.json').write_text(json.dumps(records,indent=2)+'\n');(O/'full-key-bijection.json').write_text(json.dumps(bijections,indent=2)+'\n');print(len(bijections),'exact full-key joins;',len(records),'local factors;',len(unique),'unique nonliteral local pairs')
