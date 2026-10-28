from pathlib import Path
import json,collections,hashlib
O=Path('target/wp3-units19-21/frozen-row-diagnosis'); A=O/'baseline-r2/captures';B=O/'current-metadata-r3/captures'
src=Path('target/wp3-drawing-transplant/frame-digest-boundary.py').read_text();env={};exec(src[:src.index('expected=')],env);spellings=env['spellings']
def load(p):return json.loads(p.read_text())
def norm(v):
 if isinstance(v,str):
  for a,b in spellings.items():v=v.replace(a,b)
  return v
 if isinstance(v,list):return [norm(x) for x in v]
 if isinstance(v,dict):return {k:norm(x) for k,x in v.items()}
 return v
def diffs(a,b,path=()):
 if type(a)!=type(b):yield path,a,b;return
 if isinstance(a,dict):
  for k in sorted(a.keys()|b.keys()):
   if k not in a or k not in b:yield path+(k,),a.get(k),b.get(k)
   else:yield from diffs(a[k],b[k],path+(k,))
 elif isinstance(a,list):
  if len(a)!=len(b):yield path+('length',),len(a),len(b)
  for i,(x,y) in enumerate(zip(a,b)):yield from diffs(x,y,path+(i,))
 elif a!=b:yield path,a,b
report=[]
for name in ['ym_ggg_zero_color_L3','xs_qcdqed_addx_canonical_scalar_lr_L3']:
 if not (A/name/'raw-report.txt').exists():continue
 row={'name':name,'model_differences':list(diffs(load(A/name/'model.json'),norm(load(B/name/'model.json'))))[:8],'options_equal':load(A/name/'options.json')==load(B/name/'options.json'),'grouping_equal':load(A/name/'grouping.json')==load(B/name/'grouping.json'),'stages':{}}
 for stage in ['raw','filtered','grouped']:
  if not (A/name/(stage+'-report.txt')).exists():continue
  a=load(A/name/(stage+'.json'));b=norm(load(B/name/(stage+'.json')));aa={d['id']:d for d in a['diagrams']};bb={d['id']:d for d in b['diagrams']};assert len(aa)==len(a['diagrams']);assert len(bb)==len(b['diagrams'])
  shared=sorted(aa.keys()&bb.keys());counter=collections.Counter();witnesses={};allmatches=0
  for key in shared:
   ds=list(diffs(aa[key],bb[key]));allmatches+=not ds
   for path,x,y in ds:
    coarse='/'.join(map(str,path));counter[coarse]+=1
    if coarse not in witnesses:witnesses[coarse]={'id':key,'path':path,'baseline':x,'current':y}
  ma=dict(load(A/name/(stage+'-metadata.json')));mb=dict(norm(load(B/name/(stage+'-metadata.json'))));md=list(diffs(ma,mb))
  sr={'counts':[len(aa),len(bb)],'exact_id_sets_equal':aa.keys()==bb.keys(),'only_baseline':sorted(aa.keys()-bb.keys()),'only_current':sorted(bb.keys()-aa.keys()),'full_records_equal_after_explicit_head_mapping':allmatches,'difference_paths':dict(counter),'first_witnesses':list(witnesses.values())[:10],'logical_metadata_difference_count':len(md),'first_logical_metadata_differences':md[:5],'group_payload_equal':a['groups']==b['groups'],'first_group_differences':list(diffs(a['groups'],b['groups']))[:5]}
  row['stages'][stage]=sr
 report.append(row)
(O/'stage-comparison.json').write_text(json.dumps({'normalization':'ONLY five explicitly verified presentation-head mappings, no interface omission, algebra or factor expansion','head_map':spellings,'rows':report},indent=2)+'\n')
for r in report:
 print(r['name'],'model diffs',str(r['model_differences'])[:350],'settings',r['options_equal'],r['grouping_equal'])
 for s,v in r['stages'].items():
  print(s,v['counts'],'idset=',v['exact_id_sets_equal'],'records=',v['full_records_equal_after_explicit_head_mapping'],'metadata diffs=',v['logical_metadata_difference_count'],'groups=',v['group_payload_equal']); print('paths',list(v['difference_paths'].items())[:10]);print('witness',str(v['first_witnesses'][:1])[:900])
