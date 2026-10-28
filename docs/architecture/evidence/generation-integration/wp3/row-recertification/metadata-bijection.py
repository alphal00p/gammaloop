from pathlib import Path
import json,collections
O=Path('target/wp3-units19-21/frozen-row-diagnosis');env={};s=Path('target/wp3-drawing-transplant/frame-digest-boundary.py').read_text();exec(s[:s.index('expected=')],env);heads=env['spellings']
def norm(x):
 for a,b in heads.items():x=x.replace(a,b)
 return x
result=[];witness=[]
for case in ['ym_ggg_zero_color_L3','xs_qcdqed_addx_canonical_scalar_lr_L3']:
 for stage in ['raw','filtered','grouped']:
  a=dict(json.loads((O/'baseline-decoded/captures'/case/(stage+'-metadata.json')).read_text()));b=dict(json.loads((O/'current-decoded/captures'/case/(stage+'-metadata.json')).read_text()));assert a.keys()==b.keys()
  counts=collections.Counter();bad=[]
  for id in sorted(a):
   aa={v[0]:v[1:] for v in a[id]};bb={v[0]:v[1:] for v in b[id]};assert aa.keys()==bb.keys()
   for owner in aa:
    ports,metric,composite=aa[owner];new,m2,c2=bb[owner];new=[norm(x) for x in new];counts['local_tensors']+=1
    if collections.Counter(ports)!=collections.Counter(new):bad.append({'id':id,'owner':owner,'before':ports,'after':new});continue
    assert len(set(ports))==len(ports),(case,id,owner,ports)
    permutation=[new.index(p) for p in ports];counts['port_bijections']+=1;counts['axes']+=len(ports);counts['permuted']+=ports!=new
    if (metric,composite)!=(m2,c2):counts['classification_flag_changes']+=1
    if ports!=new and len(witness)<6:witness.append({'case':case,'stage':stage,'id':id,'owner':owner,'old_ports':ports,'new_ports':new,'old_axis_to_new_axis':permutation,'old_flags':[metric,composite],'new_flags':[m2,c2]})
  result.append({'case':case,'stage':stage,'counts':dict(counts),'bad_port_sets':bad})
(O/'metadata-bijection.json').write_text(json.dumps({'method':'Decode each live historical index using its existing CookSettings::uncook owner BEFORE serialization; compare complete labelled/dual/dimensional ports, retain explicit axis permutation rather than omit interfaces. Namespaced constructor mapping only; no representation/label erasure.','rows':result,'permutation_witnesses':witness},indent=2)+'\n')
for r in result:print(r['case'],r['stage'],r['counts'],'badsets',len(r['bad_port_sets']))
