from pathlib import Path
import json,collections
O=Path('target/wp3-units19-21/frozen-row-diagnosis');env={};s=Path('target/wp3-drawing-transplant/frame-digest-boundary.py').read_text();exec(s[:s.index('expected=')],env);heads=env['spellings']
def norm(s):
 for a,b in heads.items():s=s.replace(a,b)
 return s
counts=collections.Counter()
with (O/'all-raw-axis-bijections.jsonl').open('w') as f:
 for case in ['ym_ggg_zero_color_L3','xs_qcdqed_addx_canonical_scalar_lr_L3']:
  a=dict(json.loads((O/'baseline-decoded/captures'/case/'raw-metadata.json').read_text()));b=dict(json.loads((O/'current-decoded/captures'/case/'raw-metadata.json').read_text()));pa={x['id']:x['payload'] for x in json.loads((O/'baseline-decoded/captures'/case/'raw.json').read_text())['diagrams']};pb={x['id']:x['payload'] for x in json.loads((O/'current-decoded/captures'/case/'raw.json').read_text())['diagrams']}
  for id in sorted(a):
   av={x[0]:x[1:] for x in a[id]};bv={x[0]:x[1:] for x in b[id]}
   for owner,(ports,m,c) in av.items():
    new,m2,c2=bv[owner];new=[norm(x) for x in new];assert collections.Counter(ports)==collections.Counter(new) and len(set(ports))==len(ports);assert(m,c)==(m2,c2)
    for payload,logical in [(pa[id],ports),(pb[id],new)]:
     if owner.startswith('vertex:'):decl=payload['vertices'][int(owner.split(':')[1])].get('num_interface')
     elif owner.startswith('edge:'):decl=payload['edges'][int(owner.split(':')[1])][1].get('num_interface')
     else:decl=payload.get({'prefactor':'numerator_prefactor_interface','projector':'projector_interface'}[owner])
     if decl is not None:assert [norm(x) for x in decl]==logical;counts['explicit_source_interfaces_checked']+=1
    row={'case':case,'graph_id':id,'owner':owner,'old_logical_ports':ports,'current_logical_ports_in_legacy_head_spelling':new,'old_axis_to_current_axis':[new.index(p) for p in ports],'metric':m,'composite':c};f.write(json.dumps(row)+'\n');counts['raw_tensors']+=1;counts['axes']+=len(ports);counts['permuted_owners']+=ports!=new
(O/'axis-proof-summary.json').write_text(json.dumps(dict(counts),indent=2)+'\n');print(dict(counts))
