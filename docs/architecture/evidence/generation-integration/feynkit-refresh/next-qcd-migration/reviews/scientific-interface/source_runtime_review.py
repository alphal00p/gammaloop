import ast,sys,pathlib,json,hashlib,struct,re,collections
W=pathlib.Path('/shared/lcnbr/dev/gammaloop/generation-oracle-numerator-families');D=pathlib.Path(__file__).resolve().parent;R=W/'target/raw-next-qcd-proof-preparation-20261031-r3';P=W/'target/raw-next-qcd-proof-preparation-20261031-r1';E=W/'target/raw-next-qcd-proof-native-execution-20261031-r1';C=W/'target/raw-next-qcd-preflight-20261031-r1/capsule'
sys.path.insert(0,str(R));import pure_parser as B
bindings={}
def read(p):
 p=pathlib.Path(p);b=p.read_bytes();bindings[str(p)]={'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)};return b
def j(p):return json.loads(read(p))
def sha(p):
 with pathlib.Path(p).open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(n,x):(D/n).write_text(json.dumps(x,indent=2)+'\n')
def node(p,name):return next(n for n in ast.parse(read(p)).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))and n.name==name)
functions=['owner_table','physical','shape','equal','need','defs','Reject'];sources={}
for n in functions:assert ast.dump(node(R/'strict_owners.py',n),include_attributes=False)==ast.dump(node(P/'sources/original/proof.py',n),include_attributes=False)
for n in ['render','definitions','cyclic','owners','permutation','Reject','require']:assert ast.dump(node(R/'pure_parser.py',n),include_attributes=False)==ast.dump(node(P/'sources/original/syntax.py',n),include_attributes=False)
assert read(R/'prepare.py')==read(R.parent/'raw-next-qcd-proof-preparation-20261031-r2/prepare.py')
assert read(R/'fold_plan.py')==read(R.parent/'raw-next-qcd-proof-preparation-20261031-r2/fold_plan.py')
assert read(R/'audit_native.py')==read(R.parent/'raw-next-qcd-proof-preparation-20261031-r2/audit_native.py')
for rel,meta in j(R/'unchanged-r1-prepared-data.json').items():assert sha(P/rel)==meta['sha256']
for a in j(P/'source-anchor-revalidation.json'):
 text=read(a['actual_file']).decode();assert hashlib.sha256(text.encode()).hexdigest()==a['new_source_sha256'];assert [i+1 for i,s in enumerate(text.splitlines())if a['anchor']in s]==a['actual_anchor_lines'];assert read(P/a['packet_copy'])==text.encode()
# Exact definition source constructors, not erasure of arbitrary tags.
for route in ['original_reference','current']:
 for rel in ['crates/idenso/src/color/mod.rs','crates/feynkit-graph/src/symbols.rs','crates/idenso/src/tensor/inference.rs','crates/idenso/src/tensor/mod.rs']:
  text=read(C/'source'/route/rel).decode()
  if rel.endswith('color/mod.rs'):assert 'f: tensor_symbol!("spenso::f"; Real, Antisymmetric;'in text
  if rel.endswith('symbols.rs'):
   for n in ['Q','K','P','edge','vertex']:assert '"gammalooprs::'+n+'"'in text
   assert all(s in text for s in(['spenso::tensor-label:q','spenso::tensor-label:k','spenso::tensor-label:p','spenso::index-label:e','spenso::index-label:v']if route=='original_reference'else['Notation','IndexDisplay']))
# Helper source identity and qualified ELFs, all bytes only.
typedq=j(C/'helpers/interface/3-runtime-qualification-v4.json');assert sha(P/'sources/original/interface-0-interface_probe_v4.rs')==typedq['helper_source_sha256']=='362a7fa3edce435e7957350e8f2b97da6e92a1bead4a1a450953ae6cb095afd4'
read(P/'sources/original/interface-0-interface_probe_v4.rs')
assert read(C/'helpers/c78-exact/oracle_family_exact.rs')==read(P/'sources/current8eb-exact-helper.rs')
programs=j(P/'native-tool-and-closure-bindings.json')['programs']
for q in programs.values():assert sha(q['path'])==q['sha256']
assert programs['typed8eb']['sha256']==typedq['helper_sha256']
currentq=j(W.parent/'generation-feynkit-refresh/.local/feynkit-refresh/qualification/helper-build/oracle_family_exact/qualification.json');assert currentq['sha256']==programs['8eb']['sha256']
for dep in currentq['dependencies'].values():assert sha(dep['artifact'])==dep['sha256']
# Rebuild every explicit typed parse event and first-use definition sequence.
inputs=j(P/'prepared/interface-inputs.json');schedule=j(P/'prepared/interface-definition-order.json')['graphs'];si=0;events_count=0
for cell in inputs['cells']:
 for g in cell['graphs']:
  s=schedule[si];si+=1;assert(s['case'],s['mode'],s['id'])==(cell['case'],cell['mode'],g['id'])
  for phase,side in zip(s['phases'],['current_payload','mapped_historical_payload']):
   events=[]
   def visit(v,path=()):
    if isinstance(v,dict):
     for k in sorted(v):
      if k in ('num_interface','numerator_prefactor_interface','projector_interface'):
       events.extend({'stage':'serde_interface','path':list(path+(k,i)),'text':x}for i,x in enumerate(v[k]))
      else:visit(v[k],path+(k,))
    elif isinstance(v,list):
     for i,x in enumerate(v):visit(x,path+(i,))
   p=g[side];visit(p)
   for name in ['overall_factor','numerator','numerator_prefactor','projector']:events.append({'stage':'from_serde_expression','owner':name,'text':p[name]})
   for i,v in enumerate(p['vertices']):events.append({'stage':'from_serde_expression','owner':f'vertex:{i}','text':v['numerator']})
   for i,(_,v)in enumerate(p['edges']):events.append({'stage':'from_serde_expression','owner':f'edge:{i}','text':v['numerator']})
   defs={};order=[]
   for e in events:
    for name,token in B.definitions(B.parse(e['text'])).items():
     assert name not in defs or defs[name]==token;defs[name]=token
     if token not in order:order.append(token)
   assert phase=={'payload':side,'events':events,'whole_definitions':defs,'explicit_text_first_use_order':order};events_count+=len(events)
# Parse ALL actual changed-owner lowering permutations, not just dangling multiset.
typed=j(R/'execution/records/all-six-cells--typed8eb/results.json');lower=[]
for c in typed['cells']:
 for g in c['graphs']:
  for r in g['owners']:
   if not r['affected']:continue
   for route in ['historical','current']:
    text=r[route+'_lowering']['layout'];rank=len(r[route]['ports']);identity=re.fullmatch(r'CanonicalLayout\(Identity \{ order: ([0-9]+) \}\)',text)
    if identity:assert int(identity[1])==rank;p=list(range(rank));inv=p
    else:
     m=re.fullmatch(r'CanonicalLayout\(Mapped \{ representation: RepresentationOrder\(Permutation \{ map: (\[[0-9, ]*\]), inv: (\[[0-9, ]*\]) \}\), index: IndexOrder\(Permutation \{ map: (\[[0-9, ]*\]), inv: (\[[0-9, ]*\]) \}\) \}\)',text);assert m,text
     rp,ri,ip,ii=[json.loads(x)for x in m.groups()]
     for p0,i0 in [(rp,ri),(ip,ii)]:assert sorted(p0)==sorted(i0)==list(range(rank))and all(i0[p0[k]]==k for k in range(rank))
     p=[ip[rp[k]]for k in range(rank)];inv=[ri[ii[k]]for k in range(rank)]
    assert inv==r[route]['canonical_to_logical'];lower.append({'case':c['case'],'mode':c['mode'],'id':g['id'],'owner':r['owner'],'route':route,'logical_to_canonical':p,'canonical_to_logical':inv})
assert len(lower)==420
# Static ELF metadata replay; extracting only the pure elf function avoids any runtime imports.
elfnode=node(E/'qualification.py','elf');ns={'Path':pathlib.Path,'struct':struct,'sha':sha};exec(compile(ast.Module(body=[elfnode],type_ignores=[]),'frozen-pure-elf-reader','exec'),ns)
qualification=j(E/'elf-qualification-corrected.json');objects={};cards={x['id']:x for x in j(R/'native-commands.json')['commands']};preflights=[]
for cfg in qualification['configurations']:
 assert cfg['missing']==[];card=cards[cfg['id']];assert cfg['program']==card['program'];argv=card['argv'];mounts=[[argv[i+1],argv[i+2]]for i,x in enumerate(argv)if x=='--ro-bind'];assert mounts==cfg['readonly_mounts']
 for p,meta in cfg['objects'].items():
  if p not in objects:objects[p]=ns['elf'](p)
  assert objects[p]==meta
 root=E/'preflight'/cfg['id'];probe=j(root/'work/probe.json');assert probe['readonly_mounts']==[x[1]for x in mounts]
 for q in probe['bytes']:
  path=q['path'];mapped=next(pathlib.Path(src)/path[len(dst):].lstrip('/')for src,dst in reversed(mounts)if path==dst or path.startswith(dst+'/'));assert sha(mapped)==q['sha256']
 assert probe['network_interfaces']==['lo']and probe['network_probe_errno']==101
 assert probe['private_handoff_qualified']is False and probe['scientific_helper_executed']is False
 for label in ['loader-verify','stdlib-containment']:
  life=j(root/(label+'.lifecycle.json'));cmd=j(root/(label+'.command.json'));assert life['exit']==0 and not life['timed_out']and life['process_reaped']and life['cleanup']['remaining_pids']==[]
  assert (root/(label+'.stderr')).read_bytes()==b'';assert cmd['environment']==card['environment_before_private_exec']and cmd['scientific_helper_executed']is False and cmd['license_acquired']is False
  # Prefix has identical sandbox/read-only settings, only its private work directories differ.
  divider=argv.index('--');prefix=[str(root/'work')if x==card['writable_root']else str(root/'work/home')if x==str(pathlib.Path(card['writable_root'])/'home')else x for x in argv[:divider+1]];assert cmd['argv'][:len(prefix)]==prefix
  preflights.append({'id':cfg['id'],'probe':label,'exit':life['exit'],'pid_role':'outer bwrap launcher','runtime_private_credit':False})
read(E/'elf-qualification.json');read(E/'elf-resolution.stderr');read(E/'source-and-setup-errors.json')
save('source-runtime-independent.json',{'source_bindings':bindings,'strict_original_predicates_equal':functions,'prepared_files_unchanged':87,'typed_explicit_parse_schedules':si,'typed_explicit_parse_events':events_count,'all_lowering_layouts_checked':len(lower),'lowering_layouts':lower,'static_ELF_objects':objects,'genuine_recorded_preflight_checks':preflights,'qualification_limit':'Static ELF/source checks and genuine prior unlicensed loader/stdlib checks are separate from actual21 native phase. Historical compiler/build/store trust and source/profile/observer gaps remain premises. No new loader, helper or network execution. Typed initializers/model registration are source premises, not invented global interner observations.'})
print('source/runtime reconstruction passed',si,events_count,len(lower),len(objects),flush=True)
