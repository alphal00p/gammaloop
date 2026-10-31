import sys,json,pathlib,hashlib,copy,datetime,collections,re
W=pathlib.Path('/shared/lcnbr/dev/gammaloop/generation-oracle-numerator-families'); D=pathlib.Path(__file__).resolve().parent;R=W/'target/raw-next-qcd-proof-preparation-20261031-r3';E=W/'target/raw-next-qcd-proof-native-execution-20261031-r1';P=W/'target/raw-next-qcd-proof-preparation-20261031-r1'
sys.path.insert(0,str(R));import supervisor as S

def j(p):return json.loads(pathlib.Path(p).read_bytes(),object_pairs_hook=S.A.unique)
def sha(p):
 with pathlib.Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(n,x):(D/n).write_text(json.dumps(x,indent=2)+'\n')
b=j(E/'validated-context.json');cards=j(R/'native-commands.json')['commands'];assert b['cards']==cards
contracts=j(R/'result-contracts.json');orders=j(R/'ordered-input-bindings.json');assignment=j(b['assignment']['assignment_path'])
assert [c['id'] for c in cards]==assignment['command_ids']==[x['id'] for x in j(E/'actual-invocations.json')]
assert b['all_contracts']==contracts and b['all_definition_orders']==orders
for f,k in [('supervisor.py','supervisor_source_sha256'),('supervisor_guard.py','supervisor_guard_sha256'),('native-commands.json','native_commands_sha256'),('result-contracts.json','result_contracts_sha256'),('ordered-input-bindings.json','ordered_input_bindings_sha256')]:assert sha(R/f)==assignment[k]
assert sha(R/'SEAL.json')==assignment['packet_seal_sha256']==b['supervisor_seal_sha256']
assert sha(b['assignment']['assignment_path'])==b['assignment']['assignment_sha256']
assert sha(b['assignment']['review_path'])==b['assignment']['review_sha256']==assignment['supervisor_review']['sha256']
all_rows=[];protected=j(E/'inputs-before.json');previous_end=None
for i,(card,actual,ap) in enumerate(zip(cards,j(E/'actual-invocations.json'),assignment['program_input_order_result_bindings'])):
 ident=card['id'];r=pathlib.Path(card['records_root']);receipt=j(r/'result.json');cmd=j(r/'command.json');life=receipt['lifecycle'];ct=contracts[ident]
 assert receipt==j(R/'execution/outcomes'/f'{ident}.json')
 bb={**b,'card':card,'contract':ct};expected=S._binding(bb,card,'NATIVE_ASSIGNED',b['assignment'])
 for key,val in expected.items():assert receipt[key]==val and cmd[key]==val,(ident,key)
 assert ap=={'id':ident,'program':card['program'],'input':expected['input'],'card_sha256':expected['card_sha256'],'ordered_input_binding':orders[ident],'result_contract':ct}
 assert sha(card['program']['path'])==card['program']['sha256'] and sha(card['input'])==card['input_sha256']
 assert cmd['card']==card and cmd['protected_inputs_before']==b['protected_inputs']==receipt['inputs_after']
 assert all(protected[p]['sha256']==h for p,h in b['protected_inputs'].items())
 assert receipt['deadline_seconds']==240 and receipt['retries']==0 and not receipt.get('error')
 for key,val in [('native_launched',True),('authenticated_execution',True),('validated_predicate',True),('process_reaped',True),('launch_attempted',True),('all_input_hashes_unchanged',True),('synthetic',False),('fixture_only',False),('fixture_python_launched',False)]:assert type(receipt[key])is type(val) and receipt[key]==val
 lc=S._launch_contract(card,False);assert cmd['launch_contract']==lc
 assert set(life)==set(lc)|{'exit','timed_out','pid','cleanup','seconds','process_reaped'}
 for k,v in lc.items():assert life[k]==v
 assert type(life['pid'])is int and life['pid']>0 and life['process_reaped']is True and life['timed_out']is False and life['cleanup']['remaining_pids']==[]
 assert 0<life['seconds']<240
 arts={'command.json','stdout','stderr'}|({'results.json'} if ct['classification']!='math-negative' else set())
 assert set(receipt['artifact_hashes'])==arts
 for n,h in receipt['artifact_hashes'].items():assert sha(r/n)==h
 assert actual['index']==i and actual['retries']==0 and actual['controller_exit']==actual['audit_exit']==0
 assert actual['argv'][4:]==[str(R/'supervisor.py'),'--id',ident,'--execute-licensed','--assignment',b['assignment']['assignment_path']]
 start=datetime.datetime.fromisoformat(actual['start']);end=datetime.datetime.fromisoformat(actual['end']);assert start<end
 assert previous_end is None or previous_end<start;previous_end=end
 status,audit=S.classify(bb,life,r);assert status==receipt['status']==actual['status'] and audit==receipt['audit']
 assert receipt['semantic_negative_credit']==(i>=5)
 if i>=5:assert all(x in [c['id']for c in cards[:5]]for x in ct['dependencies'])
 out=(r/'stdout').read_text();err=(r/'stderr').read_text()
 if ct['classification']=='math-baseline':assert out==f"{ct['rows']} exact factorized expression pairs\n" and err==''
 if ct['classification']=='typed-baseline':assert out==err==''
 all_rows.append({'id':ident,'kind':ct['classification'],'status':status,'actual_exit':life['exit'],'seconds':life['seconds'],'program':card['program'],'input':expected['input'],'receipt_sha256':sha(r/'result.json'),'audit':audit,'stdout':out,'stderr':err})
save('execution-independent.json',{'outcomes':all_rows,'ordered_count':len(all_rows),'baselines':5,'negative_predicates':16,'protected_inputs_per_card':len(b['protected_inputs']),'limit':'Exact cooperative command/source/input/order/result/lifecycle consistency, not hostile attestation. Replayed only pure classifiers; no native or private calls.'})
print('execution metadata and semantic replays passed',len(all_rows),flush=True)
# Authenticate EVERY structured association, and reconstruct inverse references from source obligations.
seen={};documents={};ref_count=0

def resolve(ref):
 global ref_count
 p=ref['path'];ref_count+=1
 if p not in seen:seen[p]=(sha(p),pathlib.Path(p).stat().st_size)
 assert seen[p]==(ref['sha256'],ref['bytes'])
 if 'json_pointer' not in ref:return None
 if p not in documents:documents[p]=j(p)
 v=documents[p]
 for t in ref['json_pointer'].split('/')[1:]:
  t=t.replace('~1','/').replace('~0','~');v=v[int(t)] if isinstance(v,list) else v[t]
 return v

def walk(v):
 if isinstance(v,dict):
  if {'path','sha256','bytes'}<=set(v):resolve(v)
  for a in v.values():walk(a)
 elif isinstance(v,list):
  for a in v:walk(a)
for f in ['math-member-to-native-row.json','typed-native-associations.json','prepared-data-associations.json','outcome-index.json','execution-artifacts.json']:walk(j(E/f))
inverse=j(E/'math-member-to-native-row.json');rebuilt=[]
for batch in j(P/'prepared/native-math-batches.json'):
 records=j(P/'prepared/math'/batch['batch']/'obligations.json')
 for kernel in batch['kernels']:
  native=j(R/'execution/records'/(batch['batch']+'--'+kernel)/'results.json')
  for row,(record,answer) in enumerate(zip(records,native)):
   for ai,assoc in enumerate(record['obligations']):
    for member in assoc.get('owners',[None]):rebuilt.append((batch['batch'],kernel,row,ai,member,assoc.get('obligation'),answer))
assert len(inverse)==len(rebuilt)==4268
for ref,v in zip(inverse,rebuilt):
 batch,kernel,row,ai,member,kind,answer=v
 assert (ref['batch'],ref['kernel'],ref['row'],ref['member'],ref['obligation_kind'],ref['native_predicates'])==(batch,kernel,row,member,kind,answer)
 assert ref['obligation']['json_pointer']==f'/{row}/obligations/{ai}'
 assert resolve(ref['result'])==answer
 assert member is None or member in resolve(ref['obligation'])['owners']
typed=j(R/'execution/records/all-six-cells--typed8eb/results.json');inputs=j(P/'prepared/interface-inputs.json');tas=j(E/'typed-native-associations.json');i=0
for ci,(c,nc) in enumerate(zip(inputs['cells'],typed['cells'])):
 for gi,(g,ng) in enumerate(zip(c['graphs'],nc['graphs'])):
  a=tas[i];i+=1;assert (a['case'],a['mode'],a['id'])==(c['case'],c['mode'],g['id'])
  assert resolve(a['exact_supplied_graph'])==g and resolve(a['full_native_graph'])==ng
  assert [o['owner']for o in a['owners']]==[o['owner']for o in g['owners']]
  for oi,o in enumerate(a['owners']):assert resolve(o['actual_import_export_flags_operand_maps_and_lowering'])==ng['owners'][oi]
  for k,n in [('actual_pair_plans','plans'),('actual_fold_plans_survivors_canonical_layout_and_output_axes','actual_fold'),('public_full_or_empty_numerator_of','public_numerator_of')]:assert resolve(a[k])==ng[n]
assert i==48
save('associations-independent.json',{'reference_instances_verified':ref_count,'unique_files':seen,'math_inverse_references':len(inverse),'typed_graph_associations':i,'typed_owner_associations':sum(len(g['owners'])for c in inputs['cells']for g in c['graphs'])})
print('association references passed',ref_count,flush=True)
