import sys,json,pathlib,copy,hashlib,shutil,collections
W=pathlib.Path('/shared/lcnbr/dev/gammaloop/generation-oracle-numerator-families');D=pathlib.Path(__file__).resolve().parent;R=W/'target/raw-next-qcd-proof-preparation-20261031-r3';P=W/'target/raw-next-qcd-proof-preparation-20261031-r1';E=W/'target/raw-next-qcd-proof-native-execution-20261031-r1'
sys.path.insert(0,str(R));import prepare as A;import supplement as T;import audit_native as N;import supervisor as S

def j(p):return json.loads(pathlib.Path(p).read_bytes())
def save(n,x):(D/n).write_text(json.dumps(x,indent=2)+'\n')
def diffs(a,b,p=()):
 if type(a)!=type(b):return [(p,a,b)]
 if isinstance(a,dict):
  assert set(a)==set(b);return [x for k in a for x in diffs(a[k],b[k],p+(k,))]
 if isinstance(a,list):
  assert len(a)==len(b);return [x for k,(u,v)in enumerate(zip(a,b))for x in diffs(u,v,p+(k,))]
 return [(p,a,b)]if a!=b else []
inputs=j(P/'prepared/interface-inputs.json');native=j(R/'execution/records/all-six-cells--typed8eb/results.json');native_controls=[]
for control in j(P/'prepared/native-controls.json'):
 doc=j(P/control['input'])
 if control['kind']=='typed-negative':
  cell=next(c for c in inputs['cells']if(c['case'],c['mode'])==(control['case'],control['mode']));g=next(g for g in cell['graphs']if g['id']==control['graph']);base={'cells':[copy.deepcopy(cell)]};base['cells'][0]['graphs']=[copy.deepcopy(g)];ds=diffs(base,doc);assert len(ds)==1;p,a,b=ds[0];assert a==control['before']and b==control['after']and a!=b
  if 'owner-claim'in control['id']:assert 'owners'in p and p[-1]in('is_metric','is_composite') and b is (not a)
  else:
   assert 'current_payload'in p and 'num_interface'in p and p[-1]==control['axis']
   if 'actual-dimension'in control['id']:
    slot=A.slot(a);rep,dim=slot['key'][:2];assert b==a.replace('::'+rep+'('+str(dim)+',','::'+rep+'('+str(dim+1)+',',1)
   elif 'wrong-owner'in control['id']:assert b==a.replace('::hedge(','::proof_wrong_owner(',1)
   else:assert b==a[len('spenso::{spenso::lower}::dind('):-1]
  result=j(R/'execution/records'/(control['id']+'--typed8eb')/'results.json');err=result['cells'][0]['error'];native_controls.append({'id':control['id'],'changed_path':list(p),'before':a,'after':b,'actual_graph':err['graph'],'actual_predicate':err['error'],'kind':'actual imported owner claim'if'owner-claim'in control['id']else'actual used current interface'})
 else:
  pre=j(P/'prepared/math/nonliteral-cross-and-aggregate/portable-pairs.json')[0];kind=control['id'].removeprefix('math-');a=pre['historical'];assert doc==[{'historical':{'sign':'-('+a+')','coefficient':'2*('+a+')','power':'('+a+')^2'}[kind],'current':pre['current']}]
  old=j(P/'prepared/math/nonliteral-cross-and-aggregate/obligations.json')[0];merged=j(P/'prepared/math/current-faithfulness-and-equivalence/obligations.json');idx=next(i for i,x in enumerate(merged)if(x['historical'],x['current'])==(old['historical'],old['current']))
  left,right=[A.B.parse(pre[s])for s in('historical','current')];assert left[0]==right[0]=='+'and left[2]==right[2]and left[2][0]=='symbol';z=left[2][1]
  assert z not in A.B.definitions(left[1]).values() and z not in A.B.definitions(right[1]).values()
  assert A.B.render(left[1])==A.B.render(A.B.parse(old['historical']))
  native_controls.append({'id':control['id'],'parent_equivalence_row':idx,'kernels':['c78','8eb'],'fresh_symbol':z,'actual_predicate':'assertion failed: row.parsed_equal && row.canonical_equal','nondegeneracy':'formal polynomial in z has coefficient 1; sign/coefficient changes it, square gives z-degree 2. Distinct fresh witness names in merged baseline and retained control preserve the identical base relation.'})
save('actual-negative-inputs-independent.json',native_controls)
# Fresh pure heldouts on private in-memory copies; never native receipts.
probes=[]
def probe(name,fn,stage=None):
 try:fn()
 except (A.G.Reject,A.B.Reject) as e:
  item={'name':name,'result':'REJECT','stage':e.stage,'message':str(e),'expected_stage':stage};probes.append(item);assert e.stage==stage,item
 except Exception as e:
  item={'name':name,'result':'HARNESS_ERROR','type':type(e).__name__,'message':str(e),'expected_stage':stage};probes.append(item);save('fresh-heldouts-progress.json',probes);raise
 else:
  probes.append({'name':name,'result':'ACCEPT','expected_stage':stage});assert stage is None,probes[-1]
 save('fresh-heldouts-progress.json',probes)
def typed_change(name,mut,stage):
 q=copy.deepcopy(native);mut(q);probe(name,lambda:N.audit_typed(q),stage)
probe('authentic typed baseline pure replay',lambda:N.audit_typed(native))
typed_change('last cell removed',lambda q:q['cells'].pop(),'native-inventory')
typed_change('last owner removed',lambda q:q['cells'][-1]['graphs'][-1]['owners'].pop(),'native-inventory')
typed_change('imported metric flag flipped',lambda q:q['cells'][-1]['graphs'][-1]['owners'][0]['current'].__setitem__('is_metric',not q['cells'][-1]['graphs'][-1]['owners'][0]['current']['is_metric']),'native-flags')
typed_change('last fold plan absent',lambda q:q['cells'][-1]['graphs'][-1]['actual_fold'][-1].pop('physical_plan'),'native-fold-plan')
typed_change('last fold plan kind wrong',lambda q:q['cells'][-1]['graphs'][-1]['actual_fold'][-1]['physical_plan'].__setitem__('kind','Compose'),'native-fold-kind')
def swap_cof(q):
 for c in q['cells']:
  for g in c['graphs']:
   for s in g['actual_fold']:
    ns=s['physical_plan'].get('names',[]);indices=[i for i,(a,b)in enumerate(ns)if '::cof('in a]
    if len(indices)>=2:
     i,k=indices[:2];ns[i][1],ns[k][1]=ns[k][1],ns[i][1];return
 raise AssertionError('no cof pairs')
typed_change('real fold cof right endpoints swapped',swap_cof,'native-fold-endpoint-orientation')
def duplicate(q):
 s=next(s for s in q['cells'][-1]['graphs'][-1]['actual_fold']if s['physical_plan']['kind']=='Contract');s['physical_plan']['names'].append(s['physical_plan']['names'][0])
typed_change('actual fold duplicate endpoint',duplicate,'native-fold-endpoint-duplicate')
def output_swap(q):
 r=q['cells'][-1]['graphs'][-1]['public_numerator_of'];r['canonical_to_logical'][0],r['canonical_to_logical'][1]=r['canonical_to_logical'][1],r['canonical_to_logical'][0]
typed_change('public output-only inverse permutation swapped',output_swap,'native-layout')
def survivor_loss(q):
 r=q['cells'][-1]['graphs'][-1]['actual_fold'][-1]['output'];r['ports'].pop()
typed_change('unmatched coad output omitted',survivor_loss,'native-fold')
def lowering_loss(q):
 r=next(r for c in q['cells']for g in c['graphs']for r in g['owners']if r['affected']);r['current_lowering']['dangling'].pop()
typed_change('changed-owner lowering dangling output missing',lowering_loss,'native-lowering')
math=j(R/'execution/records/current-faithfulness-and-equivalence--c78/results.json');probe('authentic math baseline replay',lambda:N.audit_math(math,118))
for field in math[-1]:
 q=copy.deepcopy(math);q[-1][field]=False;probe('last native row '+field,lambda q=q:N.audit_math(q,118),'native-math')
q=copy.deepcopy(math);q[-1]['parsed_equal']=1;probe('boolean replaced by integer',lambda:N.audit_math(q,118),'native-math')
# Source/capture guards and explicit congruence heldouts; save is replaced in memory.
A.save=lambda n,x:None
for field,val,stage in [('case','xs_qcd_gguu_L1','live-case'),('mode','None','live-mode')]:
 w={r:A.load(A.C/f'prepared/captured/xs_qcd_gggg_L1/{r}-Auto.json')for r in A.ROUTES};w['current']['live'][field]=val;probe('live '+field+' binding',lambda w=w:A.analyze('xs_qcd_gggg_L1','Auto',w),stage)
w={r:A.load(A.C/f'prepared/captured/xs_qcd_gguu_sym_L1/{r}-None.json')for r in A.ROUTES};w['current']['payload']['diagrams'][0]['payload']['symmetry_factor']=[9,1];probe('physical symmetry factor changed',lambda:A.analyze('xs_qcd_gguu_sym_L1','None',w),'cross-physical')
f='spenso::{antisymmetric,real,spenso::tensor}::f';a=A.B.parse(f+'(x::a,x::b,x::c)');b=A.B.parse(f+'(x::b,x::a,x::c)');probe('odd attributed f transposition',lambda:T.cycle_trace(a,b),'cycle-proof')
probe('whole-definition conflict',lambda:A.B.definitions(A.B.parse('spenso::{real}::f(x::a)+spenso::{antisymmetric,real}::f(x::a)')),'symbol-definition')
# Copy a real protected baseline to a labelled reader-test tree. No runtime engine is used.
root=D/'receipt-reader-copies';ident='all-six-cells--typed8eb';src=R/'execution/records'/ident;dst=root/'records'/ident;dst.mkdir(parents=True)
for f in ['command.json','result.json','results.json','stdout','stderr']:shutil.copyfile(src/f,dst/f)
(root/'purpose.json').write_text(json.dumps({'purpose':'Private copies of existing actual receipts for pure reader controls only; no new execution or authentication claim','source':str(src)}))
bundle=j(E/'validated-context.json');bundle['execution_root']=str(root);bundle['contract']={'dependencies':[ident]};probe('unchanged actual receipt copy',lambda:S._check_dependencies(bundle,False))
original=(dst/'result.json').read_bytes()
for label,mut,stage in [('missing result artifact digest',lambda r:r['artifact_hashes'].pop('results.json'),'baseline-artifact-inventory'),('missing postinput',lambda r:r['inputs_after'].pop(next(iter(r['inputs_after']))),'baseline-inputs-after'),('native receipt fixture marker',lambda r:r.__setitem__('synthetic',True),'baseline-flags')]:
 r=json.loads(original);mut(r);(dst/'result.json').write_text(json.dumps(r));probe(label,lambda:S._check_dependencies(bundle,False),stage);(dst/'result.json').write_bytes(original)
raw=(dst/'results.json').read_bytes();(dst/'results.json').write_bytes(raw+b'\n');probe('changed native result bytes without digest update',lambda:S._check_dependencies(bundle,False),'baseline-artifact-hash');(dst/'results.json').write_bytes(raw)
# Nonblocking validator limit: printable positional plan strings are diagnostic only in original audit.
q=copy.deepcopy(native);q['cells'][0]['graphs'][0]['plans'][0]['current_plan']='SYNTHETIC_DIAGNOSTIC_CORRUPTION';probe('diagnostic integer plan string ignored by packet auditor',lambda:N.audit_typed(q))
save('fresh-heldouts.json',{'probes':probes,'count':len(probes),'rejections':sum(x['result']=='REJECT'for x in probes),'native_executions':0,'validator_limit':'Packet audit checks named plans, not debug current_plan/historical_plan strings. Independent interface_review_v2.py parses both integer diagnostic plans for all actual selected pairs/folds and found no mismatch. This is not a fixed-corpus defect or an approval of a general validator.'})
print('fresh heldouts passed',len(probes),flush=True)
