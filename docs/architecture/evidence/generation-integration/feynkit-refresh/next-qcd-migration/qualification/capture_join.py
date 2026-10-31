"""Full fresh public bytes and retained exact proof associations. No scientific execution."""
import ast,collections,copy,hashlib,importlib.util,json,re,sys,types
from fractions import Fraction
from pathlib import Path
O=Path(__file__).resolve().parent;W=O.parents[2];F=W.parent/'generation-oracle-numerator-families'
P=F/'target/raw-next-qcd-migration-preparation-20261031-r1';CAP=F/'target/raw-next-qcd-captures-20261031-r2';PRE=F/'target/raw-next-qcd-preflight-20261031-r1'
R1=F/'target/raw-next-qcd-proof-preparation-20261031-r1';R3=F/'target/raw-next-qcd-proof-preparation-20261031-r3';N=F/'target/raw-next-qcd-proof-native-execution-20261031-r1'
CASES=('xs_qcd_gguu_sym_L1','xs_qcd_gguu_L1','xs_qcd_gggg_L1');MODES=('None','Auto');used={}
def need(x,*d):
 if not x:raise AssertionError(d)
def sha(data):return hashlib.sha256(data).hexdigest()
def read(p):
 p=Path(p);b=p.read_bytes();used[str(p)]={'sha256':sha(b),'bytes':len(b)};return b
def unique(pairs):
 d={}
 for k,v in pairs:need(k not in d,'uniqueJSON',k);d[k]=v
 return d
def load(p):return json.loads(read(p),object_pairs_hook=unique)
def ref(p):read(p);return {'path':str(p),**used[str(p)]}
def save(n,d):
 p=O/'joins'/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
C=types.SimpleNamespace(sha=sha,unique=unique,need=need,load=load,read=read,ref=ref,save=save,used=used,R1=R1,R3=R3,N=N,CAP=CAP,PRE=PRE,CASES=CASES,MODES=MODES)
spec=importlib.util.spec_from_file_location('approved_pure_parser',P/'scripts/pure_parser.py');B=importlib.util.module_from_spec(spec);spec.loader.exec_module(B)
ns={'collections':collections,'copy':copy,'Fraction':Fraction,'json':json,'re':re,'C':C,'B':B}
tree=ast.parse(read(P/'reconstruct.py'));defs=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in {'Digest','multiset','factor','digests','derive'}];exec(compile(ast.Module(body=defs,type_ignores=[]),str(P/'reconstruct.py'),'exec'),ns)
ns['FIELDS']=('diagrams','topologies','assignments','inverse_symmetry','weights','symmetry_factors','content','numerators','named','ids')
selfcheck='--selfcheck' in sys.argv
if not selfcheck:
 runtime=load(O/'runtime-results.json');need(runtime['passes']==runtime['actual_executions']==48 and runtime['failures']==0,'48actualpasses')
files=[];cells=[];fresh_worlds={}
for case,owners in zip(CASES,(65,117,118),strict=True):
 oldroot=PRE/f'execution/work/{case}--current-pristine/dumps'
 fresh=oldroot if selfcheck else O/f'work/pool1-{case}/dumps'
 expected={case+'.census'}|{case+f'.classes-{m}-t1'+s for m in MODES for s in ('.json','.sources.json')}
 actual={p.name for p in fresh.iterdir() if p.is_file()};need(actual==expected,'complete-exact-capture-file-domain',case,actual)
 for name in sorted(actual):
  b=read(fresh/name);need(b==read(oldroot/name),'smallest-boundary-full-byte-mismatch',case,name)
  files.append({'fresh':str(fresh/name),'approved_pristine':str(oldroot/name),'sha256':sha(b),'bytes':len(b),'byte_equal':True,'normalizations':[]})
 for mode in MODES:
  path=CAP/f'prepared/captured/{case}/current-{mode}.json';world=load(path);need((world['case'],world['mode'],world['route'])==(case,mode,'current'),'worldidentity')
  live=world['live'];need(live['case']==case and live['mode']==mode and live['request']['classes']==mode and live['grouping']=='none','case-mode-request')
  need('threads: Some(1)' in live['Run'],'retainedpool1Run')
  payload=load(fresh/(case+f'.classes-{mode}-t1.json'));source=load(fresh/(case+f'.classes-{mode}-t1.sources.json'))
  need(payload==source==world['payload'],'whole-payload-source-world')
  # Bind ID as its exact128-bit string plus entire canonical key, never ordinal or diagnostic f64.
  keyed={(g['id'],json.dumps(g['key'],sort_keys=True,separators=(',',':'))):g for g in payload['diagrams']}
  need(len(keyed)==len(payload['diagrams']),'unique-id-key');count=0;explicit=0;owner_records=[]
  raw={g['id']:g for g in world['raw']['diagrams']}
  for g in payload['diagrams']:
   need(re.fullmatch('[0-9a-f]{32}',g['id']) is not None,'128bitstring')
   public=B.owners(g['payload']);r=raw[g['id']];need(g['key']==r['key'] if 'key' in r else True,'rawfullkey')
   need(len(public)==len(r['owners']),'whole-owner-domain')
   for o in r['owners']:
    name,t=o['owner'],o['tensor'];expression,interface=public[name]
    need(expression==t['raw'],'raw-public-owner-expression')
    if interface is not None:need(interface==t['ports'],'ordered-explicit-interface');explicit+=1
    count+=1;owner_records.append({'id':g['id'],'full_key_sha256':sha(json.dumps(g['key'],sort_keys=True,separators=(',',':')).encode()),'owner':name,'expression_sha256':sha(expression.encode()),'interface_present':interface is not None,'ordered_interface':interface,'retained_decoded':t['decoded'],'retained_ports':t['ports'],'retained_flags':{'is_metric':t['is_metric'],'is_composite':t['is_composite']},'association':'Fresh full public bytes and source invariance attach the prior captured raw/cooked/decoded/flags; not a fresh observer or universal decoder claim'})
  need(count==owners,'exactcellowner-count')
  census=read(fresh/(case+'.census')).decode();report=census.splitlines()[0];need(report==live['report'],'completereport')
  row,details=ns['derive'](payload,report);expectedrow=load(P/f'rows/{case}/current-{mode}.json')['Row'];need(row==expectedrow,'reconstructed-all10fields')
  factor2=[g['id'] for g in payload['diagrams'] if 'feynkit_generator_factor::{}::NumeratorIndependentSymmetryGrouping(2)' in g['overall_factor']]
  if case==CASES[0]:need(len(factor2)==4,'four-of-five factor2')
  cells.append({'case':case,'mode':mode,'owners':count,'explicit_interfaces':explicit,'graphs':len(payload['diagrams']),'full_Row_reconstruction':row,'Row_asserted_by_native':mode=='None','Auto_report_counter_scope':'Fresh public output and success; prior mode live/report/search counters remain separately source-bound. No new Auto literal assertion.','live_inputs':live,'fresh_public':ref(fresh/(case+f'.classes-{mode}-t1.json')),'source_companion':ref(fresh/(case+f'.classes-{mode}-t1.sources.json')),'retained_world':ref(path),'owners_chain':owner_records,'factor2_ids':factor2,'per_graph_digests':details})
  fresh_worlds[str(path)]=dict(world,payload=payload)
need(sum(c['owners'] for c in cells)==600,'600 fresh current owner joins')
# Replay ONLY pure association definitions from the reviewed preparatory script, writing here.
original_load=C.load
def freshload(path):return copy.deepcopy(fresh_worlds[str(path)]) if str(path) in fresh_worlds else original_load(path)
C.load=freshload
ans={'collections':collections,'copy':copy,'json':json,'re':re,'C':C,'R':types.SimpleNamespace(B=B),'B':B,'cache':{}}
atree=ast.parse(read(P/'associate.py'));defs=[n for n in atree.body if isinstance(n,ast.FunctionDef)]
exec(compile(ast.Module(body=defs,type_ignores=[]),str(P/'associate.py'),'exec'),ans);ans['main']()
association=load(O/'joins/proof-association.json');association['status']='EXACT_FRESH_PUBLIC_SOURCE_AND_PRIOR_PROOF_ASSOCIATIONS_VERIFIED';association['approval']=False;association['new_native_or_math_credit']=0;association['scientific_review_reused']='6d892bd9263a7c093bba5fd372fb4663a5a4e81dfaad4bfc709bbba052b96d21';association['general_audit_native_approval']=False;save('proof-association.json',association)
review=F/'.local/review-next-qcd-scientific-interface-astra-20261031-5a3KbP5x'
audit=load(review/'interface-independent.json');need(audit['counts']=={'owners':600,'changed_lowerings':210,'unchanged_source_arguments':390,'selected_pairs':1890,'factor_products':96,'visited_sets':41088,'fold_endpoint_pairs':720,'fold_steps':456,'graphs':48},'accepted-exact-plan-layout-audits')
summary={'schema':'next-qcd-fresh-public-capture-and-proof-joins/1','phase':'PURE_SELFCHECK_NO_EXECUTION_CREDIT' if selfcheck else 'ACTUAL_POSTAPPLICATION_ARCHIVE','files':files,'actual_file_count':len(files),'cells':cells,'current_owner_joins':600,'historical_current_relations':600,'independent_integer_plan_layout_audit':ref(review/'interface-independent.json'),'integer_plan_layout_audit_source':ref(review/'interface.py') if (review/'interface.py').exists() else None,'native_actual_results':ref(R3/'execution/records/all-six-cells--typed8eb/results.json'),'source_invariance':ref(O/'immediate-producer-source-delta.json') if not selfcheck else 'pending actual qualified producer','raw_cooked_decoded_flags':'Retained accepted capture/observer source premises, not newly observed flags','both_operand_and_output_axis_maps_required':True,'historical_interface_absence_is_not_empty':True,'normalizations':[],'new_native_proof_runs':0,'new_observer':False,'general_audit_native_approved':False}
if not selfcheck:summary['actual_archive']=load(O/'archive-qualification.json')['archive_sha256']
save('capture-joins-selfcheck.json' if selfcheck else 'capture-joins.json',summary);save('consumed-input-bindings.json',used)
print('PASS',summary['phase'],len(files),'exact complete files;6 cells/600 fresh current owners;600 retained paired relations; both axis maps and all actual proof bindings; no native proof rerun')
