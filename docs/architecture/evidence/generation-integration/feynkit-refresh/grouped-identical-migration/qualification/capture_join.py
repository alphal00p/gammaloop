"""Fresh public payload/source/census joins only; no native imports, normalization or proof rerun."""
import ast,copy,hashlib,json,re,sys
from collections import Counter
from fractions import Fraction
from pathlib import Path
O=Path(__file__).resolve().parent;W=O.parent.parent.parent
F=W.parent/'generation-oracle-numerator-families';P=F/'target/grouped-identical-migration-preparation-20261031-r1'
G=W/'.local/feynkit-refresh/grouped-none-qualification-20261030T191919'
CASE='xs_qed_eemm_identical_L2'
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,d):(O/n).write_text(json.dumps(d,indent=2)+'\n')
def need(x,s):
 if not x:raise AssertionError(s)
def one(r,s):
 a=[x['value'] for x in r if x['stage']==s];need(len(a)==1,s);return a[0]
# Only reviewed pure algorithm definitions; no fixed-path top-level script execution.
source=P/'reconstruct.py';tree=ast.parse(source.read_bytes());allowed={'need','sha','encode','fnv','number','text','multiset','Scalar','weight','payload_digests','census'}
ns=globals().copy();defs=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in allowed]
exec(compile(ast.Module(body=defs,type_ignores=[]),str(source),'exec'),ns)
census,weight,encode=ns['census'],ns['weight'],ns['encode']

def main(dumps,selfcheck=False):
 pristine=G/'pristine-runtime/dumps';files=[]
 for suffix in ['census','classes-None-t1.json','classes-None-t1.sources.json','classes-Auto-t1.json','classes-Auto-t1.sources.json']:
  name=CASE+'.'+suffix;fresh=dumps/name;old=pristine/name
  need(fresh.read_bytes()==old.read_bytes(),'exact complete public bytes '+name)
  files.append({'file':name,'fresh':str(fresh),'approved_pristine':str(old),'bytes':fresh.stat().st_size,'sha256':sha(fresh),'byte_equal':True,'normalizations':[]})
 cells=[];source_conditions=[]
 for mode,ledger in [('None',G/'observer-runtime/captures/ledger.jsonl'),('Auto',F/'target/grouped-auto-fresh8eb-20261031-r1/observer-runtime/output/captures/ledger.jsonl')]:
  r=[json.loads(x) for x in ledger.read_bytes().splitlines()];live=one(r,'live-inputs')
  final=load(dumps/(CASE+f'.classes-{mode}-t1.json'));raw=load(dumps/(CASE+f'.classes-{mode}-t1.sources.json'))
  need(final==one(r,'final-census') and raw==one(r,'raw-census'),'complete separate mode Census')
  need(len(final['diagrams'])==29 and len(raw['diagrams'])==30,'29 finals30 sources')
  need(live['request']['classes']==mode,'mode input')
  need(live['grouping']=={'identical':{'check_canonical_numerator':False,'differentiate_particle_masses_only':True,'fully_numerical_substitution':False,'number_of_numerical_samples':5,'numerical_sample_seed':3,'symmetric_polarizations':False}},'complete grouping options')
  need(live['Run']=='Run { threads: Some(1), parallel: true, check_heredity: false, cancelled: false }','exact Run')
  need(live['options']['naming']=={'named':{'prefix':'D'}},'NamingD')
  source_conditions.append({'mode':mode,'ledger':str(ledger),'ledger_sha256':sha(ledger),'model_sha256':hashlib.sha256(encode(live['model'])).hexdigest(),'request':live['request'],'options':live['options'],'grouping':live['grouping'],'Run':live['Run'],'fresh_observation':'Public hook supplies no model/options/Run record; actual new source/filter/ELF invariance and identical mode payload/source bytes bind these accepted input premises. No fresh observer execution.'})
  owner_count=0;declared=0;records=[]
  for phase,data,graphs in [('raw',raw,one(r,'raw-diagrams')),('final',final,one(r,'final-diagrams'))]:
   need(len(data['diagrams'])==len(graphs),'all physical graphs')
   for d,g in zip(data['diagrams'],graphs,strict=True):
    need(d['payload']=={k:v for k,v in g['json'].items() if k not in ['name','id']},'whole physical payload')
    need(d['id']==g['id'] and d['name']==g['json']['name'],'full exact IDs names')
    need(d['aut']==str(d['payload']['symmetry_factor']) or d['aut']==d['payload']['symmetry_factor'],'Aut binding')
    need(d['overall_factor']==d['payload']['overall_factor'],'complete source factor')
    p=d['payload'];owners=[]
    for o in g['owners']:
     n=o['owner'];t=o['tensor']
     if ':' in n:
      kind,idx=n.split(':');obj=p['vertices'][int(idx)] if kind=='vertex' else p['edges'][int(idx)][1];field,ik='numerator','num_interface'
     else:obj=p;field=n;ik={'numerator':'num_interface','numerator_prefactor':'numerator_prefactor_interface','projector':'projector_interface'}[n]
     need(obj[field]==t['raw']==t['decoded'],'retained exact raw decoded association')
     if ik in obj:need(obj[ik]==t['ports'],'actual ordered declaration');declared+=1
     owner_count+=1;owners.append({'owner':n,'expression_sha256':hashlib.sha256(obj[field].encode()).hexdigest(),'ordered_ports':t['ports'],'is_metric':t['is_metric'],'is_composite':t['is_composite'],'flags_and_undeclared_ports_origin':'Accepted captured P3/P4 association; public hook does not independently observe flags or implicit layout'})
    records.append({'phase':phase,'id':d['id'],'name':d['name'],'full_key_sha256':hashlib.sha256(encode(d['key'])).hexdigest(),'payload_sha256':hashlib.sha256(encode(p)).hexdigest(),'owners':owners})
  seen=[];group_weights=[]
  for group,d in zip(final['groups'],final['diagrams'],strict=True):
   total=Fraction()
   for m in group['members']:
    idx=m['source_index'];s=raw['diagrams'][idx];need(m['id']==s['id'] and m['name']==s['name'] and m['overall_factor']==s['overall_factor'] and m['ratio']=='1','complete member source binding')
    total+=weight(m['overall_factor'])*Fraction(m['ratio']);seen.append(idx)
   need(total==weight(d['overall_factor']),'per group factorized scalar weight')
   group_weights.append(str(total))
  need(sorted(seen)==list(range(30)),'every source once')
  need(Counter(group_weights)==Counter({'1':17,'-1':11,'-2':1}),'exact grouped weights')
  need(sum(weight(d['overall_factor']) for d in raw['diagrams'])==sum(map(Fraction,group_weights))==4,'exact total follows membership/weights')
  report=(dumps/(CASE+'.census')).read_text().splitlines()[0] if mode=='None' else None
  reconstructed=None
  if report:
   reconstructed,per=census(final,report);expected=load(P/'reconstructed-None-row.json')['actual'];need(reconstructed==expected,'exact reviewed complete None Row');need(reconstructed['numerators']=='0xa60df38a6b8079e3','approved literal')
  cells.append({'mode':mode,'raw_graphs':30,'final_graphs':29,'owner_associations':owner_count,'explicit_declarations_raw_plus_final':declared,'physical_owners':records,'member_sources':seen,'group_weights':group_weights,'all_physical_keys_models_options_Run_factors_interfaces_bound':True,'Row':reconstructed,'report_origin':'fresh returned None census' if report and not selfcheck else 'approved original None census, selfcheck only' if report else 'No fresh Auto report or Auto literal/counter claim; original capture premise retained','new_native_layout_or_algebra_runs':0})
 need(len(cells)==2 and sum(c['owner_associations'] for c in cells)==2124,'exact fresh public association domain')
 result={'schema':'grouped-identical-public-capture-joins/1','phase':'PURE_SETUP_SELFCHECK_ONLY' if selfcheck else 'ACTUAL_NEW_ARCHIVE_RUNTIME','files':files,'cells':cells,'input_source_conditions':source_conditions,'algorithm_source':{'path':str(source),'sha256':sha(source),'definitions_only_AST':True},'byte_normalizations':[],'fresh_math':False,'source_companions_are_raw30_not_final29':True,'P1_P5_all_retained':True,'fresh_flags_or_native_layout_observation':False,'fresh_Auto_search_or_literal_assertion':False}
 if not selfcheck:
  q=load(O/'archive-qualification.json');run=load(O/'runtime-results.json');need(run['passes']==39 and run['failures']==0,'39 passes before joins');result.update(code=q['code'],archive_sha256=q['archive_sha256'],oracle_elf=q['oracle_elf'])
 save('capture-setup-selfcheck.json' if selfcheck else 'capture-joins.json',result)
 print('PASS',result['phase'],len(files),'exact files,2 mode cells,2124 retained public owner associations; Auto counters not observed')
if __name__=='__main__':main(G/'pristine-runtime/dumps' if '--selfcheck' in sys.argv else O/'work/pool1-grouped/dumps','--selfcheck' in sys.argv)
