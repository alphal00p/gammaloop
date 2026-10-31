import sys,re,ast,copy,collections,struct
from fractions import Fraction
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import common as C

class Rejected(ValueError):pass
def require(ok,label):
 if not ok:raise Rejected(label)

def encoded(x):
 # serde_json::Value uses sorted object keys in this producer; exact integer domain only here.
 def domain(v):
  require(v is None or type(v) in (bool,int,str,list,dict),'JSON-domain')
  if isinstance(v,dict):
   require(all(type(k)is str for k in v),'JSON-key')
   for val in v.values():domain(val)
  if isinstance(v,list):
   for val in v:domain(val)
 domain(x)
 return C.json.dumps(x,sort_keys=True,ensure_ascii=False,separators=(',',':'),allow_nan=False).encode()
def fnv(data):
 result=0xcbf29ce484222325
 for byte in data:result=((result^byte)*1099511628211)&0xffffffffffffffff
 return result
def u64(n):return struct.pack('<Q',n)
def text(x):return x.encode()+b'\xff'
def name(x):return b'\xfe' if x is None else text(x)
def multiset(values):return fnv(b''.join(u64(v) for v in sorted(values)))

def weights(expr):
 allowed={'AutG','CouplingsMultiplicity','InternalFermionLoopSign','ExternalFermionOrderingSign','AntiFermionSpinSumSign','NumeratorIndependentSymmetryGrouping'}
 compact=re.sub(r'feynkit_generator_factor::\{\}::([A-Za-z]+)',r'\1',expr).replace('^','**')
 tree=ast.parse(compact,mode='eval')
 def calc(n):
  if isinstance(n,ast.Constant) and type(n.value)is int:return Fraction(n.value)
  if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -calc(n.operand)
  if isinstance(n,ast.Call):
   require(isinstance(n.func,ast.Name) and n.func.id in allowed and len(n.args)==1 and not n.keywords,'bookkeeping-source-head')
   return calc(n.args[0])
  require(isinstance(n,ast.BinOp),'factor-syntax');a,b=calc(n.left),calc(n.right)
  if isinstance(n.op,ast.Add):return a+b
  if isinstance(n.op,ast.Sub):return a-b
  if isinstance(n.op,ast.Mult):return a*b
  if isinstance(n.op,ast.Div):return a/b
  require(isinstance(n.op,ast.Pow) and b.denominator==1,'factor-power');return a**int(b)
 return calc(tree.body)

def owners(payload):
 d={k:(payload[k],payload.get(f)) for k,f in [('numerator','num_interface'),('numerator_prefactor','numerator_prefactor_interface'),('projector','projector_interface')]}
 d.update({f'vertex:{i}':(v['numerator'],v.get('num_interface')) for i,v in enumerate(payload['vertices'])})
 d.update({f'edge:{i}':(e['numerator'],e.get('num_interface')) for i,(_,e) in enumerate(payload['edges'])})
 return d

def frame_numerator(payload):
 p=copy.deepcopy(payload);parts=[]
 for k in ('numerator','numerator_prefactor','projector'):parts.append((k,encoded(p.pop(k))))
 for k in ('num_interface','numerator_prefactor_interface','projector_interface'):
  if k in p:parts.append((k,encoded(p.pop(k))))
 for label,v in [(f'vertex:{i}',v) for i,v in enumerate(p['vertices'])]+[(f'edge:{i}',v) for i,(_,v) in enumerate(p['edges'])]:
  parts.append((label,encoded(v['numerator'])));v['numerator']=None
  if 'num_interface' in v:parts.append((label+'/num_interface',encoded(v.pop('num_interface'))))
 for k in ('model','symmetry_factor','overall_factor'):del p[k]
 stream=b''.join(b for _,b in parts)
 return fnv(encoded(p)),fnv(stream),[{'field':k,'sha256':C.sha(b),'bytes':len(b)} for k,b in parts]

def parse_rows(source):
 lines={};fields={}
 row_pattern=re.compile(r'    (fast|slow|heavy) (\w+): Row \{ (.*) \};')
 for line in source.splitlines(keepends=True):
  match=row_pattern.fullmatch(line.rstrip('\n'))
  if not match:continue
  tier,case,body=match.groups();require(case not in lines,'duplicate-row')
  # Translate only identifiers/borrowed slice markers; literal_eval cannot execute code.
  obj=ast.literal_eval('{'+re.sub(r'\b(\w+):',r'"\1":',body).replace('&[','[')+'}')
  require(set(obj)=={'diagrams','topologies','assignments','inverse_symmetry','weights','symmetry_factors','content','numerators','named','ids'},'Row-ten-fields')
  fields[case]=C.json.loads(C.json.dumps(obj));lines[case]=line
 require(len(lines)==149,'Row-domain')
 return fields,lines

def derive(payload,report):
 m=re.fullmatch(r'GenerationReport \{ topology_count: (\d+), interaction_assignment_count: (\d+), retained_count: (\d+), zero_numerator_count: (\d+), completed: true \}',report)
 require(m is not None,'complete-report');topology,assignments,retained,zero=map(int,m.groups())
 require(retained==len(payload['diagrams']),'retained-domain')
 content=[];numerator=[];ids=[];named=b'';hist=collections.Counter();inv=Fraction();weight=Fraction();details=[]
 for d in payload['diagrams']:
  require(re.fullmatch('[0-9a-f]{32}',d['id']) is not None,'exact-ID')
  key=fnv(encoded(d['key']));frame,num,segments=frame_numerator(d['payload']);require(type(d['aut'])is str and re.fullmatch(r'[1-9][0-9]*',d['aut']) is not None,'decimal-aut');aut=int(d['aut']);require(aut==d['payload']['symmetry_factor'],'aut')
  content.append(fnv(u64(key)+u64(aut)+text(d['overall_factor'])+u64(frame)))
  numerator.append(fnv(u64(key)+u64(num)))
  ids.append(fnv(int(d['id'],16).to_bytes(16,'little')))
  named+=name(d['name'])+u64(key);hist[aut]+=1;inv+=Fraction(1,aut);weight+=weights(d['overall_factor'])
  details.append({'id':d['id'],'full_key':d['key'],'name':d['name'],'key':f'{key:016x}','frame':f'{frame:016x}','numerator':f'{num:016x}','numerator_byte_segments':segments})
 require(len({d['id'] for d in payload['diagrams']})==retained,'unique-ID')
 require(len({encoded(d['key']) for d in payload['diagrams']})==retained,'unique-full-key')
 for group in payload['groups']:
  named+=u64(group['master'])
  for member in group['members']:named+=name(member['name'])+text(member['ratio'])+text(member['overall_factor'])
 result={'diagrams':retained,'topologies':topology,'assignments':assignments,'inverse_symmetry':[inv.numerator,inv.denominator],'weights':[weight.numerator,weight.denominator],'symmetry_factors':[list(r) for r in sorted(hist.items())],'content':multiset(content),'numerators':multiset(numerator),'named':fnv(named),'ids':multiset(ids)}
 return result,details

def validate_world(w,case,mode,era):
 require((w['case'],w['mode'],w['route'])==(case,mode,era),'world-identity')
 live=w['live'];require(live['case']==case,'live-case');require(live['mode']==mode and live['request']['classes']==mode,'live-mode')
 require(live['grouping']=='none' and 'threads: Some(1)' in live['Run'] and ('threads' not in live or live['threads']==1),'run-domain')
 p=w['payload'];raw={g['id']:g for g in w['raw']['diagrams']}
 require(w['raw']['case']==case and w['raw']['mode']==mode,'raw-identity')
 require(len(raw)==len(w['raw']['diagrams'])==len(p['diagrams'])==live['prepared_len'],'raw-domain')
 for i,d in enumerate(p['diagrams']):
  r=raw[d['id']];require(d['name']==r['name'],'raw-name')
  expected=owners(d['payload']);ro={a['owner']:a['tensor'] for a in r['owners']}
  require(len(ro)==len(r['owners']) and ro.keys()==expected.keys(),'owner-domain')
  for label,(expression,interface) in expected.items():
   tensor=ro[label];require(tensor['decoded']==expression,'cooked-decoded-public')
   if interface is not None:require(tensor['ports']==interface,'public-ports')
  require(p['groups'][i]=={'master':i,'members':[{'id':d['id'],'name':d['name'],'overall_factor':d['overall_factor'],'ratio':'1','source_index':i}]},'exact-membership')
  require(r['json']['overall_factor']==d['payload']['overall_factor']==d['overall_factor'],'source-factor')
 return raw

def main():
 before=C.read(C.P/'before.rs');old_rows,lines=parse_rows(before.decode());results=[];worlds={}
 for case in C.CASES:
  for mode in C.MODES:
   for era in ('historical','current'):
    wp=C.CAP/f'prepared/captured/{case}/{era}-{mode}.json';w=C.load(wp);validate_world(w,case,mode,era)
    roots=[C.T/'raw-next-qcd-preflight-20261031-r1/execution/work'/f'{case}--{era}-{kind}' for kind in ('pristine','observer')]
    actual=[]
    for root in roots:
     pp=root/'dumps'/f'{case}.classes-{mode}-t1.json';sp=pp.with_name(pp.stem+'.sources.json')
     require(C.read(pp)==C.read(sp),'payload-companion-bytes');p=C.load(pp);require(p==w['payload'],'exact-pristine-world')
     actual.append(C.read(pp));census=C.read(root/'dumps'/f'{case}.census').decode();require(census.splitlines()[0]==w['live']['report'],'source-report')
    require(actual[0]==actual[1],'pristine-observer-byte-identity')
    events=[C.json.loads(l,object_pairs_hook=C.unique) for l in C.read(roots[1]/'ledger.jsonl').decode().splitlines()]
    require([e['sequence'] for e in events]==list(range(len(events))),'complete-ledger')
    stages=('small-gauge-live','small-gauge-raw') if era=='historical' else ('batch-live-inputs','batch-raw-diagrams','batch-raw-census')
    seq=[]
    for stage in stages:
     es=[e for e in events if e['stage']==stage and e['value']['case']==case and e['value']['mode']==mode];require(len(es)==1,'mode-event')
     e=es[0];seq.append(e['sequence']);want=w['live'] if stage==stages[0] else w['raw'] if stage==stages[1] else None
     require(e['value']==want if want is not None else e['value']['census']==w['payload'],'actual-event-binding')
    require(seq==sorted(seq),'event-order')
    row,detail=derive(w['payload'],w['live']['report'])
    if mode=='None':
     for d in detail:require(f"{d['name']} key={d['key']} id={d['id']}" in census and f"frame={d['frame']} numerator={d['numerator']}" in census,'pristine-census-crosscheck')
     if era=='historical':require(row==old_rows[case],'historical-complete-Row')
    declared=C.load(C.P/f'rows/{case}/{era}-{mode}.json');require(row==declared['Row'],'candidate-complete-Row')
    worlds[case,mode,era]=w
    results.append({'case':case,'mode':mode,'era':era,'Row':row,'row_literal_assertion':mode=='None','diagrams':detail,'world':str(wp),'world_sha256':C.sha(C.read(wp)),'actual_event_sequence':seq,'historical_live_threads_omitted':'threads' not in w['live']})
 for case in C.CASES:
  for mode in C.MODES:
   h,c=(worlds[case,mode,e] for e in ('historical','current'))
   require(all(h['live'][k]==c['live'][k] for k in ('model','request','options','Run','grouping','prepared_len','report')),'full-live-physical-inputs')
   hr,cr=(next(r['Row'] for r in results if (r['case'],r['mode'],r['era'])==(case,mode,e)) for e in ('historical','current'))
   require([k for k in hr if hr[k]!=cr[k]]==['numerators'],'sole-numerator-difference')
 C.save('all12-Rows-independent.json',results);C.flush('rows-consumed.json')
 print('PASS12 full Rows derived independently from exact pristine serialized bytes; all ten fields match candidate. OnlyNone literal scope.')
 for r in results:
  if r['mode']=='None':print(r['case'],r['era'],hex(r['Row']['numerators']))
if __name__=='__main__':main()
