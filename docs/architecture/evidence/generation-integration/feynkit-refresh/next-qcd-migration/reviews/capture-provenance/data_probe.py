"""Independent full-domain joins and expression syntax only; no scientific imports."""
from pathlib import Path
from collections import Counter
import copy
import hashlib
import json
import re

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
P=ROOT/'target/raw-next-qcd-preflight-20261031-r1'
R=ROOT/'target/raw-next-qcd-captures-20261031-r2'
CASES={'xs_qcd_gguu_sym_L1':(5,65),'xs_qcd_gguu_L1':(9,117),'xs_qcd_gggg_L1':(10,118)}
def load(path): return json.loads(Path(path).read_text())
def need(ok, message):
    if not ok: raise ValueError(message)
def canon(value): return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False)
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def unique(rows,key,label):
    answer={key(r):r for r in rows}
    need(len(answer)==len(rows),'duplicate '+label)
    return answer
def public_owners(payload):
    rows={}
    for label,items in [('vertex',payload['vertices']),('edge',[pair[1] for pair in payload['edges']])]:
        for i,item in enumerate(items):
            rows[f'{label}:{i}']=(item['numerator'],'num_interface' in item,item.get('num_interface'))
    for field,interface in [('numerator','num_interface'),('projector','projector_interface'),('numerator_prefactor','numerator_prefactor_interface')]:
        rows[field]=(payload[field],interface in payload,payload.get(interface))
    return rows

TOKEN=re.compile(r'(?:[^\W\d]\w*::)+(?:\{[^{}]*\}::)?[^\W\d]\w*|[0-9]+𝑖|𝑖|[0-9]+|[()+*/^,\-]')
SPELLINGS={f'gammalooprs::{{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:{n.lower()}}}::{n}':f'gammalooprs::{{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}}::{n}' for n in ['Q','K','P']}
SPELLINGS.update({f'gammalooprs::{{spenso::index,spenso::index-label:{tag}}}::{name}':f'gammalooprs::{{spenso::index}}::{name}' for name,tag in [('edge','e'),('vertex','v')]})
def tokens(text,historical=False):
    result=[];cursor=0
    for m in TOKEN.finditer(text):
        need(not text[cursor:m.start()].strip(),'unsupported expression syntax')
        result.append(SPELLINGS.get(m[0],m[0]) if historical else m[0]);cursor=m.end()
    need(not text[cursor:].strip() and bool(result),'incomplete expression syntax')
    return result
def mapped_text(text,historical=False):
    tokens(text)
    return TOKEN.sub(lambda m:SPELLINGS.get(m[0],m[0]),text) if historical else text

class Expr:
    """Recursive descent, fully consumed; compares syntax, never algebra."""
    def __init__(self,text,historical=False): self.t=tokens(text,historical);self.i=0
    def peek(self): return self.t[self.i] if self.i<len(self.t) else None
    def take(self,t):
        need(self.peek()==t,'expression delimiter');self.i+=1
    def primary(self):
        t=self.peek();need(t is not None,'expression truncated');self.i+=1
        if t=='(':
            r=self.add();self.take(')');return r
        if t.isdecimal():return ('integer',t)
        if t.endswith('𝑖') and t[:-1].isdecimal():return ('*',('integer',t[:-1]),('symbol','𝑖'))
        need('::' in t or t=='𝑖','expression atom')
        if self.peek()!='(':return ('symbol',t)
        self.take('(');args=[]
        if self.peek()!=')':
            args.append(self.add())
            while self.peek()==',':self.take(',');args.append(self.add())
        self.take(')');return ('call',t,tuple(args))
    def power(self):
        lhs=self.primary()
        if self.peek()=='^':self.take('^');return ('^',lhs,self.unary())
        return lhs
    def unary(self):
        if self.peek()=='-':self.take('-');return ('neg',self.unary())
        return self.power()
    def mul(self):
        lhs=self.unary()
        while self.peek() in ['*','/']:
            op=self.peek();self.take(op);lhs=(op,lhs,self.unary())
        return lhs
    def add(self):
        lhs=self.mul()
        while self.peek() in ['+','-']:
            op=self.peek();self.take(op);lhs=(op,lhs,self.mul())
        return lhs
    def parse(self):
        r=self.add();need(self.i==len(self.t),'expression trailing tokens');return r

def definitions(text,result):
    Expr(text).parse()
    for token in tokens(text):
        if '::' not in token:continue
        name=re.sub(r'::\{[^{}]*\}::','::',token)
        need(name not in result or result[name]==token,'conflicting captured symbol tags')
        result[name]=token

def validate_world(w,expected_case,expected_era,expected_mode,table):
    case,era,mode=expected_case,expected_era,expected_mode
    count,owner_count=CASES[case]
    need((w['case'],w['route'],w['mode'])==(case,era,mode),'world identity')
    live,raw,payload=w['live'],w['raw'],w['payload']
    need(live['case']==raw['case']==case and live['mode']==raw['mode']==mode,'event identity')
    need(live['request']['classes']==mode,'requested Classes')
    need(live['grouping']=='none' and live['Run']=='Run { threads: Some(1), parallel: true, check_heredity: false, cancelled: false }','Grouping and Run')
    if era=='current':need(live['threads']==1,'pool1')
    graphs=unique(payload['diagrams'],lambda g:g['id'],'public ID')
    keys=unique(payload['diagrams'],lambda g:canon(g['key']),'full key')
    raw_graphs=unique(raw['diagrams'],lambda g:g['id'],'raw ID')
    need(len(graphs)==len(keys)==count==live['prepared_len'] and set(raw_graphs)==set(graphs),'complete graph domain')
    audit=unique(table,lambda r:(r['id'],r['owner']),'audit owner')
    stats=Counter();defs={};owner_data={}
    need(len(payload['groups'])==count,'complete groups')
    for index,g in enumerate(payload['diagrams']):
        gid=g['id'];rg=raw_graphs[gid]
        need(re.fullmatch('[0-9a-f]{32}',gid) is not None,'native exact ID')
        need(rg['name']==g['name'] and rg['json']['name']==g['name'],'native name')
        # Diagnostic serde_json::Value id is f64-rounded: NEVER use it as an ID.
        need({k:v for k,v in rg['json'].items() if k not in ['id','name']}==g['payload'],'raw JSON public payload')
        factor=g['overall_factor']
        need(factor==g['payload']['overall_factor']==rg['json']['overall_factor'],'source factor')
        group=payload['groups'][index]
        need(group=={'master':index,'members':[{'id':gid,'name':g['name'],'overall_factor':factor,'ratio':'1','source_index':index}]},'exact group ratio source factor')
        public=public_owners(g['payload']);native=unique(rg['owners'],lambda r:r['owner'],'native owner')
        need(set(native)==set(public),'complete owner domain')
        for owner,(decoded,present,interface) in public.items():
            tensor=native[owner]['tensor'];saved=audit.get((gid,owner))
            need(saved is not None and saved['full_key']==g['key'],'owner exact ID/full key')
            need(tensor['decoded']==decoded,'raw-to-public decoded owner')
            need(saved['tensor']==tensor,'verbatim tensor flags ports raw association')
            need(saved['public_interface_present']==present and ('public_interface' in saved)==present,'interface absence preservation')
            if present:need(interface==tensor['ports']==saved['public_interface'],'ordered public ports')
            stats['owners']+=1
            stats['public_interface_present']+=present;stats['public_interface_absent']+=not present
            for flag in ['is_metric','is_composite']:
                stats[flag+'_present']+=flag in tensor;stats[flag+'_absent']+=flag not in tensor
                if flag in tensor:
                    need(type(tensor[flag]) is bool,'flag boolean');stats[flag+'_true']+=tensor[flag]
            if 'ports' in tensor:stats['raw_ports']+=len(tensor['ports'])
            for text in [tensor['raw'],tensor['decoded'],*tensor.get('ports',[])]:definitions(text,defs)
            owner_data[(gid,owner)]=tensor
        definitions(factor,defs)
    need(len(audit)==stats['owners']==owner_count,'complete audit owner count')
    return graphs,owner_data,dict(stats),defs

worlds={};results=[];ledger_report=[];total_owners=0
for case,(count,_) in CASES.items():
    for era in ['historical','current']:
        pristine=P/'execution/work'/f'{case}--{era}-pristine'
        observer=P/'execution/work'/f'{case}--{era}-observer'
        census=(observer/'dumps'/f'{case}.census').read_bytes()
        need(census==(pristine/'dumps'/f'{case}.census').read_bytes(),'pristine observer census bytes')
        events=[json.loads(line) for line in (observer/'ledger.jsonl').read_text().splitlines()]
        need([e['sequence'] for e in events]==list(range(len(events))),'full ledger sequence')
        stages=['small-gauge-live','small-gauge-raw'] if era=='historical' else ['batch-live-inputs','batch-raw-diagrams','batch-raw-census']
        expected_events=2*(count+2) if era=='historical' else 6
        need(len(events)==expected_events,'complete event domain')
        ledger_report.append({'case':case,'era':era,'events':len(events),'stages':dict(Counter(e['stage'] for e in events))})
        for ordinal,mode in enumerate(['None','Auto']):
            length=count+2 if era=='historical' else 3
            segment=events[ordinal*length:(ordinal+1)*length]
            need([e['stage'] for e in segment]==([stages[0]]+['into-representative']*count+[stages[1]] if era=='historical' else stages),'source phase event order')
            live=segment[0]['value'];raw=segment[-1 if era=='historical' else 1]['value']
            name=f'{case}.classes-{mode}-t1'
            paths=[root/'dumps'/(name+suffix) for root in [pristine,observer] for suffix in ['.json','.sources.json']]
            data=[p.read_bytes() for p in paths];need(all(d==data[0] for d in data),'all pristine observer source companion bytes')
            payload=json.loads(data[0]);w=load(R/'prepared/captured'/case/(era+'-'+mode+'.json'))
            need(w['live']==live and w['raw']==raw and w['payload']==payload,'complete ledger to prepared world')
            need(w['source_companion_sha256']==hashlib.sha256(data[0]).hexdigest(),'source companion binding')
            if era=='current':need(segment[2]['value']=={'case':case,'mode':mode,'census':payload},'complete raw census event')
            else:
                pairs=[e['value'] for e in segment[1:-1]]
                need(all(p['kind']=='shared' for p in pairs),'historical ungrouped representative kind')
                native_ids={g['id'] for g in raw['diagrams']}
                need({p['comparison']['id'] for p in pairs}==native_ids and len(pairs)==len(native_ids),'whole historical intermediate domain')
            audit=R/'case-audits'/case
            graphs,owners,stats,defs=validate_world(w,case,era,mode,load(audit/(era+'-'+mode+'-owners.json')))
            need(defs==load(audit/(era+'-'+mode+'-whole-symbol-definitions.json')),'whole captured symbol table')
            declared=next(c for c in load(audit/'cells.json') if c['route']==era and c['mode']==mode)
            need(declared['counts']==stats and declared['whole_symbol_definitions']==len(defs),'exact declared counts')
            need(declared['IDs']==list(graphs) and declared['full_keys']=={i:g['key'] for i,g in graphs.items()},'declared native ID fullkeys')
            # Census exact string IDs and factors; no rounding and no ordinal-only joins.
            lines=[l for l in census.decode().splitlines() if re.match(r'D[0-9]+ key=',l)]
            need(len(lines)==count,'complete census lines')
            for line in lines:
                m=re.fullmatch(r'(D[0-9]+) key=([0-9a-f]{16}) id=([0-9a-f]{32}) S=(\d+) frame=([0-9a-f]{16}) numerator=([0-9a-f]{16}) factor=(.*) basis=(.*)',line)
                need(m is not None,'census row grammar');g=graphs[m[3]]
                need(g['name']==m[1] and str(g['aut'])==m[4] and g['overall_factor']==m[7] and g['basis']==json.loads(m[8]),'census exact identity factor basis')
            factor2=sum('::NumeratorIndependentSymmetryGrouping(2)' in g['overall_factor'] for g in graphs.values())
            need(factor2==(4 if '_sym_' in case else 0),'actual symmetry2 inventory')
            worlds[case,era,mode]=(w,graphs,owners)
            total_owners+=len(owners)
            results.append({'case':case,'era':era,'mode':mode,'diagrams':len(graphs),'counts':stats,'symbol_definitions':len(defs),'symmetry2_diagrams':factor2})

pipelines=[]
def physical(payload):
    p=copy.deepcopy(payload)
    for k in ['numerator','numerator_prefactor','projector','num_interface','numerator_prefactor_interface','projector_interface']:p.pop(k,None)
    for v in p['vertices']:v.pop('numerator');v.pop('num_interface',None)
    for _,e in p['edges']:e.pop('numerator');e.pop('num_interface',None)
    return p
for case in CASES:
    for mode in ['None','Auto']:
        hw,hg,ho=worlds[case,'historical',mode];cw,cg,co=worlds[case,'current',mode]
        need(set(hg)==set(cg) and set(ho)==set(co),'cross era complete domains')
        for key in ['model','request','options','grouping','Run','prepared_len']:need(hw['live'][key]==cw['live'][key],'cross era live '+key)
        need(hw['payload']['groups']==cw['payload']['groups'],'cross era source group ratios')
        for gid in hg:
            need({k:v for k,v in hg[gid].items() if k!='payload'}=={k:v for k,v in cg[gid].items() if k!='payload'},'cross era physical descriptor')
            need(physical(hg[gid]['payload'])==physical(cg[gid]['payload']),'cross era occurrence routing factor model')
        stats=Counter();first_raw=first_decoded=first_iface=raw_equal_decoded_changed=None
        candidates=unique(load(R/'case-audits'/case/(mode+'-owner-input-candidates.json')),lambda r:(r['id'],r['owner']),'pair input')
        need(set(candidates)==set(ho),'complete pair domain')
        for key,h in ho.items():
            c=co[key];gid,owner=key;pair=candidates[key]
            need(pair['case']==case and pair['mode']==mode and pair['full_key']==hg[gid]['key'] and pair['historical']==h and pair['current']==c,'complete owner input pair association')
            for flag in ['is_metric','is_composite']:need((flag in h)==(flag in c) and h.get(flag)==c.get(flag),'cross era flags')
            raw_equal=mapped_text(h['raw'],True)==c['raw']
            decoded_equal=Expr(h['decoded'],True).parse()==Expr(c['decoded']).parse()
            stats['owners']+=1;stats['raw_bytes_equal_after_five_source_maps']+=raw_equal;stats['decoded_syntax_equal_after_five_source_maps']+=decoded_equal
            label={'id':gid,'owner':owner}
            if not raw_equal and first_raw is None:first_raw=label
            if not decoded_equal and first_decoded is None:first_decoded=label
            if raw_equal and not decoded_equal and raw_equal_decoded_changed is None:raw_equal_decoded_changed=label
            hp=public_owners(hg[gid]['payload'])[owner];cp=public_owners(cg[gid]['payload'])[owner]
            if hp[1]!=cp[1]:
                stats['public_interface_presence_differences']+=1
                if first_iface is None:first_iface=label
        p=next(p for p in load(R/'pipeline-comparison.json') if p['case']==case and p['mode']==mode)
        need(dict(stats)==p['observations'],'pipeline counts')
        for k,value in [('first_raw_difference_after_source_maps',first_raw),('first_decoded_syntax_difference_after_source_maps',first_decoded),('first_public_interface_presence_difference',first_iface),('first_raw_equal_decoded_changed_owner',raw_equal_decoded_changed)]:need(p[k]==value,'pipeline boundary '+k)
        pipelines.append({'case':case,'mode':mode,'observations':dict(stats),'first_raw':first_raw,'first_decoded':first_decoded,'first_interface':first_iface})

controls=[]
base=worlds['xs_qcd_gguu_sym_L1','current','None'][0]
table=load(R/'case-audits/xs_qcd_gguu_sym_L1/current-None-owners.json')
mutations=[
 ('wrong-mode',lambda w:w.__setitem__('mode','Auto'),'world identity'),
 ('wrong-request-mode',lambda w:w['live']['request'].__setitem__('classes','Auto'),'requested Classes'),
 ('wrong-event-case',lambda w:w['raw'].__setitem__('case','xs_qcd_gguu_L1'),'event identity'),
 ('duplicate-public-ID',lambda w:w['payload']['diagrams'].__setitem__(1,copy.deepcopy(w['payload']['diagrams'][0])),'duplicate public ID'),
 ('duplicate-full-key',lambda w:w['payload']['diagrams'][1].__setitem__('key',w['payload']['diagrams'][0]['key']),'duplicate full key'),
 ('missing-raw-graph',lambda w:w['raw']['diagrams'].pop(),'complete graph domain'),
 ('duplicate-raw-graph',lambda w:w['raw']['diagrams'].append(copy.deepcopy(w['raw']['diagrams'][0])),'duplicate raw ID'),
 ('rounded-ID-substitution',lambda w:w['raw']['diagrams'][0].__setitem__('id',w['raw']['diagrams'][0]['json']['id']),'complete graph domain'),
 ('source-factor-change',lambda w:w['payload']['diagrams'][0].__setitem__('overall_factor','1'),'source factor'),
 ('wrong-group-ratio',lambda w:w['payload']['groups'][0]['members'][0].__setitem__('ratio','2'),'group ratio'),
 ('wrong-source-index',lambda w:w['payload']['groups'][0]['members'][0].__setitem__('source_index',1),'group ratio'),
 ('missing-owner',lambda w:w['raw']['diagrams'][0]['owners'].pop(),'complete owner domain'),
 ('duplicate-owner',lambda w:w['raw']['diagrams'][0]['owners'].append(copy.deepcopy(w['raw']['diagrams'][0]['owners'][0])),'duplicate native owner'),
 ('raw-sign-change',lambda w:w['raw']['diagrams'][0]['owners'][0]['tensor'].__setitem__('raw','-('+w['raw']['diagrams'][0]['owners'][0]['tensor']['raw']+')'),'verbatim tensor'),
 ('decoded-owner-change',lambda w:w['raw']['diagrams'][0]['owners'][0]['tensor'].__setitem__('decoded','1'),'decoded owner'),
 ('ordered-port-change',lambda w:w['raw']['diagrams'][0]['owners'][0]['tensor']['ports'].reverse(),'verbatim tensor'),
 ('flag-omission',lambda w:w['raw']['diagrams'][0]['owners'][0]['tensor'].pop('is_metric'),'verbatim tensor'),
 ('flag-change',lambda w:w['raw']['diagrams'][0]['owners'][0]['tensor'].__setitem__('is_composite',False),'verbatim tensor'),
]
for name,mutate,predicate in mutations:
    w=copy.deepcopy(base);mutate(w)
    try:validate_world(w,'xs_qcd_gguu_sym_L1','current','None',table)
    except ValueError as e:
        need(predicate in str(e),'wrong negative predicate '+name)
        controls.append({'name':name,'rejected_by':str(e),'credit':'pure capture association only'})
    else:raise AssertionError('accepted negative '+name)

report={'status':'PASS','public_mode_world_copies':24,'complete_observer_worlds':12,'case_mode_cells':6,
        'complete_owner_records_checked':total_owners,'complete_pairs_checked':sum(v[1] for v in CASES.values())*2,
        'worlds':results,'ledgers':ledger_report,'pipeline':pipelines,'negative_controls':controls,
        'negative_control_count':len(controls),'native_invocations':0,'library_imports':0,
        'scope':'Exact bytes, native string IDs/full keys, request/mode/phase, raw/decoded/ports/flags/factors and syntax only. No physical, algebra, locality or interface acceptance.'}
(OUT/'data-probe.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ['worlds','pipeline','negative_controls','ledgers']}))
