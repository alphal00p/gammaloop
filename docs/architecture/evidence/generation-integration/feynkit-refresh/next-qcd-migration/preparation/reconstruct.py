"""Reconstruct all ten fields from pristine payloads using the immutable census encoding."""
import collections
import copy
import difflib
from fractions import Fraction
import json
import re
import sys
import context as C
sys.path.insert(0,str(C.O/'scripts'))
import pure_parser as B
import strict_owners as G

FIELDS=('diagrams','topologies','assignments','inverse_symmetry','weights','symmetry_factors','content','numerators','named','ids')
OBSERVED=(0xe4b209849953f1fa,0x665e13a4c88da679,0x96e5abcebda25959)
OLD=(0xb9fdc23ed056fa0c,0x5edf8e5f7a11eeb3,0x2628c62dc373ac91)

class Digest:
    def __init__(self): self.value=0xcbf29ce484222325
    def bytes(self,b):
        for x in b: self.value=((self.value^x)*0x100000001b3)&((1<<64)-1)
        return self
    def number(self,n): return self.bytes(n.to_bytes(8,'little'))
    def text(self,s): return self.bytes(s.encode()).bytes(b'\xff')
    def name(self,s): return self.bytes(b'\xfe') if s is None else self.text(s)
    def json(self,v): return self.bytes(json.dumps(v,sort_keys=True,ensure_ascii=False,separators=(',',':'),allow_nan=False).encode())

def multiset(items):
    d=Digest()
    for x in sorted(items): d.number(x)
    return d.value

def factor(t):
    kind=t[0]
    if kind=='integer': return Fraction(int(t[1]))
    if kind=='neg': return -factor(t[1])
    if kind=='call':
        allowed={'feynkit_generator_factor::{}::'+x for x in ('AutG','CouplingsMultiplicity','InternalFermionLoopSign','ExternalFermionOrderingSign','AntiFermionSpinSumSign','NumeratorIndependentSymmetryGrouping')}
        C.need(t[1] in allowed and len(t[2])==1,'bookkeeping-factor',t)
        return factor(t[2][0])
    a,b=factor(t[1]),factor(t[2])
    if kind=='*': return a*b
    if kind=='/': return a/b
    if kind=='+': return a+b
    if kind=='-': return a-b
    C.need(kind=='^' and b.denominator==1,'bookkeeping-power',t)
    return a**int(b)

def digests(payload):
    p=copy.deepcopy(payload);n=Digest()
    for k in ('numerator','numerator_prefactor','projector'): n.json(p.pop(k))
    for k in ('num_interface','numerator_prefactor_interface','projector_interface'):
        if k in p: n.json(p.pop(k))
    for owner in [*p['vertices'],*[e[1] for e in p['edges']]]:
        n.json(owner['numerator']);owner['numerator']=None
        if 'num_interface' in owner: n.json(owner.pop('num_interface'))
    for k in ('model','symmetry_factor','overall_factor'): p.pop(k)
    return Digest().json(p).value,n.value

def derive(payload,report):
    C.need(re.fullmatch(r'GenerationReport \{ topology_count: \d+, interaction_assignment_count: \d+, retained_count: \d+, zero_numerator_count: \d+, completed: true \}',report) is not None,'exact-complete-report')
    count=lambda k:int(re.search(r'\b'+k+r': (\d+)',report)[1])
    C.need(count('retained_count')==len(payload['diagrams']),'report-count')
    content=[];num=[];ids=[];detail=[];hist=collections.Counter();inv=Fraction();weight=Fraction();named=Digest()
    for g in payload['diagrams']:
        key=Digest().json(g['key']).value;frame,n=digests(g['payload']);s=int(g['aut'])
        C.need(s>0 and s==g['payload']['symmetry_factor'],'symmetry-factor')
        hist[s]+=1;inv+=Fraction(1,s);weight+=factor(B.parse(g['overall_factor']))
        content.append(Digest().number(key).number(s).text(g['overall_factor']).number(frame).value)
        num.append(Digest().number(key).number(n).value)
        C.need(re.fullmatch(r'[0-9a-f]{32}',g['id']) is not None,'full128-id')
        ids.append(Digest().bytes(int(g['id'],16).to_bytes(16,'little')).value)
        named.name(g['name']).number(key)
        detail.append({'id':g['id'],'full_key':g['key'],'name':g['name'],'key_digest':f'{key:016x}','frame_digest':f'{frame:016x}','numerator_digest':f'{n:016x}','overall_factor':g['overall_factor']})
    for group in payload['groups']:
        named.number(group['master'])
        for member in group['members']: named.name(member['name']).text(member['ratio']).text(member['overall_factor'])
    return dict(zip(FIELDS,[len(detail),count('topology_count'),count('interaction_assignment_count'),[inv.numerator,inv.denominator],[weight.numerator,weight.denominator],[list(v) for v in sorted(hist.items())],multiset(content),multiset(num),named.value,multiset(ids)])),detail

def rows(source):
    result={}
    for m in re.finditer(r'^    (?:fast|slow|heavy) (\w+): Row \{[^\n]+\};$',source,re.M):
        name=m[1];line=m[0];C.need(name not in result,'unique-row',name)
        fields={}
        for k in ('diagrams','topologies','assignments','content','numerators','named','ids'):
            matches=re.findall(r'\b'+k+r': (0x[0-9a-f]+|\d+)',line);C.need(len(matches)==1,'row-scalar',k);fields[k]=int(matches[0],0)
        for k in ('inverse_symmetry','weights'): fields[k]=list(map(int,re.search(r'\b'+k+r': \((-?\d+), (\d+)\)',line).groups()))
        fields['symmetry_factors']=[list(map(int,v)) for v in re.findall(r'\((\d+), (\d+)\)',re.search(r'symmetry_factors: &\[([^\]]*)\]',line)[1])]
        C.need(set(fields)==set(FIELDS),'ten-row-fields')
        result[name]={'fields':{k:fields[k] for k in FIELDS},'line':line,'start':m.start(),'end':m.end(),'sha256':C.sha(line.encode())}
    C.need(len(result)==149,'149-Rows',len(result))
    return result

def world(case,mode,route):
    w=C.load(C.CAP/f'prepared/captured/{case}/{route}-{mode}.json')
    C.need((w['case'],w['mode'],w['route'])==(case,mode,route),'world-identity')
    live=w['live'];C.need(type(live['case']) is str and live['case']==case,'live-case')
    C.need(type(live['mode']) is str and live['mode']==mode,'live-mode')
    C.need(live['request']['classes']==mode,'selected-mode')
    C.need(live['grouping']=='none' and ('threads' not in live or live['threads']==1) and 'threads: Some(1)' in live['Run'],'mode-pool-grouping')
    roots={kind:C.PRE/f'execution/work/{case}--{route}-{kind}' for kind in ('pristine','observer')}
    payload=None;joins=[]
    for kind,root in roots.items():
        public=root/'dumps'/f'{case}.classes-{mode}-t1.json';source=root/'dumps'/f'{case}.classes-{mode}-t1.sources.json'
        b=C.read(public);C.need(b==C.read(source),'source-companion-bytes')
        p=C.load(public);C.need(p==w['payload'],'pristine-observer-payload',public)
        if kind=='pristine': payload=p
        census=C.read(root/'dumps'/f'{case}.census').decode()
        C.need(census.splitlines()[0]==live['report'],'pristine-report-binding')
        joins.append({'kind':kind,'payload':C.ref(public),'source_companion':C.ref(source),'census':C.ref(root/'dumps'/f'{case}.census')})
    ledger=roots['observer']/'ledger.jsonl';events=[json.loads(line,object_pairs_hook=C.unique) for line in C.read(ledger).decode().splitlines()]
    C.need([e['sequence'] for e in events]==list(range(len(events))),'ledger-global-sequence')
    stages=('small-gauge-live','small-gauge-raw') if route=='historical' else ('batch-live-inputs','batch-raw-diagrams','batch-raw-census')
    segment=[]
    for stage in stages:
        es=[e for e in events if e['stage']==stage and e['value']['case']==case and e['value']['mode']==mode]
        C.need(len(es)==1,'exact-mode-event',stage)
        segment.append(es[0]['sequence'])
        expected=w['live'] if stage==stages[0] else (w['raw'] if stage==stages[1] else None)
        C.need(es[0]['value']==expected if expected is not None else es[0]['value']['census']==payload,'exact-event-payload',stage)
    C.need(segment==sorted(segment),'mode-event-order')
    table=G.owner_table(w,route)
    row,details=derive(payload,live['report'])
    census=C.read(roots['pristine']/'dumps'/f'{case}.census').decode()
    if mode=='None':
        # Text census is the returned None reference. Auto is independently reconstructed.
        for d in details:
            C.need(f"{d['name']} key={d['key_digest']} id={d['id']}" in census and f"frame={d['frame_digest']} numerator={d['numerator_digest']}" in census,'independent-diagram-census')
    C.save(f'rows/{case}/{route}-{mode}.json',{'case':case,'mode':mode,'route':route,'Row':row,'derived_from':'Exact pristine JSON and complete pristine report; no assertion hash used as input.','native_literal_asserted':mode=='None','diagrams':details,'public_joins':joins,'live':live,'model_bindings':w['model_bindings'],'ledger':C.ref(ledger),'selected_event_sequence':segment,'owner_count':sum(len(x) for x in table.values())})
    return w,row,details

def apply_strict(before,patch):
    lines=patch.splitlines(keepends=True);C.need(lines[0]==f'diff --git a/{C.MAIN} b/{C.MAIN}\n','patch-path')
    C.need(lines[2]==f'--- a/{C.MAIN}\n' and lines[3]==f'+++ b/{C.MAIN}\n','patch-path')
    source=before.splitlines(keepends=True);out=[];cursor=0;i=4
    while i<len(lines):
        m=re.fullmatch(r'@@ -(\d+),(\d+) \+(\d+),(\d+) @@\n',lines[i]);C.need(m is not None,'patch-hunk')
        start,n,newstart,nn=map(int,m.groups());start-=1
        C.need(start>=cursor,'patch-overlap');out+=source[cursor:start];cursor=start;i+=1;oldn=newn=0
        C.need(len(out)==newstart-1,'patch-output-offset')
        while i<len(lines) and not lines[i].startswith('@@'):
            line=lines[i];tag=line[0];text=line[1:];C.need(tag in ' +-','patch-line')
            if tag in ' -': C.need(cursor<len(source) and source[cursor]==text,'patch-context-no-fuzz');cursor+=1;oldn+=1
            if tag in ' +': out.append(text);newn+=1
            i+=1
        C.need((oldn,newn)==(n,nn),'patch-count')
    out+=source[cursor:]
    return ''.join(out)

def main():
    before=C.read(C.O/'source-blobs/target'/C.MAIN);C.need(C.blob(before)=='7283198bf0152c82c1cd5f0173703cddfa7f5c0f','exact-target-blob')
    inventory=rows(before.decode());allworlds={};candidates={};physical=[]
    for ci,case in enumerate(C.CASES):
        for mode in C.MODES:
            h,hr,hd=world(case,mode,'historical');c,cr,cd=world(case,mode,'current')
            for k in ('model','request','options','Run','grouping','prepared_len','report'): C.need(h['live'][k]==c['live'][k],'physical-live-inputs',k)
            G.physical({'historical':h,'current':c},{'historical':G.owner_table(copy.deepcopy(h),'historical'),'current':G.owner_table(copy.deepcopy(c),'current')})
            C.need(len(hd)==len(cd)==(5,9,10)[ci],'case-cardinality')
            C.need(set(hr)==set(cr)==set(FIELDS),'all-ten-fields')
            C.need([k for k in FIELDS if hr[k]!=cr[k]]==['numerators'],'sole-actual-numerator-difference')
            # Exact ordered full keys and 128-bit IDs, not hash/ordinal-only matching.
            hm={(d['id'],json.dumps(d['key'],sort_keys=True)):d for d in h['payload']['diagrams']}
            cm={(d['id'],json.dumps(d['key'],sort_keys=True)):d for d in c['payload']['diagrams']}
            C.need(len(hm)==len(cm)==len(hd) and hm.keys()==cm.keys(),'full-id-full-key-correspondence')
            factor2=[d['id'] for d in c['payload']['diagrams'] if 'feynkit_generator_factor::{}::NumeratorIndependentSymmetryGrouping(2)' in d['overall_factor']]
            if case==C.CASES[0]: C.need(len(factor2)==4,'sym-factor2-four-of-five')
            physical.append({'case':case,'mode':mode,'graphs':len(hd),'owners_per_era':sum(len(B.owners(d['payload'])) for d in c['payload']['diagrams']),'model':c['live']['model'],'request':c['live']['request'],'options':c['live']['options'],'Run':c['live']['Run'],'groups':c['payload']['groups'],'factor2_ids':factor2,'physical_inputs_equal':True,'complete_worlds':{r:C.ref(C.CAP/f'prepared/captured/{case}/{r}-{mode}.json') for r in ('historical','current')},'association':'Exact full keys and128-bit IDs; complete physical fields compared by strict owner/physical predicates.'})
            allworlds[case,mode]=(h,c,hr,cr)
            if mode=='None':
                C.need(hr==inventory[case]['fields'],'historical-full-Row')
                # User observations are checked only after independent reconstruction.
                C.need(hr['numerators']==OLD[ci] and cr['numerators']==OBSERVED[ci],'observational-corroboration')
                candidates[case]={'case':case,'before_Row':hr,'after_Row':cr,'old':f"0x{hr['numerators']:016x}",'new':f"0x{cr['numerators']:016x}",'derived_from':f'rows/{case}/current-None.json','applied':False,'approval':False,'gates':['SCIENTIFIC_INTERFACE_REVIEW_PENDING','NEW_MIGRATION_REVIEW_REQUIRED']}
        nh,nc,_,_=allworlds[case,'None'];ah,ac,ahr,acr=allworlds[case,'Auto']
        candidates[case]['Auto']={'historical_Row':ahr,'current_Row':acr,'native_literal_asserted':False,'independent_mode_evidence':f'rows/{case}/current-Auto.json','None_Auto_row_differences':[k for k in FIELDS if candidates[case]['after_Row'][k]!=acr[k]],'full_id_key_correspondence':[]}
        for route,n,a in [('historical',nh,ah),('current',nc,ac)]:
            amap={(d['id'],json.dumps(d['key'],sort_keys=True)):d for d in a['payload']['diagrams']}
            for d in n['payload']['diagrams']:
                k=(d['id'],json.dumps(d['key'],sort_keys=True));C.need(k in amap,'Auto-own-domain')
                ad=amap[k];C.need(all(d[x]==ad[x] for x in ('payload','overall_factor','aut')),'Auto-own-payload')
                candidates[case]['Auto']['full_id_key_correspondence'].append({'route':route,'id':d['id'],'full_key':d['key'],'None_name':d['name'],'Auto_name':ad['name']})
    after=before.decode();scopes=[]
    for case,r in candidates.items():
        line=inventory[case]['line'];old='numerators: '+r['old'];new='numerators: '+r['new'];C.need(line.count(old)==1,'row-local-old-token')
        start=inventory[case]['start']+line.index(old)+len('numerators: ')
        scopes.append({'case':case,'field':'numerators','start_byte':start,'end_byte_exclusive':start+18,'old':r['old'],'new':r['new']})
        after=after.replace(line,line.replace(old,new),1)
    ab=after.encode();newrows=rows(after)
    C.need(len(before)==len(ab),'same-byte-length')
    changed=[k for k in inventory if inventory[k]['fields']!=newrows[k]['fields']]
    C.need(changed==[k for k in inventory if k in C.CASES],'exact-three-row-scope')
    reverse=bytearray(ab)
    for s in scopes:
        C.need(before[s['start_byte']:s['end_byte_exclusive']].decode()==s['old'] and ab[s['start_byte']:s['end_byte_exclusive']].decode()==s['new'],'exact-token-span')
        reverse[s['start_byte']:s['end_byte_exclusive']]=s['old'].encode()
    C.need(bytes(reverse)==before,'all-other-bytes-preserved')
    prior=C.load(C.T/'grouped-identical-migration-preparation-20261031-r1/preserved-prior12-Rows.json')
    prior_names=(*prior,'xs_qed_eemm_identical_L2');C.need(len(prior_names)==13,'prior13')
    preserved={}
    for case in prior_names:
        old=inventory[case];new=newrows[case];C.need(old['line']==new['line'],'prior13-full-byte-preservation',case)
        if case in prior: C.need(old['line']==prior[case]['row'],'prior12-baseline-binding',case)
        preserved[case]={'complete_Row':old['fields'],'line':old['line'],'sha256':old['sha256'],'before_after_byte_equal':True}
    patch=f'diff --git a/{C.MAIN} b/{C.MAIN}\nindex {C.blob(before)}..{C.blob(ab)} 100644\n'
    patch+=''.join(difflib.unified_diff(before.decode().splitlines(keepends=True),after.splitlines(keepends=True),fromfile='a/'+C.MAIN,tofile='b/'+C.MAIN,n=3))
    trial=apply_strict(before.decode(),patch).encode();C.need(trial==ab,'private-patch-exact-result')
    C.put('before.rs',before);C.put('after.NOT-APPLIED.rs',ab);C.put('proposed-not-applied.patch',patch.encode());C.put('source-blobs/private-patch-test/'+C.MAIN,trial)
    C.save('candidate-index.json',candidates);C.save('physical-input-bindings.json',physical);C.save('preserved-prior13-full-Rows.json',preserved)
    C.save('all149-Rows-before-after.json',{'before':inventory,'after':newrows})
    C.save('byte-scope-audit.json',{'status':'NOT_APPLIED','target_revision':C.TARGET,'path':C.MAIN,'before_blob':C.blob(before),'after_blob':C.blob(ab),'before_sha256':C.sha(before),'after_sha256':C.sha(ab),'patch_sha256':C.sha(patch.encode()),'bytes':len(before),'Rows':149,'changed_rows':changed,'unchanged_Rows':146,'preserved_prior_migrations':13,'changed_fields':['numerators'],'replacement_tokens':3,'token_spans':scopes,'actual_changed_byte_positions':[i for i,(a,b) in enumerate(zip(before,ab)) if a!=b],'all_other_bytes_equal':True,'patch_test':{'scope':'Private in-memory exact unified-diff application to frozen832 bytes; output stored only inside packet.','fuzz':0,'offset':0,'output_blob':C.blob(trial),'equals_after':True},'post_patch_runtime':'NOT_RUN','application_authorized':False})
    C.save('reconstruct-consumed-inputs.json',C.used)
    print('PASS:12 pristine mode/era Rows,3 candidate tokens,149 Rows,13 full prior Rows; after',C.blob(ab))

if __name__=='__main__': main()
