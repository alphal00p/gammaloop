"""New reviewer reconstruction. Standard library only; no candidate helper imports."""
from pathlib import Path
import ast,collections,copy,hashlib,json,re,traceback
from fractions import Fraction
O=Path(__file__).resolve().parent; W=O.parent.parent
P=W/'target/grouped-identical-migration-preparation-20261031-r1'
Q=W.parent/'generation-feynkit-refresh/.local/feynkit-refresh/grouped-none-qualification-20261030T191919'
reads={}
def read(p):
    p=Path(p);b=p.read_bytes();reads[str(p)]=hashlib.sha256(b).hexdigest();return b
def load(p):return json.loads(read(p))
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
def enc(x):return json.dumps(x,ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
def sha(b):return hashlib.sha256(b).hexdigest()
def need(x,s):
    if not x:raise AssertionError(s)
def fnv(b):
    h=14695981039346656037
    for v in b:h=((h^v)*1099511628211)%(2**64)
    return h
def u64(n):return int(n).to_bytes(8,'little')
def txt(s):return s.encode()+b'\xff'
def multi(vals):return fnv(b''.join(u64(v) for v in sorted(vals)))
def scalar(s):
    # Finite exact bookkeeping grammar, parsed as inert AST rather than evaluated code.
    t=re.sub(r'feynkit_generator_factor::\{\}::','',s).replace('^','**')
    def val(n):
        if isinstance(n,ast.Constant) and type(n.value)is int:return Fraction(n.value)
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -val(n.operand)
        if isinstance(n,ast.BinOp):
            a,b=val(n.left),val(n.right)
            if isinstance(n.op,ast.Mult):return a*b
            if isinstance(n.op,ast.Add):return a+b
            if isinstance(n.op,ast.Pow) and b.denominator==1:return a**b.numerator
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and not n.keywords:
            a=[val(v) for v in n.args]
            if n.func.id=='NumeratorDependentGrouping' and len(a)==3:return a[1]*a[2]
            if n.func.id in ['AutG','AntiFermionSpinSumSign','ExternalFermionOrderingSign','InternalFermionLoopSign'] and len(a)==1:return a[0]
        raise AssertionError('unexpected bookkeeping syntax '+ast.dump(n))
    return val(ast.parse(t,mode='eval').body)
def digests(payload):
    p=copy.deepcopy(payload);tokens=[]
    for k in ['numerator','numerator_prefactor','projector']:tokens.append(p.pop(k))
    for k in ['num_interface','numerator_prefactor_interface','projector_interface']:
        if k in p:tokens.append(p.pop(k))
    for v in p['vertices']+[e[1] for e in p['edges']]:
        tokens.append(v['numerator']);v['numerator']=None
        if 'num_interface' in v:tokens.append(v.pop('num_interface'))
    for k in ['model','symmetry_factor','overall_factor']:del p[k]
    return fnv(enc(p)),fnv(b''.join(enc(t) for t in tokens))
def census(d):
    contents=[];nums=[];ids=[];named=bytearray();records=[];sym=collections.Counter();wi=Fraction();wf=Fraction()
    for g in d['diagrams']:
        k=fnv(enc(g['key']));f,n=digests(g['payload']);a=int(g['aut'])
        need(a==g['payload']['symmetry_factor'] and g['overall_factor']==g['payload']['overall_factor'],'factor/aut binding')
        contents.append(fnv(u64(k)+u64(a)+txt(g['overall_factor'])+u64(f)))
        nums.append(fnv(u64(k)+u64(n)));ids.append(fnv(int(g['id'],16).to_bytes(16,'little')))
        named.extend((txt(g['name']) if g['name'] is not None else b'\xfe')+u64(k))
        sym[a]+=1;wi+=Fraction(1,a);wf+=scalar(g['overall_factor'])
        records.append({'id':g['id'],'key_sha256':sha(enc(g['key'])),'payload_sha256':sha(enc(g['payload'])),'key_fnv':f'{k:016x}','numerator_fnv':f'{n:016x}','keyed_numerator_fnv':f'{nums[-1]:016x}'})
    for g in d['groups']:
        named.extend(u64(g['master']))
        for m in g['members']:named.extend((txt(m['name']) if m['name']is not None else b'\xfe')+txt(m['ratio'])+txt(m['overall_factor']))
    r={'diagrams':len(records),'inverse_symmetry':[wi.numerator,wi.denominator],'weights':[wf.numerator,wf.denominator],'symmetry_factors':[[k,v] for k,v in sorted(sym.items())],'content':f'0x{multi(contents):016x}','numerators':f'0x{multi(nums):016x}','named':f'0x{fnv(named):016x}','ids':f'0x{multi(ids):016x}'}
    return r,records,sorted(nums)
def stage(rows,k):
    v=[x['value'] for x in rows if x['stage']==k];need(len(v)==1,'stage '+k);return v[0]
def events(rows,k):return [x['value'] for x in rows if x['stage']==k]
def project(g):return {k:v for k,v in g['json'].items() if k not in ['name','id']}
def owner(g):
    r={x['owner']:x['tensor'] for x in g['owners']};need(len(r)==len(g['owners'])==18,'18 owners')
    for k,t in r.items():
        if ':' in k:
            kind,i=k.split(':');p=g['json']['vertices'][int(i)] if kind=='vertex' else g['json']['edges'][int(i)][1];field,ik='numerator','num_interface'
        else:p=g['json'];field=k;ik={'numerator':'num_interface','projector':'projector_interface','numerator_prefactor':'numerator_prefactor_interface'}[k]
        need(t['raw']==t['decoded']==p[field],'owner expression')
        if ik in p:need(p[ik]==t['ports'],'declared logical order')
        need(len(t['ports'])==len(set(t['ports'])),'distinct owner labels')
    return r
def factual(rows):
    a=events(rows,'sample-zero');b=events(rows,'prepare-zero');c=events(rows,'group-cancellation');r=events(rows,'comparison-decision')
    need(len(a)==len(b)==30 and len(c)==29 and len(r)==3,'factual counts')
    need(all(type(x['sample_count'])is int and x['sample_count']==5 for x in a),'actual integer5')
    need(all(type(x['zero'])is bool and x['zero']is False for x in a+b),'actual false zero')
    need(all(type(x['cancelled'])is bool and x['cancelled']is False for x in c),'actual false cancel')
    ms=[]
    def walk(x):
        if isinstance(x,dict):
            if 'source_index'in x and 'overall_factor'in x:ms.append(x['ratio'])
            for v in x.values():walk(v)
        elif isinstance(x,list):
            for v in x:walk(v)
    walk(rows);need(len(ms)==210 and all(type(x)is str and x=='1' for x in ms),'210 ratio strings')
    need(collections.Counter(x['ratio'] for x in r)==collections.Counter([None,None,'1']),'3 comparison ratios')
    return {'samples_int5':30,'zeros_boolfalse':60,'cancellations_boolfalse':29,'member_string1':210,'ratios':[x['ratio'] for x in r]}
def main():
    d=load(Q/'pristine-runtime/dumps/xs_qed_eemm_identical_L2.classes-None-t1.json')
    rows=[json.loads(b) for b in read(Q/'observer-runtime/captures/ledger.jsonl').splitlines()]
    row,parts,sortednums=census(d)
    live=stage(rows,'live-inputs')
    # Report search counters only; no reported numerator token enters reconstruction.
    report=read(Q/'pristine-runtime/dumps/xs_qed_eemm_identical_L2.census').decode().splitlines()[0]
    row.update(topologies=int(re.search(r'topology_count: (\d+)',report)[1]),assignments=int(re.search(r'interaction_assignment_count: (\d+)',report)[1]))
    expected=load(P/'reconstructed-None-row.json')['expected']
    need([k for k in expected if expected[k]!=row[k]]==['numerators'],'nine other Row fields')
    need(row['numerators']=='0xa60df38a6b8079e3','derived token')
    need(census(dict(d,diagrams=list(reversed(d['diagrams']))))[0]['numerators']==row['numerators'],'sorted multiset order invariance')
    before=read(P/'main.before.rs');after=read(P/'main.after.NOT-APPLIED.rs');patch=read(P/'proposed-NOT-APPLIED.patch')
    lines=before.splitlines(keepends=True); matches=[x for x in lines if b'xs_qed_eemm_identical_L2: Row {' in x];need(len(matches)==1,'one row')
    orig=matches[0];repl=orig.replace(b'numerators: 0x00fac78da0e0d159',b'numerators: '+row['numerators'].encode());need(orig!=repl,'old token present')
    need(before.replace(orig,repl,1)==after and after.replace(repl,orig,1)==before,'all other bytes')
    removed=[x[1:] for x in patch.splitlines(keepends=True) if x.startswith(b'-') and not x.startswith(b'---')]
    added=[x[1:] for x in patch.splitlines(keepends=True) if x.startswith(b'+') and not x.startswith(b'+++')]
    need(removed==[orig] and added==[repl],'exact one-line patch')
    blob=lambda b:hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
    need(blob(before)=='9ef9f1ba005f20149716b22346d45c5501b62e87' and blob(after)=='7283198bf0152c82c1cd5f0173703cddfa7f5c0f','given immutable blobs')
    need(sha(before)=='6dacf61baca9f13f3e4a6d719f89eaced429fb3cb14c7904dcb41fe4c86b66df' and sha(after)=='66f01bcc855552ae888f1fb7c5bb158afbd448492b4e27a26e2bda40f09b07d2','before after pin')
    count=sum(bool(re.match(rb'\s+(?:fast|heavy) \w+: Row \{',x)) for x in lines);need(count==149,'149 rows')
    save('census-and-patch.json',{'row':row,'derived_from':str(Q/'pristine-runtime/dumps/xs_qed_eemm_identical_L2.classes-None-t1.json'),'per_diagram':parts,'sorted_keyed_numerator_multiset':[f'{x:016x}'for x in sortednums],'before_blob':blob(before),'after_unwritten_blob':blob(after),'patch_sha256':sha(patch),'Rows':count,'other_bytes_exact':True,'applied':False})
    need(stage(rows,'final-census')==d,'pristine/observer final')
    raw=stage(rows,'raw-diagrams');final=stage(rows,'final-diagrams');rc=stage(rows,'raw-census')
    approved=W/'target/xs-qed-l2/runtime/current/captures/xs_qed_eemm_L2/classes-None-t1'
    need(read(Q/'pristine-runtime/dumps/xs_qed_eemm_identical_L2.classes-None-t1.sources.json')==read(approved/'payload.json'),'complete raw payload bytes')
    for key in ['model','request','options']:need(live[key]==load(approved/(key+'.json')),'raw input '+key)
    need(live['Run']==load(approved/'capture.json')['run'],'raw Run')
    old=load(approved/'typed-owners.json');need(len(raw)==len(old)==len(rc['diagrams'])==30,'raw population')
    for g,c,pr in zip(raw,rc['diagrams'],old,strict=True):
        need(g['id']==c['id']==pr['id'] and project(g)==c['payload'],'raw id payload')
        expected_owners={a['owner']:{'raw':a['expression'],'decoded':a['producer_decoded_expression'],'ports':a['ordered_ports'],'is_metric':a['is_metric'],'is_composite':a['is_composite']} for a in pr['owners']}
        need(owner(g)==expected_owners,'540 exact raw owners')
    # Recompute every supplied source-phase axis map against original frames and accepted geometry maps.
    maps=load(W/'target/xs-qed-eemm-identical-fresh8eb-r3/maps.json')['current'];supplied=load(P/'source-phase-owner-axis-maps.json')
    pairs={p['raw_index']:p for p in events(rows,'staged-pair-final')};nperm=0
    for rec in supplied:
        idx=rec['raw_index'];side=rec['side'];m=next(x for x in maps if x['raw_index']==idx and x['side']==side);target=pairs[idx]['pair'][side]
        need(rec['full_map']==m and rec['source_native_ID']==raw[idx]['id'] and rec['target_native_ID']==target['id'],'phase exact ID/map')
        a,b=owner(raw[idx]),owner(target)
        for k,t in a.items():
            if k.startswith('vertex:'):dest='vertex:'+str(m['vertices'][int(k.split(':')[1])])
            elif k.startswith('edge:'):dest='edge:'+str(m['edges'][k.split(':')[1]][0])
            else:dest=k
            mapped=[re.sub(r'::hedge\((\d+),1\)',lambda z:'::hedge('+str(m['hedges'][z[1]])+',1)',v)for v in t['ports']]
            need(collections.Counter(mapped)==collections.Counter(b[dest]['ports']),'phase port multiset')
            perm=[b[dest]['ports'].index(v)for v in mapped];sr=rec['owners'][k]
            need(sr['target_owner']==dest and sr['source_ordered_ports']==t['ports'] and sr['mapped_ordered_ports']==mapped and sr['actual_target_ordered_ports']==b[dest]['ports'] and sr['source_axis_to_target_axis']==perm,'phase port permutation')
            nperm+=perm!=list(range(len(perm)))
    need(len(supplied)==60 and sum(len(r['owners'])for r in supplied)==1080,'60/1080')
    finalbindings=load(P/'final-physical-owner-bindings.json');chain=load(P/'membership-weight-payload-chain.json')['groups'];seen=[];weights=[]
    for i,(g,c,s,cb,cg) in enumerate(zip(final,d['diagrams'],finalbindings,chain,d['groups'],strict=True)):
        need(g['id']==c['id']==s['id']==cb['final_native_ID'] and project(g)==c['payload'],'final exact physical')
        need(owner(g)==s['owners'] and sha(enc(c['payload']))==s['physical_payload_sha256']==cb['payload_sha256'],'522 final owners')
        need(sha(enc(c['key']))==s['key_sha256']==cb['full_key_sha256'],'final full key')
        ww=[]
        for m,sm in zip(cg['members'],cb['members'],strict=True):
            idx=m['source_index'];seen.append(idx);need(m['id']==raw[idx]['id']==sm['native_ID'] and m['name']==raw[idx]['name'],'member identity')
            need(m['overall_factor']==rc['diagrams'][idx]['overall_factor']==sm['source_factor'],'source member factor')
            need(m['ratio']==sm['ratio']=='1','member ratio');ww.append(scalar(m['overall_factor']))
        w=sum(ww,Fraction());need(w==scalar(c['overall_factor']),'individual group weight');weights.append(str(w))
    need(sorted(seen)==list(range(30)),'all source members once')
    merged=[g for g in d['groups'] if len(g['members'])>1];need(len(merged)==1,'sole merge')
    need(sorted(m['id']for m in merged[0]['members'])==['02957fc4e974ffa2c5b5785e8683a949','5b4c349a74178fbed0bf4321725f53a3'],'exact unique merge IDs')
    # Auto uses its own independently approved raw data before cross-mode physical joins.
    ar=[json.loads(x)for x in read(W/'target/grouped-auto-fresh8eb-20261031-r1/observer-runtime/output/captures/ledger.jsonl').splitlines()]
    prior=[json.loads(x)for x in read(W.parent/'generation-feynkit-refresh/.local/feynkit-refresh/qualification/runtime12/observer/xs_qed_eemm_L2/world.jsonl').splitlines()]
    al=stage(ar,'live-inputs');pl=next(x['value']for x in prior if x['stage']=='batch-live-inputs' and x['value']['request']['classes']=='Auto')
    for k in ['model','request','options','Run','prepared_len']:need(al[k]==pl[k],'separate Auto '+k)
    aprior=next(x['value']['diagrams']for x in prior if x['stage']=='batch-raw-diagrams' and x['value']['mode']=='Auto');apr={g['id']:g for g in aprior}
    need(stage(ar,'raw-census')==next(x['value']['census']for x in prior if x['stage']=='batch-raw-census' and x['value']['mode']=='Auto'),'Auto raw census')
    for g in stage(ar,'raw-diagrams'):
        expected=copy.deepcopy(apr[g['id']]);expected['key']=g['key'];need(g==expected,'separate Auto owner records')
    autorec=load(P/'Auto-separate-exact-rebind.json');auto_raw=stage(ar,'raw-census')['diagrams'];none_raw=rc['diagrams']
    for r in autorec['raw_sources']:
        a,n=auto_raw[r['Auto_raw_index']],none_raw[r['None_raw_index']]
        need(a['id']==n['id']==r['exact_native_ID'] and a['key']==n['key'] and a['payload']==n['payload'],'Auto None complete raw physical')
    apairs={p['raw_index']:p for p in events(ar,'staged-pair-final')}
    for r in autorec['phase_correspondence']:
        a,n=apairs[r['Auto_raw_index']]['pair'][r['side']],pairs[r['None_raw_index']]['pair'][r['side']]
        need(project(a)==project(n) and owner(a)==owner(n),'Auto None complete phase')
    save('independent-bindings.json',{'None_raw_sources':30,'None_raw_owners':540,'Auto_raw_sources':30,'Auto_raw_owners':540,'phase_graphs_per_mode':60,'phase_owner_axis_maps':1080,'nonidentity_permutations':nperm,'final_groups':29,'final_owners':522,'group_weights':weights,'weight_histogram':dict(collections.Counter(weights)),'unique_merge_ids':sorted(m['id']for m in merged[0]['members']),'strict_facts':{'None':factual(rows),'Auto':factual(ar)},'historical_grouped_equivalence_claimed':False,'native_execution':0})
    print(json.dumps({'token':row['numerators'],'decimal':int(row['numerators'],16),'rows':count,'weights':dict(collections.Counter(weights)),'nonidentity_owner_permutations':nperm,'raw_and_Auto_bindings':'PASS'}))
try:main()
except BaseException:
    save('reconstruction-error.json',{'traceback':traceback.format_exc()});raise
finally:save('reconstruction-reads.json',reads)
