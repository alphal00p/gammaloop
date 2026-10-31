from pathlib import Path
import ast,collections,copy,hashlib,json,re,traceback
O=Path(__file__).resolve().parent;W=O.parent.parent
P=W/'target/grouped-identical-migration-preparation-20261031-r1';Q=W.parent/'generation-feynkit-refresh/.local/feynkit-refresh/grouped-none-qualification-20261030T191919';F=W.parent/'generation-feynkit-refresh'
reads={}
def read(p):
    p=Path(p);b=p.read_bytes();reads[str(p)]=hashlib.sha256(b).hexdigest();return b
def load(p):return json.loads(read(p))
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
def need(x,s):
    if not x:raise AssertionError(s)
def sha(p):return hashlib.sha256(read(p)).hexdigest()
def stages(r,s):return [x['value']for x in r if x['stage']==s]
def one(r,s):
    a=stages(r,s);need(len(a)==1,s);return a[0]
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def project(g):return {k:v for k,v in g['json'].items()if k not in ['name','id']}

def uninstrumented_sources():
    patch=read(Q/'current-observational.patch');need(hashlib.sha256(patch).hexdigest()=='5de3e290402b284f0824145a5fcbe8baa6507de4dbf3f062391a22b427bbcf92','original overlay pin')
    lines=patch.decode().splitlines(keepends=True);i=0;result={}
    while i<len(lines):
        need(lines[i].startswith('--- '),'overlay old header');old=lines[i][4:].strip();new=lines[i+1][4:].strip();name=new[2:];i+=2
        patched=read(Q/'source'/name).decode().splitlines(keepends=True);out=[];cursor=0
        while i<len(lines)and not lines[i].startswith('--- '):
            m=re.match(r'@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@',lines[i]);need(m is not None,'overlay hunk');start=int(m[3])-1;i+=1;before=[];after=[]
            while i<len(lines)and not lines[i].startswith(('@@ ','--- ')):
                line=lines[i];need(line[0]in ' +-','overlay body')
                if line[0]in ' -':before.append(line[1:])
                if line[0]in ' +':after.append(line[1:])
                i+=1
            need(patched[start:start+len(after)]==after,'exact sealed overlay hunk')
            out.extend(patched[cursor:start]);out.extend(before);cursor=start+len(after)
        out.extend(patched[cursor:]);result[name]=None if old=='/dev/null'else ''.join(out).encode()
    return result

def main():
    # Reuse our new finite checker by loading only its own function definitions.
    own=read(O/'finite_named_consumer.py');tree=ast.parse(own)
    ns=globals().copy()
    defs=ast.Module(body=[x for x in tree.body if isinstance(x,ast.FunctionDef)and x.name in ['parse','analyze']],type_ignores=[])
    exec(compile(defs,str(O/'finite_named_consumer.py'),'exec'),ns)
    analyze=ns['analyze']
    nr=[json.loads(x)for x in read(Q/'observer-runtime/captures/ledger.jsonl').splitlines()]
    ar=[json.loads(x)for x in read(W/'target/grouped-auto-fresh8eb-20261031-r1/observer-runtime/output/captures/ledger.jsonl').splitlines()]
    auto=[]
    for phase,gs in [('raw',one(ar,'raw-diagrams')),('comparison',[p['pair']['comparison']for p in stages(ar,'staged-pair-final')]),('representative',[p['pair']['representative']for p in stages(ar,'staged-pair-final')]),('final',one(ar,'final-diagrams'))]:
        for i,g in enumerate(gs):a=analyze(g);a.update(phase=phase,index=i);auto.append(a)
    nd=one(nr,'final-census')['diagrams'];ad=one(ar,'final-census')['diagrams'];nb={canon(d['key']):d for d in nd}
    for d in ad:
        n=nb[canon(d['key'])];need(d['id']==n['id']and d['payload']==n['payload'],'Auto full final physical binding')
    need(len(nb)==len(ad)==29,'29 Auto final sources')
    # Reconstruct both 58 cut/threshold unions and complete representative descent.
    chains=[]
    for mode,r in [('None',nr),('Auto',ar)]:
        collapse=one(r,'collapse-entry');pairs=one(r,'group-entry')['pairs'];plans=stages(r,'collapse-plan');final=one(r,'final-diagrams');record=[]
        for group,plan,g in zip(collapse['groups'],plans,final,strict=True):
            pi=group['master'];master=pairs[pi];rep=master['representative']if master['kind']=='distinct'else master['comparison'];changed=['cuts','topology_threshold_candidates','overall_factor','name','id','loop_momentum_basis']
            need({k:v for k,v in g['json'].items()if k not in changed}=={k:v for k,v in rep['json'].items()if k not in changed},'representative physical descent')
            need({o['owner']:o['tensor']for o in g['owners']}=={o['owner']:o['tensor']for o in rep['owners']},'representative owner descent')
            need(g['json']['loop_momentum_basis']==rep['json']['loop_momentum_basis'],'representative routing')
            for field,pk in [('cuts','cuts'),('topology_threshold_candidates','thresholds')]:
                union=[];mf=collapse['frames'][pi];mg=master['comparison']['json']
                for m in group['members']:
                    si=m['diagram'];sg=pairs[si]['comparison']['json'];sf=collapse['frames'][si]
                    for part in sg[field]:
                        sides=[]
                        for side in ['left','right']:
                            ports=part[side]['half_edges']if field=='cuts'else part[side]
                            canonical_vertices={sf['vertices'][sg['edges'][h['edge']][0][h['endpoint']]]for h in ports}
                            selected={j for j,v in enumerate(mf['vertices'])if v in canonical_vertices}
                            sides.append([{'edge':ei,'endpoint':ep}for ei,e in enumerate(mg['edges'])for ep in ['source','target']if e[0][ep]in selected])
                        if sf['left_right_swapped']!=mf['left_right_swapped']:sides=sides[::-1]
                        union.append(sides)
                need(set(map(canon,union))==set(map(canon,plan[pk]))and len(set(map(canon,plan[pk])))==len(plan[pk]),'complete partition union')
                flips=set(master.get('reversed_edges',[]))
                def normalize(hs,flip):
                    return sorted([{'edge':h['edge'],'endpoint':('target'if h['endpoint']=='source'else 'source')if flip and h['edge']in flips else h['endpoint']}for h in hs],key=lambda h:(h['edge'],h['endpoint']))
                expected=[[normalize(l,True),normalize(rr,True)]for l,rr in plan[pk]]
                actual=[[normalize(part[side]['half_edges']if field=='cuts'else part[side],False)for side in ['left','right']]for part in g['json'][field]]
                need(sorted(map(canon,expected))==sorted(map(canon,actual)),'final signed partition transport')
                record.append({'id':g['id'],'field':field,'union_cardinality':len(plan[pk])})
        chains.append({'mode':mode,'unions':record,'groups':29})
    # Full upstream seals, not merely the verdict strings.
    seal_counts={}
    for s,base in [(W/'docs/architecture/evidence/generation-xs-qed-eemm-identical-aggregation-r1/seal.json',W),(W/'docs/architecture/evidence/generation-xs-qed-eemm-identical-auto-aggregation-r1/seal.json',W/'docs/architecture/evidence/generation-xs-qed-eemm-identical-auto-aggregation-r1'),(Q/'seal.json',Q)]:
        spec=load(s)
        for name,rec in spec['files'].items():
            f=base/name;need(sha(f)==rec['sha256']and f.stat().st_size==rec['bytes'],'upstream seal '+str(f))
        seal_counts[str(s)]=len(spec['files'])
    char=load(P/'accepted-characterization-review-guard.json')
    for r in char['files']:need(sha(r['path'])==r['sha256'],'characterization review guard')
    # Source content is bound to snapshots and supplied hashes; no ref/index reads.
    originals=uninstrumented_sources()
    source_paths=[r['path']for r in load(P/'source-invariance.json')['bindings']]+['crates/idenso/src/tensor/inference.rs','crates/spenso/src/structure/ordered.rs','crates/spenso/src/structure/permuted.rs','crates/spenso/src/structure/representation.rs','crates/idenso/src/cook.rs']
    sources=[]
    for p in source_paths:
        old=originals[p] if p in originals else read(Q/'source'/p);current=read(F/p);selected=read(W/p)
        if p.endswith('generation_oracle/main.rs'):
            need(hashlib.sha256(old).hexdigest()=='aedd07670ad16ae72ae1c96970668c80ea0ecd5d6d5992c736af71de7d00999b','original main pin')
            need(current==read(P/'main.before.rs'),'current target before blob')
            project_tokens=lambda b:re.sub(rb'numerators: 0x[0-9a-f]+',b'numerators: TOKEN',b)
            need(project_tokens(old)==project_tokens(current),'only old/current numerator tokens')
        else:need(old==current==selected,'invariant source '+p)
        sources.append({'path':p,'original_source_sha256':hashlib.sha256(old).hexdigest(),'refresh_current_sha256':hashlib.sha256(current).hexdigest(),'selected_sha256':hashlib.sha256(selected).hexdigest()})
    # Every sealed original source path under the relevant code/assets is compared,
    # while original observer overlay sources remain separate and unmodified.
    domain=['crates/feynkit-generator/','crates/feynkit-graph/','crates/feynkit-model/','crates/idenso/','crates/spenso/','assets/']
    same=0;deltas=[]
    for path in (Q/'source').rglob('*'):
        if not path.is_file()or path.is_symlink():continue
        rel=path.relative_to(Q/'source').as_posix()
        if not(any(rel.startswith(x)for x in domain)or rel in ['Cargo.toml','Cargo.lock']):continue
        old=originals[rel] if rel in originals else read(path)
        if old is None:continue
        if (F/rel).is_file()and old==read(F/rel):same+=1
        else:deltas.append(rel)
    need(deltas==['crates/feynkit-generator/tests/generation_oracle/main.rs'],'finite producer-source domain deltas')
    save('supplementary-v3.json',{'Auto_graphs_independently_checked':auto,'Auto_final_complete_payloads_equal_None_by_full_key':29,'partition_chains':chains,'seals':seal_counts,'characterization_review_files':len(char['files']),'source_bindings':sources,'source_domain_unchanged_files':same,'source_domain_deltas':deltas,'source_commit_identity':'Original5ba archive5ea kernel8eb retained; d0fb/ffa remain supplied immutable applicability identities, no currentHEAD claim','VCS_invocations':0})
    print(json.dumps({'Auto_graphs':len(auto),'Auto_full_final':29,'unions_per_mode':58,'seals':seal_counts,'source_equal_files':same,'source_deltas':deltas}))
try:main()
except BaseException:save('supplementary-v3-error.json',{'traceback':traceback.format_exc()});raise
finally:save('supplementary-v3-reads.json',reads)
