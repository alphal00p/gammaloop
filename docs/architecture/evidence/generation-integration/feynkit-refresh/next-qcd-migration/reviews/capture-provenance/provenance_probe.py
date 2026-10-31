"""Pure metadata/source/record checks; no native loading, builders or secret reads."""
from pathlib import Path
import ast
import difflib
import hashlib
import json
import re
import tomllib

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
P=ROOT/'target/raw-next-qcd-preflight-20261031-r1'
R1=ROOT/'target/raw-next-qcd-captures-20261031-r1'
R2=ROOT/'target/raw-next-qcd-captures-20261031-r2'
cache=json.loads((OUT/'frozen-before-files.json').read_text());extra={}
def digest(p):
    p=Path(p).absolute();key=str(p)
    if key in cache:return cache[key]
    if key not in extra:
        with p.open('rb') as f:extra[key]=hashlib.file_digest(f,'sha256').hexdigest()
    return extra[key]
def load(p):return json.loads(Path(p).read_text())
def need(ok,text):
    if not ok:raise ValueError(text)

plans=load(P/'capture-command-plans.json')['commands']
private={'SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_FILE','SYMBOLICA_LICENSE_SIGNED'}
closure=load(P/'runtime-closure.json');tools=load(P/'runtime-tool-bindings.json')
runtime=[]
for row in plans:
    args=row['argv'];env=row['environment_before_private_exec'];producer=row['producer']
    expected_env={'PATH':'/nix/store/4168xqs4l5jbfc7xyi9zrhdhszx0pnlb-python3-3.13.16/bin','HOME':'/home/runtime','USER':'runtime','LANG':'C.UTF-8','LC_ALL':'C.UTF-8','TMPDIR':'/work/tmp','XDG_CACHE_HOME':'/work/cache','PYTHONDONTWRITEBYTECODE':'1','FEYNKIT_ORACLE_THREADS':'1','FEYNKIT_ORACLE_DUMP':'/work/dumps','INSTA_UPDATE':'no'}
    if producer.endswith('observer'):expected_env['OBS_IDENTICAL_LEDGER']='/work/ledger.jsonl'
    need(env==expected_env and not private.intersection(env),'exact allowlisted semantic environment')
    need(args[-5:]==[row['case'],'--exact','--test-threads','1','--nocapture'],'exact single libtest selector')
    need('--unshare-all' in args and '--share-net' not in args and '--die-with-parent' in args and '--new-session' in args,'isolated namespaces and lifecycle')
    rw=[args[i+1:i+3] for i,x in enumerate(args) if x=='--bind']
    need(rw==[[row['writable_root'],'/work'],[row['writable_root']+'/home','/home/runtime']],'only owned writable binds')
    ro=[args[i+1:i+3] for i,x in enumerate(args) if x=='--ro-bind']
    expected_ro=[[str(P/'capsule'),'/data'],[str(P/'capsule/runtime-input'),'/input'],[row['program'],'/input/program'],[str(P/'exec_child.py'),'/input/exec.py']]+[[p,p] for p in closure[producer]['packages']]
    need(ro==expected_ro,'exact readonly capsule input and runtime closure mounts')
    need(args[args.index('--chdir')+1]=='/input','readonly current directory')
    need(not Path(row['records_root']).is_relative_to(Path(row['writable_root'])),'protected record location')
    for source,target in ro:
        need(source not in ['/','/home','/root','/nix','/common/build-cache'] and 'socket' not in source,'no broad host or cache mount')
    for p,meta in closure[producer]['ELF'].items():
        need(digest(p)==meta['sha256'],'actual loader dependency bytes')
        need(all(not e or e.startswith('/nix/store/') for e in meta['runpath']),'nonstore RUNPATH')
    runtime.append({'id':row['id'],'exact_plan_and_environment':True,'four_private_names_absent':True,'readonly_mounts':len(ro),'writable_binds':len(rw),'runtime_ELF_objects':len(closure[producer]['ELF'])})

controls=load(R1/'pure-control-results.json')
for r in controls['results']:
    if 'exact-mount-containment' in r['name']:
        need(r['result']['exit']==0 and r['result']['process_reaped'] and r['license_access_disabled'],'original mount probe completion')
        need(r['observed']['input_read_only'] and r['observed']['capsule_read_only'] and r['observed']['four_private_names_absent'],'original observed containment')
        need(r['observed']['network_interfaces']==[[1,'lo']],'loopback-only original probe')
    elif 'exact-loader-resolution' in r['name']:
        need(r['result']['exit']==0 and r['result']['process_reaped'] and r['license_access_disabled'],'original loader probe completion')
        text=(R1/(r['name'].replace('-exact-loader-resolution','-loader.stdout'))).read_text()
        need('not found' not in text,'loader resolution text')
need(len(controls['results'])==18,'r1 control groups')
need(len(load(R1/'pure-error-results.json')['results'])==3,'final r1 simulated failure outcomes')
need(load(R2/'pure-repair-controls.json')['positive_controls']==4 and load(R2/'pure-repair-controls.json')['rejecting_controls']==31,'candidate repair counts')
before=R1/'preparation-attempts/controller-pure-fixtures-v1.py'
def fns(p):return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(p.read_text()).body if isinstance(n,ast.FunctionDef)}
old,new=fns(before),fns(R1/'controller.py')
need([k for k in old if old[k]!=new[k]]==['capture'],'original version only durability changed')
for k in ['write','subreaper','descendants','reap','bounded']:need(new[k]==fns(R2/'controller.py')[k],'r2 lifecycle owner preserved')
for oldpath,newpath,patchpath in [(R1/'controller.py',R2/'controller.py',R2/'controller-r1-r2.patch'),(P/'bind_captures.py',R2/'bind_captures.py',R2/'binder-paths-and-exact-revalidation.patch')]:
    patch=patchpath.read_text();headers=patch.splitlines()
    calculated=''.join(difflib.unified_diff(oldpath.read_text().splitlines(True),newpath.read_text().splitlines(True),fromfile=headers[0][4:],tofile=headers[1][4:]))
    need(patch==calculated,'exact full source patch')

inventory=load(P/'source-inventory.json');source_count=0
for era,files in inventory['extracted'].items():
    for rel,desc in files.items():
        data=Path(desc['copy']).read_bytes()
        need(hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==desc['blob'],'source blob '+era+' '+rel)
        need(digest(desc['copy'])==desc['sha256'],'source SHA256')
        source_count+=1
    lock=tomllib.loads((P/'capsule/source'/era/'Cargo.lock').read_text())
    relevant=[x for x in lock['package'] if x['name'] in ['symbolica','numerica','graphica']]
    need(relevant==load(P/'source-kernel-pins.json')['pins'][era],'source lock kernel package records')
for case,versions in load(P/'expected-rows.json').items():
    texts={v['text'] for v in versions.values()}
    need(len(texts)==1,'identical selected Row literals at all five revisions')
    for era,v in versions.items():
        source=P/'capsule/source'/era/'crates/feynkit-generator/tests/generation_oracle/main.rs'
        need(v['text'] in source.read_text(),'actual selected Row source text')
scope=load(P/'source-comparison.json')['current_observer_base__current']
need(scope['compared_paths']==412 and [x['path'] for x in scope['differences']]==['crates/feynkit-generator/tests/generation_oracle/main.rs'],'selected source comparison extent')
a=(P/'capsule/source/current_observer_base/crates/feynkit-generator/tests/generation_oracle/main.rs').read_text().splitlines()
b=(P/'capsule/source/current/crates/feynkit-generator/tests/generation_oracle/main.rs').read_text().splitlines()
difference=[(x,y) for x,y in zip(a,b) if x!=y]
need(len(a)==len(b) and len(difference)==4 and all(': Row {' in x and ': Row {' in y for x,y in difference),'only four unrelated Row literal differences')

# Re-bind the actual historical diagnostic library/source, not merely its test adapter.
E=ROOT/'target/xs-qed-eemm-identical'
A=ROOT/'target/small-gauge-raw-parallel-20261030-sol-01/adapter-build/historical'
aq=load(A/'qualification.json');before=load(A/'qualification-before.json');command=load(A/'compile/command.json')
need(digest(A/'qualification-before.json')==aq['before_record_sha256'] and digest(A/'compile/command.json')==aq['compile_record_sha256'],'historical adapter actual compile records')
qpath=Path(before['observer_qualification']);q=load(qpath)
need(digest(qpath)==before['observer_qualification_sha256'],'historical diagnostic qualification exact binding')
need(q['base']=='2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94' and q['source']==q['source_after'],'historical diagnostic source before/after')
historical_source_checked=0
for path,expected in q['source'].items():
    need(digest(path)==expected,'historical diagnostic source '+path)
    historical_source_checked+=1
need(command['exit']==0 and command['inputs_before']==command['inputs_after'],'historical adapter compile success and inputs')
for path,expected in command['inputs_before'].items():need(digest(path)==expected,'historical compile dependency '+path)
gen=next(r for r in q['records'] if r['name']=='feynkit_generator')
need(digest(gen['artifact'])==gen['sha256']==command['inputs_before'][gen['artifact']],'historical linked diagnostic generator artifact')
observation=E/'sources/historical-observation-v2/crates/feynkit-generator/src/identical_observation.rs'
need(digest(observation)==q['source'][str(observation)]==digest(P/'capsule/helpers/historical-observer/2-observation.rs.txt'),'actual historical tensor observer source')

# Other exact producers are already covered by the full manifests. Rejoin qualification identities.
producers=load(P/'producers.json')['producers'];producer_report={}
for name,prod in producers.items():
    need(digest(prod['path'])==prod['sha256'],'actual producer ELF')
    qual=load(prod['qualification'])
    if name=='historical-pristine':need(qual['elf_sha256']==prod['sha256'] and qual['profile']=='dev-optim','historical pristine actual qualification')
    if name=='historical-observer':need(qual['elf_sha256']==prod['sha256'] and qual['base']==prod['code'],'historical observer actual qualification')
    if name=='current-pristine':need(qual['code']==prod['code'] and qual['archive_sha256']==load(P/'archive-binding.json')['sha256'],'current archive actual qualification')
    if name=='current-observer':need(qual['sha256']==prod['sha256'] and qual['code']==prod['code'] and qual['kernel']==prod['kernel'],'current observer actual qualification')
    producer_report[name]=prod

# Global inventories must describe exactly the original physical trees and local case tables.
capture=load(R2/'capture-record-manifest.json');actual=set()
for row in plans:
    for key in ['writable_root','records_root']:
        actual.update(str(p) for p in Path(row[key]).rglob('*') if p.is_file())
need(set(capture)==actual and len(capture)==120,'complete physical capture manifest')
worlds=load(R2/'complete-world-bindings.json')
need(set(worlds)=={str(p) for p in (R2/'prepared/captured').rglob('*') if p.is_file()} and len(worlds)==15,'complete world manifest')
for mapping in [capture,worlds]:
    for path,desc in mapping.items():need(digest(path)==desc['sha256'] and Path(path).stat().st_size==desc['bytes'],'exact inventory hash and size')
cells=load(R2/'six-cells.json');need(len(cells)==6,'six actual cells')
for cell in cells:
    local=load(R2/'case-audits'/cell['case']/'cells.json')
    need(cell['routes']=={r['route']:r for r in local if r['mode']==cell['mode']},'global/local cell identity')
    need(cell['pool']==1 and cell['grouping']=='None' and cell['historical_None_Row_native_pass']==(cell['mode']=='None') and not cell['Auto_literal_asserted'] and not cell['current_native_pass'],'honest cell native credit')
    for route,path in cell['worlds'].items():need(path==str(R2/'prepared/captured'/cell['case']/(route+'-'+cell['mode']+'.json')),'actual cell world path')
for row in plans:
    origin=R1 if row in plans[:3] else R2
    need(digest(Path(row['records_root'])/'result.json')==digest(origin/'outcomes'/(row['id']+'.json')),'original external protected result and outcome bytes')
need(len(list((R1/'outcomes').glob('*.json')))==3 and len(list((R2/'outcomes').glob('*.json')))==9,'3 old plus 9 fresh durable outcomes')

refs=load(R1/'review-read-bindings.json')
char=load(P/'capsule/qualification/accepted-characterization-review.json')
bounded=char['exact_inputs']['reused_bounded12_review'];refs[bounded['path']]=bounded['sha256']
refs[str(P/'capsule/qualification/accepted-characterization-review.json')]='fb688585bc7efd0b494fef7648df5abad398f88bfe8c439e8b42c0dd520eaa37'
for path,expected in refs.items():need(digest(path)==expected,'exact prior review identity')
report={'status':'PASS','runtime_plans':runtime,'r1_original_control_groups':18,'r1_original_reexecuted_negative_controls':23,
        'r1_final_source_simulated_errors':3,'original_unlicensed_mount_probes':4,'original_unlicensed_loader_probes':4,
        'r2_candidate_repair_controls':{'positive':4,'negative':31},'source_snapshot_blobs_verified':source_count,
        'historical_observer_actual_source_files_rehashed':historical_source_checked,
        'historical_observer_actual_generator_artifact_sha256':gen['sha256'],'historical_observer_actual_tensor_source_sha256':digest(observation),
        'producer_qualification_joins':producer_report,'actual_physical_manifest_files':120,'complete_world_manifest_files':15,
        'prior_reviews_hash_bound':refs,'source_invariance':{'paths_inspected':412,'different_paths':1,'different_unrelated_Row_literals':4},
        'archive_probe':{'command':'tar --zstd -xOf ARCHIVE target/ci-optim/deps/generation_oracle-3892f55f31822180 | sha256sum','observed_member_sha256':'5d7275d26392474f2490963dd0503f4f46ccebc30572a4f31f7a66228303b30f','extracted_files':0},
        'reused_qualification_limits':'Source/ELF/compiled closure facts accepted only for the exact retained artifacts and reviewed records. Prior Nix/store/build honesty assumptions remain. No rebuilt compiler/kernel theorem, no profile equivalence, and no repair of unrecorded historical grouped provenance.',
        'new_native_invocations':0,'new_loader_probes':0,'new_license_fixture_runs':0,'candidate_library_imports':0}
(OUT/'provenance-probe.json').write_text(json.dumps(report,indent=2)+'\n')
(OUT/'supplemental-bindings.json').write_text(json.dumps(extra,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ['runtime_plans','producer_qualification_joins','prior_reviews_hash_bound']}))
