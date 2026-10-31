"""Seal the authorized candidate; pure analysis, no Nix or license access."""
from pathlib import Path
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess
import tomllib
import xml.etree.ElementTree as ET

assert not any(k in os.environ for k in ['SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED'])
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh');O=Path(__file__).resolve().parent
E=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/sm-eeww-migration'
BASE='48d0b98b97b368b2a2f58abe67cb19dddec623df';CODE='605178a65620b104805b012f1edc299710e8a176'
ENV=dict(os.environ,GIT_OPTIONAL_LOCKS='0',GIT_CONFIG_NOSYSTEM='1',GIT_CONFIG_GLOBAL='/dev/null')

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def write(p,data):Path(p).write_text(json.dumps(data,indent=2)+'\n')
def load(p):return json.loads(Path(p).read_text())
def git(*args):return subprocess.check_output(['git','--git-dir='+str(W.parent/'main/.git'),*args],env=ENV)
def blob(b):return hashlib.sha1(b'blob %d\0'%len(b)+b).hexdigest()

def scopes(data,reader):
    out=[]
    for sec in tomllib.loads(data.decode())['section']:
        for n in sec.get('note',[]):
            for s in n.get('scope',[]):
                b=reader(s['path']);t=b.decode();a=s.get('anchor');z=s.get('end_anchor')
                assert not a or t.count(a)==1;assert not z or t.count(z)==1
                r=t[t.index(a):t.index(z)].encode() if z else b
                assert blob(r)==s['digest'],(n['id'],s)
                out.append({'note':n['id'],**s,'actual_blob':blob(r)})
    return out

q=load(O/'archive-qualification.json');r=load(O/'runtime-results.json');assert q['code']==r['code']==CODE
assert r['actual_executions']==r['passes']==27 and r['failures']==r['skip_credit']==0
expected=set(x['name'] for x in load(O/'selected-inventory.json')['actual_compiled_selected']);assert len(expected)==9
for p in [1,8,32]:
    cmd=load(O/f'pool{p}/command.json');assert cmd['exit']==0 and cmd['inputs_before']==cmd['inputs_after']
    assert cmd['environment_before_exec']['FEYNKIT_ORACLE_THREADS']==str(p)
    assert cmd['environment_before_exec']['NEXTEST_TEST_THREADS']=='16'
    tests=ET.parse(O/f'pool{p}/nextest/ci_gammaloop/junit.xml').findall('.//testcase');assert len(tests)==9 and {c.attrib['name'] for c in tests}==expected
    assert all(c.find('failure') is None and c.find('error') is None and c.find('skipped') is None for c in tests)
inv=load(O/'inventory/stdout');suite=inv['rust-suites']['feynkit-generator::generation_oracle'];assert suite['binary-path']=='/archive/target/ci-optim/deps/'+Path(q['oracle_elf']['path']).name
assert sha(q['oracle_elf']['path'])==q['oracle_elf']['sha256'] and sha(q['archive'])==q['archive_sha256']
for p,info in load(O/'runtime-closure.json')['ELF_RUNPATH_NEEDED'].items():assert sha(p)==info['sha256']
for p,oid in load(O/'fetched-tree-check.json')['blobs'].items():
    f=Path(q['source'])/p;b=os.fsencode(os.readlink(f)) if f.is_symlink() else f.read_bytes();assert blob(b)==oid
path='crates/feynkit-generator/tests/generation_oracle/main.rs';before=git('show',BASE+':'+path);after=(W/path).read_bytes()
assert after==before.replace(b'numerators: 0xec3019645ae8d3e5',b'numerators: 0x0098a2a731e2a6fb',1)
assert blob(before)=='d61a4beebfc1b0d650cd17b001a1fdc6433b60d5' and blob(after)=='e6af4ec644688e1333064c50e3b157ae631473ee'
assert git('diff-tree','--no-commit-id','--name-only','-r',BASE,CODE).decode().splitlines()==[path]
oldreg=git('show',BASE+':docs/developers.toml');newreg=(W/'docs/developers.toml').read_bytes();oldsc=scopes(oldreg,lambda p:git('show',BASE+':'+p));newsc=scopes(newreg,lambda p:(W/p).read_bytes());assert len(oldsc)==118 and oldsc==newsc
oldnotes=[n for s in tomllib.loads(oldreg.decode())['section'] for n in s.get('note',[])];newnotes=[n for s in tomllib.loads(newreg.decode())['section'] for n in s.get('note',[])];assert [n for n in newnotes if n['id']!='generation-sm-eeww-migration']==oldnotes
note=W/'docs/architecture/generation-sm-eeww-migration.typ';title='= SM eeww one-token migration qualification';assert note.read_text().splitlines()[0]==title and sum(l.startswith('= ') for l in note.read_text().splitlines())==1
newnote=next(n for n in newnotes if n['id']=='generation-sm-eeww-migration');assert newnote['title']==title[2:] and 'review_ref' not in newnote and 'reviewed_at' not in newnote
ci=(W/'nix-ci.nix').read_text();assert re.search(r'^  enable = false;$',ci,re.M) and (W/'nix-ci.nix').read_bytes()==git('show',BASE+':nix-ci.nix')
start=load(O/'start.json');g=W.parent/'main/.git';main={p:sha(g/p) if (g/p).exists() else None for p in start['main_git']};assert main==start['main_git']
refs=lambda data:{r:h for h,r in (l.split() for l in data.splitlines())};oldrefs=refs(start['refs']);newrefs=refs(git('show-ref').decode());assert all(newrefs.get(k)==v for k,v in oldrefs.items())
integration=subprocess.check_output(['jj','--ignore-working-copy','log','--no-graph','--no-pager','-r','generation-integration','-T','commit_id'],cwd=W,env=ENV,text=True).strip();assert integration=='fdff0d18953348b84a0f1b7e6ab5ef703a5090de'
guard={'base_preserved':BASE,'code_preserved':CODE,'main_HEAD_index_packed_refs_equal':True,'preexisting_git_refs_equal':True,'new_git_refs':{k:v for k,v in newrefs.items() if k not in oldrefs},'generation_integration':integration,'CI':False,'all118_ordered_scopes_and_digests_unchanged':True,'scope_anchors_loosened':False,'preexisting_notes_unchanged':True,'new_note_approval_claim':False,'full_source_and_runtime_ELF_libraries_rehashed':True,'runtime_27_results_reparsed':True}
write(O/'final-guard.json',guard)
# Recompile only the final standalone note. No Nix/network/license operation.
typst=Path('/nix/store/wzpxqpc3v2xxah69rqifq5kac67hx36g-typst-0.15.1/bin/typst');args=[str(typst),'compile','--root',str(W),str(note),str(O/'canonical-note-final.pdf')];p=subprocess.run(args,cwd=W,env=ENV,capture_output=True,text=True);assert p.returncode==0,p.stderr
write(O/'typst-final.json',{'argv':args,'exit':p.returncode,'stdout':p.stdout,'stderr':p.stderr,'tool_sha256':sha(typst),'source_sha256':sha(note),'output_sha256':sha(O/'canonical-note-final.pdf'),'scope':'only new standalone canonical note; final PDF is ignored review artifact'})
names=['start.json','source-check.json','format.json','copied-evidence.json','flake.json','archive-derivation.json','producer-delta.json','producer-chain.json','selected-source-delta.json','fetched-tree-check.json','archive-qualification.json','compiler.stdout','selected-drv.txt','nix-complete.json','runtime-closure.json','elf-parser-binding.json','network-preflight.json','selected-inventory.json','runtime-results.json','runtime-driver.stdout','runtime-driver.stderr','runtime-driver-attempt1.stdout','runtime-driver-attempt1.stderr','qualify.stdout','qualify.stderr','qualify.stdout.premature-attempt1','qualify.stderr.premature-attempt1','patch-check.json','package-build-log.stdout','package-build-log.stderr','archive-build-log.stdout','archive-build-log.stderr','build.stdout','build.stderr','final-guard.json','typst-check.json','typst-final.json','initial-tool-setup.json','pure.py','qualify.py','runtime.py','runtime-attempt1.py','elf_helpers.py','finish.py']
(E/'qualification').mkdir(exist_ok=True)
for name in names:shutil.copyfile(O/name,E/'qualification'/name)
for f in O.glob('*.command.json'):shutil.copyfile(f,E/'qualification'/f.name)
for root in ['preflight','inventory','pool1','pool8','pool32']:
    dest=E/'runtime'/root;dest.mkdir(parents=True,exist_ok=True)
    for name in ['command.json','stdout','stderr']:shutil.copyfile(O/root/name,dest/name)
    if root.startswith('pool'):shutil.copyfile(O/root/'nextest/ci_gammaloop/junit.xml',dest/'junit.xml')
# Hash every retained local packet file; omit this manifest itself to avoid a cycle.
live={'local-artifacts.json','finish.stdout','finish.stderr'}
local={str(f.relative_to(O)):sha(f) for f in sorted(O.rglob('*')) if f.is_file() and str(f.relative_to(O)) not in live}
write(O/'local-artifacts.json',{'root':str(O),'files':local,'exclusions':sorted(live),'exclusion_reason':'Manifest itself and live seal-driver streams; completed evidence receipt/manifest is sealed separately.','retention':'Archive and extracted immutable producer, fingerprints, commands, failed setup attempts and outputs retained for independent review; no owned process remains.'})
write(E/'qualification/local-artifacts.json',load(O/'local-artifacts.json'))
receipt={'schema':'sm-eeww-one-token-integration-candidate/1','status':'NEW_INDEPENDENT_INTEGRATION_REVIEW_REQUIRED','created':datetime.datetime.now(datetime.UTC).isoformat(),'author':'Lucien Huber <im@lcnbr.ch>','workspace':str(W),'base_metadata_commit':BASE,'new_code_commit':CODE,'code_bookmark':'generation-sm-eeww-migration-code-20261031','decision':{'verdict':'APPROVE_BOUNDED_SM_ONE_TOKEN_MIGRATION','review_sha256':sha(E/'approval/review.json'),'attribution':'Independent Astra migration reviewer; root expressly adopted and authorized exact token plus routine qualification. This candidate is not independently integration-reviewed.','source_blob_before':blob(before),'source_blob_after':blob(after),'source_sha256_after':sha(W/path),'every_other_byte_equal':True,'other_nine_Row_fields_and_prior_eight_preserved':True,'proposed_patch_sha256':sha(E/'preparation/proposed-not-applied.patch')},'new_producer':{'uri':q['uri'],'archive_drv':q['archive_drv'],'archive':q['archive'],'archive_sha256':q['archive_sha256'],'oracle_elf':q['oracle_elf'],'source_paths':q['source_paths'],'selected_source_files':670,'dependency_units':452,'unchanged_dependency_layers':q['compiled_layers'][1:],'new_terminal_layer':q['compiled_layers'][0],'compiler':q['compiler']['rustc'],'profile':'ci-optim','kernel':'8eb7c3729d6abefc916decc1bdf0c06536f8af80','prior_proof_source':'5ba04864ca506e259d9eef4a54f65acdd8323fc2','prior_built_source':'5ea149ffd88ae3afc4072a38c0871c1329be8424','one_selected_Nix_realization':True,'Nix_pure_license_keys_removed':['SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED'],'Nix_calls_after_licensed_phase':0,'cache_policy':'Normal configured substitution; unchanged dependency producers reused. No cache transplant, relink, pin/feature/budget change.'},'runtime':r,'coverage':{'compiled_inventory':9,'oracle_suite_listed':183,'other_suites':'skipped by binary/package filter; zero execution credit','complete_Row_assertions':27,'class_modes':'Actual classes::check_request loops None and Auto; Named/Minimal, count, payload and source checks run; only returned None reference asserts historical Row. No duplicate Auto Row claim.','generation_pools':[1,8,32],'profile':'ci_gammaloop','timeout_seconds':240,'workers':16,'retries':0,'record_update_flags':'ABSENT','math_scope':'Finite SM proof remains pool1. Pools8/32 are runtime evidence only. No new algebra, observer or old-capture rerun.'},'source_scope_policy':{'all118_preserved':True,'anchors_and_digests_changed':0,'owner_test_expectation_bindings':'All existing main.rs owner regions exclude this SM row and remain byte-identical; no anchor change or blanket approval.','checks':['whole physical base tree before edit','exact reviewed patch and after blob','all149 fast/heavy literal Row source entries','all9 non-numerator fields','prior8 and every other row byte','all118 scope bindings','pinned rustfmt','one pure archive build','actual nine compiled identities','27 actual JUnit outcomes','standalone final Typst compile'],'additional_Clippy_full_docs_suite_credit':False,'Rust_policy_departures':'None; one lowercase hexadecimal token only.'},'carried_qualifications':{'preparation_global_VCS_guard':'FAILED during separate root-authorized grouped metadata; verbatim disclosure retained. No retroactive global guard pass.','formal_premises':'Exact original-input invariance, definitions and finite capture/decoder premises, not digest equality. No universal decoder/consumer theorem.','native_reuse':'32 SM substage batches/1022 obligations within72 total/4110; 8 SM baseline batches/230 obligations, 115 per mode; other24 successful control substages; zero new native batches.','controls':'15 implementation pure controls and8 additional reviewer predicates remain separate. Original inverse-QCD label perturbed a vertex; additive same-rank edge/vertex controls, journal completion/archived-version and duplicate4-graph qualifications preserved.','definition_limit':'Native subset-tag parsing is not attribute/tag identity evidence; exact definition tables remain premises.'},'failures_retained':[{'kind':'qualification infrastructure','result':'premature call KeyError before build completion; no new build or test credit','path':'qualification/qualify.stderr.premature-attempt1'},{'kind':'pure loader preflight','result':'libstd visited before libgcc carrier discovery; queue ordering corrected; no license delivered or tests executed','path':'qualification/runtime-driver-attempt1.stderr'}],'build_warnings':'Existing Clarabel unused import and elided-lifetime warnings retained; no lint fixes or test weakening.','remaining_gates':['NEW independent integration review of this candidate','other3 raw rows: declared interfaces and ordered consumers','groupedAuto and full371 lane','M-rebuild and combined/WP3/whole-workspace integration approval','installed-host/browser/stub/full-docs/performance/distribution/publication'],'guard':guard,'canonical_note':{'path':str(note.relative_to(W)),'sha256':sha(note),'standalone_compile':'PASS; final PDF ignored, no full-site/doc-suite claim'},'process':{'delegations':0,'owned_background_processes_remaining':[],'publication_push_fetch_PR_upload_label':False,'old_history_bookmarks_main_index_rewritten':False,'integration_ref_advanced':False,'CI':False,'budget_policy_changed':False,'twenty_GiB_proposal_adopted':False},'retention':{'tracked':'Scoped verbatim approval/bindings, actual producer qualification, commands, JUnit/logs and receipt.','ignored_holder':str(O),'ignored_purpose':'Actual archive/extraction, fingerprints, complete queries/failed setup artifacts and standalone PDF retained for independent review; shared caches and host-task productions untouched.'}}
receipt['failures_retained'].append({'kind':'initial tool lookup','result':'Unqualified python3 lookup exited127 before the guard script ran; corrected to explicit existing pinned interpreter. No test or license delivery.','path':'qualification/initial-tool-setup.json'})
write(E/'receipt.json',receipt)
files={str(f.relative_to(E)):sha(f) for f in sorted(E.rglob('*')) if f.is_file() and str(f.relative_to(E)) not in ['SHA256SUMS','SEAL.json']}
(E/'SHA256SUMS').write_text(''.join(h+'  '+p+'\n' for p,h in files.items()))
write(E/'SEAL.json',{'schema':'sm-eeww-migration-evidence-seal/1','status':'NEW_INDEPENDENT_INTEGRATION_REVIEW_REQUIRED','manifest_sha256':sha(E/'SHA256SUMS'),'files':len(files),'receipt_sha256':sha(E/'receipt.json'),'code_commit':CODE,'canonical_note_sha256':sha(note),'exclusions':['SHA256SUMS','SEAL.json'],'no_combined_approval':True})
print('sealed',len(files),'evidence files; receipt',sha(E/'receipt.json'),'27 actual passes; all118 scopes unchanged',flush=True)
