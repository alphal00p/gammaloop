"""Additive durable evidence; no Nix/compiler/native operations."""
import datetime,gzip,hashlib,json,os,re,shutil,subprocess,tomllib
from pathlib import Path
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh');O=Path(__file__).resolve().parent
E=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/grouped-identical-migration';NOTE=W/'docs/architecture/generation-grouped-identical-migration.typ'
STATUS='NEW_INDEPENDENT_INTEGRATION_REVIEW_REQUIRED'
def load(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def blob(b):return hashlib.sha1(b'blob %d\0'%len(b)+b).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,indent=2)+'\n')
Q=load(O/'archive-qualification.json');RUN=load(O/'runtime-results.json');JOIN=load(O/'capture-joins.json');START=load(O/'start.json');INV=load(O/'source-invariance.json');DEP=load(O/'dependency-invariance-probes.json')
assert RUN['passes']==RUN['actual_executions']==39 and RUN['failures']==0
cases=INV['prior12']+['xs_qed_eemm_identical_L2'];assert len(cases)==13
for pool in [1,8,32]:
 rs=[r for r in RUN['results'] if r['pool']==pool];names=[t['case'] for r in rs for t in r['outcomes']]
 assert len(names)==13 and set(names)==set(cases) and all(r['exit']==0 for r in rs)
for p,r in START['base_tree'].items():
 if p in ['docs/developers.toml','crates/feynkit-generator/tests/generation_oracle/main.rs']:continue
 f=W/p;b=os.fsencode(os.readlink(f)) if r['mode']=='120000' else f.read_bytes();assert blob(b)==r['blob'],p
assert sha(W/'crates/feynkit-generator/tests/generation_oracle/main.rs')=='66f01bcc855552ae888f1fb7c5bb158afbd448492b4e27a26e2bda40f09b07d2'
REG=W/'docs/developers.toml';reg=tomllib.loads(REG.read_text());notes=[n for s in reg['section'] for n in s.get('note',[])];new=[n for n in notes if n['id']=='generation-grouped-identical-migration'];assert len(new)==1 and new[0]['status']==STATUS and 'reviewed_at' not in new[0]
sc=[]
for n in notes:
 for s in n.get('scope',[]):
  b=(W/s['path']).read_bytes();t=b.decode();a,z=s.get('anchor'),s.get('end_anchor');assert not a or t.count(a)==1;assert not z or t.count(z)==1
  x=t[t.index(a):t.index(z)].encode() if z else b;assert blob(x)==s['digest'];sc.append({'note':n['id'],**s})
assert sc==INV['scope_bindings'] and len(sc)==118
# After deleting only our one additive entry, the complete old registry semantics remain exact.
oldreg=subprocess.check_output(['git','--git-dir='+str(W.parent/'main/.git'),'show',START['base']+':docs/developers.toml'],env={k:v for k,v in os.environ.items() if k not in ['SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE']})
projected=tomllib.loads(REG.read_text())
for s in projected['section']:s['note']=[n for n in s.get('note',[]) if n['id']!='generation-grouped-identical-migration']
assert projected==tomllib.loads(oldreg.decode())
text=NOTE.read_text();assert text.count('ACTUAL_QUALIFICATION_PENDING') in [0,1]
chain=load(O/'producer-chain.json')
newbody=f'''The full immutable Git URI for code `{Q['code']}` supplies the actual
new package archive. One pure offline realization, with import-from-derivation
disabled, built only the selected generator terminal layer and archive. Normal
cached dependency outputs were reused. All four private licensing variables were
absent from Nix/build environments; no private license input was read for Nix.
The exact graph binds{Q['source_paths']} full source paths,{len(load(O/'selected-source-delta.json')['new_files'])} selected source files,
{len(chain['derivations'])} retained derivations and{Q['closure_units']} actual fingerprint units, compiler
rustc1.98.1, default+native and ci-optim. Kernel
`8eb7c3729d6abefc916decc1bdf0c06536f8af80`, dependencies, features, filters,
locks and policy remain unchanged. The actual terminal library payloads receive
new qualification:24 dependency artifacts checked; rebuilt rlib object payloads
and lib.rmeta are byte-identical, while compiler-session member-name changes
are explicitly bijected in ar/link metadata. Whole-rlib identity is not claimed.

Archive SHA-256: `{Q['archive_sha256']}`.
Actual executed oracle ELF SHA-256:
`{Q['oracle_elf']['sha256']}`.
The compiled ELF contains the new expected token and no old-token occurrence.
No old ELF was copied into this producer.

All13 exact complete migrated Rows passed at pool1, pool8 and pool32:39 actual
executions, zero failures, retries, skips or empty selections. Pool1 ran the new
grouped case alone with public dumps, followed by the exact prior12 without
dumps; pools8/32 ran all13 without dumps. Existing ci_gammaloop240-second
policy,16-worker cap and resource groups remain intact. Pools8/32 supply runtime
evidence only. No additional native capture test or full371 lane was run.

Pure source/ELF/dependency, no-network namespace, read-only-mount, loader and
compiled-identity qualification all precede licensing. Empty-root bwrap
capsules expose only read-only actual archive/source/helper/exact runtime
libraries and launcher stdlib; fresh work alone is writable. Command/stdout/
stderr records stay outside child-writable mounts. The controller preserves an
exported license or reads current-user XDG privately, then supplies an anonymous
pipe to immediate native execution. No license value/path/hash/file/argv dump
is retained. No Nix/compiler runs follow licensed start.

All five fresh grouped public files are byte-identical to the reviewed pristine
captures: None/Auto final29 payloads, their raw30 source companions and the
returned None census. Complete mode-specific physical keys/models/source
factors/member identities/ordered interfaces and2124 raw/final owner associations
join the accepted proof inputs. Model/options/Run conditions come from the
accepted per-mode inputs bound by exact actual source/filter/ELF invariance;
the public hook does not independently emit those inputs or flags. None Census
reconstruction yields the approved full Row/token. Auto search counters and its
literal remain original capture premises; no fresh Auto report assertion is
invented. There is no normalization, new decoder/observer/native algebra or
native layout run.'''
NOTE.write_text(text.replace('ACTUAL_QUALIFICATION_PENDING',newbody))
# Unlicensed structural checks; explicitly no Typst rendering after licensed start.
text=NOTE.read_text();titles=re.findall(r'^= (.+)$',text,re.M);assert titles==[new[0]['title']]
# Evidence link is created below, checked on the final packet pass.
for p in E.rglob('*.json'):json.loads(p.read_text())
commands=[load(p) for p in O.glob('*.command.json')];nix=[r for r in commands if r['argv'][0]=='nix'];boundary=load(O/'nix-complete.json')['completed'];licensed=load(O/'licensed-phase-begin.json')['time']
assert all(r['exit']==0 and r['end']<=boundary<licensed for r in nix)
assert sum('build' in r['argv'] for r in nix)==1
save(O/'pure-chronology.json',{'Nix_commands':len(nix),'selected_build_calls':1,'all_exit_zero_before_boundary':True,'boundary':boundary,'licensed_start':licensed,'Nix_calls_after_licensed_start':0,'compiler_calls_after_licensed_start':0,'four_private_variables_removed':load(O/'nix-complete.json')['license_keys_removed'],'realization_script_environment':'Four exclusions; no private input read. Later graph/qualification queries use explicit environment allowlist.'})
save(O/'canonical-structure-check.json',{'single_matching_title':True,'all118_scopes_unchanged':True,'old_registry_semantics_unchanged_after_removing_one_new_entry':True,'old_base_files_preserved_except_main_literal_and_additive_registry':True,'new_review_ref_pending':True,'new_approval_date':False,'Typst_rendered':False,'reason':'No Typst/compiler after licensed phase; source explicitly UNRENDERED','full_docs_Clippy_claimed':False,'applicable_Rust_policy':'Single lowercase hexadecimal literal only; pinned rustfmt passed; no design-policy departure.'})
# Copy own actual evidence, exact bytes or deterministic lossless gzip within1MiB.
bindings=[]
def cp(source,dest):
 source=Path(source);dest=E/dest;dest.parent.mkdir(parents=True,exist_ok=True);b=source.read_bytes();large=len(b)>1048576
 encoding='identity';parts=[]
 if large:
  compressed=gzip.compress(b,compresslevel=9,mtime=0);assert gzip.decompress(compressed)==b
  gz=Path(str(dest)+'.gz')
  if len(compressed)<=1048576:dest=gz;dest.write_bytes(compressed);encoding='gzip'
  else:
   encoding='gzip-chunks'
   for i,start in enumerate(range(0,len(compressed),524288)):
    f=Path(str(gz)+f'.part{i:03}');f.write_bytes(compressed[start:start+524288]);parts.append({'path':str(f.relative_to(W)),'sha256':sha(f),'bytes':f.stat().st_size})
   dest=Path(str(gz)+'.index.json');save(dest,{'encoding':encoding,'uncompressed_sha256':hashlib.sha256(b).hexdigest(),'uncompressed_bytes':len(b),'parts':parts})
   if gz.exists():gz.unlink()
 else:dest.write_bytes(b)
 assert dest.stat().st_size<=1048576
 bindings.append({'external_path':str(source),'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b),'durable_copy':str(dest.relative_to(W)),'copy_sha256':sha(dest),'copy_bytes':dest.stat().st_size,'encoding':encoding,'parts':parts,'durable_reference':'generation-grouped-identical-migration-evidence-20261031'})

for p in sorted(O.iterdir()):
 if p.is_file() and p.name not in ['retained-actual-archive.tar.zst','finalize.stdout','finalize.stderr']:cp(p,'qualification/'+p.name)
for p in sorted((O/'capsule-input').rglob('*')):
 if p.is_file():cp(p,'qualification/capsule-input/'+str(p.relative_to(O/'capsule-input')))
for p in sorted((O/'runs').rglob('*')):
 if p.is_file():cp(p,'runtime/'+str(p.relative_to(O/'runs')))
for p in sorted((O/'work/pool1-grouped/dumps').iterdir()):cp(p,'captures/'+p.name)
# Retain fresh terminal fingerprints and markers; no binary target/cache transplant.
for p in sorted((O/'new-fingerprints').rglob('*')):
 if p.is_file():cp(p,'fingerprints/'+str(p.relative_to(O/'new-fingerprints')))
save(E/'large-evidence-bindings.json',{'tracking_limit_bytes':1048576,'unchanged':True,'large_records':[b for b in bindings if b['encoding']!='identity'],'exact_original_bytes_available':True,'all_own_evidence_bindings':bindings,'large_native_archive':load(O/'retained-artifacts.json'),'source_ref':load(O/'code-freeze.json')['bookmark'],'external_kernel_dependency_store_inputs_retained':True})
review=load(E/'reviews/migration/review.json')
receipt={'schema':'grouped-identical-applied-migration-candidate/1','status':STATUS,'recorded_utc':datetime.datetime.now(datetime.UTC).isoformat(),'author':'Lucien Huber <im@lcnbr.ch>','workspace':str(W),'accepted_host_base':START['base'],'actual_code':Q['code'],'code_bookmark':'generation-grouped-identical-migration-code-20261031','evidence_bookmark':'generation-grouped-identical-migration-evidence-20261031','candidate_bookmark':'generation-grouped-identical-migration-candidate-20261031','literal_approval':{'review_sha256':sha(E/'reviews/migration/review.json'),'seal_sha256':sha(E/'reviews/migration/seal.json'),'root_acceptance':str(E/'roots/root-grouped-identical-literal-acceptance-20261031.json'),'case':'xs_qed_eemm_identical_L2','old':'0x00fac78da0e0d159','new':'0xa60df38a6b8079e3','before_blob':'9ef9f1ba005f20149716b22346d45c5501b62e87','after_blob':'7283198bf0152c82c1cd5f0173703cddfa7f5c0f','patch_sha256':'2450f56adfe23ea3b157b4d687c4ad6edd39866ff841f6ba8e85fa6f6d864fb7','all149_Rows_prior12_other9fields_every_other_byte_preserved':True},'source_invariance':INV,'actual_producer':{'code':Q['code'],'uri':Q['uri'],'source_paths':Q['source_paths'],'selected_source':Q['selected_source'],'selected_source_files':len(load(O/'selected-source-delta.json')['new_files']),'archive_drv':Q['archive_drv'],'archive':Q['archive'],'archive_sha256':Q['archive_sha256'],'archive_nar':Q['archive_nar'],'oracle_elf':Q['oracle_elf'],'compiler':Q['compiler'],'features':Q['actual_root_fingerprint']['data']['features'],'profile':Q['profile'],'kernel_packages':Q['kernels'],'actual_fingerprint_units':Q['closure_units'],'derivations':len(chain['derivations']),'terminal_dependency_artifacts':DEP['count'],'terminal_metadata_qualification':'Object payloads and lib.rmeta byte-exact; whole-rlib session member-name metadata separately qualified','compiled_layers':Q['compiled_layers'],'one_realization':True,'new_source_archive_ELF':True,'pure_chronology':load(O/'pure-chronology.json')},'actual_runtime':RUN,'fresh_public_joins':{'record':'qualification/capture-joins.json.gz','sha256_uncompressed':sha(O/'capture-joins.json'),'files':JOIN['files'],'mode_cells':2,'raw_graphs_per_mode':30,'final_graphs_per_mode':29,'retained_owner_associations':2124,'ordered_interface_and_flags_scope':'Explicit serialized ordered declarations freshly byte-joined; flags and implicit ports remain accepted P3/P4 association premises','None_complete_Row_reconstructed':True,'fresh_Auto_search_counters_or_literal_assertion':False,'normalizations':[],'new_observer_decoder_native_layout_algebra_runs':0},'scientific_record':{'domain':review['domain'],'WP9_argument':review['WP9_rationale'],'current_interface_source_argument':review['interface_disposition'],'accepted_premises':review['accepted_premises'],'accepted_None_candidate':'ad7ccbe7cee2373ed54fe16d17971491eb581e80','accepted_Auto_candidate':'113608b3917d903d0fcc730aa3760f114fbc3f04','separate_roots_and_full_reviews':'roots/ and reviews/accepted-{None,Auto}/','raw_three_and_SM_scopes_preserved':True,'old_raw_decoded_association_observer_definition_tag_subset_limits_retained':True},'retained_failures':load(O/'setup-failures.json'),'prior_full_lane':{'passed':302,'failed':69,'ignored':49,'fresh':False,'new303_68_49_claim':False},'host_scope':'Accepted separate private RustRed arity1/tracingfeature and ordinary absenceIBPFamily; no host rerun/newhostscope','validation':load(O/'canonical-structure-check.json'),'retention':{'local':str(O),'reason':'Actual new archive, extracted ELF, fingerprints, raw commands/dumps, original failures and external dependencies retained for pending independent review','live_processes':[],'compressed_durable_copies':'large-evidence-bindings.json','retained_input_bindings':'input-bindings.json'},'remaining_gates':['NEW independent applied integration review','All other rows/groupings and remaining full-lane failures','WP3/WP9/M-rebuild/global integration','Full docs/site/Clippy and unrendered standalone Typst','Performance/memory policy','Distribution/publication'],'global_integration':'fdff0d18953348b84a0f1b7e6ab5ef703a5090de','global_advance':False,'CI':False,'resource_policy_changed':False,'publication':False,'self_approval':False}
save(E/'receipt.json',receipt)
# Manifest binds complete committed packet + canonical note + registry; excluded seal is noncircular.
files=sorted([p for p in E.rglob('*') if p.is_file() and p not in [E/'SHA256SUMS',E/'SEAL.json']]+[NOTE,REG]);assert all(p.stat().st_size<=1048576 for p in files)
(E/'SHA256SUMS').write_text(''.join(sha(p)+'  '+str(p.relative_to(W))+'\n' for p in files))
save(E/'SEAL.json',{'schema':'grouped-identical-migration-evidence-seal/1','status':STATUS,'code':Q['code'],'accepted_base':START['base'],'manifest_sha256':sha(E/'SHA256SUMS'),'receipt_sha256':sha(E/'receipt.json'),'canonical_note_sha256':sha(NOTE),'files':len(files),'bytes':sum(p.stat().st_size for p in files),'actual_selected_executions':39,'actual_failures_retries_skips':0,'self_or_integration_approval':False,'declared_large_input_binding':'large-evidence-bindings.json','exclusions':['SHA256SUMS (bound by seal)','SEAL.json (exact hash returned separately)']})
for target in re.findall(r'#link\("([^"#]+)(?:#[^"\)]*)?"\)',NOTE.read_text()):
 if not re.match(r'^[A-Za-z]+:',target):assert (NOTE.parent/target).exists(),target
for p in E.rglob('*.json'):json.loads(p.read_text())
for b in bindings:
 f=W/b['durable_copy'];packed=b''.join((W[x['path']]).read_bytes() for x in b['parts']) if b['encoding']=='gzip-chunks' else f.read_bytes();raw=gzip.decompress(packed) if b['encoding']!='identity' else packed;assert hashlib.sha256(raw).hexdigest()==b['sha256'] and len(raw)==b['bytes']
save(O/'committed-seal.json',{'seal_sha256':sha(E/'SEAL.json'),'manifest_sha256':sha(E/'SHA256SUMS'),'receipt_sha256':sha(E/'receipt.json'),'packet_files':len(files)-2,'all_max_new_file_size_bytes':max(p.stat().st_size for p in files),'status':STATUS})
print('PASS evidence files',len(files),'all<=1MiB;39 real passes;5 exact joins;118 old scopes;unrendered Typst;seal',sha(E/'SEAL.json'),flush=True)
