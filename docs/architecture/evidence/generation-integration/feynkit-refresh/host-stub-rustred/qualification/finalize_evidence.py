from pathlib import Path
import hashlib,json,os,shutil,subprocess,tomllib,ast,re,datetime
P=Path(__file__).resolve().parent;W=P.parents[2];D=Path(json.loads((P/'evidence-directory.json').read_text())['path']);A=json.loads((P/'artifact-qualification.json').read_text());b=A['source_bindings'];r=json.loads((P/'runtime/results.json').read_text());B=P.parent/'host-browser-correction-20261031-r1'
ENV={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')};ENV['GIT_OPTIONAL_LOCKS']='0'
def sha(p):return hashlib.file_digest(Path(p).open('rb'),'sha256').hexdigest()
def save(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n')
def copy(p,q):q=D/q;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
assert len(r)==8 and all(x['passed'] and x['exit_code']==0 and not x['secret_redaction_required'] for x in r)
qcd=json.loads((P/'runtime/private/work/qcd/qcd-boundary.json').read_text());assert qcd['complete_script_passed'] and qcd['targets']==[[2],[3]] and qcd['residuals']==[[1]] and qcd['trees_count']==2 and qcd['loop_groups']==6
assert qcd['sha256']=='2622f4f8f53624e177e6ae059eb12953f1b16f6cf08ebe85b70ff6155423dfc6';assert qcd['native_sha256']==A['records']['private']['native']['sha256']
script=Path(b['full_source'])/'crates/feynkit-py/tests/installed_feyncalc_qcd_ghost_vertex.py';assert sha(script)==qcd['sha256'];tree=ast.parse(script.read_text());assert len([n for n in ast.walk(tree) if isinstance(n,ast.Assert)])==23;assert 'max_depth=2' in script.read_text()
browser=json.loads((P/'runtime/private/work/browser/product.json').read_text());assert browser['complete'] and browser['tests_run']==1 and browser['failures']==browser['errors']==browser['skips']==0;assert browser['drawings']==[[0,2],[22,2],[0,42],[64,2]];assert [e['engine'] for e in browser['engines']]==['chromium','webkit'];assert [e['version'] for e in browser['engines']]==['153.0.8010.36','26.5'];assert '/webkit-2311/' in browser['engines'][1]['default_executable'];assert len(browser['pages'])==2 and all(x['dpr']==1 for x in browser['pages']);assert all(e['private_variables_absent'] for e in browser['engines'])
# Retain source-only correction and exact final private manifest/lock/registration patch.
patch=[]
for rel in ('Cargo.toml','Cargo.lock','src/lib.rs'):
 import difflib
 f='examples/notebooks/symbolica-host/'+rel
 patch.extend(difflib.unified_diff((Path(b['ordinary_source'])/f).read_text().splitlines(keepends=True),(Path(b['private_source'])/f).read_text().splitlines(keepends=True),fromfile='a/'+f,tofile='b/'+f))
(P/'private-host-final.patch').write_text(''.join(patch))
# Exact final source proof and qualification records. Large binaries and store inputs remain in place.
for f in sorted(P.iterdir()):
 if not f.is_file():continue
 if f.name in ('prepare-inputs.stdout','prepare-inputs.stderr','runtime-licensed.stderr') or f.suffix in ('.py','.json','.patch','.nix','.lock') or f.name.endswith(('.stdout','.stderr')):
  if f.stat().st_size<2_000_000:copy(f,Path('qualification')/f.name)
for prefix in ('pure','proposal','proposal-attempt1','qualification','artifacts','runtime'):
 for f in sorted((P/prefix).rglob('*')):
  if not f.is_file() or f.is_symlink() or f.stat().st_size>=2_000_000:continue
  rel=f.relative_to(P)
  if prefix=='pure' and ('graph.json' in f.name or f.name=='stdout' and f.parent.name.endswith('-graph')):continue
  if prefix=='artifacts' and f.suffix!='.json':continue
  if prefix=='runtime' and '/installed/' in '/'+str(rel)+'/':continue
  if prefix=='runtime' and '/work/' in '/'+str(rel)+'/' and f.suffix!='.json':continue
  copy(f,Path('qualification')/rel)
# Installed inventories are sufficient; private complete environments retained only locally for review.
reuse={'ordinary':{'native_comparison':'qualification/ordinary-native-comparison.json','prior_review':'reviews/host-browser-stub/review.json','actual_new_checks':['fresh complete-wheel install/RECORD/stub','ELF loader and containment','native import identity','canonical three-target exporter --check'],'prior_passes_reused_not_fresh':['three selected Rust embedding tests','installed UFO crossing','squared display','two sector signatures','two-engine native browser method'],'ordinary_QCD_failure':'Original missing-IBPFamily failure remains in exact prior packet; not rerun or erased. Current ordinary import still verifies absence of IBPFamily.'},'exporter':{'sha256':A['exporter']['sha256'],'owner_source_proof':'qualification/exporter-reuse-source-proof.json','fresh_check':True,'new_build':False},'private':{'new_composition':True,'actual_native_hash':A['records']['private']['native']['sha256'],'old_minimal_passes_not_transferred':True}}
save(D/'qualified-reuse.json',reuse)
retained={'workspace':str(W),'ignored_packet':str(P),'code':b['new_code'],'sources':{k:b[k] for k in ('full_source','full_source_nar','ordinary_source','ordinary_nar','private_source','private_nar','rustred_source','rustred_nar','proposal_source','proposal_nar')},'outputs':{n:{k:v[k] for k in ('drv','output','output_nar','wheel')} for n,v in A['records'].items()},'compiler':'qualification/compiler-binding.json','large_graphs':[str(p) for p in (P/'pure').glob('*graph.json')],'native_extractions':[v['native']['path'] for v in A['records'].values()],'fresh_private_installations':[str(P/'runtime'/n/'installed/venv') for n in A['records']],'purpose':'Retain exact store sources, two built distributions, final native/fingerprint/depfile provenance, fresh whole-wheel environments and failures for independent review. No shared cache cleanup or donor/workspace removal.'};save(D/'retained-artifacts.json',retained)
receipt={'schema_version':1,'kind':'corrected-stub-distributions-and-private-rustred-validation-host','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'NEW_INDEPENDENT_SOURCE_PACKAGING_PRIVATE_HOST_RUNTIME_REVIEW_REQUIRED','independent_approval_claimed':False,'base':'d0fb4e906140a6a0fc8dd6401ca839cab0c2734d','base_acceptance':'Root accepted bounded twelve-row step only under reviews/bounded-twelve-row; all restrictions retained','code':b['new_code'],'code_bookmark':'generation-host-stub-code-20261031','author':'Lucien Huber <im@lcnbr.ch>','code_change':{'only':'examples/notebooks/symbolica-host/python/symbolica/community/hepkit/__init__.pyi','approved_patch_sha256':'0666eba15326649c0e1fb3dc035bfe0bc8cbc1238ec5fcec99622ecbf2d58e47','old_sha256':'bf9d4ce3f4b5f3871c39164d3c11db248123eecd8a1939a845ad7569da4de92f','new_sha256':'5e2dfa31a10d9c545d38f4257b006daa7b1b53bcf0d96aafb060ee71a68f77e5','all149_literal_Rows_and12_migrations_and118_scope_bindings_preserved':True},'pure_build':{'first_coordinated_attempt':'Private --locked metadata exit101; ordinary output not realized. No native compilation in failing private lane. Preserved build-completion/logs.','private_lock_wiring_correction':'Only prune rustred-app/rustred-python and unused optional bridge edge; no version/source/feature change. Exact reviewed proposal and corrected final patch both retained.','final_coordinated_attempt':json.loads((P/'build-r2-completion.json').read_text()),'successful_wheel_outputs':2,'extra_archive_exporter_site_builds':0,'new_builds_used_old_stub':False,'no_semantic_RustRed_or_kernel_modification':True,'kernel':'8eb7c3729d6abefc916decc1bdf0c06536f8af80','rustred':'7c1ed03722b8c05daf60c89ba4ecc79457ed2ada','source_bindings':'qualification/source-bindings-r2.json','metadata_and_final_fingerprints':'qualification/artifacts/','ordinary_vs_private':'Ordinary remains minimal production composition. RustRed is only private validation composition; no production dependency/lock/registration change.','private_scope':'Single extension, single pinned kernel/PyO3 graph, features[], no campaign/app/python/reconstruction, upstream RUSTRED_RUNTIME_ARITIES=1.'},'artifacts':{n:{'drv':v['drv'],'output':v['output'],'output_nar':v['output_nar'],'wheel':v['wheel'],'native_sha256':v['native']['sha256'],'stub_sha256':v['stub_sha256']} for n,v in A['records'].items()},'runtime':{'fresh_complete_wheel_installations':2,'unlicensed_loader_probes':2,'containment_probes':2,'licensed_checks':8,'passes':8,'failures':0,'fresh_native_import_identities':2,'canonical_exporter_targets_checked':3,'QCD':{'script_sha256':qcd['sha256'],'assertion_sites_preserved':23,'targets':qcd['targets'],'max_depth':2,'residuals':qcd['residuals'],'trees':qcd['trees_count'],'loop_groups':qcd['loop_groups'],'completion_marker_present':True,'original_source_fixture_and_expectations':True,'observer_no_scientific_modification':True},'UFO':{'EW_crossings':36,'QCD_crossings':3,'Lorentz_linearity_probes':15,'original_marker':True},'squared_display':True,'sector_native_methods':2,'sector_three_stub_checks':3,'browser':browser,'budgets_timeouts_retries':'unchanged; no new cap or retry','compiler_or_nix_after_licensed_start':False,'private_delivery':'Existing exported license else current-user XDG config, privately in controller memory -> anonymous pipe -> immediate native exec. All four private variables excluded from pure children and browser descendants. No value/path/hash/file logged.','network_compiler_credentials_exclusion':'Fresh empty-root bwrap, read-only exact inputs, network namespace loopback only/no routes, no host root/config/Nixsocket/compiler/GPU. Stdout/stderr and controller receipts outside child writable root.','runtime_commands':'qualification/runtime/'},'reused':'qualified-reuse.json','original_failures_preserved':['Original minimal-host QCD IBPFamily absence','Old browser/stub setup and runtime failures under exact prior packets/reviews','Original proposed private metadata lock failure','Guard Git object-store discovery setup failure and Nix-v4 graph JSON adapter setup failure'],'mandatory_prior_distinctions':{'ASTRA-N1':'Historical combined_community .symtab/.strtab four hidden LLVM suffix changes remain distinct from allocated-content equality. Current ordinary comparison is separate and does not generalize generator loader-only deltas.','ASTRA-N2':'Historical workspace-root formatter cache write, ownership/cleanup and --no-cache correction remain disclosed; no perfect historical no-outside-write assertion.','prior_twelve_row_review':'All mandatory retained premises remain binding; no generator/Auto/math/native obligation reruns in this step.'},'pending_gates':['NEW independent source/packaging/private-host/runtime combined review','No full Community/campaign/multi-denominator certification','No broad generator/math/native obligation or full module/docs/site/Clippy execution credit','No WP3/M-rebuild/global integration/performance/distribution/publication approval'],'retention':'retained-artifacts.json','cleanup':{'owned_build_and_runtime_processes_reaped':True,'dispensable_lock_work_removed':True,'store_inputs_scratch_and_environments_retained_for_review':True,'no_publication_push_PR_fetch_upload':True,'no_CI_budget_or_global_integration_change':True}}
save(D/'receipt.json',receipt)
# Canonical note/registry are additive and deliberately carry no approval date/ref.
note=W/'docs/architecture/generation-host-stub-rustred.typ';note.write_text('''= Corrected notebook stub and private RustRed host qualification
<generation-host-stub-rustred>

This is a candidate for NEW independent source, packaging, private-host and
runtime review. Code `2cfb52b28d99f92305a9d571bbe3911cf04b8a19` applies only the
independently approved generated notebook stub patch to bounded base
`d0fb4e906140a6a0fc8dd6401ca839cab0c2734d`. The notebook stub changes from
`bf9d4ce3…` to `5e2dfa31…`; all three canonical targets now agree. All 149 literal
Rows, twelve migrations, every other production byte and all 118 existing
scope bindings remain unchanged. The prior twelve-row review retains all its
restrictions; this step grants no global integration or WP3 approval.

== Exact compositions and producers

The ordinary production notebook host remains minimal. The private validation
composition adds only a host manifest, standalone lock and registration delta,
using unchanged RustRed `7c1ed03722b8c05daf60c89ba4ecc79457ed2ada` and kernel
`8eb7c3729d6abefc916decc1bdf0c06536f8af80`. It is one native extension, with one
PyO3 0.28.3 and one Symbolica/Numerica/Graphica 3.0.1 graph. RustRed features
remain empty; campaign, Python/application adapters and reconstruction remain
disabled. Upstream `RUSTRED_RUNTIME_ARITIES=1` bounds this validation to one
denominator. No production RustRed dependency, registration or lock is added.

Both immutable distribution sources derive from the new code. Ordinary source
selection differs from the reviewed f47 host source only in the notebook stub;
the private source adds only those three host files. Pure offline evaluation
with IFD disabled passes. The first coordinated build stops at the required
locked metadata gate: the proposed lock still lists two disabled optional
adapters. A retained offline resolver reconciliation prunes those unused nodes
and the unused bridge edge, changing no package version, source or feature.
The corrected locked gates and one coordinated realization of the two required
wheel outputs pass. No archive, exporter or full-site producer is built.

The receipt binds actual wheel and output NARs, RECORD/METADATA, native ELF,
compiler, metadata, fingerprints, depfiles, source filters and all source deltas.
The old-stub proposal and old binaries remain separate historical inputs.

== Actual bounded runtime

Two separate fresh private environments install their complete exact wheels
offline. Packaged and installed notebook stubs both hash to `5e2dfa31…`.
Both actual loader and containment probes pass before activation. Eight licensed
checks pass, with zero failures: two native import identities, the canonical
three-target exporter check, and private-host QCD, UFO, squared display, sector
signatures and browser execution. The unchanged exporter binary is reused under
exact Rust-generator source invariance; expectations are never regenerated.

The original QCD script retains all 23 assertion sites and its exact source and
fixture layout. It completes with targets `[[2], [3]]`, depth 2, residual
`[[1]]`, two trees and six loop groups, printing its original completion marker.
UFO crossing passes 36 EW and three QCD crossings plus 15 Lorentz-linearity
probes. Squared display and both twelve-keyword sector signatures pass, including
the distinct diagram/family powers annotations.

The unchanged selected native browser method produces four drawings and passes
both Chromium 153.0.8010.36 and WebKit 26.5/revision 2311: two launches, two pages,
DPR 1 and zero skips, failures or errors. The exact reviewed Fontconfig and
software EGL inputs are retained; no tolerance or product/test source changes.

The ordinary native ELF is newly built. Every section except `.dynamic` and
`.dynstr` matches the reviewed old native byte-for-byte; only its output RUNPATH
changes, with new exact resolution and loader probes. Its prior bounded
coexistence/browser and three Rust embedding passes are explicitly reused,
not fresh. The old minimal-host QCD failure stays historical; the ordinary
host still has no IBPFamily. Only the private composition receives fresh QCD
and coexistence credit.

All pure queries, builds, qualification and installation precede the licensed
phase. Empty-root bwrap capsules mount exact read-only inputs and fresh writable
work, with loopback-only isolated networking and no routes, host credentials,
compiler, Nix socket or host driver tree. Records and stdout/stderr stay outside
the child writable mount. Private activation uses controller memory, an anonymous
pipe and immediate native exec; browsers receive no private variables. Budgets,
timeouts and retry policy remain unchanged. Concurrent work means build timing
is not isolated performance evidence.

== Evidence, distinctions and remaining gates

The additive packet is
`evidence/generation-integration/feynkit-refresh/host-stub-rustred/`.
Its receipt, retained-artifact index, exact reviews/seal bindings and checksums
keep every original packet, failure and old/current distinction. ASTRA-N1's
historical hidden static-symbol suffix changes and ASTRA-N2's formatter-cache
correction remain separate. There is no retrospective perfect no-write claim.

NEW independent review is required. Full Community, multi-denominator/campaign
coverage, broad generator or mathematical obligations, full modules/docs/site,
Clippy, WP3/M-rebuild, global integration and publication remain outside this
candidate. Only structural validation of this new note/registry is claimed;
no compiler or Nix command runs after the licensed boundary. Main and its index,
CI, resource budgets and the global integration bookmark are preserved. Exact
store inputs, native extractions and complete private environments remain
retained locally for review, as listed in `retained-artifacts.json`.
''')
registry=W/'docs/developers.toml';text=registry.read_text();marker='[[section.note]]\nid = "generation-raw-three-migration"';assert text.count(marker)==1
entry='''[[section.note]]
id = "generation-host-stub-rustred"
title = "Corrected notebook stub and private RustRed host qualification"
summary = "Two corrected distributions, a private one-denominator RustRed composition, and eight exact bounded runtime passes."
source = "docs/architecture/generation-host-stub-rustred.typ"
kind = "Execution record"
status = "New independent source, packaging, private-host and runtime review required"
lifecycle = "investigation"
owner = "unassigned"
products = ["feynkit", "gammaloop"]
topics = ["python", "testing", "provenance", "integration"]
review_triggers = ["examples/notebooks/symbolica-host/**", "docs/architecture/generation-host-stub-rustred.typ", "docs/architecture/evidence/generation-integration/feynkit-refresh/host-stub-rustred/**"]

''';registry.write_text(text.replace(marker,entry+marker))
oldreg=subprocess.check_output(['git','--git-dir='+str(W.parent/'main/.git'),'show',receipt['base']+':docs/developers.toml'],env=ENV);old=[n for s in tomllib.loads(oldreg.decode())['section'] for n in s.get('note',[])];new=[n for s in tomllib.loads(registry.read_text())['section'] for n in s.get('note',[])];assert [n for n in new if n['id']!='generation-host-stub-rustred']==old
newnote=next(n for n in new if n['id']=='generation-host-stub-rustred');assert not any(k in newnote for k in ('reviewed_at','review_ref','scope'));assert note.read_text().splitlines()[0]=='= '+newnote['title'];assert len([l for l in note.read_text().splitlines() if l.startswith('= ')])==1
save(P/'canonical-structure-check.json',{'registry_TOML':True,'all_prior_notes_unchanged':True,'one_matching_Typst_title':True,'invented_approval_fields':False,'full_docs_Clippy_or_Typst_compile':False,'no_compiler_or_Nix_after_licensed_boundary':True,'source_sha256':sha(note)})
# Hash/re-read every copied evidence file, without invoking anything from licensed runtime.
for f in (D/'qualification').rglob('*.json'):json.loads(f.read_text())
print(json.dumps({'durable_evidence':str(D),'checks':len(r),'passes':8,'private_QCD':True,'browser_drawings':4,'engines':2,'pages':2,'skips':0,'canonical_note':str(note)}))
