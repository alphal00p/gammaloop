"""Complete the ignored packet from frozen inputs; never apply the proposal."""
from pathlib import Path
exec(Path(__file__).with_name('evidence_io.py').read_text())
import difflib,shlex,tomllib

CODE='2cfb52b28d99f92305a9d571bbe3911cf04b8a19';EVIDENCE='457e361d7d39e0a3daae801651f8c1dfabbf47aa';GLOBAL='fdff0d18953348b84a0f1b7e6ab5ef703a5090de'
URI='git+file:///shared/lcnbr/dev/gammaloop/main?rev='+CODE+'&shallow=1&allRefs=1'
def emit(name,value):(P/name).write_text(json.dumps(value,indent=2)+'\n')
def own(name):return json.loads((P/name).read_text())
later=own('later-grouped-applicability.json');dec=own('decision-index.json')
dec.append({'acceptance':'.local/feynkit-refresh/root-grouped-identical-literal-acceptance-20261031.json','review':later['review'],'expected_sha256':later['review_sha256'],'actual_sha256':later['review_sha256'],'hash_matches':True,'verdict':'APPROVE','scope':'Exact single-token migration PREPARATION under pinned None/Auto premises; not application/runtime/integration approval','seal_sha256':later['review_seal_sha256']})
emit('decision-index.json',dec)
matrix=own('status-matrix.json')
for row in matrix:
 if row['gate']=='grouped None / Auto and finite runner':
  row['status']='accepted conditional aggregates/finite R3 AND single-token migration preparation; application integration review separate'
  row['limits']='At immutable457 the grouped literal is still old. The later token00fac78da0e0d159 -> a60df38a6b8079e3 is independently approved for the assigned writer, with NEW archive/39bounded passes/fresh grouped public joins and NEW integration review required. No postapplication result observed. None/Auto premises, operand AND output-axis maps and finite R3 restrictions remain. No need to wait for unrelated QCD proofs to review457.'
  row['evidence'].append('.local/feynkit-refresh/root-grouped-identical-literal-acceptance-20261031.json')
emit('status-matrix.json',matrix)

# Re-read only immutable objects for the proposal base. No moving working-copy chase.
files=['docs/architecture/evidence/generation-integration/execution.json','docs/architecture/generation-feynkit-community-refresh.typ','docs/developers.toml']
frozen={}
git=['git','--no-optional-locks','--git-dir='+str(W.parent/'main/.git')]
for i,f in enumerate(files):
 r=run('immutable457-metadata-'+str(i),git+['show',EVIDENCE+':'+f]);assert r['exit_code']==0;frozen[f]=r['stdout']
guard=own('metadata-proposal-guards.json')
for f in files:
 assert hashlib.sha256(frozen[f].encode()).hexdigest()==guard['base_file_sha256'][f], 'proposal base differs from457: '+f
guard['base_inputs_reverified_from_immutable_git_objects']=True
guard['later_grouped_preparation']=later['review_sha256']
guard['live_writer_movement']='Expected. Never apply this patch blindly to a later source/registry head; reconcile on a separately frozen metadata candidate and preserve writer changes.'
emit('metadata-proposal-guards.json',guard)

proposed=own('execution.NOT-APPLIED.json');current=proposed['current_status_reconciliation_20261031']
current['accepted_reviews']=[{'review':x['review'],'sha256':x['actual_sha256'],'scope':x.get('scope'),'verdict':x.get('verdict')} for x in dec]
current['grouped_scope']='Frozen457 retains old literal. Subsequent single-token preparation approved under7d96505d; assigned writer owns application, newarchive/39passes/joins and NEW integration review. No postapplication pass claimed.'
current['later_applicability']='Single test-literal delta conditionally preserves13 inspected non-generator compile source graphs; re-attest exact drv identities and runtime-input invariance on the writer frozen commit. Metadata/full-docs runtime closure treated separately.'
emit('execution.NOT-APPLIED.json',proposed)
patch=own('execution-reconciliation.NOT-APPLIED.json-patch');patch[-1]['value']=current;emit('execution-reconciliation.NOT-APPLIED.json-patch',patch)
addition=(P/'canonical-note-addition.NOT-APPLIED.typ').read_text()
addition=addition.replace('grouped literal and active migration review are outside frozen457.', 'grouped literal remains outside frozen457; subsequent preparation is approved.')
addition+='''
Subsequent preparation review7d96505dcb700eeae42b5a4e87676551ba910cb660df586cafe34d6b90e944f1
approves only token00fac78da0e0d159 to a60df38a6b8079e3 for xs_qed_eemm_identical_L2.
The separate writer owns its application on457, new generator archive,
39bounded row passes and fresh grouped joins, followed by NEW integration
review. This preflight claims no postapplication execution. The sole literal
is absent from thirteen inspected non-generator compile input graphs; those
producers may transfer only after exact derivation/runtime-input invariance
is re-attested on the new frozen head. Docs runtime and registry inputs differ
from compiled docs-library inputs and must be treated separately.
'''
(P/'canonical-note-addition.NOT-APPLIED.typ').write_text(addition)
note=frozen[files[1]];registry=frozen[files[2]];start=registry.index('id = "generation-feynkit-community-refresh"');end=registry.index('[[section.note]]',start);section=registry[start:end]
section=section.replace('summary = "Non-destructive source reconciliation and separately bounded provenance, finite-application and certificate reviews; aggregation remains unreviewed."','summary = "Separately accepted source, provenance, finite migration, pool1 and host evidence; new combined refresh review remains required."')
section=section.replace('status = "Bounded reviews recorded; metadata candidate and independent combined review pending"','status = "Bounded constituents accepted; NEW combined review and broader WP3 gates remain open"').replace('captured_at = "2026-10-30"','captured_at = "2026-10-31"').replace('evidence_revision = "cacd5182704082df929a066548e3501be22c9533"','evidence_revision = "'+EVIDENCE+'"')
newregistry=registry[:start]+section+registry[end:]
def scopes(text):return [(n['id'],s) for a in tomllib.loads(text)['section'] for n in a.get('note',[]) for s in n.get('scope',[])]
assert scopes(registry)==scopes(newregistry) and len(scopes(registry))==118
diff=''.join(difflib.unified_diff(note.splitlines(True),(note+addition).splitlines(True),fromfile='a/'+files[1],tofile='b/'+files[1]))+''.join(difflib.unified_diff(registry.splitlines(True),newregistry.splitlines(True),fromfile='a/'+files[2],tofile='b/'+files[2]))
(P/'metadata.NOT-APPLIED.patch').write_text(diff)

commands=own('commands.json');historical=own('historical-command-index.json')
commands['later_WP3_example_and_original_oracle_commands']={
 'classification':'Reference recipes for R/F completion; not selected or executed for this bounded preflight. Compilation and native execution must be split at the private-license boundary.',
 'pure_four_loop_build_argv':['env','-u','SYMBOLICA_LICENSE','-u','SYMBOLICA_LICENSE_PATH','-u','SYMBOLICA_LICENSE_SIGNED','-u','SYMBOLICA_LICENSE_FILE','CARGO_TARGET_DIR=/work/pure-four-loop-target','nix','--offline','develop',URI,'-c','cargo','build','--locked','-p','feynkit-generator','--example','four_loop_numerator','--profile','dev-optim','--message-format=json'],
 'pure_build_cwd':'Disposable byte-exact code2cf source copy; /work/pure-four-loop-target is a declared future output mount, not an artifact already present.',
 'runtime_four_loop_argv_each_pool':[['/work/pure-four-loop-target/dev-optim/examples/four_loop_numerator',str(n),'-','false'] for n in [1,8,32]],
 'runtime_four_loop_count_contract':'Actual4664 diagrams/306 zero numerators at each pool; also preserve original filtered true controls. Dev-optim example is a separately qualified artifact, not the ci-optim archive ELF.',
 'four_loop_perturbation_original':next((x for x in historical if x['path'].endswith('/four-loop-perturbation.argv')),None),
 'original_fast_oracle_source_recipe':['nix','--offline','develop',URI,'-c','just','generation-oracle','--profile','ci_gammaloop'],
 'oracle_environment':'FEYNKIT_ORACLE_THREADS=1,8,32; existing record/update flags absent. Recipe covers linnet, feynkit-generator and feynkit-py; accepted generator371 is not that whole recipe. Actual inventory mandatory; no364 or full8/32 claim from a pool1 package result.',
 'heavy':'WP9 specifies final-tip qcd_gse_L5@32 and four32-only rows at1/8/32 outside nextest. No new heavy run or guessed command selection is introduced here.'}
emit('commands.json',commands)

source_inputs={}
immutable=Path('/nix/store/3sj3nw06s3vhl3j12wfp5gc3r1h6f0k3-source')
for f in ['Cargo.toml','Cargo.lock','flake.nix','flake.lock','nix/rust-workspace.nix','nix/ci.nix','nix/documentation.nix','.config/nextest.toml','.config/hakari.toml','just/tests.just','examples/notebooks/symbolica-host/Cargo.lock','crates/tydenso/Cargo.lock']:
 q=immutable/f
 source_inputs[f]={'exists':q.exists(),'sha256':hashlib.sha256(q.read_bytes()).hexdigest() if q.exists() else None}
emit('immutable-source-input-hashes.json',{'source':str(immutable),'code':CODE,'inputs':source_inputs,'compiled_unit_features_and_phase_inputs':'producer-roots.json and demand-frontier.json','compiler_hashes':'tool-hashes.json','actual_native_qualification_this_task':False})

front=own('demand-frontier.json')
report={
 'schema':'generation-combined-integration-preflight/1','status':'READ_ONLY_RECONCILIATION_AND_PURE_PREFLIGHT_COMPLETE; NO_APPROVAL',
 'workspace':str(W),'packet':str(P),'frozen_code':CODE,'frozen_evidence':EVIDENCE,'global_integration':GLOBAL,
 'kernel':'8eb7c3729d6abefc916decc1bdf0c06536f8af80','source_scope':'Immutable457/2cf only; later grouped preparation additive applicability, no live-source chase.',
 'main_findings':[
  'Accepted plan requires NEW independent review of each bounded conjunction, not a prior migration of every outstanding numerator row.',
  'Constituent acceptances cannot supply combined approval. Globalfdff has not moved in this task.',
  'Source46 cumulative approval79dfeaab is real; old pending phrases and exclusions do not require repeating implementation.',
  'Exact2cf FeynKit group and package-specific generator derivations equal qualifiedffa inputs;12-row/full371 reruns or host rebuilds are unnecessary for metadata/stub alone.',
  'Current8eb production core/API/evaluator and selected library/docs test outputs are absent in the canonical local store. Their old1a executions cannot become current runtime credit.',
  'No additional generator/host artifact is missing for the already accepted bounded457 scopes. There ARE missing artifacts for extending the review to current production-core/evaluator qualification.',
  'Later grouped preparation accepted separately; no postapplication outcome claimed.13 non-generator compile source graphs omit its testfile; exact future derivation/runtime inputs must still be re-attested.',
  'All118 current frozen scope bindings validate, and proposed reconciliation preserves every ordered binding. Metadata proposal is NOT APPLIED.'],
 'minimum_new_combined_review_input_set':[
  'Exact frozen source/evidence and current truthful metadata delta with14 rehashed independent review identities; no selfapproval.',
  'SourceR4 + sharedScopeA + metadata48d0 within their own scopes; prior-eight restrictions and independent finite representation/migration approvals retained.',
  'Accepted twelve-row integration36passes15joins, full pool1 characterization420/371/302/69/49, applied-stub/ordinary/private-host acceptance8/8 with exact feature/loader/wheel limits.',
  'Separate None/Auto aggregate and finite R3 dispositions; operand AND output-axis maps, historical/cook/definition/observer/counter premises, no generic runner claim.',
  'Explicit unresolved current production-core/evaluator, full fast/oracle8/32 and WP3/M-rebuild/later gates. If reviewer scope includes fresh core claims, qualify the six-producer/core-evaluator menu first; those claims are currently unsupported.',
  'If reviewing later grouped application instead of457, use assigned writer frozen source, new archive39passes/joins and new integration review; retain this packet as conditional non-generator input qualification.'],
 'strict_gate_interpretation':'The accepted plan does not enumerate a standalone refresh-core/API/evaluator/library checklist before any bounded bookmark advance. Original R/F exit and historical bounded selections are recovered exactly in original-selection-records and historical-command-index; commands marks proposed impact selections separately. Only a NEW combined review/root disposition can decide the actual bounded advance. No new all69 or fullsite prerequisite is invented.',
 'resource_frontier':{'core_evaluator':commands['coordinated_build_core_and_evaluator']['frontier'],'all_recovered_consumers':commands['coordinated_build_all_recovered_consumers']['frontier'],'qualified_generator_host_reuse':'accepted-local-output-path-info.json; selected-producer-invariance.json','all_selected_outputs_missing':13,'unknown_output_paths':len(front['unknown_required_outputs']),'network_cache_queried':False,'measured_time_memory_disk':False,'builds_run':0},
 'gaps':[
  'No current production-core/evaluator native executions or actual discovery counts. Source-derived17 evaluator names are a future discovery assertion, not a pass.',
  'Broad current library/source46 colour qualification remains old-kernel/source-bound where no fresh exact run is recorded; code approval stays accepted.',
  'Docs builder runtime executable/catalog/exporter closure and final metadata inputs need resolution before strict licensed full docs execution; final site/WASM/SVG/packaging belongs later gates.',
  'Remote cache availability unqueried. Local demand counts are not compile-unit/time/memory estimates; offline build may expose further unavailable source outputs.',
  'Later grouped candidate commit and application results intentionally unknown. Do not transfer generator result credit across its changed test source.',
  'Preflight evidence depends on retained store/local packets; not a standalone publication bundle.'],
 'root_next_steps':[
  'Verify SHA256SUMS/SEAL and consume this packet as preflight data only.',
  'Coordinate at most one selected producer union if extending qualification: six core/evaluator roots or thirteen recovered consumer roots, never overlapping independent builds. Preserve current normal pinned owners/features/source filters.',
  'After pure artifact/loader/closure qualification, perform exact native discovery and only selected missing tests in the existing private no-network runtime boundary. Reuse accepted input-invariant row/host outputs; do not duplicate QCD capture/proofs.',
  'Have the metadata/source owner reconcile the NOT-APPLIED proposal on a separately frozen candidate; preserve expected grouped-writer changes and all118 scope records.',
  'Obtain NEW independent combined review and root scoped disposition before any refreshed integration bookmark advance. Keep WP3/M-rebuild and publication gates open.'],
 'files':{'status':'status-matrix.json','decisions':'decision-index.json','commands':'commands.json','roots':'producer-roots.json','frontier':'demand-frontier.json','filtered_sources':'source-filter-bindings.json','later_grouped':'later-grouped-applicability.json','metadata_patch':'metadata.NOT-APPLIED.patch','execution_patch':'execution-reconciliation.NOT-APPLIED.json-patch','proposal_guards':'metadata-proposal-guards.json','source_guards':'code2cf-store-parity.json','reads':'read-inventory.json','tools':'tool-hashes.json','errors':'errors.json'},
 'runtime_credit':0,'compiler_executions':0,'test_discovery_executions':0,'network_fetch_upload_publication':0,'delegation':0,'tracked_edits':0,'JJ_snapshot_ref_index_mutation':0,'license_config_reads':0,
 'retention':'This ignored packet and its source/store metadata retained for root continuation; no workspace removal or cleanup of others. No owned background process is intended to remain.'}
emit('preflight.json',report)

typ='''= NEXT combined refresh integration preflight

Frozen evidence457e361d7d39e0a3daae801651f8c1dfabbf47aa; code2cfb52b28d99f92305a9d571bbe3911cf04b8a19; kernel8eb7c3729d6abefc916decc1bdf0c06536f8af80.
Read-only reconciliation and pure preflight only. No build, compiler, native
discovery, runtime, license/config read, network, delegation, approval or ref
movement. All proposals below are NOT APPLIED.

== Exact next boundary

The accepted plan requires NEW independent review of each bounded conjunction
(plan:4–13,137–140). It does not require all69 old numerator rows to migrate
before every bounded refresh integration. Constituent verdicts do not approve
their conjunction. Global generation-integration remains fdff in this task.
Root must freeze truthful metadata and obtain NEW combined review with explicit
remaining production-core/evaluator and WP3 limits before any scoped advance.

No further generator/host artifact is missing for the already accepted457
bounded scopes. Actual selected generator drv b3psgc7y equals the qualifiedffa
producer; full FeynKit group hk3sccsq is also invariant. All4970 code2cf store
entries match Git. The only non-doc delta fromffa is the approved notebook stub;
457 adds docs only. Reuse12rows/full371/host evidence when actual runtime inputs
remain invariant. Matching69 Summary pairs are not physics or migration proof.

Current production-core/API/evaluator archives are absent and unqualified on8eb.
Old1a core2214, inference/BRST3, tensor202 plus17doctests and docs116 are historical
source/profile-bound credit. Source46 was cumulatively APPROVED at79dfeaab
(final Opus review ee6e37ee…); old pending wording is stale, not a reason to repeat
implementation. New8eb consumer qualification remains separate.

== Accepted constituent scopes

- SourceR4 and sharedScopeA; metadata48d0 attribution only.
- Twelve literal None rows:36 exact passes at1/8/32 and15 fresh joins. Auto
  physical captures retain historical counter/representation premises.
- Full pool1:420 identities,371selected,302pass69fail49ignored. No complete8/32 claim.
- Corrected ordinary wheel/stub plus separate private RustRed7c1 validation host:
  reviewer8/8. Ordinary lacks IBPFamily; private native arities[1],23assertions,
  targets2/3,depth2,residual1. Its tracing_max_level_info differs from ordinary
  host features. No production RustRed or full Community/campaign/multiple-denominator credit.
- None and Auto conditional grouped aggregates and exact finite runnerR3
  accepted separately. Operand AND output-axis maps mandatory; no generic runner.

The14 review JSON hashes are reverified in decision-index.json. Exact paths,
hashes, source identities, limits and per-gate citations are in status-matrix.json.

== Later grouped applicability

Subsequent preparation7d96505d… approves only token00fac78da0e0d159 to a60df38a6b8079e3.
The assigned writer owns application on457, new generator archive,39 bounded
passes/fresh grouped joins and NEW integration review. No postapplication result
is claimed. The changed test file occurs in none of13 non-generator compile
source graphs inspected, including required uncached ancestors. Re-attest exact
drv and runtime-input invariance at the writer frozen head. Docs runtime/registry
inputs differ from compiled docs-library inputs. Expected live writer movement
must be preserved; this packet stays on457 and does not wait for QCD proofs.

== Concrete qualification menu and resource frontier

commands.json supplies exact2cf Nix argv, immutable archive paths, selectors and
discovery-count contracts. Original selections, features, profiles and all80
historical argv records are retained separately. Focused core3, API topology1,
proposed evaluator16+1, ghost1, original NLO/CoolQM3, B2six, tensor6 and source46
controls are distinguished from the full canonical core lane and later R/F gates.
Historical counts are not substituted for actual current discovery.

The six-producer core/evaluator menu has58 missing required outputs and133
reused outputs over188 derivations. The optional13-producer library/docs union
has86 missing and159 reused over242. These are dependency counts, not measured
compile units, time, disk or memory. No remote cache query occurred. Select one
coordinated union if the intended review needs these fresh production claims;
these menus do not invent universal prerequisites for every bounded review.

Pure Nix children omit LICENSE, LICENSE_PATH, LICENSE_SIGNED and LICENSE_FILE
under the SYMBOLICA prefix; IFD is disabled and evaluation/store queries offline.
Build/qualify artifacts first with normal pinned owners. Only then perform
native listing/runs in the private no-network capsule, source at/build/source,
writable work/JUnit/state roots, required runtime tools prebound. No Nix/build
command receives the private license. Canonical240s,zero retries,16workers and
two memory-intensive slots remain; no20GiB cap is adopted.

== Metadata and later gates

NOT-APPLIED execution JSON Patch and canonical-note/registry diff preserve all118
ordered scope bindings, all verbatim reviews and original runtime identities.
The proposal adds current attribution outside historical records and claims no
combined approval. Its base hashes are verified against immutable457; reconcile
it with the source writer on a new frozen candidate before use.

WP3/M-rebuild stays incomplete: full fast/oracle matrix,3150,4664/306,
five fresh a_ddx processes, BRST/NLO/CoolQM, R-stage docs/assets and milestone
evidence remain distinct from bounded integration. Names, route, Classes,
correctness/WP9, performance, memory and final packaging are later work. Measure
memory before proposing a limit; no1.2x ceiling. Lucien Huber is settled author.
Publication/fetch/push/PR/labels and historical cleanup retain explicit user-go
gates (plan:665–680), with no general authorization here.

preflight.json lists remaining gaps and root next steps. errors.json retains
all read/probe/setup failures, including bounded Nix timeouts, an unlocked-store
reference rejection and a corrected Nix lambda. No test/build failure or pass
was produced. The ignored sealed packet is retained for continuation.
'''
(P/'preflight.typ').write_text(typ)
print(json.dumps({'review_hashes':len(dec),'scope_bindings_preserved':118,'NOT_APPLIED':True,'root_next_step_ready':True}))
