"""Seal only this independent review directory; never modify candidate inputs."""
from pathlib import Path
import datetime
import hashlib
import json

O=Path(__file__).resolve().parent;W=O.parents[1]
P=W/'target/raw-next-qcd-preflight-20261031-r1'
R1=W/'target/raw-next-qcd-captures-20261031-r1'
R2=W/'target/raw-next-qcd-captures-20261031-r2'
def read(path):return json.loads(Path(path).read_text())
def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
frozen=read(O/'frozen-after.json');parser=read(O/'parser-probe.json');data=read(O/'data-probe.json');prov=read(O/'provenance-probe.json')
assert frozen['pass'] and frozen['original_before_after_equal']
assert parser['status']==data['status']==prov['status']=='PASS'
rows=[]
for case in ['xs_qcd_gguu_sym_L1','xs_qcd_gguu_L1','xs_qcd_gggg_L1']:
    native=[r for r in parser['record_replays'] if r['id'].startswith(case+'--')]
    current=[r for r in native if r['actual'] is not None]
    assert len(current)==2 and current[0]['actual']==current[1]['actual']
    rows.append({'case':case,'None_diagrams':current[0]['actual']['diagrams'],
                 'native_invocations':4,'historical_native_passes':2,'current_native_failures_exit101':2,
                 'current_only_different_field':'numerators','old_numerators_hex':hex(current[0]['original_Row']['numerators']),
                 'current_numerators_hex':hex(current[0]['actual']['numerators']),
                 'both_modes_observed':True,'Auto_literal_asserted':False})
receipt=read(R2/'invocation3-revalidation.json')
report={
 'schema':'independent-next-qcd-capture-provenance-review/1',
 'reviewer':{'model':'Astra','role':'independent NEW review sub-agent','fallback':'User-authorized Astra fallback','delegations':0},
 'completed_utc':datetime.datetime.now(datetime.UTC).isoformat(),
 'decision':'approve',
 'verdict':'APPROVE_EXACT_NEXT_THREE_QCD_CAPTURE_PROVENANCE_ASSOCIATION_AND_NARROW_CONTROLLER_REPAIR',
 'blocking_findings':[],
 'workspace':str(W),'review_directory':str(O),
 'scope':{
   'cases':['xs_qcd_gguu_sym_L1','xs_qcd_gguu_L1','xs_qcd_gggg_L1'],
   'Classes':['None','Auto'],'Grouping':'None','pool':1,
   'meaning':'Only factual provenance, retained capture associations, exact saved libtest outcome classification, and the r2 Summary framing/additive prerequisite repair for this completed finite stage. No scientific, migration, finite-interface, general-runner or global integration approval.',
   'authorized_write_scope':'This new ignored review directory only; all main/source/candidate/frozen inputs read-only.'},
 'frozen_input_verification':frozen,
 'candidate_actual_counts':{
   'planned_native_invocations':12,'r1_actual_invocations':3,'r2_fresh_invocations':9,'r1_invocations_reused_without_rerun':3,
   'combined_actual_native_invocations':12,'native_reruns':0,'native_retries':0,
   'historical_passed_tests':6,'current_failed_tests_exit101':6,'tests_selected_per_invocation':1,
   'public_mode_world_copies':24,'complete_observer_worlds':12,'case_mode_cells':6,
   'algebra_or_interface_native_invocations':0,
   'mode_contract':'Each exact test internally executes both requested modes. Only the returned None Summary receives the original literal Row assertion; Auto is not assigned a literal Row pass.'},
 'per_case':rows,
 'original_failure_and_repair':{
   'r1_original_result':receipt['source_result_path'],'r1_original_result_sha256':receipt['source_result_sha256'],
   'native_exit':101,'native_test_pass':False,'original_controller_validation':False,
   'original_controller_error':{'phase':'assertion_boundary','type':'ValueError','predicate':'ten Summary fields'},
   'reproduced_cause':'Original flat regex yields 12 entries by including nested Fraction denominator fields. It fails before validating the saved ten-field Summary; that controller error is not native success and is not a semantic rejection.',
   'original_missing_inputs_after_preserved':True,
   'repair_sources':{str(R2/name):sha(R2/name) for name in ['controller.py','row_boundary.py','revalidation.py','controller-r1-r2.patch','binder-paths-and-exact-revalidation.patch']},
   'parser':'Balanced delimiters, exactly ten unique top-level Summary fields, exact Fraction fields and integer syntax, valid symmetry tuples; then the existing parsed_summary/expected_row owner. A separate reviewer recursive grammar independently agrees on all six saved failing Summary pairs and the immutable three Rows.',
   'additive_receipt':{'path':str(R2/'invocation3-revalidation.json'),'sha256':sha(R2/'invocation3-revalidation.json'),'original_artifact_bindings':len(receipt['bindings']),'original_before_and_current_pure_after_inputs':len(receipt['inputs_before'])},
   'receipt_limit':'The after hash set is a new pure rehash, corroborated by immutable r1 integrity evidence and read-only mounts. No original result.inputs_after record is invented.',
   'prerequisites':'Only the exact xs_qcd_gguu_sym_L1--current-pristine ID can use the fixed receipt path and pinned receipt hash; it rechecks the original result/command/streams/artifacts, immutable expected Row and parser hash. All other prerequisites resolve their own protected records. Original two historical results remain bound to r1. First-three indices are prohibited in the r2 capture path.',
   'unchanged_runtime_owners':['write','subreaper','descendants','reap','bounded'],
   'unchanged_guards':'Exact frozen plans/producer hashes, sequence and matching-pristine gates, fresh resolved work/record paths, protected results, 240-second deadline, zero retries, private pipe and reaping remain. Source patches reproduce the complete old-to-new files.',
   'r2_failed_pure_attempt_preserved':'preparation-attempts/v1-* retains the tuple-versus-JSON-list receipt mismatch and its source/error; it preceded the corrected pure controls and the nine native captures.'},
 'capture_associations':{
   'public_joins':'All 12 case/era/mode joins independently check pristine JSON == observer JSON == both source companions as exact bytes. Six case/era census pairs are byte-identical. Census lines independently bind exact string IDs, names, aut, source factors and momentum bases.',
   'manifest_inventory':'The 120 physical capture/record files exactly equal the twelve actual record/work tree inventories. The 15 complete-world bindings exactly equal twelve prepared worlds plus three binding files. All sizes/hashes and root six-cell/per-case tables agree.',
   'whole_ledgers':data['ledgers'],
   'phase_and_mode_binding':'Historical complete segments are live, one into-representative event per actual diagram, then raw; current complete segments are live, raw, census. None then Auto are each selected by actual case/mode and request.classes; sequences are contiguous and complete. Null OrderClasses is never used as Auto evidence. Older sample/prefix events are not substituted for the final complete diagram domain.',
   'identity':'Joins use native 32-hex string IDs and full exact key objects, with duplicate/missing domains rejected. Diagnostic rg.json.id is an f64-rounded serde_json::Value field and is explicitly excluded from identity decisions; ordinal-only mapping is not used.',
   'owners':'65/117/118 owners per era and mode: 1,200 complete raw/decoded/ports/flags records and 600 historical/current input pairs checked. Every retained tensor, exact full key, owner label, ordered port list, present/absent interface field and flag matches its direct source world.',
   'definitions':'Independent full-consumption syntax/token scans reproduce every complete captured symbol table, including tags/attributes, for raw tensors, decoded expressions, ordered ports and source factors. These tables are captured spellings only, not uncaptured callbacks or user data.',
   'world_inventory':data['worlds'],
   'shared_inputs':'For each cell, full live model, request, options, grouping, Run and prepared length agree between eras. Exact diagram physical descriptors, incidence/routing payloads, full keys, factors and ungrouped source membership/ratios agree after removing only the compared numerator/interface fields.',
   'symmetry_factor':'Exactly four of the five symmetrized diagrams retain NumeratorIndependentSymmetryGrouping(2); the fifth lacks that head. Plain cases have no such factor. No symmetry-key overlap or old one-PI overlap earns equivalence/native credit.',
   'old_onePI_native_or_proof_reuse':'None; fresh complete owner-input association does not inherit old native mathematics/interface proof.'},
 'pipeline_observations':{
   'cells':data['pipeline'],
   'meaning':'The first differing captured representation is at the raw cooked tensor spelling boundary. This is distinct from decoded syntax and later public interface serialization, and is not necessarily the first internal production divergence.',
   'decoded_syntax':'After the five exact source-qualified spelling maps, all 65 sym gguu and all 117 plain gguu decoded owner trees match per mode; gggg has 107 of 118 matches, leaving 11 syntax differences per mode.',
   'serialization':'Historical explicit public interfaces are absent for all owners. Current explicit interfaces occur on 34/63/8 owners per case and mode; absence is not rewritten as an empty list. Raw port order and flags remain preserved separately.',
   'attribution_limit':'No deeper algebraic cause follows from names, equal counts/digests, or matching decoded trees. No physical correctness, numerator equivalence or finite interface policy is approved.'},
 'actual_producer_binding':{
   'producers':prov['producer_qualification_joins'],
   'historical_source':'Actual retained original source 2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94; selected literals introduced at 8cf6cecaf22be4d5c62d467b8e6827e80c6ba66d. All five source snapshots retain those same three Row texts; fresh historical selected tests are the actual new pass evidence.',
   'historical_kernel':'Symbolica/Numerica git c78cffc20f318df1fce50539a587cbe8ffc717c0; Graphica registry 3.0.1. Original dev-optim and actual diagnostic compiler/features/flags remain distinct from current CI qualifications; no all-git-c78 or profile-equivalence claim.',
   'historical_observer_supplement':{'source_files_rehashed':3551,'actual_linked_generator_sha256':prov['historical_observer_actual_generator_artifact_sha256'],'actual_observer_source_sha256':prov['historical_observer_actual_tensor_source_sha256'],'method':'Actual adapter qualification-before hash -> historical manual-build qualification -> equal source/source_after map and linked generator artifact -> adapter compile before/after hashes -> retained adapter ELF; all relevant current file bytes rechecked. The original tensor observer equals the sealed capsule source.'},
   'current_pristine':'Actual ffa76ec420642c03349ab756d2d7720ab93c0057 / archive8993e6c8346b84f6d91c7e65c4c4fc5ed83cf6d932e969bb7fc55cbdced20740 / kernel8eb7c3729d6abefc916decc1bdf0c06536f8af80. Pure streaming extraction of its selected archive member reproduces executed ELF5d7275d26392474f2490963dd0503f4f46ccebc30572a4f31f7a66228303b30f.',
   'current_observer':'Actual built5ea149ffd88ae3afc4072a38c0871c1329be8424 plus qualified additive observer overlay, kernel8eb. Proof source5ba, actual built5ea and evidencec344 remain distinct. Inspected 412-path invariance to ffa has only one changed file with four unrelated Row literals; no broader equivalence claimed.',
   'source_snapshot_blobs_verified':759,
   'prior_review_identities':prov['prior_reviews_hash_bound'],
   'reused_scope':'Shared ScopeA565caeace0112a30fb9118f94b2bb759e4267d2b802395c314e485b89b4b2b70 supplies its exact producer/closure provenance only. Corrected gauge disposition6e4410dd16ef6c2366b08793af0e835c21e55526bf124a2d3581283c32f1a9a2 supersedes the report blockers of initial1591408c within its stated finite scope. Bounded12 review535b7926 and full-pool1 characterizationfb688585 remain separate; neither supplies next-three scientific proof.',
   'limits':'Build/Nix/store honesty and qualified lock/vendor/compiler provenance remain inherited premises. This review does not rebuild, verify a general compiler theorem or repair old historical grouped/source gaps. eae/1a remains the failing baseline, never the original reference.'},
 'capsule_and_controller':{
   'actual_records':'All twelve actual command plans equal the prepared argv/environment; direct exact selectors, one libtest thread, oracle pool1, 240-second per invocation deadline, zero retries and no new 20-GiB cap.',
   'mounts':'Fresh work/record destinations; input parent, program, exec bridge and capsule read-only; only owned work/home binds writable, plus private tmp/proc/dev mounts. Exact runtime store packages only; no host-root/home/config/credential, compiler, donor-cache or Nix-socket mount. Network unshared; original four unlicensed mount probes show only loopback.',
   'original_validation':'18 r1 pure control groups, including the original 23 negative finite guards and four mount/four loader probes, plus three final-r1-source simulated durable failures. These are retained evidence, not newly run by this reviewer. Initial/final/r2 lifecycle function ASTs agree.',
   'private_boundary':'All four private environment names are absent from every actual native preexec environment and pure subprocess plan. Controller alone obtains the license, preserving a nonempty exported value or using its current-user config fallback; anonymous stdin pipe goes to the immediate native exec bridge. No secret/config lookup, secret hashing or secret-value scan was performed in this review. Source/record inspection grants no adversarial secrecy theorem.',
   'fixture_disclosure':'Original r1 pure_errors.py uses a public stand-in environment value only with native launch replaced. This review ran no license fixture and imported no controller/module top-level code.',
   'durability':'Original outcomes and streams remain protected outside writable child mounts, including the frozen failed invocation3 without an original after record. Every actual result reports no remaining descendants; candidate cleanup records are retained.',
   'limit':'Cooperative lifecycle and finite capture containment only; no general future-runner or adversarial TOCTOU approval.'},
 'fresh_vs_candidate_validation':{
   'candidate_r2_repair_controls':{'positive':4,'negative':31,'native':0},
   'reviewer_saved_invocation_replays':12,'reviewer_independent_failed_Summary_pairs':6,'reviewer_independent_original_Rows':3,
   'reviewer_repaired_parser_receipt_negative_controls':43,'reviewer_original_parser_bug_reproductions':1,'reviewer_data_association_negative_controls':18,
   'reviewer_complete_world_joins':12,'reviewer_owner_records':1200,'reviewer_input_pairs':600,
   'reviewer_native_invocations':0,'reviewer_loader_mount_probes':0,'reviewer_builds':0,'reviewer_license_config_reads':0,
   'control_credit':'All reviewer rejection controls reach their named ValueError predicates. These are pure syntax/schema/identity/hash/association checks; no parser/setup failure is counted as scientific or native semantic rejection.',
   'method':'Standard-library scripts and explicitly extracted named pure AST functions only; no candidate module import, former helper entry point, scientific library, native capture, build/compiler, Nix, network or delegation.'},
 'reviewer_probe_history':[
   {'probe':'provenance_probe.py v1','error':'Reviewer initially required historical_None_Row_native_pass=true even on an Auto cell. Candidate correctly sets it false for Auto. Reviewer predicate corrected to (mode==None); v1 script/stdout/stderr retained. No candidate defect and no semantic credit.'},
   {'probe':'read-only source lookup','error':'An initial guessed target/xs-qed-eemm-identical/manual-source path did not exist; followed the actual qualification to sources/historical-observation-v2. No producer execution or mutation.'}],
 'remaining_gates':[
   'Independent factor-preserving exact mathematics for these complete new owner domains, including the 11 gggg decoded syntax differences per mode.',
   'Native raw/decoded admission, full captured symbol and port/flag finite-interface obligations, with separate approval before adopting interface explicitness/order or Row migration.',
   'Physical locality remains required: edge factors use only their own signed momentum and vertex factors incident momenta. No universal locality, denominator/integrand/UV or physical-equivalence theorem is supplied.',
   'No migration/literal change, general decoder, fullpool8/32, grouped, WP3, global integration, source46 expansion or publication approval.',
   'Old historical grouped/source gaps stay unapproved; fresh case association cannot reconstruct unrecorded past provenance.'],
 'process_and_retention':{
   'writes':str(O),'ignored_by':str(W/'.gitignore')+' (.local/)',
   'main_source_candidates_unchanged':True,'jj_ref_index_snapshot_operations':0,'network_operations':0,'delegations':0,
   'retained':'Review scripts, independent probe results, full before/after hash inventories, supplemental source bindings and the reviewer schema-error history are retained for audit. All frozen candidate inputs preserved.',
   'owned_processes_remaining':0,'workspace_or_ref_cleanup':'No workspace/ref removed or modified; this review is a new ignored artifact only.'}
}
with (O/'review.json').open('x') as f:json.dump(report,f,indent=2,ensure_ascii=False);f.write('\n')
files=sorted(p for p in O.rglob('*') if p.is_file() and p.name not in ['SHA256SUMS','SEAL.json'])
assert all(not p.is_symlink() for p in files)
manifest=''.join(sha(p)+'  '+str(p.relative_to(O))+'\n' for p in files)
with (O/'SHA256SUMS').open('x') as f:f.write(manifest)
seal={'schema':'independent-next-qcd-capture-provenance-review-seal/1','decision':report['decision'],
      'verdict':report['verdict'],'files':len(files),'review_sha256':sha(O/'review.json'),
      'manifest_sha256':sha(O/'SHA256SUMS'),'candidate_r2_handoff_sha256':sha(R2/'handoff.json'),
      'candidate_r2_SEAL_sha256':sha(R2/'SEAL.json'),'candidate_r2_manifest_sha256':sha(R2/'SHA256SUMS'),
      'native_invocations':0,'owned_processes_remaining':0}
with (O/'SEAL.json').open('x') as f:json.dump(seal,f,indent=2);f.write('\n')
for line in manifest.splitlines():
    expected,rel=line.split('  ',1);assert sha(O/rel)==expected
for p in O.rglob('*'):
    if p.is_file():p.chmod(0o444)
O.chmod(0o555)
print(json.dumps({'review':str(O/'review.json'),'review_sha256':seal['review_sha256'],'seal':str(O/'SEAL.json'),'seal_sha256':sha(O/'SEAL.json'),'manifest_sha256':seal['manifest_sha256'],'files':len(files),'decision':'approve'}))
