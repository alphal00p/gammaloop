import sys,datetime
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import common as C
rows=C.load(C.O/'all12-Rows-independent.json');proof=C.load(C.O/'proof-review-independent.json');patch=C.load(C.O/'patch-scope-independent.json');integrity=C.load(C.O/'integrity-final.json');initial=C.load(C.O/'integrity-before.json');refs=C.load(C.O/'refs-final-disclosure.json');controls=C.load(C.O/'fresh-controls.json');science=C.load(C.S/'review.json');handoff=C.load(C.P/'handoff.json')
per=[]
for case in C.CASES:
 h,c=(next(x['Row'] for x in rows if (x['case'],x['mode'],x['era'])==(case,'None',era)) for era in ('historical','current'))
 accepted=next(x for x in science['per_case_mode_verdicts'] if (x['case'],x['Classes'])==(case,'None'))
 per.append({'case':case,'field':'numerators','decision':'APPROVE','old':f"0x{h['numerators']:016x}",'new':f"0x{c['numerators']:016x}",'before_Row':h,'after_Row':c,'Classes':'None','Grouping':'None','captured_pool':1,'physical_argument':{'complete_owners':accepted['owner_relations'],'whole_tree_identities':accepted['complete_tree_identities'],'even_f_cycle_relations':accepted['even_f_cycle_relations'],'changed_owner_lowerings':accepted['changed_owner_lowerings'],'unchanged_imported_input_source_arguments':accepted['unchanged_identical_input_source_arguments'],'basis':'Accepted independent full factorized expression correctness and named-axis transport attached to exact pristine serialized digest owners, with full physical inputs, source factors, memberships, definitions and raw/cooked/decoded capture premises. All other nine fields independently reconstructed unchanged.'},'remaining_preapplication_gap':None,'postapplication_test_status':'NOT_RUN; future separately authorized owner required'})
conditions=[
'Only the three named numerators tokens in the exact832 source blob and exact proposed patch are approved. All149 Rows,13 prior full migrated Rows,118 ordered source scopes and every other byte must remain exact. xs_qcd_gggg_sym_L1 is excluded despite its shared old hash.',
'The accepted scientific premise is finite: the six separately captured case/None-or-Auto cells, GroupingNone, pool1. OnlyNone receives the original literal Row assertion. Auto remains independently captured/reconstructed and supplies no separate literal assertion.',
'Full exact model, request, options, Run, graph128-bit string IDs and full canonical keys, group source indices/members/ratios, source factors, cuts, signs, routing, incidence/half-edge slots, particle/mass/rule/coupling ownership remain bound. NumeratorIndependentSymmetryGrouping(2) applies to exactly4of5 sym gguu diagrams, not the entire case.',
'Both entering-operand and surviving-output maps are essential. The current typed8eb consumer imports both eras; this is not new c78 typed execution. Historical interface absence is not an empty interface. No historical positional identity or stale cached-plan interchangeability is authorized.',
'All600 full owner relations,96 unexpanded products, complete166/118 native pair partitions per kernel and exact definition/side orders remain dependencies.4268 inverse references are associations, not new obligations. Native568 pair-result rows and5baseline/16specific-negative evidence remain historical retained execution, never rerun here.',
'For210 changed owners retain both actual lowerings and the accepted420 layout audits;390 unchanged owners use exact imported expression/ports/flags plus the qualified deterministic source. Actual flags belong to reconstructed owners, not serialized-flag admission. Rank8 numerator plus rank4 projector leaves FOURcoad8 outputs. No scalar closure or component evaluation is inferred.',
'Historical2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94 uses Symbolica/Numerica c78cffc20f318df1fce50539a587cbe8ffc717c0 and registryGraphica3.0.1. Current pristine actualffa76ec420642c03349ab756d2d7720ab93c0057/archive8993e6c8346b84f6d91c7e65c4c4fc5ed83cf6d932e969bb7fc55cbdced20740/kernel8eb7c3729d6abefc916decc1bdf0c06536f8af80 is distinct from observerbuilt5ea149ffd88ae3afc4072a38c0871c1329be8424+qualified overlay and proofsource5ba04864ca506e259d9eef4a54f65acdd8323fc2. No producer is relabelled832/6f.',
'Whole symbol spellings/tags/attributes, uncaptured callback/userdata semantics, constructors/model/module initialization, original profile/compiler/kernel/store honesty and qualified observer/raw-cooked-decoded provenance remain explicit premises. Historicaldev-optim/currentCI/direct-rustc qualifications are not declared equivalent. eae/1a remains a failing baseline, never the original reference.',
'Finite locality912 records/1096momentum calls certify their exact dependencies only.41088visited factor subsets are finite pure enumeration, not native runs;456actual folds/720endpoint pairs are the actual retained typed observations. No general decoder/locality/UV/full-integrand/denominator-evenness/trace/component/evaluator/cross-section theorem is approved.',
'Cooperative protected-record/source/command/input/result consistency is the execution premise; no hostile or cryptographic attestation, new global interner observation or inner-PID/private-value evidence is invented.',
'INFO-1 remains nonblocking solely because the accepted independent scientific audit checked the exact actual integer plans, dispatch kinds and420lowering layouts. This migration review additionally checks all actual selected and fold integer plans and output maps. No general audit_native validator reuse is approved.26builtin bis no-ops earn zero rejection credit.',
'All unnamed rows, grouped modes, new sources/models/options/definitions, pool8/32, full371/WP3/M-rebuild/main/publication/global integration remain outside this approval. Metadata acceptance at832 is not a fullWP3 pass. Latest actual full-package result remains302passed/69failed/49ignored; no subtraction-based updated count is claimed.',
'External capture, source, helper, proof, actual execution and review dependencies remain required at their recorded paths/hashes; this review packet is not self-contained.'
]
findings={'decision':'APPROVE','blocking_findings':[],'smallest_remaining_preapplication_gap':None,'nonblocking_findings':[{'id':'MIG-INFO-1','title':'General auditor limitations remain excluded','disposition':'Accepted scientific review independently closed integer-plan/dispatch/lowering-layout checks for this exact sealed corpus. No additional native execution is needed. No broader validator qualification.'},{'id':'MIG-INFO-2','title':'Postapplication execution is a separate gate','disposition':'No tracked patch applied or runtime claimed. Future owner needs exact new producer archive, actual full-Row tests at1/8/32, fresh capture/proof joins where required and new integration review.'},{'id':'MIG-INFO-3','title':'Historical pending and failed-preparation records remain unchanged','disposition':'Scientific gate consumed through additive root acceptance225b2c..., while original preparation pending strings, eight authoring fixes, PREP1/PREP2 and SUP1/SUP2/SUP3 history, and original nested-Fraction capture error remain historical records.'}]}
review={
'schema':'independent-next-three-raw-qcd-numerator-migration-review/1','decision':'APPROVE','accepted_scope':'EXACT_THREE_NUMERATOR_LITERAL_MIGRATIONS_ONLY_PREAPPLICATION_WP9_DECISION','summary':'APPROVE all three exact numerator replacements. Full pristine ten-field Row reconstruction, exact proof-to-digest-owner bindings and immutable target applicability close the preapplication WP9 justification. No blocker or additional native proof requirement. This is not application, postapplication runtime or broader integration approval.',
'reviewer':{'role':'NEW independent migration review sub-agent','model':'Astra','fallback':'User-authorized after Opus rate limits','candidate_author':False,'delegated':False,'distinct_from_scientific_and_prior_row_reviews':True},
'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'workspace':str(C.W),'review_directory':str(C.O),'author_identity_already_settled':'Lucien Huber <im@lcnbr.ch>; no VCS identity/change operation',
'per_row_verdicts':per,'blocking_findings':[],'smallest_remaining_preapplication_gap':None,'additional_native_proof_requirement':None,
'candidate':{'path':str(C.P),'SEAL_sha256':'b4e75df006b68563219947ecbdef5cb8fd289f9b4c7ba840d79c0ce599377324','SHA256SUMS_sha256':'7c8c48341828a16d00b4f692f6ce3941d35b5679e73eac2fd4176f139942b503','manifest_files':122,'handoff_sha256':'1b69607db0184093de3b16d32541d6f7a5a3f00968b85694685a6d17fb1506b3','complete_input_manifest_sha256':'600ced69f0172b82a4043aa3d2a255e270c9782a79726f4d745ab93d3cc40e6e','patch_sha256':patch['patch_sha256'],'target':handoff['target'],'before_sha256':'66f01bcc855552ae888f1fb7c5bb158afbd448492b4e27a26e2bda40f09b07d2','after_sha256':'17e0abe46c76c6b6922e5b33782925020e1412e1a6ed8cd8dc7e23de438d4aa2','applied':False},
'additive_scientific_join':{'root_acceptance':str(C.F/'root-next-qcd-scientific-interface-acceptance-20261031.json'),'root_acceptance_sha256':'225b2c244ba02f16ca3ee73a59d23df3fff12204deef44f3aa405397247a0ea7','review':str(C.S/'review.json'),'review_sha256':'6d892bd9263a7c093bba5fd372fb4663a5a4e81dfaad4bfc709bbba052b96d21','SEAL_sha256':'888b82c478697e917cb243587653cc16d09a484f299a671c5ffd8981e58c1d99','manifest_files':57,'disposition':'Completed exact six-cell finite mathematics and named-consumer acceptance consumed additively; no original PENDING string altered or treated as a remaining blocker. The scientific review itself did not grant migration approval; this NEW decision does.','all_original_source_capture_definition_premises':science['source_capture_and_definition_premises'],'original_scientific_restrictions_verbatim_with_original_scope':science['restrictions']},
'WP9':{'path':science['WP9']['path'],'sha256':science['WP9']['sha256'],'rule':science['WP9']['rule'],'reasoning':'The historicalNone fullRows match original literals; pristine currentNone fullRows differ solely in numerators. Full physical graph domains, source factors and memberships are equal, and independently accepted whole-expression congruence/even-f-cycle mathematics plus finite named operand/output-axis covariance attach to every exact raw/cooked-decoded/public owner serialized into the digest. Exact source invariance transports this finite scientific justification to832 without claiming a new832 runtime. Thus this is a justified representation-sensitive expectation change, not re-recording a failed actual hash. Commit-only-passing-expectations remains a future postapplication gate, not asserted complete here.','global_completion':False},
'pipeline_comparison':science['pipeline_comparison'],
'independent_reconstruction':{'full_ten_field_Rows':12,'eras':['historical','current'],'modes':['None','Auto'],'literal_assertion_mode':'None','implementation':'Fresh stdlib digest/factor/Row decoder implemented from immutable Rust census; no candidate derive() execution and no failed assertion hash as input. Pure owner/parser helpers are inspected private byte-identical copies with stdlib imports only.','encoding':'FNV1a wrappingu64; exact LEu64/u128; canonical serialized JSON, text0xff/nameNone0xfe; full key/frame/numerator/interface vertex-then-edge encoding; ordered named/group folds; sorted digest multisets; exact rational source-factor evaluation.','complete_capture_joins':'All twelve full worlds directly joined to pristine and observer public JSON/source companions, exact census report, complete mode-event ledgers and raw owner data. Historical missing live.threads remains omitted and actual Run binds Some(1). No ordinal-only or rounded-ID matching.','report':'all12-Rows-independent.json'},
'proof_payload_binding':{'report':'proof-review-independent.json','owner_bindings':'proof-owner-join-independent.json','fold_bindings':'fold-output-join-independent.json','method':'All600 native input/actual owners rejoined to exact pristine serialized expression and explicit-or-absent interface segments; complete raw/decoded/port/flag data, model/request/options/Run, full keys/IDs, source factor/rule/coupling membership and full definition maps checked. Both actual operand and canonical output maps satisfy the complete commuting square.','local_relations':{'whole_tree_literal':578,'nonliteral_even_f_cycle':22,'unexpanded_products':96},'native_math_partitions':proof['native_math_partitions'],'native_rows_total':568,'inverse_references':4268,'selected_integer_plan_pairs':1890,'actual_fold_steps':456,'actual_fold_endpoint_pairs':720,'changed_owner_lowerings':210,'accepted_lowering_layouts_bound_to_same_actual_owners':420,'unchanged_equal_input_source_arguments':390,'new_native_credit':0,'scientific_proof_status':'Accepted independent finite physical mathematics retained, not inferred from counts/weights/hash consistency.'},
'source_applicability':{'actual_pristine':C.ACT,'target':C.REV,'relevant_code':C.CODE,'immutable_tree':'c9de5c3fd962bd5d7c9375ab1f86f75699e8a69b','source_objects_each_revision':1521,'method':'Read and independently hash every immutable object across complete crates/assets/root manifests/lock domain; match mode/type/blob and file inventories for actualffa/832/6f. Compare all repository deltas separately. Recheck16 frozen source files directly against Git objects.','only_generator_delta':'Exactly one unrelated previously approved xs_qed_eemm_identical_L2.numerators token; replacing that complete row in target with actualffa row reproduces every main.rs byte. All other1520 objects identical.832 and6f full relevant domains identical.','outside_domain_deltas':'1193 docs/evidence paths and one notebook .pyi stub; complete paths and object IDs retained. All other root source/test/build configuration unchanged.','producer_relabelled':False,'new_target_compile_or_runtime':False,'reports':['source-objects-independent.json','all-target-deltas-independent.json','ordered-scopes-independent.json']},
'patch_scope':patch,
'validation':{'fresh_controls':{'total':controls['total'],'positive':controls['accepted'],'intended_rejections':controls['rejected'],'report':'fresh-controls.json','credit':'Pure byte/identity/definition/map controls only; no native or scientific negative credit.'},'before_unique_files':initial['unique_files'],'final_unique_files':integrity['unique_files'],'all_initial_inputs_unchanged':True,'integrity_failures':[],'candidate_external_files_rehashed':3869,'previously_unread_source_configs_rehashed':integrity['source_configs_freshly_rehashed'],'supplemental_scientific_source_and_association_files_rehashed':221,'Git_refs_before_after_equal':refs['equal'],'refs_action':'Read-only disclosure; never restore or mutate concurrent work.','reviewer_errors':'errors.json; interpreter lookup, initial integer-aut assumption and overbroad schema discovery retained, no semantic credit.','candidate_authoring_history':'Eight errors/fixes retained unchanged in candidate/errors.json, with original failed scripts/streams; no native failure rewritten.','Typst':'UNRENDERED; no compiler invoked'},
'historical_reviews':'accepted-source-audits-rehashed.json retains exact capture/R1 PREP/R2 PREP+supervisor/R3 supervisor review hashes and original decisions; later native execution is separately authorized by650f597b... and accepted scientifically by6d892bd9... . No retroactive overwrite.',
'conditions_and_restrictions':conditions,'nonblocking_findings':findings['nonblocking_findings'],
'future_authorized_owner_boundary':{'preapplication_missing_proof':None,'required_later':['Apply only the exact three-token patch to its exact immutable before bytes in the separately authorized implementation workspace.','Create a new postapplication producer archive/source qualification and run actual full-Row tests for all three cases at pools1,8,32.','Attach fresh capture/proof joins wherever the ensuing source, build, payload or interface boundary requires them; preserve all inherited finite/source conditions.','Obtain a NEW integration review of the actual postapplication evidence before any broader checkpoint or commit/publication claim.'],'approved_runtime_passes_now':False,'proof_native21_rerun_required':False},
'latest_full_package_result':{'passed':302,'failed':69,'ignored':49,'updated_full_run':False,'inventory_note':'All unnamed rows excluded regardless of the historical other60 phrase; no inferred remaining-failure count.'},
'activity':{'new_native_or_loader_executions':0,'scientific_library_imports':0,'private_user_config_or_license_reads':0,'Nix_build_compiler_network_commands':0,'tracked_edits':0,'VCS_index_ref_snapshot_mutations':0,'writes_outside_new_review_directory':0,'delegations':0,'publication':False,'patch_application':'Only frozen private copy within review directory','owned_processes_remaining':0},
'guidance':{'read':['main/AGENTS.md','main/CONTRIBUTING.typ','selected AGENTS.md','selected CONTRIBUTING.typ','accepted integration plan WP9','applicable source and prior review guidance'],'local':'No deeper applicable AGENTS under selected generator/source/review paths; .local/ ignore rule read directly. Explicit inspection-only task overrides implementation/build/license/VCS workflows. No Rust source implementation edit or style departure.'},
'retention':'Retain this NEW ignored review directory, independent scripts/controls, original failed reviewer probes, source copies, exact private patch result and hash manifests for root continuation. All original candidate/capture/proof/review inputs retained at original paths. No new ref created under explicit prohibition; no workspace removed; no owned process remains. External dependencies required; not self-contained.'}
C.save('findings.json',findings);C.save('review.json',review)
C.save('external-input-bindings.json',{'files':integrity['checks'],'count':integrity['unique_files'],'all_initial_inputs_unchanged':True,'immutable_source_objects':'source-objects-independent.json','self_contained':False,'original_candidate_and_scientific_seals':integrity['pinned_files']})
C.put('findings.typ',b'''= Exact three raw QCD numerator migration review

Decision: *APPROVE* each of the three exact numerator literal replacements in
`review.json`. No blocking finding or remaining preapplication proof gap was
found. This is a new migration review, distinct from the accepted scientific
review and earlier row reviews. No tracked source was changed.

== Decision and WP9 basis

- xs_qcd_gguu_sym_L1: b9fdc23ed056fa0c to e4b209849953f1fa.
- xs_qcd_gguu_L1: 5edf8e5f7a11eeb3 to 665e13a4c88da679.
- xs_qcd_gggg_L1: 2628c62dc373ac91 to 96e5abcebda25959.

The accepted independent scientific review supplies full finite physical
correctness and named-consumer compatibility under its explicit premises. This
review attaches that argument to every exact pristine payload owner feeding the
serialized digest, reconstructs all ten Row fields, and independently verifies
applicability to immutable832. The conjunction satisfies WP9 for these three
representation-sensitive expectation changes. Counts, weights and observed
failed hashes alone do not establish the decision.

All12 Rows were reconstructed independently from pristine None/Auto payloads
for both eras using the actual immutable Rust census encoding. Each historical
None Row matches its original literal; each current None Row differs solely in
numerators. Auto has its own live/request/event/source/owner associations and is
not granted a literal Row assertion. Explicit historical interface absence is
preserved. Source factor NumeratorIndependentSymmetryGrouping(2) belongs to
four of five sym gguu diagrams only.

The full600 owner associations retain model/request/options/Run, exact string
IDs and complete keys, routing/incidence/half-edge slots, memberships, source
factors, masses, couplings and whole symbol definitions. Both operand AND output
maps are checked. All1890 selected integer plans and456 actual folds with720
endpoint pairs agree with the accepted data. Full166/118 native partitions per
kernel, formal witnesses and definition-side order reconstruct from the actual
owners;4268 inverse references supply no extra obligation credit. All96 products
remain unexpanded. The scientific proof covers578 literal relations and22 even
f-cycle relations. Rank8 numerator plus rank4 projector leaves four coad8 axes.

== Exact target and patch

Target83276dc047815d1dcc7f003aeb512c204b26758e has tree
c9de5c3fd962bd5d7c9375ab1f86f75699e8a69b. The sole patch path is
crates/feynkit-generator/tests/generation_oracle/main.rs, before blob
7283198bf0152c82c1cd5f0173703cddfa7f5c0f and after blob
cc4505a6c4f6fc1687c075ed79e154a8e19f05a6.

Every one of1521 immutable source objects in the complete crates/assets/root
manifests/lock domain was read and checked. Actual pristineffa and target832
differ solely by the unrelated already approved grouped eemm numerator token;
all other1520 objects match.832 and relevantcode6f match throughout that domain.
Complete outside-domain changes are1193 docs/evidence paths and one notebook
stub. This is source applicability, not a relabelled producer or new832 runtime.

The candidate patch applies to a reviewer-private frozen copy with zero offset
and zero fuzz and produces the exact after blob. Exactly3 numerator fields
change;146 complete Rows and all other source bytes remain unchanged. All13
prior full migrated Rows and118 ordered source bindings are preserved. The
same-old-hash xs_qcd_gggg_sym_L1 sibling remains untouched.

== Preserved premises and boundaries

Historical2b/c78 Symbolica-Numerica with registryGraphica3.0.1, current
pristineffa/archive8993/kernel8eb, observerbuilt5ea plus overlay and proofsource5ba
retain their separate identities. Source/profile/compiler/store, initialization,
whole-definition/callback/userdata and raw/cooked/decoded capture premises stay
explicit.210 changed owners retain actual lowerings on both representations;
390 unchanged owners use exact imported-input equality and deterministic source.
Actual owner flags do not constitute serialized-flag admission. Historical
positional identity and stale cached-plan interchangeability remain excluded.

The original auditor's integer-plan/layout limitations are nonblocking for this
fixed corpus because the accepted scientific review independently checked the
actual data. No general validator reuse is approved.26 built-in bis no-ops earn
no rejection credit. Pure controls and this review run no native mathematics.
All historical preparation errors, PREP and SUP findings, the original capture
Fraction framing error and additive revalidation remain distinct, unchanged
records. Original PENDING text is resolved additively through the root scientific
acceptance, not rewritten.

All26 fresh controls pass:4 valid positives and22 intended rejections, covering
wrong case/mode/live/model/source/owner/full-definition/operand/output/plan data,
wrong literals, sibling/prior rows and extra patch scope. The reviewer initial
aut-schema mistake and interpreter/discovery errors are retained without credit.
Final input rehash passes; exact details and ref disclosure are in review.json.
Two repository source configurations were freshly hashed; no private user
configuration or license was read. No owned process remains.

A future separately authorized owner still needs a new postapplication archive,
actual full-Row tests at pools1/8/32, fresh capture/proof joins where required,
and new integration review. No additional native21 proof rerun is requested.
No postapplication runtime, globalWP3/full371, other row, grouped, broader pool,
main or publication pass is approved. The latest actual full package result
remains302passed/69failed/49ignored, with no new full run. The finite scientific
argument is not a universal decoder, locality/UV, full-integrand or component
physics theorem. External dependencies remain required; this is not a
self-contained packet. Retain the ignored review artifacts for root review.
''')
print('Wrote APPROVE per-row migration verdicts and retained limitations.')
