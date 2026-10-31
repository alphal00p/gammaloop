"""Write reviewable data/proposals only into this ignored packet."""
from pathlib import Path
exec(Path(__file__).with_name('evidence_io.py').read_text())
import difflib, re, tomllib

CODE='2cfb52b28d99f92305a9d571bbe3911cf04b8a19'
EVIDENCE='457e361d7d39e0a3daae801651f8c1dfabbf47aa'
GLOBAL='fdff0d18953348b84a0f1b7e6ab5ef703a5090de'
KERNEL='8eb7c3729d6abefc916decc1bdf0c06536f8af80'
URI='git+file:///shared/lcnbr/dev/gammaloop/main?rev='+CODE+'&shallow=1&allRefs=1'
SOURCE='/nix/store/3sj3nw06s3vhl3j12wfp5gc3r1h6f0k3-source'
PLAN='docs/architecture/generation-integration-plan.typ'
EXECUTION='docs/architecture/evidence/generation-integration/execution.json'
NOTE='docs/architecture/generation-feynkit-community-refresh.typ'
REGISTRY='docs/developers.toml'
def emit(name,value): (P/name).write_text(json.dumps(value,indent=2)+'\n')
def gate(name,status,requirement,evidence,limits):
 return dict(gate=name,status=status,accepted_requirement=requirement,evidence=evidence,limits=limits)

decisions=data(P/'decision-index.json')
matrix=[
 gate('NEXT bounded combined refresh integration','NEW combined review unresolved',
      'Each bounded step needs a NEW independent review; the new conjunction must be reviewed on an exact frozen candidate. Exit criteria are necessary, not sufficient; reviewers may add checks.',
      [PLAN+':4-13',PLAN+':137-140',NOTE+':247-251','decision-index.json'],
      'Constituent APPROVEs do not approve their conjunction. No accepted text makes all 69 old numerator rows a prerequisite for every bounded refresh integration. The next review must state its bounded scope and retain production-core/evaluator, full WP3 and publication exclusions. No bookmark movement is proposed here.'),
 gate('WP0','accepted','Surviving relevant evidence and refs durable; unavailable raw evidence classified rather than silently certified.',[PLAN+':142-154',EXECUTION+':packages.WP0'],'Archiving is not runtime or physics certification.'),
 gate('WP1 / B2','accepted finite original-pair theorem; R/F rerun obligation remains','Original GL028/GL034 exact 4D SU(3) N28=(7/2)N34 certificate, zero-sample and sampled production controls; explicitly rerun ignored original-pair regression on rebuilt and final tips.',[PLAN+':155-189','docs/architecture/evidence/generation-integration/wp3/exit-matrix/receipt.json:matrix.b2'],'Old b9 six-pass replay is 1a, not a current8eb execution. Full generator pool1 does not execute the Rust-ignored original-pair test. No new mathematics or bisect warranted.'),
 gate('WP2 / refresh dependency base','accepted historical base; bounded refresh source accepted','Exact base/dependency manifest reviewed; source replay and new kernel pinned explicitly.',[PLAN+':191-206',EXECUTION+':packages.WP2','decision-index.json:source R4'],
      '3f0df2a8 is historical accepted base. Refresh targets 9102f885722e856e30227ed46e685e5c465905f7, kernel8eb. R4 source approval is not whole-workspace runtime compatibility.'),
 gate('WP3 / M-rebuild','incomplete','Conflict-free formatted chain, deny-warnings Clippy, green fast lane; all fast oracle rows at pools1/8/32, SM3150 at1/8, four-loop4664/306 plus perturbation at1/8/32, a_ddx5 fresh processes, BRST, epem NLO, CoolQM; R-stage validation matrix and milestone recording.',[PLAN+':468-474',PLAN+':620-646',EXECUTION+':active_wp3_exit_matrix'],'Full current pool1 is 302pass69fail49ignored, not green. Twelve exact rows pass36 at1/8/32; this is not complete oracle coverage. No M-rebuild claimed. Full fast-green exit and WP9 failure classification must remain distinct from a bounded integration approval.'),
 gate('source46 integration','accepted cumulative at79dfeaab; stale pending wording','Approved colour-cycle owner integrated; exact producer/profile and consumer limits separately bound.',[EXECUTION+':active_slice46','../generation-integration-plan/.local/review-main46-round4-opus-rZyrnOS8/review.json#sha256=ee6e37ee207747d4f65c6a45b619412738af3248a65301b396768f27f5227e22'],
      '192 profile+2 default-filter controls,16 networks17484 permutations were old1a observations. The implementation need not be repeated. Host/full-generator review exclusions of source46 do not revoke its approval or supply fresh8eb standalone colour-test credit.'),
 gate('tensor layout / core repairs / XS-QED original integration','accepted bounded historical steps','Preserve reviewed tensor restoration, BTreeSet inference fix, API compile-time fixture, generated-evidence owner and eight prior migrations.',[EXECUTION+':active_tensor_layout',EXECUTION+':active_core_test_inference',EXECUTION+':active_core_docs','decision-index.json:original XS scopes A/B'],
      'Tensor202pass3ignore plus17doctests, core3exact and later2214fast, docs116 library tests are bound to old source/profile/kernel. Original b9 core failure remains immutable history. XS A and B each scoped; current globalfdff is the accepted original integration boundary, not approval of the8eb refresh.'),
 gate('current12 rows / full generator','accepted bounded current8eb evidence','Reuse accepted ffa producer,36 exact passes and15 capture joins, full420inventory371selected302pass69fail49ignored at pool1.', ['.local/feynkit-refresh/root-three-integration-host-acceptance-and-next-20261031.json','.local/feynkit-refresh/root-full-generator-pool1-acceptance-20261031.json','candidate-code-delta.json','code-ffa-delta.json'],
      'None alone asserts literal full Row; Auto payload fresh but historical counter premises remain.69 matching ten-field Summary pairs classify assertion history only; no physical theorem or migration authorization. No repeat merely for metadata or host stub outside generator inputs.'),
 gate('ordinary host / private RustRed','accepted bounded current8eb composition','Exact applied generated stub, corrected ordinary wheel installation/import/parity; separate private7c1 single-extension one-denominator validation host;8/8 fresh reviewer checks.', ['.local/feynkit-refresh/root-host-stub-rustred-acceptance-20261031.json','decision-index.json:host review'],
      'Ordinary host lacks IBPFamily: no ordinary QCD pass. Private RustRed adds symbolica/tracing_max_level_info, native runtime arities[1], genuine reduction targets[[2],[3]],depth2,residual[[1]],23 assertions. No production RustRed, full Community, campaign/reconstruction or multiple-denominator adoption; no full docs/widget/site claim. Wheel native bytes, loose fixed DSO and loader sections remain distinct.'),
 gate('grouped None / Auto and finite runner','accepted exact conditional aggregates/finite R3; migration separate','Preserve independently pinned None and Auto premises and exact finite runner maps.', ['.local/feynkit-refresh/root-grouped-none-acceptance-20261031.json','.local/feynkit-refresh/root-grouped-auto-acceptance-20261031.json','.local/feynkit-refresh/root-interface-runner-r3-acceptance-20261031.json'],
      'Grouped Identical literal remains old at457; migration review active elsewhere. No waiting for it is required to review frozen457. No generic runner/adapter, scalar/sign, other rows/pools, raw-cook, positional/full-integrand or maximality approval. Operand AND output-axis maps mandatory; None and Auto stay separate.'),
 gate('WP4 names','later; no exit claimed','Canonical GL names and name@id pins, drop Request classes/serde fallback, settled57/8 normalization; names stable across threads/modes/grouping and card selections/pin errors/cost validated.',[PLAN+':476-490'],'Do not add to the immediate bounded refresh checklist.'),
 gate('WP5 route','later; no exit claimed','RoutedSet and numerator stage, select before numerators; build-equivalence oracle, routed peak/Ir, signed gammabench adapter.',[PLAN+':492-499'],'No exit established by the bounded host or row reviews.'),
 gate('WP6 Classes::Auto','later; no exit claimed','Control hash variance; Auto/None <=1.005 Ir generation/full on named workloads; DHAT bytes<=+5%; heavy tier both modes.',[PLAN+':501-513'],'Default remains Auto. Gate failure calls for investigation or measured user decision, not weakening.'),
 gate('WP7 performance','later','Same-output-contract stage table versus stage1/v4; measure serial sampling/pools/cache hosting; ablation justifies kept optimizations.',[PLAN+':515-529'],'No pre-agreed target ratio. Wall time needs spread; no isolated timing credit from shared characterization.'),
 gate('WP8 memory','later; measurement before limit','Measure DHAT holders and peak ratios versus stage1 at tip, then propose justified limit.',[PLAN+':531-543',PLAN+':113-115'],'No20GiB cap adopted and no preauthorized1.2x ceiling.'),
 gate('WP9 correctness / snapshots / heavy','unresolved later closure','Every failure inherited/new/justified expectation change; recertification requires membership+weights or full expression equivalence plus independent correctness argument; heavy final tip.',[PLAN+':545-563'],'Equal counts/summed weights/Summary pairs are insufficient. Active threeQCD capture and grouped migration owners continue separately; no duplicate proof work.'),
 gate('WP10 milestones / ablation','later','M-rebuild and decision-stage/M-final registrations; >=21 balanced rounds, paired spread/counters; keep/drop table and local performance script.',[PLAN+':565-579'],'No milestone or performance approval in this packet.'),
 gate('WP11 notices / drawing','notices accepted historically; shipping gate persists','Keep Clarabel and required Clarabel/NetworkX notices in every shipped artifact; approved drawing layers only.',[PLAN+':581-589',EXECUTION+':packages.WP11-notices'],'Parallel lane does not block generation delivery; nothing ships without notices.'),
 gate('WP12 cleanup','later','Bounded listed minors and repository-owned leftovers.',[PLAN+':591-596'],'No unrelated implementation authorized by this preflight.'),
 gate('WP13 final','later','Final-tip full matrix, independent adversarial whole-stack review, just docs-check, crates.io release suitability check before dropping patch.',[PLAN+':598-602'],'Full site/WASM/SVG/wheels/Pyodide/publication packaging are not inferred from bounded installed-host checks. No network now.'),
 gate('WP14 / publication / historical cleanup','requires explicit user go','No push/fetch/PR/label without user go; before push every commit conflict-free and authored Lucien Huber, B1/B2 closed, WP13 complete, notices; CI enable/validate/upload/push/final-review sequence. Cleanup selective and authorized.',[PLAN+':665-680'],'Lucien Huber <im@lcnbr.ch> already settled. No general publication authorization. Preserve sources/bookmarks/caches and do not touch /common/build-cache.'),
]
emit('status-matrix.json',matrix)

# Recover precise historical selections as data rather than silently rewriting their source/profile.
historical=[]
for file in ['docs/architecture/evidence/generation-integration/wp3/prerequisites08-14.json','docs/architecture/evidence/generation-integration/wp3/units27-29/validation.json','docs/architecture/evidence/generation-integration/wp3/staged-combined-corrections.json']:
 v=data(file)
 def visit(x,path=''):
  if isinstance(x,dict):
   for k,y in x.items():
    if k in ['command','argv','commands']: historical.append({'source':file,'json_pointer':path+'/'+k,'value':y})
    elif k not in ['artifact_hashes','artifact_sha256','source_hashes']:visit(y,path+'/'+k)
  elif isinstance(x,list):
   for i,y in enumerate(x):visit(y,path+'/'+str(i))
 visit(v)
emit('original-selection-records.json',historical)

# A metadata proposal does not update any historical runtime identity or assert combined approval.
accepted=[{'review':x['review'],'sha256':x['actual_sha256'],'scope':x.get('scope'),'verdict':x.get('verdict')} for x in decisions]
current={
 'schema':'generation-refresh-current-status-proposal/1','status':'CONSTITUENT_BOUNDED_ACCEPTANCES; NEW_COMBINED_REVIEW_REQUIRED',
 'original_approved_global_integration':GLOBAL,'frozen_bounded_refresh_code':CODE,'frozen_bounded_refresh_evidence':EVIDENCE,
 'source_upstream':'9102f885722e856e30227ed46e685e5c465905f7','kernel':KERNEL,
 'author':'Lucien Huber <im@lcnbr.ch>','CI_enabled':False,'publication_authorized':False,
 'accepted_reviews':accepted,'full_generator_pool1':{'inventory':420,'selected':371,'passed':302,'failed':69,'ignored':49,'classification':'numerator-only asserted Summary history; not physics approval'},
 'row_scope':'12 literal None rows;36 exact passes at1/8/32;15 fresh capture joins; Auto physical/counter premises retained separately',
 'grouped_scope':'None and Auto conditional aggregates separately accepted; grouped literal unchanged at457; active migration review excluded',
 'host_scope':'ordinary corrected wheel plus private RustRed7c1/kernel8eb one-denominator validation host; no production RustRed/full Community',
 'finite_runner':'Exact r3 accepted separately; no general runner/native/math cross-credit',
 'current8eb_production_core_evaluator_status':'not yet runtime-qualified by the accepted evidence read here; exact core archive outputs absent in local store',
 'WP3':'incomplete','M_rebuild':'not claimed','scope_bindings_preserved':118,
 'historical_records':'Verbatim reviews and original receipts remain unchanged, including their then-pending language; current attribution is additive outside them.'}
old_execution=read(EXECUTION);new_execution=json.loads(old_execution)
updates={
 '/packages/WP3/status':'Sources through46 and original tensor/core/docs/XS integration accepted in their recorded old-kernel scopes; bounded8eb refresh constituents separately accepted. NEW combined refresh review required; WP3/M-rebuild incomplete.',
 '/packages/WP3/current_slice_review':'Historical source30 pending description superseded by recorded later bounded approvals. Current frozen refresh457/code2cf awaits NEW independent combined review.',
 '/packages/WP3/review':'See current_status_reconciliation_20261031 for exact bounded attributions; original receipts remain immutable. No current8eb whole-WP3 approval.',
 '/packages/WP3/parallel_lanes/color46/status':'Original source46 cumulatively APPROVED at79dfeaab by final independent Opus round4; no repeat implementation required. Current8eb standalone colour/evaluator qualification remains distinct.',
 '/active_xs_qed_l2/status':'Original fdff/eae bounded integration approved through independent scopes A/B and accepted globalfdff boundary;298/73/49 remains old1a history. Current8eb ffa pool1 is302/69/49; no whole-WP3 approval.',
}
patch=[]
for pointer,value in updates.items():
 node=new_execution;parts=pointer.strip('/').split('/')
 for part in parts[:-1]:node=node[part]
 before=node[parts[-1]];patch += [{'op':'test','path':pointer,'value':before},{'op':'replace','path':pointer,'value':value}];node[parts[-1]]=value
patch.append({'op':'add','path':'/current_status_reconciliation_20261031','value':current});new_execution['current_status_reconciliation_20261031']=current
emit('execution-reconciliation.NOT-APPLIED.json-patch',patch)
emit('execution.NOT-APPLIED.json',new_execution)
addition='''
== Current bounded dispositions reconciled 2026-10-31

This additive status proposal remains unreviewed. The original six verbatim
reviews and their then-pending wording are immutable history. Metadata48d0 was
independently accepted for its attribution only (review SHA256567fba9cfa9287478d5c1531fd46eac49c79f4cd6fbcc878c50bd364cf279342).
No constituent review approves the current conjunction.

Frozen code2cfb52b28d99f92305a9d571bbe3911cf04b8a19 and evidence457e361d7d39e0a3daae801651f8c1dfabbf47aa contain twelve bounded
literal migrations and the one generated notebook stub. SourceR4/sharedScopeA,
the twelve-row integration, full pool1 characterization, corrected ordinary
distribution and private RustRed validation host have separately bounded
acceptances. Root records under .local/feynkit-refresh/ bind the exact reviews;
the preflight decision index records their rehashed immutable identities.
Full characterization is420 inventory,371 selected,302 passes,69 numerator-only
assertion failures and49 ignored at pool1. Matching Summary pairs supply no
physics or migration credit. The36 exact row passes at1/8/32 assert only None
literal full Rows; Auto captures and historical counter premises stay separate.

None and Auto grouped conditional aggregates are separately accepted; the old
grouped literal and active migration review are outside frozen457. Exact finite
runnerR3 is accepted separately without generic runner or adapter approval.
The corrected ordinary host lacks IBPFamily. Private RustRed7c1 is a separate
single-extension, one-denominator validation host; its symbolica tracing_max_level_info
feature differs from the ordinary host. This grants no production RustRed,
full Community, campaign, reconstruction or multiple-denominator adoption.

Original source46 is cumulatively APPROVED at79dfeaab77958064bf802bc88db23eb86ac2f8ce
(final review ee6e37ee207747d4f65c6a45b619412738af3248a65301b396768f27f5227e22).
Old receipt exclusions do not require repeating its implementation. Original
tensor/core/docs/ghost-weight/XS approvals retain their source/kernel limits;
they supply no fresh8eb production-core/evaluator execution. Current pure scope
inspection validates all118 bindings but is not compiled/docs-runtime credit.

A NEW independent combined review on an exact frozen metadata candidate is
still required before any refreshed integration bookmark advance. The accepted
bounded-step workflow does not require all69 old rows to migrate before every
bounded refresh integration. WP3's full exit, M-rebuild and later names, route,
Classes, correctness, performance, memory, WP9, final packaging and publication
gates remain separately open. Measure memory before proposing a limit; no20GiB
cap or1.2x ceiling was adopted. Author Lucien Huber is settled. Publication and
destructive historical cleanup still need the user's explicit go.
'''
note=read(NOTE);proposed_note=note+addition
registry=read(REGISTRY);start=registry.index('id = "generation-feynkit-community-refresh"');end=registry.index('[[section.note]]',start)
section=registry[start:end]
section=re.sub(r'^summary = .*$', 'summary = "Separately accepted source, provenance, finite migration, pool1 and host evidence; new combined refresh review remains required."',section,flags=re.M)
section=re.sub(r'^status = .*$', 'status = "Bounded constituents accepted; NEW combined review and broader WP3 gates remain open"',section,flags=re.M)
section=section.replace('captured_at = "2026-10-30"','captured_at = "2026-10-31"').replace('evidence_revision = "cacd5182704082df929a066548e3501be22c9533"','evidence_revision = "'+EVIDENCE+'"')
proposed_registry=registry[:start]+section+registry[end:]
def scopes(v):return [(n['id'],s) for a in v['section'] for n in a.get('note',[]) for s in n.get('scope',[])]
assert scopes(tomllib.loads(registry))==scopes(tomllib.loads(proposed_registry))
assert 'review_ref =' not in section and 'reviewed_at =' not in section
(P/'canonical-note-addition.NOT-APPLIED.typ').write_text(addition)
diff=''.join(difflib.unified_diff(note.splitlines(True),proposed_note.splitlines(True),fromfile='a/'+NOTE,tofile='b/'+NOTE))+''.join(difflib.unified_diff(registry.splitlines(True),proposed_registry.splitlines(True),fromfile='a/'+REGISTRY,tofile='b/'+REGISTRY))
(P/'metadata.NOT-APPLIED.patch').write_text(diff)
emit('metadata-proposal-guards.json',{'NOT_APPLIED':True,'base_evidence':EVIDENCE,'base_file_sha256':{f:hashlib.sha256(read(f).encode()).hexdigest() for f in [EXECUTION,NOTE,REGISTRY]},'scope_bindings_before':118,'scope_bindings_after':118,'scope_records_order_and_values_identical':True,'no_new_approval_fields':True,'no_combined_approval_claim':True,'publication_gate_unchanged':True,'apply_only_after_coordinator_freezes_separate_metadata_candidate':True,'rebase_if_concurrent_state_changes':True})
print(json.dumps({'status_rows':len(matrix),'original_selection_records':len(historical),'metadata_proposal_applied':False,'scope_bindings_preserved':118}))
