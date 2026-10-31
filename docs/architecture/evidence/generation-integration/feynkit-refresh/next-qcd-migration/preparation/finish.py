"""Close pure audits, preserve error accounting, rehash inputs, then seal closed files."""
import ast
import collections
import json
from pathlib import Path
import context as C
import reconstruct as R
import controls

def main():
    before=C.read(C.O/'before.rs');after=C.read(C.O/'after.NOT-APPLIED.rs')
    candidate=C.load(C.O/'candidate-index.json');audit=C.load(C.O/'byte-scope-audit.json')
    controls.scope(before,after,candidate)
    patch=C.read(C.O/'proposed-not-applied.patch')
    C.need(R.apply_strict(before.decode(),patch.decode()).encode()==after,'final-private-patch-application')
    # The real owner hashes complete files unless an explicit end anchor scopes them.
    ordered=C.load(C.O/'ordered118-source-bindings.json');after_scopes=[]
    for scope in ordered['ordered_scopes']:
        text=after.decode() if scope['path']==C.MAIN else C.git('show',C.TARGET+':'+scope['path']).decode()
        C.need(text.count(scope['anchor'])==1,'final-scope-start')
        if 'end_anchor' in scope:
            C.need(text.count(scope['end_anchor'])==1,'final-scope-end')
            text=text[text.index(scope['anchor']):text.index(scope['end_anchor'])]
        actual=C.blob(text.encode());C.need(actual==scope['digest']==scope['actual_digest'],'all118-after-scope-digests',scope)
        after_scopes.append({'position':scope['position'],'note':scope['note'],'path':scope['path'],'before_digest':scope['digest'],'after_digest':actual,'byte_scope_preserved':True})
    C.save('all118-scopes-after-candidate.json',{'count':len(after_scopes),'all_preserved':True,'ordered_bindings':after_scopes,'registry_edited':False,'docs_rendered':False})
    py=[]
    for p in sorted(C.O.rglob('*.py')):
        ast.parse(p.read_text(),filename=str(p));py.append({'path':str(p.relative_to(C.O)),'sha256':C.sha(p.read_bytes()),'parse':'PASS; pure AST parse, no compiler/native helper'})
    C.save('script-source-validation.json',py)
    errors=[
        {'stage':'read-only discovery','error':'Guessed capture meta/handoff.json and capsule current_pristine paths did not exist. Correct actual paths used; no data changed.','evidence':'Original terminal output; no output files created by those reads.'},
        {'stage':'interpreter discovery','error':'python3 was absent from shell PATH. Existing stdlib-only Python3.13.12 executable selected directly; no Nix/build/install command.','evidence':'Original terminal output.'},
        {'stage':'read-only schema inspection','error':'Prepared capture world has no model_bindings key. Exact owner/model predicates reconstruct that field separately; original schema unchanged.','evidence':'Original terminal KeyError.'},
        {'stage':'bootstrap attempt1','error':'Mistook100 architecture notes for the118 ordered code-scope bindings. Corrected to nested note.scope owner algorithm; all118 original and candidate digests pass.','evidence':'bootstrap.attempt1.py + bootstrap.attempt1.stderr'},
        {'stage':'reconstruction attempt1','error':'Historical live capture omits threads. Preserved omission; actual captured Run(threads Some1) binds pool. No absent-to-empty inference.','evidence':'reconstruct.attempt1.py + reconstruct.attempt1.stderr'},
        {'stage':'applicability attempt1','error':'Assumed four prior tokens differed from actualffa. Immutable comparison proves exactlyONE prior grouped token fromffa to832. Four-token observer/proof source qualification remains separately attributed. Corrected exact one-token fence passes.','evidence':'applicability.attempt1.py + applicability.attempt1.stderr'},
        {'stage':'association attempt1','error':'Unexpanded aggregate obligations use direct case/mode/graph identity and memberNone, unlike local owner lists. Added exact branch; all4268 retained member associations now bind.','evidence':'associate.attempt1.py + associate.attempt1.stderr'},
        {'stage':'association expansion attempt2','error':'All new owner-chain checks passed; exclusive existing consumed-input output refused a changed dependency inventory. Original inventory retained, new complete-chain inventory written under a distinct name.','evidence':'associate.attempt2.py + associate.attempt2.stderr'},
    ]
    C.save('errors.json',{'analysis_and_preparation_only':errors,'native_or_scientific_library_runs':0,'scientific_negative_credit':0,'scope_mismatch_in_final_candidate':False,'original_capture_controller_error':'Original nested Fraction controller error retained in original records; additive invocation3 revalidation stays separate. Neither rewritten nor upgraded to native pass.','external_native_setup_errors':C.ref(C.N/'source-and-setup-errors.json')})
    # Every external file read by the pure reconstruction is included. All old
    # manifest bindings are retained even where source-config reads are forbidden.
    originals=C.load(C.O/'original-inputs-before.json');consumed={}
    for p in sorted(C.O.glob('*consumed-inputs*.json')):
        for name,r in C.load(p).items():
            if not Path(name).is_relative_to(C.O):
                C.need(name not in consumed or consumed[name]==r,'complete-input-consistency',name)
                consumed[name]=r
    for name,r in C.used.items():
        if not Path(name).is_relative_to(C.O): consumed[name]=r
    expected=dict(originals)
    for name,r in consumed.items():
        C.need(name not in expected or expected[name]==r,'consumed-original-binding',name)
        expected[name]=r
    actual={name:C.fingerprint(name) for name in sorted(expected)}
    C.need(actual==expected,'original-input-before-after-invariance')
    C.save('original-inputs-after.json',{name:actual[name] for name in originals})
    C.save('complete-input-manifest.json',{'external_files':expected,'external_file_count':len(expected),'fresh_rehash_before_after':True,'original_input_files':len(originals),'original_inputs_preserved':True,'source_config_references_not_read':C.load(C.O/'recorded-unread-config-dependencies.json'),'immutable_Git_sources':'frozen-source-bindings.json + complete-source-object-bindings.json + ordered118-source-bindings.json','internal_helper_sources':'copied-pure-source-bindings.json; all scripts and failed versions included in packet seal','self_contained':False,'retention':'Original capture/proof/preparation/review/native/helper/runtime/source Git dependencies stay at their exact paths; packet does not copy or substitute their execution.'})
    now={'jj_ignore_working_copy':C.jj_state(),'git_refs':C.git('for-each-ref','--format=%(refname) %(objectname)').decode()}
    earlier=C.load(C.O/'refs-before.json')
    C.save('refs-after-and-disclosure.json',{'before':earlier,'after':now,'concurrent_ref_change_observed':earlier!=now,'changes_restored_or_reverted':False,'this_task_VCS_mutations':0,'target_binding':'Always immutable832; live bookmark and live@ never select patch base.'})
    warnings=[
        'SCIENTIFIC_INTERFACE_REVIEW_PENDING: running Astra review not accessed, consumed, waited for or polled; no scientific verdict or root scientific acceptance invented.',
        'NEW_MIGRATION_REVIEW_REQUIRED: candidate supplies no migration approval or application authorization. Root must join the later completed scientific verdict to this exact seal, then assign a NEW independent migration review.',
        'NOT_APPLIED: no source/expectation/ref/index/JJ/snapshot edits. Private patch test only. No new target producer build, native runtime pass, full371, WP3/M-rebuild/main/publication or broader combined approval.',
        'Typed8eb imports both eras; no new historical c78 consumer. Both operand AND output maps are mandatory; numerator rank8 plus projector rank4 leaves four coad8 outputs.',
        '210 changed-owner lowerings have retained actual outputs;390 unchanged cases retain identical-input/deterministic-source conditions.26builtin bis noops earn zero credit. No positional identity or cached-plan theorem.',
        'Historical2b/c78 Symbolica/Numerica and registryGraphica3.0.1; current pristine actualffa/archive8993/kernel8eb; observer built5ea+qualified overlay and proofsource5ba remain distinct. No profile/compiler/observer/callback/userdata/source qualification upgrade.',
        'Raw/cooked/decoded association, whole definitions, symbol callbacks/userdata and qualified actual source/feature/registration state remain finite captured premises. Locality912 records do not prove UV subtraction, decoder universality or arbitrary component physics.',
        'OnlyNone receives original Row assertion. Auto has its own captured worlds and separately reconstructed Rows; equal finite values do not invent an Auto literal assertion. Symfactor2 applies to4of5 sym graphs.',
        'External dependencies are retained honestly.3861 original files and all additionally consumed external files are rehashed. Two source configuration dependencies remain recorded-only and unread under the explicit no-config-read boundary; their earlier native/root verification is not a fresh check here.',
        'Cooperative native receipt consistency is retained evidence, not hostile/cryptographic attestation. Pure controls earn zero native or scientific proof credit.',
    ]
    proof=C.load(C.O/'proof-association.json');checks=C.load(C.O/'focused-controls.json')
    handoff={'schema':'next-three-raw-qcd-migration-preparation/1','status':'NOT_APPLIED_PREPARATORY_CANDIDATE_READY_FOR_ROOT_JOIN_AFTER_SCIENTIFIC_VERDICT','gates':['SCIENTIFIC_INTERFACE_REVIEW_PENDING','NEW_MIGRATION_REVIEW_REQUIRED'],'author':'Lucien Huber <im@lcnbr.ch>','workspace':str(C.W),'packet':str(C.O),'target':{'revision':C.TARGET,'tree':'c9de5c3fd962bd5d7c9375ab1f86f75699e8a69b','relevant_code':C.CODE,'path':C.MAIN,'before_blob':audit['before_blob'],'after_blob':audit['after_blob']},'patch':C.ref(C.O/'proposed-not-applied.patch'),'before':C.ref(C.O/'before.rs'),'after_not_applied':C.ref(C.O/'after.NOT-APPLIED.rs'),'candidates':{case:{k:v[k] for k in ('old','new','before_Row','after_Row')} for case,v in candidate.items()},'validation':{'ten_field_pristine_Rows_reconstructed':12,'None_literal_candidates':3,'changed_numerator_tokens':3,'all_Rows':149,'unchanged_Rows':146,'prior_migrated_full_Rows_preserved':13,'ordered_source_scopes_preserved':118,'source_objects_bound':1521,'actual_pristine_to_target_generator_delta':'ONE unrelated prior grouped numerator token; all remaining crate/asset/manifests/lock objects identical.','proof_counts':proof['counts'],'retained_math_member_associations':4268,'focused_controls':{k:checks[k] for k in ('total','accepted','rejected')},'original_files_rehashed':len(originals),'total_external_files_rehashed':len(expected),'unread_config_dependencies':2,'originals_unchanged':True,'post_patch_runtime':'NOT_RUN','Typst':'UNRENDERED; no compiler invoked','owned_processes':0},'warnings_and_conditions':warnings,'scientific_review':{'task_name':'generation-next-qcd-scientific-interface-astra-review-20261031-r1','status':'PENDING_NOT_CONSUMED','verdict':None,'root_acceptance':None},'approval':False,'applied':False,'new_native_scientific_credit':0,'input_manifest':'complete-input-manifest.json','byte_scope_audit':'byte-scope-audit.json','scope_after_audit':'all118-scopes-after-candidate.json','physical_argument':'wp9-draft.typ','errors':'errors.json','retention':'Retain entire NEW ignored packet, before/after blobs, scripts, failed authoring versions, hashes and external dependencies for root/new migration review. No owned processes; no cleanup of other work; no new durable VCS reference created under explicit no-ref-mutation scope.'}
    C.save('handoff.json',handoff)
    text='''= Next three raw QCD numerator migration preparation

Author: Lucien Huber <im@lcnbr.ch>. NOT APPLIED. This is a scoped WP9 draft,
with SCIENTIFIC_INTERFACE_REVIEW_PENDING and NEW_MIGRATION_REVIEW_REQUIRED.
No scientific verdict, root acceptance or application authorization is supplied.

== Exact candidate

The sole proposed file is crates/feynkit-generator/tests/generation_oracle/main.rs,
frozen from immutable83276dc047815d1dcc7f003aeb512c204b26758e,
tree c9de5c3fd962bd5d7c9375ab1f86f75699e8a69b, relevant code6f41807b.
Before blob7283198bf0152c82c1cd5f0173703cddfa7f5c0f; after
blobcc4505a6c4f6fc1687c075ed79e154a8e19f05a6. Exactly three numerator tokens:

- xs_qcd_gguu_sym_L1: b9fdc23ed056fa0c → e4b209849953f1fa.
- xs_qcd_gguu_L1: 5edf8e5f7a11eeb3 → 665e13a4c88da679.
- xs_qcd_gggg_L1: 2628c62dc373ac91 → 96e5abcebda25959.

All149 Rows remain;146 are byte-identical. All13 prior migrated full Rows,
all other fields/bytes and all118 ordered source-scope bindings remain exact.
The patch applies without offset or fuzz to the frozen private bytes and yields
the exact after blob. It has not been applied to any repository checkout.

== Independent reconstruction and finite argument

Twelve full ten-field Rows are reconstructed from sealed pristine JSON for
three cases, historical/current and None/Auto separately. FNV1a byte updates,
little-endian u64/u128, canonical JSON, text/name terminators, numerator/interface
field removal, vertex/edge ordering, frame hashing, exact rational bookkeeping,
multiset folds and ordered named/group folds follow the frozen real census owner.
The source report is independently joined to the pristine census; public payload,
source companion and observer equality are exact. Failed actual assertion hashes
are not inputs to the reconstruction. The supplied observed hashes corroborate
the derived values only after all fields have been calculated.

Each None historical Row equals the immutable expectation. Each current None Row
differs solely in numerators. Auto uses its own event segment, request, payload,
source and physical domain; no Auto literal assertion is invented. Original
nested-Fraction controller failure and additive revalidation remain separate.

The complete physical domains use exact128-bit IDs AND full canonical keys,
not rounded IDs, counts or ordinal correspondence. Captured model, process,
request/options/Run, signed routes, incidences, cut memberships, masses, all owner
factors and singleton groups are bound. Cases have5/9/10 diagrams and65/117/118
owners per era/mode. Historical missing interfaces remain missing; current
explicit counts are34/63/8. The symmetrized factor2 belongs to exactly4of5 graphs.

The retained candidate has600 local relations:578 full-tree literal identities
and22 nonliteral even-f-cycle relations with144 explicit rewrites, plus96
unexpanded products and912 finite locality records. Completed native evidence
has5 baseline passes and16 specific semantic rejections, zero infrastructure
outcomes/retries. Whole-definition/ordered-side native partitions are166/118
unique pairs per c78/8eb. This packet binds4268 inverse member references to
their actual result rows and exact definition order; repeated members earn no
additional proof credit. Original native outputs, helper/source bindings,
physical and locality obligations remain external, completely indexed inputs.

All600 actual typed owner imports bind exact raw/cooked-decoded/ports/flags to
their own captured case/mode/graph/owner. If pi maps historical logical axes to
current logical axes, and p_h/p_c map canonical axes to each logical interface,
the retained maps satisfy p_h = inverse(pi) composed with p_c. Actual operand
transport, canonical output maps and each surviving fold output axis are bound.
All1890 selected pair inputs and456 actual fold steps retain exact oriented
named contractions. Both operand AND output-axis maps are necessary for a
component-covariance argument; positions alone are not interchangeable.

The210 changed-owner lowerings bind actual output dangling/layout records.
The390 unchanged owners retain deterministic identical-input source conditions;
26builtin bis noops receive no rejection credit. Native typed8eb imports both
eras; no new c78 consumer ran. Rank8 numerator plus rank4 projector leaves FOUR
coad8 axes, not a scalar. Projector/full-output association here is a finite
named contraction argument, not newly executed component mathematics.

Conditional WP9 claim: if the independent scientific/interface review accepts
the exact retained finite identities, complete factors, model/routing/definition
premises and commuting operand/output maps, the three complete factorized graph
domains preserve the source-defined physical tensor expression under those maps.
The representation-sensitive numerator digests may then be considered by a NEW
independent migration review. Counts, weights, hashes and this preparer's
association checks do not themselves establish independent physical correctness.

== Applicability and retained limits

All1521 crate/asset/manifests/lock source objects are compared independently
between actual pristineffa76ec and target832. The only difference is one unrelated
prior grouped numerator token in main.rs; reverting it yields exact main bytes.
Target832 and relevant code6f have identical expectation bytes. Complete path
deltas isolate docs/evidence and notebook stub changes. Actual capture producer
is never relabelled832/6f. This proves source applicability, not a new832 build
or runtime pass. Historical2b/c78 Symbolica/Numerica with registryGraphica3.0.1,
currentffa/archive8993/kernel8eb, observer built5ea plus qualified overlay and
proofsource5ba remain distinct. Profile/compiler/feature/source qualification,
raw/cooked/decoded capture, callback/userdata and registration premises remain.

The argument is finite: three cases, separate None/Auto, GroupingNone, pool1.
It grants no general decoder/free-index/locality/UV theorem, historical positional
identity, cached-plan theorem, grouped/other-row/pool8/32 upgrade or full integrand
and component-physics approval. Cooperative receipt consistency is not hostile
or cryptographic attestation. Native report restrictions are carried verbatim
in proof-association.json. PREP and supervisor approvals confer neither math nor
migration approval. Pending scientific review is never consumed or polled.

== Validation, retention and next gate

Thirty-four meaningful pure controls accept the exact candidate/actual maps and
reject wrong case/mode/owner/source hash/literal/definitions/operand map/output
axes/endpoints, changed prior Rows and extra patch bytes. No new native/helper,
scientific import, compiler/Nix/build/network, license/config read, publication,
delegation or VCS mutation occurs. Typst is UNRENDERED. Failed authoring/schema
assumptions and their corrected source versions are retained in errors.json;
they earn no scientific negative credit.

Originals and every additionally consumed external file are rehashed before/after.
Two previously bound source-configuration files stay unread under the explicit
no-config-read boundary; their recorded hashes are retained without claiming a
fresh verification. The packet is NOT self-contained: exact capture/proof/native,
review/helper/source dependencies must remain at the indexed external locations.
No process is retained. Any concurrent ref change is disclosed without restoration.

Root must first attach the completed scientific/interface verdict and its exact
conditions, then assign NEW independent migration review of this sealed candidate.
No application or new producer qualification is authorized by this preparation.
WP3/full371/main/M-rebuild/publication remain outside scope. Retain this ignored
packet and original external evidence for review; no other work is cleaned up.
'''
    C.put('wp9-draft.typ',text.encode())
    C.put('handoff.typ',('= Unapplied next QCD migration handoff\n\nAuthor: Lucien Huber <im@lcnbr.ch>.\n\nSCIENTIFIC_INTERFACE_REVIEW_PENDING; NEW_MIGRATION_REVIEW_REQUIRED.\n\nSee handoff.json for exact hashes, before/after ten-field Rows, counts and all warnings;\nwp9-draft.typ contains the scoped conditional physical argument. The candidate changes\nonly three numerator tokens, preserving149 Rows, all13 prior full Rows and118 scopes.\nNo application, new native execution, scientific verdict or migration approval.\nAll artifacts are retained in this new ignored packet; no owned process remains.\n').encode())
    C.save('final-validation.json',{'status':'PASS_PURE_PREPARATION_ONLY','byte_scope':audit,'focused_controls':{k:checks[k] for k in ('total','accepted','rejected')},'originals_before_after_equal':True,'ordered118_after_equal':True,'external_dependencies':len(expected),'warnings_and_conditions':warnings,'refs_changed':earlier!=now,'owned_processes':0,'approval':False,'application':False})
    # All process streams below were closed by earlier completed commands.
    for p in C.O.rglob('*.json'): json.loads(p.read_bytes(),object_pairs_hook=C.unique)
    entries=[]
    for p in sorted(C.O.rglob('*')):
        if p.is_file() and p.name not in ('SEAL.json','SEAL.sha256','SHA256SUMS'):
            entries.append((str(p.relative_to(C.O)),C.fingerprint(p)))
    manifest=''.join(r['sha256']+'  '+name+'\n' for name,r in entries).encode()
    C.put('SHA256SUMS',manifest)
    seal={'schema':'next-qcd-not-applied-migration-preparation-seal/1','status':handoff['status'],'gates':handoff['gates'],'manifest_sha256':C.sha(manifest),'files':len(entries),'handoff_sha256':C.sha((C.O/'handoff.json').read_bytes()),'patch_sha256':C.sha(patch),'before_blob':C.blob(before),'after_blob':C.blob(after),'external_input_manifest_sha256':C.sha((C.O/'complete-input-manifest.json').read_bytes()),'self_contained':False,'applied':False,'approval':False,'owned_processes':0}
    C.save('SEAL.json',seal);sealhash=C.sha((C.O/'SEAL.json').read_bytes());C.put('SEAL.sha256',(sealhash+'  SEAL.json\n').encode())
    print(json.dumps({'packet':str(C.O),'seal_sha256':sealhash,'files':len(entries),'patch_sha256':C.sha(patch),'handoff_sha256':seal['handoff_sha256'],'external_inputs':len(expected),'unread_config_references':2,'refs_changed':earlier!=now,'status':handoff['status']}))

if __name__=='__main__': main()
