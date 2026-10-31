"""Produce command data; this program runs no producer, test, or discovery."""
from pathlib import Path
import json, shlex
P=Path(__file__).parent
CODE='2cfb52b28d99f92305a9d571bbe3911cf04b8a19'
URI='git+file:///shared/lcnbr/dev/gammaloop/main?rev='+CODE+'&shallow=1&allRefs=1'
SOURCE='/nix/store/3sj3nw06s3vhl3j12wfp5gc3r1h6f0k3-source'
KERNEL='8eb7c3729d6abefc916decc1bdf0c06536f8af80'
roots=json.loads((P/'producer-roots.json').read_text())
front=json.loads((P/'demand-frontier.json').read_text())
old=json.loads((P/'historical-command-index.json').read_text())
license_vars=['SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE']
unset=['env']+[x for k in license_vars for x in ['-u',k]]
nix=['nix','--offline','--option','allow-import-from-derivation','false','--option','accept-flake-config','false']
nextest='/nix/store/kwm3g9vh36d9zd171v0nz64fyan1hvfz-cargo-nextest-0.9.143/bin/cargo-nextest'
generator='/nix/store/qg5nc1vcxjbdb5zzjigilx3wngy8xyfd-gammaloop-nextest-binaries-feynkit-feynkit-generator-0.1.0/gammaloop-nextest-feynkit-feynkit-generator.tar.zst'
def archive(key):
 group,package=key.split('/')
 return roots[key]['output']+'/gammaloop-nextest-'+group+'-'+package+'.tar.zst'
def frontier_for(keys):
 wanted={roots[k]['drv']:{'out'} for k in keys};todo=list(wanted);seen=set();missing=[];reused=[]
 while todo:
  drv=todo.pop(); requested=tuple(sorted(wanted[drv])); marker=(drv,requested)
  if marker in seen:continue
  seen.add(marker);d=front['graph'][drv];absent=False
  for out in requested:
   path=d['outputs'][out]
   # Match the observation recorded in the sealed frontier, not live store changes.
   if any(r['drv']==drv and r['output']==out for r in front['missing_required_outputs']):missing.append((drv,out));absent=True
   else:reused.append((drv,out))
  if absent:
   for dep,outs in d['inputDrvs'].items():
    prior=wanted.setdefault(dep,set())
    if set(outs)-prior:prior.update(outs);todo.append(dep)
 return {'roots':len(keys),'missing_required_outputs':len(set(missing)),'reused_frontier_outputs':len(set(reused)),'examined_derivations':len(wanted)}
core_keys=[k for k in roots if k.startswith('core/')]+['integration/gammaloop-integration-tests']
all_keys=list(roots)
def build(keys):
 return unset+nix+['build','--no-link','--print-out-paths','--print-build-logs','--max-jobs','1','--cores','2']+[roots[k]['drv']+'^out' for k in keys]

cards=[]
def runtime(id,key,filterset=None,count=None,when='current8eb production consumer qualification; not an extra universal bounded-review gate',ignored=False,bypass=False,historical_count=None):
 path=generator if key=='qualified-generator' else archive(key)
 common=['--archive-file',path,'--extract-to','/work/'+id,'--extract-overwrite','--workspace-remap','/build/source','--profile','ci_gammaloop']
 if bypass:common+=['--ignore-default-filter']
 if filterset:common+=['-E',filterset]
 if ignored:common+=['--run-ignored','all']
 cards.append({'id':id,'phase':'FUTURE_NATIVE_DISCOVERY_AND_RUNTIME_ONLY; NOT EXECUTED','source_code':CODE,'source_mount':SOURCE+' -> /build/source read-only','producer':key,'archive':path,'when':when,
  'list_argv':[nextest,'nextest','list',*common,'--message-format','json'],
  'run_argv':[nextest,'nextest','run',*common,'--test-threads','16','--retries','0','--no-fail-fast','--final-status-level','all','--no-tests=fail'],
  'expected_selected_exact':count,'historical_selected_not_a_current_result':historical_count,
  'count_contract':'Before run, retain actual binary+test IDs and runnable/ignored/filtered inventory; require this exact count when specified, otherwise exact frozen discovered IDs and a nonempty selection. After run require same IDs/count, zero failures/errors/skips/retries and successful command plus JUnit. Do not infer a count from source or an old receipt.',
  'deduplication':'Skip this extra execution if an already accepted current-source/input-invariant run, or the selected canonical core run, includes every exact test and retains its result; no double credit.'})
focus=['integrands::process::tests::lmb_sampling_detects_only_oppositely_stored_parallel_edges','integrands::process::tests::graph_derived_two_loop_reflections_preserve_local_odd_numerator_and_partition','processes::cross_section::gauge_sewing_tests::massive_tau_covariant_cuts_satisfy_the_brst_quartet_identity']
runtime('core-inference-brst','core/gammalooprs',' or '.join('test(='+x+')' for x in focus),3,historical_count=3)
runtime('api-topology','core/gammaloop-api','test(=commands::generate::tests::topology_card_projector_crosses_the_current_cli_boundary)',1,historical_count=1)
runtime('evaluator-api-command','integration/gammaloop-integration-tests','binary(=test_evaluation_api) or binary(=test_evaluation_command)',17,
        when='Proposed focused source-impact selection of the existing16 API and1 command tests; not a separately enumerated accepted-plan gate. Confirm17 by actual discovery before execution.')
runtime('ghost-statistics','integration/gammaloop-integration-tests','binary(=test_feyngen) and test(=closed_ghost_loops_preserve_statistics_weights)',1)
runtime('nlo-coolqm','integration/gammaloop-integration-tests',' or '.join('test(='+x+')' for x in ['epem_a_ddx::epem_a_ddx_nlo_is_alpha_s_over_pi_times_lo_in_all_local_uv_routes','epem_a_ttx::epem_a_ttx_msbar_nlo_matches_the_published_inclusive_ratio_in_all_local_uv_routes','cool_qm_eos::cool_qm_eos_nlo_matches_example_target']),3,
        when='Original R/F WP3 acceptance selection; not a new requirement to close every bounded integration step.',historical_count=3)
runtime('sm3150','integration/gammaloop-integration-tests','test(=test_vacuum_amplitude_generation_full_sm_three_loop_filtered)',1,
        when='Original WP3/R and final gate; body must actually check3150 and pools1/8 with unchanged weights/zero controls.',bypass=True,historical_count=1)
runtime('a-ddx-one-process','integration/gammaloop-integration-tests','test(=test_generate_sm_a_ddx)',1,
        when='Original WP3/R and final gate: invoke in FIVE fresh processes with five distinct work/JUnit directories; one selection/pass each. No five-process claim from one run.',bypass=True,historical_count=1)
b2=['grouping::tests::vacuum_group28_original_pair_identity','grouping::tests::sampled_quotient_groups_physical_local_factors','grouping::tests::sampled_quotient_requires_scalar_interfaces_and_all_samples','grouping::tests::sampled_quotient_rejects_undefined_ratios','grouping::tests::sampled_quotient_does_not_group_nonfinite_operands','grouping::tests::sampled_quotient_checks_nonfinite_operands_in_every_sample']
runtime('b2-six','qualified-generator',' or '.join('test(='+x+')' for x in b2),6,
        when='Current8eb R/F B2 qualifier still distinct from accepted pool1 package run because original-pair test is ignored. Existing qualified archive suffices; no new generator build.',ignored=True,historical_count=6)
runtime('source46-colour','spenso/idenso','test(color::)',historical_count=192,
        when='Optional current8eb requalification of the approved source46 implementation. Do not reimplement it; actual profile count must be discovered and reconciled, not assumed192.')
runtime('source46-two-controls','spenso/idenso','test(=color::test::failing::test_color_matrix_element) or test(=color::test::failing::test_color_matrix_element_two)',2,bypass=True,historical_count=2)
tensor=['reduction::blocks::tests::loop_vector_sums_are_retained_under_a_rank_two_projector','reduction::blocks::tests::mixed_loop_and_fixed_sums_keep_the_fixed_rank_two_term','reduction::external::tests::one_external_direction_rank_one_and_mixed_rank_two','reduction::tests::graph_tensor_reduce_selects_internal_edges_and_promotes_lorentz_slots','reduction::tests::reduction_keeps_declared_spectator_axes_and_zero_interfaces','reduction::tests::reduction_materializes_positional_spectators_with_actual_witnesses']
runtime('tensor-six-controls','feynkit/feynkit-tensor',' or '.join('test(='+x+')' for x in tensor),6,historical_count=6)
for key in ['feynkit/feynkit-graph','feynkit/feynkit-model','feynkit/feynkit-tensor','spenso/spenso','spenso/spenso-hep-lib']:
 filt='test(finalization::tests::) or test(numerator_interface::tests::)' if key.endswith('feynkit-graph') else 'binary(='+key.split('/')[1].replace('-','_')+')'
 runtime('library-'+key.split('/')[1],key,filt,when='Recovered graph/interface or library consumer selection; broader current8eb qualification menu. Historical tensor202+17doctests is a different dev-optim/native contract, not an expected current archive count.')
runtime('docs-builder-lib','docs/alphal00p-docs-builder','binary(=alphal00p_docs_builder)',historical_count=116,
        when='Recovered docs-builder library selection; count116 historical, actual current inventory required. Compilation inputs alone do not qualify the final metadata runtime sources.')
# Carry exact accepted original prerequisites and scoped sampling selectors, with package mapping explicit.
original=json.loads((P/'original-selection-records.json').read_text())
selectors=[]
for r in original:
 if isinstance(r['value'],list) and '-E' in r['value']:
  argv=r['value'];selectors.append({'source':r['source'],'json_pointer':r['json_pointer'],'filter':argv[argv.index('-E')+1],'original_argv':argv,'profile_and_features':'Preserve original in reference record; use matching normal pinned archive context only after actual feature/profile qualification. No claim old dev-optim and ci-optim ELFs identical.'})
canonical_core=[]
for key in roots:
 if not key.startswith('core/'):continue
 canonical_core.append({'producer':key,'archive':archive(key),'list_argv':[nextest,'nextest','list','--archive-file',archive(key),'--extract-to','/work/canonical-'+key.split('/')[1],'--workspace-remap','/build/source','--profile','ci_gammaloop','--message-format','json'],
 'run_argv':[nextest,'nextest','run','--archive-file',archive(key),'--extract-to','/work/canonical-'+key.split('/')[1],'--extract-overwrite','--workspace-remap','/build/source','--profile','ci_gammaloop','--test-threads','16','--retries','0','--no-fail-fast','--no-tests=pass'],
 'count_contract':'Historical total2249 inventory/2214 selected is a comparison, not a current count. Record exact current inventory and selected IDs per package. Zero-test macros package may be recorded as zero with no pass credit; other packages and group total must be nonempty. No default-filter/ignored credit.'})
canonical_static=[
 {'id':'format','phase':'future pure static check','argv':unset+nix+['build','--no-link','--print-out-paths',URI+'#checks.x86_64-linux.gammaloop-fmt']},
 {'id':'ci-graph','phase':'future pure static check','argv':unset+nix+['build','--no-link','--print-out-paths',URI+'#checks.x86_64-linux.gammaloop-guppy-workspace-graph']},
 {'id':'clippy','phase':'future pure compilation/static check; no runtime credit','argv':unset+nix+['build','--no-link','--print-out-paths','--max-jobs','1','--cores','2',URI+'#checks.x86_64-linux.gammaloop-clippy']},
 {'id':'hakari-diff','phase':'future pure source/static check in exact disposable immutable source copy; no writes to selected checkout','argv':unset+['nix','--offline','develop',URI,'-c','cargo','hakari','generate','--diff']},
 {'id':'hakari-verify','phase':'future pure source/static check in exact source copy','argv':unset+['nix','--offline','develop',URI,'-c','cargo','hakari','verify']},
 {'id':'feature-isolation','phase':'future pure source/static check in exact source copy','argv':unset+['nix','--offline','develop',URI,'-c','just','check-symbolica-feature-isolation']},
]
docs={
 'immediate_bounded_scope':'Current structural/TOML/title/118 scope checks and faithful evidence attribution; no full documentation-runtime or site approval inferred.',
 'original_commands_pinned_to_code2cf':[
  ['nix','--offline','develop',URI,'-c','cargo','test','--locked','-p','alphal00p-docs-builder','--lib','--','--test-threads=1'],
  ['nix','--offline','develop',URI,'-c','cargo','run','--locked','-p','alphal00p-docs-builder','--','check'],
  ['nix','--offline','develop',URI,'-c','just','docs-check']],
 'classification':'Original workflow command references only, NOT commands to run after license handoff. They can build and/or invoke other native exporters. For the strict two-phase boundary, prepare normal pinned producers and all compiler/catalog/exporter artifacts with four license vars omitted; then execute their exact resolved runtime commands in a fresh no-network capsule. Full final docs belongs WP13 and final metadata must be frozen first.',
 'pure_artifact_owner_argv':unset+nix+['build','--no-link','--print-out-paths','--max-jobs','1','--cores','2',URI+'#packages.x86_64-linux.alphal00p-docs-cargo-artifacts'],
 'canonical_check_owner_eval_argv':unset+nix+['eval','--raw',URI+'#checks.x86_64-linux.alphal00p-docs.drvPath'],
 'remaining_exact_runtime_gap':'Builder/exporter/catalog executable member paths, runtime source closure and final metadata inputs must be resolved from this normal owner before strict licensed execution. No executable path or full-site result invented here.'}

obj={'schema':'generation-combined-qualification-command-data/1','NOT_EXECUTED':True,'immutable_code':CODE,'immutable_evidence':'457e361d7d39e0a3daae801651f8c1dfabbf47aa','kernel':KERNEL,'source':SOURCE,'source_nar':'sha256-RCGp7HErhoxa4nMbFUDXneWTb4992VxGY1T7eB57Rmk=',
 'decision_boundary':{'mandatory_next':'Freeze truthful additive metadata and obtain a NEW independent combined review with exact scoped evidence and unresolved production-core/WP3 exclusions. Do not infer approval or advance global bookmark.','no_extra_generator_or_host_build_for_frozen457':'Existing qualified ffa generator producer is EXACT same2cf selected derivation; approved ordinary/private host compositions already reviewed. New metadata alone does not require their reruns.','production_qualification_missing':'The selected current8eb core/API/evaluator archives and library/docs test archives in producer-roots are absent; needed for new current production-consumer or WP3 claims, not automatically prerequisites for every bounded review.','checklist_ambiguity':'Accepted plan specifies bounded review and R/F matrix, not an exhaustive standalone refresh-core/API/evaluator/library list before any bookmark advance. Proposed impact selections below are labeled; reviewer/root may require them for the scope they actually approve. No permission question is reopened.'},
 'pure_before_licensed':{'pure':['normal immutable Nix evaluation','local store/NAR/derivation/source membership queries','authorized pure realization/static checks after this task'], 'omitted_environment':license_vars,'IFD':False,'substitution':'No network queried here. Offline build argv use cached local dependencies; no remote availability claim. Parent may coordinate normal substitution only under its separate authorized execution scope.','artifact_qualification':'After a build, bind actual compiler argv, synthetic manifests, source/dependency NARs, kernel triple, complete features/profile, retained layered fingerprints and archive/member/loader identities before runtime. A successful producer is compilation only.','runtime':'Fresh minimal bwrap --unshare-all capsule, source/store inputs read-only, writable /work and /tmp, network namespace distinct, no host home/config/Git/Nix socket. Source at/build/source; TESTS_GAMMALOOP_STATE_PATH=/work/test-states. Obtain loader closure and runtime tools from pinned Nix owners before handoff. Some evaluator tests invoke FORM/GCC/runtime compilation: keep those exact required runtime tools available inside the capsule and attribute generated compiler work; no build/Nix command with the private license in its environment.','license':'No license/config read here. Parent preserves an exported license and supplies it privately only to the immediate no-network runtime child through existing anonymous-pipe handoff. Never include values/hashes/files in commands, store or logs.','runtime_environment':{'TESTS_GAMMALOOP_STATE_PATH':'/work/test-states','TMPDIR':'/tmp','NEXTEST_TEST_THREADS':'16','NEXTEST_RETRIES':'0','INSTA_UPDATE':'no','FEYNKIT_ORACLE_THREADS':'explicit by selected gate; no guessed8/32 credit'},'runtime_unset':['INSTA_FORCE_PASS','INSTA_UPDATE_FAST','FEYNKIT_ORACLE_RECORD','FEYNKIT_ORACLE_DUMP'],'limits':'Canonical ci_gammaloop240s,zero retries,16workers,two memory-intensive slots unchanged. No20GiB cap. Discovery counts and successful JUnit completion are mandatory.'},
 'coordinated_build_core_and_evaluator':{'select_only_if_claim_needed':True,'producer_keys':core_keys,'frontier':frontier_for(core_keys),'argv':build(core_keys)},
 'coordinated_build_all_recovered_consumers':{'optional_union_not_a_new_gate':True,'producer_keys':all_keys,'frontier':frontier_for(all_keys),'argv':build(all_keys),'prefer_single_union':'If root selects these consumers, realize the union once; do not also launch the smaller overlapping build. Existing exact dependency outputs are reused.'},
 'runtime_cards':cards,'canonical_core_fast':canonical_core,'original_scoped_selectors':selectors,'static_cards':canonical_static,'documentation':docs,
 'later_grouped_delta':'later-grouped-applicability.json;13 non-generator compile filters and required source graphs omit main.rs. Sole token delta needs generator rebuild/39passes/joins by assigned writer. Re-attest exact non-generator drv identities on writer frozen commit before transfer; docs runtime/registry inputs separately change. No postapplication pass observed.'}
(P/'commands.json').write_text(json.dumps(obj,indent=2)+'\n')
print(json.dumps({'runtime_cards':len(cards),'core_frontier':obj['coordinated_build_core_and_evaluator']['frontier'],'union_frontier':obj['coordinated_build_all_recovered_consumers']['frontier']}))
