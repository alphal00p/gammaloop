# Spenso Network Execution Test Baseline

Baseline captured on 2026-05-05 in jj change `lwvmvnsn`
(`profile raw full-depth Spenso parsing and isolate execution orchestration`).

Command:

```bash
just test-ci
```

Result:

```text
Nextest run ID: ce9f0e8a-4eb3-45f3-aafc-dc3251c4ae05
Summary: 1107 tests run: 993 passed, 112 failed, 2 timed out, 128 skipped
Exit code: 100
JUnit artifact: target/nextest/ci_gammaloop/target/nextest/ci/junit.xml
```

The two timeouts are the large Spenso execution diagnostics added for the raw
large input. Both parsed into symbolic networks in seconds before timing out in
execution:

- `idenso::large_spenso_inputs::smallest_root_input_symbolic_network_parse_default`
  parsed in `1.554s`, graph nodes `69975`, graph edges `123483`, tensors
  `30184`, scalars `11663`.
- `idenso::large_spenso_inputs::smallest_root_input_symbolic_network_parse_without_scalar_precontraction`
  parsed in `2.774s`, graph nodes `200789`, graph edges `254297`, tensors
  `30184`, scalars `89427`.

The exact non-passing records from the JUnit artifact are:

```text
gammaloop-api::commands::set::test::set_process_updates_are_transactional_across_multiple_targets
gammaloop-api::commands::set::test::set_process_without_integrand_applies_to_all_generated_integrands_in_one_process
gammaloop-api::commands::set::test::set_process_rejects_integrand_without_process_target
gammaloop-api::commands::set::test::set_process_without_targets_applies_to_all_generated_integrands
gammaloop-api::state::tests::activate_loaded_integrand_backends_falls_back_to_eager_when_external_artifacts_are_missing
gammaloop-integration-tests::test_cli::cli_stateful_workflow_behaviors
gammalooprs::uv::hedge_poset::tests::sunrise
gammalooprs::uv::hedge_poset::tests::fourloop_b
gammalooprs::uv::hedge_poset::tests::mercedes
gammalooprs::uv::hedge_poset::tests::triple_tadpole
gammalooprs::uv::hedge_poset::tests::dumbells
gammalooprs::processes::amplitude::test::generation_orientation_pattern_filters_evaluator_orientations
gammalooprs::uv::hedge_poset::tests::dotted_sunrise
gammalooprs::uv::hedge_poset::tests::basketball
gammalooprs::uv::tests::scalars_profile_new
gammalooprs::uv::hedge_poset::tests::bugblatter
gammalooprs::uv::hedge_poset::tests::four_loop_a
gammalooprs::uv::hedge_poset::tests::dotted
gammalooprs::uv::hedge_poset::tests::spectacles
gammalooprs::uv::tests::spinney_partial_cmp_is_equal_for_identical_subgraphs
gammalooprs::uv::tests::scalars_profile
idenso::metric::test::permute
idenso::color::test::test_color_structures
idenso::tensor::contract::test::parse
idenso::color::test::minus_sign
idenso::tensor::test::parse
idenso::chain::tests::collect_oppositely_oriented_chain
idenso::dirac::test::feyncalc_reference::dirac_simplify_id6_projectors_move_through_gamma
idenso::color::test::form_reference::su_f4f4
idenso::color::test::form_reference::chain_two_generator_trace_normalizes
idenso::color::test::form_reference::tloop_g14
idenso::color::test::form_reference::chain_one_generator_trace_vanishes
idenso::color::test::form_reference::tloop_qgloop_size_3
idenso::color::test::form_reference::adjacent_generator_casimir_chain
idenso::tensor::tests::parsing::parse_scalar_tensors_step_by
idenso::color::test::form_reference::separated_generator_casimir_trace
idenso::color::test::form_reference::su_fnfn_n3
idenso::color::test::form_reference::tloop_fiveq
idenso::color::test::form_reference::su_a4a4
idenso::color::test::form_reference::three_generator_trace_terminal
idenso::dirac::test::form_reference::gamma_five_symmetric_12_term_count
idenso::dirac::test::gamma_chain_canonical_ordering_is_opt_in
idenso::color::test::form_reference::su_f4a4
idenso::tensor::tests::parsing::dot_derivative
idenso::color::test::form_reference::four_generator_trace_terminal
idenso::tensor::tests::parsing::parse_val
idenso::color::test::form_reference::simpli_contracts_projected_f_pair
idenso::color::test::form_reference::symmetric_invariant_d33_partial_contraction
idenso::tensor::tests::parsing::infinite_execution
idenso::color::test::form_reference::two_f_loop_contracts_to_ca_metric
idenso::color::test::form_reference::su_f3f3
idenso::color::test::form_reference::tloop_qloop_size_3
idenso::metric::test::dots
idenso::color::test::form_reference::mixed_trace_structure_contraction
idenso::color::test::form_reference::separated_generator_casimir_shortcut
idenso::dirac::test::form_reference::gamma_five_regular_12_term_count
idenso::color::test::form_reference::three_f_loop_contracts_to_ca_f
gammaloop-integration-tests::test_runs::examples::test_scalar_bubble_example_cli
gammaloop-integration-tests::test_runs::inspect::test_mass_approach_threshold_subtraction
gammaloop-integration-tests::test_runs::smoke::oak
gammaloop-integration-tests::test_runs::smoke::addbar
gammaloop-integration-tests::test_runs::inspect::test_mass_approach_scalar_self_energy
gammaloop-integration-tests::test_runs::inspect::inspect_x_space_reports_invalid_coordinate_count_cleanly
gammaloop-integration-tests::test_runs::test_integrated_uv_cts::important::scalar_bubble_inspect
gammaloop-integration-tests::test_runs::test_integrated_uv_cts::scalar_bubble_integrated
gammaloop-integration-tests::test_runs::multi_integrand::test_inspect_uses_per_integrand_model_parameters
gammaloop-integration-tests::test_runs::multi_integrand::test_multi_integrand_batching_preserves_results
gammaloop-integration-tests::test_runs::profile_bulk::massless_triangle_bulk_profile_passes
gammaloop-integration-tests::test_runs::inspect::inspect_momentum_space_graph_id_is_sampling_agnostic
gammaloop-integration-tests::test_runs::events::amplitude_events_surface_threshold_counterterms_and_reproduce_weight
gammaloop-integration-tests::test_runs::events::amplitude_observables_accumulate_without_returning_events
gammaloop-integration-tests::test_runs::multi_integrand::test_multi_integrand_with_local_model_parameters
gammaloop-integration-tests::test_runs::inspect::inspect_x_space_reports_missing_discrete_dimensions_cleanly
gammaloop-integration-tests::test_runs::events::amplitude_selectors_generate_internal_events_without_surfacing_them
gammaloop-integration-tests::test_runs::aa_aa::important::aa_aa_local_inspect_precisions_and_backends
gammaloop-integration-tests::test_runs::multi_integrand::test_multi_integrand
gammaloop-integration-tests::test_runs::spin_sums::cross_section_vector_spin_sum_matches_explicit_incoming_helicity_average
gammaloop-integration-tests::test_runs::differential::lu_save_dot_silently_overwrites_existing_files
gammaloop-integration-tests::test_runs::spin_sums::important::cross_section_fermion_spin_sum_matches_explicit_incoming_helicity_average
gammaloop-integration-tests::test_runs::differential::lu_differential_integration_hwu_output_is_optional_and_single_file
gammaloop-integration-tests::test_runs::differential::lu_differential_integration_cli_flag_writes_iteration_observables
gammaloop-integration-tests::test_runs::differential::lu_differential_observables_without_selectors_still_fill_histograms
gammaloop-integration-tests::test_runs::integrations::v_diag
gammaloop-integration-tests::test_runs::integrations::scalar_bubble
gammaloop-integration-tests::test_runs::multi_integrand::test_integration_workspace_model_mismatch_requires_restart
gammaloop-integration-tests::test_runs::multi_integrand::test_integration_workspace_resume_preserves_current_effective_model_override
idenso::network_informed::network_informed
gammalooprs::test_renormalization::scalar_pole_part
gammalooprs::test_renormalization::finit_part_ghlo
gammalooprs::test_renormalization::finite_part_quark_lo
gammalooprs::test_renormalization::finite_part_ghost_2loop
spenso::structure::slot::shadowing_tests::to_symbolic
spenso-hep-lib::gamma_algebra_validate::validate
idenso::large_spenso_inputs::smallest_root_input_symbolic_network_parse_default
idenso::large_spenso_inputs::smallest_root_input_symbolic_network_parse_without_scalar_precontraction
gammaloop-integration-tests::test_feyngen::cross_section_standalone_export_writes_archive_and_loader
gammaloop-integration-tests::test_feyngen::cp_fix_from_symbolica
gammaloop-integration-tests::test_feyngen::test_vacuum_amplitude_kaapo
gammaloop-integration-tests::uv::dod1_bubble_uv
gammaloop-integration-tests::uv::dod2_bubble_uv
gammaloop-integration-tests::uv::dod0_bubble_uv
gammaloop-integration-tests::uv::epem_a_bbx_amp_uv
gammaloop-integration-tests::test_evaluation_api::lu_rust_xspace_integrator_weights_align_per_sample_and_scale_events
gammaloop-integration-tests::test_evaluation_api::lu_rust_evaluate_samples_respect_event_generation_and_observables
gammaloop-integration-tests::test_evaluation_api::lu_rust_precise_evaluate_samples_returns_precision_tagged_results
gammaloop-integration-tests::test_evaluation_api::lu_rust_momentum_space_evaluate_sample_reports_no_parameterization_jacobian
gammaloop-integration-tests::test_evaluation_api::lu_rust_momentum_space_integrator_weights_scale_events_without_jacobian
gammaloop-integration-tests::test_evaluation_api::lu_rust_default_clustered_pdgs_match_explicit_massless_qcd_list
gammaloop-integration-tests::test_evaluation_api::lu_rust_generated_events_follow_graph_grouping_and_cut_ids
gammaloop-integration-tests::test_evaluation_api::lu_rust_get_integrand_info_reports_groups_orientations_lmbs_and_cuts
gammaloop-integration-tests::test_evaluation_api::lu_rust_xspace_integrator_weights_scale_observable_histograms
gammaloop-integration-tests::test_evaluation_api::lu_rust_explicit_lmb_multichanneling_groups_channel_events_and_tags_metadata
gammaloop-integration-tests::test_evaluation_api::lu_rust_evaluate_samples_can_override_returned_events_without_forcing_internal_generation
gammaloop-integration-tests::test_evaluation_api::lu_rust_minimal_output_keeps_events_and_observables
```
