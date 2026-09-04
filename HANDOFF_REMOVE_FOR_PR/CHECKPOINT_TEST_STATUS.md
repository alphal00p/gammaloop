# Handoff checkpoint test status

This status belongs to the temporary handoff checkpoint prepared on
2026-09-04. It is intentionally committed under `HANDOFF_REMOVE_FOR_PR/` and
must be removed with the rest of that directory before the merge-ready PR.

## Exhaustive curated result

The primary command was:

```bash
export SYMBOLICA_LICENSE='dc3cebc3#6c76981c#e3f8341d-402d-522d-87fa-6415319915a3'
export RUST_MIN_STACK=67108864
nix develop --command just test_gammaloop --no-fail-fast
```

It compiled with warnings denied and ran 1,955 of the 1,967 tests selected by
the curated profile before nextest reached its failure cap:

```text
1,926 passed
28 failed
1 timed out
12 not run
295 excluded/skipped by the curated profile
```

The 12 withheld tests were then selected exactly and run with the same
warnings-as-errors build configuration, no retries, and no fail-fast. Eleven
passed and one failed. Combining the two non-overlapping runs gives the
exhaustive status of all 1,967 selected tests:

```text
1,937 passed
29 failed
1 timed out
0 unclassified
```

The primary log was `/tmp/raised-energy-handoff-test.log`; the residual log was
`/tmp/raised-energy-handoff-residual.log`. Those machine-local logs are not
committed; all non-passing tests are recorded below.

## Failing tests

1. `gammalooprs cff::tests::stored_production_cff_uses_typed_energy_factor_ownership`
2. `gammalooprs cff::tests::exact_unexpanded_cff_matches_ordinary_for_cross_loop_factorized_energy`
3. `gammalooprs graph::three_d_source::tests::exact_sub_lmb_taylor_vacuum_keeps_dotted_bubble_owner_incidence`
4. `gammalooprs cff::tests::exact_cff_lu_residue_factorizes_from_quadratic_cubic_spectator`
5. `gammalooprs graph::lmb::test::shrunken_disconnected_subgraph`
6. `gammalooprs uv::approx::local_3d::tests::projected_source_sum_falls_back_when_the_cut_excludes_compatible_hosts`
7. `gammalooprs uv::approx::tests::expanded_4d_setting_does_not_change_the_empty_forest_root`
8. `gammalooprs uv::hedge_poset::tests::union_terms_replay_component_paths_from_typed_roots`
9. `gammalooprs uv::approx::projected_4d::tests::typed_taylor_next_component_reuses_one_topology_for_prior_residue_states`
10. `gammaloop-integration-tests::test_runs smoke::addbar`
11. `gammalooprs uv::tests::scalars_profile_new`
12. `gammalooprs uv::tests::scalars_profile`
13. `gammaloop-integration-tests::test_runs test_3d_reps::cff_cli_validate_and_build_use_gammaloop_graph_state`
14. `gammalooprs::test_renormalization finite_part_ghost_2loop`
15. `gammaloop-integration-tests::test_runs test_integrated_uv_cts::important::scalar_bubble_inspect`
16. `gammaloop-integration-tests::test_runs repeated_masses::repeated_bubble_split_mass_limit_without_threshold_subtraction`
17. `gammaloop-integration-tests::test_runs repeated_masses::repeated_bubble_model_mass_limit_without_threshold_subtraction`
18. `gammaloop-integration-tests::test_runs inspect::bare_raised_scalar_self_energy_has_the_reported_uv_degree`
19. `gammaloop-integration-tests::uv scalar_spectacles_integrated_uv_factorizes_over_bridge`
20. `gammaloop-integration-tests::test_runs scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::quadratic_energy_numerators::scalar_3l_cross_section_gl16_quadratic_energy_inspects_match::q7_squared`
21. `gammaloop-integration-tests::test_runs scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl04_right_temporal_product_inspects_match`
22. `gammaloop-integration-tests::test_runs scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl04_q5_temporal_square_inspects_match`
23. `gammaloop-integration-tests::test_runs smoke::oak`
24. `gammaloop-integration-tests::test_runs scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl04_right_owned_dot_inspects_match`
25. `gammaloop-integration-tests::test_runs scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::quadratic_energy_numerators::scalar_3l_cross_section_gl00_quadratic_energy_inspects_match::q7_squared`
26. `gammaloop-integration-tests::test_runs scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl04_inspects_match`
27. `gammaloop-integration-tests::uv sunrise_scalar_1_uv`
28. `gammalooprs uv::approx::tests::gl24_direct_3d_modes_preserve_orientation_selector_contracts`
29. `gammaloop-integration-tests::test_runs integrations::v_diag`

The additional residual failure, `integrations::v_diag`, fails during generation
at `crates/gammalooprs/src/uv/approx/direct_3d/forest.rs:155` with:

```text
orientation pattern selects no production residue maps
```

## Timed-out test

- `gammaloop-integration-tests::test_evaluation_api gl20_multichannel_local_inspect_event_snapshot`

It timed out on both attempts at the configured 600-second boundary while
generating the selected GL20 NLO graph.

## The exact residual set

The following 12 tests were withheld by the primary run. Their checkpoint
status is recorded explicitly so none is mistaken for untested coverage:

| Test | Status |
|---|---|
| `gammalooprs integrate::tests::integration_result_always_contains_first_non_trivial_discrete_breakdown` | pass |
| `gammalooprs integrate::tests::integration_workspace_manifest_accepts_only_the_current_version` | pass |
| `gammaloop-integration-tests::test_runs differential::lu_differential_integration_cli_flag_writes_iteration_observables` | pass |
| `gammaloop-integration-tests::test_runs differential::lu_differential_integration_hwu_output_is_optional_and_single_file` | pass |
| `gammaloop-integration-tests::test_runs integrations::scalar_bubble` | pass |
| `gammaloop-integration-tests::test_runs integrations::v_diag` | **fail** |
| `gammaloop-integration-tests::test_runs multi_integrand::test_integration_workspace_model_mismatch_requires_restart` | pass |
| `gammaloop-integration-tests::test_runs multi_integrand::test_integration_workspace_resume_preserves_current_effective_model_override` | pass |
| `gammaloop-integration-tests::test_epem_a_ddx_nlo_acceptance epem_a_ddx_nlo_is_alpha_s_over_pi_times_lo_in_all_local_uv_routes` | pass |
| `gammaloop-integration-tests::test_epem_a_ttx_nlo_acceptance epem_a_ttx_msbar_nlo_matches_the_published_inclusive_ratio_in_all_local_uv_routes` | pass |
| `gammaloop-integration-tests::test_gamma_star_ddx_nlo_acceptance gamma_star_ddx_msbar_nlo_matches_the_published_graph_cross_sections` | pass |
| `gammaloop-integration-tests::test_gamma_star_ttx_nlo_acceptance gamma_star_ttx_msbar_nlo_matches_the_published_absolute_cross_sections` | pass |
