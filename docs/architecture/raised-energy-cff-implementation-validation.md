> Historical record for the revisions and dates stated below. Its test counts,
> API descriptions and continuation instructions do not certify or govern the
> reconstructed stack. See the [stack review](raised-energy-cff-stack-review.md)
> and [fresh validation record](raised-energy-cff-stack-review-validation.md),
> together with [current architecture](architecture-current.md) and
> [CONTRIBUTING.md](../../CONTRIBUTING.md). The original body is preserved.

# Raised-energy CFF implementation validation

Production validation ran on 2026-09-08; approved test-oracle corrections were verified on
2026-09-09 for implementation change `wuuyylluomunwsqoonmvmqttszsnynvr`.
The reviewed tip was `91142139e1fde44cc4709a35ad090631d27ca451`; the untouched
pre-implementation baseline used for the independent pole check was
`55db15c1c14ae34d1bdae8566c8e1c3c50a494e2`.

Commands ran from the repository root in its Nix development environment with the supplied
Symbolica license available to the processes and `RUST_MIN_STACK=67108864`. The license is
not included in this record. Independent tests may run concurrently with this license. The
final main suites used Nextest profile `test_gammaloop`, four test threads and no retries.
Existing tests retained their own execution-order annotations; data serialization and its
behavioral tests remain in scope. Timings below exclude compilation and build-lock waits.

## Coefficient migration: 2026-09-09

The follow-up change `wuqmzzux` migrates the four `LinearEnergyExpr` coefficient fields and
`CFFVariant::prefactor` to `Rational`, including constructors, consumers and native serde/bincode
support. All **805 distinct selected tests are verified passing** across the migrated
production-code runs and the approved fixture's focused rerun. The original selection passed
804 tests; the one failed fixture is now corrected and passes. Repeated executions are counted
once, and no automatic retries were used. The two unaffected Spenso tests from the earlier
807-test selection were not rerun. This is an affected-scope validation, not a full-workspace pass.

| Suite | Executed | Passed | Failed | Seconds | Session log |
| --- | --- | --- | --- | --- | --- |
| shared | 106 | 106 | 0 | 10.579 | `/tmp/raised-rational-shared-tests.log` |
| core before fixture correction | 276 | 275 | 1, resolved below | 115.068 | `/tmp/raised-rational-core-tests.log` |
| corrected affine fixture | 1 | 1 | 0 | 0.634 | `/tmp/raised-rational-approved-affine-tests.log` |
| api | 375 | 375 | 0 | 13.231 | `/tmp/raised-rational-api-tests.log` |
| integrations | 48 | 48 | 0 | 801.824 | `/tmp/raised-rational-integrations-tests.log` |

Formatting and all-target checks/Clippy passed for the affected shared/core/API/integration
packages. The shared crate passed all-target checks with both all features and default features
disabled. The latter retains pre-existing unused imports/helpers in eval-gated test modules;
`proc-macro-error2` also retains its existing future-compatibility notice. No new lint warnings
remain in the all-feature shared or consumer checks.

The existing persistence test now exercises negative fractions and coefficients exceeding i64,
then compares complete public expressions and energy maps after binary and JSON round-trips.
Source-mapping fixtures express their arbitrary symbolic witnesses through existing energy
substitutions; every original assertion, GL04 common-chart certificate and tolerance is
unchanged. The corrected fixture substitution is recorded below. Native Rational serialization
changes coefficient encodings; old-format compatibility is not retained under repository policy.

The first shared check found two Rust ownership/type errors, which were corrected before tests
ran. A now-unused production import was moved into its test module before the final checks.
The initial run failed `exact_taylor_vacuum_affine_circulation_preserves_energy_mapping`:
its second, independently constructed pinched mapper lacked the newly introduced symbolic
witness substitution. The expressions differed only by `OSE(7)` versus `affine_pole`.
With explicit user approval, the fixture now shares that witness substitution between both
mappers, retaining their own physical replacements and every original assertion. The focused
rerun passed. Formatting, the core all-target check and Clippy also passed after the correction.
Only this test fixture changed; the passing production-code results remain applicable.

The driver ran in the existing licensed Nix development shell with `RUST_MIN_STACK=67108864`
and `INSTA_UPDATE=no`. Independent Nextest suites ran concurrently; Cargo's build-directory lock
coordinated compilation. The commands and driver result are recorded below; no secret is included.

```sh
set -o pipefail
cargo fmt --all -- --check > /tmp/raised-rational-fmt.log 2>&1 || exit
cargo check -p three-dimensional-reps --all-features --all-targets --locked > /tmp/raised-rational-shared-check.log 2>&1 || exit
cargo check -p three-dimensional-reps --no-default-features --all-targets --locked > /tmp/raised-rational-minimal-check.log 2>&1 || exit
cargo check -p gammalooprs -p gammaloop-api -p gammaloop-integration-tests --all-targets --locked > /tmp/raised-rational-consumers-check.log 2>&1 || exit
cargo clippy -p three-dimensional-reps --all-features --all-targets --locked > /tmp/raised-rational-shared-clippy.log 2>&1 || exit
cargo clippy -p gammalooprs -p gammaloop-api -p gammaloop-integration-tests --all-targets --locked > /tmp/raised-rational-consumers-clippy.log 2>&1 || exit
cargo nextest run -p three-dimensional-reps --all-features --lib --locked -P test_gammaloop --test-threads 4 --retries 0 --no-fail-fast > /tmp/raised-rational-shared-tests.log 2>&1 &
raised_rational_shared_pid=$!
cargo nextest run -p gammalooprs --lib --test test_renormalization --locked -P test_gammaloop -E 'test(cff::) or test(graph::three_d_source::) or test(graph::lmb::) or test(uv::approx::) or test(uv::hedge_poset::) or test(integrands::process::) or test(integrands::evaluation::) or test(numerator::energy_degree::) or test(processes::process::tests::computed_uv_forest_process_export) or test(=scalar_pole_part) or test(=finite_part_ghost_2loop)' --test-threads 4 --retries 0 --no-fail-fast > /tmp/raised-rational-core-tests.log 2>&1 &
raised_rational_core_pid=$!
cargo nextest run -p gammaloop-api --lib --locked -P test_gammaloop --test-threads 4 --retries 0 --no-fail-fast > /tmp/raised-rational-api-tests.log 2>&1 &
raised_rational_api_pid=$!
cargo nextest run -p gammaloop-integration-tests --test test_cli --test test_runs --test test_evaluation_api --test test_gamma_star_ttx_nlo_acceptance --test uv --locked -P test_gammaloop -E '((binary(=test_cli) or test(lu_rust_get_integrand_info_reports_groups_orientations_lmbs_and_cuts) or test(test_3d_reps::) or test(repeated_masses::) or test(raised_cut_numerator_cancels_one_propagator_in_both_orientation_modes) or test(bench_cli_profiles_fixed_scalar_triangle_point_and_restores_settings) or test(integrate_writes_numerical_stability_histograms_for_scalar_triangle) or test(comparisons_reject_nonfinite_and_large_relative_errors) or test(default_scalar_3l_cross_section_inspects::) or binary(=test_gamma_star_ttx_nlo_acceptance) or (binary(=uv) and (test(=scalar_amplitudes_match_across_local_uv_routes) or test(=scalar_spectacles_integrated_matches_across_local_uv_routes))))) or (test(scalar_3l_cross_section_gl04_q5_spatial_square_inspects_match) or test(scalar_3l_cross_section_gl04_q5_temporal_square_inspects_match) or test(scalar_3l_cross_section_gl04_right_temporal_product_inspects_match) or test(scalar_3l_cross_section_gl04_right_owned_dot_inspects_match) or test(scalar_3l_cross_section_gl00_quadratic_energy_inspects_match::q7_squared) or test(=scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl04_inspects_match))' --test-threads 4 --retries 0 --no-fail-fast > /tmp/raised-rational-integrations-tests.log 2>&1 &
raised_rational_integrations_pid=$!
wait "$raised_rational_shared_pid"
raised_rational_shared_status=$?
wait "$raised_rational_core_pid"
raised_rational_core_status=$?
wait "$raised_rational_api_pid"
raised_rational_api_status=$?
wait "$raised_rational_integrations_pid"
raised_rational_integrations_status=$?
printf 'SHARED=%s CORE=%s API=%s INTEGRATIONS=%s\n' "$raised_rational_shared_status" "$raised_rational_core_status" "$raised_rational_api_status" "$raised_rational_integrations_status" > /tmp/raised-rational-status.log
exit $((raised_rational_shared_status || raised_rational_core_status || raised_rational_api_status || raised_rational_integrations_status))
```

Initial driver result: `SHARED=0 CORE=100 API=0 INTEGRATIONS=0`.
The approved fixture correction was then verified in the licensed Nix environment with the
same stack setting and zero retries:

```sh
cargo fmt --all
cargo fmt --all -- --check > /tmp/raised-rational-approved-affine-fmt.log 2>&1
cargo check -p gammalooprs --all-targets --locked > /tmp/raised-rational-approved-affine-check.log 2>&1
cargo clippy -p gammalooprs --all-targets --locked > /tmp/raised-rational-approved-affine-clippy.log 2>&1
cargo nextest run -p gammalooprs --lib --locked -P test_gammaloop -E 'test(=graph::three_d_source::tests::exact_taylor_vacuum_affine_circulation_preserves_energy_mapping)' --test-threads 4 --retries 0 --no-fail-fast > /tmp/raised-rational-approved-affine-tests.log 2>&1
```

Correction driver result: `RATIONAL_APPROVED_FIX_EXIT=0`.


## Earlier implementation results

Before the coefficient migration, all 807 distinct selected tests passed after the approved test-only corrections: 106 shared,
276 core, 375 API, 48 integrations and two Spenso tests. The earlier core run passed 275 of
276; its sole cache-oracle failure is resolved by the passing focused rerun below, which also
checks the analogous triangle case. The complete shared suite was rerun. Production code did
not change in that test-only completion; the coefficient migration subsequently reran
the affected shared, core, API and integration selections. Repeated tests and the untouched baseline are excluded from the total.
A full workspace run is not claimed.

| Run | Executed | Passed | Failed | Seconds | Session log |
| --- | --- | --- | --- | --- | --- |
| Shared crate after approved corrections | 106 | 106 | 0 | 10.240 | `/tmp/raised-impl-approved-shared-tests.log` |
| Core and forest selection before oracle correction | 276 | 275 | 1, subsequently resolved | 101.439 | `/tmp/raised-impl-owned-mass-core.log` |
| Corrected cache oracle and analogous triangle | 2 | 2 | 0 | 2.497 | `/tmp/raised-impl-approved-core-tests.log` |
| API library | 375 | 375 | 0 | 14.945 | `/tmp/raised-impl-owned-mass-api.log` |
| Six GL00/GL04 regressions | 6 | 6 | 0 | 445.424 | `/tmp/raised-impl-owned-mass-regressions.log` |
| Complementary integration selection | 42 | 42 | 0 | 466.762 | `/tmp/raised-impl-owned-mass-integrations.log` |
| Two affected Spenso tests | 2 | 2 | 0 | 5.035 | `/tmp/raised-impl-spenso-tests2.log` |
| Contour cleanup rerun | 1 | 1 | 0 | 0.024 | `/tmp/raised-impl-final-contour.log` |
| Positive-sign representation cleanup rerun | 1 | 1 | 0 | 0.058 | `/tmp/raised-impl-final-positive-sign.log` |
| Untouched baseline two-loop pole | 1 | 1 | 0 | 69.121 | `/tmp/raised-impl-baseline-ghost.log` |

Formatting and affected all-target Cargo checks/Clippy passed. The test-only completion
reran shared/core all-target checks and Clippy; its lint logs are
`/tmp/raised-impl-approved-shared-clippy.log` and `/tmp/raised-impl-approved-core-clippy.log`.
The production-code lint logs include
`/tmp/raised-impl-owned-mass-clippy.log`, `/tmp/raised-impl-final-shared-clippy.log`,
`/tmp/raised-impl-final-positive-sign-clippy.log` and `/tmp/raised-impl-spenso-clippy-final.log`.
The only reported notice is the dependency `proc-macro-error2` future-compatibility warning.

## Validation commands

The final-driver commands are transcribed below; equivalent shared/Spenso reproduction commands
follow. Long nextest expressions preserve the selected scope. The shell drivers ran independent
suites concurrently; the commands below can also be run individually.

### Scoped UV correction: check, Clippy, core and six physical regressions

```sh
cargo fmt --all -- --check

cargo check -p gammalooprs -p gammaloop-api -p gammaloop-integration-tests --all-targets --locked

cargo clippy -p gammalooprs -p gammaloop-api -p gammaloop-integration-tests --all-targets --locked

cargo nextest run -p gammalooprs --lib --test test_renormalization --locked -P test_gammaloop -E 'test(cff::) or test(graph::three_d_source::) or test(graph::lmb::) or test(uv::approx::) or test(uv::hedge_poset::) or test(integrands::process::) or test(integrands::evaluation::) or test(numerator::energy_degree::) or test(processes::process::tests::computed_uv_forest_process_export) or test(=scalar_pole_part) or test(=finite_part_ghost_2loop)' --test-threads 4 --retries 0 --no-fail-fast

cargo nextest run -p gammaloop-integration-tests --test test_runs --locked -P test_gammaloop -E 'test(scalar_3l_cross_section_gl04_q5_spatial_square_inspects_match) or test(scalar_3l_cross_section_gl04_q5_temporal_square_inspects_match) or test(scalar_3l_cross_section_gl04_right_temporal_product_inspects_match) or test(scalar_3l_cross_section_gl04_right_owned_dot_inspects_match) or test(scalar_3l_cross_section_gl00_quadratic_energy_inspects_match::q7_squared) or test(=scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl04_inspects_match)' --test-threads 4 --retries 0 --no-fail-fast
```

### Final production code: complementary API and integration selections

```sh
cargo nextest run -p gammaloop-api --lib --locked -P test_gammaloop --test-threads 4 --retries 0 --no-fail-fast

cargo nextest run -p gammaloop-integration-tests --test test_cli --test test_runs --test test_evaluation_api --test test_gamma_star_ttx_nlo_acceptance --test uv --locked -P test_gammaloop -E '(binary(=test_cli) or test(lu_rust_get_integrand_info_reports_groups_orientations_lmbs_and_cuts) or test(test_3d_reps::) or test(repeated_masses::) or test(raised_cut_numerator_cancels_one_propagator_in_both_orientation_modes) or test(bench_cli_profiles_fixed_scalar_triangle_point_and_restores_settings) or test(integrate_writes_numerical_stability_histograms_for_scalar_triangle) or test(comparisons_reject_nonfinite_and_large_relative_errors) or test(default_scalar_3l_cross_section_inspects::) or binary(=test_gamma_star_ttx_nlo_acceptance) or (binary(=uv) and (test(=scalar_amplitudes_match_across_local_uv_routes) or test(=scalar_spectacles_integrated_matches_across_local_uv_routes)))) and not (test(scalar_3l_cross_section_gl04_q5_spatial_square_inspects_match) or test(scalar_3l_cross_section_gl04_q5_temporal_square_inspects_match) or test(scalar_3l_cross_section_gl04_right_temporal_product_inspects_match) or test(scalar_3l_cross_section_gl04_right_owned_dot_inspects_match) or test(scalar_3l_cross_section_gl00_quadratic_energy_inspects_match::q7_squared) or test(=scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl04_inspects_match))' --test-threads 4 --retries 0 --no-fail-fast
```

### Final shared test diagnostic cleanup

```sh
cargo fmt --all -- --check

cargo check -p three-dimensional-reps --all-features --all-targets --locked

cargo clippy -p three-dimensional-reps --all-features --all-targets --locked

cargo nextest run -p three-dimensional-reps --all-features --lib --locked -P test_gammaloop -E 'test(ordinary_terminal_sources_match_signed_contours)' --test-threads 4 --retries 0
```

### Final positive-sign test representation cleanup

```sh
cargo fmt --all -- --check

cargo check -p gammalooprs --all-targets --locked

cargo clippy -p gammalooprs --all-targets --locked

cargo nextest run -p gammalooprs --lib --locked -P test_gammaloop -E 'test(=cff::generation::tests::projected_denominator_chains_preserve_occurrence_multiplicity)' --test-threads 4 --retries 0
```

### Shared suite and affected Spenso tests

```sh
cargo nextest run -p three-dimensional-reps --all-features --lib --locked -P test_gammaloop --test-threads 4 --retries 0 --no-fail-fast

cargo nextest run -p spenso --features shadowing --lib --locked -P test_gammaloop -E 'test(large_scaled_tensor_sum_contracts_without_broadcasting_coefficients) or test(grouped_atom_multi_contraction_preserves_metric_permutations_and_factors)' --test-threads 4 --retries 0

cargo clippy -p spenso --features shadowing --all-targets --locked
```

### Schema export and document build

```sh
GAMMALOOP_SCHEMA_PATH="$PWD/assets/schemas" cargo run -p gammaloop-api --bin gammaloop --locked -- -s /tmp/raised-impl-schema -n save schema

typst compile docs/architecture/raised-energy-cff-implementation.typ docs/architecture/raised-energy-cff-implementation.pdf

git diff --check
```

The schema export completed successfully and all 205 local schema references resolve.
Only `runhistory.json` changed, reflecting the template payload and refreshing existing
UV-profile schema drift. Snapshot audit compared all 198 archived blobs to the reviewed
originals; every file is identical. The audit ledger is `/tmp/raised-final-snapshot-move-audit.json`.

## Approved test-oracle corrections

The user authorized the prepared corrections on 2026-09-09 after the repository-required
approval request. All affected checks now pass. The three original failures were mistakes in
the replacement oracles; numerical tolerances remain unchanged.

| Test | Corrected outward contract |
| --- | --- |
| `cff::generation::tests::exact_cff_batch_reuses_owner_relabelled_sub_lmb_topology` | Compare complete production CFF values with all convention and on-shell factors resolved at two exact common shell points. Require a nonzero scalar residue. |
| `uv::approx::local_4d::tests::dod_one_triangle_keeps_separate_denominator_topologies` | Apply the same full production conversion to each Taylor term. Require at least one nonzero residue, allowing legitimate zero terms and cancellation in the sum. |
| `generation::causal_generation_tests::finite_pole_sampling_preserves_reversed_repeated_channel_moments` | Reverse endpoints along with the complete momentum signature and require graph validation. Keep the `1e-11` tolerance. |
| `generation::cff_tests::cff_generation_assigns_distinct_labels_to_distinct_energy_maps` | Distinct canonical residue maps require distinct complete public labels. Require a distinct-map witness without fixing label format, count, density or ordering. |

The following commands reran after the approved changes:

```sh
cargo fmt --all
cargo check -p three-dimensional-reps --all-features --all-targets --locked
cargo check -p gammalooprs --all-targets --locked
cargo clippy -p three-dimensional-reps --all-features --all-targets --locked
cargo clippy -p gammalooprs --all-targets --locked
cargo nextest run -p three-dimensional-reps --all-features --lib --locked -P test_gammaloop --test-threads 4 --retries 0 --no-fail-fast
cargo nextest run -p gammalooprs --lib --locked -P test_gammaloop -E 'test(exact_cff_batch_reuses_owner_relabelled_sub_lmb_topology) or test(dod_one_triangle_keeps_separate_denominator_topologies)' --test-threads 4 --retries 0 --no-fail-fast
```

No production code changed during this completion. Final formatting, whitespace, document
compilation and visual PDF checks also passed.

## Independent diagnosis of the mixed UV correction

Unchanged GL04 strict Arb comparisons initially disagreed by roughly `8.5749e-21` in the
imaginary part. Captured complete evaluator inputs and normalized function maps were identical.
Independent 90, 180 and 330 digit rational-program replays reproduced the discrepancy;
coefficient extraction isolated it to a one-loop integrated localization factor. Varying its
localization scale showed that the small default-scale error was systematic, not an Arb
precision limit.

The direct 3D active coefficient lacked the vacuum-mass hard scaling already present in the
four-dimensional forest Taylor operation when that operation encloses the coefficient owner.
Normalized localization kernels stay frozen; a disjoint coefficient stays fixed for that
particular Taylor operation but may become active for a later enclosing one. The final direct
correction retains the owner of each connected coefficient using transient
`mUV(owner)` applications. The Taylor operation scales contained owners, leaves disjoint ones
fixed, and rejects partial overlaps. Output copies restore ordinary `mUV`; stored sectors
retain their owners for later enclosing steps. Four-dimensional Taylor and integrated
coefficient code retain their established mass scaling.

A blanket direct mass substitution also failed GL00 q7-squared and the GL04 quadratic-pair
numerator, both of which contain disjoint divergent components. Fresh GL00 capture isolated
the mismatch to the bare cut contribution: threshold terms agreed, all paired Arb inputs
were identical, and all 43 normalized function definitions matched. Independent program
replay reproduced the generated-function difference. The scoped correction preserves the
component information needed at this boundary. Individual LU-order coefficient comparisons
were treated as diagnostics; they were not substituted for the complete-cut oracle.

The independent two-loop renormalization test was run unchanged on an untouched archive of
the pre-implementation baseline with `--locked --offline` and `INSTA_UPDATE=no`; it passed.
It also passes on the final code, together with depth three and all six discriminating
GL00/GL04 cases. Incorrect four-dimensional alternatives were discarded. Their interrupted broad runs are not included
in the successful validation totals, and no snapshots were accepted or updated.
