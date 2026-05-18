# UV scalar profile benchmark notes

Date: 2026-05-18

This note summarizes the `dev-optim` benchmark and sampling pass for the tests that had previously shown stack-sensitive behavior in the default cargo `test` profile.

## Scope

The tests checked were:

- `uv::tests::scalars_profile`
- `uv::tests::scalars_profile_new`
- `uv::tests::spinney_partial_cmp_is_equal_for_identical_subgraphs`
- `processes::amplitude::test::generation_orientation_pattern_filters_evaluator_orientations`

These tests passed in `dev-optim`; the investigation here was about where the slow runtime is spent.

## Benchmark commands

Representative commands used:

- `cargo nextest run --profile test_gammaloop --cargo-profile dev-optim -p gammalooprs -E 'test(=uv::tests::scalars_profile)' --retries 0`
- `cargo nextest run --profile test_gammaloop --cargo-profile dev-optim -p gammalooprs -E 'test(=uv::tests::scalars_profile_new)' --retries 0`
- `cargo nextest run --profile test_gammaloop --cargo-profile dev-optim -p gammalooprs -E 'test(=uv::tests::spinney_partial_cmp_is_equal_for_identical_subgraphs)' --retries 0`
- `cargo nextest run --profile test_gammaloop --cargo-profile dev-optim -p gammalooprs -E 'test(=processes::amplitude::test::generation_orientation_pattern_filters_evaluator_orientations)' --retries 0`

Sampling was done with macOS `sample` against the `dev-optim` test binary while running targeted tests directly.

## Benchmark results

| Test | Result | Observed runtime |
| --- | ---: | ---: |
| `uv::tests::scalars_profile` | passed | ~101s in a clean timing run; a later rerun under noisier load took ~123s |
| `uv::tests::scalars_profile_new` | passed | ~99.5s |
| `uv::tests::spinney_partial_cmp_is_equal_for_identical_subgraphs` | passed | ~84.7s |
| `processes::amplitude::test::generation_orientation_pattern_filters_evaluator_orientations` | passed | ~0.31s |

The slow set is therefore the three scalar UV tests. The orientation-pattern test is not materially slow in `dev-optim`.

## High-level timing breakdown

A timestamped debug trace of `uv::tests::scalars_profile_new` showed the following approximate breakdown:

| Phase | Approximate time | Notes |
| --- | ---: | --- |
| UV expression / forest preprocessing | ~67–69s | Dominant cost |
| Evaluator construction | ~17–20s | Second-largest cost |
| Actual UV profile evaluation loop | ~14–19s | Smaller than construction |
| `Spinney::partial_cmp` itself | negligible | The spinney test is slow because it builds the same UV scalar amplitude first |

## Main hotspot

The largest single gap in both `scalars_profile` and `scalars_profile_new` occurs during UV final-integrand construction:

`CutForests::compute` → `Approximation::compute` → `Approximation::final_integrand`

The trace repeatedly showed a large gap between:

- `Numerator constructed, numerator: 1`
- `Integrand before parsing for S_GV⊛GS*Top(S_GS⊛44*Top(S_44⊛0)) for dod4: ...`

That single gap was about:

- ~50.6s in `scalars_profile_new`
- ~52.7s in `scalars_profile`

The sampled stack during that gap was dominated by symbolic simplification:

- `Approximation::final_integrand`
- `simplify_color`
- `idenso::color::color_simplify_impl`
- `collect_chains`
- Symbolica `collect_symbol` / polynomial conversion / zero testing
- recursive `symbolica::evaluate_impl`
- `symbolica::zero_test_impl`
- `symbolica::to_polynomial_in_vars_impl`
- GMP / `rug` rational allocation and conversion routines

Conclusion: the main cost is symbolic final-integrand construction and simplification, especially color simplification / chain collection and Symbolica polynomial/zero-test recursion for the `dod4` `S_GV⊛GS*Top(S_GS⊛44*Top(...))` case.

## Extracting color-simplification inputs

`Approximation::final_integrand` now supports a focused dump mode for the exact expressions passed to the hot `simplify_color()` call. Set `GAMMALOOP_DUMP_UV_COLOR_SIMPLIFY_INPUTS` to an output directory and rerun the scalar profile test. Use an absolute path when you want the files in a predictable location, because test runners may use a package-local working directory.

For example:

- `GAMMALOOP_DUMP_UV_COLOR_SIMPLIFY_INPUTS=<repo>/target/uv-color-simplify-inputs cargo nextest run --profile test_gammaloop --cargo-profile dev-optim -p gammalooprs -E 'test(=uv::tests::scalars_profile_new)' --retries 0`

For each final-integrand term, this writes:

- `*.sym`: the parser-oriented `Atom::to_plain_string()` expression after `GS.dim -> 4` and immediately before `simplify_color()`.
- `*.meta.txt`: graph name, UV topological order, `dod`, term index, approximation label, and the expression file path.

The target hotspot should be identifiable in the metadata as `dod = 4` with an approximation containing `S_GV⊛GS*Top(S_GS⊛44*Top(S_44⊛0))`. In the extraction run from this note, the matching expression files were:

| Expression file | Size |
| --- | ---: |
| `target/uv-color-simplify-inputs/uv-final-integrand-color-input-banana-topo29-dod4-term0.sym` | 517,558 bytes |
| `target/uv-color-simplify-inputs/uv-final-integrand-color-input-banana-topo36-dod4-term0.sym` | 116,635 bytes |
| `target/uv-color-simplify-inputs/uv-final-integrand-color-input-banana-topo39-dod4-term0.sym` | 30,203 bytes |
| `target/uv-color-simplify-inputs/uv-final-integrand-color-input-banana-topo45-dod4-term0.sym` | 30,203 bytes |

The largest hotspot expression is also checked in as `crates/gammalooprs/tests/resources/uv_color_simplify/banana_topo29_dod4_term0.sym`. Run the ignored regression/hotspot test with:

- `cargo test --profile dev-optim -p gammalooprs --test test_uv_color_simplify_dump -- --ignored --exact loads_uv_scalar_profile_hotspot_dump_and_simplifies_color --nocapture`

That test initializes GammaLoop/Symbolica, loads the fixture from disk, parses it, and calls `simplify_color()`.

## Evaluator-build hotspot

The next-largest cost occurs after `Orientation Parametric integrand` logging, during evaluator construction.

The sampled stack included:

- `Amplitude::build_integrand`
- `AmplitudeGraph::generate_term_for_graph`
- `EvaluatorStack::new_with_timings`
- `EvaluatorStack::new_iterative`
- `EvaluatorStack::new_single_parametric`
- Symbolica `to_evaluator`
- `optimize_horner_scheme_multiple`
- `ExpressionEvaluator::remove_common_pairs`
- `ExpressionEvaluator::merge`

Conclusion: after UV expressions are built, the next bottleneck is Symbolica evaluator generation, especially Horner optimization and common-pair removal/merging.

## Debug logging impact

The direct traced run of one scalar profile test produced roughly 307 MB of debug output and about 50k timestamped debug events.

`init_test_tracing()` configures debug logging for tests, and some debug paths build expensive previews before logging. In particular, `Approximation::final_integrand` constructs a debug preview using operations such as replacement, factor collection, numerator collection, and `log_print(...)` in the same hot symbolic region.

This does not fully explain the runtime—the CPU sample still shows genuine Symbolica simplification work—but debug preview construction likely amplifies both runtime and memory pressure.

## Summary

The slow scalar tests are slow because they all call `build_uv_scalars_amplitude(...)`.

For `spinney_partial_cmp_is_equal_for_identical_subgraphs`, the actual comparison is not the expensive part; almost all runtime is setup.

The primary bottleneck is `Approximation::final_integrand` symbolic simplification, especially `simplify_color` / chain collection / Symbolica polynomial and zero-test recursion.

The secondary bottleneck is evaluator construction in `EvaluatorStack::new_with_timings`, mostly Symbolica `to_evaluator`, Horner optimization, and common-subexpression/common-pair processing.
