# Multi-Way Comparison Debug Status

This note tracks the scalar cross-section rich local inspect parity effort
across:

1. CFF with local UV counterterms from 3D expansions,
2. CFF with local UV counterterms from expanded 4D integrands,
3. LTD with local UV counterterms from expanded 4D integrands.

The invariant being protected is event-level parity of `Original`, threshold
counterterms, full multiplicative factors, cut ids, event grouping, and final
weights. Overall common phases are less important than relative parity between
the three modes, but graph-specific sign modifiers and loop-number fudge
factors are not allowed.

## Working Rule

The active debugging discipline is deliberately conservative:

1. Start from the latest pushed checkpoint with the largest confirmed passing
   set, especially the full 143-test slow scalar cross-section sweep.
2. Pick exactly one related failing bucket.
3. Apply only the minimal first-principles fix for that bucket.
4. As soon as the bucket is fixed, rerun the protected regression gates.
5. If the gates do not regress, document the fix, commit, and push. If they do
   regress, keep the patch small enough to understand and repair or replace it
   before moving on.

This file should stay in sync with that process so each commit has a clear
known-good test matrix and the next pass does not mix unrelated failures.

## Current Checkpoint

Status as of 2026-05-19: the latest pushed checkpoint before this pass is
`3e2fe0d4f6606ab34d28e0a5e72d0b6477b900ff`, which had the protected scalar
cross-section sweep green and only two selected-suite failures left. This pass
fixes exactly that remaining 3Drep diagnostic bucket:

- `generation::ltd_tests::ltd_equal_signature_vacuum_hexagon_matches_cff_through_quintic_numerator`
- `generation::ltd_tests::ltd_gl06_forward_with_initial_state_cut_square_energy_numerator_matches_cff`

The code fix is limited to
`crates/three-dimensional-reps/src/generation.rs`. In the LTD-style
finite-pole contact branch, the lower sector is generated as an LTD object, so
the contact quotient enters with its algebraic sign. The extra one-edge pinch
sign is specific to lifting a lower CFF denominator tree back into a pure-CFF
source/sink orientation and is not part of the LTD lower-sector comparison.

The corresponding derivation note was added to
`docs/3Dreps/generalised_ltd.tex`. The implementation remains free of
graph-name branches, loop-count sign exceptions, production use of
`graph_from_signature`, and temporary `GAMMALOOP_TRACE_LTD_LU_SIGNS`
diagnostics.

## Validation Matrix

Targeted two-test 3Drep diagnostic bucket:

```text
env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 \
  cargo nextest run -p three-dimensional-reps --features eval \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(=generation::ltd_tests::ltd_equal_signature_vacuum_hexagon_matches_cff_through_quintic_numerator) | test(=generation::ltd_tests::ltd_gl06_forward_with_initial_state_cut_square_energy_numerator_matches_cff)' \
  --no-capture --retries 0

2 tests run: 2 passed
```

Broad 3Drep eval sweep:

```text
env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 \
  cargo nextest run -p three-dimensional-reps --features eval \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  --no-capture --retries 0

60 tests run: 60 passed
```

Focused scalar cross-section regression bucket that caught the accidental
pure-CFF sign edit during this pass:

```text
env EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 RUSTFLAGS=-L/opt/local/lib/libgcc \
  INSTA_FORCE_PASS=1 \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(=scalar_3l_cross_section_inspects::slow::quadratic_energy_numerators::scalar_3l_cross_section_gl22_quadratic_energy_inspects_match::q1_squared) | test(=scalar_3l_cross_section_inspects::slow::quadratic_energy_numerators::scalar_3l_cross_section_gl22_quadratic_energy_inspects_match::q7_squared) | test(=scalar_3l_cross_section_inspects::slow::quadratic_energy_numerators::scalar_3l_cross_section_gl24_quadratic_energy_inspects_match::q1_squared) | test(=scalar_3l_cross_section_inspects::slow::quadratic_energy_numerators::scalar_3l_cross_section_gl24_quadratic_energy_inspects_match::q7_squared) | test(=scalar_3l_cross_section_inspects::slow::quadratic_energy_numerators::scalar_3l_cross_section_gl35_quadratic_energy_inspects_match::q1_squared) | test(=scalar_3l_cross_section_inspects::slow::quadratic_energy_numerators::scalar_3l_cross_section_gl40_quadratic_energy_inspects_match::q1_squared) | test(=scalar_3l_cross_section_inspects::slow::quadratic_energy_numerators::scalar_3l_cross_section_gl40_quadratic_energy_inspects_match::q7_squared) | test(=scalar_3l_cross_section_inspects::slow::quadratic_energy_numerators::scalar_3l_cross_section_gl48_quadratic_energy_inspects_match::q1_squared) | test(=scalar_3l_cross_section_inspects::slow::quadratic_energy_numerators::scalar_3l_cross_section_gl48_quadratic_energy_inspects_match::q7_squared)' \
  --no-capture --retries 0

9 tests run: 9 passed
```

Protected scalar cross-section sweep:

```text
env EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 RUSTFLAGS=-L/opt/local/lib/libgcc \
  INSTA_FORCE_PASS=1 \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(/scalar_3l_cross_section_inspects::slow/)' \
  --no-capture --retries 0

143 tests run: 143 passed, 124 skipped
```

Formatting, whitespace, and compilation:

```text
cargo fmt
git diff --check
env EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 cargo check

passed
```

Full selected GammaLoop suite:

```text
env EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 just test_gammaloop

1131 tests run: 1131 passed, 271 skipped
```

As of this checkpoint, there are no known remaining multi-way comparison
failures in the selected GammaLoop suite or the protected slow scalar
cross-section sweep.

## Sensitive Construction Points

- `Graph::active_three_d_repeated_channels()` is the single repeated-channel
  descriptor used for global 3D sign excess, residue bridges, source duplicate
  bridges, and raised-edge support normalization.
- `CrossSectionGraph::residue_selector_for_raised_cut_group` is the central
  cross-section entry point for LU residue planning. All simple, raised,
  repeated, threshold, and expanded-4D UV residue selectors must flow through
  this graph/LMB-derived data instead of separate per-case sign algorithms.
- `CrossSectionGraph::lu_cut_edge_sets_with_cutkosky_signs` maps raised LU
  groups to positive-energy Cutkosky edge support and signs.
- `CrossSectionGraph::simple_ltd_lu_cut_local_coordinate_signs` computes the
  generated selected `E`-surface sign, local-series denominator signs, and
  LU surface-family bridge from the resolved coordinate determinant.
- `ResidueSelector` stores selected signs, local-series signs, ordinary
  residue prefactors, direct-original prefactors, and threshold prefactors.
- `cff/expression.rs::ltd_lu_local_series_coefficients_from_parametric_atom`
  consumes the LTD LU local-series coordinates and signs. The expressions
  passed to it for LTD production are LTD expressions.
- `uv/approx/mod.rs` and `uv/approx/expanded_4d.rs` dispatch threshold residue
  selection by representation and must remain LTD-expression-facing for LTD
  output.
- `Graph::cff_source_has_repeated_active_denominators()` is the CFF-side guard
  for enabling confluent source evaluation on actual repeated active
  denominators. It uses the same graph source, initial-state cut handling, and
  preserved 4D denominator normalization as production 3D-expression
  generation.
- `three-dimensional-reps::generate_cff_ltd_comparison_expression` is the
  diagnostic CFF/LTD comparison entry point. It is intentionally not the
  GammaLoop production CFF local-UV generator.

The relevant theory note for the repeated/confluent basis distinction is
`docs/3dreps/generalised_ltd.tex`, especially the diagnostic normal-form and
embedding discussion. Repeated active denominators require a confluent CFF
limit for ordinary expanded sources; generated LU local-series sources are
already coordinates of a selected local residue and must not be conflated with
that ordinary source basis.

## Remaining Plan

No multi-way comparison bucket is currently known to be failing. Future changes
to CFF/LTD signs, repeated-channel handling, LU cut selection, threshold
residue selection, or finite-pole contact completion should keep the same
workflow:

1. Start from this green checkpoint.
2. Change one related bucket at a time.
3. Rerun the targeted bucket, the protected 143-test slow scalar sweep, and
   `just test_gammaloop` before committing.
4. Update this file with exact commands and counts before pushing a new
   checkpoint.
