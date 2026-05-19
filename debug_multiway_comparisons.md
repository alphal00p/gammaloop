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

Status as of 2026-05-19: the latest pushed checkpoint is
`d9d2b81f3e8cf4224b155aa503a4bc479f5c4be8` and is the scalar cross-section
anchor. The current working tree contains one additional targeted fix for the
amplitude-like divergent scalar-bubble CFF expanded-4D local-UV bucket.

The fix is in `crates/gammalooprs/src/uv/approx/expanded_4d.rs` and changes
when CFF expanded sources with repeated active denominators are evaluated in
the confluent equal-energy basis. Ordinary CFF expanded sources still take the
confluent limit before residue selection. UV-leading LU local-series sources
do not, because they are already Laurent coefficients in the selected
Cutkosky coordinate and receive their CFF/LTD bridge from the LU residue
construction.

This is intended to repair the CFF local expanded-4D behavior for repeated
active denominators without touching the scalar cross-section LU bridge. It is
not a per-graph modifier and does not depend on loop count or graph name.

The current working tree has no graph-name branches, no loop-count sign
exceptions, no production use of `graph_from_signature`, and no temporary
`GAMMALOOP_TRACE_LTD_LU_SIGNS` diagnostics in Rust code.

## Validation Matrix

Targeted bucket fixed by this patch:

```text
env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(=inspect::divergent_bubble_local_uv_from_expanded_4d_inspects_match_cff_and_ltd) | test(=inspect::divergent_bubble_integrated_uv_from_expanded_4d_inspects_match_cff_and_ltd) | test(=inspect::cff_local_uv_from_expanded_4d_works_without_explicit_orientation_sum)' \
  --no-capture --retries 0

3 tests run: 3 passed, 264 skipped
```

Formatting and compilation:

```text
cargo fmt --check

passed
```

```text
env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 cargo check

passed
```

Protected scalar cross-section anchor after the targeted fix:

```text
env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(/scalar_3l_cross_section_inspects::slow/)' \
  --no-capture --retries 0

143 tests run: 143 passed (1 slow), 124 skipped
```

Full selected GammaLoop suite after the targeted fix:

```text
env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T just test_gammaloop

1131 tests run: 1126 passed, 5 failed, 271 skipped
```

The five failing tests are the remaining known buckets, not scalar
cross-section regressions:

- `profile_bulk::dotted_bubble_amp_bulk_profile_passes`
- `profile_bulk::double_dotted_bubble_amp_bulk_profile_passes`
- `profile_bulk::scalar_self_energy_amp_bulk_profile_passes`
- `generation::ltd_tests::ltd_equal_signature_vacuum_hexagon_matches_cff_through_quintic_numerator`
- `generation::ltd_tests::ltd_gl06_forward_with_initial_state_cut_square_energy_numerator_matches_cff`

The profile failures are amplitude-like IR/profile tests. They show sparse
zero/infinite-value samples in the individual original/CT fits and a failing
sum fit. They must be handled as their own bucket, not mixed into the
expanded-4D CFF local-UV fix.

The 3Drep failures are diagnostic CFF/LTD comparison mismatches:

- equal-signature vacuum hexagon:
  `cff=-0.026196853417673083`,
  `ltd=0.009207375767346093`;
- GL06 forward initial-state cut square numerator:
  `cff=-66.79480952028291`,
  `ltd=548.3436360239327`.

These diagnostic failures were also present at the clean pushed scalar-sweep
anchor and are not introduced by the current targeted fix.

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

Do not start from a broad edit. Pick one bucket and keep it isolated.

1. Commit and push the current targeted fix only if the staged diff is limited
   to the CFF repeated-source guard and this status note.
2. Next bucket candidate: the three profile-bulk amplitude-like tests. First
   determine whether the zero/infinite samples are caused by original/CT
   cancellation being probed term-by-term, by the profile evaluator choosing a
   singular sample sequence for repeated bubbles, or by a real missing local
   limit in the amplitude-like subtraction. Fix only the shared profile/local
   limit construction, then rerun:
   - the three profile-bulk tests plus `massless_triangle`,
   - the 3-test divergent-bubble inspect bucket,
   - the full 143 slow scalar cross-section sweep,
   - `just test_gammaloop`.
3. Separate bucket: the two 3Drep diagnostic CFF/LTD mismatches. They should be
   fixed through the diagnostic comparison entry point or shared graph/LMB
   coordinate construction. Do not redirect GammaLoop production CFF local-UV
   generation and do not use `graph_from_signature` in production logic.
4. After each fixed bucket, update this file with the exact commands and
   counts, then commit and push before moving to the next bucket.
5. Final acceptance after all buckets are green:
   - `cargo fmt`
   - `env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T cargo check`
   - `env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T just test_gammaloop`
   - full slow scalar cross-section sweep without `INSTA_FORCE_PASS`
   - broad 3Drep integration sweep
   - no generated `.snap.new` files unless deliberately accepted
   - commit, push, and update the PR body with the implemented changes and
     validation matrix.
