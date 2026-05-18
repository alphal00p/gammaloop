# Multi-Way Comparison Debug Status

This note tracks the scalar cross-section rich local inspect parity effort
across:

1. CFF with local UV counterterms from 3D expansions,
2. CFF with local UV counterterms from expanded 4D integrands,
3. LTD with local UV counterterms from expanded 4D integrands.

The goal is event-level parity, including `Original`, threshold counterterms,
full multiplicative factors, cut ids, event grouping, and weights. Overall
common phases are less important than relative parity between the three modes,
but no graph-specific or unexplained sign modifiers are allowed.

## Current Checkpoint

The scalar cross-section local inspect matrix is green at this checkpoint and
has been pushed at commit `c714ef9b`. That commit deliberately restores the
pure-CFF production path for true repeated denominator channels, while keeping
the denominator-edge-coordinate source map for pure-CFF numerator localization.
The previous automatic redirection of all repeated CFF graphs through the
confluent/LTD-like bounded path made the 3Drep repeated comparison matrix
greener, but it spoiled the slow scalar cross-section local-inspect sweep. The
scalar sweep is therefore the guardrail: repeated-channel CFF/LTD diagnostic
comparisons must be repaired through an explicit, centralized comparison normal
form, not by changing the GammaLoop production CFF local-UV basis.

The current working revision after `c714ef9b` only fixes the macOS build/link
environment coverage. Crates that can link Symbolica/GMP/MPFR-backed tests now
honor `EXTRA_MACOS_LIBS_FOR_GNU_GCC=T` themselves, including subcrates that
previously had no `build.rs`. This removes the need for manual `RUSTFLAGS=-L...`
or `-l...` workarounds when running `just test_gammaloop`.

The 3Drep CFF/LTD diagnostic path now does exactly that.  The centralized
entry point is `generate_cff_ltd_comparison_expression`.  For non-repeated
graphs, and for non-CFF requested representations, it falls through to ordinary
`generate_3d_expression`.  For true repeated CFF channels it uses the common
derivative-free confluent normal form used by the LTD-like comparison basis and
disables the pure-CFF duplicate-signature excess convention.  That duplicate
sign belongs to the independent pure-CFF residue/source basis; it is not part
of the shared CFF/LTD diagnostic projection.  This affects `3Drep
test-cff-ltd`, `eval::compare_cff_ltd`, and the old Python probe matrix only.
It does not redirect GammaLoop production CFF local-UV generation.

The already-committed scalar sign work centralizes active repeated-channel
metadata in `Graph::active_three_d_repeated_channels()` and derives all
repeated-channel bridges from that single source:

- duplicate-signature excess for the global 3D sign exponent;
- repeated-pole residue bridge signs;
- normalized raised-edge support signs;
- the ordinary full-source same-routing duplicate bridge used by simple LTD LU
  direct-original residues.

The extra direct-original bridge is deliberately narrow. It applies only to
unselected active repeated channels whose copies have the same canonical
routing, where the full-source LTD residue must still be compared to the CFF
projection convention of the same source. Mixed-routing channels already carry
their source orientation in the generated repeated-LTD routing. When such a
channel is selected, the repeated-pole derivative bridge and normalized
edge-support bridge carry the remaining orientation.

This keeps the scalar cross-section sign logic graph/LMB-derived. The current
checkpoint still has no GL-specific branches, no loop-number fudge factors, no
production use of `graph_from_signature`, and no temporary
`GAMMALOOP_TRACE_LTD_LU_SIGNS` diagnostics in the touched Rust code.

The `143/143` scalar sweep is not a diagnostic self-comparison. The GammaLoop
cross-section tests still generate the three production modes separately:
pure-CFF local UV from 3D expansions, pure-CFF local UV from expanded 4D
integrands, and LTD local UV from expanded 4D integrands. The helper
`generate_cff_ltd_comparison_expression` is wired only into 3Drep diagnostics
(`test-cff-ltd`, `eval::compare_cff_ltd`, and the old Python probe matrix), not
into GammaLoop production CFF local-UV generation. `INSTA_FORCE_PASS=1` was used
only to prevent snapshot write gating from hiding numerical mismatches; the
generated `.snap.new` files were deleted and not accepted.

## Validation Matrix

Commands below were run with `INSTA_FORCE_PASS=1` so snapshot drift did not
hide numerical parity failures. Generated `.snap.new` files were removed and
were not accepted as new references.

Focused repeated-channel guardrail:

```text
env EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 RUST_MIN_STACK=33554432 \
  RUSTFLAGS=-L/opt/local/lib/libgcc INSTA_FORCE_PASS=1 \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E '<21 focused scalar inspect tests>' --no-capture --retries 0

21 tests run: 21 passed, 246 skipped
```

Focused crossing and mixed repeated-channel set:

```text
env EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 RUST_MIN_STACK=33554432 \
  RUSTFLAGS=-L/opt/local/lib/libgcc INSTA_FORCE_PASS=1 \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E '<GL21, GL38, GL40, GL46 focused scalar inspect tests>' \
  --no-capture --retries 0

8 tests run: 8 passed, 259 skipped
```

Full slow scalar cross-section local inspect sweep:

```text
env EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 RUST_MIN_STACK=33554432 \
  RUSTFLAGS=-L/opt/local/lib/libgcc INSTA_FORCE_PASS=1 \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(/scalar_3l_cross_section_inspects::slow/)' \
  --no-capture --retries 0

143 tests run: 143 passed, 124 skipped
```

The checkpoint was revalidated on 2026-05-18 after the pure-CFF restoration and
the denominator-edge-coordinate source map change:

```text
find tests/tests/test_runs/snapshots -name '*.snap.new' -delete

env EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 RUST_MIN_STACK=33554432 \
  RUSTFLAGS=-L/opt/local/lib/libgcc INSTA_FORCE_PASS=1 \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(/scalar_3l_cross_section_inspects::slow/)' \
  --no-capture --retries 0

143 tests run: 143 passed (2 slow), 124 skipped
```

The previous sign-only GL38/GL46 failures are now green. Earlier GL24 and GL35
local-UV drift failures also remain green.

Most recent full selected GammaLoop suite status before the 3Drep diagnostic
fix:

```text
just test_gammaloop

1131 tests run: 1089 passed, 42 failed, 271 skipped
```

Those failures are outside the now-green slow scalar cross-section matrix. The
important buckets are:

- forward-scattering LU generation tests where the local-series external
  coordinate selector finds two equivalent initial-state candidates with the
  same external signature;
- the divergent/amplitude-like scalar bubble local/integrated UV inspect tests;
- generated forward cross-section smoke tests that currently hit the same
  coordinate ambiguity;
- amplitude profile-bulk tests;
- spin-sum cross-section generation tests that also stop at the coordinate
  ambiguity;
- one default quartic scalar snapshot drift when `INSTA_FORCE_PASS` is not set.

The checkpoint should therefore be read as a scalar cross-section multi-way
local-inspect parity checkpoint, not as a full-suite clean point.

Most recent full selected GammaLoop suite status after the macOS build-link
cleanup:

```text
env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T \
  just test_gammaloop -- --no-fail-fast --final-status-level=fail \
  --status-level=fail

1131 tests run: 1101 passed, 30 failed, 271 skipped
```

The current 30 failures group as follows:

- `test_evaluation_api` LU event/integration metadata tests: 12 failures. These
  all stop during generation with
  `Cannot determine a unique LTD LU cut local-series external coordinate`.
  The ambiguous candidates are the two incoming external half-edges of a
  forward-scattering cut carrying the same external energy shift.
- `test_runs::differential` LU differential tests: 4 failures. Same
  forward-scattering local-series coordinate ambiguity.
- `test_runs::spin_sums`: 2 failures. Same generation ambiguity in the
  spin-summed cross-section construction.
- Generated forward cross-section inspect tests: 3 aborts:
  `ltd_generated_forward_cross_section_quartic_numerator_matches_cff_and_energy_trade`,
  `ltd_generated_forward_cross_section_threshold_inspects_match_cff`, and
  `ltd_generated_forward_cross_section_threshold_and_4d_uv_inspects_match_cff`.
  These are expected to share the same forward-scattering coordinate root until
  proven otherwise by a targeted trace.
- Divergent/amplitude-like scalar bubble inspect/profile tests: 6 failures:
  `cff_local_uv_from_expanded_4d_works_without_explicit_orientation_sum`,
  `divergent_bubble_local_uv_from_expanded_4d_inspects_match_cff_and_ltd`,
  `divergent_bubble_integrated_uv_from_expanded_4d_inspects_match_cff_and_ltd`,
  `dotted_bubble_amp_bulk_profile_passes`,
  `double_dotted_bubble_amp_bulk_profile_passes`, and
  `scalar_self_energy_amp_bulk_profile_passes`.
- One scalar cross-section default snapshot drift:
  `scalar_3l_cross_section_gl24_quartic_energy_inspects_match::q1_quartic`.
  The visible drift is numerical, e.g. `8.39312238214002094e-12` to
  `8.38188823800568133e-12`, and must not be accepted until its source is
  understood.
- 3Drep unit diagnostics: 2 failures:
  `ltd_equal_signature_vacuum_hexagon_matches_cff_through_quintic_numerator`
  and
  `ltd_gl06_forward_with_initial_state_cut_square_energy_numerator_matches_cff`.
  These are not part of the GammaLoop production scalar sweep, but they are
  useful probes for the remaining CFF/LTD diagnostic and forward-cut theory.

The build-link fix was validated independently with:

```text
env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T \
  RUST_MIN_STACK=33554432 \
  cargo test -p idenso --profile dev-optim --no-run

Finished dev-optim profile
```

## 3Drep Status At This Checkpoint

The pure-CFF source-map correction fixes the broad non-repeated and
initial-state-cut 3Drep mismatches without changing the scalar sweep:

```text
env EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 RUST_MIN_STACK=33554432 \
  RUSTFLAGS=-L/opt/local/lib/libgcc \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(=test_3d_reps::cli_aa_aa_box_and_double_box_multibackend_threedrep_comparisons) | test(=test_3d_reps::cli_physical_hexagon_high_power_energy_sectors_match_ltd)' \
  --no-capture --retries 0

2 tests run: 2 passed, 265 skipped
```

The imported `box_pow3` repeated-channel fixture now has direct CFF/LTD
diagnostic parity through the centralized comparison normal form: both CFF and
LTD serialize to the derivative-free repeated-channel comparison expression
with `17` orientations and `24` nodes.  The fixture also still asserts the
actual five-variant LTD repeated-channel structure.

```text
env EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 RUST_MIN_STACK=33554432 \
  RUSTFLAGS=-L/opt/local/lib/libgcc \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(=test_3d_reps::cli_imported_box_pow3_3drep_test_uses_gammaloop_graph_path)' \
  --no-capture --retries 0

1 test run: 1 passed, 266 skipped
```

The old Python 3Drep case matrix is now green:

```text
old Python 3Drep case matrix: 310 ok, 0 failed, 310 total
```

The broad 3Drep integration sweep is also green after the diagnostic entry
point change:

```text
env EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 RUST_MIN_STACK=33554432 \
  RUSTFLAGS=-L/opt/local/lib/libgcc \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(/test_3d_reps/)' --no-capture --retries 0

17 tests run: 17 passed (8 slow), 250 skipped
```

## Theory Anchor

`docs/3Dreps/generalised_ltd.tex` states that local unitarity cut selection and
threshold counterterm construction should use a representation-neutral
canonical causal `E`-surface catalogue. CFF exposes this catalogue directly.
LTD may contain additional `H`-surfaces, but those are not physical threshold
surfaces after the orientation sum.

The implementation identifies LU cuts and threshold residues once from
graph/LMB/linnet and canonical `E`-surface data, then maps the selected residue
variables back to the requested representation. Simple, raised, confluent,
threshold, and expanded-4D UV sources should not use separate ad-hoc sign
algorithms.

Edges with `is_cut` are forward-scattering cut edges. For graph manipulation
they should behave as two external half-edges with the same momentum, one
incoming and one outgoing. LU/Cutkosky and threshold support normalization must
therefore treat them as cut/external boundaries rather than ordinary internal
propagators.

## Sensitive Code Paths

- `Graph::active_three_d_repeated_channels()` is the single repeated-channel
  descriptor used for global 3D sign excess, residue bridges, source duplicate
  bridges, and raised-edge support normalization.
- `CrossSectionGraph::residue_selector_for_raised_cut_group` builds the LU
  residue plan and is the main cross-section entry point to keep centralized.
- `CrossSectionGraph::lu_cut_edge_sets_with_cutkosky_signs` maps raised LU
  groups to positive-energy Cutkosky edge support and signs.
- `CrossSectionGraph::simple_ltd_lu_cut_local_coordinate_signs` computes the
  generated selected `E`-surface sign, the local-series denominator sign, and
  the surface-family bridge from the resolved LMB coordinate determinant.
- `CrossSectionGraph::simple_lu_cut_has_mixed_repeated_channel_contact`
  identifies the support-local case where a simple cut bridges multiple
  repeated channels.
- `ResidueSelector` stores selected signs, local-series signs, ordinary residue
  prefactors, direct-original prefactors, and threshold prefactors.
- `cff/expression.rs::ltd_lu_local_series_coefficients_from_parametric_atom`
  consumes the LTD LU local-series coordinates and signs.
- `uv/approx/mod.rs` and `uv/approx/expanded_4d.rs` dispatch threshold residue
  selection by representation and must remain LTD-expression-facing for LTD
  output.
- `three-dimensional-reps::generate_cff_ltd_comparison_expression` is the only
  diagnostic CFF/LTD comparison entry point.  It is intentionally not the
  GammaLoop production CFF local-UV generator.

## Remaining Follow-Up

The scalar cross-section slow matrix remains the protected anchor. The remaining
work must proceed in dependency order, with explicit regression gates after
each principled change rather than by toggling signs graph by graph.

1. Freeze the anchor before any more algorithmic edits.
   - Run the slow scalar cross-section sweep on the current build-link revision
     after deleting `.snap.new` files.
   - Inspect the multi-way harness once more to confirm it still compares the
     three independent GammaLoop production modes and is not routed through the
     3Drep diagnostic common projection.
   - Keep `INSTA_FORCE_PASS=1` only as a snapshot gating bypass; never stage
     `.snap.new` files unless a drift is deliberately accepted.
2. Repair the forward-scattering LU coordinate ambiguity first.
   - Treat every `is_cut` edge as two external half-edges with the same momentum,
     one incoming and one outgoing.
   - Build one canonical external half-edge orbit for equal-energy
     forward-scattering candidates, using graph/LMB/linnet data rather than
     external edge index order.
   - Make the local coordinate construction choose that canonical orbit instead
     of requiring a unique raw candidate.
   - Targeted gate:
     `lu_rust_default_clustered_pdgs_match_explicit_massless_qcd_list` plus one
     differential LU test and one spin-sum test.
   - Regression gate immediately after the targeted fix: full slow scalar
     `143/143`.
3. Re-run all currently coordinate-blocked tests.
   - The expected bucket is the 12 evaluation API tests, 4 differential tests,
     2 spin-sum tests, and likely the 3 generated forward inspect tests.
   - If any generated forward inspect test remains after coordinate selection,
     trace threshold-residue selection and Cutkosky support normalization before
     touching signs.
   - Regression gate: full slow scalar `143/143` again if any code touched
     cross-section generation, residue selectors, LU coordinates, or threshold
     selection.
4. Address amplitude-like scalar bubble UV behavior separately.
   - First evaluate the original divergent bubble in CFF and LTD at the same
     sample and confirm the original graph reference agrees.
   - Then compare CFF 3D local UV, CFF expanded-4D local UV, and LTD expanded-4D
     local UV before event assembly, term by term.
   - Classify the mismatch as source-basis projection, duplicate/contact
     channel collapse, local/integrated UV normalization, or numerator
     localization.
   - Fix only the shared mathematical source of the mismatch. Do not modify the
     now-green scalar LU residue bridge unless the bubble exposes the same
     graph/LMB invariant.
   - Regression gate after any UV-generation fix: full slow scalar `143/143`,
     then the divergent-bubble inspect/profile bucket.
5. Revisit the two remaining 3Drep unit diagnostics after the production
   forward-cut and bubble issues are understood.
   - The equal-signature vacuum hexagon likely probes diagnostic CFF/LTD
     high-power numerator canonicalization.
   - The GL06 initial-state-cut square numerator likely probes the same
     forward-cut half-edge coordinate/orientation theory as the GammaLoop
     coordinate ambiguity, but in the standalone 3Drep diagnostic path.
   - Fix through the centralized diagnostic entry point or the shared graph/LMB
     coordinate construction only; do not route GammaLoop production CFF through
     diagnostic projection.
6. Final acceptance sequence after every bucket is green:
   - `cargo fmt`
   - `env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T cargo check`
   - `env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T just test_gammaloop`
   - full slow scalar cross-section sweep `143/143`
   - broad 3Drep integration sweep
   - remove all generated `.snap.new` files unless deliberately accepted
   - update this markdown and the theory note for any mathematical construction
     that changed.
