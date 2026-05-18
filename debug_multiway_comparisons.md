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

The scalar cross-section local inspect matrix is green at this checkpoint. The
last source change centralizes active repeated-channel metadata in
`Graph::active_three_d_repeated_channels()` and derives all repeated-channel
bridges from that single source:

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

This keeps the scalar cross-section sign logic graph/LMB-derived. There are no
GL-specific branches, no loop-number fudge factors, no production use of
`graph_from_signature`, and no temporary `GAMMALOOP_TRACE_LTD_LU_SIGNS`
diagnostics in the touched Rust code.

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

The previous sign-only GL38/GL46 failures are now green. Earlier GL24 and GL35
local-UV drift failures also remain green.

Full selected GammaLoop suite status before the 3Drep fixture refresh:

```text
just test_gammaloop

1131 tests run: 1096 passed, 35 failed, 271 skipped
```

Those failures are outside the now-green slow scalar cross-section matrix. The
important remaining buckets are:

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

## 3Drep Broad Revisit

The broad 3Drep tests were rerun after the CFF repeated/high-power fixes. The
remaining failure was not numerical parity: the imported `box_pow3` fixture
still expected the old pure-CFF shape (`62` CFF orientations). The fixed
generator now routes repeated high-power CFF through the same derivative-free
bounded repeated-channel expression as LTD for that case, giving `17` CFF
orientations, `17` LTD orientations, and identical CFF/LTD expression JSON.
The fixture now asserts that shared expression explicitly.

Validated 3Drep commands:

```text
cargo test -p three-dimensional-reps --features eval -- --nocapture

60 passed, 0 failed
```

```text
env EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 RUST_MIN_STACK=33554432 \
  RUSTFLAGS=-L/opt/local/lib/libgcc \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(/test_3d_reps/)' --no-capture --retries 0

17 tests run: 17 passed, 250 skipped
```

The old Python matrix remains green inside that run:

```text
old Python 3Drep case matrix: 310 ok, 0 failed, 310 total
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

## Remaining Follow-Up

The scalar cross-section matrix no longer has a known multi-way local-inspect
parity failure. The remaining work should proceed in dependency order rather
than by toggling signs graph by graph.

1. Fix the LU coordinate ambiguity generically. When multiple external
   half-edges represent the same forward-scattering energy shift, the selected
   coordinate should be chosen from a canonical external half-edge orbit, not by
   failing uniqueness. This belongs in the local coordinate construction, using
   graph/LMB data and the `is_cut` half-edge interpretation.
2. Re-run the failing LU/evaluation/spin-sum smoke tests and then the scalar
   cross-section guardrails to ensure the ambiguity fix does not perturb the
   now-green scalar inspect matrix.
3. Address the amplitude-like divergent scalar bubble CFF local-UV behavior
   with a small first-principles comparison:

   - evaluate the original bubble graph in CFF and LTD at the same sample and
   confirm the common original-integrand reference.
   - extract the CFF 3D local UV, CFF expanded-4D local UV, and LTD expanded-4D
   local UV limits term by term before event assembly.
   - classify any mismatch as one of: source-basis projection sign, repeated or
   duplicate-channel collapse, integrated/local UV normalization, or numerator
   localization.
   - fix only the shared source-basis construction if the mismatch is in CFF
   local UV generation. Do not change the cross-section LU residue bridge unless
   the bubble reduction exposes the same graph/LMB invariant.
4. Re-run `cargo fmt`, `cargo check`, `just test_gammaloop`, and the slow
   scalar cross-section sweep after each principled fix. The repeated-channel
   descriptor is shared and should remain a single implementation rather than a
   special-case branch.
