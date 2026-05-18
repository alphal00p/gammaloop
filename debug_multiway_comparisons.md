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

The current implementation keeps LU residue metadata centralized in
`CrossSectionGraph::residue_selector_for_raised_cut_group`. The plan now stores
separate signs for:

- the selected generated LTD `E`-surface denominator;
- the LU local-series denominator used for Laurent extraction;
- the CFF/LTD surface-family projection bridge;
- the positive-energy Cutkosky edge-flow orientation;
- the direct original-integrand prefactor, which is distinct from the ordinary
  multi-residue prefactor used by expanded-source residues.

The direct original-integrand path is now stable on the previously oscillating
simple/repeated guardrails. Simple mixed contacts with repeated channel support
use the generated selected-surface orientation and Cutkosky edge flow. Simple
cuts in graphs with repeated LU groups use the full-graph projection bridge and
the generated-LTD repeated-channel routing bridge, because the direct
`Original` extraction has not selected the repeated-channel residue yet.
Repeated/confluent cuts use the repeated-channel local Laurent wedge parity.
This is still a single graph/LMB-derived construction; there are no GL-specific
branches.

Threshold residues are now representation-aware:

- CFF threshold selection uses canonical CFF surface-family residue selection.
- LTD threshold selection consumes the generated LTD threshold-surface
  denominator directly through `select_esurface_residue_in_generated_basis`.
- Cross-section threshold `ResidueSelector` prefactors are neutral for
  left/right threshold counterterms because the LU threshold evaluator supplies
  the subtractive Cauchy orientation after threshold localization.

This checkpoint is useful because it fixes the GL03/GL12 pure sign family that
remained after the previous 132/143 full-sweep checkpoint, without temporary
`GAMMALOOP_TRACE_LTD_LU_SIGNS` diagnostics. It is not the final parity state:
the remaining known failures are the GL24/GL35 non-sign local-UV drifts.

## Validation Snapshot

Commands were run with `INSTA_FORCE_PASS=1` so generated snapshot drift did not
hide numerical parity failures.

Passing at the previous full-sweep checkpoint:

- Focused threshold cluster: GL28, GL29, GL32, GL38, GL39, GL41.
- Focused direct-original cluster: GL08, GL21, GL37, GL38, GL44, GL45, GL46,
  GL47.
- Original simple/repeated LU guardrail: GL07, GL47, GL06 q1, GL21 q1, GL27
  q1, GL40.
- Full slow scalar cross-section sweep before the latest GL03/GL12 fix:

```text
143 tests run: 132 passed, 11 failed, 124 skipped
```

Additional targeted validation after the latest direct-original prefactor fix:

```text
GL03 no numerator, GL03 q1, GL03 q7,
GL12 no numerator, GL12 q1, GL12 q7,
GL07, GL08, GL21 q1, GL47: 10/10 passed
```

Known remaining slow sweep failures:

```text
GL24 q1, q7, no_numerator
GL35 q1, no_numerator
```

The failures split into two clean categories:

- GL03 and GL12 were pure sign mismatches in `Original`. GL12 threshold
  counterterms followed the same overall sign as `Original`, so the immediate
  issue was not threshold selection. They are now covered by the targeted
  simple/repeated LU guardrail.
- GL24 and GL35 are small non-sign `Original` drifts with matching event/cut
  metadata and matching full multiplicative factors. Threshold subtraction is
  disabled in these tests; this points at the expanded-4D local UV/original
  extraction path rather than threshold normalization.

Generated `.snap.new` files were produced by forced-pass runs and should not be
staged unless accepted deliberately as new references.

## Theory Anchor

`docs/3dreps/generalised_ltd.tex` states that local unitarity cut selection and
threshold counterterm construction should use a representation-neutral
canonical causal `E`-surface catalogue. CFF exposes this catalogue directly. LTD
may contain additional `H`-surfaces, but those are not physical threshold
surfaces after the orientation sum.

The implementation should identify LU cuts and threshold residues once from
graph/LMB/linnet and canonical `E`-surface data, then map the selected residue
variables back to the requested representation. Simple, raised, confluent,
threshold, and expanded-4D UV sources should not use separate sign algorithms.

Edges with `is_cut` are forward-scattering cut edges. For graph manipulation
they should behave as two external half-edges with the same momentum, one
incoming and one outgoing. LU/Cutkosky and threshold support normalization must
therefore treat them as cut/external boundaries rather than ordinary internal
propagators.

## Sensitive Code Paths

- `CrossSectionGraph::residue_selector_for_raised_cut_group` builds the LU
  residue plan and is the main entry point to keep centralized.
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
- `uv/approx/mod.rs` and `uv/approx/expanded_4d.rs` now dispatch threshold
  residue selection by representation and must remain LTD-expression-facing for
  LTD output.

## Next Debugging Order

1. Keep GL03/GL12 in the guardrail set. Their fix came from the same generated
   LTD repeated-channel bridge used for repeated-channel residues, applied to
   simple direct-original cuts in graphs whose LU surface family contains an
   unselected repeated channel.
2. Re-run the guardrails after any sign or routing change: GL06 q1, GL07, GL21 q1, GL27
   q1, GL40, GL47, plus GL08/GL37/GL38/GL44/GL45/GL46.
3. Debug GL24/GL35 as non-sign drift.
   Start from event selection, selected support, representation choice,
   numerator completion, and expanded-4D local UV/original separation before
   touching any normalization.
4. Run the full slow scalar sweep after each accepted principled fix. A change
   is not accepted if it reopens an already green guardrail.
5. Before final handoff, run `cargo fmt`, `cargo check`, `just test_gammaloop`,
   and the full slow scalar cross-section sweep. Do not stage generated
   `.snap.new` files unless they are deliberately accepted as references.
