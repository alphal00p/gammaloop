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
`CrossSectionGraph::residue_selector_for_raised_cut_group`. The LU residue
plan separates:

- selected generated LTD `E`-surface denominator signs;
- local-series denominator signs used for Laurent extraction;
- CFF/LTD surface-family projection bridges;
- positive-energy Cutkosky edge-flow orientation signs;
- repeated-channel source-basis bridges;
- direct-original prefactors, which are distinct from expanded-source residue
  prefactors.

Two sign conclusions are now encoded in production code:

- Expanded-4D LTD UV-leading local sources use the same repeated-channel
  source-basis bridge as finite cograph residues. The bridge is a property of
  the reduced source basis, and UV rescaling can collapse distinct original
  denominators onto the same equal-energy channel.
- Ordinary simple direct-original LU cuts use their own resolved
  Laurent-coordinate Jacobian. An unrelated repeated channel elsewhere in the
  graph belongs to that repeated residue, not to this simple local series. The
  repeated-channel bridge remains active for repeated residues and for explicit
  simple cuts whose support spans repeated-channel supports.

This is still a graph/LMB-derived construction. There are no GL-specific
branches, no loop-number fudge factors, and no temporary
`GAMMALOOP_TRACE_LTD_LU_SIGNS` diagnostics in the touched Rust code.

## Validation Matrix

Commands below were run with `INSTA_FORCE_PASS=1` so snapshot drift did not
hide numerical parity failures. Generated `.snap.new` files were not accepted
or staged at this checkpoint.

Passing guardrail after the current direct-original sign cleanup:

```text
GL03 no numerator, GL03 q1, GL03 q7
GL12 no numerator, GL12 q1, GL12 q7
GL07, GL08, GL21 q1, GL23 no numerator, GL23 q1, GL23 q7
GL24 no numerator, GL24 q1, GL24 q7
GL35 no numerator, GL35 q1
GL47

18/18 passed
```

Interrupted broad slow scalar cross-section sweep:

```text
122/143 tests run before manual cancellation
117 passed
4 real failures
1 SIGTERM from cancellation
```

Real failures observed in that sweep:

```text
quadratic_energy_numerators::GL38 q1_squared
quadratic_energy_numerators::GL38 q7_squared
quadratic_energy_numerators::GL46 q1_squared
quadratic_energy_numerators::GL46 q7_squared
```

The GL46 q1 failure was inspected in the run output. It is a pure relative sign
mismatch in `Original` for `LTD local-4D vs CFF local-3D`:

```text
group 0 event 1, cut_id 1
FullMultiplicativeFactor matches
ThresholdCounterterm { subset_index: 0 } is zero in both modes
Original actual   = +2.6279245830775977e-9
Original expected = -2.6279245830775977e-9
event weight differs by the same sign
```

GL24 and GL35, which were the previous non-sign local-UV drift failures, are
now green in the current guardrail. The remaining known failures are sign-only
quadratic numerator cases in the GL38/GL46 family.

## Theory Anchor

`docs/3dreps/generalised_ltd.tex` states that local unitarity cut selection and
threshold counterterm construction should use a representation-neutral
canonical causal `E`-surface catalogue. CFF exposes this catalogue directly.
LTD may contain additional `H`-surfaces, but those are not physical threshold
surfaces after the orientation sum.

The implementation should identify LU cuts and threshold residues once from
graph/LMB/linnet and canonical `E`-surface data, then map the selected residue
variables back to the requested representation. Simple, raised, confluent,
threshold, and expanded-4D UV sources should not use separate ad-hoc sign
algorithms.

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
- `uv/approx/mod.rs` and `uv/approx/expanded_4d.rs` dispatch threshold residue
  selection by representation and must remain LTD-expression-facing for LTD
  output.

## Next Debugging Plan

1. Re-run a narrow GL38/GL46 diagnostic set with sign tracing restored only
   locally and remove the trace before committing. Capture, for each failing
   event, the residue selector signs, selected denominator signs, Cutkosky
   orientation signs, local-series prefactor, direct-original prefactor, and
   source-basis bridge.
2. Classify GL38/GL46 by the same graph/LMB invariants used by the green
   guardrails: simple vs repeated cut, mixed repeated-channel support, raised
   support normalization, threshold enabled state, and numerator support. Do
   not introduce graph-name conditionals.
3. Compare GL38/GL46 against the green neighbors GL37/GL40/GL45/GL47. The
   question is whether the remaining sign is attached to the numerator-induced
   expanded-4D source basis, the direct-original simple local-series Jacobian,
   or a Cutkosky orientation convention for a specific support type.
4. Accept a change only if it can be phrased as a single construction in the
   centralized residue plan. After any change, re-run the 18-test guardrail
   above first, then the GL38/GL46 focused set, then the full slow scalar sweep.
5. Before final handoff, run `cargo fmt`, `cargo check`, `just test_gammaloop`,
   and the full slow scalar cross-section sweep. Do not stage generated
   `.snap.new` files unless they are deliberately accepted as references after
   the implementation is fixed.
