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

The latest implementation centralizes the LU sign pieces for each raised LU
group in one plan component before constructing `ResidueSelector`. The direct
LTD original-integrand prefactor is now stored explicitly in the LU residue
plan, separate from the local-series extraction prefactor and the ordinary
multi-residue prefactor.

The simple-cut oscillation between GL07, GL21, and GL47 was resolved by making
the simple-cut bridge depend on normalized edge-support contact with repeated
channels:

- a simple LU cut that touches multiple distinct repeated channel supports
  through multiple distinct cut edges uses the generated selected `E`-surface
  orientation;
- an ordinary simple LU cut uses the local-series coordinate orientation;
- raised/confluent LU cuts use the repeated-channel residue bridge.

This is still an intermediate checkpoint. It is not yet commit-ready for the
final parity goal because the full slow scalar sweep still has failures, but the
change is useful because it removes the previous sign regression loop.

Temporary `GAMMALOOP_TRACE_LTD_LU_SIGNS` diagnostics were removed from the
checkpointed code.

## Validation Snapshot

Commands were run with `INSTA_FORCE_PASS=1` so generated snapshot drift did not
hide numerical parity failures.

Passing:

- `cargo check -p gammalooprs --profile dev-optim`
- Focused simple/repeated LU guardrail:
  GL07 no numerator, GL47 no numerator, GL06 q1, GL21 q1, GL27 q1, GL40 no
  numerator.
- Slow scalar cross-section sweep improved from the previous 91/143 pass count
  to:

```text
143 tests run: 109 passed, 34 failed, 124 skipped
```

Slow sweep failures:

```text
GL08 q7, no_numerator
GL24 q1, q7, no_numerator
GL28 q7, no_numerator
GL29 q1, no_numerator
GL32 q7, no_numerator
GL35 q1, no_numerator
GL37 q1, q7, no_numerator
GL38 q1, q7, no_numerator
GL39 q1, q7, no_numerator
GL41 q1, q7, no_numerator
GL44 q1, q7, no_numerator
GL45 q1, q7, no_numerator
GL46 q1, q7, no_numerator
```

The failures currently split into three categories:

- Direct `Original` pure sign mismatches with matching event/cut metadata and
  matching full multiplicative factors: GL08, GL37, GL38, GL44, GL45, GL46.
- Threshold counterterm sign or threshold normalization mismatches while
  `Original` usually matches: GL28, GL29, GL32, GL39, GL41.
- Small non-sign `Original` drift with matching event/cut metadata and matching
  full multiplicative factors: GL24 and GL35.

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

Edges with `is_cut` are forward-scattering cut edges: for graph manipulation
they should behave as two external half-edges with the same momentum, one
incoming and one outgoing. LU/Cutkosky and threshold support normalization must
therefore treat them as cut/external boundaries rather than ordinary internal
propagators.

## Sensitive Code Paths

- `CrossSectionGraph::residue_selector_for_raised_cut_group`
  builds the LU residue plan and is the main entry point to keep centralized.
- `CrossSectionGraph::lu_cut_edge_sets_with_cutkosky_signs`
  maps raised LU groups to positive-energy Cutkosky edge support and signs.
- `CrossSectionGraph::simple_ltd_lu_cut_local_coordinate_signs`
  computes the generated selected `E`-surface sign, the local-series
  denominator sign, and the surface-family bridge from the resolved LMB
  coordinate determinant.
- `CrossSectionGraph::simple_lu_cut_has_mixed_repeated_channel_contact`
  identifies the support-local case where a simple cut must use generated
  selected `E`-surface orientation because it bridges multiple repeated
  channels.
- `ResidueSelector` stores selected signs, local-series signs, ordinary residue
  prefactors, direct-original prefactors, and threshold prefactors.
- `cff/expression.rs::ltd_lu_local_series_coefficients_from_parametric_atom`
  consumes the LTD LU local-series coordinates and signs.
- `uv/approx/mod.rs` and `uv/approx/expanded_4d.rs` consume the selector for
  expanded-4D UV counterterms and must remain LTD-expression-facing.

## Next Debugging Order

1. Lock the direct original-integrand signs before touching thresholds.
   Use representatives GL08, GL37, GL38, GL44, GL45, and GL46. These are clean
   pure sign mismatches with matching event metadata, so any fix should be a
   determinant/Jacobian correction in the centralized LU plan, not a
   threshold/UV change.
2. Once direct `Original` parity is stable, isolate threshold counterterm
   residues with GL28, GL29, GL32, GL39, and GL41. The likely target is the
   threshold residue prefactor path or threshold support normalization, not the
   already-matching original residue.
3. After sign issues are stable, inspect GL24 and GL35 as non-sign `Original`
   drift. Start by checking event selection, selected support, representation
   choice, numerator completion, and local UV/threshold contribution separation
   before changing any normalization.
4. Re-run the guardrail set after each principled fix, then the slow scalar
   sweep. A fix is not accepted if it reopens GL06, GL07, GL21, GL27, GL40, or
   GL47.
5. Before final handoff, run `cargo fmt`, `cargo check`, `just test_gammaloop`,
   and the full slow scalar cross-section sweep without committing generated
   `.snap.new` files unless they are explicitly accepted.
