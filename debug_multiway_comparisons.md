# Multi-Way Comparison Debug Status

This note tracks the unfinished effort to make scalar cross-section rich local
inspect output agree across:

1. CFF with local UV counterterms from 3D expansions,
2. CFF with local UV counterterms from expanded 4D integrands,
3. LTD with local UV counterterms from expanded 4D integrands.

The goal is event-level parity, including `Original`, threshold counterterms,
full multiplicative factors, cut ids, event grouping, and weights. Overall
common phases are less important than relative parity between the three modes,
but no graph-specific or unexplained sign modifiers are allowed.

## Current Baseline

Commands were run with `INSTA_FORCE_PASS=1` to avoid treating snapshot drift as
the primary signal.

Passing:

- `cargo fmt`
- `cargo check`
- Focused GL06/GL21/GL27/GL40 LU sign control set.
- Dedicated GL48 generated UV/threshold repeated inspect check.
- Local UV inspect checks, including divergent scalar bubble:
  CFF expanded-4D local UV agrees with CFF local-3D, and LTD expanded-4D local
  UV agrees with CFF for these inspect fixtures.
- Repeated/mass-shift bubble limit checks:
  split-mass and equal-mass repeated bubble tests pass.

Failing:

- Focused generated no-UV LTD-vs-CFF matrix: 6 passed, 10 failed.
  Failing cases include GL00 triple repeated, GL10 simple/raised/repeated
  quartic, GL12 simple/raised, GL13 simple, and the amplitude-like bulk
  profiles.
- Slow scalar cross-section sweep:
  143 tests run, 91 passed, 52 failed, 124 skipped.
  Failures are primarily `LTD local-4D vs CFF local-3D`; many are `Original`,
  while at least GL28 also flips a threshold counterterm and GL48 shows a small
  non-sign drift.
- Amplitude bulk profiles:
  dotted bubble, double-dotted bubble, and scalar self-energy all fail their IR
  profile checks. These are not the direct local-UV inspect parity failures; the
  direct divergent bubble inspect comparison currently passes.

Do not stage generated `.snap.new` files unless the implementation is fixed and
the snapshots are deliberately accepted as new references.

## Theory Anchor

`docs/3dreps/generalised_ltd.tex` states that local unitarity cut selection and
threshold counterterm construction should use a representation-neutral
canonical causal `E`-surface catalogue. CFF exposes this catalogue directly. LTD
may contain additional `H`-surfaces, but those are not physical threshold
surfaces after the orientation sum.

The implementation should therefore identify LU cuts and threshold residues once
from graph/LMB/linnet and canonical `E`-surface data, then map the selected
residue variables back to the requested representation. Simple, raised,
confluent, threshold, and expanded-4D UV sources should not use separate sign
algorithms.

## Planned Central Abstraction

The scattered pieces should be centralized in a single residue-coordinate plan
owned by `ResidueSelector`:

```text
ResiduePlan
  LU cut group, if any
  left/right threshold groups, if any
  normalized Cutkosky edge support
  generated E-surface ids
  canonical/projection E-surface ids when needed
  localized external-energy coordinate
  positive-energy Cutkosky orientation
  determinant/Jacobian sign of the resolved LMB coordinate map
  selected-variable signs for each basis
  ordinary-residue and direct-local-series prefactors as derived accessors
```

Callers should request one of a small set of modes:

- select in canonical/CFF `E`-surface basis,
- select in generated LTD `E`-surface basis,
- select in positive-energy Cutkosky basis,
- extract direct LTD LU Laurent coefficients,
- select threshold residues and then localize threshold variables.

No raising is `max_occurrence = 1`; raised/confluent cuts use the same plan with
higher multiplicity. Threshold E-surfaces use the same plan with
`ResidueKind::Threshold`, not a threshold-specific sign bridge.

## Validation Ladder

Work should proceed in layers and lock each layer before moving on:

1. Behavior-preserving refactor into `ResiduePlan`.
   The focused pass/fail set should not change.
2. Original generated no-UV parity, simple cuts first:
   GL10/GL12/GL13 simple cuts. Magnitude mismatches here indicate residue
   support or coordinate selection bugs, not UV issues.
3. Original generated no-UV parity for raised/confluent cuts:
   GL00 triple repeated, GL10/GL12 raised, GL10 quartic energy-trade.
   Repeated/mass-shift tests must stay green.
4. Threshold normalization in isolation:
   compare threshold weights separately from `Original`, especially GL28.
5. CFF local UV:
   keep CFF expanded-4D vs CFF local-3D green, including divergent bubbles.
6. LTD expanded-4D UV:
   only after original residues and threshold selection are stable.
7. Full validation:
   focused scalar controls, slow scalar cross-section sweep, `cargo fmt`,
   `cargo check`, and `just test_gammaloop`.

Temporary sign tracing such as `GAMMALOOP_TRACE_LTD_LU_SIGNS` must be removed
before the work is considered commit-ready.
