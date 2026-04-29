# three-dimensional-reps

`three-dimensional-reps` builds generalized three-dimensional loop-energy
representations used by GammaLoop. It owns the expression data model that was
previously local to `gammalooprs::cff`. Diagnostic tooling is integrated in the
GammaLoop CLI through `gammaloop-api`.

## Library Path

The GammaLoop-facing API is intentionally small:

- `generate_3d_expression(...)`
- `Generate3DExpressionOptions`
- `RepresentationMode::{Ltd, Cff}`
- the serializable `ThreeDExpression<OrientationID>` data model

The production GammaLoop path calls the graph-first generator through traits
implemented for GammaLoop's already-parsed graph representation. This crate does
not parse user DOT graphs directly; DOT import remains owned by the GammaLoop
CLI/state layer.

`old_cff` exists only behind the `old-cff` feature as a migration/parity anchor
for the legacy affine CFF implementation. It rejects higher energy-numerator
powers.

## Features

- `serde` (default): serializable expression data model.
- `display`: colored/table expression summaries.
- `diagnostics`: exposes diagnostic-only modes such as `pure_ltd`.
- `eval`: eager f64 evaluator used by CLI diagnostics.
- `old-cff`: exposes the legacy CFF comparison mode.
- `test-support`: exposes internal diagnostic modes such as `pure_ltd`.

The compiled Symbolica runtime evaluator from the Python prototype is not
ported yet. The current `eval` feature is an eager f64 diagnostic evaluator, not
GammaLoop's production evaluator path.

## CLI

Diagnostic commands are integrated in `gammaloop-api` under the `3Drep`
command. They operate on graphs already present in the active GammaLoop state,
except for the orthogonal `3Drep graph-from-signatures` helper, which emits a
GammaLoop-format DOT graph from a propagator-signature expression.

## Current Limits

The generalized CFF port covers the validated classes that are currently wired
into tests: affine CFF, one-loop bounded higher powers, repeated one-loop
channels, selected repeated spectator cases, and a conservative multiloop
single-quadratic lower-sector completion. Broader multiloop high-power CFF
completion is intentionally rejected until the remaining product/lower-sector
assembly cases are ported and validated.
