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

`old_cff` exists only behind the `old_cff` feature as a migration/parity anchor
for the legacy affine CFF implementation. It rejects higher energy-numerator
powers.

## Features

- `serde` (default): serializable expression data model.
- `display`: colored/table expression summaries.
- `diagnostics`: exposes diagnostic-only modes such as `pure_ltd`.
- `eval`: eager f64 evaluator used by CLI diagnostics.
- `old_cff`: exposes the legacy CFF comparison mode.
- `test-support`: exposes internal diagnostic modes such as `pure_ltd`.

The compiled Symbolica runtime evaluator from the Python prototype is not part
of this library crate. The current `eval` feature is an eager f64 diagnostic
evaluator used by the GammaLoop integration tests. The GammaLoop CLI-side
`3Drep evaluate` command reuses GammaLoop's evaluator construction machinery to
materialize the Symbolica expression and parameter builder artifacts.

## CLI

Diagnostic commands are integrated in `gammaloop-api` under the `3Drep`
command. They operate on graphs already present in the active GammaLoop state,
except for the orthogonal `3Drep graph-from-signatures` helper, which emits a
GammaLoop-format DOT graph from a propagator-signature expression.
`3Drep evaluate` can either load an explicit oriented-expression JSON with
`--json-in` or resolve a cached build artifact from
`-p/-i/-g --representation ...`, with the old latest-artifact pointer kept only
as an ad-hoc fallback.

## Current Limits

The expression-generation and eager-evaluation coverage mirrors the validated
old Python prototype surface currently ported into GammaLoop tests: affine CFF,
LTD with repeated propagators, bounded-energy CFF completions for the repeated
and multiloop high-power cases in the old matrix, uniform numerator sampling
scale modes, graph-from-signatures reconstruction, and the slow five-loop
ultimate-basis diagnostic. `pure_ltd` remains diagnostic-only for repeated
propagators because derivative-free pure LTD is not a valid numerical equality
test in that case.
