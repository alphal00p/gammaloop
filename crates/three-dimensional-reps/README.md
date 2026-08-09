# three-dimensional-reps

`three-dimensional-reps` builds generalized CFF loop-energy representations
used by GammaLoop. It owns the expression data model that was previously local
to `gammalooprs::cff`.

## Library Path

The GammaLoop-facing API is intentionally small:

- `generate_3d_expression(...)`
- `Generate3DExpressionOptions`
- `GeneratedThreeDExpression` and its `CffEnergyFactorOwnership`
- the serializable `ThreeDExpression<OrientationID>` data model

The production GammaLoop path calls the graph-first generator through traits
implemented for GammaLoop's already-parsed graph representation. This crate does
not parse user DOT graphs directly; DOT import remains owned by the GammaLoop
CLI/state layer.

## Features

- `serde` (default): serializable expression data model.
- `display`: colored/table expression summaries.
- `eval`: eager f64 evaluator used by integration tests.

The compiled Symbolica runtime evaluator from the Python prototype is not part
of this library crate. The current `eval` feature is an eager f64 diagnostic
evaluator used by the GammaLoop integration tests. The GammaLoop CLI-side
`3Drep build` command validates, renders, and optionally writes the oriented
expression; production evaluator construction remains in GammaLoop.

## Current Limits

The production boundary covers affine CFF, bounded-energy CFF with
repeated-channel normal form, multiloop high-power numerators, and uniform
numerator sampling-scale modes. GammaLoop owns both UV orchestration orders:
direct local 3D subtraction on CFF expressions and 4D-local UV followed by an
explicit-orientation CFF sum. `RepresentationMode::Ltd` is reserved for a
future proper LTD backend and currently returns an explicit not-implemented
error.
