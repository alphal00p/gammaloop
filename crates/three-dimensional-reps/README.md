# three-dimensional-reps

`three-dimensional-reps` builds generalized CFF loop-energy representations
used by GammaLoop. It owns the expression data model that was previously local
to `gammalooprs::cff`.

## Library Path

The GammaLoop-facing API is intentionally small:

- `generate_3d_expression(...)`
- `Generate3DExpressionOptions`
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
`3Drep evaluate` command reuses GammaLoop's evaluator construction machinery to
materialize the Symbolica expression and parameter builder artifacts.

## Current Limits

The initial production boundary covers affine CFF, bounded-energy CFF
completions for repeated and multiloop high-power cases, and uniform numerator
sampling-scale modes. Four-dimensional UV support and a public LTD engine are
staged separately so each layer has one clear owner.
