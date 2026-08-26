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

- `serde` (default): retained as downstream compatibility vocabulary. Serde
  dependencies and derives are unconditional, so disabling this no-op feature
  does not change the crate's serialization API.
- `display`: colored/table expression summaries.
- `eval`: optional eager f64 diagnostic evaluator and oracle tests.

The compiled Symbolica runtime evaluator from the Python prototype is not part
of this library crate. Neither GammaLoop nor the normal workspace test suite
enables `eval`; run its diagnostic inventory manually and serially with:

```bash
cargo nextest run -p three-dimensional-reps --features eval \
  --no-fail-fast --test-threads 1
```

The GammaLoop CLI-side `3Drep build` command is diagnostic-only: it validates,
renders, and optionally writes the oriented expression, but its input and
expression preparation are not a GammaLoop production contract. Production
evaluator construction remains in GammaLoop.

## Current Limits

The production boundary covers affine CFF, bounded-energy CFF with
repeated-channel normal form, multiloop high-power numerators, and uniform
numerator sampling-scale modes. GammaLoop owns both UV orchestration orders:
direct local 3D subtraction on CFF expressions and 4D-local UV followed by an
explicit-orientation CFF sum. `RepresentationMode::Ltd` is reserved for a
future proper LTD backend and currently returns an explicit not-implemented
error.
