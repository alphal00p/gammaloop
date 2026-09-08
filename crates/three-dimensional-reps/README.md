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
of this library crate. The normal workspace test suite enables `eval` through
its integration-test dependency. Run its diagnostic inventory directly and
serially with:

```bash
cargo nextest run -p three-dimensional-reps --features eval \
  --no-fail-fast --test-threads 1
```

Thermal generation uses `Generate3DExpressionOptions::medium_mode`. Finite-
temperature and zero-temperature equilibrium use the same structural recursion,
including cyclic orientations and distribution derivatives. The serializable
`ThermalWeight` on each variant retains symbolic edge distributions separately
from rational CFF coefficients. GammaLoop owns their expansion using particle
statistics, chemical potentials, and temperature. The standalone `eval` feature
does not accept these model-dependent distribution inputs and rejects thermal
expressions explicitly. Thermal generation retains on-shell numerator maps and
does not support uniform numerator sampling scales.

The GammaLoop CLI-side `3Drep build` command is diagnostic-only: it validates,
renders, and optionally writes the oriented expression. It forwards the configured
medium to the shared generator and records it in the JSON output, but its input and
expression preparation are not a GammaLoop production contract. Production
evaluator construction remains in GammaLoop.

## Current Limits

The vacuum production boundary covers affine CFF, bounded-energy CFF with
repeated-channel normal form, multiloop high-power numerators, and uniform
numerator sampling-scale modes. GammaLoop owns both UV orchestration orders
in vacuum: direct local 3D subtraction on CFF expressions and 4D-local UV
followed by an explicit-orientation CFF sum. Medium modes and optional vacuum
subtraction support only direct local 3D subtraction; the local-4D option is
rejected. `RepresentationMode::Ltd` is reserved for a
future proper LTD backend and currently returns an explicit not-implemented
error.
