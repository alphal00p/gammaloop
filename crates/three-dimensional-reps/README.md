# three-dimensional-reps

`three-dimensional-reps` builds generalized CFF and LTD loop-energy representations
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
its integration-test dependency. Run its diagnostic inventory directly with:

```bash
cargo nextest run -p three-dimensional-reps --features eval \
  --no-fail-fast --retries 0
```

The GammaLoop CLI-side `3Drep build` command is diagnostic-only: it validates,
renders, and optionally writes the oriented expression, but its input and
expression preparation are not a GammaLoop production contract. Production
evaluator construction remains in GammaLoop.

## Current Limits

The production boundary covers affine CFF, bounded-energy CFF with
independent-occurrence normal form, multiloop high-power numerators, and uniform
numerator sampling-scale modes. GammaLoop owns both UV orchestration orders:
direct local 3D subtraction on CFF expressions and 4D-local UV followed by an
explicit-orientation CFF sum.

The shared generator also supports `RepresentationMode::Ltd`: signed energy
residues with affine numerator maps, explicit E/H surfaces, and term-local
on-shell energy factors. Simple poles do not need numerator-degree bounds;
repeated poles require explicit bounds for exact derivative-free polynomial
sampling (`Some([])` selects an energy-independent numerator). All residue
terms contribute to an evaluation, including the cancellation of H surfaces.
Physical-pole selection retains the factorized numerator maps and requests
polynomial bounds when a raised pole needs them.

The LTD core and its tests are ported from `ltd_support` commits `6fac4ffaf4`
and `6b35f5a2a4`. GammaLoop's production LTD orchestration and Python residue-map
bindings are separate integrations.
