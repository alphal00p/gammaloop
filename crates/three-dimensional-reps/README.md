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

Serde serialization is always enabled.

- `display`: colored/table expression summaries.
- `eval`: optional eager f64 diagnostic evaluator and oracle tests.
- `json_schema`: JSON schema support for the representation and medium enums.
- `cli`: representation values for the existing GammaLoop CLI.
- `python`: Python bindings for the shared representation and medium enums.

The compiled Symbolica runtime evaluator from the Python prototype is not part
of this library crate. The normal workspace test suite enables `eval` through
its integration-test dependency. Run its diagnostic inventory directly with:

```bash
cargo nextest run -p three-dimensional-reps -p gammaloop-workspace-hack \
  --features three-dimensional-reps/eval \
  --no-fail-fast --retries 0
```

Including the workspace feature-unification package selects Symbolica's integer
and floating-point backends, as in GammaLoop builds.

Thermal CFF and LTD generation use `Generate3DExpressionOptions::medium_mode`
and support uniform numerator sampling scales in both equilibrium modes. GammaLoop supplies
the model-dependent distributions; the standalone `eval` feature rejects thermal
expressions because it does not accept those inputs. See the
[CFF architecture](../../docs/architecture/architecture-current.typ#cff-production-and-numerator-energy-ownership)
for thermal weights, numerator reduction, and ownership boundaries.

The GammaLoop CLI-side `3Drep build` command is diagnostic-only: it validates,
renders, and optionally writes the oriented expression. It forwards the configured
medium to the shared generator and records it in the JSON output, but its input and
expression preparation are not a GammaLoop production contract. Production
evaluator construction remains in GammaLoop.
Select `--representation ltd` to display affine residue keys and explicit E/H
surfaces, and `--show-details-for-residue <id>` to inspect a complete map row.
There is no separate CLI evaluation subcommand.

## Current Limits

The vacuum production boundary covers affine CFF, bounded-energy CFF with
independent-occurrence normal form, multiloop high-power numerators, and uniform
numerator sampling-scale modes. GammaLoop owns both UV orchestration orders:
direct local 3D subtraction on CFF expressions and 4D-local UV followed by a
complete CFF or LTD residue sum. Medium modes and optional vacuum subtraction
support only CFF with direct local 3D subtraction; production and diagnostic CLI
LTD still reject both equilibrium modes. The shared LTD generator supports them
by reflecting each vacuum residue's cut poles and attaching thermal weights,
including analytic distribution derivatives for raised poles.
`RepresentationMode::Ltd` retains the signed energy-residue sum and its
explicit E/H surfaces. Simple LTD poles do not need numerator-degree
bounds; repeated poles require explicit bounds for exact derivative-free
polynomial sampling (`Some([])` selects a scalar numerator).
LTD residue terms must be summed together for their H-surface cancellation.
Physical-pole localization uses the same factorized affine maps. When
localization creates a raised pole, its Laurent coefficients request polynomial
bounds lazily; ordinary simple-pole generation remains independent of them.
