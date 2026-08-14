#import "../../shared.typ": boundary, catalog-contract, source-link

#let api = [
= Rust and Python API boundaries

#catalog-contract(
  rust-scope: "vakint",
  python-scope: "symbolica.community.vakint",
)

== Rust package

The main Rust types make the workflow state explicit:

- `Vakint` owns the generated topology library and matching/evaluation operations;
- `VakintSettings` owns symbols, precision, normalization, tool paths, validation, and backend
  order;
- `VakintExpression` and `VakintTerm` separate numerators from topology atoms;
- `EvaluationOrder` and `EvaluationMethod` select AlphaLoop, MATAD, FMFT, or pySecDec;
- `NumericalEvaluationResult` stores Laurent coefficients in the dimensional regulator.

For a dependency-free canonicalization setup, select an empty backend order before validating:

```rust
use vakint::{EvaluationOrder, Vakint, VakintSettings};

let engine = Vakint::new()?;
let settings = VakintSettings {
    evaluation_order: EvaluationOrder::empty(),
    ..VakintSettings::default()
};
engine.validate_settings(&settings)?;
```

Parsing macros and helpers build Vakint-namespaced Symbolica atoms. The generated catalog owns
their exact grammar and error types. Do not interpolate untrusted input into a macro example;
use the runtime parser path when the expression arrives dynamically.

== Python community module

#boundary("Community module and feature pairing", [
  Python users import `symbolica.community.vakint`. It is not a standalone Vakint wheel. The
  Rust build must include `symbolica_community_module`; `python_stubgen` generates type metadata
  but does not, by itself, enable the community module in the current feature graph.
])

Install a Symbolica Python distribution only when its release includes Vakint and verify it with
`python -c "import symbolica.community.vakint"`; Vakint is not installed by a separate wheel in
this repository. Source embedders add the crate to the external
#link("https://github.com/benruijl/symbolica-community")[symbolica-community] assembly, enable
`symbolica_community_module`, and invoke `VakintWrapper`'s `SymbolicaCommunityModule`
registration while building the Symbolica extension. Stub generation describes that mounted
surface but never mounts it into an existing Python installation.

The module exports four classes: `Vakint`, `VakintExpression`, `VakintEvaluationMethod`, and
`VakintNumericalResult`. This example constructs a matching-only engine without probing external
backends:

```python
from symbolica.community.vakint import Vakint

vakint = Vakint(evaluation_order=[])
# Supply a Symbolica Expression in the Vakint namespace to `to_canonical`.
canonical = vakint.to_canonical(expr, short_form=True)
```

`Vakint` exposes canonicalization, tensor reduction, integral evaluation, complete evaluation,
and numerical conversion. `VakintEvaluationMethod` constructs configured backend choices.
`VakintNumericalResult` converts Laurent coefficients to Python tuples and compares results with
thresholds and optional errors.

Construct one engine and reuse it. When enabling pySecDec, give explicit numerical masses and
external momenta and decide whether generated output may be reused. Do not run those heavy,
external-tool examples as part of an ordinary documentation render.

Source starting points are #source-link("crates/vakint/src/lib.rs", label: "the Rust API") and
#source-link("crates/vakint/src/symbolica_community_module/mod.rs", label: "the Python module").
]
