#import "../../shared.typ": boundary, catalog-contract, source-link

#let api = [
= Rust and Python APIs

#catalog-contract(
  rust-scope: "vakint",
  python-scope: "symbolica.community.vakint",
)

== Rust package

The main Rust types make the workflow state explicit:

- `Vakint` owns the topology library and matching/evaluation operations;
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

`vakint_parse!` and the related helpers parse Symbolica expressions in Vakint's namespace. Handle
parse failures before passing a user-supplied expression to matching or evaluation.

The #link("reference/rust/")[curated Rust API] presents the supported entry points before the
full Rustdoc, while the generated #link("reference/topologies/")[topology reference] records the
runtime topology library and supported external-tool versions for this release.

== Python community module

#boundary("Python availability", [
  Python users import `symbolica.community.vakint`; Vakint is not distributed as a standalone
  Python package. Check that the installed Symbolica distribution includes the Vakint community
  module before relying on this interface.
])

Verify availability with `python -c "import symbolica.community.vakint"`. Custom Symbolica builds
can add Vakint to the
#link("https://github.com/benruijl/symbolica-community")[symbolica-community] assembly, enable the
`symbolica_community_module` feature, and register `VakintWrapper` while building the extension.

The structured #link("reference/python/")[Python reference] covers the four exported classes:
`Vakint`, `VakintExpression`, `VakintEvaluationMethod`, and
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

Construct one engine and reuse it. When enabling pySecDec, provide numerical masses and external
momenta where the topology requires them. Reuse existing pySecDec output only while debugging,
and remove it before evaluating changed inputs.

Source starting points are #source-link("crates/vakint/src/lib.rs", label: "the Rust API") and
#source-link("crates/vakint/src/symbolica_community_module/mod.rs", label: "the Python module").
]
