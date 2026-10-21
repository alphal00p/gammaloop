#import "../../shared.typ": boundary, catalog-contract, source-link

#let api = [
= Rust and Python APIs

#catalog-contract(
  rust-scope: "vakint",
  python-scope: "symbolica.community.hepkit.vakint",
)

== Rust package

The main Rust types make the workflow state explicit:

- `Vakint` owns the topology library and matching/evaluation operations;
- `VakintSettings` owns symbols, precision, normalization, tool paths, validation, and backend
  order;
- `VakintExpression` and `VakintTerm` separate numerators from topology atoms;
- `EvaluationOrder` and `EvaluationMethod` select RustRed, AlphaLoop, MATAD, FMFT, or pySecDec;
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

== External RustRed data

`Vakint::load_official_external_data(max_loop_count, &options)` installs and loads
all loop packs from one through the requested maximum. `OfficialExternalDataOptions`
defaults to version `latest`, the current working directory, and no refresh.
The matching Python method accepts `max_loop_count`, `version`,
`local_artifacts_path`, and `refresh`. A path containing a checked-out external-data
manifest is validated and loaded without network access or checkout modification.
Installed releases are reused offline; explicit refresh checks for updates.

`RustRedMasterMode` chooses `Retain`, `Symbolic`, or `Numerical` (default). Python
selects the same behavior with `new_rustred_method(master_mode="symbolic")`.
Symbolic mode leaves numerical-only Laurent coefficients as named symbols, which
the same engine can resolve later. Numerical mode reports finite-source precision limits,
and rejects missing epsilon coefficients or unresolved masters with their exact identities.
All three RustRed modes preserve inexact input coefficients for later tracked evaluation.
Passing retained uncertainty sources to pySecDec returns an explicit error before input generation.
Use #source-link("crates/vakint/external_data.typ", label: "the external-data guide")
for campaign coverage, cache bindings, and release conventions.

The #link("reference/rust/")[Rust orientation] identifies the compiled crates, versions, and
feature matrix and links directly to native Rustdoc. The generated
#link("reference/topologies/")[topology reference] records the runtime topology library and
supported external-tool versions for this release.

== Python community module

#boundary("Python availability", [
  Python users import `symbolica.community.hepkit.vakint`; Vakint is not distributed as a standalone
  Python package. Check that the installed Symbolica distribution includes the Vakint community
  module before relying on this interface.
])

Install and verify the published assembly with
`python -m pip install --upgrade symbolica` and
`python -c "import symbolica.community.hepkit.vakint"`. Custom Symbolica builds can add Vakint to the
#link("https://github.com/symbolica-dev/symbolica-community")[symbolica-community] assembly, enable the
`symbolica_community_module` feature, and register `VakintWrapper` while building the extension.

The structured #link("reference/python/")[Python reference] covers the engine, expression,
evaluation-method and numerical-result classes, plus typed constant/master identities,
inspection records, precision metadata and warnings. This example constructs a matching-only
engine without probing external backends:

```python
from symbolica import E
from symbolica.community.hepkit.vakint import Vakint

vakint = Vakint(evaluation_order=[])
expr = E(
    "topo(prop(18,edge(7,7),k(99),muvsq,1))",
    default_namespace="vakint",
)
canonical = vakint.to_canonical(expr, short_form=True)
assert "I1L" in str(canonical)
```

`Vakint` exposes canonicalization, tensor reduction, integral evaluation, complete evaluation,
and numerical conversion. `VakintEvaluationMethod` constructs configured backend choices.
`VakintNumericalResult.coefficients` maps epsilon powers to read-only `NumericalCoefficient`
objects. Their `numerical_value` is a native Symbolica `ComplexFloat`, and `precision.real`
and `precision.imaginary` expose achieved significant digits, absolute error, exact-zero status
and error quality. Absolute errors and real constant/master-coefficient values use Symbolica's
existing arbitrary-precision `Float` class. `requested_decimal_precision` and `warnings` belong
to the complete result.
`numerical_evaluation` continues to return `(central_result, optional_statistical_error)`.

Parameter values and external four-vector components accept native `Float` and `ComplexFloat`
objects without changing their precision. Python `float` and `complex` inputs carry their
binary64 precision, even when more working digits are requested. Python integers and numerical
Symbolica expressions are substituted symbolically before tracked evaluation, so exact rational
masses and exact cancellations in momentum contractions are preserved. Supply high-precision
data as, for example, `Float("1.234567890123456789", decimal_digits=200)`, rather than first
converting it to a Python float.
Floating coefficients embedded in an `Expression` retain separate source identities for
different parameters. Equal central values therefore do not turn independent uncertain inputs
into a certified exact cancellation.
Parameter mappings cannot override registered constants or reserved data-source identities;
qualified ordinary parameter names remain supported.

`to_expression()` preserves independently achieved component precision and raises an error for
uncertain zeros or values without reliably encodable digits. `to_list()` intentionally rounds
to Python floats and drops accuracy metadata. `compare_to()` retains its relative threshold
and optional statistical-error arguments.

```python
from symbolica.community.hepkit.vakint import Vakint, ConstantSymbol

engine = Vakint(evaluation_order=[])
constant = engine.get_constant(ConstantSymbol("z3"), prec=40)
print(constant.available_precision.kind)
print(constant.precision.real.significant_digits)

master = engine.get_master(engine.master_symbols()[0], prec=32)
for power, coefficient in master.coefficients.items():
    print(power, coefficient.numerical_value, coefficient.symbolic_expression.kind)
```

`get_constant` and `get_master` query loaded or vendored RustRed data without downloading packs.
Both accept typed identifiers and default to 16 decimal digits. Discovery returns qualified
identifiers; ambiguous friendly names report the matching alternatives. `VakintConstant`,
`VakintMaster`, `MasterCoefficient`, and the shared precision/provenance records have structured
terminal representations and `to_html()` / `_repr_html_()` notebook views. Native Symbolica
rendering preserves long formulas and numbers behind expandable sections, and visible warnings
explain precision shortfalls.
Loaded finite constants such as `K5Moment(1)` report their validated source cap.
Master formulas retain their original snapshot internally; notebook and terminal views show the
mathematical names, while the expression getter preserves the binding needed for evaluation.
Symbolic reductions preserve the same finite-period bindings when their formulas are evaluated
later; loading another data snapshot does not silently replace those sources.

`master.graph` is the normal FeynKit `FeynmanDiagram`, with its existing DOT, SVG and HTML
interfaces. In Marimo, display the master directly or use `mo.Html(master.to_html())`.
Graph previews require FeynKit's optional `typst-py` renderer; metadata and coefficient displays
remain available without it. Rendering does not reevaluate numbers or change the data snapshot.

Construct one engine and reuse it. When enabling pySecDec, provide numerical masses and external
momenta where the topology requires them. Reuse existing pySecDec output only while debugging,
and remove it before evaluating changed inputs.

Source starting points are #source-link("crates/vakint/src/lib.rs", label: "the Rust API") and
#source-link("crates/vakint/src/symbolica_community_module/mod.rs", label: "the Python module").
]
