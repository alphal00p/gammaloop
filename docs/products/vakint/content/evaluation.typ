#import "../../shared.typ": callout, boundary, source-link

#let evaluation = [
= Topology matching, reduction, and evaluation

Vakint recognizes an integral only after its propagators, masses, powers, and momentum routing
can be mapped to a supported topology pattern. The topology reference lists the patterns available
in the selected version.

== Parsing and topology matching

A `VakintExpression` is a sum of numerator/topology pairs. Parsing checks the Vakint namespace;
matching then searches allowed substitutions and momentum shifts. Unknown-topology handling is
a setting: preserving an unknown term is useful for partial simplification, but does not make
that term evaluable.

Canonical long form is explicit and suitable for auditing. Canonical short form is convenient
for display and library lookup. Both represent the same matched topology only after
canonicalization has succeeded.

== Normalization and canonicalization

Normalization conventions cover loop-measure factors, mass scales, epsilon powers, and the
requested expansion depth. Set them before comparing results from different backends. Momentum
canonicalization may rename loop variables and reorder propagators; compare canonical
expressions rather than input spelling when testing equivalence.

#callout("Matching does not require an evaluation backend", [
  Parsing, matching, normalization, and canonical routing work without FORM or pySecDec when the
  evaluation order is empty. Use this configuration when you only need a canonical topology.
])

== Tensor reduction

Tensor numerators are reduced to scalar integrals before analytic evaluation where required.
Reduction depends on the topology, Lorentz rank, dimension convention, and scalar-product
normalization. Tensor reduction requires FORM. When diagnosing a mismatch, keep the FORM input
and Vakint temporary directory so that the failing reduction can be inspected.

== Evaluation order and backends

The evaluation order is an ordered selection policy, not a request to combine backend results.
Vakint skips methods whose capability declaration does not match the canonical topology, then
uses the first matching method. An error from that selected backend is returned immediately; it
does not trigger fallback to the next method. AlphaLoop, MATAD, and FMFT are analytic FORM-backed
paths. pySecDec is numerical and additionally needs numeric parameters and external momenta where
the topology requires them.

For reproducible comparisons record:

- the canonical topology and normalization convention;
- epsilon expansion depth and decimal precision;
- the ordered backend list and backend-specific settings;
- FORM, pySecDec, MATAD, and FMFT versions where applicable;
- parameter substitutions and any retained temporary-output directory.

== Configure one analytic tensor reduction

This program makes normalization, precision, and backend order explicit before reducing a
rank-two one-loop numerator. It compiles without running external tools in the documentation
harness; running it requires a supported FORM installation for the AlphaLoop path.

// docs-example: compile vakint-backend-policy
```rust
use std::collections::HashMap;

use vakint::{
    EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings, vakint_parse,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let settings = VakintSettings {
        allow_unknown_integrals: false,
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        number_of_terms_in_epsilon_expansion: 2,
        run_time_decimal_precision: 32,
        evaluation_order: EvaluationOrder::alphaloop_only(),
        ..VakintSettings::default()
    };
    let vakint = Vakint::new()?;
    let input = vakint_parse!(
        "(k(1,1)*k(1,2)+k(1,3)*p(1,3))*topo(prop(1,edge(1,1),k(1),muvsq,1))"
    )?;
    let canonical = vakint.to_canonical(&settings, input.as_view(), true)?;
    let scalar = vakint.tensor_reduce(&settings, canonical.as_view())?;
    let evaluated = vakint.evaluate_integral(&settings, scalar.as_view())?;

    let parameter_values = HashMap::from([
        ("muvsq".to_owned(), 1.0),
        (settings.mu_r_sq_symbol.clone(), 1.0),
    ]);
    let real_parameters = vakint.params_from_f64(&settings, &parameter_values);
    let complex_parameters = HashMap::default();
    let (numerical, numerical_error) = vakint.numerical_evaluation(
        &settings,
        evaluated.as_view(),
        &real_parameters,
        &complex_parameters,
        None,
    )?;

    for (epsilon_power, coefficient) in numerical.get_epsilon_coefficients() {
        println!("epsilon^{epsilon_power}: {coefficient}");
    }
    if let Some(error) = numerical_error {
        println!("backend uncertainty: {error}");
    }
    Ok(())
}
```

The result invariant is an MS-bar Laurent series in the dimensional regulator: the reduced
expression no longer contains the original loop-momentum numerator, and only the first
successful method in `alphaloop_only()` contributes. Inspect the exact
#link("reference/rust/vakint/struct.VakintSettings.html")[`VakintSettings`],
#link("reference/rust/vakint/struct.VakintExpression.html#method.tensor_reduce")[tensor
reduction], and
#link("reference/rust/vakint/struct.VakintExpression.html#method.evaluate_integral")[evaluation]
Rustdoc before changing normalization or backend order.

== Numerical substitution and result interpretation

Backend evaluation returns a symbolic Laurent series. Numerical substitution is a separate
boundary: supply every remaining real/complex parameter and, when the numerator contains
external scalar products, the external four-momenta. `params_from_f64` and
`externals_from_f64` convert ordinary inputs while retaining at most their original binary64
precision. Increasing the evaluation precision does not add source digits to these inputs.

`NumericalEvaluationResult` stores an ordered map from epsilon powers to named
`NumericalCoefficient` records. Negative powers are poles, power zero is the finite term, and
positive powers are higher orders. Each record has a complex central value and separate real
and imaginary precision information: achieved significant decimal digits, an arbitrary-precision
absolute error, exact-zero status, and whether the error is bounded or estimated. Requested
precision and source availability are distinct from achieved accuracy.

For native Rust callers, `MATAD::substitute_masters_directly` and
`FMFT::substitute_masters` prepare finite-source identities and exact Gamma expressions.
Expand and normalize the complete Laurent expression before passing it to numerical evaluation.
The complete Vakint evaluation entry points already perform these steps in order.

Evaluation makes one pass at the requested precision plus the fixed guard allowance. It does
not retry after cancellation. Structured warnings identify limiting constants or master
coefficients, epsilon orders, requested and achieved precision, and the data version. Increase
the requested precision explicitly when cancellation leaves too few useful digits; a finite
source limit still cannot be overcome by adding working precision.

The optional second result remains the backend-reported statistical uncertainty series. It is
separate from numerical rounding and source-data uncertainty. Compare corresponding powers;
do not collapse the series to one complex number or silently discard either uncertainty.

Rust `to_atom` and Python `to_expression` are checked conversions. They round each nonzero
component to its achieved precision, preserve exact zeros, and reject uncertain zeros or
components without reliably encodable digits. The raw central values and precision records
remain inspectable, including in notebook displays, when conversion fails. Python `to_list()`
is an explicitly lossy conversion to ordinary Python floats; use named coefficient records
to retain precision information.

The #source-link("crates/vakint/precision_validation.typ", label: "precision validation report")
records numerical, graph and interface checks separately from pending scientific and release gates.

== Inspecting registered constants and masters

`get_constant(symbol, prec)` and `get_master(symbol, prec)` inspect the active RustRed registry;
they do not inspect FORM's private tables or trigger downloads. Rust takes an explicit decimal
precision, while Python defaults to 16. Discover qualified identities with `constant_symbols()`
and `master_symbols()`. Friendly names must be unambiguous within the snapshot.

Constant records expose `AvailablePrecision::ArbitraryPrecision` or
`AvailablePrecision::FixedPrecision(maximum_digits)`, a numerical value, achieved precision,
warnings, provenance, and `definition`: the exact defining multiple integral, written with the
inert `Integral`, `LoopIntegral` and `SeriesCoefficient` heads and explicit normalization
factors described in the masters documentation. An exact expression depending on a finite period inherits that
period's limit. Master records contain every available Laurent order in their recorded
unit-mass normalization, independent of the engine's usual expansion depth. Each coefficient
has its own availability and either `MasterSymbolicExpression::Expression` or
`MasterSymbolicExpression::NumericalOnly`. Declared uncomputed terminals have empty coefficient
maps; missing orders are never supplied as zeros.

The four-loop pack supplies finite `K5Moment(1)`, `K5Moment(3)`,
`K5LogMoment(1)` and `K5LogMoment(3)` values with certified source precision.
Exact PR4-family formulas using these symbols inherit their finite limits;
the function names alone do not imply arbitrary-precision evaluation.

Lookups cap numerical values at source availability and report the limitation. Returned records
retain the original immutable snapshot after another release is loaded. Their FeynKit graphs
preserve the authenticated master topology, powers, numerator and momentum routing. Structured
terminal and notebook representations use only these stored records; displaying an object
does not download data or reevaluate it.

Use `partial_numerical_evaluation` when you intentionally want to substitute known constants and
retain unresolved symbolic factors for inspection. Use `numerical_evaluation` for the final
boundary: it reports an error when tensor structure or a required symbol remains. This distinction
prevents a partially substituted expression from being mistaken for a fully numerical result.

#callout("Interpret failures at the owning stage", [
  An engine-construction error means the topology library could not be initialized. A
  canonicalization error means the propagators did not match a supported topology. A reduction
  error points to the numerator, Lorentz rank, or FORM translation; retain the temporary directory
  in that case. An evaluation error after successful reduction means settings validation, the
  selected backend's capabilities, or the backend invocation failed.
])

#boundary("Python is an embedded community module", [
  `symbolica.community.hepkit.vakint` is registered into a Symbolica installation; it is not a
  standalone PyPI package. Constructing its `Vakint` class validates the configured backends.
  Pass an empty evaluation order for pure matching work on machines without FORM/pySecDec.
])

See the supported-topology and external-dependency reference for the patterns and minimum tool
versions available here. Their Rust definitions begin in
#source-link("crates/vakint/src/topologies.rs", label: "Vakint's topology module").

== Methods and software to cite

Vakint combines distinct methods rather than treating every backend as interchangeable. Cite the
software version and the method actually selected for the reported result:

- #link("https://arxiv.org/abs/1203.6543")[FORM] for the symbolic reduction engine;
- #link("https://arxiv.org/abs/hep-ph/0009029")[MATAD] when its massive tadpole tables are used;
- #link("https://arxiv.org/abs/1707.01710")[FMFT] for the four-loop fully massive
  tadpole path; and
- #link("https://arxiv.org/abs/1703.09692")[pySecDec] for numerical sector
  decomposition.

Vakint also uses #link("https://symbolica.io/")[Symbolica] for expression manipulation. Record the
Vakint revision, normalization, epsilon depth, precision, selected backend and dependency
versions with the result; a generic citation to the package does not encode those choices.
]
