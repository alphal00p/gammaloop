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

Tensor reduction validates Lorentz index domains before vacuum projection. It passes universal loop-momentum tensor kernels to FORM while retaining graph coefficients as factorized expressions. Analytic integration uses temporary coefficient symbols and restores their values afterward; numerical integration combines the scalar sectors into one estimator. Numerator and loop-normalization poles in epsilon increase the required integral order. Positive epsilon prefactors never reduce the conservative requested depth. Rational dependence and scalar log/exp/sin/cos series are supported; retained Laurent coefficients must be epsilon-independent. Opaque epsilon-dependent functions produce an error.

Callers with additional spin tensors declare their Lorentz slots using `vakint::tensor(body, slot, ...)`. The loop-independent body stays opaque; metrics rename complete slots without modifying spin or color indices. Projection can introduce `vakint::tensor_index(family, loop, ordinal)` dummy slots, which callers translate to their own index representation before unwrapping tensor bodies. GammaLoop performs both conversions automatically and completes any new d-dimensional Dirac contractions before scalar backend evaluation truncates epsilon. Vakint does not evaluate opaque spin algebra.

`VakintSettings.project_onto_tensor_integrals` defaults to `true`. Setting it to `false` sends each complete numerator through FORM tensor reduction and the scalar backend, instead of using universal tensor kernels and coefficient aliases. Distinct denominator topologies remain separate in both modes. GammaLoop exposes this choice as `global.generation.uv.project_integrated_uv_cts_onto_tensor_integrals`; both choices share its d-dimensional algebra and gamma5 checks.

`PySecDecOptions.numerical_parameters` stores mass and numerator parameters as complex values. Real external momentum components retain their signs. This adapter requires real pole masses and real external momenta; complex numerator coefficients are supported. Physical PySecDec integration tests are manual only: they are ignored and excluded from automatic local/CI test selection. Missing opt-in or missing PySecDec causes an explicit manual invocation to fail instead of reporting an unevaluated test as passed. Pure Rust adapter tests remain in the normal suite because they do not launch PySecDec.

Run the four fast one-loop numerical checks explicitly when changing that backend:

```bash
RUN_PYSECDEC_TESTS=1 cargo nextest run -p vakint \
  --cargo-profile dev-optim -P test_gammaloop --test integral_evaluation_pysecdec_tests \
  --run-ignored only --ignore-default-filter \
  -E 'test(/^test_integrate_1l_(simple|complex_and_signed_parameters|cross_product|cross_product_with_additional_symbols_numerator)$/)'
```

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
`externals_from_f64` convert ordinary inputs to the configured decimal precision before the
series is evaluated.

`NumericalEvaluationResult` stores sorted `(epsilon power, complex coefficient)` pairs. Negative
powers are poles, power zero is the finite term, and positive powers are higher orders. The
optional second result carries a backend-reported uncertainty series when the evaluated
expression contains one. Compare corresponding powers; do not collapse the series to one complex
number or silently discard the uncertainty.

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
  `symbolica.community.vakint` is registered into a Symbolica installation; it is not a
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
