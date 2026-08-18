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

The evaluation order is an ordered fallback policy, not a request to combine backend results.
Vakint tries each configured method whose capability matches the canonical topology and stops
at the first successful evaluation. AlphaLoop, MATAD, and FMFT are analytic FORM-backed paths.
pySecDec is numerical and additionally needs numeric parameters and external momenta where the
topology requires them.

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

// docs-example: compile
```rust
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

    println!("{evaluated}");
    Ok(())
}
```

The result invariant is an MS-bar Laurent series in the dimensional regulator: the reduced
expression no longer contains the original loop-momentum numerator, and only the first
successful method in `alphaloop_only()` contributes. Inspect the exact
#link("reference/rust/supported/vakint/#supported-vakintsettings")[`VakintSettings`],
#link("reference/rust/supported/vakint/#supported-vakintexpression-tensor-reduce")[tensor
reduction], and
#link("reference/rust/supported/vakint/#supported-vakintexpression-evaluate-integral")[evaluation]
references before changing normalization or backend order.

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
]
