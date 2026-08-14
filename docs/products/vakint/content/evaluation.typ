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

#boundary("Python is an embedded community module", [
  `symbolica.community.vakint` is registered into a Symbolica installation; it is not a
  standalone PyPI package. Constructing its `Vakint` class validates the configured backends.
  Pass an empty evaluation order for pure matching work on machines without FORM/pySecDec.
])

See the supported-topology and external-dependency reference for the patterns and minimum tool
versions available here. Their Rust definitions begin in
#source-link("crates/vakint/src/topologies.rs", label: "Vakint's topology module").
]
