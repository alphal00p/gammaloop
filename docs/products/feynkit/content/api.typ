#import "../../shared.typ": catalog-contract

#let api = [
= Rust and Python interfaces

#catalog-contract(
  rust-scope: "feynkit, feynkit-model, feynkit-ufo, feynkit-kinematics, feynkit-graph, feynkit-generator, feynkit-cff, feynkit-tensor",
  python-scope: "symbolica.community.feynkit",
)

== Rust ownership map

- #link("reference/rust/feynkit/")[`feynkit`] re-exports the focused crates through feature gates.
- #link("reference/rust/feynkit_model/")[`feynkit-model`] owns validated models, parameter cards,
  typed IDs, and the explicit recomputation protocol.
- #link("reference/rust/feynkit_ufo/")[`feynkit-ufo`] loads UFO models through a caller-owned Python
  interpreter and `ufo_model_loader`.
- #link("reference/rust/feynkit_kinematics/")[`feynkit-kinematics`] owns momenta, transformations,
  signatures, helicities, and generalized-kt clustering.
- #link("reference/rust/feynkit_graph/")[`feynkit-graph`] owns finalized diagrams, cuts, routing,
  serialization, and Linnest source output.
- #link("reference/rust/feynkit_generator/")[`feynkit-generator`] owns process selectors,
  generation options, progress/cancellation, and diagram grouping.
- #link("reference/rust/feynkit_cff/")[`feynkit-cff`] owns CFF expressions, surface arenas,
  orientations, residues, and topology conversion.
- #link("reference/rust/feynkit_tensor/")[`feynkit-tensor`] owns vacuum projection, contraction
  orbits, and exact coefficient tables.

The facade's default features are `cff`, `generator`, `graph`, `kinematics`, `model`, and
`tensor`. `ufo` is opt-in. Native Rustdoc for this product includes that optional interface;
applications choose their features explicitly. When using `feynkit-model` in isolation, enable
its `native` feature to select the GMP integer and MPFR floating-point backends. The isolated
kinematics reference selects Numerica's `integer-gmp` and `float-mpfr` features explicitly.
The APIs shown here correspond to the source
revision of this manual, so pin matching dependencies when reproducing a calculation.

== Python ownership map

The #link("reference/python/feynkit-community/")[generated Python reference] covers native
classes, properties, signatures, examples, and error types. The primary owners are `Model`,
`Process`, `Generator`, `SnailFilterOptions`, `NumeratorGrouping`, `FeynmanDiagram`, `CffGenerator`, `TensorReducer`,
`JetDefinition`, and the optional `UfoLoader`. Expressions cross the boundary as Symbolica values.

Use the #link("quickstart/python/")[quickstart] for an installed host or the
#link("guides/community-host/")[community-host guide] when packaging the module. The Python
adapter does not publish a standalone wheel. Exception classes keep model, generation, graph,
CFF, kinematics, and tensor-reduction failures distinguishable.
]
