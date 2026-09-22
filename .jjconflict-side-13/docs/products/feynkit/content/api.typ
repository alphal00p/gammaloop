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
  signatures, helicities, symbolic scalar products, Mandelstam substitutions, and generalized-kt clustering.
- #link("reference/rust/feynkit_graph/")[`feynkit-graph`] owns finalized diagrams, cuts, routing,
  serialization, integral-family bases, and Linnest source output.
- #link("reference/rust/feynkit_generator/")[`feynkit-generator`] owns process selectors,
  generation options, progress/cancellation, diagram grouping, and shared external-state spin sums.
- #link("reference/rust/feynkit_cff/")[`feynkit-cff`] owns CFF expressions, surface arenas,
  orientations, residues, and topology conversion.
- #link("reference/rust/feynkit_tensor/")[`feynkit-tensor`] owns covariant tensor reduction, contraction
  orbits, and exact coefficient tables.

The facade's default features are `cff`, `generator`, `graph`, `kinematics`, `model`, and
`tensor`. `ufo` is opt-in. Native Rustdoc for this product includes that optional interface;
applications choose their features explicitly. When using `feynkit-model` in isolation, enable
its `native` feature to select the GMP integer and MPFR floating-point backends. The kinematics crate selects those backends through its default `native` feature;
use `default-features = false` and `wasm` for the portable backends.
The APIs shown here correspond to the source
revision of this manual, so pin matching dependencies when reproducing a calculation.

== Python ownership map

The #link("reference/python/feynkit-community/")[generated Python reference] covers native
classes, properties, signatures, examples, and error types. The primary owners are `Model`,
`Process`, `Generator`, `SnailFilterOptions`, `NumeratorGrouping`, `FeynmanDiagram`, `Subgraph`, `CffGenerator`, `TensorReducer`,
`Kinematics`, `IntegralFamily`, `JetDefinition`, and the optional `UfoLoader`. Expressions cross the boundary as Symbolica values.

Use the #link("quickstart/python/")[quickstart] for an installed host or the
#link("guides/community-host/")[community-host guide] when packaging the module. The Python
adapter does not publish a standalone wheel. Exception classes keep model, generation, graph,
CFF, kinematics, and tensor-reduction failures distinguishable.

== CFF coefficients and generalized cuts

`diagram.subgraph(selection).build_cff()` constructs the selected region's CFF;
`CffGenerator.generate` accepts the diagram or its `Subgraph` view.
`result.to_expression()` returns the same
canonical surface placeholders used by GammaLoop. Set `expand_surfaces=True`
to substitute energies, or `normalized=True` to include the internal
propagator energy factors and GammaLoop's loop measure. Numerators, projectors,
and global weights remain separate.

`result.raised_surface_groups(edge_representatives)` identifies equivalent
surfaces after mapping repeated propagators to their representative edge.
`result.pole_coefficients(group)` returns coefficients in increasing pole order;
these are not yet analytic residues. `result.residue` additionally requires the
integration variable, its root, the surface in that variable, and the complete
remaining coefficient. Surface objects and groups belong to the result that
created them; mixing results is rejected.

```python
q0 = Expression.symbol("q0")
cut = CutPropagator(q0, 2, power=3, normalization=1)
assert cut.apply(q0**2, q0) == -Expression.num(1) / 128
```

`CutPropagator` follows the generalized distributional cutting rules in
#link("https://arxiv.org/abs/2203.11038")[Local Unitarity, section 2.2].
Its default factor is $-2 pi i$ for a reciprocal propagator. The orientation
chooses the positive or negative energy root; the prescription specifies the
sign of the infinitesimal imaginary part. A raised cut differentiates both
the remaining integrand and the uncut energy factor before imposing the shell.
The energy-space symbol `Delta(n,x)` denotes the normalized divided-difference
distribution whose action is the derivative of order `n-1` divided by `(n-1)!`.
Route to an independent energy coordinate before calling `apply`.

]
