#import "../../shared.typ": boundary, developer-link, product-link

#let overview = [
= Overview

FeynKit is a standalone particle-physics toolkit for validated models, Feynman diagrams,
deterministic generation, Cross-Free Family (CFF) expressions, relativistic kinematics, and
vacuum tensor reduction. Rust clients use the focused crates or the `feynkit` facade. Python
clients use `symbolica.community.feynkit` inside a shared Symbolica kernel.

== Choose a task

- Generate and inspect a small diagram with the #link("quickstart/rust/")[Rust] or
  #link("quickstart/python/")[Python] quickstart.
- Load a UFO model, configure generation, or retain symbolic couplings with the
  #link("tutorial/")[model-to-diagram tutorial].
- Reduce vacuum numerators while keeping external momenta distinct with the
  #link("guides/tensor-reduction/")[tensor-reduction guide].
- Render diagrams and inspect CFF results with the #link("guides/notebooks/")[notebook guide].
- Embed the Python module in a Symbolica distribution using the
  #link("guides/community-host/")[community-host integration guide].

#boundary("Where FeynKit ends", [
  FeynKit returns reusable models, finalized diagrams, and symbolic results. It does not own a
  GammaLoop state, numerical integration, event accumulation, or a complete scattering
  calculation. #product-link("gammaloop", label: "GammaLoop") supplies those application layers.
])

== Shared representations

A `FeynmanDiagram` retains its model, particle and interaction assignments, numerator,
external-state projector, scalar prefactor, loop routing, and physical cuts. Generation builds
these once; downstream applications add their own numerical caches. CFF surface IDs resolve
through the arena returned with the expression, including when related expressions share it.

FeynKit builds on #product-link("linnet", label: "Linnet") graphs and uses
#product-link("spenso", label: "Spenso") tensor notation with
#product-link("idenso", label: "Idenso") identities. Its vacuum projector is also the default
numerator-reduction backend in #product-link("vakint", label: "Vakint"); scalar-integral
matching and evaluation remain Vakint's responsibility.

The #link("reference/interfaces/")[interface guide] links all eight Rust components and the
Python reference. Contributors can follow the
#developer-link("gammaloop-architecture", "architecture-current.typ", "ownership architecture")
and #developer-link("cff-surface-cache-ownership", "cff-surface-cache-ownership.typ", "CFF arena invariants").
]
