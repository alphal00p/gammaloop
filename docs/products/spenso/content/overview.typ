#import "../../shared.typ": callout, boundary, product-link, source-link

#let overview = [
= Overview

Spenso represents and evaluates tensors with abstract indices. It separates a tensor's
structure—representations, slots, dimensions, names, and index order—from its stored data.
Dense and sparse storage, symbolic and parametric values, automatic index matching, and
network execution build on that separation.

#callout("Two questions, two layers", [
  Tensor structure answers “which indices and representations does this object carry?” Tensor
  data answers “which values occupy that structure?” Contraction first matches compatible
  slots, then combines data. Keep these phases distinct when diagnosing a shape, duality, or
  contraction error.
])

== From tensors to networks

The core package supports three related levels of work:

- tensor structures describe slots and abstract indices;
- dense, sparse, numeric, symbolic, and parametric tensors attach data to those structures;
- tensor networks represent a collection of tensors and choose an execution or contraction
  strategy.

Networks can be assembled directly. With the `shadowing` feature, they can also be inferred
from Symbolica expressions whose functions carry recognizable tensor structure. Parsing is not
mere string conversion: representation initialization, index variance, names, and the tensor
library determine the inferred network.

```toml
[dependencies]
spenso = "0.6"
```

Enable `shadowing` only when Symbolica-backed parsing or symbolic operations are needed. It
pulls in the Symbolica integration and corresponding Linnet support. The default core remains
useful for typed tensors and contractions without that layer.

== Companion packages

#boundary("Component ownership", [
  `spenso` owns the tensor model and execution engine. `spenso-macros` owns the
  `SimpleRepresentation` derive used to declare custom index representations.
  `spenso-hep-lib` owns concrete high-energy-physics tensor-library data such as gamma matrices
  and projectors. The `spynso3` adapter owns the Python community module. Their versions and
  feature gates are tracked independently in the product registry.
])

Spenso uses #product-link("linnet", label: "Linnet") for the underlying network graph.
#product-link("idenso", label: "Idenso") is the symbolic-identity layer for Spenso-formatted
Symbolica expressions. #product-link("gammaloop", label: "GammaLoop") consumes these components
inside a broader collider workflow.

Start with the #source-link("crates/spenso/README.md", label: "Spenso README") for the package
summary and use the generated catalog for the exact compiled surface.
]
