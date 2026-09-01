#import "../../shared.typ": callout, boundary, developer-link, product-link

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

== Choose a task

- To install Symbolica and construct a typed tensor with the bundled Python module, use the
  #link("quickstart/python/")[Python guide].
- To verify slot duality with one numerical Rust contraction, follow the
  #link("quickstart/rust/")[Rust guide], then continue with the
  #link("tutorial/")[first-contraction tutorial] and the exact
  #link("reference/rust/spenso/contraction/trait.Contract.html")[`Contract` reference].
- To assemble and execute several tensor/scalar nodes, use the
  #link("guides/networks/")[tensor-network guide] with the
  #link("reference/rust/spenso/network/struct.Network.html")[`Network` reference].
- To construct tensors or execute a network from Python, follow the
  #link("guides/python/")[Python tensor-workflow guide] with the exact
  #link("reference/python/spynso3/TensorNetwork/")[`TensorNetwork` reference].
- To parse symbolic tensor functions or apply Dirac/color identities, review the
  #link("reference/interfaces/")[feature and interface guide], then continue with
  #product-link("idenso", label: "Idenso") for representation-aware rewrites.

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

// docs-example: syntax
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
  feature gates are independent, so select and upgrade each component explicitly.
])

Spenso uses #product-link("linnet", label: "Linnet") for the underlying network graph.
#product-link("idenso", label: "Idenso") is the symbolic-identity layer for Spenso-formatted
Symbolica expressions. #product-link("gammaloop", label: "GammaLoop") consumes these components
inside a broader collider workflow.

The #link("reference/interfaces/")[interface guide], native Rustdoc, and generated Python
reference list the public modules, traits, feature gates, and import paths.

Contributors changing structures, contraction, parsing, storage, or network scheduling should
also read the source-audited
#developer-link("spenso-architecture", "spenso-architecture.typ", "Spenso implementation architecture").
]
