#import "../../shared.typ": callout, developer-link, product-link

#let algebra = [
= Syntax, indices, and algebra passes

Idenso consumes the Symbolica function form emitted for Spenso structures. Function heads carry
representation meaning and their arguments carry abstract indices or tensor data. Import the
Idenso community module before parsing so Symbolica attributes, Spenso tags, and dual
representations are registered in one deterministic order.

== Indices and cooking

Free indices describe the result; repeated compatible indices describe contractions. Before
combining independently built expressions, wrap or rename dummy-index namespaces so equal
printed names do not create accidental contractions. Cooking temporarily replaces selected
index payloads or function subexpressions with compact symbols. Retain the generated map and
settings whenever uncooking is required.

Cooking is especially useful before expensive canonicalization, but it changes what downstream
matchers can see. Uncook before applying an identity whose pattern depends on the hidden tensor
head or index structure.

== Keep independent dummy namespaces independent

Two factors may legitimately print the same local dummy name before they are multiplied. Wrap
each factor with a distinct header first, while leaving its external indices untouched:

// docs-example: compile idenso-dummy-namespaces
```python
import symbolica as sp
from symbolica.community.idenso import list_dangling, wrap_dummies
from symbolica.community.spenso import Representation, TensorName

rep = Representation.euc(3)
mu = rep("mu")
nu = rep("nu")
rho = rep("rho")
g = TensorName.g()
p = TensorName("p")
q = TensorName("q")

left = g(mu, nu) * p(mu)
right = g(mu, rho) * q(mu)
safe_product = wrap_dummies(left, sp.S("lhs")) * wrap_dummies(right, sp.S("rhs"))

assert len(list_dangling(safe_product)) == 2
```

The stable invariant is two free indices, `nu` and `rho`; the two occurrences of local `mu`
belong to separate contractions after wrapping. The generated
#link("reference/python/idenso-community/wrap_dummies-function/")[`wrap_dummies` reference] records the
Python signature, while the exact
#link("reference/rust/idenso/trait.IndexTooling.html")[`IndexTooling` Rustdoc] covers
the underlying Rust boundary.

#callout("Interpret index failures before simplifying", [
  More or fewer than two dangling indices means a name collided or a slot's representation or
  duality differs from the intended one. An unchanged plain Symbolica function means it was not
  constructed through Spenso tensor names/representations, or the Idenso module was imported
  only after parsing.
  Correct those structural issues before metric, Dirac, or color simplification.
])

== Metric and epsilon operations

Metric contraction raises, lowers, or identifies compatible Lorentz indices according to the
registered representation. Epsilon identities depend on dimension, ordering, and sign
conventions; expand them only when the next pass benefits from the larger expression. Keep
Minkowski expansion selective to avoid distributing unrelated scalar factors.

== Dirac and color algebra

Dirac passes simplify gamma chains, traces, slashes, and spinor-compatible contractions. Color
passes handle fundamental/adjoint deltas, generators, structure constants, and registered
group parameters. Apply one algebra family at a time and inspect the intermediate expression;
an all-at-once fixed-point loop can obscure which convention produced a sign or normalization.

#callout("Canonical does not mean physically equivalent by itself", [
  Canonical ordering makes structurally equivalent expressions comparable under the registered
  rules. On-shell relations, gauge choices, dimension-specific identities, and model parameter
  substitutions must still be requested explicitly.
])

== Schoonschip network parsing

The Schoonschip path parses large gamma-chain expressions into the same Spenso-compatible
network model before contraction. It resolves aliases and normalizes the expression before
constructing the network. Spenso then plans and executes contractions, while Idenso supplies
the representation and algebra rules. Avoid substituting scalar parameters too early: doing so
can substantially increase intermediate expression size even when the resulting tensor network
is unchanged.

Concrete syntax and rewrite cases are documented in the
#developer-link(
  "spenso-symbolica-syntax-and-rewrites",
  "spenso-symbolica-syntax-and-rewrites.typ",
  "Spenso/Symbolica syntax note",
)
and the rendered #link("reference/form-color-dirac/")[shipped color and Dirac convention reference].
The #developer-link(
  "schoonschip-network",
  "schoonschip-net-parsing.typ",
  "Schoonschip parsing guide",
)
shows how normalization, network construction, and contraction fit together.

Tensor storage and execution remain owned by #product-link("spenso", label: "Spenso").
]
