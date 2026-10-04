#import "../../shared.typ": callout, developer-link, product-link

#let algebra = [
= Syntax, indices, and algebra passes

Idenso consumes the Symbolica function form emitted for Spenso structures. Function heads carry
representation meaning and their arguments carry abstract indices or tensor data. Import the
Spenso community module before parsing so Symbolica attributes, Spenso tags, and dual
representations are registered in one deterministic order.

== Indices and cooking

Free indices describe the result; repeated compatible indices describe contractions. Before
combining independently built expressions, wrap or rename dummy-index namespaces so equal
printed names do not create accidental contractions. Cooking temporarily replaces selected
index payloads or function subexpressions with compact symbols. Rust `CookSettings::uncook`
restores reversible encodings. Python uses `intern="indices"` for reversible index encoding or
`intern="flattened"` for readable index names when constructing or indexing a tensor.
The default `intern=None` leaves indices unchanged; neither mode hides tensor heads.

Cooking changes what downstream matchers can see. Retain the original expression when a later
identity needs the hidden function head or index payload. Typed dummy scoping below does not
hide those tensor structures.

== Keep independent dummy namespaces independent

Two factors may legitimately print the same local dummy name before they are multiplied. Wrap
each factor with a distinct header first, while leaving its external indices untouched:

// docs-example: compile idenso-dummy-namespaces
```python
import symbolica as sp
from symbolica.community.tensor import (
    Representation,
    TensorExpression,
    TensorName,
)

rep = Representation.euc(3)
mu = rep("mu")
nu = rep("nu")
rho = rep("rho")
g = TensorExpression.g(rep)
p = TensorName.vector("p")
q = TensorName.vector("q")

left = TensorExpression(g(mu, nu).to_expression() * p(mu).to_expression())
right = TensorExpression(g(mu, rho).to_expression() * q(mu).to_expression())
safe_product = (
    left.wrap_indices(sp.S("lhs"), dummies_only=True)
    * right.wrap_indices(sp.S("rhs"), dummies_only=True)
)

assert len(safe_product.list_dangling()) == 2
assert (
    safe_product.contract(rank_one=False).to_expression()
    == p(nu).to_expression() * q(rho).to_expression()
)
```

The stable invariant is two free indices, `nu` and `rho`; the two occurrences of local `mu`
belong to separate contractions after wrapping. `wrap_indices(..., dummies_only=True)`
returns a typed tensor with scoped dummy indices and the original external ports. Its default
scopes every explicit index; unresolved ports remain unresolved in either mode. No cooking
roundtrip is required. The generated
#link("reference/python/spynso3/TensorExpression/#exports-tensorexpression-wrap-indices-method")[`wrap_indices` reference]
records the signature; the shared Rust `SymbolicTensor` owns the transformation.

#callout("Interpret index failures before simplifying", [
  More or fewer than two dangling indices means a name collided or a slot's representation or
  duality differs from the intended one. An unchanged plain Symbolica function means it was not
  constructed through Spenso tensor names/representations, or the Spenso module was imported
  only after parsing.
  Correct those structural issues before metric, Dirac, or color simplification.
])

Bracketed products with explicit indices participate directly in metric, Dirac,
epsilon, and color simplification. The rules see contractions across the bracket's
arguments; scalar results and zero lose their wrappers. Unresolved contracted
products retain a bracket when needed to keep powers and denominators atomic.
This normalization preserves unrelated scalar sums and their factorization.

== Metric and epsilon operations

Metric contraction raises, lowers, or identifies compatible Lorentz indices according to the
registered representation. Epsilon identities depend on dimension, ordering, and sign
conventions; expand them only when the next pass benefits from the larger expression. Keep
Minkowski expansion selective to avoid distributing unrelated scalar factors.

== Dirac and color algebra

Dirac passes simplify gamma chains, traces, slashes, and spinor-compatible contractions. Color
passes handle fundamental/adjoint deltas, generators, structure constants, and registered
group parameters. Select the required families explicitly with `simplify_algebra(gamma=True, color=True, epsilon=True)`; the shared planner chooses
their order and revisits only affected domains. Inspect intermediate output when comparing
conventions, without prescribing a global pass pipeline.

For explicit color generators, `dirac_adjoint` exchanges
fundamental and antifundamental slots and transposes the generator ports. This
uses the Hermiticity of the SU(N) generators. Scalar representation labels in
Casimir and index invariants are preserved. Real momenta and couplings still
need explicit assumptions or substitutions for unevaluated conjugations.

Compact fundamental chains and ordered traces of these Hermitian generators
use the same conjugation: reverse the generator sequence and, for open chains,
exchange and dualize the endpoints. Color simplification with
`simplify_algebra(color=True, color_evaluate_traces=False, color_expand_fierz=False)`
collects explicit generator words while retaining compact chains and traces;
this commutes with conjugating an explicit SU(3) network. If both operands
remain explicit, scope their internal dummies before multiplying a word by
its adjoint, as above. A trace of three generators is not assumed real. The
existing trace builder retains cyclic normalization. Symmetric, antisymmetric
and cyclic groups are traversed recursively, keeping the projectors compact.
Reversing an antisymmetric group
produces its permutation sign through Spenso's existing normalization. Numeric
coefficients and scalar symbols, including their sums and products, are
conjugated within groups.
For three generators, the symmetric trace is real and the antisymmetric trace
is purely imaginary. These rules recognize generator words; they do not assign
Hermiticity to arbitrary matrix-valued functions.

Color simplification contracts metrics inside collected chains and traces before
applying trace identities. A symmetric trace through degree four with a repeated
adjoint pair is expanded using Spenso's normalized projector and reduced by the
existing Casimir rules. Open symmetric invariants and higher-degree projectors
remain symbolic; this does not claim a general higher-rank invariant reduction.

#callout("Canonical does not mean physically equivalent by itself", [
  Canonical ordering makes structurally equivalent expressions comparable under the registered
  rules. On-shell relations, gauge choices, dimension-specific identities, and model parameter
  substitutions must still be requested explicitly.
])

== Shared structural contraction

`contract()` reads the borrowed Spenso factor graph and applies Idenso’s shared component
reducer. It retains generated sums in factored form and uses the ordered substitution owner
when callbacks require it. Scalar coefficients stay factored. Concrete tensor-network execution
remains a separate Spenso operation when component data are required. Avoid substituting scalar
parameters too early: this can increase expression size without changing tensor incidence.

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
