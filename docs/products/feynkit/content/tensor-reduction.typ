#import "../../shared.typ": callout, product-link

#let tensor-reduction = [
= Vacuum tensor reduction and momentum selectors

FeynKit reduces Lorentz-covariant vacuum tensor numerators to scalar Spenso invariants, retaining
metric tensors when free Lorentz indices remain. Vectors carry `spenso::mink(D,index)` slots;
compact vectors omit the index, and scalar contractions are represented by `spenso::dot`.
The reducer's dimension must match every input slot exactly.

== Select what is integrated

`TensorReducer.feynkit(dimension)` selects every `FeynKit::Momentum` tensor. Use it only for a
pure vacuum numerator whose momenta are all integrated. A generated scattering diagram still
contains external momenta; calling it a vacuum diagram does not change that. Construct
`TensorReducer(dimension)` and add exact compact vectors with `with_integrated_vector(...)`
when integrated and external momenta share a head.

This self-contained example integrates `k` and leaves `p` external:

// docs-example: compile feynkit-tensor-selectors
```python
from symbolica import S
import symbolica.community.feynkit as fk

D = S("feynkit_docs::D")
mu = S("feynkit_docs::mu")
nu = S("feynkit_docs::nu")
k = S("feynkit_docs::k")
p = S("feynkit_docs::p")
mink = S("spenso::mink")
dot = S("spenso::dot")
k_compact = k(mink(D))
p_compact = p(mink(D))
numerator = (
    k(mink(D, mu)) * k(mink(D, nu))
    * p(mink(D, mu)) * p(mink(D, nu))
)
reducer = fk.TensorReducer(D).with_integrated_vector(k_compact)
scalar = reducer.reduce(numerator)
expected = dot(k_compact, k_compact) * dot(p_compact, p_compact) / D
assert scalar == expected
```

For native generated momenta, select the corresponding exact
`FeynKit::Momentum(edge_id,spenso::mink(D))` vectors from your vacuum/routing construction.
Selecting an entire head with `with_integrated_head(...)` is appropriate only when every vector
under that qualified name is integrated. Ordinary generated FeynKit rules carry dimension `4`;
a Taylor-expanded expression whose slots carry symbolic `D` requires that same `D` instead.

== Reduce a diagram or a standalone expression

`reducer.reduce(expression)` transforms a numerator without a graph.
`diagram.reduce_tensor_numerator(reducer)` multiplies the stored numerator by its external-state
projector before reducing. It returns a Symbolica expression, including residual metrics.
`diagram.reduce_tensor_graphs(reducer)` requires a fully contracted result and returns one
scalar graph per compact term. It consumes and resets the projector, preserves the scalar
numerator prefactor and topology ID, and assigns deterministic `name.tensor[index]` names.

The Rust facade exposes `TensorReducer` and `FeynmanDiagramTensorExt` with the same ownership
boundary. Import the extension trait to call diagram reduction methods. FORM is not required
for this native projection; #product-link("vakint", page: "guides/evaluation/", label: "Vakint")
adds vacuum-topology matching and scalar-integral evaluation with separately selected backends.

== Rank, symmetry, and limits

Exact orthogonal-Weingarten coefficients support even ranks through 20. They depend on
integer-partition classes rather than every labeled pairing: rank 20 needs 42 coefficient
classes. Fully contracted numerators with repeated vectors use smaller contraction-orbit
systems; unsymmetrized free-index output can still contain factorially many terms.
The pairing, pairing-product, and output-term budgets make that limit explicit.

#callout("Keep mixed high-rank dimensions symbolic", [
  Fixed low integer dimensions can make the universal metric basis singular. Retain `D` or a
  dimensional-regulator expression through mixed high-rank reduction and substitute afterward.
  The all-equal isotropic fast path does not require that matrix inverse.
])

Consult the #link("reference/rust/feynkit_tensor/")[Rust tensor reference] for coefficient
engines and error variants, and the #link("reference/python/feynkit-community/TensorReducer/")[Python
reducer reference] for selector and budget methods.
]
