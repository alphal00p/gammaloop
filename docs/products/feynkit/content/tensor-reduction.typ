#import "../../shared.typ": callout, product-link

#let tensor-reduction = [
= Covariant tensor reduction and momentum selectors

FeynKit reduces Lorentz-covariant tensor numerators to scalar Spenso invariants, retaining
metric tensors when free Lorentz indices remain. Vectors carry `spenso::mink(D,index)` slots;
compact vectors omit the index, and scalar contractions are represented by `spenso::dot`.
The reducer's dimension must match every Lorentz slot; native encoded dimensions
are interpreted through their scalar source payloads.
Without an external basis, the denominator must have vacuum rotational symmetry.

== Select what is integrated

`TensorReducer.feynkit(dimension)` selects every `gammalooprs::Q` tensor. Use it only for a
pure vacuum numerator whose momenta are all integrated. A generated scattering diagram still
contains external momenta; calling it a vacuum diagram does not change that. Construct
`TensorReducer(dimension, integrated=[...])` with exact compact vectors
when integrated and external momenta share a head.

This self-contained example integrates `k` and leaves `p` external:

// docs-example: compile feynkit-tensor-selectors
```python
from symbolica import S
import symbolica.community.hepkit as fk

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
reducer = fk.TensorReducer(D, integrated=[k_compact])
scalar = reducer.reduce(numerator)
expected = dot(k_compact, k_compact) * dot(p_compact, p_compact) / D
assert scalar.rank == 0 and scalar.axes == ()
assert scalar.to_expression() == expected
```

Compact contractions can be passed directly: the same reducer gives zero for
`dot(k_compact, p_compact)` and gives `expected` for
`dot(k_compact, p_compact)**2`. A product mixing compact dots and free-index
vectors is handled in the same projection, without an `undo_dots()` step.
Dots between two integrated vectors remain scalar weights, as do products with
declared external basis vectors. Nonnegative integer powers of spectator
contractions are polynomial tensor numerators; negative or noninteger powers
require a different integral treatment and are rejected. The existing rank
limit applies before repeated contraction legs are allocated.

For native generated momenta, select the corresponding exact
`gammalooprs::Q(edge_id,spenso::mink(D))` vectors from your vacuum/routing construction.
Selecting an entire head with `TensorReducer(D, integrated=[k])` is appropriate only when every vector
under that symbol is integrated. Pass a `TensorName` or a bare Symbolica `Expression`, such as
`S("Loop::k")`. Compact vectors may be `Expression` or `TensorExpression` objects.
Ordinary generated FeynKit rules carry dimension `4`;
a Taylor-expanded expression whose slots carry symbolic `D` requires that same `D` instead.

Exact compact selectors are source inputs. The constructor admits each selector
once into the same native index and dimension encoding used by retained
numerators. This matters for a compound dimension such as `4 - 2*epsilon`:
`Q(0,mink(4 - 2*epsilon))` must match the native vector for edge zero while
leaving edge one under the same `Q` head unselected. Comparing the raw selector
to an encoded numerator by Atom spelling would silently miss that selection.
This admission does not decode or re-admit the numerator passed to `reduce`.

== Denominators with external momenta

Supply every independent external direction in the denominator through the constructor
keyword `external=[p_compact]`. Both `Expression` and `TensorExpression` entries are accepted.
The reducer decomposes loop vectors into
components parallel and perpendicular to that span. Symbolica inverts the
external Gram matrix, and the existing vacuum projector averages the transverse
components in dimension `D - n`, where `n` is the basis size. Odd total ranks can
then be nonzero.

// docs-example: compile feynkit-tensor-external-basis
```python
from symbolica import S
from symbolica.community.tensor import TensorExpression
import symbolica.community.hepkit as fk

D, k, p, mu = S("external_docs::D", "external_docs::k", "external_docs::p", "external_docs::mu")
mink, dot = S("spenso::mink", "spenso::dot")
kc, pc = k(mink(D)), p(mink(D))
reducer = fk.TensorReducer(D, integrated=[kc], external=[pc])
result = reducer.reduce(k(mink(D, mu)))
expected = TensorExpression(dot(kc, pc) / dot(pc, pc) * p(mink(D, mu)))
assert result.rank == 1 and result.axes == expected.axes
assert (result - expected).together().is_zero()
```

The Gram determinant must remain nonzero after imposing kinematics. For a
lightlike one-vector basis, supply an additional independent auxiliary direction
before reducing, then impose its scalar products on the result. The reducer
works with symbolic scalar products and cannot infer an unstated on-shell or
momentum-conservation relation. It leaves loop scalar products in the numerator;
it does not perform momentum shifts, partial fractions or IBP reduction.
External basis vectors take precedence over integrated-head selectors.

== Reduce a diagram or a standalone expression

`diagram.tensor_reduce(D)` selects the graph's internal edge momenta automatically, promotes
four-dimensional Lorentz slots in the numerator and projector to `D`, and returns a
`TensorExpression`. Dimension promotion uses the retained tensor's dimension
mapper, including for a compound `D`, preserving logical axis order and ranked
zero. The scalar dimension in projector coefficients remains the supplied
expression. External edge momenta remain projector vectors; they do not automatically define
the denominator's external basis. Use a configured `TensorReducer` for a
non-vacuum denominator. The scalar numerator
prefactor remains separate, and at least one internal edge is required.

Use `diagram.tensor_reduce(D, expression=prepared)` after contractions or a UV expansion.
The supplied expression replaces the stored numerator and projector and must use the graph's
`gammalooprs::Q(edge_id, ...)` names. It cannot be combined with an explicit `projector`.
Internal edge momenta are still selected from the diagram. For a `Subgraph`, call
`region.tensor_reduce(D, expression=prepared)` to integrate only its selected internal
momenta. An explicit projector is required when reducing a proper region's stored numerator.

Tensor reduction keeps loop-dependent vertex sums factored. It separates fixed
Lorentz tensors from the loop-dependent blocks and projects the open indices of
those blocks, rather than distributing all their internal contractions. The result
can therefore retain contracted indices inside sums. Metric and vector contractions
that preserve the sums are performed automatically; expanding those sums into scalar
dot-product monomials is a separate, potentially expensive operation.

For standalone tensors, `TensorExpression(expression)` performs the same conversion as
`as_tensor(expression)`. Algebraic rearrangements return tensor expressions with their ordered
interfaces and data identities preserved: `expand`, `expand_num`, `factor`, `collect`,
`collect_symbol`, `collect_num`, `collect_factors`, `collect_by_coefficient`, `collect_horner`,
`together`, `cancel`, and `apart`. Python's `copy.copy` also preserves tensor metadata.
Trusted algebraic rearrangements retain the established interface; collection callbacks
require checking the resulting interface and reject incompatible changes.
Tensor equality compares the ordered interface as well as the expression, so a raw
Symbolica expression is not equal to a `TensorExpression`. For a scalar oracle, check
that the result has no open axes and compare `result.to_expression()` with the raw
expected expression. For tensor-valued oracles, construct a typed expected value and
compare its axes; use `is_zero()` to test the exact difference.

`TensorExpression` subclasses Symbolica's `Expression`: inspection, matching, conversion,
and evaluation remain available directly. The tensor overrides of `replace`,
`replace_multiple`, `map`, and `derivative` preserve the ordered interface and reject
incompatible changes. Other inherited methods keep Symbolica's return types.
Use `to_expression()` to deliberately discard the tensor interface, for example before a
rank-changing rewrite. It returns an ordinary `Expression` without tensor metadata.
The expression can still contain interned index or dimension symbols;
`to_expression()` does not decode those source labels.
Use `tensor.with_lorentz_dimension(D)` before contracting indices to promote four-dimensional
Lorentz slots, including compact representations. Spinor/color dimensions, scalar coefficients,
and Lorentz slots already in other dimensions remain unchanged.

`reducer.reduce(expression)` transforms a numerator without a graph.
`diagram.reduce_tensor_numerator(reducer)` multiplies the stored numerator by its external-state
projector before reducing. Both return a `TensorExpression`, retaining ordered
ports and residual metrics, including the interface of a zero result.
`diagram.reduce_tensor_graphs(reducer)` requires a fully contracted result and returns one
scalar graph per compact term. It consumes and resets the projector, preserves the scalar
numerator prefactor and topology ID, and assigns deterministic `name.tensor[index]` names.

The Rust facade exposes `TensorReducer` and `FeynmanDiagramTensorExt` with the same ownership
boundary. Import the extension trait to call diagram reduction methods. FORM is not required
for this native projection; #product-link("vakint", page: "guides/evaluation/", label: "Vakint")
adds vacuum-topology matching and scalar-integral evaluation with separately selected backends.

`TensorReducer::admit` admits an external expression already converted to native
Spenso syntax and preserves the reducer's index and power diagnostics on rejected
input. It does not cook source notation or perform reduction. Existing tensor
owners call `reduce` directly; no re-admission is needed.

The Rust source-selector builders `with_integrated_vector(Atom)` and
`with_external_vector(Atom)` return `Result<Self, TensorReductionError>` because
their one-time source admission can fail. Propagate the result, for example
`TensorReducer::new(D).with_integrated_vector(k_compact)?`. Head selectors retain
their existing behavior. This is a separate boundary from `TensorReducer::admit`:
the latter expects an already converted native numerator, while selector builders
accept compact source notation. A raw compound dimension is not itself a Spenso
`Dimension` value; graph promotion uses the existing checked dimension mapper
rather than requiring callers to flatten that expression to a symbol.

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
