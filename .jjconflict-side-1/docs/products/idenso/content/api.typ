#import "../../shared.typ": boundary, source-link

#let api = [
= Rust and Python APIs

== Rust package

The `idenso` crate exposes representation types and syntax macros together with several
rewrite families:

- `IndexTooling` covers canonicalization, conjugation, index wrapping, and dangling-index
  inspection for Symbolica atoms;
- `Cookable`, `CookSettings`, and the cook filters control reversible or flattening encodings;
- `SymbolicTensor::collect` retains selected sectors in the shared alias registry;
  `coefficient_list` returns typed selected sectors and their factored coefficients;
- `dirac`, `color`, `epsilon`, and shorthand modules implement algebra-specific rewrites;
- `representations::initialize` installs the standard representation and tensor symbols.

Representation helper macros such as `bis!`, `cof!`, and `coad!` construct the symbolic forms
expected by Spenso and Idenso. The #link("reference/rust/")[Rust orientation] leads to the
revision-specific Rustdoc for their accepted forms, return types, feature gates, and source
locations. APIs behind `bincode` and `reference-cases` are available only when the matching Cargo
feature is enabled. Python bindings belong to `spynso3`.

== Python community module

#boundary("Part of Symbolica community", [
  Import `TensorExpression` and representation types from `symbolica.community.tensor`.
  The `spynso3` Cargo package supplies the unified Python bindings for Spenso and Idenso;
  Idenso remains the Rust algebra implementation. There is no separate `idenso` Python module.
])

Install a Symbolica community build with the unified `spynso3` bindings, then verify it with
`python -c "import symbolica.community.tensor"`. There is no `pip install idenso` fallback.
Source embedders add `spynso3` to their
#link("https://github.com/symbolica-dev/symbolica-community")[symbolica-community] assembly and
register `SpensoModule` through `SymbolicaCommunityModule`. Building the Rust crate alone does
not add the community module to an already installed Symbolica package.

The generated #link("reference/python/")[Python API] records exact signatures and defaults. Its
operations cover:

- setup: importing the community module registers its symbols;
- collection and materialization: `collect` preserves factored sectors; `expand` explicitly
  materializes an expanded result;
- index preparation: typed `wrap_indices`, `list_dangling`, and the explicit
  `intern="indices"` or `intern="flattened"` construction/indexing policy;
- reduction: `contract(...)` and `simplify_algebra(...)`;
- notation: `to_dots`, `undo_dots`, `undo_chain`, and `undo_trace`;
- canonical labels: `canonize`;
- conjugation: `dirac_adjoint`.

```python
from symbolica.community.tensor import Representation, TensorExpression, TensorName

minkowski = Representation.mink(4)
mu = minkowski("mu")
nu = minkowski("nu")
metric = TensorExpression.g(minkowski)
momentum = TensorName.vector("p")
expression = metric(mu, nu) * momentum(mu)

external_indices = expression.list_dangling()
reduced = expression.contract().to_expression()
assert len(external_indices) == 1
assert len(reduced.list_dangling()) == 1
assert reduced == momentum(nu)
```

`contract()` returns an `AliasedTensorExpression` and reports whether contraction is complete.
Use `contract(rank_one=False, collect_chains=False, collect_traces=False)` for metrics only. `contract_ports(rhs, left=..., right=...)`
is the separate binary operation on selected logical ports. `to_dots()` changes notation;
request `contract()` first when repeated explicit indices need contraction.
Resolve aliases with `to_expression()`, or explicitly distribute the result with `expand()`.

`simplify_algebra(...)` selects gamma, color, and epsilon identities
independently. Its default enables no identity family; `simplify_algebra(gamma=True, color=True, epsilon=True)` enables
all three. Flat keyword arguments such as `gamma_ordering="canonical"`,
`gamma_evaluate_traces=False`, and `color_expand_fierz=False` specify family
conventions and requested output. Supplying a family option requires explicitly
enabling that family. Dimensions belong to representations. Enabled identities authorize their necessary
structural contractions; emitting an epsilon does not itself enable epsilon identities.
The default `contract="fully"` also performs structural contractions, including
closed colour deltas, without enabling colour generator identities. No trailing
`contract()` is needed to obtain the colour dimension. `contract="dots"` requests
the same work plus canonical dot notation. `contract="selected"` requires
`representations=[...]` and permits additional connections only in those
representations; an empty list permits none. `contract="none"` performs only the
structural prerequisites of enabled identities. Representation filters never
suppress those prerequisites or undo construction's intrinsic normalization.
Trace collection does not evaluate the trace. Unrelated sums remain factored.
All modes use the same Rust domain planner and contraction implementation; dot
notation is scheduled with the other operations and obeys the same work budget.
Rust uses `AlgebraSettings.contract: AlgebraContraction` with `Fully` (default),
`Dots`, `Selected(Vec<LibraryRep>)`, and `None`. No nested contraction settings
are required.

`contract()` independently controls `metrics`, `rank_one`, `representations`,
`collect_chains`, `collect_traces`, and the work budget `max_steps_per_domain`. Representation filters select connections,
including a Lorentz connection on a gamma matrix whose spinor ports survive. `None`
permits every representation; `[]` permits none. Trace collection preserves an unevaluated
trace. Disabling it retains indexed factors or closed chains without creating trace notation.
Neither structural contraction nor color-delta contraction applies generator or Fierz identities.

Construction already normalizes a metric on equal compatible indices to its dimension,
substitutes compact vectors into metric ports, and orders the symmetric `dot` arguments.
Filters cannot suppress those intrinsic constructor rules. Structural contraction of
an indexed product uses the existing compact metric spelling `g(p(R), q(R))` for its
scalar product. `to_dots()` converts this surviving notation to `dot(p(R), q(R))`;
the two forms share the scalar-product display. Typed port composition can already
produce `dot` intrinsically. The conversion remains meaningful for surviving compact
metrics, and never searches repeated indices or contracts indexed products.
The nested spelling `p(q(R))` is not accepted by typed tensor admission: a registered
rank-one tensor requires its final structural port.

`undo_dots()` introduces compatible dummy indices for scalar products. `undo_trace()`
opens a trace as a closed indexed chain without evaluation; `undo_chain()` then exposes
its ordered factors. These operations preserve tensor interfaces, aliases, and unrelated
factorization. They unfold notation currently present and cannot recover an expression's
history after algebraic reduction. There is no aggregate `undo_all` operation.

Results report `reduction_status` as `Complete`, `Deferred`, or `Capped`.
Completion is relative to enabled capabilities and output requirements. Budget exhaustion
retains the exact unfinished tensor and its aliases; rerunning may make further progress.
A zero `max_passes` budget preserves eligible work and reports `Capped`.
Completed reruns are stable. Global arithmetic distribution remains an explicit `expand()` request.
Color dimension invariants are selected with
`simplify_algebra(color=True, color_substitute_cof_dimension_invariants=True)`.

The construction and indexing surface is intentionally retained: `flat`, `sigma`, and `f`
construct tensors; `rank` and `shape` inspect their interface; `index`, `reindex`, and
`rename_indices` assign or change indices; `permute_axes` changes logical axis order;
`compose`, `outer`, and `trace` build tensor operations. `format_tensor` and `to_latex`
control presentation. These are supported tensor operations even when a current application
has no call site. Arbitrary unchecked construction and post-hoc `reinfer` are not public
escape hatches: reconstruct through the checked `TensorExpression` constructor.

Idenso does not define a second parser syntax: the example constructs a Spenso-compatible
`TensorExpression` and then applies one Idenso transformation through its methods. Select the required identity families together; the shared planner schedules their prerequisites
while retaining unrelated factorization.

For implementation details, start with
#source-link("crates/idenso/src/lib.rs", label: "the Rust API") and
#source-link("crates/spynso3/src/expression.rs", label: "the Python binding").
]
