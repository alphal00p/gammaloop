#import "../../shared.typ": boundary, source-link

#let api = [
= Rust and Python APIs

== Rust package

The `idenso` crate exposes representation types and syntax macros together with several
rewrite families:

- `IndexTooling` covers canonicalization, conjugation, index wrapping, and dangling-index
  inspection for Symbolica atoms;
- `Cookable`, `CookSettings`, and the cook filters control reversible or flattening encodings;
- `SymbolicTensor::collect` collects selected sectors while retaining factorization;
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
- arithmetic expansion: `expand` explicitly materializes an expanded result;
- index preparation: typed `wrap_indices`, `list_dangling`, and the explicit
  `intern="indices"` or `intern="flattened"` construction/indexing policy;
- reduction: `contract(...)` and `simplify_algebra(...)`;
- notation: `to_dots`, `undo_dots`, `undo_chain`, and `undo_trace`;
- canonical labels: `canonize`;
- conjugation: `dirac_adjoint`.

== Choose a reduction operation
<tensor-reduction-guide>

Use `contract()` for structural work: metrics and identity tensors, dimension
factors, vector substitutions, scalar products, ordered matrix chains and
unevaluated traces. Use `simplify_algebra()` when the calculation also requires
Dirac, colour or epsilon identities. Both return a `TensorExpression`, preserve
its logical external ports and share the Rust planner and contraction engine.

Keep the expression typed throughout the calculation. Construct tensor names and
representations once, label their ports, then call the operation on the product.
The Python examples below run in order in one session:

```python
from symbolica.community.tensor import (
    ReductionStatus, Representation, TensorExpression, TensorName,
)

lorentz = Representation.mink(4)
metric = TensorExpression.g(lorentz)
p = TensorName.vector("p")(lorentz)
expression = metric("mu", "nu") * p("mu")

reduced = expression.contract()
assert reduced == p("nu")
assert reduced.reduction_status == ReductionStatus.Complete
```

The constructors check tensor structure. Index labels identify axes, rather than
numerical components; after admission they need not be simple symbols. Matching
labels and representation names alone are insufficient: dimensions, duality,
orientation and dummy scope must also be compatible. Use `contract_ports(rhs,
left=..., right=...)` when explicitly connecting selected logical ports of two
tensors. The unary `contract()` instead discovers permitted connections in an
expression.

Use `to_expression()` when leaving the tensor API for ordinary Symbolica algebra.
Do not insert raw-expression round trips between tensor operations. To admit
external raw syntax, use the checked `TensorExpression(raw)` constructor.

=== Structural contraction settings

All arguments are keyword-only. Their defaults permit all supported structural
work, with no work limit:

#table(
  columns: (1.4fr, .8fr, 3fr),
  [Keyword], [Default], [Meaning],
  [`representations`], [`None`], [Allowed connecting representation families; `None` permits all and `[]` permits none.],
  [`metrics`], [`True`], [Substitute metric/identity ports and form dimension factors.],
  [`rank_one`], [`True`], [Substitute vectors into compatible tensor or vector ports, including gamma slashes and scalar products.],
  [`collect_chains`], [`True`], [Collect connected matrix factors with their order and orientation intact.],
  [`collect_traces`], [`True`], [Represent closed compatible chains as traces without evaluating them.],
  [`expand`], [`True`], [Permit distribution needed by a selected contraction. `False` retains products of independent sum alternatives; it still permits boundary substitutions and recursion into existing branches.],
  [`order`], [`None`], [Optional zero-based permutation of normalized top-level factors for the initial contraction; normally let the planner choose.],
  [`max_steps_per_domain`], [`None`], [Optional limit on successful planner transformations; see completion below.],
)

`expand=True` permits necessary local distribution; it does not request global
polynomial expansion. For metrics only, use
`expression.contract(rank_one=False, collect_chains=False, collect_traces=False)`.
With `collect_traces=False`, chain collection does not create traces indirectly;
an exact indexed or closed-chain form remains. Collecting a trace never evaluates
it. Contracting a colour delta never enables generator or Fierz identities.

Filters apply to connections, not whole factors. A Lorentz connection can be
contracted into a gamma matrix while its spinor ports remain external. A filter
accepts `Representation` or `RepresentationName` values and selects the base
family, including its dual: `Representation.mink(4)` also permits compatible
symbolic-D Minkowski connections. It does not make unequal dimensions compatible.
Chain assembly and trace closure need permission for their connecting family.

```python
from symbolica import S

D = S("D")
lorentz_D = Representation.mink(D)
spin = Representation.bis(4)
gamma = TensorExpression.dirac_gamma(D)
p_D = TensorName.vector("p_D")(lorentz_D)
mixed = gamma("a", "b", "mu") * p_D("mu")

slashed = mixed.contract(representations=[Representation.mink(4)])
assert slashed.rank == 2
assert slashed == mixed.contract(representations=[lorentz_D.name])
assert expression.contract(representations=[]) == expression
```

The remaining ports in `slashed` are the original ordered spinor ports $a,b$.
`None` and `[]` are different permissions, rather than two spellings of a default.

=== Algebra identities and contraction modes

Python enables `gamma=True`, `color=True` and `epsilon=False` by default.
For a gamma-only calculation, write `simplify_algebra(color=False)`; for colour
only, write `simplify_algebra(gamma=False)`. Epsilon-only reduction requires
`simplify_algebra(gamma=False, color=False, epsilon=True)`.

The `contract` keyword controls structural work in addition to prerequisites of
the enabled identities. These are the supported modes:

#table(
  columns: (1fr, 3.5fr),
  [Mode], [Requested structural work],
  [`"fully"`], [Default: all supported contractions, chain collection and trace collection.],
  [`"dots"`], [The same work as `"fully"`, plus canonical index-free dot notation.],
  [`"selected"`], [Additional connections in `representations=[...]`, using the same base-family filter as `contract()`. An empty list allows no additional connections.],
  [`"minimal"`], [Structural boundary substitutions and work inside existing branches, without distributing products of independent sum alternatives. Enabled identity prerequisites remain available.],
  [`"none"`], [Only connections required by enabled identities; no additional structural contraction.],
)

`representations` is required for `"selected"` and is rejected for the other
modes. The filter never suppresses an identity's prerequisites. For example,
the contracted Clifford pair reduces even with no additional structural work:

```python
pair = gamma("a", "b", "mu") * gamma("b", "c", "mu")
identity = D * spin.g("a", "c")

assert pair.simplify_algebra(color=False, contract="none") == identity
assert pair.simplify_algebra(
    color=False, contract="selected", representations=[]
) == identity
assert pair.simplify_algebra(
    gamma=False, color=False, contract="none"
) == pair
```

Necessary Lorentz and spinor contractions are part of the selected gamma
identity. This permission does not extend to every identity in a connected
component. In particular, a gamma rewrite that emits epsilon does not enable
epsilon identities. The default `"fully"` contracts closed colour deltas even
with `color=False`, so no trailing `contract()` is needed to obtain their
representation dimension. Python's separate `color=True` default is what
authorizes colour generator identities.

`"minimal"` requests structural work without forming cross-products of separate
sum alternatives. It is more permissive than `"none"`: a metric can still relabel
the exposed ports of a sum, and work inside its existing branches can finish.
Standalone `contract(expand=False)` uses the same structural policy. This does
not restrict the contractions required by enabled algebra identities:

```python
q = TensorName.vector("q")(lorentz)
A = TensorName("A")(lorentz, lorentz)
product = (metric("mu", "nu") + A("mu", "nu")) * (p("mu") + q("mu"))
minimal = product.simplify_algebra(
    gamma=False, color=False, contract="minimal"
)
assert minimal == product
assert minimal.reduction_status == ReductionStatus.Complete
assert product.contract(expand=False) == product
assert product.contract().to_expression() != product.to_expression()

boundary = metric("mu", "nu") * (p("mu") + q("mu"))
assert boundary.contract(expand=False) == p("nu") + q("nu")
```

Here `Complete` means no work remains under the chosen distribution permission.
Use `contract()` or `contract="fully"` when the remaining sum-product connections
should also be contracted. Minimal structural work is not a promise that enabled
gamma or colour identities cannot introduce sums.

`collect_coefficients=True` collects exact coefficients produced in contracted
Dirac traces and their bound vector factors, combining equal terms and exposing
cancellations. It includes Symbolica's `expand_num()` numerical distribution
within those generated coefficients, such as $2(x+y) = 2x + 2y$.
Set `collect_coefficients=False` to retain the nested output of
those identities. This option does not expand unrelated input sums, scalar
prefactors, or free Dirac traces. Other identity kernels retain their existing
coefficient normalization; disabling collection cannot undo earlier reduction.
For example, `numerator.simplify_algebra(collect_coefficients=False)` requests
the nested form without changing the enabled identity families or contraction
mode. Global arithmetic expansion remains an explicit `expand()` call.

=== Gamma output and conventions

Optional family keywords default to `None`, meaning the native family default
shown below. Supplying any gamma option while `gamma=False` is an error, even
when supplying a value equal to its default.

#table(
  columns: (2.3fr, 3.8fr),
  [Keyword], [Default and effect],
  [`gamma_output`], [Default `"reduced"` applies Clifford identities; `"chains"` collects Dirac words without Clifford reduction or trace evaluation.],
  [`gamma_ordering`], [Default `"repeated_pairs"` moves repeated gammas toward their partners. `"canonical"` orders open chains through Clifford swaps and may introduce sums.],
  [`gamma_evaluate_traces`], [Default `True`: evaluates supported closed Dirac traces in reduced output; inactive for `"chains"`.],
  [`gamma0`], [Default `False`: enables additional gamma-zero preparation and factoring. It is not an off switch for every supported gamma-zero identity.],
  [`gamma_conjugate`], [Default `False`: prepares explicitly conjugated Dirac matrices for reduction; does not conjugate the whole expression.],
  [`gamma_expand_three_gamma_epsilon`], [Default `False`: allows the four-dimensional three-gamma expansion in the gamma5/epsilon basis, without enabling epsilon identities.],
)

`gamma0` and `gamma_conjugate` preparation can run with chain output. Ordering and
the three-gamma epsilon expansion act on reduced output. The ordinary Clifford
algebra and traces support a symbolic Lorentz dimension. Gamma-five, axial,
Chisholm and three-gamma epsilon rules require their supported four-dimensional
conventions; the operation does not invent a dimensional-regularization scheme
for gamma-five. Spinor trace normalization belongs to the spinor representation,
independently of the Lorentz dimension. `dirac_gamma(D)` uses spinor dimension four.

Trace collection and evaluation are separate requests:

```python
closed = gamma("a", "b", "mu") * gamma("b", "a", "nu")
collected = closed.contract()
evaluated = closed.simplify_algebra(color=False)

assert collected != evaluated
assert evaluated == 4 * lorentz_D.g("mu", "nu")
assert closed.simplify_algebra(
    color=False, gamma_evaluate_traces=False
) == collected
```

=== Colour output and conventions

The three colour options require `color=True`. `color_evaluate_traces=True`
evaluates supported closed colour traces; disabling it retains traces, though
enabled Fierz operations can still join them. `color_expand_fierz=True` enables
cross-chain Fierz identities. Disabling it does not disable other colour
identities or guarantee a sum-free result.

`color_substitute_cof_dimension_invariants=False` retains symbolic group
invariants. Set it to `True` to express supported SU($N$) invariants for recognized
fundamental or adjoint dimensions using $T_R=1/2$, $C_F=(N^2-1)/(2N)$ and $C_A=N$.
It is not a blanket SU(3) substitution: the representation dimensions choose $N$.
Compare independent results using the same invariant convention.

Closed products of structure constants can require quartic invariants. The
scalar `gram(4,coad(N_A),coad(N_A))` denotes the full contraction
$d_A^(a b c d) d_A^(a b c d)$, without division by $N_A$. The symmetric trace
defining each $d_A$ uses normalized symmetrization. Retaining this invariant is a
completed colour reduction, not an unresolved tensor contraction. With explicit dimension
substitution, recognized adjoint dimensions $N_A=N^2-1$ give
$d_A^(a b c d) d_A^(a b c d)=N^2(N^2+36)(N^2-1)/24$.
For example, the closed colour factor of the four-loop ghost graph FK0032 is
$N_A C_A^4/24 - d_A^(a b c d) d_A^(a b c d)$; for SU(3) the two contributions are
$27$ and $135$, giving $-108$. Its independent Lorentz momentum factors remain
factorized under both output conventions.

```python
fundamental = Representation.cof(3)
color_identity = TensorExpression.g(fundamental, fundamental.dual())
generator = TensorExpression.color_t(8, 3)
casimir = generator("a", "i", "j") * generator("a", "j", "k")

symbolic = casimir.simplify_algebra(gamma=False)
su3 = casimir.simplify_algebra(
    gamma=False, color_substitute_cof_dimension_invariants=True
)
assert su3 == color_identity("i", "k") * 4 / 3

delta_loop = color_identity("i", "j") * color_identity("j", "i")
assert delta_loop.contract() == TensorExpression(3)
assert delta_loop.simplify_algebra(gamma=False, color=False) == TensorExpression(3)
```

Here `symbolic` retains $C_F$; `su3` substitutes $4/3$. The delta loop requires
only structural contraction. Fundamental and dual ports are written explicitly
so their orientation is preserved.

=== Factorization, completion and work limits

Reduction works outside-in: the planner starts from a shallow view of factors
and their external interfaces, then opens a selected region when an authorized
contraction or identity needs its interior. It also processes internal work in
scalar-interface leaves. An unrelated sum or scalar coefficient stays factored.
A selected operation can require local distribution, for example when contracting
into a tensor sum; this is different from expanding the entire numerator.
Use `expand()` for an explicitly requested polynomial output.

Choose the required families and output form together. The planner can contract
before an identity, process its result, and revisit newly eligible work without
a user-written sequence of whole-expression passes. Settings specify permitted
work and output; `order` is an optional structural hint, not a family pipeline.
Reusable Python settings are ordinary dictionaries passed with `**options`.

`max_steps_per_domain=None` is unlimited for both operations. There are no
automatic term-count, workspace or recursion-depth work caps. An explicit
nonnegative integer limits successful planner transformations, not terms,
wall-clock time or individual inner-kernel steps. Inspect `reduction_status`:

- `Complete`: no eligible work remains under the requested identities and output.
  Disabled identities can still be present; this is not global canonicality.
- `Deferred`: the exact result contains eligible work not performed, for example
  an unsupported case or a detected rewrite cycle. Repeating identical settings
  need not resolve it.
- `Capped`: the caller's explicit budget stopped work. The entire unfinished
  expression remains available for a subsequent call.

```python
capped = closed.simplify_algebra(color=False, max_steps_per_domain=0)
assert capped.reduction_status == ReductionStatus.Capped
assert capped == closed

finished = capped.simplify_algebra(color=False)
assert finished.reduction_status == ReductionStatus.Complete
assert finished.simplify_algebra(color=False) == finished
```

Zero reports `Capped` only when eligible work remains; an already complete input
can report `Complete`. Completed reruns with the same settings are stable.
`contraction_complete` specifically reports certified structural completion,
rather than completion of every algebra family. No terms are discarded when work
is unfinished.

=== Explicit notation conversions

The notation methods preserve the represented tensor and its external ports.
They do not evaluate traces or other identities, and retain unrelated sums and
scalar factorization:

- `to_dots()` changes surviving index-free scalar-product notation to canonical
  dots. It does not search for repeated indices or contract indexed products.
- `undo_dots()` unfolds scalar products into vector components with fresh,
  compatible dummy indices.
- `undo_trace()` unfolds a trace into a closed indexed chain without evaluation.
- `undo_chain()` exposes the chain's ordered factors, preserving external ports
  and orientation. Apply it after `undo_trace()` to expose a trace's factors.

```python
q = TensorName.vector("q")(lorentz)
indexed = p("mu") * q("mu")
assert indexed.to_dots() == indexed

dots = indexed.simplify_algebra(gamma=False, color=False, contract="dots")
assert dots == indexed.contract().to_dots()
assert dots.undo_dots().contract().to_dots() == dots

explicit = collected.undo_trace().undo_chain()
assert explicit.contract() == collected
```

Structural contraction can leave the compact metric spelling `g(p(R), q(R))`;
`to_dots()` converts it to `dot(p(R), q(R))`. Both share the scalar-product display.
Construction already performs intrinsic normalization, including compatible
metric self-traces, vector substitutions in compact metrics and symmetric dot
ordering. Typed composition may already produce a dot. Filters cannot suppress
or undo those constructor rules. The nested spelling `p(q(R))` is not supported
typed syntax; a registered vector needs its structural port.

The undo methods allocate collision-free scoped indices, and their explicit
results survive construction. They unfold the notation currently present;
they cannot recover historical factors after algebraic evaluation. There is no
aggregate `undo_all`. Construction, index renaming, arithmetic expansion,
canonicalization, inspection and tensor conjugation remain separate operations.

=== Rust settings and defaults

Rust calls the same operations with `ContractSettings` and `AlgebraSettings`.
`ContractSettings::default()` enables structural work without filters or a limit.
`AlgebraSettings::default()` instead has `gamma=None`, `color=None`,
`epsilon=false`, `contract=AlgebraContraction::Fully`, `collect_coefficients=true`
and `max_passes=None`.
This differs from Python's enabled gamma/colour default. Enable a native family
with `Some(GammaSimplifySettings::default())` or
`Some(ColorSimplifySettings::default())`; `AlgebraSettings::hep()` enables all
three families, including epsilon.

`AlgebraContraction` has `Fully`, `Dots`, `Selected(Vec<LibraryRep>)`, `Minimal`
and `None`.
Standalone contraction uses `representations: Option<&[LibraryRep]>`, where
`None` is unrestricted and `Some(&[])` permits no connections. Both native
settings use `max_passes: Option<usize>` for the opt-in budget.
`ContractSettings.expand=false` requests the same structural distribution policy
as `AlgebraContraction::Minimal`.
`AlgebraSettings` has no nested `ContractSettings`: it enables algebraic
capabilities on the same planner. See the generated references for complete
signatures and native family types.

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
