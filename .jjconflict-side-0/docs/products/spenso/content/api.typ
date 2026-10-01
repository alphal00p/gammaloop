#import "../../shared.typ": boundary, source-link

#let api = [
= Rust, macros, and Python APIs

== Core Rust package

The `spenso` crate organizes its public surface into `structure`, `tensors`, `contraction`,
`network`, `iterators`, and `algebra`. The important abstractions form a progression:

- `TensorStructure` and related traits describe slots, names, and contraction compatibility;
- `DenseTensor`, `SparseTensor`, and the heterogeneous tensor enums own storage;
- contraction traits perform pairwise or multi-tensor operations;
- network stores and libraries bind symbolic tensor names to concrete data and execute a graph.

Trait implementations and generic constraints determine which combinations of structure and
data support an operation. Consult the #link("reference/rust/")[Rust orientation] to choose the
relevant crate, then use its revision-specific Rustdoc when a method is unavailable for a
particular tensor type. The `shadowing` API and symbolic parallelism controls require their
corresponding Cargo features.

== Proc macros and HEP data

`spenso-macros` is a separate proc-macro crate. Its `SimpleRepresentation` derive generates the
representation and duality boilerplate used by Spenso index types. A declaration supplies a
symbolic name and chooses either `self_dual` or a `dual_name`:

```rust
use spenso_macros::SimpleRepresentation;
use spenso::structure::representation::RepName;

#[derive(SimpleRepresentation)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
#[representation(name = "flavor", dual_name = "AntiFlavor")]
struct Flavor {}
```

The #link("reference/rust/spenso_macros/derive.SimpleRepresentation.html")[`SimpleRepresentation`
Rustdoc] lists the derive's helper attributes and allowed targets. Macro expansion happens at
compile time and produces ordinary Rust implementations.

`spenso-hep-lib` supplies domain data and tensor-library construction for high-energy physics.
It is intentionally separate from the generic core. Users who need gamma matrices or physics
projectors add that package; generic Spenso users do not inherit those conventions implicitly.

== Python community module

#boundary("An adapter, not a Spenso wheel", [
  Python users import `symbolica.community.tensor`. The implementation is the `spynso3` Rust
  adapter and is distributed through the Symbolica community-module mechanism. Enabling
  Spenso's Rust `python` feature only enables conversion interoperability; it does not create a
  standalone importable `spenso` Python distribution.
])

Install the published Symbolica assembly with `python -m pip install --upgrade symbolica`, then verify
`python -c "import symbolica.community.tensor"`. Module availability follows the Symbolica
assembly version; a local Spenso source checkout does not add the module to an installed wheel.
Source embedders must add `spynso3`
to the external #link("https://github.com/symbolica-dev/symbolica-community")[symbolica-community]
assembly and invoke its `SymbolicaCommunityModule` registration while building that extension;
building the Rust crate alone does not inject the module into another Symbolica wheel.

The same module exports Idenso's algebra settings and exceptions. Its simplification
operations are `TensorExpression` methods, so symbolic pipelines can chain
`expression.simplify_algebra(gamma=True, color=True, epsilon=True)` while retaining the
tensor interface and literal aliases. Resolve the result with `to_expression()`;
`expand()` explicitly requests polynomial materialization.

Colour reduction accepts arbitrarily long traces in the supported fundamental
and adjoint representations, and structure-constant networks of arbitrary size.
Ordered generator traces are reduced to symmetric traces and
commutator terms; closed structure-constant cycles use the same trace reduction
in the adjoint representation. There is no fixed maximum trace length or cycle
size. The existing short-word identities remain useful shortcuts.

The result stays parametric in the group and representation. Higher symmetric
traces and their scalar contractions, such as `gram(rank, rep1, rep2)`, can remain
as exact invariants. A completed result therefore need not be a polynomial only
in `Ca`, `Cf`, and `Na`; higher-rank invariant relations are a separate algebraic
question. This also applies to repeated indices within a higher symmetric
trace. General networks of several symmetric invariants, words with multiple
symmetric blocks, and arbitrary matrix insertions can also remain exact; this
is not a complete canonical basis for every colour tensor. Disabling trace
evaluation preserves ordered trace notation.

Local colour identities can produce combinatorial term growth as their input
size increases. They use the shared reduction planner and preserve unrelated
sums and scalar factors; they do not expand the whole numerator. An explicit
work budget can return an exact capped result that can be resumed. Step limits
count kernel transformations; they do not impose a time limit within one identity.
Completion is relative to the enabled identities and requested output, while deferred
structural work retains its exact remaining expression.

`TensorNetwork` retains component data and library bindings while composing and executing a
graph. Its binary positional contraction is `contract_ports(rhs, left=..., right=...)`,
matching `TensorExpression.contract_ports`; `TensorExpression.contract()` instead performs
unary symbolic index contraction. Graph composition, indexing, and axis permutation remain
available because converting a network to a symbolic descriptor would lose stored components.
Execution and rendering use `execute`, `step`, `status`, `result_tensor`, `result_scalar`,
`expression`, `to_dot`, `to_html`, `render`, and `to_linnest`.

Index syntax is checked when a tensor structure is admitted. Subsequent
contraction, algebra, relabelling, and reconstruction use the admitted label as
an opaque identity: named and scoped labels have the same capabilities as plain
symbols. A label named `in`, `out`, or `mu_` remains literal index data; it does
not become a chain endpoint or pattern wildcard during algebra. Explicit pattern
APIs retain their pattern semantics. Numeric admission must preserve a label
exactly or reject it; it must never narrow distinct labels into the same index.
Representation compatibility, index incidence, and port order still govern
valid tensor operations. Component coordinates and representation dimensions
have their own input requirements; they are not abstract labels.

Component evaluation uses exact HEP data by default: four-dimensional Dirac
matrices in the Weyl basis, SU(3) generators and structure constants, and
dimension-dependent identity and metric factories. `TensorLibrary.hep_lib_atom()`
constructs an independent library with that same exact convention. Its entries
are Symbolica `Expression` values, including rational coefficients and algebraic
constants such as $sqrt(3)/2$. `TensorLibrary.hep_lib()` explicitly selects the
corresponding double-precision real or complex data. Selecting a component
library does not enable symbolic colour or gamma identities.

`Tensor.dense(descriptor, values)` preserves explicitly supplied Expressions;
it does not convert `E("1/3")` into a floating-point approximation. Integer-only
sequences are exact too. Sequences of floats or complex numbers retain numerical
storage. `Tensor.sparse(descriptor, Expression)` stores exact symbolic components
and accepts exact integer assignments. Libraries retain the component storage
chosen during registration.

```python
from symbolica import E
from symbolica.community.tensor import Representation, Tensor, TensorLibrary, TensorName

p = TensorName.vector("p")(Representation.mink(2))
library = TensorLibrary.hep_lib_atom()
library.register(Tensor.dense(p, [E("1/3"), E("2/7")]))
(p("mu") * p("mu")).to_tensor(library).scalar()  # 13/441 exactly
```

`Tensor`, `TensorNetwork`, and `TensorExpression` expose an immutable
`TensorStructure` through `.structure`. It is the canonical free-axis signature
of the whole expression treated as one opaque tensor. It contains `.axes`,
`.rank`, and `.shape`, without a tensor name, scalar arguments, or a layout
permutation. Thus `A(i,j)`, `A(j,i)`, and `A(i,j) + A(j,i)` have equal structures.
Construction sorts the supplied free axes by Spenso's representation/index order;
this is independent of function arguments, multiplication order, and summand order.

`structure.axes` is the canonical tuple of `Slot` and `Representation` objects.
`structure.slots()` returns a homogeneous list of indexed `Slot` objects, raising
`ValueError` at the first unresolved axis in canonical order.
`structure.representations()` returns a homogeneous list of unresolved
`Representation` objects, raising `ValueError` at the first indexed axis.
Neither method invents indices, discards labels, or filters out axes. Different
spaces and dimensions are allowed; scalars return empty lists. Unresolved axes
retain their multiplicity, with occurrence-local identifiers excluded from
signature equality. Arithmetic still checks positional compatibility of unresolved
axes before operating: equal opaque signatures do not imply equal wiring.

The expression retains its actual tensor calls and their argument order.
`expression.name` and `expression.arguments` describe its optional data identity.
`tensor.expression()` returns the symbolic descriptor, and `network.expression()`
returns its source computation. A signature alone cannot reconstruct a call:
`A(j,i)` has canonical free slots but still refers to the transposed occurrence.
`TensorExpression(raw, structure=signature)` validates the free axes without
reordering the expression; for a scalar zero, the signature supplies its missing rank.

For component coordinates and positional operations, use `value.axes` and
`value.shape`. Those describe the current view, following construction and
`value.permute_axes(...)`. Permuting explicitly indexed axes preserves the
canonical `.structure`. Permuting unresolved expression or network ports gives
them fresh dummy labels so their occurrences remain distinct; the resulting
signature includes those labels.
Network parsing and evaluation retain the argument-to-component correspondence;
there is no second mapping stored on the public signature.

The metadata displays separate presentation labels from exact index expressions.
`rep.name` is a `RepresentationName` with dimension-independent identity and
duality; its `metric_sign(i)` queries the canonical contraction sign. Dualizable
representations display a dual pairing. Structure displays show their canonical
free axes in selectable boxes; expression displays retain mathematical notation.

For indexed tensors, `expression.structure.slots()` exposes the external `Slot` objects.
The read-only `slot.representation` property returns their typed `Representation`,
including dimension and duality. Filter slots with `slot.representation == rep`
when constructing a projector for a particular representation.

Tensor multiplication contracts matching explicit indices together, including several
pairs in one product. Unresolved ports are contracted when their representations and
remaining labels determine a unique maximum pairing; competing pairings require an
explicit `contract` call. When explicit labels determine every contraction, their
pairing takes precedence over matrix-channel inference, including when a rank-two
expression is a product of vectors. Unresolved matrix channels retain ordered chain
composition. Symbolic contractions remain compact dots, chains, or brackets until their
evaluation is requested.

```python
from symbolica.community.tensor import (
    Representation,
    Tensor,
    TensorName,
)

rep = Representation.euc(2)
structure = TensorName("I")(rep("i"), rep("j"))
identity = Tensor.dense(structure, [1.0, 0.0, 0.0, 1.0])
```

The #link("guides/python/")[Python tensor-workflow guide] connects construction, libraries,
network execution, evaluators, and symbolic-parallelism policy. Use it with the structured
#link("reference/python/")[Python reference], whose exact signatures can differ from the generic
Rust API because `spynso3` provides Python-specific conversions and defaults.

Source starting points are #source-link("crates/spenso/src/spenso.rs", label: "the core crate"),
#source-link("crates/spenso-macros/src/lib.rs", label: "the derive crate"),
#source-link("crates/spenso-hep-lib/src/lib.rs", label: "the HEP library"), and
#source-link("crates/spynso3/src/lib.rs", label: "the Python adapter").
]
