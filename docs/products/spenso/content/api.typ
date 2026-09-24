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
  Python users import `symbolica.community.spenso`. The implementation is the `spynso3` Rust
  adapter and is distributed through the Symbolica community-module mechanism. Enabling
  Spenso's Rust `python` feature only enables conversion interoperability; it does not create a
  standalone importable `spenso` Python distribution.
])

Install the published Symbolica assembly with `python -m pip install --upgrade symbolica`, then verify
`python -c "import symbolica.community.spenso"`. Module availability follows the Symbolica
assembly version; a local Spenso source checkout does not add the module to an installed wheel.
Source embedders must add `spynso3`
to the external #link("https://github.com/symbolica-dev/symbolica-community")[symbolica-community]
assembly and invoke its `SymbolicaCommunityModule` registration while building that extension;
building the Rust crate alone does not inject the module into another Symbolica wheel.

The same module exports Idenso's algebra settings and exceptions. Its simplification
operations are `TensorExpression` methods, so symbolic pipelines can chain
`expression.simplify_gamma().simplify_color().simplify_metrics()` while retaining the
tensor interface.

`Tensor`, `TensorNetwork`, and `TensorExpression` expose an immutable
`TensorStructure` through the `.structure` property. It records an optional
`TensorName`, scalar `.arguments`, and the ordered `.slots`, `.rank`, and `.shape`.
Shape entries are Python integers or symbolic expressions. Unresolved ports stay
as `Representation` objects; assigning indices produces new metadata rather than
mutating an earlier structure snapshot.

`tensor.expression()` returns its symbolic descriptor, independently of stored
component values. `network.expression()` returns its source computation, which
execution does not replace with the result. These replace the former expression-returning
`structure()` methods; the former `.interface` tuple is now `.structure.slots`.

The metadata displays separate presentation labels from exact index expressions.
`rep.name` is a `RepresentationName` with dimension-independent identity and
duality; its `metric_sign(i)` queries the canonical contraction sign. Dualizable
representations display a dual pairing rather than claiming a metric signature on
one space. Structure displays put the named head and ordered ports in selectable
boxes. Ordinary `TensorExpression` outputs retain their mathematical rendering.

For indexed tensors, `expression.structure.slots` exposes the external `Slot` objects.
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
from symbolica.community.spenso import (
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
