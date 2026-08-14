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
data support an operation. Consult the Rust API reference when a method is unavailable for a
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

The proc-macro API reference lists the derive's helper attributes, allowed targets, and feature
conditions. Macro expansion happens at compile time and produces ordinary Rust
implementations.

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

For a released Python assembly that includes Spenso, install the provider's Symbolica
distribution with `python -m pip install symbolica`, then verify
`python -c "import symbolica.community.spenso"`. Module availability follows the Symbolica
assembly version; a local Spenso source checkout does not add the module to an installed wheel.
Source embedders must add `spynso3`
to the external #link("https://github.com/benruijl/symbolica-community")[symbolica-community]
assembly and invoke its `SymbolicaCommunityModule` registration while building that extension;
building the Rust crate alone does not inject the module into another Symbolica wheel.

```python
from symbolica.community.spenso import (
    Representation,
    Tensor,
    TensorIndices,
)

rep = Representation.euc(2)
indices = TensorIndices(rep("i"), rep("j"))
identity = Tensor.dense(indices, [1.0, 0.0, 0.0, 1.0])
```

The Python surface includes tensor and structure types, tensor networks, execution modes,
tensor libraries, expression evaluators, and symbolic-parallelism controls. Its signatures can
differ from the generic Rust API because `spynso3` provides the Python-specific conversions and
defaults.

Source starting points are #source-link("crates/spenso/src/spenso.rs", label: "the core crate"),
#source-link("crates/spenso-macros/src/lib.rs", label: "the derive crate"),
#source-link("crates/spenso-hep-lib/src/lib.rs", label: "the HEP library"), and
#source-link("crates/spynso3/src/lib.rs", label: "the Python adapter").
]
