#import "../../shared.typ": callout

#let tutorial = [
= Tutorial

This tutorial performs one concrete tensor contraction in Rust. The example makes the tensor
structure, dual index matching, storage, and numerical result visible—the same layers that a
larger Spenso network coordinates automatically.

== Prerequisites

Use Rust 1.85 or newer (Spenso uses Rust 2024 edition), then create a binary project:

// docs-example: syntax
```sh
cargo new spenso-first-contraction
cd spenso-first-contraction
cargo add spenso@0.6
```

This first contraction uses dense integer tensors and needs neither Symbolica nor the
`shadowing` feature.

== Contract one shared index

Replace `src/main.rs` with:

// docs-example: run
```rust
use spenso::{
    contraction::Contract,
    structure::{
        OrderedStructure, PermutedStructure,
        representation::{LibraryRep, Lorentz, RepName},
        slot::{DualSlotTo, IsAbstractSlot},
    },
    tensors::data::{DenseTensor, GetTensorData, SetTensorData},
};

fn main() {
    let rep = Lorentz {};
    let left_free = rep.new_slot(2, 0).to_lib();
    let right_free = rep.new_slot(2, 2).to_lib();
    let shared = rep.new_slot(2, 10).to_lib();

    let left_structure: OrderedStructure<LibraryRep> =
        PermutedStructure::from_iter([left_free, shared]).structure;
    let right_structure: OrderedStructure<LibraryRep> =
        PermutedStructure::from_iter([right_free, shared.dual()]).structure;

    let mut left = DenseTensor::<i32, _>::zero(left_structure);
    let mut right = DenseTensor::<i32, _>::zero(right_structure);
    left.set(&[0, 1], 2).unwrap();
    right.set(&[0, 1], 5).unwrap();

    let result = left.contract(&right).unwrap();
    assert_eq!(*result.get_ref([0, 0]).unwrap(), 10);
    println!("result structure: {}", result.structure);
    println!("result data: {:?}", result.data);
}
```

Run `cargo run`. Success means the assertion passes and the printed result contains `10` at
the component selected by the two remaining free indices. Spenso matched `shared` with
`shared.dual()`, summed that axis, and preserved `left_free` and `right_free` in the output
structure. The contraction follows index identity and duality, not merely equal dimensions.

#callout("Verification scope and cost", [
  The docs harness compiles and runs this Rust program; it syntax-checks the setup commands
  without creating a project or using the network. Success is the asserted component value
  `10` with two free output indices. The contraction is tiny; a clean external Cargo build can
  take minutes, while the program itself should finish in well under a second.
])

#callout("Structure is part of the value", [
  The data vectors alone do not say which axes may contract. Preserve the `OrderedStructure`
  with the tensor, and inspect the resulting structure before interpreting a flat component
  offset. A dimension match with unrelated abstract indices produces an exterior product, not
  the contraction shown above.
])

== From two tensors to a network

Pairwise `Contract` is the clearest first success. For a larger expression, place tensors in a
`Network`, inspect its graph, choose or benchmark a contraction strategy, and execute it. Use
`DataTensor` when individual nodes may be dense or sparse. Enable `shadowing` only when the
network must mix concrete data with Symbolica-backed parametric tensors.

== Troubleshooting and next steps

- If `contract` leaves both shared-looking axes free, print their slots and check that one is
  the dual of the other; dimensions alone are insufficient.
- If `set` rejects a component, compare the coordinate count and each coordinate with the
  structure's rank and dimensions.
- If a sparse contraction reports an empty input, distinguish an intentionally zero tensor
  from a tensor whose nonzero entries were never populated.
- Continue with the tensor-structures and networks manual before selecting an automatic
  contraction strategy.

]
