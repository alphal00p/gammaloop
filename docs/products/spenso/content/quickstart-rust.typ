#import "../../shared.typ": callout

#let quickstart-rust = [
= Using Spenso from Rust

This example constructs two dense tensors, contracts one pair of dual indices, and checks one
output component. It uses a two-dimensional Euclidean representation as a compact software
exercise rather than a Lorentz-space physics example.

== Create a Rust project

Use Rust 1.85 or newer:

// docs-example: syntax
```sh
cargo new spenso-quickstart
cd spenso-quickstart
cargo add spenso@0.6.0
```

Replace `src/main.rs` with:

// docs-example: run spenso-rust-quickstart
```rust
use spenso::{
    contraction::Contract,
    structure::{
        OrderedStructure, PermutedStructure,
        representation::{Euclidean, LibraryRep, RepName},
        slot::{DualSlotTo, IsAbstractSlot},
    },
    tensors::data::{DenseTensor, GetTensorData, SetTensorData},
};

fn main() {
    let rep = Euclidean {};
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

Run `cargo run`. Success means the assertion passes and the printed result contains `10` at the
component selected by the two remaining free indices. Spenso matches index identity and duality,
not merely equal dimensions.

#callout("The docs run this example", [
  The documentation test compiles and executes this program against the current workspace. The
  contraction itself needs neither Symbolica nor Spenso's optional `shadowing` feature.
])

Continue with the #link("tutorial/")[first tensor workflow] for a more deliberate explanation of
structure and storage, then use the #link("guides/networks/")[tensor-network guide] when more than
two tensors should be planned and executed together.
]
