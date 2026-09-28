#import "../../shared.typ": callout

#let quickstart-rust = [
= Using Idenso from Rust

This route constructs the same indexed metric contraction as the Python example, applies
Idenso’s shared `SymbolicTensor::contract`, and checks the exact reduced Symbolica expression.

#callout("Use the current Git API until the next crate release", [
  The published `idenso 0.3.0` predates the current module layout and uses an older Symbolica
  dependency. Mixing its API with the current manual can create incompatible `Atom` types. The
  example below uses a local checkout carrying this API and its pinned Symbolica revision.
  Use that checkout’s dependency versions when embedding it in another workspace.
])

== Create a Rust project

Use Rust 1.89 or newer:

// docs-example: syntax
```sh
cargo new idenso-quickstart
cd idenso-quickstart
cargo add idenso --path /path/to/gammaloop/crates/idenso
cargo add spenso --path /path/to/gammaloop/crates/spenso
```

Copy the checkout’s `[patch.crates-io]` dependency overrides into your `Cargo.toml`
so all tensor crates share one Symbolica build. Replace `src/main.rs` with:

// docs-example: compile idenso-rust-quickstart
```rust
use idenso::tensor::SymbolicTensor;
use spenso::{g, mink, vector};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    idenso::representations::initialize();

    let expression = g!(mink!(4, 0), mink!(4, 1)) * vector!(p, mink!(4, 1));
    let reduced = SymbolicTensor::infer(expression)?
        .contract(Default::default())?
        .resolved()?
        .into_expression();

    assert_eq!(reduced, vector!(p, mink!(4, 0)));
    println!("{reduced}");
    Ok(())
}
```

Run `cargo run`. Success means the assertion passes and the vector retains index `0` after the
metric contracts index `1`. Symbolica's license terms apply to the native crate as well as the
Python module.

Continue with the #link("guides/algebra/")[algebra guide] before composing multiple identity
families or multiplying expressions with independently allocated dummy indices.
]
