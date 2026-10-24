#import "../../shared.typ": callout

#let quickstart-rust = [
= Using Idenso from Rust

This route constructs the same indexed metric contraction as the Python example, applies
Idenso's native tensor contraction, and checks the exact reduced Symbolica expression.

#callout("Use the current Git API until the next crate release", [
  The published `idenso 0.3.0` predates the current module layout and uses an older Symbolica
  dependency. Mixing its API with the current manual can create incompatible `Atom` types. The
  pinned Git dependency below matches this documentation; switch to the next published Idenso
  version once it carries the same API.
])

== Create a Rust project

Use Rust 1.96 or newer:

// docs-example: syntax
```sh
cargo new idenso-quickstart
cd idenso-quickstart
```

Replace `Cargo.toml` with the following pinned dependencies. Keep these patches in your
project's root manifest: dependency manifests do not pass their patches to downstream projects.
They select one shared Symbolica kernel, matching numerical and graph primitives, and the
Spenso Typst bridge.

// docs-example: syntax
```toml
[package]
name = "idenso-quickstart"
version = "0.1.0"
edition = "2024"
rust-version = "1.96"

[dependencies]
idenso = { git = "https://github.com/alphal00p/gammaloop.git", rev = "8ef03d7fd7557772e1b531b71b71c08e3419ff54" }
symbolica = { version = "=3.0.1", default-features = false }

[patch.crates-io]
symbolica = { git = "https://github.com/symbolica-dev/symbolica", rev = "9cd57678d8817ec33c05cacba4925162f1c86b0b" }
numerica = { git = "https://github.com/symbolica-dev/symbolica", rev = "9cd57678d8817ec33c05cacba4925162f1c86b0b" }
graphica = { git = "https://github.com/symbolica-dev/symbolica", rev = "9cd57678d8817ec33c05cacba4925162f1c86b0b" }
symbolica-typst-plugin = { git = "https://github.com/symbolica-dev/symbolica-typst-plugin", rev = "57b7455b37c30f4f87353c0073e49f5225e937ab" }
```

Replace `src/main.rs` with:

Declare `p` as a rank-one tensor so native inference can check its representation and index.

// docs-example: compile idenso-rust-quickstart
```rust
use idenso::tensor::{ContractSettings, SymbolicTensor};
use symbolica::{parse_lit, symbol};

fn main() {
    idenso::representations::initialize();
    symbol!("p", tags = ["spenso::tensor", "spenso::rank1"]);

    let expression =
        parse_lit!(spenso::g(spenso::mink(4, 0), spenso::mink(4, 1)) * p(spenso::mink(4, 1)));
    let reduced = SymbolicTensor::infer(expression)
        .unwrap()
        .contract(ContractSettings::default().without_rank_one_tensors())
        .unwrap()
        .into_expression();

    assert_eq!(reduced, parse_lit!(p(spenso::mink(4, 0))));
    println!("{reduced}");
}
```

Run `cargo run`. Success means the assertion passes and the vector retains index `0` after the
metric contracts index `1`. Symbolica's license terms apply to the native crate as well as the
Python module.

Continue with the #link("guides/algebra/")[algebra guide] before composing multiple identity
families or multiplying expressions with independently allocated dummy indices.
]
