#import "../../shared.typ": callout

#let quickstart-rust = [
= Using Idenso from Rust

This route constructs the same indexed metric contraction as the Python example, applies
Idenso's native `MetricSimplifier`, and checks the exact reduced Symbolica expression.

#callout("Use the current Git API until the next crate release", [
  The published `idenso 0.3.0` predates the current module layout and uses an older Symbolica
  dependency. Mixing its API with the current manual can create incompatible `Atom` types. The
  pinned Git dependency below matches this documentation; switch to the next published Idenso
  version once it carries the same API.
])

== Create a Rust project

Use Rust 1.85 or newer:

// docs-example: syntax
```sh
cargo new idenso-quickstart
cd idenso-quickstart
cargo add idenso \
  --git https://github.com/alphal00p/gammaloop.git \
  --rev 41e0eddb39c7c668074af483ac1d566e46c247a8
cargo add symbolica@2.2.0 --no-default-features
```

Replace `src/main.rs` with:

// docs-example: compile idenso-rust-quickstart
```rust
use idenso::shorthands::metric::MetricSimplifier;
use symbolica::parse_lit;

fn main() {
    idenso::representations::initialize();

    let expression = parse_lit!(
        spenso::g(spenso::mink(4, 0), spenso::mink(4, 1))
            * p(spenso::mink(4, 1))
    );
    let reduced = expression.simplify_metrics();

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
