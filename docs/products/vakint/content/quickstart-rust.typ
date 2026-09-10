#import "../../shared.typ": callout

#let quickstart-rust = [
= Using Vakint from Rust

This route parses and canonicalizes the same one-loop integral as the Python example while
keeping matching policy and expression ownership in Rust.

#callout("Use the current Git API until the next crate release", [
  The published `vakint 0.1.2` API differs from the current source while carrying the same version
  number. The pinned Git dependency below matches this manual. Switch to the next published
  Vakint version only after its Rustdoc exposes the `Vakint::new()` and explicit-settings API used
  here.
])

== Create a Rust project

Use Rust 1.85 or newer:

// docs-example: syntax
```sh
cargo new vakint-quickstart
cd vakint-quickstart
cargo add vakint \
  --git https://github.com/alphal00p/gammaloop.git \
  --rev 41e0eddb39c7c668074af483ac1d566e46c247a8
cargo add symbolica@2.2.0 --no-default-features
```

Replace `src/main.rs` with:

// docs-example: compile vakint-rust-quickstart
```rust
use vakint::{
    vakint_parse, EvaluationOrder, Vakint, VakintExpression, VakintSettings,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let settings = VakintSettings {
        allow_unknown_integrals: false,
        evaluation_order: EvaluationOrder::empty(),
        ..VakintSettings::default()
    };
    let engine = Vakint::new()?;
    engine.validate_settings(&settings)?;

    let input = vakint_parse!(
        "topo(prop(18,edge(7,7),k(99),muvsq,1))"
    )?;
    let canonical = engine.to_canonical(&settings, input.as_view(), true)?;
    let canonical = VakintExpression::try_from(canonical)?;

    assert!(canonical.to_string().contains("I1L"));
    println!("{canonical}");
    Ok(())
}
```

Run `cargo run`. Success means Vakint recognizes the input as its canonical one-loop topology and
exits without probing an evaluation executable. Symbolica's license terms apply to this native
route as well as to the Python community module.

Continue with the #link("guides/evaluation/")[matching and evaluation guide] before adding tensor
numerators or enabling analytic and numerical backends.
]
