#import "../../shared.typ": boundary, callout, source-link

#let tutorial = [
= Tutorial

This tutorial parses and canonicalizes one vacuum integral in Rust. It deliberately disables
all evaluation backends, so the first success exercises Vakint's topology library and momentum
routing without requiring FORM, MATAD, FMFT, or pySecDec.

== Prerequisites

Use Rust 1.85 or newer. Clone the workspace, create a sibling binary project, and add the current
Vakint crate by path:

```sh
git clone https://github.com/alphal00p/gammaloop.git
cargo new vakint-first-match
cd vakint-first-match
cargo add vakint --path ../gammaloop/crates/vakint
```

Vakint uses Symbolica; ensure your use satisfies its license terms. External integration tools
are unnecessary for the matching-only settings below. For a reproducible checkout, use the
source revision shown in this page's footer.

#callout("Crates.io 0.1.2 uses the previous settings API", [
  The published `0.1.2` crate stores settings in `Vakint::new(Some(settings))` and calls
  `to_canonical(input.as_view(), true)`. This manual follows the current workspace API, where
  settings are passed explicitly to each operation; the versioned reference is authoritative for
  whichever source revision you install.
])

== Match a one-loop topology

Replace `src/main.rs` with:

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
    let vakint = Vakint::new()?;
    vakint.validate_settings(&settings)?;

    let input = vakint_parse!(
        "topo(prop(1,edge(1,1),k(1),muvsq,1))"
    )?;
    let canonical = vakint.to_canonical(
        &settings,
        input.as_view(),
        true,
    )?;

    println!("input: {}", VakintExpression::try_from(input)?);
    println!("canonical: {}", VakintExpression::try_from(canonical)?);
    Ok(())
}
```

Run `cargo run`. Success means Vakint recognizes the propagator as a supported one-loop vacuum
topology, prints a canonical short form, and exits without probing an external executable. The
canonical spelling may reorder identifiers or momentum routing; that normalized result, not
the input's labels, is the comparison boundary.

#callout("An empty evaluation order is intentional", [
  `EvaluationOrder::empty()` makes this a topology-matching tutorial. Default settings include
  analytic and numerical backends, and validating those settings checks their external
  dependencies even if you only planned to call `to_canonical`.
])

== Move from matching to evaluation

Once the canonical form is understood, choose one supported backend explicitly, validate its
dependencies, then apply the stages in order: tensor reduction, integral evaluation, and
numerical substitution. Record the normalization factor, epsilon depth, decimal precision,
backend settings, and tool versions with the result. Reuse one `Vakint` engine because topology
construction has a cost.

== Troubleshooting and next steps

- An unrecognized-integral error means the propagators, powers, masses, or momentum routing did
  not match a registered topology and `allow_unknown_integrals` is false. Print the long form
  and compare it with the generated topology table.
- A FORM or pySecDec version error means a nonempty evaluation order slipped into the settings;
  return to `EvaluationOrder::empty()` for pure canonicalization.
- If two inputs look different after matching, compare their canonical forms rather than their
  original propagator numbering.
- Continue with the topology-and-evaluation manual before enabling a backend or adding tensor
  numerators.

#boundary("Ground truth for this walkthrough", [
  The workflow follows #source-link("crates/vakint/README.md", label: "Vakint's maintained Rust example")
  and the matching path in #source-link("crates/vakint/src/lib.rs", label: "Vakint::to_canonical").
])
]
