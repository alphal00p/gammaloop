#import "../../shared.typ": callout

#let tutorial = [
= Tutorial

This tutorial parses and canonicalizes one vacuum integral in Rust. It deliberately disables
all evaluation backends, so the first success exercises Vakint's topology library and momentum
routing without requiring FORM, MATAD, FMFT, or pySecDec.

== Prerequisites

Use Rust 1.85 or newer. Create an independent binary project and pin the published Vakint release
and the direct Symbolica dependency required by its parsing macro:

// docs-example: syntax
```sh
cargo new vakint-first-match
cd vakint-first-match
cargo add vakint@0.1.2 symbolica@1.2
```

Vakint uses Symbolica; ensure your use satisfies its license terms. External integration tools
are unnecessary for the matching-only settings below. This example targets the crates.io
`vakint` 0.1.2 API rather than an unversioned repository checkout.

== Match a one-loop topology

Replace `src/main.rs` with:

// docs-example: syntax
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
    let vakint = Vakint::new(Some(settings))?;

    let input = vakint_parse!(
        "topo(prop(1,edge(1,1),k(1),muvsq,1))"
    )?;
    let canonical = vakint.to_canonical(input.as_view(), true)?;

    println!("input: {}", VakintExpression::try_from(input)?);
    println!("canonical: {}", VakintExpression::try_from(canonical)?);
    Ok(())
}
```

Run `cargo run`. Success means Vakint recognizes the propagator as a supported one-loop vacuum
topology, the `canonical:` line contains the short form `topo(I1L(muvsq,1))`, and the program
exits without probing an external executable. Other canonical spellings may reorder identifiers
or momentum routing; the normalized result, not the input's labels, is the comparison boundary.

#callout("Verification scope and cost", [
  The docs harness syntax-checks this release-targeted Rust source and the setup commands. A clean
  external project must compile and run it to verify that both printed forms are produced and the
  canonical form is recognized. Expect a dependency-heavy cold build to take minutes and the
  matching operation itself to take seconds. Symbolica may refuse a second concurrent unlicensed
  process; stop the other instance or use a licensed environment before retrying.
])

#callout("An empty evaluation order is intentional", [
  `EvaluationOrder::empty()` makes this a topology-matching tutorial. Default settings include
  analytic and numerical backends. Passing these settings into `Vakint::new` avoids selecting
  those backends when you only plan to call `to_canonical`.
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
  and compare it with the supported-topology reference.
- A FORM or pySecDec version error means a nonempty evaluation order slipped into the settings;
  return to `EvaluationOrder::empty()` for pure canonicalization.
- A Symbolica concurrency error means another unlicensed Symbolica process is active; stop it or
  retry in an environment licensed for concurrent instances.
- If two inputs look different after matching, compare their canonical forms rather than their
  original propagator numbering.
- Continue with the topology-and-evaluation manual before enabling a backend or adding tensor
  numerators.
]
