#import "../../shared.typ": callout

#let tutorial = [
= Tutorial

This tutorial parses and canonicalizes one vacuum integral in Rust. It deliberately disables
all evaluation backends, so the first success exercises Vakint's topology library and momentum
routing without requiring FORM, MATAD, FMFT, or pySecDec.

== Prerequisites

Use Rust 1.85 or newer and work from the GammaLoop source revision identified in this page's
footer. The `latest` documentation follows that checkout's API; use an immutable documentation
snapshot when working from a released package. Enter the repository development environment:

// docs-example: syntax
```sh
nix develop
cargo test --locked -p alphal00p-docs-examples
```

Vakint uses Symbolica; ensure your use satisfies its license terms. External integration tools
are unnecessary for the matching-only settings below. The documentation test compiles this
example against the same workspace revision as the generated Rust reference.

== Match a one-loop topology

The following is a complete program. The documentation test compiles it directly; to inspect
its output, place the same program in a scratch binary target in this checkout and run that
target.

// docs-example: compile vakint-first-match
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
    let canonical = vakint.to_canonical(&settings, input.as_view(), true)?;

    println!("input: {}", VakintExpression::try_from(input)?);
    println!("canonical: {}", VakintExpression::try_from(canonical)?);
    Ok(())
}
```

Runtime success means Vakint recognizes the propagator as a supported one-loop vacuum topology,
the `canonical:` line contains the short form `topo(I1L(muvsq,1))`, and the program exits without
probing an external executable. Other canonical spellings may reorder identifiers or momentum
routing; the normalized result, not the input's labels, is the comparison boundary.

#callout("Verification scope and cost", [
  The docs harness compiles this source against the current workspace and syntax-checks the setup
  command. Run the program in a provisioned checkout to verify that both printed forms are
  produced and the canonical form is recognized. Expect a dependency-heavy cold build to take
  minutes and the matching operation itself to take seconds. Symbolica may refuse a second
  concurrent unlicensed process; stop the other instance or use a licensed environment before
  retrying.
])

#callout("An empty evaluation order is intentional", [
  `EvaluationOrder::empty()` makes this a topology-matching tutorial. Default settings include
  analytic and numerical backends. Validating the explicit empty order before `to_canonical`
  prevents an accidental external-backend probe in a matching-only workflow.
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
- Continue with the topology-and-evaluation guide before enabling a backend or adding tensor
  numerators.
]
