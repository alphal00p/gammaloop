# Symbolica external-constant merge reproducer

This independent crate uses only Symbolica and its own graphica/numerica crates.
It compares four rational expressions with builtin i and pi, evaluated
independently, built jointly, and merged from separate evaluators. Factorized
and expanded forms are checked against analytical numerical values. No
GammaLoop application or OEM-license initialization is involved.

The recorded before/after validation compiled this exact source against the
already built unpatched and patched Symbolica libraries: exit codes 101 and 0,
respectively. The standalone Cargo manifest and lockfile provide the independent
reproduction route below.

Run the unpatched pin from this directory:

```sh
cargo run --locked
```

It exits 101 because merged output 1 is wrong. The expected value is
−0.023271056693257727. The factorized merged value is −0.04334315115609011;
the expanded merged value is 0.011635528346628864 + 0.1096622711232151 i.
Independent and jointly built evaluators match the analytical oracle.

To test the patch, copy this entire directory to a scratch location and, from
that copy, run:

```sh
git clone https://github.com/alphal00p/symbolica symbolica-patched
git -C symbolica-patched checkout 4d0a833eb8e059d1f95bdae5abed2559830b235f
git -C symbolica-patched apply ../merge-constant-indices.patch
```

Replace the copy's three `[patch.crates-io]` entries with:

```toml
symbolica = { path = "symbolica-patched" }
graphica = { path = "symbolica-patched/lib/graphica" }
numerica = { path = "symbolica-patched/lib/numerica" }
```

Then `cargo run` updates dependency selection and passes the unchanged
expressions and assertions. Supply an authorized `SYMBOLICA_LICENSE` through
the environment if required; no key is included in this reproducer.

The patch keeps the incoming evaluator's external-constant descriptors intact
while scanning its constant table, remapping only the descriptor copied into
the destination evaluator. Mutating the incoming index early can incorrectly
classify a later literal as that external constant. Only
`src/evaluate/optimize.rs` changes.

The equivalent fix is published upstream as
[ba4958cb](https://github.com/symbolica-dev/symbolica/commit/ba4958cb43b9370bfe9230b7f541df8ed47b5968).
This small backport deliberately retains the existing 2.2 API; adopting that
newer upstream revision requires a separate dependency/API migration.
