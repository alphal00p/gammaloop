# Full GL262 replay against official Symbolica main

This standalone Cargo package consumes the complete six-file builder capture
from `gl262_faithful_20260914/GL262_capture_r1/capture/build_1602207_0`.
It has only Symbolica and serialization dependencies. No GammaLoop or Spenso.
The Rust main.rs loads the unchanged raw expression and ordered parameters,
restores the original symbol IDs, adapts the new default inlining policy in the
small FunctionMap, and certifies every preexisting definition field/raw Atom.
It then invokes the native evaluator `.build()` with the original H1/CPE5 settings.
No generated evaluator compilation, expression reduction or numerator expansion.

Dependency and checks: official main pinned at
[ba3737137c2a2ccd7bb39f0441837d38ec867e78](https://github.com/symbolica-dev/symbolica/tree/ba3737137c2a2ccd7bb39f0441837d38ec867e78).
Formatting, cargo check, optimized Rust build and clippy all pass; debug assertions
remain enabled. Both the original and main source declare package version 2.2.0.

Set `SYMBOLICA_LICENSE` to a valid S-format key, then:

```sh
cargo check --locked --profile dev-optim
cargo build --locked --profile dev-optim
target/dev-optim/symbolica-upstream-exact-builder --exact-builder /path/to/capture
```

The complete GL00 control builds and matches all four pinned operation counts.
With the corrected license, the full GL262 case reproduces the exact original
`advance out of bounds` panic in 186.136 s, sampled peak tree RSS 94.133 GB.
The 500 GB guard did not intervene. Final logs/receipts are at
`../GL00_exact_r3/` and `../GL262_exact_r3/`; `../immutable_exact_r1/` contains
frozen executable/source/lock identities. The code compiled before the final
license correction is identical; the key is a runtime setting.

The portable sharing folder is
`../../symbolica_gl262_reproducer_20260914`, containing both dependency revisions,
the complete input and small control, exact source, lockfiles and evidence.
See that folder's README.md and RESULTS.md for commands, stack and qualifications.
