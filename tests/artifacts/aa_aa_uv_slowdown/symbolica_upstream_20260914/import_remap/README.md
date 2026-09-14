# Nonsymmetric function import check against upstream Symbolica

This standalone Cargo package uses only Symbolica, bincode and serde_json. The
Rust reproducer in `src/bin/import_remap.rs` is copied unchanged from the pinned
fork test. The function `f` has no symmetry attributes.

Dependency: [official upstream main at ba3737137c2a2ccd7bb39f0441837d38ec867e78](https://github.com/symbolica-dev/symbolica/tree/ba3737137c2a2ccd7bb39f0441837d38ec867e78),
observed on 2026-09-14. Its package version remains `2.2.0`; the exact commit
identifies the tested implementation. `Cargo.lock` pins all dependencies.

Set `SYMBOLICA_LICENSE` to a valid key for current upstream, then run from this
folder:

```sh
cargo run --locked --profile dev-optim --bin import_remap -- write tiny.symbolica
cargo run --locked --profile dev-optim --bin import_remap -- read tiny.symbolica
```

Use a new output path for another writer run: existing files are deliberately
not overwritten. Each command starts a separate process. The writer registers
`x,y,f` and exports `f(x,y)`; the reader registers `y,x,f` and requires exact
symbolic equality with `f(x,y)`. The expected ordered arguments are unchanged.

The old pinned fork at `4d0a833eb8e059d1f95bdae5abed2559830b235f` incorrectly
imported `f(y,x)`. Upstream now traverses the arguments of a function whose own
symbol ID did not change, in `AtomView::rename_no_norm`.

Build and runtime receipts are retained alongside the source. This small import
check is independent of the large GL262 evaluator-construction reproducer.

## Verified result on 2026-09-14

PASS with the original Rust source unchanged: a fresh writer exported `f(x,y)`;
a fresh reader imported `f(x,y)` and passed exact equality. `cargo check`, the
optimized native Rust build, and `cargo clippy -- -D warnings` all passed.
`run_r3/receipt.json`, `write.log` and `read.log` record the final run with the corrected S-format license. The full
binary identity is in `run_r1/receipt.json`; the immutable local binary is retained
at `run_r1/import_remap` but should not be committed or required for sharing.

Earlier attempts used a hash-delimited license key rejected by this upstream main:
[`validate_license_key`](https://github.com/symbolica-dev/symbolica/blob/ba3737137c2a2ccd7bb39f0441837d38ec867e78/src/license.rs#L435)
explicitly rejects keys containing `#`. The test nevertheless completes under
Symbolica's restricted single-thread mode. `run_r2` assigns the new key after the
local launcher, excluding inherited environment overrides as the explanation.
The license key itself is not recorded; its SHA256 identifies the supplied value.

The corrected S-format key supplied later by the user is accepted without the
license warning. `run_r3` repeats both fresh processes with that key and passes.
The environment variable is `SYMBOLICA_LICENSE` (with the final A), as read by
upstream. No license value is embedded in this package or its receipts.
