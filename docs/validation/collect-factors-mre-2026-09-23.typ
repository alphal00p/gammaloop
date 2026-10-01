= Standalone Symbolica factor-collection stress reproducer

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.

#link("../../examples/collect_factors_mre.rs")[`examples/collect_factors_mre.rs`]
is a single executable Rust script with an embedded Cargo manifest pinning
Symbolica to `=3.0.0`. It needs no GammaLoop code, input files or workspace.
Supply the usual `SYMBOLICA_LICENSE` environment variable.

```sh
rust-script examples/collect_factors_mre.rs
rust-script examples/collect_factors_mre.rs 2048 128
rust-script examples/collect_factors_mre.rs 4096 256
```

Arguments are outer summand count and inner-sum width; defaults are 1024 and 64.
Each invocation builds `sum_i c(i) * (sum_j a(i,j))`, prints input bytes and RSS,
calls `collect_factors()` exactly once, reports elapsed collection time and
Linux peak RSS, and asserts structural equality with the original expression.
There is no `.expand()` or repeated factorization. Other platforms can use an
external memory profiler instead of Linux's `/proc/self/status`.

Each inner sum becomes an owned factor key. At the outer sum, the collector
fills absent keys into every term's factor map even when the exponent is zero.
There are `2*n*n` entries after that phase, plus repeated copies of the inner
sums. Increasing the width amplifies the copied-key cost independently of the
outer term count. Input size grows linearly in each parameter.

#table(
  columns: (auto, auto, auto, auto, auto),
  [Terms], [Width], [Input MiB], [Peak RSS GiB], [Result],
  [256], [64], [0.224], [0.090], [Unchanged],
  [512], [64], [0.464], [0.310], [Unchanged],
  [1024], [2], [0.050], [0.264], [Unchanged],
  [1024], [64], [0.944], [1.219], [Unchanged],
  [2048], [128], [3.763], [8.847], [Unchanged in 11.188 s],
  [4096], [256], [15.031], [12.027], [Stopped at 12-GiB bound],
)

These are fresh release-mode processes on a shared Linux host. Collection times
exclude expression construction and compilation; peak RSS includes construction.
The measurements use an external monitor with 20-ms polling and a 12-GiB RSS /
300-second wall bound. The standalone script itself has no resource cap: its
`4096 256` case deliberately continues allocating unless externally stopped.
The final row is a lower bound on memory required to finish, not a completed run.

The adjacent evidence directory records the source/toolchain manifest, raw
measurement outputs, and command-line checks. Build and execution logs plus the
standalone generated Cargo package and lockfile remain under
`/tmp/collect-factors-rustscript-2026-09-23/`.
