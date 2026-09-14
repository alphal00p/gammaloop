# Confirmed full GL262 standalone Symbolica failure — 2026-09-14

Both the live GammaLoop generation and the independent Rust replay fail in
Symbolica 2.2.0 at [4d0a833eb8e059d1f95bdae5abed2559830b235f](https://github.com/alphal00p/symbolica/tree/4d0a833eb8e059d1f95bdae5abed2559830b235f).
The replay has no GammaLoop or Spenso dependency and does not compile generated
evaluator code. Debug assertions are enabled; H1, CPE5, one core, direct translation.

```text
advance out of bounds: the len is 1721315294 but advancing by 10311249876
EvaluatorBuilder::build                evaluate/function_map.rs:343
AtomView::to_evaluator                 evaluate/tree.rs:258
optimize_horner_scheme_multiple        evaluate/tree.rs:411
horner_scheme                         collect.rs:1421
horner_scheme_impl                    collect.rs:1449–1453
normalize                             normalize.rs:1540
normalize (product iteration)         normalize.rs:662
ListIterator::next                    atom/representation.rs
bytes::Buf::advance                   panic
```

This detects an invalid child length during normalization of Horner output.
It does not yet identify the operation that first creates the invalid length.
Both runs use the same complete raw expression (19,982,016,710 bytes), ordered
220 parameters, function map, registry and settings. Identity registry and exact
raw FunctionMap certificates pass before the standalone build.
Expression SHA256: `f3f36f029ae973d75a8238b4654be4a30c8a1e1210d0441f9368bc8b7c31b3a6`.

| Run | Total wall time | Memory before builder | Peak memory | Exit |
|---|---:|---:|---:|---:|
| GammaLoop live, full generation | 3640.592 s | 46.674 GB RSS | 121.099 GB VmHWM | 101 |
| Standalone, complete input | 179.911 s | 19.988 GB RSS | 93.856 GB sampled tree RSS | 101 |

The standalone reaches `.build()` at 29.703 s. The remaining interval includes
panic reporting and process teardown; it is not a successful evaluator build time.
Both watchdogs enforce 500 GB and complete normally without intervention.
Live pre-builder peak VmHWM was 69.168 GB. Capture costs 7.479 s.

Local artifacts relative to `tests/artifacts/aa_aa_uv_slowdown`:

- `gl262_faithful_20260914/GL262_capture_r1/capture/build_1602207_0/`: all six input files.
- `gl262_faithful_20260914/GL262_standalone_r1/`: exact command, hashes, full backtrace and watchdog.
- `gl262_faithful_20260914/immutable_replay_r5/`: frozen validated Rust executable/source/lock.
- `gl262_faithful_20260914/GL262_capture_r1/`: live run logs, card, stage memory and receipts.

From the standalone package, with its pinned dependencies and configured license:

```sh
cargo check --locked --profile dev-optim
cargo build --locked --profile dev-optim
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export RUST_MIN_STACK=67108864 RUST_BACKTRACE=full
target/dev-optim/symbolica-evaluator-mre --exact-builder /path/to/build_1602207_0
```

The captured application callbacks are inventoried but cannot be serialized.
Native Symbolica special-function callbacks are restored by their native initializer.
This limitation does not prevent reproducing this failure; the GL00 physical control
also matches all four live evaluator operation counts.
