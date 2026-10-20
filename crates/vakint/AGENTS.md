- If running into a linking issue complaining about missing `__emul...` symbol, then it may be a macos specific issue, that can be fixed by building with the following environment variable change `EXTRA_MACOS_LIBS_FOR_GNU_GCC=T` (see impact of this in file `build.rs`.
- A quick way to run all tests except pySecDec related ones (which can be slow) is the following command:

VAKINT_SKIP_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --package vakint --no-default-features

## Reduction benchmark maintenance

Always keep `docs/benchmark/benchmark.pdf`, its canonical Typst source
`docs/benchmark/benchmark.typ`, and `docs/benchmark/data/results.json` up to date
when reduction code, caches, external data, or relevant dependencies change.
The reproducible Rust collector is `examples/benchmark_reduction.rs`. From the
repository workspace root, refresh with:

```sh
cargo run --locked --release -p vakint --example benchmark_reduction -- \
  --form /path/to/form \
  --work-dir /path/to/workspace/vakint-benchmark \
  --output crates/vakint/docs/benchmark/data/results.json
typst compile crates/vakint/docs/benchmark/benchmark.typ \
  crates/vakint/docs/benchmark/benchmark.pdf
```

Use a single CPU affinity for comparable timings (for example, prefix the Cargo
command with `taskset -c CPU` on Linux); children inherit it. Keep build outputs
and temporary files under the development workspace. `--artifacts PATH`
optionally selects an existing external-data checkout or installation directory;
the default official `latest` loader downloads missing packs and reuses validated
local data. The collector resolves one data version for the whole run.

Preserve the exact-agreement and cache-hit checks, report the measured source,
compiler, FORM and data versions, and keep setup outside the reduction timer.
Never hand-edit timings or carry old measurements forward as new results.
Review the rendered one-page PDF and commit its source, measured JSON and PDF
together. Keep the large coloured comparison cells and the explicit
RustRed/FMFT ratio convention when changing the layout.
