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
taskset -c CPU cargo run --locked --release -p vakint --example benchmark_reduction \
  --features feynkit-ingress,experimental-rustred,weighted-lowering-diagnostics,symbolica_community_module -- \
  --form /path/to/form \
  --work-dir /path/to/workspace/vakint-benchmark \
  --output crates/vakint/docs/benchmark/data/results.json \
  --build-profile release --select-current
typst compile crates/vakint/docs/benchmark/benchmark.typ \
  crates/vakint/docs/benchmark/benchmark.pdf
```

Replace `CPU` with one permitted CPU on Linux, or use the equivalent affinity
tool on another platform; children inherit it. Keep build outputs
and temporary files under the development workspace. `--artifacts PATH`
optionally selects an existing external-data checkout or installation directory;
the default official `latest` loader downloads missing packs and reuses validated
local data. The collector resolves one data version for the whole run.

Completed, validated runs are appended under `reduction_runs`; original measurements,
ingress records, raw samples and protocol metadata remain unchanged. `--select-current`
selects the new run for the report. Omit it to record a run while retaining the
existing selection. An initially empty output selects its first run. The selected
identifier is recorded explicitly as `selected_reduction_run` in the JSON.
Declare the actual Cargo profile with `--build-profile`; this label is not
compiler introspection. Final speed comparisons require the documented `--release`
build without profile overrides. Diagnostic `ci-optim` runs must declare that
profile and retain their distinct provenance. Optimization/LTO settings are never
inferred from the absence of debug assertions.

Preserve the exact-agreement and cache-hit checks, report the measured source,
compiler, FORM and data versions, and keep setup outside the reduction timer.
Never hand-edit timings or carry old measurements forward as new results.
Review the rendered PDF and commit its source, measured JSON and PDF
together. Keep the large coloured comparison cells and the explicit
RustRed/FMFT ratio convention when changing the layout.

For five-loop report integration, follow
[the benchmark maintenance contract](docs/benchmark/maintenance.typ). Use the
existing report adapter to append hashed cache/suite/primary evidence into a new
workspace JSON. Functional evidence never supplies performance medians: final
comparisons require matched fresh samples and an explicitly accepted serialized
window. Preserve original cohorts, full result/guard checks, and unmeasured null
phases. The adapter's synthetic file-contract tests do not run reductions.
