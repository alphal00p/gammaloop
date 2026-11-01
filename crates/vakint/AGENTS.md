- If running into a linking issue complaining about missing `__emul...` symbol, then it may be a macos specific issue, that can be fixed by building with the following environment variable change `EXTRA_MACOS_LIBS_FOR_GNU_GCC=T` (see impact of this in file `build.rs`.
- A quick way to run all tests except pySecDec related ones (which can be slow) is the following command:

VAKINT_SKIP_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --package vakint --no-default-features

## Reduction benchmark records

`examples/benchmark_reduction.rs` measures fresh-process first-call timings
(FMFT against RustRed's ordered scheduler) on the four-loop cases of
`tests/fixtures/benchmark_reduction.json` and appends each completed, validated
run to `docs/benchmark/data/results.json`. From the repository workspace root:

```sh
taskset -c CPU cargo run --locked --release -p vakint --example benchmark_reduction \
  --features feynkit-ingress,experimental-rustred,weighted-lowering-diagnostics,symbolica_community_module -- \
  --form /path/to/form \
  --work-dir /path/to/workspace/vakint-benchmark \
  --output crates/vakint/docs/benchmark/data/results.json \
  --build-profile release --select-current
```

Replace `CPU` with one permitted CPU on Linux. `--artifacts PATH` selects an
existing external-data checkout; otherwise the official `latest` loader is
used. The JSON is append-only: earlier runs, including the historical
depth-first and cache cohorts, remain unchanged records and are never presented
as current measurements. Never hand-edit timings. Declare the actual Cargo
profile with `--build-profile`; final speed comparisons need a `--release` build.

Vakint reduces only with RustRed's ordered scheduler. Depth-first reduction and
its prepared/running caches were removed.

For physical rung measurements use `examples/ordered_reduction_rungs.rs` (one
engine per fresh process, `--no-monitor` for timings). All rung cases use MSbar.
The historical depth-first results are exact references (the frozen rung
fixtures and the saved five-loop result), never new timing samples. Reduction
timing uses `reduction_math_seconds` (IBP plus exact weighted
assembly/grouping, requires `experimental-rustred`), kept separate from the
broader public-call timer. See `docs/ordered_reduction.typ` for the phase and
grouping-counter semantics.

The five-loop physical collector (`examples/five_loop_qcd_ladder.rs`) saves its
complete scalar input as `prepared-input.json` plus authenticated
input/condition payloads; `--prepared-input PATH` shares it across fresh
processes. Keep native preparation timing separate from reduction timing.

The Typst benchmark report, its maintenance contract, the Python merge
adapters and the provenance archives that used to live under `docs/benchmark/`
were agent-generated evidence and are no longer part of the repository.
