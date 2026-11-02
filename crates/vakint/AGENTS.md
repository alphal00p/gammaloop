- If running into a linking issue complaining about missing `__emul...` symbol, then it may be a macos specific issue, that can be fixed by building with the following environment variable change `EXTRA_MACOS_LIBS_FOR_GNU_GCC=T` (see impact of this in file `build.rs`.
- pySecDec tests are opt-in (`RUN_PYSECDEC_TESTS=1`). Run the fast local gate with
  `just test_vakint` and the full CI selection with `just test_vakint_full`
  (see "Testing" in `README.md`). Keep slow cases out of the gate by adding them
  to the `test_vakint` default-filter in `.config/nextest.toml`, never by
  wrapping them in `slow::` modules: CI skips those.
- Integration tests share their setup: `tests/test_utils.rs` (evaluation lanes,
  `compare_evaluations`, `compare_exact_peers` at a `ReferencePoint`),
  `tests/support/rustred.rs` (FORM-free RustRed settings, the pinned
  external-data release) and `tests/support/reference_table.rs` (Laurent tables).
  Each test should read as one short call of such a helper.

## Readable physics tests

- `tests/reference_tables.rs` (fast): one test per named one- to four-loop
  integral, written inline in `topo(prop(...))` notation. vakint's default
  engine (RustRed, Ordered scheduler, numerical masters, no FORM) is compared
  coefficient by coefficient with the frozen MATAD/alphaLoop/FMFT (and pySecDec)
  values of the original tests. Four-loop cases need the official external
  data (`VAKINT_TEST_EXTERNAL_DATA`). Run with `--nocapture` to see the tables.
- `tests/ordered_reduction_rungs.rs`: the physical gluon-rung ladder. One to
  four loops (fast) check exact equality with the frozen depth-first results
  and exact = modular weights, then print a table against the stored
  high-precision golden. The five-loop fixtures (`#[ignore]`, slow) keep their
  exact/sha256 gates and then print a master-product summary: the leading
  Laurent terms of every master-product coefficient and a master legend.
- New tables reuse `tests/support/reference_table.rs`. Never loosen an
  existing exact or digest assertion to make a table pass.

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
