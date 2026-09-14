# Main CI and native Cargo measurements, 9 September 2026

The measured candidate is `5181661ec340ebfb181a0045dac79fcec4f35525`, based on
main `395610143576507503fd2c785db3ba62340f4277`. Its manifests, lockfile, crates and
tests are unchanged from that main revision. It includes the CI improvements,
including the narrower Python runtime dependencies, without the FeynKit split.

## Local unchanged warm run

Both paths executed the same 1,633 selected tests, including the five Python API
tests; all passed and 155 were skipped. Inventories match exactly by package,
binary, test name, kind, ignored flag and selection. Neither path recompiled
anything. This comparison covers the seven runtime groups, excluding Clippy and
doctests.

| Operation | Native Cargo | Nix |
|---|---:|---:|
| Execute all selected tests with artifacts available | 187.00s | 171.71s |
| Nix evaluation and disposable-store definition copy | — | 12.39s |
| Execution including that setup | 187.00s | 184.10s |
| Compilation-only warm check (`cargo test --no-run`) | 0.62s | — |
| Reuse an already successful Nix result, without evaluation | — | 0.26s |

The three-second end-to-end difference is too small to establish an advantage.
The 0.26s Nix repeat runs no tests: it reuses their successful result. It also
excludes flake evaluation, so it is not the elapsed time of a full flake command.
The execution comparison deliberately changes only a runtime benchmark tag to
make Nix execute every test group again while retaining compiled artifacts.

Native execution uses `cargo nextest run`, the project's test runner, with the
same `ci-optim` build profile and `ci_gammaloop` test profile as Nix:

```sh
cargo nextest run --offline --locked --workspace \
  --features gammaloop-integration-tests/python-api-tests \
  --cargo-profile ci-optim --profile ci_gammaloop \
  -E 'not package(=spynso3) & not package(=gammaloop-workspace-hack)'
```

Literal `cargo test` execution has different filtering and scheduling; the
compilation-only measurement uses `cargo test --no-run` with the matching build
profile and requested features. Native workspace resolution enables every
CI-selected feature, but unifies dependencies across the workspace; Nix retains
several package-specific Cargo contexts. This is a comparison of the two current
workflows with matching test coverage, not identical Cargo unit graphs. Cargo's
incremental compilation is enabled for the native persistent target. Both routes use pinned Rust 1.97.0, nextest 0.9.140 and the
same eight assigned CPUs (8–15). Runs are sequential on a shared host; Nix may
run four builders within those eight CPUs. A single pair is indicative.

The same prebuilt Python extension is supplied to both routes. Building the
native Rust workspace with `python-api-tests` enables the integration harness;
it does not build the Python-enabled extension. Extension preparation is outside
this warm comparison. Native compilation preparation took 885.41s; Nix had some
seeded project artifacts and also prepared the extension. Those unequal starting
conditions do not constitute a controlled clean-build comparison.

Raw evidence is under
`/tmp/gammaloop-ci-validation/main-native-comparison/local/`, including
`warm-plan.json`, `warm-phases.json`, `warm-timings.csv`,
`warm-comparison.json`, `warm-inventory-comparison.json`, command logs and the
17 per-package Nix archive inventories.

## Local small-source-edit run

The second comparison makes one equivalent implementation edit in
`clear_interrupt_request`: it directly performs the same relaxed atomic store
as the existing setter. Both copies start from the warmed application revision
above. The six ordinary Rust groups execute 1,628 tests; all passed, with 155
skipped and exactly matching inventories. The five Python API tests and their
separately built extension are excluded from this source-edit experiment.

| Operation | Native Cargo | Current Nix configuration |
|---|---:|---:|
| Compile tests / prepare test archives | 44.18s | 1,381.95s (23m02s) |
| Execute the six Rust groups | 182.37s | 177.42s |
| Nix evaluation and definition copy | — | 12.75s |
| Combined test workflow | 226.54s (3m47s) | 1,572.12s (26m12s) |
| Cargo compilation messages during preparation | 3 | 42 |

Formatting and `cargo check` ran before native compilation, as required by the
contribution guidance. They took 3.46s and 245.35s respectively, separately from
the timed test workflow. Including those preliminary checks, the native sequence
took 475.35s (7m55s). No test was changed. Both temporary source edits were restored
and their original hashes verified after the measurements and inventory checks.

The package boundaries remain useful: 13 of 16 ordinary Rust package archives
are identical. Only GammaLoop, its API and the integration-test package archives
change; the clinnet, linnet, spenso and vakint group archives remain identical.
However, the affected Cargo contexts still do excessive preparation for this
edit. Native compiles GammaLoop, its API and integration tests once each. Nix
logs GammaLoop five times, the API five times and integration tests three times,
plus upstream rebuilding including Symbolica. Its 42 messages also include six
generated dependency-anchor crates. No unchanged baseline Cargo derivation
compiled locally; the changed source invalidates the consumer contexts, which
then rebuild these Cargo units. This differs from NixCI rebuilding an identical
published derivation on another worker.

The code explicitly rebuilds workspace libraries against each artifact merge's
dependency metadata. Together with the different Cargo unit graphs and native
incremental state, this explains why per-package archive reuse does not produce
native incremental-build latency. The Nix number includes compilation and archive
preparation. It is a local-store experiment with NixCI uploads disabled, so the
remaining local build design needs improvement independently of remote behavior.

These are single sequential samples on the same eight assigned CPUs. This
comparison does not isolate the effect of incremental compilation alone, compare
old and new Nix implementations after the same edit, or establish cold-build
performance. It demonstrates a large gap in the current edit-and-test workflows.

Exact commands, phase timings and CSV are in `local/native-edit-phases.json`,
`local/nix-edit-phases.json` and `local/edit-timings.csv` under the measurement
root. `local/edit-comparison.json`, `local/edit-compilation-details.json`,
`local/nix-edit-package-archive-boundaries.json` and
`local/edit-inventories/comparison.json` retain the timing, compilation, cache
identity and complete coverage evidence.

## Main NixCI follow-up

The additional main suite was explicitly requested after the original 14-suite
campaign. It was pushed at 15:35:55.941 UTC and completed at 16:27:21 UTC.
Both revisions report all 11 required checks successful, including Clippy,
doctests, formatting, graph validation and the seven runtime groups.

| Observed metric | Previous main, 983fc214f | Current main, 5181661ec |
|---|---:|---:|
| Suite start to required checks | 4h03m05s | 50m28s |
| Push start to required checks | 4h03m34s | 51m00s |
| Whole suite, from suite start | 4h03m08s | 50m53s |
| Push start to first selected worker, including replaced attempts | 70m27s | 21m17s |
| Observed worker minutes, lower bounds | 330.82 | 96.77 |
| Reported intermediate restored content | 16.90 GiB | 23.78 GiB |
| Publication intervals summed across workers | 2h24m48s | 30m11s |
| Cargo `Compiling` messages | 2,015 | 897 |
| Exact Cargo artifact derivations built in multiple jobs | 38 | 5 |
| Final success after the last required check | 3s | 25s |

The [previous main run](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/983fc214fc27a78b08c79f31cccd14380910d3a7)
and [new main run](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/5181661ec340ebfb181a0045dac79fcec4f35525)
show faster observed completion. Queues and remote cache state differ, and the
baseline warmed artifacts used by the candidate. Replaced worker logs also make
resource accounting incomplete. These numbers cannot isolate the speedup caused
by the Python dependency change or establish a reliable expected saving.

The narrower Python preparation does work on main: its
[module job](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/5181661ec340ebfb181a0045dac79fcec4f35525/77e4cec7-aa67-4672-8c8c-d3862183289d)
completed in 9.17s with no Cargo compilation and without a separately scheduled
production GammaLoop producer. The integration group was the final required
check. The new Clippy and doctest cache experiments were not in this push.

Seven runtime/doctest workers have execution summaries matching their baseline
counts. The eighth, core, has a successful result but its final nested-Nix log
contains only abbreviated copying/evaluation progress and `Run succeeded.`.
Execution versus reuse cannot be established from that log, so the reporter
retains `unknown`; this is not evidence that all groups reexecuted. The separate
local full-inventory comparison above verifies all 1,633 selected runtime tests.

The remote warm goals remain unmet: substantial compilation and restoration
continue. Reported intermediate content increased in this sample. There is no
measurement of network bytes or CHF savings, and no controlled project-cold
comparison in this follow-up.

### Unexpected provider behavior retained for investigation

- The workspace-hack job exposed three worker starts under one URL, and the
  doctest job exposed two. Every observed stream is retained. Their retry or
  replacement causes are unknown.
- [Linnet-py's producer](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/5181661ec340ebfb181a0045dac79fcec4f35525/9bb8c315-9195-49d7-a40b-ed0d777f055e)
  reported uploading the `gnjpz47gk5hrbrkh9mwf65qd063lcsc5` test artifact at
  16:14:34 UTC. The later
  [Linnet runtime worker](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/5181661ec340ebfb181a0045dac79fcec4f35525/40e47d17-d48e-4b2a-bdc6-f23904ec757d)
  rebuilt that identical derivation at 16:19:08 UTC. This is an exact identity
  match; cache visibility, retention and substitution need investigation.
- The [graph-validation worker](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/5181661ec340ebfb181a0045dac79fcec4f35525/f3037a10-8dd5-46b5-a8a7-e09bf15f1334)
  took 8m03s with zero Cargo compilation. It rebuilt 310 vendored source outputs
  and had 6m03s of publication activity. This is not time uploading test binaries.

The paired JSON/CSV report and retained logs are under
`/tmp/gammaloop-ci-validation/main-native-comparison/remote/reports/watched-suites-2026-09-09T16-27-51.135Z/`.
A compact summary, per-group outcomes and exact repeated artifact links are in
`remote/main-comparison-summary.json` under the same measurement root.
