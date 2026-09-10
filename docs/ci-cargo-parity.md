# Closing the gap with native Cargo

The objective is comparable **edit-to-test completion**, including Nix evaluation,
artifact restoration and export, compilation, and execution of the same tests.
A cached successful result that does not execute tests is a separate measurement.
The current change improves main, but **does not achieve Cargo parity**.

## Integrated cache implementation, 10 September 2026

The follow-up combines four changes:

- Update Crane to upstream revision `eb35abda9f232cc6610b1d1e3200d15c49b7ac54`
  (3 September, after v0.24.0); other existing input pins stay unchanged.
- Apply the existing target/build-script regex feature alignment to every context
  already consuming Symbolica. Cargo unit-graph inspection finds matching Symbolica
  and upstream library units, without changing GammaLoop's resolved graph.
- Enable incremental compilation consistently in per-crate library/test producers
  and their nextest archive consumers.
- Reuse compiler state from a pinned, compatible application revision. Normal
  dependency artifacts remain separate from compiler state, so downstream crates
  do not restore their parents' query caches.

`ci-cache-base` is a source-only flake input. The current CI implementation and
Rust toolchain prepare that revision's per-context caches. A manifest or lockfile
change disables reuse; an identical full derivation identity reuses its existing
outputs. Changed producers restore only compiler state and retain current Cargo
source fingerprints and dependency outputs. The library/test split is preserved.

The seed has separate `out` and `incremental` outputs. A changed producer exports
only its compiled artifacts: its newly generated query cache would be unused by
subsequent builds from the fixed baseline. NixCI producer wrappers publish the
current artifacts and reusable baseline state; ordinary consumers use raw `out`.
The initial baseline should be prepared before comparing changed revisions.

Refresh deliberately after a green application revision:

```sh
just ci-cache-base FULL_COMMIT_HASH
```

This prepares a new seed once using the current CI definition. It does **not**
automatically promote the latest successful commit's compiler state. Each edit
reuses the pinned baseline independently, avoiding a history of cache archives.
FeynKit needs its own compatible pin after rebasing and regenerating its CI graph.

Crane's normal compressed deltas remain enabled for compiled dependency artifacts.
Compiler-state snapshots are self-contained: a filesystem probe demonstrates that
rustc-style new-session hardlinks with epoch-1 timestamps disappear from Crane's
mtime-only deltas. Touching new paths fixes inclusion but does not represent
removed files. Full state snapshots reconstruct exactly. Latest Crane also fixes
inherited `.cargo-build-lock` and `.cargo-artifact-lock` handling; two successive
restorations passed the standalone regression probe.

All seven FeynKit invalidation scenarios match the previous implementation.
A GammaLoop edit selects six affected library/test producers for compiler-state
reuse; the extracted libraries remain unchanged. Python selects two, CFF twelve,
model twenty-four, and a fixture six. Manifest and lockfile edits select zero.
The existing Python consumer-expectation discrepancy described below is unchanged.

Publishing the baseline state adds that archive to the producer wrapper's closure.
Warm NixCI measurements must include this transfer, even though downstream build
and runtime outputs do not retain compiler state. Intermediate multi-output cache
publication and skipped-wrapper behavior still need remote verification.

The integrated main source-edit run finishes in **6m35s**, versus **17m47s** for
the earlier candidate and **3m47s** for native Cargo. It executes the same 1,628
tests, with 155 skipped. Names, target kinds, ignored flags and filter selections
match exactly. The full prepared revision also passes all 1,633 tests, including
Python. These source-edit timings exclude the five Python tests, Clippy and
doctests, consistently with the native comparison.

| Same main source edit | Artifact/compilation preparation | Tests | Total including setup |
|---|---:|---:|---:|
| Native Cargo, persistent incremental target | 44.18 s | 182.37 s | **3m47s** |
| Earlier Nix candidate | 843.06 s | 203.99 s | **17m47s** |
| Integrated Nix candidate | 185.90 s | 182.31 s | **6m35s** |

The integrated total includes 26.47 seconds of evaluation and derivation setup.
It is 63% faster than the earlier candidate, but remains 74% slower than Cargo.
These are single sequential samples on a shared host with the same eight CPUs,
Rust 1.97, nextest 0.9.140 and test profiles; they are not controlled cold runs.

Only GammaLoop, its API, integration tests and their synthetic feature anchors
compile after the edit: 12 Cargo compilation messages across six producers.
Their Cargo stages take 13.63/28.74, 13.92/24.78 and 13.88/20.78 seconds respectively
(library/test). Symbolica and the unchanged workspace libraries are reused.
All final package archive steps compile nothing. The serial API and integration
input merges still take 38.40 and 30.34 seconds; another 13.43-second merge prepares
the integration test producer. Repeated restoration and export remain a bottleneck.
These overlapping builder intervals must not be summed as whole-suite wall time.

Compiler state is absent from all 64 materialized raw-output closures inspected
across the prepared and edited revisions. Four optional producer entries were not
materialized. The 30 unique prepared compiler-state outputs total **2.24 GiB of
compressed files**; GammaLoop's test state is about 420 MiB and its library state
268 MiB. This is storage footprint, not measured network traffic. A seeded
GammaLoop negative control fails on an intentional compile error even when the
source timestamp is reset to epoch zero.

Preparing caches with the new Crane revision takes 28m35s for package archives.
The subsequent full runtime phase takes 15m14s, but includes 65 compilation
messages from the separate Python-module preparation. It must not be described
as 15 minutes of test execution or compared with the warm edit benchmark. The
current Python extension enables additional Python/ABI features, so the ordinary
API test library cannot substitute for it unchanged.

Clippy, formatting and generated graph/configuration checks pass. Doctests pass
43 tests, with 20 ignored. This separate static run takes 10m14s after the Crane
update; its workspace-wide compilation is not included in the edit timing.
FeynKit validation also passes on its isolated application revision `07d0efc91`.
The full prepared revision runs **1,834 tests**, with 127 skipped, including all
five Python tests. Exact inventories match across all 26 package archives.
The independently prepared GammaLoop edit runs **1,829 ordinary tests** in
**7m19s**: 37.20 seconds of setup, 218.09 seconds of artifact preparation and
183.87 seconds of test execution. Its 25 package inventories also match exactly.
Only GammaLoop, its API, integration tests and six synthetic feature-anchor
compilations occur; all nine extracted libraries and Symbolica remain cached.
Final archive steps compile nothing. This confirms reuse on FeynKit; it does not
measure a speedup against FeynKit's earlier implementation or native Cargo.

FeynKit cache preparation takes 33m49s for archives, then 18m33s for the runtime
phase including its separate Python-module preparation (90 compilation messages).
These are preparation timings on the existing local store, not cold results.
The main-only NixCI seed, unchanged and source-edit runs are being collected
separately; the historical experiments below describe earlier implementations.

## Separate Python artifact-reuse experiment

A scratch-only main experiment prepares the Python module from the existing API
test-library artifacts. Its Python features, ABI, interpreter, fixtures and test
selection stay unchanged. The old production preparation graph is absent from
both the candidate module and its runtime check.

Module preparation takes **245.37 seconds** with 14 Cargo compilation messages;
Symbolica and vakint are reused. The remaining Python-specific variants still
compile. Packaging compiles nothing. The runtime check takes **1.77 seconds**,
including **0.70 seconds executing the same five passing tests**, with an exact
inventory match. There is no previous-revision incremental state for this Python
variant yet.

This experiment is not included in commit `5b87a850e` or its three approved NixCI
runs. The older 914.17-second runtime phase included ordinary tests in parallel,
so it cannot supply an isolated Python speedup percentage. A paired preparation
measurement and the FeynKit Python layout still need validation before integration.
Evidence: `python-graph-probe/timing-v1/{phases,result}.json` under the main
measurement directory.

## Main source-edit measurement, 9 September 2026

Application revision: `395610143576507503fd2c785db3ba62340f4277`.
CI baseline: `08ffdaf61`, with the implementation previously pushed as `5181661`.
Candidate implementation: `8043bd6ab`.
Each variant starts from its prepared clean application revision and replaces
`set_interrupted(false)` inside `clear_interrupt_request` with the identical
atomic store used by the setter. The temporary source edit is not part of the change.

| Local run | Compilation/artifact preparation | Test execution | Total |
|---|---:|---:|---:|
| Native Cargo/nextest, persistent incremental target | 44.18 s | 182.37 s | **3m47s** |
| Previous Nix implementation | 1381.95 s | 177.42 s | **26m12s** |
| This Nix candidate | 843.06 s | 203.99 s | **17m47s** |

Nix totals include evaluation and derivation setup: 12.75 seconds before and
20.17 seconds after. Artifact preparation improves by 39%; total completion
improves by 32%. There are 19 Cargo compilation messages instead of 42, including
three synthetic feature anchors instead of six. GammaLoop and API each compile
twice instead of five times; integration tests compile twice instead of three times.

These are single sequential samples on a shared host, using the same eight CPUs,
Rust 1.97.0, nextest 0.9.140, `ci-optim` and `ci_gammaloop`. The native preflight
`cargo check` took another 245 seconds and is reported separately. Cache preparation
histories differ; these are not controlled project-cold measurements. The runtime
variation is another reason to avoid treating the measured percentages as promises.

All three source-edit runs execute the same **1,628 tests**, with **155 skipped**.
This timing excludes Clippy, doctests, and the five Python API tests. A separate full
candidate run passes **1,633 tests**, including Python, with **155 skipped**. Exact
inventories match, including names, target kinds, ignored tests, and filter selection.

## Earlier narrow implementation

The existing feature anchors align the regex dependency flags used by GammaLoop,
its API and integration tests, on both the target and build-script sides. Cargo's
full recursive GammaLoop compilation-unit signatures now match in those contexts.
Without the build-script alignment, the Symbolica dependency graph still differs.
In that earlier implementation, only contexts that depend on GammaLoop received the additional anchor dependencies;
independent FeynKit crates retain their existing feature sets.

The dependency producer now runs its existing test-compatible library build once;
the unused argument helper for the removed command is also removed.
The preceding library build used a different feature context and created work that
could not be reused by that build. Normal libraries remain separate from test
binaries, preserving fixture-only invalidation boundaries. Source filters, profiles,
archive commands, output names and scheduling are unchanged.

In the source-edit measurement, 13 of 16 package test archives stay unchanged.
The remaining work is still substantial: GammaLoop's dependency stage recompiles
Symbolica and several unchanged workspace libraries, and each affected package
still has separate library and test compilation. This change does not provide
previous-revision incremental compiler state.

## FeynKit and rejected approaches

On the isolated FeynKit integration revision `07d0efc91`, the candidate has exactly
the baseline's invalidation sets for GammaLoop, Python-binding, CFF, model, embedded
fixture, manifest and lockfile edits. All nine extracted libraries remain stable
for a GammaLoop-only edit. These are derivation inspections, not FeynKit runtime
or speed measurements.

The existing probe expects only one production package to change for a Python
binding edit, but both baseline and candidate change `feynkit-py` and its API
consumer. Its existing assertion remains unchanged and reports that discrepancy.

A trial that made downstream libraries inherit their parents' combined library/test
artifacts was rejected: an embedded fixture edit invalidated 13 FeynKit test-cache
entries instead of seven. The smaller candidate preserves seven.

## Incremental-state experiment

A separate, unintegrated GammaLoop-only probe explicitly gives the edited Nix
producer the preceding revision's compiler cache. It rebuilds GammaLoop in
**22.41 seconds of Cargo time**, or **55.80 seconds including local Nix restoration
and artifact export**. The same combined producer without incremental state took
4m14s inside Cargo. This is a one-crate diagnostic; it is not comparable to the
44-second native compilation of all three affected packages as a whole-suite result.

The resulting test archive requires no additional compilation and passes all
456 selected GammaLoop tests, with 61 skipped and an identical 517-test inventory.
A deliberately inserted compile error fails even with source timestamps reset to
the Nix epoch, demonstrating that saved binaries do not mask source changes.
Tests run in the existing Nix fixture/snapshot environment; an initial attempt to
execute the archive outside that environment failed and is not a correctness result.

The full incremental archive grows from **1.74 GiB compressed / 5.74 GiB expanded**
to **2.28 GiB / 7.03 GiB** after the edit, retaining two compiler-cache sessions.
That transfer cost makes blindly exporting incremental caches unsuitable for rollout.
Crane's existing directory/symlink mode does not improve this result: converting
that seed takes 94.57 seconds, and the first edit takes **216.52 seconds** despite
only 29.62 seconds of Cargo compilation. Its output closure has 11.61 GiB of NAR
content versus 2.90 GiB for the compressed edited output; these are restored-content
sizes, not network bytes. A second reuse fails because `.cargo-build-lock` has become
a symlink into the read-only store. This storage mode is rejected; neither the slow
export nor the lock-file failure is evidence of a NixCI service fault. The former pinned
Crane inheritance hook removes `.cargo-lock`, but does not remove the
`.cargo-build-lock` used by this Cargo version.

## Earlier validation and evidence

In the earlier narrow implementation, Clippy, doctest, formatting and workspace-graph derivations are unchanged from the
baseline. **Clippy, formatting and generated graph/configuration checks passed
locally.** The unchanged doctest check was not rerun in this experiment.
Removing the unused helper also leaves all seven measured test-archive roots
identical to their validated derivations.

The static checks were built with:

```sh
nix build --no-link \
  .#checks.x86_64-linux.gammaloop-clippy \
  .#checks.x86_64-linux.gammaloop-fmt \
  .#checks.x86_64-linux.gammaloop-guppy-workspace-graph
```

No additional NixCI suite has been pushed for this change, and no cold-build
or remote worker-time improvement is claimed.

Local evidence is under `/tmp/gammaloop-ci-validation/cargo-parity-main/`:

- `library-comparison.json`, `edit-libraries-v1/phases.json` and
  `edit-libraries-v1/compilation-details.json`: timings and compilation contexts.
- `prepare-libraries-v1/inventories/comparison.json` and
  `edit-libraries-v1/inventories/comparison.json`: exact full and edited inventories.
- `library-package-boundaries.json`, `library-boundary-comparison.json` and
  `feynkit-boundaries-{baseline,libraries}/report.json`: cache invalidation.
- `aligned-v2-unit-signatures.json`: matching Cargo compilation graphs.
- `incremental-state/{phases,archive-sizes,archive-contents,validation,runtime-nix,inventory-comparison}.json`:
  the incremental-state timing, footprint and freshness controls.
- `incremental-directory/{phases,sizes,decision}.json`: storage-mode timing,
  closure size, and the failed second reuse.
- `measure.py`, `measure-incremental.py`, `incremental-probe.nix` and
  `feynkit-boundaries.py`: the external measurement harnesses.

The previous Cargo and Nix measurements are documented in
[the native comparison](ci-main-native-comparison.md). The integrated implementation above applies the resulting compiler-state experiment
separately to each affected library and test producer. Remote transfer cost, cold
completion and the remaining artifact-merge overhead still require measurement.

The portable [timing and validation summary](ci-cargo-parity-results.json) is
committed alongside this document. Detailed integrated evidence in the same directory:

- `edit-integrated-v2/phases.json` and `builder-timings.json`: full timing and
  overlapping per-builder phase intervals.
- `integrated-validation/{prepare,edit}-integrated-v2/inventories/comparison.json`:
  exact test inventories.
- `integrated-validation/{closure-validation,compiler-state-sizes,negative-control-result}.json`:
  state separation, storage footprint and stale-source rejection.
- `integrated-validation/feynkit-boundary-comparison.json`: unchanged invalidation sets.
- `crane-update-probes/latest.json`: repeated lock cleanup and archive regression probes.

The CI log reporter's 23 existing tests pass. The first of three approved main
NixCI suites is [running the committed integrated candidate](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-cargo-parity-main-seed/5b87a850e77900029a487b3b8c77d1c7138fa8da).
Its unchanged follow-up and independent source-edit commits have the exact
locally validated archive roots. Remote conclusions await those completed runs.
FeynKit timing and inventory evidence is under
`/tmp/ci-dependency-reuse-audit/integrated/benchmark/`, including
`integrated-result.json` and `inventory-comparison-v1.json`.
