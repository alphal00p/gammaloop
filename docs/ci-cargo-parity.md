# Closing the gap with native Cargo

The objective is comparable **edit-to-test completion**, including Nix evaluation,
artifact restoration and export, compilation, and execution of the same tests.
A cached successful result that does not execute tests is a separate measurement.
The current change improves local main source-edit completion, but **does not
achieve Cargo parity**.
Local correctness and cold-latency gates pass. Remote zero-recompilation and
final-success goals remain unmet; the integrated compiler-state candidate is
not accepted for merge as a remote performance improvement.

The latest [artifact-preparation follow-up](ci-artifact-preparation.md) removes
redundant merges and cache references, shares Python/static-check dependencies,
and makes publication stable. Its fresh main source-edit pair, including Python,
finishes in **4m53s versus 10m12s**. The measurements below describe the preceding
implementation; they remain useful baseline evidence.

A subsequent [matched cold Cargo/Nix comparison](ci-cold-cargo-comparison.md)
includes Python preparation and all 1,633 runtime tests in both routes: **21m46s
with direct Cargo versus 28m48s through Nix**. It excludes Clippy/doctests and is
separate from the source-edit measurements below.

## Earlier integrated cache implementation, 10 September 2026

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
subsequent builds from the fixed baseline. Normal NixCI producer wrappers now publish only the
current artifacts; ordinary consumers use raw `out`. The optional
`ci-compiler-state` output exposes the pinned seeds for local/cache preparation.
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
model twenty-four, and a fixture six. Manifest and lockfile edits disable
prior-state reuse and retain their broad derivation invalidation.
The existing Python consumer-expectation discrepancy described below is unchanged.

The earlier NixCI measurements below include baseline-state publication in each
producer wrapper, even though downstream build and runtime outputs do not retain
compiler state. That version restores state on unchanged commits. The latest
follow-up removes those links and the revision salt; reliable remote retrieval
of stable artifacts and separately published seeds remains to verify.

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
FeynKit's boundary checker and its self-test pass; all nine extracted libraries
remain independent of GammaLoop and its private workspace-hack dependency.
Its Python-stub workflow is preserved, but stub generation/staleness checking is
not part of these cold builds and was not rerun. That workflow builds a separate
`python_stubgen` executable with additional optional dependencies and a different
profile; the measured Python module cannot substitute for that check.
The main-only NixCI measurements below assess the integrated implementation;
the historical experiments later in this document describe earlier implementations.

## Controlled local cold follow-up

All four builds pass the complete measured scope: seven/eight runtime groups,
Python, Clippy, licensed doctests, formatting and generated graph/configuration
checks. The baseline is the earlier narrow CI candidate, before the four changes
above, **not original main CI** (`b9d6071fc` versus `5b87a850e` on main).
FeynKit uses the same application revision within its own pair.

| Layout | Baseline | Candidate | Completion change | Passing tests / doctests per variant |
|---|---:|---:|---:|---:|
| Main | 46m26s | 44m01s | −5.19% | 1,633 / 43 |
| FeynKit | 48m15s | 48m04s | −0.38% | 1,834 / 62 |

Exact archive inventories match across all 17/26 packages, including filters and
ignored cases. All selected nextest cases execute successfully; exact doctest
names and outcomes match, with 20 ignored per variant. Archive steps compile
nothing. All 19/21 selected roots pass. Exact repeats take 0.70–0.77 seconds with
zero builders; those repeats reuse results and do not execute tests again.
Both pairs meet the no-more-than-10%-slower cold goal. The FeynKit difference is
too small to treat as a meaningful improvement in one shared-host observation.

Combined builder occupancy increases: **123.32 → 126.57 minutes** on main (+2.64%)
and **162.64 → 163.70 minutes** on FeynKit (+0.65%). These are sums of overlapping
builder wall intervals, not CPU or billed time. Main has 807 → 794 Cargo compilation
messages; FeynKit has 896 in each variant. The new compiler-state outputs increase
total produced project NAR content from 29.88 → 31.54 GiB and 44.20 → 46.69 GiB
respectively, despite slightly smaller ordinary artifact totals. This is local
output storage, not restored content or network traffic. There is **no demonstrated
cold resource saving**.

Completion improvements are uneven. Main core/Python/integration groups finish
roughly 2½–3 minutes earlier. On FeynKit, extracted-crate tests finish 1m56s earlier
and core 1m42s earlier, while Clippy finishes 1m10s later, doctests 1m33s later and
Python 1m02s later. Integration finishes 24 seconds earlier; Python is last.

Waiting matters independently of build speed. FeynKit Clinnet waits **505 → 302
seconds** after all direct inputs are ready, with all four build slots occupied
during over 98% of that wait. Its archive becomes ready 237 seconds later, so its
runtime group still finishes 35 seconds later. These are observed local delays,
not a NixCI anomaly or guaranteed removable critical-path time. The portable JSON
records direct-input readiness, waiting and completion for every runtime group.

Each variant starts in a fresh store with the same external seed within its
layout, no project outputs or compiler state, CPU affinity 8–15, eight cores and
four Nix jobs. Remote substitution, remote builders and upload hooks are disabled.
Evaluation/graph inspection and setup copies are separately recorded. Store
validation and other preflight work are excluded; those fields do not sum to
launch-to-finish time. Memory, storage and other host workloads remain shared.

Evidence is under `/tmp/ci-dependency-reuse-audit/cold-followup/`, including each
layout's `comparison.json` and `detailed-summary.json`. After all four builds,
postprocessing needed a missing extraction directory created for nextest. The
original collector and failed logs are preserved alongside the corrected copy.
The scratch summary parser also needed to handle nextest's `(3 slow)` annotation
before `67 skipped`; corrected counts match 127 skipped on both FeynKit variants.
Original summaries are retained. No timed build or test assertion changed or was rerun.

## Integrated NixCI measurements, 10 September 2026

All three approved runs pass, but **the unchanged run misses the zero-Cargo-
recompilation and final-success goals**. These scenarios use the same implementation;
the seed is not an implementation baseline and cannot supply a causal speedup
percentage. It also uses shared caches, so it is not a controlled cold run.
The [9 September remote improvement](ci-measurement-results.md) predates this
implementation. These three runs do not isolate an additional NixCI speedup.

| Main scenario | Required checks | Whole suite | Push to required | Observed worker minutes | Reported restored content |
|---|---:|---:|---:|---:|---:|
| [Prepare new caches](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-cargo-parity-main-seed/5b87a850e77900029a487b3b8c77d1c7138fa8da) | 1h41m49s | 1h42m09s | 1h47m21s | 238.99 | 31.97 GiB |
| [Unchanged commit](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-cargo-parity-main-warm/f1a086ead8a176a6e63c3d686b9317c841d373e8) | 8m56s | 13m31s | 9m02s | at least 41.25 | at least 27.74 GiB |
| [GammaLoop edit, with retries](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-cargo-parity-main-edit/3da4b1d3da993cfb23d8f0bf70f07c2aa307af4c) | 19m10s | 19m21s | 19m18s | at least 66.17 | at least 20.26 GiB |

Worker minutes are observed log spans, not billed compute. Content sizes are not
measured network bytes. The seed observer had a 36-minute gap; transient attempts
inside it may be missing even though the collector retrieved all available logs.
The unchanged run has a verified replaced worker log, making resource totals
lower bounds. Its earlier observed 33.24 seconds and 199.4 MiB are kept separately
because overlap with the replacement is unknown.

The seed executes all 1,633 nextest cases and 43 doctests. In the unchanged run,
all six ordinary group checks finish by **2m16s**, but only 390 executed cases are
visible (Clinnet and Spenso). Four ordinary groups and Python have successful jobs
with insufficient evidence to distinguish execution from cached result reuse.
The 43 doctests execute and determine required completion. Therefore 2m16s is
not a measurement of executing all tests and cannot be compared with native Cargo's
3m47s edit-and-test run. Local exact-inventory validation remains the coverage proof.

Two distinct problems remain:

- **Our scheduling:** all 33 selected publication wrappers change with each commit.
  The unchanged run restores 30 compiler-state outputs, totaling **2.24 GiB** of
  reported content and 131.41 seconds of summed transfer intervals. Every named
  state download matches its own producer wrapper by its unique logged name;
  remote archive hashes are not available. These are not inherited parent states
  in ordinary artifact consumers. Ten preparation jobs finish after all
  required checks, extending final success by **4m35s**. The aggregate itself
  finishes six seconds after the final producer. This preparation did not help
  that run's test completion.
- **Unexplained remote reuse failures:** the unchanged run logs 162 Cargo compilation
  messages (73 in the doctest job), including builds of exact outputs uploaded by
  the seed. The same workspace-hack dependency is built in three seed jobs, including
  starts five and sixteen minutes after its first successful upload. The unchanged
  doctest output is also rebuilt after its seed upload. Successful upload messages
  do not establish later cache availability; the logs do not identify the root cause.

Three package archives (`gammaloop-tracing-filter-macros`, `symbolica-utils` and
`linnet-py`) list zero test cases in both layouts. Their unchanged-run binary
workers consume 100.44 observed worker-seconds and restore 1.12 GiB, compiling
nothing. A separate `linnet-py` dependency producer takes 64.92 seconds and logs
32 Cargo compilation messages, although no other workspace package depends on
`linnet-py`. Its own test binaries still consume that library. The other two
libraries are needed by real test consumers. These are additional preparation
costs, but no completion-time saving is established: the binary workers start
after the ordinary groups finish. Empty inventories alone do not establish that
compile-only coverage is redundant, so no coverage removal is included.

A build-free scheduling probe removes only the per-commit identity change. Across
the controlled GammaLoop edit, **7 wrappers change instead of 33**, exactly matching
the six affected library/test producers and Python module; 26 unrelated wrappers
stay stable. This is a follow-up hypothesis, not an integrated change or measured
saving. The original salt worked around NixCI memoizing top-level jobs without
realizing their closures, so dependable remote cache access still matters.

The seed's [Python module worker](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-cargo-parity-main-seed/5b87a850e77900029a487b3b8c77d1c7138fa8da/640af92b-e576-4bfd-abe0-13be90552a12)
spends 49m12s in its observed worker span. Upload activity covers 31m48s; Cargo
commands cover about 6m53s. Many short preparation builds are followed by lengthy
upload intervals. These intervals overlap and include unspecified processing, so
they do not establish pure network time. Nix documents that a
[post-build hook blocks the build loop](https://releases.nixos.org/nix/nix-2.34.8/manual/advanced-topics/post-build-hook.html).
The worker's public closure links curl 8.20.0, so the upstream Nix fix specifically
for curl 8.21 and newer is not an established explanation for these delays.

In the source-edit suite, the six ordinary groups finish in **17m28s**. Visible
summaries prove 1,174 nextest executions across the suite, including five Python
cases; Spenso and Vakint succeed with unknown execution versus reuse. The affected GammaLoop, API and integration Cargo commands
take **7–12 seconds each**, confirming that compiler-state reuse works remotely.
Whole dedicated workers take 66–291 seconds, including preparation and transfers.
This is not a clean speed comparison with the 6m35s local run: hardware differs,
not every test is proven to reexecute, and the remote run has retries and duplicates.

The [dedicated GammaLoop producer](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-cargo-parity-main-edit/3da4b1d3da993cfb23d8f0bf70f07c2aa307af4c/7f86ebce-a314-4984-bc99-0f20fa54320e)
finishes uploading at 11:56:57.195. The [integration producer](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-cargo-parity-main-edit/3da4b1d3da993cfb23d8f0bf70f07c2aa307af4c/bbfdb96d-84ea-46bd-80f6-74da94ce8c71)
starts later and rebuilds the identical derivation at 12:02:50.156. It also rebuilds
four exact seed artifact derivations, compiling unchanged libraries. That repeated
work is distinct from the expected rebuilding of changed application code.

Fourteen HTTP/2 download errors from `cache.nix-ci.com` hit six attempts within
6.18 seconds. Each later succeeds under another job URL. GitHub updates the original
check but keeps its old URL; the reporter now uses validated reciprocal NixCI retry
links to associate its clocks with the replacement. Original failed-attempt check
clocks remain unknown. A replacement Python worker also has an overwritten stream;
its earlier 12.20 seconds and 12.76 MiB remain separate because overlap is unresolved.
All original logs and reports are preserved. No fourth suite was submitted.

The replacement Python module worker takes **11m47s**, including about **4m39s** of
Cargo command time and **6m41s** of upload activity. These intervals overlap. It
still prepares production/Python feature variants separately; the scratch reuse
experiment below is not part of these runs. The final five Python tests execute
in 0.47 seconds, while their check completes last at 19m10s. Final success follows
11 seconds later, meeting the final-tail goal for this scenario.

The controlled local cold follow-up is reported above. It does not establish
a causal remote implementation speedup.

Detailed JSON, CSV, raw logs, replaced streams and snapshots are retained under
`/tmp/gammaloop-ci-validation/cargo-parity-main/nixci-approved-three/`. The portable
JSON beside this document includes all three final scenarios, per-check completion
and execution labels, and their limitations.

## Syd's suggestions and the remaining cache question

Synchronous dependency discovery is already enabled in this candidate. We retain
our explicit generated graph, so this does not adopt all of [PR #104](https://github.com/alphal00p/gammaloop/pull/104).
[PR #105](https://github.com/alphal00p/gammaloop/pull/105) was absent from these
measurements. Its local cache-upload setup is now incorporated in a separate
follow-up change, adapted to the refactored flake. Benchmark runs still disable
upload hooks so transfers do not contaminate the build clock. Setup and limitations
are documented in [CONTRIBUTING.md](../CONTRIBUTING.md#nixci-cache).

The incorporated hook also corrects one conditional issue: the
[pinned upstream hook](https://github.com/alphal00p/gammaloop/blob/13d6826ddee00d7cb4cbd733a837d23c248658f3/flake.nix#L2517-L2533)
changes only `XDG_CACHE_HOME`, but Nix 2.34.7
[uses `NIX_CACHE_HOME` first when inherited](https://github.com/NixOS/nix/blob/2.34.7/src/libutil/users.cc#L15-L24).
Setting `NIX_CACHE_HOME="$XDG_CACHE_HOME/nix"` as well preserves the hook's intended
temporary cache isolation.

Local validation on 2026-09-11 with Nix 2.34.8 confirmed:

- All nine dev shells evaluate, and activation preserves existing Nix settings
  both with and without a readable netrc file.
- The Nix-built hook forwards multiple outputs, isolates both cache directories,
  cleans up, and returns success with a warning when the upload command fails.
- A real build in a disposable local store succeeds when the hook's upload command
  is deliberately made to fail; repeating that cached build does not invoke the
  hook. These use a stub upload command, not an authenticated remote transfer.
- All 53 selected CI derivations, 15 runtime/Python roots and 35 package producers
  retain their derivation paths. Generated `nix-ci.nix` matches, and the flake
  passes parsing and Alejandra checks. This addition does not invalidate them.

The current shell has no cache-token file and reports `trusted: false` for its
Nix daemon. Authenticated upload and signed restore remain unverified here, and
automatic uploads require that machine setup. No NixCI suite was started for this
change. Scripts and raw results are under `/tmp/ci-pr105-integration/`.

Nix also [caches missing paths for one hour by default](https://github.com/NixOS/nix/blob/2.34.7/src/libstore/include/nix/store/globals.hh#L55-L73),
and [checks process memory before its disk cache](https://github.com/NixOS/nix/blob/2.34.7/src/libstore/store-api.cc#L608-L650).
A separate upload cache does not clear another consumer's earlier missing-path
entry. This is a source-supported mechanism to investigate with Syd, **not the
established cause of our duplicate builds**: hosted worker environments, cache
sharing and actual expiry settings are unknown. The source audit and proposed
local-hook patch are under `/tmp/ci-dependency-reuse-audit/narinfo-cache-audit/`.

## Separate Python artifact-reuse experiment

A scratch-only main experiment prepares the Python module from the existing API
test-library artifacts. Its Python features, ABI, interpreter, fixtures and test
selection stay unchanged. The old production preparation graph is absent from
both the candidate module and its runtime check.

A controlled pair starts from two fresh stores seeded with exactly the same
ordinary test artifacts and external prerequisites (16.39 GiB of logical NAR
content). Python-specific/production preparation and compiler-state outputs are
absent from both seeds. Both variants use CPUs 8–15, eight cores, four Nix jobs,
no remote builders/substitution and no upload hook. This measures **Python
preparation after ordinary artifacts are available**, not a project-cold build
or a source edit. Seed copies and validation are outside the timing.

| Measurement | Existing production preparation | Reuse API test-library artifacts |
|---|---:|---:|
| Python module preparation | 14m36s | 3m26s |
| Same five-test check | 1.62s | 1.66s |
| Preparation plus test check | 14m38s | 3m27s |
| Cargo compilation messages | 288 | 14 |
| Preparation builders | 108 | 3 |
| Module builder occupancy (minutes) | 18.72 | 3.42 |

Module preparation improves by **76.53%** in this pair; preparation plus the
five-test check improves by **76.38%**. Exact executed test names match and all
five pass. Runtime checks compile nothing; exact unchanged repeats start zero
builders. Module runtime-closure NAR content is slightly smaller, 188.46 → 186.88
MB across seven paths, so this does not expand that final closure. These sizes
are local logical content, not transferred bytes. Builder occupancy sums local
builder wall intervals, not CPU or billed time. No prior-revision incremental
state is used for this Python variant.

The earlier unpaired candidate observation was 245.37 seconds. It is retained
as history, not combined with this pair. The older 914.17-second mixed runtime
phase also ran ordinary tests and is not an isolated Python baseline.

**This experiment remains unintegrated.** It is not part of `5b87a850e` or the three
approved NixCI runs, and does not establish full CI completion savings. In the
measured main cold candidate, the Python module was already ready roughly 170
seconds before the final integration group completed; an isolated module gain
cannot simply be subtracted from that total. Full required-check latency and
FeynKit's Python layout/cache boundaries still need validation before integration.

Evidence: `/tmp/ci-dependency-reuse-audit/cold-followup/python-pair/`, including
`plan.json`, `comparison.json`, exact five-test inventories and redacted phase
logs. The prior observation remains under `python-graph-probe/timing-v1/` in the
main measurement directory.

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
That earlier candidate still rebuilt Symbolica and several unchanged workspace
libraries and lacked previous-revision incremental compiler state. Each affected
package retained separate library and test compilation.

## Earlier FeynKit and rejected approaches

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

## Earlier incremental-state experiment

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

At the time of the earlier narrow experiment, no additional NixCI suite had
been pushed and no cold-build or remote worker-time improvement was claimed.
The integrated follow-up measurements are reported above.

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
separately to each affected library and test producer. The integrated results
above measure remote transfer overhead. The completed cold pairs are reported
above; eliminating remaining artifact-merge overhead is a separate experiment.

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

The CI log reporter's 24 tests pass, including validated retry-check association. All three approved main NixCI
suites have been submitted. The unchanged and independent source-edit commits
have the exact locally validated archive roots; their remote evidence is described
above.
FeynKit timing and inventory evidence is under
`/tmp/ci-dependency-reuse-audit/integrated/benchmark/`, including
`integrated-result.json` and `inventory-comparison-v1.json`.
