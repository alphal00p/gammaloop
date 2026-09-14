# CI artifact preparation

## Implemented follow-up, 10 September 2026

The changes are on `codex/ci-cargo-parity-libraries`, with FeynKit validated in the
isolated `codex/ci-artifact-feynkit-local` integration branch. The shared FeynKit
branch and application code are unchanged. These are local measurements; no new
NixCI suite has been launched.

The implementation removes repeated archive work while retaining each package's
compilation boundary. A single artifact input is inherited directly; multiple
inputs retain only independent dependency contexts. FeynKit's separate UFO input
is preserved. Compilation-artifact producers no longer embed nonexistent output
`/lib` paths in their helper binaries, which previously retained obsolete archives.
Actual native-library dependencies are preserved and runtime checks pass.

NixCI publication wrappers now depend on their artifacts, without a commit marker
or eager compiler-state links. All 34 selected wrappers, including the new static
check dependency producer, have identical derivations across a revision-only
change. None requests an `incremental` output. Optional `ci-compiler-state` exposes
only the scheduled producers' pinned seeds, independently of the automatic jobs:

```sh
nix build .#ci-compiler-state --out-link result-ci-compiler-state
```

The result link keeps the prepared seeds available locally. Publishing them still
requires the configured cache upload mechanism; the command alone does not upload
them. An existing local build verified that the collection adds just one symlink
farm, with no otherwise-unused producer compilation. Remote seed availability,
substitution, worker minutes and final-success timing still need a NixCI check.
The earlier 4m35s tail has **not** been remeasured remotely.

Python preparation inherits the API's test-library artifacts, retaining its
Python-specific features, ABI, interpreter and optional FeynKit dependencies.
Its own pinned compiler state is reusable after an edit. Packaging inherits these
artifacts without compiling again. Cargo check, Clippy and doctests now share a
dependency producer built with their exact workspace features and compilation
modes. Ordinary test producers retain their narrower inputs. This does not make
workspace Clippy/doctest source inputs granular.

### Correctness and an archive bug found during validation

The complete prepared revisions pass 1,633 tests and 43 doctests on main, and
1,834 tests and 62 doctests on FeynKit. There are 155/127 skipped nextest cases and
20 ignored doctests per layout. Package inventories, names, target kinds, ignored
flags and filters match the earlier measurements exactly. Python contributes five
passing cases in each layout. Cargo check, Clippy, formatting and generated
workspace/CI configuration checks pass. FeynKit's boundary guard and self-test pass;
its eighth group and nine extracted crates remain intact. Its existing Python-stub
workflow is preserved, but stub generation was not rerun in these local suites.

The first FeynKit changed-source candidate compiled the Python extension twice.
A Cargo fingerprint trace identified missing `.rmeta` files: rustc restored newly
named metadata from its query cache with epoch-1 timestamps, so Crane's normal
mtime-based delta omitted them. The existing seeded-artifact helper now records
inherited filenames and touches newly created files before publication. It retains
compact deltas; it does not repack every inherited file or combine query caches.
Both layouts now perform zero Cargo compilation during final Python packaging.
A negative control still rejects a real compile error with source mtime set to
zero, confirming that restored state does not bypass changed source.

The derivation-boundary probes retain GammaLoop-only isolation for all nine
extracted libraries, CFF/model reverse dependencies, fixture-only test changes,
and broad manifest/lockfile invalidation. One pre-existing assertion still expects
a `feynkit-py` edit to affect only that package; the actual graph also invalidates
its API consumer. That discrepancy predates these changes. The assertion was not
changed and its failure remains in the evidence.

An earlier main baseline timing attempt had two intermittent state-file failures.
The identical tests passed on an unchanged retry. Several feyngen cases share and
clean a workspace directory; the logs preserve this possible race for separate
investigation. No failing tests, filters or retry policy were changed. The failed
attempt is not presented as a successful full-suite measurement.

### Local timings

| Same-layout source edit, all runtime groups | Before | After | Change |
|---|---:|---:|---:|
| Main: completion including evaluation/setup | 10m12s | 4m53s | −52.0% |
| FeynKit: completion including evaluation/setup | 9m58s | 5m00s | −49.8% |
| Main: artifact preparation alone | 138.32s | 93.84s | −32.2% |
| FeynKit: artifact preparation alone | 155.72s | 97.58s | −37.3% |

The main pair sums 15.63 versus 8.68 builder-minutes (−44.5%); FeynKit sums
14.98 versus 8.70 minutes (−41.9%). Main's Python module
is ready 449 versus 101 seconds into the runtime phase; the final candidate Python
packaging step compiles nothing. The candidate's integration tests finish last.
Python preparation still has 14/20 Cargo compilation messages on main/FeynKit,
using its distinct feature context and restored query caches. This does not
invalidate the ordinary per-crate FeynKit test artifacts.

The fresh main artifact-preparation breakdown is:

| Work on the completion path | Before | After |
|---|---:|---:|
| Cargo commands | 49.35s | 50.20s |
| Archive merging | 60.19s | 12.89s |
| Other builder work, including restore/export | 27.64s | 29.76s |
| Waiting and command completion | 1.14s | 0.99s |

This follows the last completed direct input through the build graph and does not
sum overlapping branches. Most of this preparation gain removes filesystem work,
not compiler work. The older 44.18s native Cargo / 185.90s Nix observations below
are historical, not the baseline of this fresh pair.

Warm reruns force test execution while reusing the same artifacts. Main passes
1,633 cases with zero Cargo compilation: 0.17s artifact preparation and 167.30s
runtime phase. FeynKit passes 1,834 cases with zero Cargo compilation: 0.18s and
162.38s respectively. Evaluation/setup adds approximately 17s/26s. These warm
results execute tests; the cold experiment's exact repeats instead reuse successful
results and run zero builders.

These are single sequential observations on a shared host with eight selected
CPUs, eight Nix cores and four jobs. Each layout compares the same application
source between implementations; source edits start independently from their clean
pinned revision. Test execution is forced for the warm measurement. A cached
successful result that executes no tests is reported separately. The changed-source
rows include all runtime groups and Python, but exclude Clippy/doctests.

| Project-cold full scope | Earlier controlled baseline | New candidate |
|---|---:|---:|
| Main: all required checks | 44m01s | 38m21s (−12.9%) |
| Main: combined builder occupancy | 126.57 minutes | 106.18 minutes (−16.1%) |
| FeynKit: all required checks | 48m04s | 42m13s (−12.2%) |
| FeynKit: combined builder occupancy | 163.70 minutes | 141.86 minutes (−13.3%) |

All seven main runtime groups finish earlier. The static checks have a cold-start
trade-off: Clippy finishes at 31m20s instead of 21m40s, and doctests at 37m13s
instead of 29m20s, after waiting for the shared dependency producer. Both still
finish before the final integration group. A separate source-edit pair runs Clippy and doctests concurrently, with the
dependency caches already prepared. Completion improves from 9m00s to 8m14s
(8.6%), and combined builder time falls 11.0%. Clippy itself takes 372s versus
318s; doctests take 540s versus 493s. All 43 doctests have identical names and
outcomes. This is a modest warm benefit with a cold-start trade-off; these checks
still compile the broad workspace inputs and do not reuse per-crate query caches.
The static timings exclude the ordinary runtime groups.

Main's exact repeat takes 0.26s with no builders. All 19 roots succeed; 1,633 test
cases and 43 doctests match the baseline, including exact names and statuses.
No archive step compiles anything. FeynKit passes all 21 roots, with the same
1,834 test cases and 62 doctests, including exact names and statuses. Its exact
repeat takes 0.27s with no builders. All eight FeynKit runtime groups, Clippy and
doctests finish earlier than in its baseline; the static-check cold regression
above is specific to the main observation.

The main cold completion path contains 24m30s of Cargo commands, 1m28s of archive
merges, 1m51s of other builder work, 3m22s in the final runtime group, and 7m10s of
waiting/command completion. The earlier path spent 6m25s in merges. These are
non-overlapping wall intervals along one completion path, not aggregate CPU usage;
other required checks execute in parallel.

The cold candidate starts in a disposable store containing no project outputs or
compiler state. Its complete external seed closure and baseline derivation
identities match the earlier same-day controlled measurement. The comparison uses
that historical baseline, not a fresh contemporaneous baseline run. Both use CPUs
8–15, eight cores, four jobs, sandboxing, and no substituters, remote builders or
upload hooks. Shared memory/storage and host activity are not isolated.

### Stored closures

| Main output | Previous closure | New closure |
|---|---:|---:|
| GammaLoop test-library artifacts | 3.44 GiB | 1.11 GiB |
| GammaLoop test-binary artifacts | 1.78 GiB | 1.22 GiB |
| API test-library artifacts | 1.70 GiB | 1.13 GiB |
| Integration test-binary artifacts | 1.94 GiB | 1.38 GiB |

The GammaLoop library closure falls 67.6%. Its own compressed archive remains
about 36.5 MB; most savings remove retained ancestor copies. Compiler-artifact
content within that closure falls 82.2%, but this is **not** a whole-suite 80%
transfer reduction. These are NAR storage/closure sizes, not network bytes or CHF
cost estimates. Inspected ordinary closures contain no compiler-state outputs.

### Evidence and rollout

Machine-readable observations are in
[ci-artifact-preparation-results.json](ci-artifact-preparation-results.json).
Raw phases, redacted activity logs, inventories, derivations and negative controls
are under `/tmp/ci-artifact-implementation`; frozen cold inputs are under its
`cold/` directory. The final source-edit snapshots are independent of the working
application checkout. Commands and source/CI hashes are recorded with each run.

The existing report collector retains the earlier NixCI anomaly evidence. No new
remote observations have been made, and the local missing-metadata bug and flaky
state-file failures are not attributed to NixCI. Stable publication is locally
verified, but remote cache retrieval must be tested before claiming the late-job
and transfer goals are met. Synchronous dependency discovery remains integrated;
Syd's local post-build upload hook still needs deployment-specific setup. No new
paid suite is authorized by this follow-up.

The investigation below records the earlier implementation and diagnostics; its
"current" references describe that historical snapshot.

## Earlier investigation

Follow-up to the controlled Cargo/Nix and NixCI measurements in
[ci-cargo-parity.md](ci-cargo-parity.md). Application and CI source files were not
changed in this investigation. Experiments ran in the existing disposable local
store with remote builders, substituters and upload hooks disabled. No NixCI suite
was launched. Detailed scratch evidence is in /tmp/ci-artifact-followup.

## The 44-second versus 186-second preparation gap

The native Cargo run took 44.1777 seconds for compilation/artifact preparation.
Following the last completed direct input recursively through the recorded Nix
build gives this non-overlapping partition of its 185.8983 seconds:

| Work on the completion path | Seconds |
| --- | ---: |
| Cargo commands | 62.55 |
| Archive merges | 82.17 |
| Other builder work | 39.02 |
| Ready waiting and command completion | 2.16 |

The Cargo chain is GammaLoop library 13.63s, API library 13.92s, integration library
13.88s, integration test binaries 20.78s and nextest archive command 0.34s.
Parallel GammaLoop/API test-binary builders are not added to this chain.

The 39.02s includes 20.51s in patch phases (including Crane's artifact restoration),
8.81s in build phases outside reported Cargo time (including compiler-state
restoration), 4.87s installing outputs, 3.47s fixing up outputs and other small
phases. These are wall intervals, not CPU time. The separate 26.47s evaluation and
definition setup and 182.31s test execution are outside the 185.90s build phase.

Removing the observed 82.17s of merges with every other duration fixed would leave
about 103.73s. This is an accounting counterfactual, not a measured candidate
speedup; it still exceeds native Cargo. A changed build graph can alter contention,
restoration costs, cache closures and Cargo freshness.

## Full-content checks of redundant merges

The existing binary merge combines the common dependency cache with a library
archive that already contains those dependencies. Full restoration comparisons
checked SHA-256 for every regular file, path/type/mode/symlink metadata, writable
permissions and epoch-1 file timestamps, with Cargo locks cleaned on both sides.

| Layout and consumer | Comparison | Result |
| --- | --- | --- |
| Main: GammaLoop, API, integration test binaries | Current merge vs its existing library archive | Equal: 4,641 / 4,686 / 4,733 files |
| Main: API library | Current merge vs GammaLoop library archive | Equal: 4,641 files |
| Main: integration library | Current merge vs API library archive | Equal: 4,686 files |
| FeynKit: GammaLoop, API, integration test binaries | Current merge vs its existing library archive | Equal: 4,878 / 4,954 / 5,001 files |
| FeynKit: integration library | Current merge vs API library archive | Equal: 4,954 files |
| FeynKit: API library | Current merge vs GammaLoop library alone | Missing 31 UFO/feature-anchor files |
| FeynKit: API library | Current merge vs UFO then GammaLoop library archives | Equal: 4,909 files |

The FeynKit difference includes the feynkit-ufo rlib/rmeta, fingerprints and its
feature anchor. The optimization must retain that independent dependency.
A fixed “use GammaLoop for every API input” rule is incorrect.

All three binary merges are redundant in the inspected snapshots of both layouts.
Their recorded main durations total 41.53 worker-seconds; only the integration
binary merge's 13.43s lies on the original final completion path. The two main
library merges contribute another 68.74s on that path. The restoration comparisons
are correctness probes, not an end-to-end speed benchmark; their timing order
favors the second restoration through page-cache warming. General source-edit,
feature, manifest and fixture boundaries still need validation after implementation.

## New finding: unused ELF search paths retain large cache closures

Blindly using the existing GammaLoop library archive in place of its binary-input
merge would increase that input's stored closure from 1,797,614,944 to
3,688,328,960 NAR bytes. The library delta references a merged ancestor archive
whose compressed contents accidentally retain older artifact outputs.

The exact pinned stdenv adds an output's own lib directory to the linker search
path unless NIX_NO_SELF_RPATH=1:
`/nix/store/sg7xjcivmh51ckjdr7bjik3lghd8kjax-stdenv-linux/setup:513`.
These dependency outputs contain target.tar.zst and optional .prev links, with no
lib directory. Nevertheless compiled build scripts and shared objects contain
RUNPATH entries such as an old artifact output followed by /lib. Some store hashes
remain visible in compressed data, and Nix records those references.

The diagnostic examined this existing GammaLoop input:
`/nix/store/m7ljznqc675vmbj8i1966gsgyfx9pz3n-gammaloop-crate-test-dependencies-be74a11b6e47d216-inputs`.

It removed only five known nonexistent artifact /lib search paths from 19 ELF
files, checked that DT_NEEDED entries did not change, and repacked the archive
with the existing Crane hook and compression settings.

| Metric | Original | Diagnostic |
| --- | ---: | ---: |
| Compressed target archive | 491,472,970 bytes | 491,471,928 bytes |
| Total stored closure | 3,651,829,808 bytes | 1,205,876,816 bytes |

The closure reduction is approximately 67%, while archive size barely changes.
This is a local diagnostic, not a measured network transfer or full-suite speedup.
It identifies removable dependency retention before further compression tuning.

The appropriate implementation experiment is to prevent self-RPATHs specifically
in compilation-artifact producers, then regenerate compatible seeds and inspect
their output closures. Preserve actual native-library dependencies explicitly:
compression can also hide legitimate references from Nix's scanner. The current
diagnostic does not establish all standalone runtime closure requirements.
Do not disable reference checking globally or remove all RPATH entries.
Old pinned artifacts may retain old paths until rebuilt.

## Publication and final success

A new same-source, different-revision evaluation confirms:

- Current selected publication wrappers: 33; all 33 change on an empty commit.
- Underlying compilation derivations: zero changes.
- Wrappers with the revision marker removed: zero changes.
- Removing compiler-state links from publication wrapper commands reduces requested
  incremental outputs from 30 to zero, while retaining their ordinary outputs.
- The earlier GammaLoop-edit evaluation changes only seven stable publication
  wrappers, preserving the other 26.

These are derivation-identity results. NixCI has not run either publication variant.
The current wrappers were introduced because memoized top-level success did not
always ensure that downstream workers could retrieve the closure. Stable wrappers
and separate seed publication therefore need a remote cache-availability check.
Compiler state should be available when a changed build actually needs its pinned
seed; dropping all seed publication without replacement could move rebuilding into
those workers.

The current configuration selects 52 jobs. In the measured warm suite, ten
preparation jobs finish after required completion at 8m56; final deployment reports
success at 13m31, a 275-second tail. Those jobs have 496.08 total worker-seconds and
9.29 GB of reported restored content across their complete lifetimes. Some work
precedes the cutoff, so those totals are not pure tail costs or guaranteed savings.

Removing the unnecessary selected work is the useful change. Reporting final
success sooner while leaving the same background jobs running would not reduce
that work. Preserve every required check and the prerequisites used by impure
runtime runners. Syd's [PR 104](https://github.com/alphal00p/gammaloop/pull/104)
provides the direction of automatic discovery with only otherwise invisible
ordering edges; retain our explicit selection of intended outputs when evaluating it.

## Python, Clippy and doctests

The existing controlled main Python pair already showed preparation falling from
875.96s to 205.63s by reusing ordinary test-library artifacts. Builders fell from
108 to three, and Cargo compilation messages from 288 to 14, with the same five
passing test cases and unchanged Python features/ABI/interpreter.

The remaining 14 messages include pyo3 components and several workspace crates.
Sharing an archive does not make every Python compilation context identical.
This pair started with ordinary test artifacts available; it is not a full-suite
or changed-source timing. The current main cold suite had the Python module ready
about 170 seconds before final completion. FeynKit's Python behavior, changed-source
reuse and full-suite impact remain to validate before integration.

Static checks have a different mismatch: the common prebuild invokes Cargo build
through a dependency anchor excluding linnet-py/spynso3, whereas Clippy and doctests
use --workspace and all workspace sources. The complete workspace adds dependencies
and feature combinations that the common input did not prepare.

In the main cold run, Clippy took 519.72s with 45 Compiling messages; doctests took
835.96s with 73. All 45 Clippy crate names also occur in doctest compilation, but
matching names do not establish identical compiler work. Clippy checking metadata
and rustdoc/build artifacts must be distinguished. Both checks discard newly
produced Cargo artifacts, so they cannot share those new results with one another.

Prepare dependency inputs for their exact features and compilation modes, then
measure whether sharing or splitting checks advances required completion.
[Crane supports inherited Cargo artifacts for Clippy](https://crane.dev/API.html#cranelibcargoclippy),
but merely pointing it at a test archive does not prove compatible fingerprints.
Keep all existing lint targets, doctest cases and Python/stub coverage.

## Next implementation and validation

1. Prevent unnecessary artifact self-RPATHs and inspect references/native runtimes.
2. Bypass the proven redundant binary merges; reduce library inputs according to
   the actual dependency contexts, retaining FeynKit's UFO branch.
3. Measure independent changed-source builds from newly prepared seeds on both
   layouts, preserving inventories and cache-boundary probes.
4. Validate Python reuse in FeynKit and in complete/changed-source suites.
5. Trial stable publication and seed handling remotely, preserving required-check
   coverage and watching for cache misses, duplicate builds and late preparation.
6. Prepare matching static-check inputs as a separately measured experiment.

The original remote run allowance is exhausted. No additional push or paid suite
is included in this investigation. No production CI changes or project tests were
modified. One initial diagnostic used an uncached default Python version and failed
during dependency realization; it was rerun with the already-cached Python 3.13.
The failed log is retained and no project compilation timing includes that attempt.
