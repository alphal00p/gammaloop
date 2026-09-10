# Artifact preparation investigation

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
