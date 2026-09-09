# Closing the gap with native Cargo

The objective is comparable **edit-to-test completion**, including Nix evaluation,
artifact restoration and export, compilation, and execution of the same tests.
A cached successful result that does not execute tests is a separate measurement.
The current change improves main, but **does not achieve Cargo parity**.

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

## What changes

The existing feature anchors align the regex dependency flags used by GammaLoop,
its API and integration tests, on both the target and build-script sides. Cargo's
full recursive GammaLoop compilation-unit signatures now match in those contexts.
Without the build-script alignment, the Symbolica dependency graph still differs.
Only contexts that depend on GammaLoop receive the additional anchor dependencies;
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
export nor the lock-file failure is evidence of a NixCI service fault. The pinned
Crane inheritance hook removes `.cargo-lock`, but does not remove the
`.cargo-build-lock` used by this Cargo version.

## Validation and evidence

Clippy, doctest, formatting and workspace-graph derivations are unchanged from the
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
[the native comparison](ci-main-native-comparison.md). The next experiment must preserve previous-revision compiler state separately for
each affected library and test producer, while returning identical artifacts for
unchanged packages. It must combine that state with the current dependency outputs,
keep fixture-only changes out of production libraries, and reject incompatible
toolchain, manifest, lockfile or feature contexts. The transfer and export cost must
be included in the same full edit-to-test comparison. Simply enabling incremental
compilation inside a fresh Nix build directory does not provide this state.
