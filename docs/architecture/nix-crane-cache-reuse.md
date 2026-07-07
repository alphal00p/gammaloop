# Nix/Crane Cargo artifact reuse

This note records the cache-reuse audit for the NixCI Rust build graph.

## Sources checked

- Crane FAQ: <https://crane.dev/faq/constant-rebuilds.html>
- Crane quick-start example: <https://crane.dev/examples/quick-start.html>
- Crane API reference: <https://crane.dev/API.html>
- Crane issue 180 and comments: <https://github.com/ipetkov/crane/issues/180>

The important Crane contracts are:

- Nix rebuilds a derivation when any direct input changes, including files in
  `src`. Source filtering and narrow source sets prevent unrelated invalidation.
- `cargoArtifacts` is an existing Cargo `target` directory. Crane inherits it
  after patching and before build hooks, so Cargo can mark already-built units as
  fresh.
- The inherited artifact only avoids compilation when Cargo asks for the same
  unit: same profile, target kind, package graph, and feature set.
- `doNotLinkInheritedArtifacts = true` does not disable reuse. It asks Crane to
  deep-copy inherited artifacts instead of symlinking reusable `.rlib`/`.rmeta`
  files. That can be slower, but it is not the reason Cargo recompiles.
- Crane issue 180 has no built-in `buildWorkspace` answer. The maintainer points
  at a Guppy/topological chain of package derivations with filtered sources and
  `cargoArtifacts` passed from dependencies. Later comments describe the same
  pattern: per-crate dependency artifacts trade more derivations for better
  code-change reuse.

## Tested MWEs

The MWE workspace was created under `/tmp/crane-reuse-mwe` with crates `a` and
`b`, where `b` depends on `a`. Commands used the repository's pinned Rust:

```bash
PATH=/nix/store/qqq8y0xm1n3y9wy8a7hb05clyr9dbcb9-rust-stable-with-components-2026-04-16/bin:$PATH
```

The source files were:

```toml
# Cargo.toml
[workspace]
members = ["a", "b"]
resolver = "2"

[profile.ci]
inherits = "dev"
debug = false
```

```toml
# a/Cargo.toml
[package]
name = "a"
version = "0.1.0"
edition = "2021"

[features]
extra = []

[lib]
path = "src/lib.rs"
```

```rust
// a/src/lib.rs
pub fn value() -> u32 {
    1
}
```

```toml
# b/Cargo.toml
[package]
name = "b"
version = "0.1.0"
edition = "2021"

[dependencies]
a = { path = "../a" }

[lib]
path = "src/lib.rs"
```

```rust
// b/src/lib.rs
pub fn value() -> u32 {
    a::value() + 1
}
```

### Source path alone was not enough to force rebuilds

Build `a` in one copy of the workspace, copy `target`, then build `b` in an
identical workspace at another path:

```bash
CARGO_INCREMENTAL=0 cargo build \
  --manifest-path /tmp/crane-reuse-mwe/source-a/Cargo.toml \
  --target-dir /tmp/crane-reuse-mwe/target-a \
  --profile ci -p a -vv

cp -a /tmp/crane-reuse-mwe/target-a /tmp/crane-reuse-mwe/target-b

CARGO_INCREMENTAL=0 CARGO_LOG=cargo::core::compiler::fingerprint=info cargo build \
  --manifest-path /tmp/crane-reuse-mwe/source-b/Cargo.toml \
  --target-dir /tmp/crane-reuse-mwe/target-b \
  --profile ci -p b -vv
```

Observed result: Cargo printed `Fresh a` and only compiled `b`. The absolute
workspace path alone did not explain the reuse failure.

### Feature mismatch does force a second crate variant

Add feature `extra` to crate `a`, build `a` with that feature, then consume the
copied target from `b` without the feature:

```bash
CARGO_INCREMENTAL=0 cargo build \
  --manifest-path /tmp/crane-reuse-mwe/source-a/Cargo.toml \
  --target-dir /tmp/crane-reuse-mwe/target-feature-a \
  --profile ci -p a --features a/extra -vv

cp -a /tmp/crane-reuse-mwe/target-feature-a /tmp/crane-reuse-mwe/target-feature-b-narrow

CARGO_INCREMENTAL=0 CARGO_LOG=cargo::core::compiler::fingerprint=info cargo build \
  --manifest-path /tmp/crane-reuse-mwe/source-c/Cargo.toml \
  --target-dir /tmp/crane-reuse-mwe/target-feature-b-narrow \
  --profile ci -p b -vv
```

Observed result: Cargo compiled `a` again. The log looked for a different
fingerprint (`a-0a0f...`) than the already-built feature-enabled artifact
(`a-5b84...`).

Running the same consumer with matching qualified features:

```bash
cp -a /tmp/crane-reuse-mwe/target-feature-a /tmp/crane-reuse-mwe/target-feature-b-wide

CARGO_INCREMENTAL=0 CARGO_LOG=cargo::core::compiler::fingerprint=info cargo build \
  --manifest-path /tmp/crane-reuse-mwe/source-c/Cargo.toml \
  --target-dir /tmp/crane-reuse-mwe/target-feature-b-wide \
  --profile ci -p b --features a/extra -vv
```

Observed result: Cargo printed `Fresh a` and only compiled `b`.

## Local audit findings

The repo uses Hakari and Guppy to create a cacheable workspace graph, but Cargo
artifact reuse depends on more than the Nix graph shape. Cargo reuses an
inherited artifact only if the requested unit has the same package, profile,
target kind, source freshness metadata, feature set, and dependency hashes.

### Artifact merge and freshness fixes

The artifact merge layer had to preserve Crane's previous-artifact chain and
keep inherited targets writable:

```bash
if [ -e "$artifact.prev" ] || [ -L "$artifact.prev" ]; then
  unpack_artifact "$(realpath "$artifact.prev")"
fi
zstd -d "$artifact" --stdout | tar --no-same-permissions -x -C "$out/target"
rsync -a --chmod=u+w "$artifact/" "$out/target/"
```

The merge layer now materializes even a single inherited artifact. This is
intentional: Crane's single-artifact path does not recursively unpack the
`target.tar.zst.prev` chain. Returning the single artifact directly made later
derivations lose the root prebuild and rebuild the heavy dependency stack.

Dependency-only archives also strip dummy workspace artifacts after compiling
third-party dependencies. The strip logic intentionally ignores binary target
names: `clinnet` has a binary named `linnet`, and stripping by every target name
would delete `liblinnet-*`.

The Hakari package needs a real build artifact, not just a dependency-only
archive. Otherwise downstream Cargo fingerprints point at a different
`gammaloop-workspace-hack` unit. The flake now builds:

- `workspaceHackDependencyArtifacts`: third-party seed for the hack.
- `workspaceHackBuildArtifacts`: real `cargoBuild` of the hack from that seed.

The hack build script also needed stable source freshness metadata. Cargo records
one local-file marker for build-script freshness; across filtered Nix sources the
marker changed between `Cargo.toml`, `build.rs`, and `src/lib.rs`, which made
downstream fingerprints stale. The flake normalizes those timestamps:

```bash
touch -d @0 crates/gammaloop-workspace-hack/Cargo.toml
touch -d @0 crates/gammaloop-workspace-hack/src/lib.rs
touch -d @1 crates/gammaloop-workspace-hack/build.rs
```

One bug in this layer was that scripts passed to Crane's `mkDummySrc` wrote
paths such as `crates/foo/src/lib.rs` directly. `mkDummySrc` builds its output
under `$out`, so those writes missed the dummy source. Dummy-source edits now
write to `$out/...`; post-patch edits in real package builds still use
workspace-relative paths.

Dependency-only sources now restore real source directories for already-built
workspace dependencies, while keeping the current package's targets dummy. This
lets Cargo reuse upstream workspace artifacts in final package builds, and keeps
the dependency archive stable across ordinary source edits to the current
package.

The Cargo CLI feature set is intentionally narrower than the source/artifact
closure. Nix derivations keep resolved workspace dependency sources and
artifacts so Cargo can validate fingerprints, but `--features` only names the
selected package plus direct workspace dependencies. Cargo rejects qualified
features for transitive packages that are not selected and are not direct
dependencies. The `65e5d86e` NixCI run failed this way for
`crate-deps-gammaloop-tracing-filter` and `crate-deps-spenso-hep-lib`.

### What is fixed locally

`crate-spenso` is the clean sentinel for normal target-library reuse because it
consumes `linnet` and the workspace hack. After the fixes above:

```text
cargo build --profile ci-optim --locked -p gammaloop-workspace-hack -p spenso \
  --features linnet/bincode,linnet/serde,linnet/symbolica,spenso-macros/shadowing,spenso/shadowing,symbolica/tracing_max_level_info

Compiling spenso v0.6.0
```

The final `crate-spenso-build` phase no longer compiles `linnet` or
`gammaloop-workspace-hack`, and it no longer compiles `spenso-macros` in the
final package phase. The proc-macro unit is produced by the dependency archive
for the consumer context.

The same final-build reuse pattern was observed for narrower consumers before
the `gammalooprs` check:

```text
crate-idenso-build: compiled idenso only
crate-spenso-hep-lib-build: compiled spenso-hep-lib only
crate-gammaloop-tracing-filter-build: compiled gammaloop-tracing-filter only
```

The two feature-selection failures from the `65e5d86e` NixCI run now build
locally. Their dependency artifacts compile the required upstream workspace
crates in the consumer context, and their final package derivations compile only
the selected package:

```text
crate-deps-gammaloop-tracing-filter:
  Compiling gammaloop-tracing-filter
  Compiling gammaloop-tracing-filter-macros
  Compiling linnet
  Compiling spenso

crate-gammaloop-tracing-filter-build:
  Compiling gammaloop-tracing-filter

crate-deps-spenso-hep-lib:
  Compiling spenso-hep-lib
  Compiling linnet
  Compiling spenso
  Compiling idenso

crate-spenso-hep-lib-build:
  Compiling spenso-hep-lib
```

The root `crate-gammalooprs-build` phase is the broadest normal target-library
reuse check. After adding the final feature anchors, a local build produced:

```text
cargo build --profile ci-optim --locked -p gammaloop-workspace-hack -p gammalooprs \
  --features gammaloop-tracing-filter/clap,gammaloop-tracing-filter/symbolica,idenso/bincode,idenso/reference-cases,linnet/bincode,linnet/serde,linnet/symbolica,spenso/shadowing,symbolica/tracing_max_level_info

Compiling gammalooprs v0.3.3
Finished `ci-optim` profile [optimized] target(s) in 1m 39s
```

That final root phase did not compile `vakint`, `linnet`, `spenso`, `idenso`,
`spenso-hep-lib`, `gammaloop-tracing-filter`, or
`gammaloop-workspace-hack`. Those crates were consumed from earlier workspace
artifact derivations.

The source-filter cache boundary was checked with temporary code-only comments:

```text
spenso code edit:
  crate-deps-spenso stayed /nix/store/207kj3kab5j87n07i8kims0is6wkv3rq-gammaloop-crate-spenso-deps-0.1.0.drv
  crate-spenso changed to /nix/store/vs6a6mckr8d9yd8w6r8mjli8c8a5xdsb-gammaloop-crate-spenso-build-0.1.0.drv

gammalooprs code edit:
  crate-deps-gammalooprs stayed /nix/store/xs8jqfllizzc6bbjyv3qj2li05r7zs2m-gammaloop-crate-gammalooprs-deps-0.1.0.drv
  crate-gammalooprs changed to /nix/store/0dmnbsxd0cyqirjxr4vjs7z397sv5pnk-gammaloop-crate-gammalooprs-build-0.1.0.drv
```

## Follow-up audit: workspace crate reuse is context-sensitive

The later `linnest -> linnet` investigation found an additional Cargo
constraint: a workspace crate artifact is not universally reusable just because
the package name and feature list look the same. Cargo's unit hash also changes
with the selected package graph.

Observed locally:

```text
crate-linnet-build / standalone linnet context:
  target/ci-optim/.fingerprint/linnet-f122c64f71757531

crate-linnest-deps / linnet as dependency of linnest:
  target/ci-optim/.fingerprint/linnet-58bb076a94f2462f
```

Trying to feed the standalone `crate-linnet` artifact into `crate-linnest`
therefore still made Cargo compile `linnet`. The fix is to make
`crate-deps-<consumer>` the cache boundary for workspace dependencies in that
consumer's graph, and to preserve those workspace dependency artifacts for the
final package derivation.

The final package input also has to go through the custom artifact merge step.
Crane's normal single-archive unpack path does not recursively materialize the
`.prev` chain created by `buildDepsOnlyWithArtifacts`; passing the bare
`crate-deps-linnest` archive caused the final package build to lose inherited
root artifacts and recompile third-party crates. The final package derivations
now receive a merged target directory so the `.prev` chain is expanded before
Crane copies artifacts.

The hack timestamp normalization also has to run after unpack in
`buildDepsOnlyWithArtifacts`. Running it only while constructing a dummy source
is insufficient because Nix store normalization resets file mtimes. Without the
post-unpack normalization, Cargo marked `gammaloop-workspace-hack` dirty with:

```text
dirty: PrecalculatedComponentsChanged { old: "1.000000000s (src/lib.rs)", new: "1.000000000s (build.rs)" }
```

That hack rebuild then made dependent workspace crates, such as `linnet`, dirty
through `UnitDependencyInfoChanged`.

### Current focused verification

The focused `linnest` package path now has this shape:

```text
crate-linnet-deps:
  Compiling cgmath
  Compiling linnet

crate-linnest-deps:
  Compiling linnest   # dummy target package
  Compiling linnet    # dependency in the linnest consumer context

crate-linnest-build:
  Compiling linnest   # final real package only
```

The final `crate-linnest-build` phase no longer compiles `linnet`,
`gammaloop-workspace-hack`, Symbolica, or the rest of the external graph.

The focused test-binary path for `linnet-py` also no longer recompiles `linnet`
in the final test-binary derivation:

```text
crate-test-binaries-linnet-py:
  Compiling gammaloop-workspace-hack
  Compiling linnet-py
  Executable unittests src/lib.rs
  Executable unittests src/bin/stubgen.rs
```

This still compiles the current test package and the Hakari anchor in that test
context. It does not compile upstream `linnet` again.

The nextest archive layer is now a two-step cache boundary:

```text
gammaloop-nextest-binaries-linnet-test-deps:
  cargo test --profile ci-optim --no-run --locked \
    -p clinnet -p gammaloop-workspace-hack -p linnest -p linnet -p linnet-py ...
  Compiling iai-callgrind-runner
  Compiling gammaloop-workspace-hack
  Compiling linnet
  Compiling linnest
  Compiling linnet-py

gammaloop-nextest-binaries-linnet:
  cargo nextest archive ...
  Finished `ci-optim` profile [optimized] target(s) in 0.29s
```

Selecting `gammaloop-workspace-hack` in the nextest package set is required for
Symbolica feature anchoring. Without that selected package, the linnet archive
path requested a narrower Symbolica unit and rebuilt Symbolica. With the anchor,
the test-binary artifact still compiles the selected workspace test units, but
it does not compile Symbolica, numerica, or graphica.

Before the later per-package nextest split, the spenso group used the same
group-level shape:

```text
gammaloop-nextest-binaries-spenso-test-deps:
  cargo test --profile ci-optim --no-run --locked \
    -p gammaloop-workspace-hack -p idenso -p spenso -p spenso-hep-lib \
    -p spenso-macros -p spynso3 ...
  Compiling spenso-macros
  Compiling gammaloop-workspace-hack
  Compiling linnet
  Compiling spenso
  Compiling idenso
  Compiling spenso-hep-lib
  Compiling spynso3

gammaloop-nextest-binaries-spenso:
  cargo nextest archive ...
  Finished `ci-optim` profile [optimized] target(s) in 0.30s
```

An intermediate experiment fed each package's dependency-mode artifact into the
nextest group artifact. That was worse: `spynso3`'s dependency setup compiled
`linnet`, `spenso`, `idenso`, and `spenso-hep-lib` in a separate consumer
context before the spenso test-binary group compiled them again. The nextest
groups were later changed again; see "Follow-up audit: nextest package
boundaries" below for the current per-package shape and the removal of
`spynso3` from the fast `spenso` archive.

The archived nextest runner also has to recreate the source at `/build/source`.
The binaries are compiled with `CARGO_MANIFEST_DIR=/build/source/...`; copying
the source to the runCommand's default working directory made insta miss
snapshots. The runner additionally installs dummy targets for workspace members
outside the narrowed nextest source so `cargo metadata` succeeds.

Local runner checks:

```text
checks.x86_64-linux.gammaloop-nextest-linnet:
  150 tests run: 150 passed

checks.x86_64-linux.gammaloop-nextest-spenso:
  350 tests run: 350 passed, 27 skipped
```

`checks.x86_64-linux.gammaloop-nextest-core` also runs from the archive, but in
the local environment it fails test initialization with `OEM license key not
set`. That is a compile-time license input issue, not a cache-reuse or archive
layout issue.

### Remaining caveats

The root prebuild is intentionally broad and still cold-builds a large unified
external graph, including Python/stubgen-related units and
`numerica`/`graphica`. This is a cache-root cost, not downstream double work.
If that cost is too high, the next optimization should split Python/stubgen
prebuilds from the normal Rust test prebuild rather than trying to reuse
standalone workspace package artifacts across incompatible Cargo unit contexts.

The nextest archive step is no longer a compile boundary. It receives a
group-matched test-binary artifact and only packages fresh binaries.

### Feature-context anchors

Before the final fix, `crate-gammalooprs-build` had already stopped recompiling
most of the workspace stack, but still compiled `vakint`:

```text
gammalooprs
vakint
```

The inherited merged target contained a `vakint` artifact, so this was not a
missing Nix edge. Comparing the two `vakint` fingerprints showed that every
direct dependency hash matched except `eyre`:

```text
vakint archive: eyre dependency hash 1421407417017949601
root archive:   eyre dependency hash 566811204752360852
```

Following the dependency hash down showed that `eyre` was built against a
different `indenter` unit in the root graph:

```text
vakint archive: indenter features ["default"]
root archive:   indenter features ["default", "std"]
```

`arrayvec` had the same shape earlier: the root graph enabled its normal
`default,std` feature set through `vakint`, but the workspace hack did not carry
that feature context while `vakint` stayed outside Hakari traversal.

The workspace hack therefore has manual normal-dependency anchors outside
Hakari's managed block:

```toml
[dependencies.arrayvec]
version = "0.7"
features = ["std"]

[dependencies.indenter]
version = "0.3"
features = ["std"]
```

It also has matching build-dependency anchors because Hakari verification checks
host-side proc-macro feature contexts separately:

```toml
[build-dependencies.arrayvec]
version = "0.7"
features = ["std"]

[build-dependencies.indenter]
version = "0.3"
features = ["std"]
```

After those anchors, the `vakint` artifact and root graph request the same
third-party units, so the root build reuses `vakint`. The build-dependency
anchors also make `cargo hakari verify` pass for the host platform.

Proc-macro packages using Symbolica also need the workspace-hack anchor. Without
that, `crate-deps-spenso-macros` inherited the hack artifact but requested a
narrower `symbolica` unit and rebuilt `symbolica`, `numerica`, and `graphica`.
After allowing proc-macro packages to select the hack anchor, the local
`crate-deps-spenso-macros` build was:

```text
cargo build --profile ci-optim --locked \
  -p gammaloop-workspace-hack -p spenso-macros \
  --features spenso-macros/shadowing,symbolica/tracing_max_level_info

Compiling spenso-macros v0.3.2
Compiling gammaloop-workspace-hack v0.1.0
Finished `ci-optim` profile [optimized] target(s) in 0.78s
```

It no longer rebuilt the Symbolica stack.

### Rejected experiment

Adding `gammalooprs` as a selected feature-context anchor for every dependency
crate did make the Cargo commands closer to the final graph:

```text
cargo build --profile ci-optim --locked \
  -p gammaloop-workspace-hack -p gammalooprs -p linnet \
  --features gammaloop-tracing-filter/clap,gammaloop-tracing-filter/symbolica,idenso/bincode,idenso/reference-cases,linnet/bincode,linnet/serde,linnet/symbolica,spenso/shadowing,symbolica/tracing_max_level_info
```

But it also made every producer compile the dummy `gammalooprs` closure. For
example, `crate-linnet-build` compiled `vakint`, `gammalooprs`,
`spenso-macros`, `spenso`, `idenso`, `gammaloop-tracing-filter`, and
`spenso-hep-lib` in addition to `linnet`. That is worse for CI latency and was
backed out.

## Follow-up audit: archive and dependency-closure reuse

The later integration archive audit found two more places where the Nix graph
looked right but Cargo still rebuilt workspace crates.

### Stripped archives must be flattened after Crane installs them

`buildDepsOnlyWithArtifacts` originally stripped dummy workspace artifacts in
`postBuild`, then tried to recursively materialize the inherited
`target.tar.zst.prev` chain in `postInstall`. That was too early for Crane's
`installCargoArtifactsHook`: the hook writes `target.tar.zst` and `.prev` while
`runHook postInstall` is running, after the shell `postInstall` body has already
checked for the files.

The flatten-and-strip pass now runs in `preFixup`, after Crane has installed the
artifact archive. For stripped dependency archives this leaves a single
self-contained `target.tar.zst` and removes `target.tar.zst.prev`.

The strip script also removes unhashed dependency outputs such as:

```text
target/ci-optim/deps/libgammaloop_api.rlib
target/ci-optim/deps/libgammaloop_api.so
target/ci-optim/deps/gammaloop_api.d
```

Those files were the source of the stale dummy `gammaloop-api` failure. The
merged integration dependency input had previously contained a 15 KiB dummy
`libgammaloop_api.so`, so `gammaloop-integration-tests` failed with missing
`gammaloop_api::commands`, `session`, and `state`. After stripping unhashed
outputs and pre-stripping the current package before package-specific builds,
the same merged input contains the real `libgammaloop_api.so` artifact
(about 28 MiB locally), and the integration dependency artifact compiles.

### Synthetic consumers need artifacts for the same closure they select

The dependency-mode `*-deps-deps` derivations generate a synthetic consumer that
depends on the resolved workspace dependency closure. Before this fix, their
`cargoArtifacts` input only merged direct dependency artifacts. That mismatch
made Cargo rebuild transitive workspace crates in the synthetic consumer layer:

```text
gammaloop-crate-gammaloop-tracing-filter-dependency-deps-deps:
  Compiling linnet
  Compiling spenso

gammaloop-crate-idenso-dependency-deps-deps:
  Compiling gammaloop-workspace-hack
  Compiling linnet
  Compiling spenso
```

The merge inputs now use the same resolved dependency closure as the synthetic
consumer source. Re-running the same targets locally produced only the
synthetic consumer in the deps-of-deps layer:

```text
gammaloop-crate-gammaloop-tracing-filter-dependency-deps-deps:
  Compiling gammaloop-ci-consumer-gammaloop-tracing-filter

gammaloop-crate-idenso-dependency-deps-deps:
  Compiling gammaloop-ci-consumer-idenso
```

The real package layers then compiled only their own package plus the synthetic
consumer:

```text
gammaloop-crate-gammaloop-tracing-filter-dependency-deps:
  Compiling gammaloop-tracing-filter
  Compiling gammaloop-ci-consumer-gammaloop-tracing-filter

gammaloop-crate-idenso-dependency-deps:
  Compiling idenso
  Compiling gammaloop-ci-consumer-idenso
```

The same closure-merge shape was verified for `spenso-hep-lib`, `gammalooprs`,
`gammaloop-api`, and `gammaloop-integration-tests`.

### Plain integration must not depend on the Python module

The plain integration archive was incorrectly using the same Python environment
setup as the Python API archive. That made
`gammaloop-nextest-binaries-integration` depend on `gammaloop-python-module`,
which in turn rebuilt Python-feature workspace crates before ordinary
integration tests could even be archived.

The nextest graph now has a separate `nextestUsesPythonModule` predicate based
on the `python-api-tests` feature. The plain integration archive no longer
builds `gammaloop-python-module`; the Python module is only an input for the
`python-api` archive/check.

### Nextest needs a generated feature anchor

The nextest test-binary layer inherited the right per-crate artifacts but still
rebuilt upstream workspace crates because its Cargo command only enabled
features for selected packages and direct dependencies:

```text
cargo test --profile ci-optim --no-run --locked \
  -p gammaloop-integration-tests -p gammaloop-workspace-hack \
  --features spenso/shadowing,symbolica/tracing_max_level_info

Compiling gammaloop-workspace-hack
Compiling linnet
Compiling spenso
Compiling idenso
Compiling gammaloop-tracing-filter
Compiling spenso-hep-lib
Compiling gammalooprs
```

Cargo rejects qualified features for transitive packages unless they are
selected or direct dependencies of a selected package. Selecting the full
workspace closure would compile too much, so nextest archive builds now add a
generated package under `crates/gammaloop-ci-nextest-<group>-features`. The
anchor depends on the resolved source package closure with the common feature
set and is selected with `-p`. The archive/test commands use `--offline` because
the generated path-only workspace member is not in the checked-in lockfile.

With that anchor, the plain integration test-binary layer became:

```text
cargo test --profile ci-optim --no-run --offline \
  -p gammaloop-ci-nextest-integration-features \
  -p gammaloop-integration-tests \
  -p gammaloop-workspace-hack \
  --features gammaloop-tracing-filter/clap,gammaloop-tracing-filter/symbolica,...

Compiling gammaloop-workspace-hack
Compiling gammaloop-integration-tests
Compiling gammaloop-ci-nextest-integration-features
Finished `ci-optim` profile [optimized] target(s) in 40.39s
```

The archive layer then compiled only the generated feature anchor and packaged
already-built test binaries:

```text
gammaloop-nextest-binaries-integration:
  Compiling gammaloop-ci-nextest-integration-features
  Finished `ci-optim` profile [optimized] target(s) in 0.74s
  Archiving 9 binaries ...
```

Local verification:

```text
nix build --impure .#checks.x86_64-linux.gammaloop-nextest-binaries-integration --no-link -L --print-out-paths
  /nix/store/khpd1v0xp2jrbl0x6ida5275yxyavda1-gammaloop-nextest-binaries-integration-0.1.0

cached rerun:
  /nix/store/khpd1v0xp2jrbl0x6ida5275yxyavda1-gammaloop-nextest-binaries-integration-0.1.0
```

The Python API archive has the same nextest shape after the separate Python
module build:

```text
gammaloop-nextest-binaries-python-api-test-deps:
  Compiling gammaloop-integration-tests
  Compiling gammaloop-workspace-hack
  Compiling gammaloop-ci-nextest-python-api-features

gammaloop-nextest-binaries-python-api:
  Compiling gammaloop-ci-nextest-python-api-features
  Finished `ci-optim` profile [optimized] target(s) in 0.72s
  Archiving 1 binary ...
```

Local verification:

```text
nix build --impure .#checks.x86_64-linux.gammaloop-nextest-binaries-python-api --no-link -L --print-out-paths
  /nix/store/hiv5ils5196i8pis4697353sjzh2bci0-gammaloop-nextest-binaries-python-api-0.1.0

cached rerun:
  /nix/store/hiv5ils5196i8pis4697353sjzh2bci0-gammaloop-nextest-binaries-python-api-0.1.0
```

### Python remains a separate feature family

The Python module path cannot reuse the normal `gammalooprs` and
`gammaloop-tracing-filter` artifacts byte-for-byte. Enabling
`gammaloop-api/python_abi` and `gammaloop-api/pyo3-extension-module` activates
Python/PyO3 features on those crates, so Cargo requests different units.

Observed local Python dependency layer:

```text
gammaloop-api-python-deps:
  Compiling gammalooprs
  Compiling pyo3-build-config
  Compiling pyo3-ffi
  Compiling pyo3-macros-backend
  Compiling pyo3
  Compiling numpy
  Compiling gammaloop-api      # dummy target, stripped from the dependency artifact
  Compiling pyo3-macros
  Compiling gammaloop-tracing-filter
```

The following Python build/package/archive layers reused that Python-feature
artifact family: the real Python build compiled only `gammaloop-api`, the
package install did not print a `Compiling ...` line, and nextest compiled only
the Python API test package plus the anchor. If this remaining Python-feature
workspace compile is too coarse, the follow-up is to split Python-feature
per-crate artifacts rather than moving PyO3 features into the common prebuild.

## Follow-up audit: nextest package boundaries

The nextest archive graph was changed from one test-binary derivation per group
to one test-binary derivation per package, with the final group archive merging
those package artifacts. This makes the Nix graph expose package-sized cache
boundaries instead of hiding all test compilation inside
`gammaloop-nextest-binaries-<group>-test-deps`.

The checked Guppy graph now also records workspace dependency edge features:

```json
"spynso3": {
  "spenso": [
    "python",
    "shadowing"
  ]
}
```

Those incoming features are folded into the common Crane feature set, except for
Python wrapper packages that are intentionally outside the fast Rust cache root:

```nix
workspaceFeatureUnificationExcludedPackages = [
  "linnet-py"
  "spynso3"
];
```

`spynso3` was removed from the fast `spenso` nextest archive. Keeping it there
forces `spenso/python` and `symbolica/python_export` into the same resolver
context as normal Rust `spenso` tests. A local experiment added
`symbolica/python_export` to the workspace-hack anchor; that made the hack
producer cold-build PyO3/numpy/stubgen-related Symbolica dependencies and still
did not stop `spynso3` from rebuilding workspace crates because PyO3 features
were not anchored. The measured cold hack producer was:

```text
gammaloop-crate-gammaloop-workspace-hack-deps:
  cargo build ... --features symbolica/python_export,symbolica/tracing_max_level_info
  Finished in 5m33s
  target archive: 818 MiB -> 218 MiB compressed
```

That route was backed out. The next step for `spynso3`, if it must be tested in
NixCI, should be a separate Python-wrapper graph with its own PyO3/Symbolica
feature family. Putting it back into the fast `spenso` graph would either
reintroduce workspace recompilation or make the root prebuild much heavier.

Local verification after the per-package nextest split:

```text
nix build --impure .#checks.x86_64-linux.gammaloop-guppy-workspace-graph --no-link -L --print-out-paths
  /nix/store/lkqg567a9ig75jjh01l6yh7wqfkvx1fq-gammaloop-guppy-workspace-graph-check

nix build --impure .#checks.x86_64-linux.gammaloop-nextest-binaries-spenso --no-link -L --print-out-paths
  /nix/store/k60aaha4crsvnk1l570lr963w0yaqr4q-gammaloop-nextest-binaries-spenso-0.1.0

nix build --impure .#checks.x86_64-linux.gammaloop-nextest-binaries-integration --no-link -L --print-out-paths
  /nix/store/0imvvji66p0ww16xzlqfr60pbwlfvahm-gammaloop-nextest-binaries-integration-0.1.0

nix build --impure .#checks.x86_64-linux.gammaloop-nextest-binaries-python-api --no-link -L --print-out-paths
  /nix/store/2jyb556w3za039kryg6rg30h534x693l-gammaloop-nextest-binaries-python-api-0.1.0
```

Observed compile shape:

```text
gammaloop-nextest-binaries-integration-gammaloop-integration-tests-test-deps:
  Compiling gammaloop-workspace-hack
  Compiling gammaloop-integration-tests
  Compiling gammaloop-ci-nextest-integration-gammaloop-integration-tests-features

gammaloop-nextest-binaries-python-api-gammaloop-integration-tests-test-deps:
  Compiling gammaloop-integration-tests
  Compiling gammaloop-workspace-hack
  Compiling gammaloop-ci-nextest-python-api-gammaloop-integration-tests-features

gammaloop-nextest-binaries-spenso:
  Compiling gammaloop-ci-nextest-spenso-features
  Archiving 14 binaries
```

So the final archive layers are no longer compile boundaries; they only compile
their generated feature anchor and package already-built binaries.

The per-package nextest prebuilds now consume shared test-support artifacts
instead of normal package artifacts. This matters because `cargo build` and
`cargo test --no-run` produce different Cargo units even under the same named
profile. A normal `crate-deps-*` artifact cannot satisfy test-harness units.

The public `crate-test-binaries-*` package attrs also alias the shared
test-support artifact for that package's test SCC. For example, these two evals
now return the same derivation:

```text
nix eval --impure --raw .#packages.x86_64-linux.crate-test-binaries-spenso-macros.drvPath
nix eval --impure --raw .#packages.x86_64-linux.crate-test-support-spenso.drvPath
  /nix/store/d3lvg9h3pa9yl4yi6pxrcwf7gsiny5my-gammaloop-crate-test-support-spenso-deps-0.1.0.drv
```

That removes the separate exposed per-package test-binary build that could make
NixCI schedule package-specific duplicate work.

The downstream nextest package prebuilds are now clean. They inherit
`crate-test-support-spenso`, then compile only their generated feature-anchor
package before archiving/running existing binaries:

```text
gammaloop-nextest-binaries-spenso-spenso-macros-test-deps:
  cargo test --profile ci-optim --no-run --offline \
    -p gammaloop-ci-nextest-spenso-spenso-macros-features \
    -p gammaloop-workspace-hack -p spenso-macros ...
  Compiling gammaloop-ci-nextest-spenso-spenso-macros-features
  Finished in 0.74s

gammaloop-nextest-binaries-spenso-spenso-hep-lib-test-deps:
  Compiling gammaloop-ci-nextest-spenso-spenso-hep-lib-features
  Finished in 0.70s

gammaloop-nextest-binaries-spenso-idenso-test-deps:
  Compiling gammaloop-ci-nextest-spenso-idenso-features
  Finished in 0.75s

gammaloop-nextest-binaries-spenso-spenso-test-deps:
  Compiling gammaloop-ci-nextest-spenso-spenso-features
  Finished in 0.72s
```

The final archive layer remains clean:

```text
gammaloop-nextest-binaries-spenso:
  Compiling gammaloop-ci-nextest-spenso-features
  Finished in 0.66s
  Archiving 14 binaries
```

The next workspace-crate recompilation failure was one layer earlier, between
shared test-support SCCs. Before the synthetic-consumer fix,
`crate-test-support-linnest` compiled `linnet` for the `linnet/linnest` test
component:

```text
crate-test-support-linnest:
  cargo test --profile ci-optim --no-run --offline \
    -p gammaloop-ci-test-support-linnest-features \
    -p gammaloop-workspace-hack -p linnest -p linnet ...
  Compiling gammaloop-workspace-hack
  Compiling linnet
  Compiling linnest
  Compiling gammaloop-ci-test-support-linnest-features
```

Then `crate-test-support-spenso`, which depended on that artifact, still
compiled `linnet` again in the `spenso/spenso-macros` test component:

```text
crate-test-support-spenso:
  cargo test --profile ci-optim --no-run --offline \
    -p gammaloop-ci-test-support-spenso-features \
    -p gammaloop-workspace-hack -p spenso -p spenso-macros ...
  Compiling spenso-macros
  Compiling linnet
  Compiling spenso
  Compiling gammaloop-ci-test-support-spenso-features
```

That was a failure under the strict "no workspace crate recompilation" goal. It
was no longer hidden downstream in package/archive derivations, but the
fine-grained test-support SCC graph still asked Cargo for a different enough
test context that the inherited `linnet` unit was not reused.

A focused fingerprint check showed that the two `linnet` units have the same
visible `linnet` features:

```text
["bincode", "default", "drawing", "nodestore-vec", "serde", "symbolica"]
```

One difference was the `cgmath` dependency unit: the `linnest` support context
enabled `cgmath/serde`, while the `spenso` support context reached `cgmath`
through `linnet/drawing`. A temporary manual `cgmath/serde` workspace-hack
anchor was tested and then backed out. It made the common root broader but
`crate-test-support-spenso` still compiled `linnet`, so `cgmath` feature
mismatch was not the whole cause. The manifest-side fix is to make
`linnet/drawing` consistently enable `cgmath/serde`, which removes the separate
no-`serde` `cgmath` unit.

The remaining difference was Cargo context: `spenso` wanted a `lib-linnet` unit
for "dependency of a selected test package", while the earlier `linnest`
support derivation had only built the standalone/component context. The support
graph now prebuilds dummy next-consumer contexts and strips those dummy consumer
artifacts afterward. This keeps downstream source code out of the upstream
support source hash, but gives Cargo the exact upstream workspace units later
components request.

Local verification after narrowing the synthetic consumers:

```text
crate-test-support-linnest:
  cargo test ... \
    -p gammaloop-ci-test-support-linnest-features \
    -p gammaloop-workspace-hack -p linnest -p linnet -p spenso -p spenso-macros ...
  Compiling linnet
  Compiling spenso-macros
  Compiling spenso
  Compiling linnest
  target archive: 220 MiB => 66.4 MiB

crate-test-support-spenso:
  cargo test ... \
    -p gammaloop-ci-test-support-spenso-features \
    -p gammaloop-tracing-filter -p gammaloop-workspace-hack -p idenso -p spenso -p spenso-macros ...
  Compiling spenso-macros
  Compiling spenso
  Compiling idenso
  Compiling gammaloop-tracing-filter
  no Compiling linnet
  target archive: 136 MiB => 42.4 MiB

crate-test-support-idenso:
  cargo test ... \
    -p gammaloop-ci-test-support-idenso-features \
    -p gammaloop-workspace-hack -p idenso -p spenso-hep-lib ...
  Compiling idenso
  Compiling spenso-hep-lib
  no Compiling spenso
  no Compiling linnet
```

The exported support archives were checked to ensure dummy consumer artifacts
are stripped. For example, `crate-test-support-spenso` exports `linnet`,
`spenso`, and `spenso-macros` fingerprints, but not `idenso` or
`gammaloop-tracing-filter` fingerprints.

The current conclusion is:

- Normal package outputs reuse workspace crate artifacts in their final package
  layer.
- Test-support SCCs prebuild the next dummy consumer context, so later support
  SCCs reuse upstream workspace crates instead of recompiling them.
- Downstream nextest package prebuilds and archive layers reuse the shared
  test-support artifacts and no longer recompile workspace crates.
- The tradeoff is a larger earlier support artifact. The `linnest` support
  archive grew from roughly `152 MiB => 45 MiB` to `220 MiB => 66 MiB`, but the
  work moved to a NixCI-cacheable dependency layer.

The previous failed shape is kept here for comparison. Feeding `spenso`'s
normal package artifact into `spenso-macros` was insufficient because Cargo
rebuilt `linnet` and `spenso` in the `spenso-macros` test context:

```text
gammaloop-nextest-binaries-spenso-spenso-macros-test-deps:
  Compiling gammaloop-workspace-hack
  Compiling linnet
  Compiling gammaloop-ci-nextest-spenso-spenso-macros-features
  Compiling spenso-macros
  Compiling spenso
```

A follow-up experiment added a test-specific dependency producer with the same
generated anchor name. That producer compiled `linnet`, `spenso`, and
`spenso-macros`, but the final `spenso-macros` test derivation still compiled
the same workspace units again. This was backed out because it added work
without improving reuse. The current hypothesis is that Cargo's test-unit
fingerprint includes enough context from the real current package graph that
the dummy/dependency producer cannot satisfy the final package test request.

## Follow-up audit: package artifacts preserve workspace crates

The `65e5d86e` NixCI run failed before the nextest archive nodes could prove
anything about reuse. The hard failures were:

```text
packages.x86_64-linux.crate-deps-gammaloop-tracing-filter
packages.x86_64-linux.crate-deps-spenso-hep-lib
```

Both failed because the generated command passed transitive package feature
flags to a package that did not directly own those feature names:

```text
cargo build ... -p spenso-hep-lib \
  --features idenso/bincode,idenso/reference-cases,linnet/bincode,linnet/serde,linnet/symbolica,...

error: the package 'spenso-hep-lib' does not contain these features:
  linnet/bincode, linnet/serde, linnet/symbolica
```

The deeper bug was that public `crate-deps-*` outputs used the older
dummy-current-package artifact path. That path compiled with the current
package target dummied and then stripped the current package artifact, so a
downstream package received dependency archives but not the compiled workspace
crate it needed to reuse. The public `crate-deps-*` outputs now point at the
dependency-mode artifacts instead. Those build a generated consumer package and
preserve the real package artifact.

The dependency-mode artifacts now also preserve the package's resolved
workspace dependencies. Without that, the final `crate-*` package derivation
received the current package artifact but not upstream workspace crates, and
Cargo rebuilt `linnet`/`spenso` in the final package build.

Final package commands now request the full resolved workspace feature set and
select the non-proc-macro feature owners. Proc-macro crates are kept as
dependencies, not selected roots, because selecting `spenso-macros` as a root
created a separate Cargo unit and rebuilt it.

Local verification after this change:

```text
nix build --impure .#packages.x86_64-linux.crate-gammalooprs --no-link -L --print-out-paths
  crate-deps layer:
    Compiling gammalooprs
    Compiling gammaloop-ci-consumer-gammalooprs
  final crate layer:
    Finished in 0.32s, no Compiling lines

nix build --impure .#packages.x86_64-linux.crate-spenso-hep-lib --no-link -L --print-out-paths
  crate-deps layer:
    Compiling spenso-hep-lib
    Compiling gammaloop-ci-consumer-spenso-hep-lib
  final crate layer:
    Finished in 0.29s, no Compiling lines

nix build --impure .#packages.x86_64-linux.crate-gammaloop-tracing-filter --no-link -L --print-out-paths
  final crate layer:
    cargo build ... -p gammaloop-tracing-filter -p gammaloop-workspace-hack -p linnet -p spenso ...
    Finished in 0.29s, no Compiling lines

nix build --impure .#checks.x86_64-linux.gammaloop-guppy-workspace-graph --no-link --print-out-paths
  /nix/store/lkqg567a9ig75jjh01l6yh7wqfkvx1fq-gammaloop-guppy-workspace-graph-check
```

The previous remaining double work was in the shared test-support SCCs:
`crate-test-support-spenso` compiled `linnet` after inheriting
`crate-test-support-linnest`. Synthetic next-consumer contexts now move that
work into the earlier cacheable support layer. The public `crate-test-binaries-*`
attrs point at those shared support artifacts, so they do not add another
per-package test-prebuild layer on top of that.

## Current strict reuse audit

The current success criterion is stricter than "Nix derivations are ordered
correctly": once a workspace crate has been compiled in the cacheable support
boundary, later package-specific test-prebuilds and nextest archive derivations
must reuse that artifact. Recompiling a generated feature-anchor crate is
allowed; recompiling `linnet`, `spenso`, `idenso`, `spenso-hep-lib`,
`gammalooprs`, `gammaloop-api`, or similar real workspace crates downstream is a
failure.

Two details were needed after the package-artifact fix:

- `clinnet` is its own nextest group. It has no `symbolica` dependency, but the
  old grouped `linnet` archive selected `clinnet` together with
  `linnet`/`linnest`/`linnet-py`. That broader group changed Cargo fingerprints
  through shared dependencies such as `indicatif`, so the final archive rebuilt
  `clinnet`.
- Test support builds binary targets only for the real component packages in
  that support SCC. Selecting binaries from synthetic consumers made support
  layers such as `crate-test-support-gammalooprs` build the
  `gammaloop-api` binary, which moved high-level work into a lower-level cache
  boundary.

The aggregate local check now completes:

```text
nix build --impure .#checks.x86_64-linux.gammaloop-nextest-binaries --no-link -L --print-out-paths
  /nix/store/x7ymfpn5j4q2nhxjhfk2c911bfikp7hj-gammaloop-nextest-binaries
```

During that run, changed package-specific prebuilds compiled only generated
anchor crates. Examples:

```text
gammaloop-nextest-binaries-spenso-spenso-test-deps:
  Compiling gammaloop-ci-nextest-spenso-spenso-features
  no Compiling linnet
  no Compiling spenso

gammaloop-nextest-binaries-linnet-linnet-test-deps:
  Compiling gammaloop-ci-nextest-linnet-linnet-features
  no Compiling linnet

gammaloop-nextest-binaries-core-gammalooprs-test-deps:
  Compiling gammaloop-ci-nextest-core-gammalooprs-features
  no Compiling gammalooprs
  no Compiling gammaloop-api
```

The final archive derivations also compiled only their generated group anchors:

```text
gammaloop-nextest-binaries-clinnet:
  Compiling gammaloop-ci-nextest-clinnet-features

gammaloop-nextest-binaries-linnet:
  Compiling gammaloop-ci-nextest-linnet-features

gammaloop-nextest-binaries-spenso:
  Compiling gammaloop-ci-nextest-spenso-features

gammaloop-nextest-binaries-vakint:
  Compiling gammaloop-ci-nextest-vakint-features

gammaloop-nextest-binaries-core:
  Compiling gammaloop-ci-nextest-core-features
```

Normal integration and Python API remain distinct checks. Normal integration no
longer enables `python-api-tests` and reuses the ordinary Rust support
artifacts. The Python API check intentionally compiles
`gammaloop-integration-tests` once with `python-api-tests`, and the Python module
path intentionally compiles a separate PyO3 feature family for `gammaloop-api`
and its dependents. That is not reusable with the normal Rust feature family.

## Follow-up audit: crate-level check splits

The per-crate package and nextest cache graph should not be copied directly to
every Cargo mode. I tested crate-level `clippy`, `doc`, and doctest splits and
rejected them because the actual logs showed downstream rebuilds of already
compiled dependency work.

The `clippy` experiment used a dependency-mode layer first, then a final real
package lint layer. The dependency-mode nodes improved reuse, but linting the
real selected package still changed Cargo's unit fingerprints and rebuilt the
expensive dependency graph:

```text
crate-clippy-deps-spenso:
  no Checking symbolica
  no Checking linnet

crate-clippy-gammalooprs:
  Compiling symbolica
  Checking linnet
  Checking spenso
  Checking idenso
```

The `doc` experiment failed in the same way. A narrow package could sometimes
reuse the inputs, but another package in the same dependency family re-entered
the Symbolica graph:

```text
crate-doc-linnet:
  no Checking symbolica

crate-doc-spenso-macros:
  Checking numerica
  Checking graphica
  Checking symbolica
```

Even the single workspace doc artifact fed from the merged workspace check
artifact does doc-mode work that is not reusable from the normal check/test
artifacts:

```text
cargo doc --profile ci-optim --locked --workspace ... --no-deps
  Compiling rug
  Compiling numerica
  Checking graphica
  Checking symbolica
```

The doctest split looked viable for one component, but failed the broader
representative path. The `spenso` component reused cleanly:

```text
crate-doctest-spenso:
  cargo test --profile ci-optim --doc ... -p spenso -p spenso-macros ...
  Finished in 0.29s
  no Compiling symbolica
  no Compiling linnet
  no Compiling spenso
```

When the `gammalooprs` doctest path pulled in downstream components, the split
either produced illegal transitive `pkg/feature` flags or, after narrowing those
flags to direct packages, recompiled real workspace crates:

```text
crate-doctest-gammaloop-tracing-filter:
  Compiling linnet
  Compiling spenso
  Compiling gammaloop-tracing-filter

crate-doctest-spenso-hep-lib:
  Compiling linnet
  Compiling spenso
  Compiling idenso
  Compiling spenso-hep-lib
```

The workspace aggregate artifact approach also proved too expensive as a cache
object. In the `000c330d7caae29bf4d4fd7136e4929c5e6f1f19` NixCI run,
`workspaceCheckCargoArtifacts`, `workspaceClippyCargoArtifacts`, and
`workspaceDocCargoArtifacts` spent minutes on cache download/upload and
compression even when they were already cached. They are invalidated by ordinary
workspace crate edits, and no later high-value test step consumes their exact
Cargo mode output. The replacement is a terminal `gammaloop-check` job that runs
`cargo check --all-targets` and writes only a tiny success marker.

So the current rule is: only expose crate-level CI attrs for a Cargo mode after
the persisted logs prove that downstream nodes do not compile real workspace
crates or the Symbolica/numerica/graphica graph. For now, crate-level reuse is
kept for package and nextest/test-binary artifacts. Workspace `check`, `clippy`,
`doc`, and doctest are terminal checks rather than misleading per-crate DAGs or
reusable artifact producers. They consume the base dependency artifact to avoid
starting completely cold, but they do not publish merged workspace target trees
back to the cache.

## Remaining caveats

- Dependency-only derivations can still compile generated feature-anchor or
  consumer crates and root-specific third-party crates. The public
  `crate-deps-*` package outputs now preserve the real current package artifact,
  and `crate-deps-gammalooprs` does not compile `linnet`, `spenso`, `idenso`,
  `vakint`, `spenso-hep-lib`, or `gammaloop-tracing-filter`.
- The generated nextest feature anchor is a build-source-only workspace member.
  Because it is not in `Cargo.lock`, nextest archive/test-binary Cargo commands
  use `--offline` rather than `--locked`.
- Test harness binaries are different Cargo units from normal library/package
  builds. The shared `crate-test-support-*` derivations are therefore the cache
  boundary for test binaries. The nextest package/archive derivations should
  only add generated feature-anchor packages and package fresh binaries.
- The fine-grained test-support SCC graph avoids downstream workspace
  recompilation by compiling dummy next-consumer contexts in earlier support
  derivations and stripping those dummy consumer artifacts afterward. This keeps
  ordinary downstream source edits from invalidating the upstream support
  derivation, but it intentionally increases the amount of work and archive size
  in earlier cacheable layers.
- Python ABI builds intentionally use different `gammaloop-api` features
  (`python_abi`, `pyo3-extension-module`) and remain a separate artifact family.
  Normal Rust test archives should not depend on that family.
- `doNotLinkInheritedArtifacts = true` still deep-copies inherited artifacts in
  some derivations. That should not cause recompilation, but it can add wall
  time and larger outputs. If compile reuse is fixed but CI still spends a lot
  of time copying artifact trees, the next follow-up is to remove that override
  and use Crane's default symlink-heavy inheritance mode where safe.
