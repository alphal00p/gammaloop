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

The repo already used Hakari and Guppy to create a cacheable workspace graph, and
the offline Cargo feature graph shows the workspace hack pulling Symbolica with
the same main third-party feature set (`bincode`, `gmp`, `serde`) through the
Rust workspace.

The local flaw was that the Nix derivations did not always ask Cargo for the
same feature set as the artifact they inherited:

- `cargoArtifacts` was built with the full workspace feature superset.
- `gammaloop-clippy`, docs, doctests, and some product/package derivations used
  narrower package-local feature arguments.
- Per-crate package builds inherited dependency package artifacts, but their
  consuming Cargo command only enabled features for the selected package. For
  example, `crate-spenso-hep-lib` could inherit `crate-idenso`, but then ask
  Cargo for an `idenso` variant without the same extra features used by
  `crate-idenso`. Cargo correctly compiled a second `idenso` unit.

This matches the Crane FAQ warning: differing package selections or feature
flags between `buildDepsOnly` and consumers invalidate artifact reuse.

## Fix applied

The flake now constructs Cargo args for the relevant workspace dependency
closure, not just the selected package:

- `cargoPackageCiArgsFor package` enables CI features for `package` and its
  normal workspace dependency closure.
- `cargoPackageTestArgsFor package` does the same for test dependency closure.
- `ciArgs.cargoExtraArgs` now includes the full workspace feature superset, so
  clippy/doc/doctest consume the same kind of dependency artifact as
  `cargoArtifacts`.
- nextest archive derivations enable features for their archive input closure,
  so the cacheable archive build can reuse the normal package artifacts it
  receives.

Representative inspected commands after the fix:

```text
cargoWithProfile build --locked -p spenso-hep-lib \
  --features idenso/bincode,idenso/reference-cases,linnet/bincode,linnet/serde,linnet/symbolica,spenso-macros/shadowing,spenso/shadowing
```

```text
cargoWithProfile clippy --locked --workspace \
  --features gammaloop-integration-tests/python-api-tests,gammaloop-tracing-filter/clap,gammaloop-tracing-filter/symbolica,idenso/bincode,idenso/reference-cases,linnet/bincode,linnet/serde,linnet/symbolica,spenso-macros/shadowing,spenso/shadowing \
  --all-targets -- --deny warnings
```

```text
cargo nextest archive --cargo-profile ci-optim \
  --locked -p idenso -p spenso -p spenso-hep-lib -p spenso-macros -p spynso3 \
  --features idenso/bincode,idenso/reference-cases,linnet/bincode,linnet/serde,linnet/symbolica,spenso-macros/shadowing,spenso/shadowing
```

## Verified per-crate feature sets

"Common feature set" here means the common qualified feature set for a crate's
workspace dependency closure, not one global feature set forced onto every
crate. `clinnet` and `vakint` intentionally stay empty because they do not have
feature-bearing Symbolica workspace dependencies.

The table below was extracted from the generated derivations with
`nix derivation show`. For every workspace package, the extracted feature string
for `packages.x86_64-linux.crate-deps-<package>` was compared with
`packages.x86_64-linux.crate-<package>`.

```text
matched feature strings for 16 workspace packages
```

| Workspace crate | Qualified features used by both `crate-deps-*` and `crate-*` |
| --- | --- |
| `clinnet` | none |
| `gammaloop-api` | `gammaloop-tracing-filter/clap,gammaloop-tracing-filter/symbolica,idenso/bincode,idenso/reference-cases,linnet/bincode,linnet/serde,linnet/symbolica,spenso-macros/shadowing,spenso/shadowing` |
| `gammaloop-integration-tests` | `gammaloop-tracing-filter/clap,gammaloop-tracing-filter/symbolica,idenso/bincode,idenso/reference-cases,linnet/bincode,linnet/serde,linnet/symbolica,spenso-macros/shadowing,spenso/shadowing` |
| `gammaloop-tracing-filter` | `gammaloop-tracing-filter/clap,gammaloop-tracing-filter/symbolica,linnet/bincode,linnet/serde,linnet/symbolica,spenso-macros/shadowing,spenso/shadowing` |
| `gammaloop-tracing-filter-macros` | none |
| `gammaloop-workspace-hack` | none |
| `gammalooprs` | `gammaloop-tracing-filter/clap,gammaloop-tracing-filter/symbolica,idenso/bincode,idenso/reference-cases,linnet/bincode,linnet/serde,linnet/symbolica,spenso-macros/shadowing,spenso/shadowing` |
| `idenso` | `idenso/bincode,idenso/reference-cases,linnet/bincode,linnet/serde,linnet/symbolica,spenso-macros/shadowing,spenso/shadowing` |
| `linnest` | `linnet/bincode,linnet/serde,linnet/symbolica` |
| `linnet` | `linnet/bincode,linnet/serde,linnet/symbolica` |
| `linnet-py` | `linnet/bincode,linnet/serde,linnet/symbolica` |
| `spenso` | `linnet/bincode,linnet/serde,linnet/symbolica,spenso-macros/shadowing,spenso/shadowing` |
| `spenso-hep-lib` | `idenso/bincode,idenso/reference-cases,linnet/bincode,linnet/serde,linnet/symbolica,spenso-macros/shadowing,spenso/shadowing` |
| `spenso-macros` | `spenso-macros/shadowing` |
| `spynso3` | `idenso/bincode,idenso/reference-cases,linnet/bincode,linnet/serde,linnet/symbolica,spenso-macros/shadowing,spenso/shadowing` |
| `vakint` | none |

## Remaining caveats

- Test harness binaries are different Cargo units from normal library/package
  builds. The nextest archive job is therefore the cache boundary for test
  binaries; the in-repo NixCI test runner should only execute the archive.
- Python ABI builds intentionally use different `gammaloop-api` features
  (`python_abi`, `pyo3-extension-module`) and remain a separate artifact family.
- `doNotLinkInheritedArtifacts = true` still deep-copies inherited artifacts in
  some derivations. That should not cause recompilation, but it can add wall
  time and larger outputs. If compile reuse is fixed but CI still spends a lot
  of time copying artifact trees, the next follow-up is to remove that override
  and use Crane's default symlink-heavy inheritance mode where safe.
