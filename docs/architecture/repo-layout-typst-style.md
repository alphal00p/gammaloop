# Proposed Typst-Style Repository Layout

This document proposes an incremental migration toward a `typst/typst`-style workspace layout with a `crates/` directory and a cleaner top-level.

## Target Top-Level Shape

```text
.
├── crates/
│   ├── gammalooprs/          # current root crate (`src/`)
│   ├── gammaloop-api/        # current `gammaloop-api/`
│   └── integration-tests/    # current `integration-tests/`
├── docs/
├── tests/                    # repo-level tests + fixtures
├── assets/
├── models/
├── examples/
├── tools/                    # helper scripts (new)
├── Cargo.toml                # workspace root only
├── Cargo.lock
└── justfile
```

## Incremental Migration Plan

1. Move all Rust crates under `crates/` while keeping package names unchanged.
2. Turn root `Cargo.toml` into a workspace manifest and keep shared lint/profile config there.
3. Keep integration fixtures in `tests/` and crate-specific fixtures under each crate.
4. Move ad-hoc scripts into `tools/` and reserve root for project entry points.
5. Keep local run outputs under `.local/scratch/` and avoid creating new top-level working directories.

## Module Foldering Strategy Inside `gammalooprs`

To reduce a flat `src/`, move high-churn top-level modules into directory modules first, then split internals:

1. Already moved in this change:
- `src/momentum.rs` -> `src/momentum/mod.rs`
- `src/model.rs` -> `src/model/mod.rs`
- `src/integrate.rs` -> `src/integrate/mod.rs`
- `src/observables.rs` -> `src/observables/mod.rs`

2. Next low-risk foldering candidates:
- `src/signature.rs` -> `src/signature/mod.rs`
- `src/momentum_sample.rs` -> `src/momentum_sample/mod.rs`
- `src/improve_ps.rs` -> `src/improve_ps/mod.rs`
- `src/initialisation.rs` -> `src/initialisation/mod.rs`
- `src/inspect.rs` -> `src/inspect/mod.rs`

3. Then split internals by concern while preserving public paths with re-exports:
- `momentum/{energy,rotation,polarization,spinor}.rs`
- `model/{particle,parameter,coupling,lorentz,ufo}.rs`
- `integrate/{engine,state,stats,format}.rs`

## Compatibility Rules During Migration

1. Keep existing `crate::module` paths stable while splitting.
2. Use `pub use` in `mod.rs` to preserve public API.
3. Migrate one module at a time and run `cargo check` between steps.
