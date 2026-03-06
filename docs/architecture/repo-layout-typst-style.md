# Typst-Style Repository Layout

This document tracks the `typst/typst`-style workspace direction with a `crates/` directory and a cleaner top-level.

## Target Top-Level Shape

```text
.
├── crates/
│   ├── gammalooprs/          # former top-level crate (`src/`)
│   ├── gammaloop-api/
├── tests/                    # integration-test crate + fixtures
├── docs/
├── assets/                   # includes `assets/models/`
├── examples/
├── tools/                    # helper scripts (new)
├── Cargo.toml                # workspace root only
├── Cargo.lock
└── justfile
```

## Migration Status

1. Rust crates are under `crates/`, with integration tests promoted to top-level `tests/`.
2. Root `Cargo.toml` is now a workspace manifest.
3. Integration fixtures live under `tests/resources/`.
4. Local run outputs should stay outside the top-level (for example `.local/scratch/`).

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
