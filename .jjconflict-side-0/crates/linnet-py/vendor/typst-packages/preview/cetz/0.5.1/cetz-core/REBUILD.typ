= Rebuilding the bundled CeTZ core

The bundled binary is built from CeTZ `v0.5.1`, commit
`100ea3c52342c8bf866e8f154589ce8cd309265c`, with
`batched-bounds.patch` and `content-placement.patch` applied. The source remains under CeTZ's
LGPL-3.0-or-later license; see `../LICENSE`.

The bounds patch adds one export that calls the existing cubic-extrema and AABB
routines. Their coefficient rounding, arithmetic and traversal order are
unchanged. The original exports remain available. The extra `finite` result
allows callers to aggregate a finite path's two bounds corners while retaining
the full point stream for nonfinite coordinates.

The content patch adds a batch export for prepared rectangular content. It
preserves the scalar transform, corner and frame arithmetic and line order.
Nonfinite geometry is flagged so the caller retains Typst's scalar comparison
behavior. Styles, measurements, wrapping and drawable templates remain in Typst.

The checked-in binary uses Rust 1.98.0 (`88d9e12ae`, LLVM 22.1.8), the upstream
`cetz-core/Cargo.lock`, and the upstream release profile with `CARGO_PROFILE_RELEASE_OPT_LEVEL=3`
for speed. Run from the repository
root with that toolchain and its `wasm32-unknown-unknown` target installed.
The build runs from the temporary checkout, outside this repository's Cargo
configuration, with no custom Rust flags or target directory:

```sh
cetz_source=$(mktemp -d)
git clone --depth 1 --branch v0.5.1 https://github.com/cetz-package/cetz.git "$cetz_source"
test "$(git -C "$cetz_source" rev-parse HEAD)" = 100ea3c52342c8bf866e8f154589ce8cd309265c
git -C "$cetz_source" apply "$PWD/crates/linnet-py/vendor/typst-packages/preview/cetz/0.5.1/cetz-core/batched-bounds.patch"
git -C "$cetz_source" apply "$PWD/crates/linnet-py/vendor/typst-packages/preview/cetz/0.5.1/cetz-core/content-placement.patch"
(
  cd "$cetz_source" || exit
  unset RUSTFLAGS CARGO_ENCODED_RUSTFLAGS CARGO_BUILD_RUSTFLAGS CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_RUSTFLAGS CARGO_TARGET_DIR
  CARGO_PROFILE_RELEASE_OPT_LEVEL=3 cargo +1.98.0 build --locked --release --target wasm32-unknown-unknown --manifest-path cetz-core/Cargo.toml
)
cp "$cetz_source/cetz-core/target/wasm32-unknown-unknown/release/cetz_core.wasm" crates/linnet-py/vendor/typst-packages/preview/cetz/0.5.1/cetz-core/cetz_core.wasm
```

The resulting binary's SHA-256 is
`fd5197c55f3160d260473c581f75eff7b1e0fb7f40559ce143c1a5661501a6e2`.
The bounds behavior fixture compares batched results against the scalar CeTZ
exports, including degenerate curves, mixed dimensions, signed zeros and
nonfinite aggregation. It runs as part of the Clinnet public Typst fixture. The identity-target fixture
also compares batched content against scalar placement across styles, transforms,
coordinate hooks, degenerate frames and overflowing coordinates.
