= Rebuilding the bundled CeTZ core

The bundled binary is built from CeTZ `v0.5.1`, commit
`100ea3c52342c8bf866e8f154589ce8cd309265c`, with
`batched-bounds.patch`, `content-placement.patch`, `sampled-curves.patch`,
`marker-geometry.patch`, and `marker-footprints.patch` applied in that order. The source remains under CeTZ's
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

The sampled-curve patch moves CeTZ's existing fixed-sample cubic arc-length
and distance-to-parameter loops into the native core. It retains the same
sample counts, Bernstein evaluation order, vector dimensions, endpoint behavior
and integer/float return types. It does not use adaptive curve measurements.

The marker patch batches resolved marker contacts, alignment matrices, and
shaft trimming. It reuses the fixed-sample cubic owner and preserves the
reference transform and subcurve arithmetic. Mark shapes and styling remain
in Typst. Invalid dimensions and nonfinite geometry are rejected before
indexing or painting; empty unmarked carriers remain valid. The patch also
adds the pinned `libm` dependency used by Typst's trigonometric functions.

The footprint patch batches candidate arrow collision bounds from shared numeric
marker templates. It resolves each carrier's relative stations, normalizes and
joins its paths, and uses the same marker contact/trim and cubic-extrema owners
as painting. Shaft bounds remain per command; each marker drawable retains one
whole-path bound. Content keeps its original dimensions after position
transforms. Integer and floating-point CBOR coordinates are both accepted;
invalid dimensions, nonfinite values, and overflowing geometry are rejected.
Style resolution and context-dependent custom shape evaluation remain in Typst.
Canonical numeric carrier packets may remain encoded between plugins. The same
footprint input owner expands them and selects the resolved single- or multi-segment
shaft style before entering the shared painter; bounds can also remain encoded.

The checked-in binary uses Rust 1.98.1 (`48a229cea`, LLVM 21.1.8), the patched
`cetz-core/Cargo.lock`, and the upstream release profile with
`CARGO_PROFILE_RELEASE_OPT_LEVEL=3` for speed. The build uses the cached
`rust-lld` recorded below: LLD 22.1.8 (LLVM revision
`52ed14fcd56afc30f9cccd8ca8ce237c2eef7e04`), from a separate Rust toolchain. Run from
the repository root with Rust 1.98.1 and its `wasm32-unknown-unknown` target
installed. The build runs in the temporary checkout, outside this repository's
Cargo configuration, with no custom Rust flags or target directory:

```sh
cetz_source=$(mktemp -d)
git clone --depth 1 --branch v0.5.1 https://github.com/cetz-package/cetz.git "$cetz_source"
test "$(git -C "$cetz_source" rev-parse HEAD)" = 100ea3c52342c8bf866e8f154589ce8cd309265c
git -C "$cetz_source" apply "$PWD/crates/linnet-py/vendor/typst-packages/preview/cetz/0.5.1/cetz-core/batched-bounds.patch"
git -C "$cetz_source" apply "$PWD/crates/linnet-py/vendor/typst-packages/preview/cetz/0.5.1/cetz-core/content-placement.patch"
git -C "$cetz_source" apply "$PWD/crates/linnet-py/vendor/typst-packages/preview/cetz/0.5.1/cetz-core/sampled-curves.patch"
git -C "$cetz_source" apply "$PWD/crates/linnet-py/vendor/typst-packages/preview/cetz/0.5.1/cetz-core/marker-geometry.patch"
git -C "$cetz_source" apply "$PWD/crates/linnet-py/vendor/typst-packages/preview/cetz/0.5.1/cetz-core/marker-footprints.patch"
(
  cd "$cetz_source" || exit
  unset RUSTC RUSTDOC RUSTFLAGS CARGO_ENCODED_RUSTFLAGS CARGO_BUILD_RUSTFLAGS CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_RUSTFLAGS CARGO_TARGET_DIR
  export CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_LINKER=/nix/store/c3bpv60f5wgp4pmgsifaykxk214h6463-rust-stable-with-components-2026-09-03/lib/rustlib/x86_64-unknown-linux-gnu/bin/rust-lld
  CARGO_PROFILE_RELEASE_OPT_LEVEL=3 cargo +1.98.1 build --locked --release --target wasm32-unknown-unknown --manifest-path cetz-core/Cargo.toml
)
cp "$cetz_source/cetz-core/target/wasm32-unknown-unknown/release/cetz_core.wasm" crates/linnet-py/vendor/typst-packages/preview/cetz/0.5.1/cetz-core/cetz_core.wasm
```

The resulting binary's SHA-256 is
`a84fde7399862dd8af3c503721532fb46107c2d5edf0daa9a5d55eeed570aaec`.
The bounds behavior fixture compares batched results against the scalar CeTZ
exports, including degenerate curves, mixed dimensions, signed zeros and
nonfinite aggregation. It runs as part of the Clinnet public Typst fixture. The identity-target fixture
also compares batched content against scalar placement across styles, transforms,
coordinate hooks, degenerate frames and overflowing coordinates.

The sampled-curve fixture compares exact CBOR results against the scalar loops
across sample counts, signed zeros, mixed dimensions and forward/reverse distances.

The marker fixture compares native placement with the original Typst contact
and trimming algorithm, including short carriers, reversed/slanted markers,
zero-length marks, and repeated controls. Native tests also cover malformed
plans and empty unmarked carriers. Marker numeric output uses floating-point
coordinates; mathematical integer values need not retain CBOR integer tags.

Candidate-footprint tests compare exact CBOR bounds with the ordinary painted
owner across relative and repeated marks, curved and short carriers, custom
content, hidden drawables, transforms, and per-path stroke styles. Native tests
also cover trimming, relative stations, singular transforms, and invalid inputs.
Packet tests require byte-identical bounds to expanded numeric input, preserve
candidate order across packets, and reject invalid transport and geometry.
