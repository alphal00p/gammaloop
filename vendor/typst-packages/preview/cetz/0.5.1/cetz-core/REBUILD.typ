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
`CARGO_PROFILE_RELEASE_OPT_LEVEL=3` for speed. It links with the `rust-lld` of a
separate Rust toolchain, the one in this repository's pinned development shell:
LLD 22.1.8 (LLVM revision `52ed14fcd56afc30f9cccd8ca8ce237c2eef7e04`).

Rustc stores dependency source paths for panic locations, so the bytes depend on the
`CARGO_HOME` path, and nothing else about the build machine enters them. The checked-in
binary was built with `CARGO_HOME=/build/rust/.cargo`. It ran in a sandbox whose root holds only
the Nix store, the checkout at `/build/source` and that Cargo home, without network after
`cargo fetch --locked`. Another Cargo home gives the same code with different path strings.
The compiler must be a Rust 1.98.1 build linked against LLVM 21.1.8, such as nixpkgs' `rustc`
1.98.1. The official 1.98.1 distribution, which bundles LLVM 22.1.8, produces different code.

To reproduce the checked-in bytes, run the recipe below from the repository root. It needs
network access for the clone and `cargo fetch`, and `bwrap` (bubblewrap) on `PATH`. The whole
recipe is one script run inside the pinned development shell, which supplies the linker and the
host C toolchain used by build scripts. The compiler and Cargo come from the nixpkgs revision
pinned in `flake.lock`; the script checks both compiler versions and the linker version. Any
failed step stops the script, and its last command compares the result with the recorded digest.
It leaves `vendor/` unchanged:

```sh
nix develop --command bash -euo pipefail -s <<'EOF'
expected=d1250cd165874a289b28bc410f27ea449182bfd5daf5ef3cbdc11a8bbbcebd85
patches="$PWD/vendor/typst-packages/preview/cetz/0.5.1/cetz-core"
rust_lld="$(rustc --print sysroot)/lib/rustlib/x86_64-unknown-linux-gnu/bin/rust-lld"
"$rust_lld" -flavor wasm --version | grep -F 'LLD 22.1.8' | grep -F 52ed14fcd56afc30f9cccd8ca8ce237c2eef7e04
rustc_bin="$(nix build --no-link --print-out-paths --inputs-from . 'nixpkgs#rustc^out')/bin/rustc"
cargo_bin="$(nix build --no-link --print-out-paths --inputs-from . 'nixpkgs#cargo^out')/bin/cargo"
"$rustc_bin" -vV | grep -Fx 'release: 1.98.1'
"$rustc_bin" -vV | grep -Fx 'LLVM version: 21.1.8'
root=$(mktemp -d)
git clone --depth 1 --branch v0.5.1 https://github.com/cetz-package/cetz.git "$root/source"
test "$(git -C "$root/source" rev-parse HEAD)" = 100ea3c52342c8bf866e8f154589ce8cd309265c
for patch in batched-bounds content-placement sampled-curves marker-geometry marker-footprints; do
  git -C "$root/source" apply "$patches/$patch.patch"
done
mkdir -p "$root/rust/.cargo"
CARGO_HOME="$root/rust/.cargo" "$cargo_bin" fetch --locked --manifest-path "$root/source/cetz-core/Cargo.toml"
bwrap --unshare-net --ro-bind /nix /nix --ro-bind /usr /usr --ro-bind /bin /bin --ro-bind /etc /etc \
  --proc /proc --dev /dev --bind "$root" /build --chdir /build/source --clearenv \
  --setenv PATH "$PATH" --setenv HOME /build --setenv CARGO_HOME /build/rust/.cargo \
  --setenv RUSTC "$rustc_bin" --setenv CARGO_NET_OFFLINE true \
  --setenv CARGO_PROFILE_RELEASE_OPT_LEVEL 3 \
  --setenv CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_LINKER "$rust_lld" \
  "$cargo_bin" build --locked --release --target wasm32-unknown-unknown --manifest-path cetz-core/Cargo.toml
actual=$(sha256sum "$root/source/cetz-core/target/wasm32-unknown-unknown/release/cetz_core.wasm" | cut -d ' ' -f 1)
echo "built $actual, recorded $expected; output in $root"
test "$actual" = "$expected"
EOF
```

The resulting binary's SHA-256 is
`d1250cd165874a289b28bc410f27ea449182bfd5daf5ef3cbdc11a8bbbcebd85`.

For iterating on the patches, an optional rustup build is quicker. It needs Rust 1.98.1 with its
`wasm32-unknown-unknown` target installed through rustup, run outside the development shell,
and it still takes `rust-lld` from the pinned shell. Its official compiler produces a
working binary whose bytes differ from the checked-in one, so the recipe writes to a scratch
directory and never replaces `cetz_core.wasm`. The build runs in the temporary checkout,
outside this repository's Cargo configuration, with no custom Rust flags or target directory.
Run it from the repository root; any failed step stops the subshell:

```sh
(
  set -euo pipefail
  rust_lld="$(nix develop --command rustc --print sysroot)/lib/rustlib/x86_64-unknown-linux-gnu/bin/rust-lld"
  "$rust_lld" -flavor wasm --version | grep -F 'LLD 22.1.8' | grep -F 52ed14fcd56afc30f9cccd8ca8ce237c2eef7e04
  patches="$PWD/vendor/typst-packages/preview/cetz/0.5.1/cetz-core"
  cetz_source=$(mktemp -d)
  out=$(mktemp -d)
  git clone --depth 1 --branch v0.5.1 https://github.com/cetz-package/cetz.git "$cetz_source"
  test "$(git -C "$cetz_source" rev-parse HEAD)" = 100ea3c52342c8bf866e8f154589ce8cd309265c
  for patch in batched-bounds content-placement sampled-curves marker-geometry marker-footprints; do
    git -C "$cetz_source" apply "$patches/$patch.patch"
  done
  cd "$cetz_source"
  unset RUSTC RUSTDOC RUSTFLAGS CARGO_ENCODED_RUSTFLAGS CARGO_BUILD_RUSTFLAGS CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_RUSTFLAGS CARGO_TARGET_DIR
  export CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_LINKER="$rust_lld"
  CARGO_PROFILE_RELEASE_OPT_LEVEL=3 cargo +1.98.1 build --locked --release --target wasm32-unknown-unknown --manifest-path cetz-core/Cargo.toml
  cp cetz-core/target/wasm32-unknown-unknown/release/cetz_core.wasm "$out/cetz_core.wasm"
  sha256sum "$out/cetz_core.wasm"
)
```

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
