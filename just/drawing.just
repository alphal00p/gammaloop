# Sync drawing Typst templates and Nix-built WASM bundles for development.
sync-drawing-assets:
    #!/usr/bin/env bash
    set -euo pipefail

    root="{{ justfile_directory() }}"
    cd "$root"

    install_asset() {
      local src="$1"
      local dest="$2"
      rm -f "$dest"
      install -m 0644 "$src" "$dest"
    }

    install_assets() {
      local target="$1"
      shift
      local src
      for src in "$@"; do
        install_asset "$src" "$target/$(basename "$src")"
      done
    }

    copy_gammaloop_templates() {
      local target="$1"
      mkdir -p "$target"
      install_assets "$target" assets/embedded/drawing/templates/*.typ
    }

    clean_old_flat_bundle() {
      local target="$1"
      rm -f \
        "$target/curve.typ" \
        "$target/draw.typ" \
        "$target/graph.typ" \
        "$target/linnest.typ" \
        "$target/physics-edge-style.typ" \
        "$target/subgraph.typ" \
        "$target/kurvst.wasm" \
        "$target/linnest.wasm" \
        "$target/crates/linnest/typst/kurvst.wasm"
    }

    copy_package_sources() {
      local target="$1"
      mkdir -p \
        "$target/crates/linnest/typst/src" \
        "$target/crates/kurvst/typst/src"
      install_asset crates/linnest/typst/typst.toml "$target/crates/linnest/typst/typst.toml"
      cp -R crates/linnest/typst/src/. "$target/crates/linnest/typst/src/"
      install_asset crates/kurvst/typst/typst.toml "$target/crates/kurvst/typst/typst.toml"
      cp -R crates/kurvst/typst/src/. "$target/crates/kurvst/typst/src/"
    }

    copy_wasm_bundles() {
      local target="$1"
      mkdir -p \
        "$target/crates/linnest/typst" \
        "$target/crates/kurvst/typst"
      install_asset result/linnest.wasm "$target/crates/linnest/typst/linnest.wasm"
      install_asset result/kurvst.wasm "$target/crates/kurvst/typst/kurvst.wasm"
    }

    if [ -d gammaloop_state/drawings/templates ]; then
      clean_old_flat_bundle gammaloop_state/drawings/templates
      copy_gammaloop_templates gammaloop_state/drawings/templates
      install_asset assets/embedded/drawing/justfile gammaloop_state/justfile
    fi

    nix build .#linnest-wasm --print-build-logs

    install_asset result/kurvst.wasm crates/kurvst/typst/kurvst.wasm
    install_asset result/linnest.wasm crates/linnest/typst/linnest.wasm

    if [ -d gammaloop_state/drawings/templates ]; then
      copy_package_sources gammaloop_state/drawings/templates
      copy_wasm_bundles gammaloop_state/drawings/templates
    fi

    cmp result/kurvst.wasm crates/kurvst/typst/kurvst.wasm
    cmp result/linnest.wasm crates/linnest/typst/linnest.wasm

# Sync drawing Typst templates and Cargo-built WASM bundles for non-Nix development.
sync-drawing-assets-cargo:
    #!/usr/bin/env bash
    set -euo pipefail

    root="{{ justfile_directory() }}"
    cd "$root"

    wasm_target="wasm32-unknown-unknown"
    wasm_profile="release"
    wasm_dir="target/$wasm_target/$wasm_profile"

    install_asset() {
      local src="$1"
      local dest="$2"
      rm -f "$dest"
      install -m 0644 "$src" "$dest"
    }

    install_assets() {
      local target="$1"
      shift
      local src
      for src in "$@"; do
        install_asset "$src" "$target/$(basename "$src")"
      done
    }

    copy_gammaloop_templates() {
      local target="$1"
      mkdir -p "$target"
      install_assets "$target" assets/embedded/drawing/templates/*.typ
    }

    clean_old_flat_bundle() {
      local target="$1"
      rm -f \
        "$target/curve.typ" \
        "$target/draw.typ" \
        "$target/graph.typ" \
        "$target/linnest.typ" \
        "$target/physics-edge-style.typ" \
        "$target/subgraph.typ" \
        "$target/kurvst.wasm" \
        "$target/linnest.wasm" \
        "$target/crates/linnest/typst/kurvst.wasm"
    }

    copy_package_sources() {
      local target="$1"
      mkdir -p \
        "$target/crates/linnest/typst/src" \
        "$target/crates/kurvst/typst/src"
      install_asset crates/linnest/typst/typst.toml "$target/crates/linnest/typst/typst.toml"
      cp -R crates/linnest/typst/src/. "$target/crates/linnest/typst/src/"
      install_asset crates/kurvst/typst/typst.toml "$target/crates/kurvst/typst/typst.toml"
      cp -R crates/kurvst/typst/src/. "$target/crates/kurvst/typst/src/"
    }

    copy_wasm_bundles() {
      local target="$1"
      mkdir -p \
        "$target/crates/linnest/typst" \
        "$target/crates/kurvst/typst"
      install_asset "$wasm_dir/linnest.wasm" "$target/crates/linnest/typst/linnest.wasm"
      install_asset "$wasm_dir/kurvst.wasm" "$target/crates/kurvst/typst/kurvst.wasm"
    }

    if command -v rustup >/dev/null 2>&1; then
      rustup target add "$wasm_target"
    fi

    cargo build --release -p linnest -p kurvst --features linnest/custom --target "$wasm_target"

    install_asset "$wasm_dir/kurvst.wasm" crates/kurvst/typst/kurvst.wasm
    install_asset "$wasm_dir/linnest.wasm" crates/linnest/typst/linnest.wasm

    if [ -d gammaloop_state/drawings/templates ]; then
      clean_old_flat_bundle gammaloop_state/drawings/templates
      copy_gammaloop_templates gammaloop_state/drawings/templates
      install_asset assets/embedded/drawing/justfile gammaloop_state/justfile
      copy_package_sources gammaloop_state/drawings/templates
      copy_wasm_bundles gammaloop_state/drawings/templates
    fi

    cmp "$wasm_dir/kurvst.wasm" crates/kurvst/typst/kurvst.wasm
    cmp "$wasm_dir/linnest.wasm" crates/linnest/typst/linnest.wasm
