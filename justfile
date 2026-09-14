# Gammaloop build and development commands

ci_cargo_profile := "ci-optim"

mod linnet 'crates/linnet/Justfile'

# Imported recipes keep the same top-level commands and repository working directory.
import 'just/ci.just'
import 'just/tests.just'
import 'just/drawing.just'

alias sync-draw-assets := sync-drawing-assets

drawing_asset_sync := "none"

_sync-drawing-assets:
    #!/usr/bin/env bash
    set -euo pipefail
    case "{{ drawing_asset_sync }}" in
      nix)
        just sync-drawing-assets
        ;;
      cargo)
        just sync-drawing-assets-cargo
        ;;
      none|off|false)
        ;;
      *)
        echo "unknown drawing_asset_sync={{ drawing_asset_sync }}; expected nix, cargo, or none" >&2
        exit 2
        ;;
    esac

# Build gammaloop Python CLI with UFO support and dev-optim profile
build-cli: _sync-drawing-assets
    cargo build -p gammaloop-api --bin gammaloop --features ufo_support --profile dev-optim

# Build gammaloop CLI without Python (no pyo3)
build-cli-no-pyo3: _sync-drawing-assets
    cargo build -p gammaloop-api --bin gammaloop --no-default-features --features cli,no_pyo3 --profile dev-optim

# Build gammaloop Python CLI with UFO support and stable ABI (dev-optim profile)
build-cli-abi: _sync-drawing-assets
    cargo build -p gammaloop-api --bin gammaloop --features ufo_support,python_abi --profile dev-optim

# Build gammaloop Python CLI in release mode
build-cli-release: _sync-drawing-assets
    cargo build -p gammaloop-api --bin gammaloop --features ufo_support --release

# Build gammaloop Python CLI in release mode with stable ABI
build-cli-release-abi: _sync-drawing-assets
    cargo build -p gammaloop-api --bin gammaloop --features ufo_support,python_abi --release

# Build gammaloop Python API with UFO support and dev-optim profile
build-api: _sync-drawing-assets
    maturin develop -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_api --profile=dev-optim

# Build gammaloop Python API with UFO support and stable ABI (dev-optim profile)
build-api-abi: _sync-drawing-assets
    maturin develop -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_abi --profile=dev-optim

# Build gammaloop Python API with UFO support and dev-optim profile
build-api-wheel: _sync-drawing-assets
    maturin build -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_api --profile=dev-optim

# Build gammaloop Python API with UFO support and stable ABI (dev-optim profile)
build-api-abi-weel: _sync-drawing-assets
    maturin build -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_abi --profile=dev-optim

# Build gammaloop Python API with UFO support and release profile
build-api-release: _sync-drawing-assets
    maturin develop -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_api --profile=release

# Build gammaloop Python API with UFO support and stable ABI (release profile)
build-api-abi-release: _sync-drawing-assets
    maturin develop -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_abi --profile=release

# Build gammaloop Python API wheel with UFO support and release profile
build-api-release-wheel: _sync-drawing-assets
    maturin build -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_api --profile=release

# Build gammaloop Python API weel with UFO support and stable ABI (release profile)
build-api-abi-release-wheel: _sync-drawing-assets
    maturin build -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_abi --profile=release

# Build all workspace packages
build-all:
    cargo build --workspace

# Clean build artifacts
clean:
    cargo clean

# Check code without building
check:
    cargo check --workspace --all-targets --locked

doc:
    cargo doc --workspace --no-deps --locked --profile {{ ci_cargo_profile }}

# Format code
fmt *lint_args:
    cargo fmt --all {{ lint_args }}

# Run clippy linter
clippy *lint_args:
    cargo clippy --workspace --all-targets --locked --profile {{ ci_cargo_profile }} {{ lint_args }}

# Build everything in release mode for maximum performance
build-release-all:
    cargo build --workspace --release

# Quick development cycle: build deps, then build in release
dev-release: build-deps-nix build-release-all

# Run gammaloop
run *ARGS:
    cargo run -p gammaloop-api --bin gammaloop --features ufo_support --profile dev-optim -- {{ ARGS }}
