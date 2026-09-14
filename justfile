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

# Regenerate every SVG used by the website from its committed Typst source.
docs-svg-assets:
    scripts/render-docs-svg-assets.sh

# Fail when any checked-in website SVG is stale relative to its Typst source,
# Linnest renderer, palette, or real process/test graph data.
docs-svg-assets-check:
    #!/usr/bin/env bash
    set -euo pipefail

    temp_base=${TMPDIR:-/tmp}
    check_root=$(mktemp -d "$temp_base/alphal00p-docs-svg.XXXXXX")
    cleanup() {
        case "$check_root" in
            "$temp_base"/alphal00p-docs-svg.*) rm -rf -- "$check_root" ;;
        esac
    }
    trap cleanup EXIT

    scripts/render-docs-svg-assets.sh "$check_root"
    checked_assets=(
        docs/assets/graphs/portal-*.svg
        docs/assets/local-unitarity-*.svg
        docs/assets/spensologo.svg
        assets/gammalooplogo*.svg
    )
    generated_assets=(
        "$check_root"/docs/assets/graphs/portal-*.svg
        "$check_root"/docs/assets/local-unitarity-*.svg
        "$check_root"/docs/assets/spensologo.svg
        "$check_root"/assets/gammalooplogo*.svg
    )
    [ "${#checked_assets[@]}" -eq 28 ]
    [ "${#generated_assets[@]}" -eq 28 ]
    for checked in "${checked_assets[@]}"; do
        cmp "$checked" "$check_root/$checked"
    done

# Build one product documentation site, or all five sites.
docs-site PRODUCT="all" CHANNEL="latest" SNAPSHOT_TAG="" OUTPUT="target/alphal00p-docs":
    #!/usr/bin/env bash
    set -euo pipefail

    channel={{ quote(CHANNEL) }}
    snapshot_tag={{ quote(SNAPSHOT_TAG) }}
    args=(
        build
        --product {{ quote(PRODUCT) }}
        --channel "$channel"
        --output {{ quote(OUTPUT) }}
    )

    case "$channel" in
        latest)
            if [ -n "$snapshot_tag" ]; then
                echo "SNAPSHOT_TAG is only valid when CHANNEL=snapshot" >&2
                exit 2
            fi
            ;;
        snapshot)
            if [ -z "$snapshot_tag" ]; then
                echo "SNAPSHOT_TAG is required when CHANNEL=snapshot" >&2
                exit 2
            fi
            args+=(--snapshot-tag "$snapshot_tag")
            ;;
        *)
            echo "CHANNEL must be latest or snapshot" >&2
            exit 2
            ;;
    esac

    unset GITHUB_REF_NAME
    just docs-check
    cargo run --locked -p alphal00p-docs-builder -- "${args[@]}"

# Validate the five-product documentation registry and generated inputs.
docs-check:
    cargo run --locked -p alphal00p-docs-catalogs --features gammaloop-reference --bin alphal00p-docs-gammaloop-reference -- --check
    cargo run --locked -p alphal00p-docs-catalogs --features vakint-reference --bin alphal00p-docs-vakint-reference -- --check
    cargo run --locked -p alphal00p-docs-python-exporter --features gammaloop -- gammaloop-python docs/api/python/gammaloop-python.pyi --check
    cargo run --locked -p alphal00p-docs-python-exporter --features linnet -- linnet-py docs/api/python/linnet-py.pyi --check
    cargo run --locked -p alphal00p-docs-python-exporter --features spenso -- spynso3 docs/api/python/spynso3.pyi --check
    cargo run --locked -p alphal00p-docs-python-exporter --features idenso -- idenso-community docs/api/python/idenso-community.pyi --check
    cargo run --locked -p alphal00p-docs-python-exporter --features vakint -- vakint-community docs/api/python/vakint-community.pyi --check
    cargo test --locked -p alphal00p-docs-python-exporter
    cargo test --locked -p alphal00p-docs-python-exporter --features gammaloop gammaloop_runtime_surface_and_signatures_match_the_docs_stub
    cargo test --locked -p alphal00p-docs-examples
    cargo run --locked -p alphal00p-docs-builder -- check
    just docs-svg-assets-check
    just docs-linnet-python-check

# Build Linnet's extension and compare its real import surface with the checked-in stub.
docs-linnet-python-check:
    #!/usr/bin/env bash
    set -euo pipefail

    temp_base=${TMPDIR:-/tmp}
    test_root=$(mktemp -d "$temp_base/alphal00p-docs-linnet-python.XXXXXX")
    cleanup() {
        case "$test_root" in
            "$temp_base"/alphal00p-docs-linnet-python.*)
                rm -rf -- "$test_root"
                ;;
        esac
    }
    trap cleanup EXIT

    python_bin=${PYTHON_BIN_PATH:-python3}
    uv venv "$test_root/venv" --python "$python_bin"
    VIRTUAL_ENV="$test_root/venv" maturin develop --uv --locked --manifest-path crates/linnet-py/Cargo.toml --features extension-module,abi3-py310
    "$test_root/venv/bin/python" -m unittest crates/linnet-py/tests/test_basic.py

# Regenerate source-backed CLI/settings and topology/dependency snapshots.
docs-generated:
    cargo run --locked -p alphal00p-docs-catalogs --features gammaloop-reference --bin alphal00p-docs-gammaloop-reference
    cargo run --locked -p alphal00p-docs-catalogs --features vakint-reference --bin alphal00p-docs-vakint-reference

# Run the provisioned FORM-backed Vakint documentation examples without pySecDec.
docs-vakint-form-check:
    VAKINT_SKIP_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test tensor_reduction_tests
    VAKINT_SKIP_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test integral_evaluation_freeform_tests test_integrate_1l_decorated_indices_matad -- --exact
    VAKINT_SKIP_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test integral_evaluation_freeform_tests test_integrate_1l_decorated_indices_fmft -- --exact
    VAKINT_SKIP_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test integral_alphaloop_vs_matad_tests

# Run one numerical pySecDec smoke comparison. The caller must provision pySecDec.
docs-vakint-pysecdec-smoke-check:
    RUN_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test integral_evaluation_pysecdec_tests test_integrate_1l_simple -- --exact

# Run the complete opt-in pySecDec and analytic-backend comparison tier.
docs-vakint-pysecdec-full-check:
    RUN_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test integral_evaluation_pysecdec_tests
    RUN_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test integral_comparison_vs_pysecdec_tests

# Continuously build and serve one product with native browser reloads.
docs-watch PRODUCT PORT="8000" OUTPUT="target/alphal00p-docs-watch":
    #!/usr/bin/env bash
    set -euo pipefail

    exec cargo run --locked --profile docs-watch -p alphal00p-docs-builder \
        --features persistent-typst -- \
        watch \
        --product {{ quote(PRODUCT) }} \
        --port {{ quote(PORT) }} \
        --output {{ quote(OUTPUT) }}

# Backward-compatible name for the native live watcher.
docs-serve PRODUCT PORT="8000" OUTPUT="target/alphal00p-docs-serve":
    just docs-watch {{ quote(PRODUCT) }} {{ quote(PORT) }} {{ quote(OUTPUT) }}

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
