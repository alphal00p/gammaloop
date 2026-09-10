# Gammaloop build and development commands

ci_cargo_profile := "ci-optim"

mod linnet 'crates/linnet/Justfile'

alias sync-draw-assets := sync-drawing-assets

drawing_asset_sync := "none"

_sync-drawing-assets:
    #!/usr/bin/env bash
    set -euo pipefail
    case "{{drawing_asset_sync}}" in
      nix)
        just sync-drawing-assets
        ;;
      cargo)
        just sync-drawing-assets-cargo
        ;;
      none|off|false)
        ;;
      *)
        echo "unknown drawing_asset_sync={{drawing_asset_sync}}; expected nix, cargo, or none" >&2
        exit 2
        ;;
    esac

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
    #!/usr/bin/env bash
    if [ -n "{{ lint_args }}" ]; then
        cargo fmt --all {{ lint_args }}
    else
        cargo fmt --all
    fi

# Run clippy linter
clippy *lint_args:
    #!/usr/bin/env bash
    if [ -n "{{ lint_args }}" ]; then
        cargo clippy --workspace --all-targets --locked --profile {{ ci_cargo_profile }} {{ lint_args }}
    else
        cargo clippy --workspace --all-targets --locked --profile {{ ci_cargo_profile }}
    fi

# Run workspace clippy via Nix (same as CI)
clippy-nix:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-clippy

# Check formatting via Nix (same as CI)
fmt-check-nix:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-fmt

# Run security audit via Nix (same as CI)
audit-nix:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-audit

# Run license check via Nix (same as CI)
deny-nix:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-deny

# Build documentation via Nix (same as CI)
doc-nix:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-doc

# Run doctests via Nix (same as CI)
doctest-nix:
    nix build --impure .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-doctest

test:
    cargo nextest run --workspace --cargo-profile dev-optim -P local_test

test_gammaloop *args:
    #!/usr/bin/env bash
    set -euo pipefail

    enforce_warnings_as_errors=1
    gammaloop_packages=(
        gammaloop-api
        gammalooprs
        idenso
        linnest
        linnet
        spenso
        spenso-hep-lib
        spenso-macros
        vakint
        gammaloop-integration-tests
    )
    gammaloop_module_classes=(important slow failing)
    gammaloop_base_excluded_classes=(slow failing)
    raw_args=({{ args }})
    selected_classes=()
    nextest_args=()
    passthrough_mode=0

    contains_value() {
        local needle="$1"
        shift
        local value
        for value in "$@"; do
            if [ "$value" = "$needle" ]; then
                return 0
            fi
        done
        return 1
    }

    for arg in "${raw_args[@]}"; do
        if [ "$passthrough_mode" -eq 1 ]; then
            nextest_args+=("$arg")
            continue
        fi

        case "$arg" in
            --)
                passthrough_mode=1
                nextest_args+=("--")
                ;;
            --allow-warnings|--no-warnings-as-errors)
                enforce_warnings_as_errors=0
                ;;
            --fail-fast|--ff|--no-fail-fast|--nff)
                nextest_args+=("$arg")
                ;;
            *)
                selected_classes+=("$arg")
                ;;
        esac
    done

    if [ ${#selected_classes[@]} -eq 0 ]; then
        selected_classes=(base)
    fi

    want_base=0
    selected_modules=()
    declare -A seen_modules=()

    for class in "${selected_classes[@]}"; do
        if [ "$class" = "base" ]; then
            want_base=1
        elif contains_value "$class" "${gammaloop_module_classes[@]}"; then
            if [ -z "${seen_modules[$class]+x}" ]; then
                selected_modules+=("$class")
                seen_modules[$class]=1
            fi
        else
            echo "unknown test_gammaloop class: $class" >&2
            echo "available classes: base ${gammaloop_module_classes[*]}" >&2
            exit 1
        fi
    done

    filter_terms=()
    if [ "$want_base" -eq 1 ]; then
        remaining_excluded=()
        for class in "${gammaloop_base_excluded_classes[@]}"; do
            if [ -z "${seen_modules[$class]+x}" ]; then
                remaining_excluded+=("$class")
            fi
        done

        if [ ${#remaining_excluded[@]} -gt 0 ]; then
            excluded_regex="$(printf '%s|' "${remaining_excluded[@]}")"
            excluded_regex="${excluded_regex%|}"
            filter_terms+=("not test(/(^|::)(${excluded_regex})::/)")
        fi

        if [ ${#remaining_excluded[@]} -eq 0 ] && [ ${#selected_modules[@]} -gt 0 ]; then
            filter_terms=("all()")
        fi
    fi

    if ! { [ "$want_base" -eq 1 ] && [ ${#filter_terms[@]} -eq 0 ]; }; then
        for class in "${selected_modules[@]}"; do
            filter_terms+=("test(/(^|::)${class}::/)")
        done
    fi

    run_ignored=0
    for class in "${selected_modules[@]}"; do
        case "$class" in
            slow|failing)
                run_ignored=1
                ;;
        esac
    done

    existing_rustflags="${RUSTFLAGS-}"

    if [ "$enforce_warnings_as_errors" -eq 1 ]; then
        compile_rustflags="${existing_rustflags:+$existing_rustflags }-Dwarnings"
        cmd=(
            env
            "RUSTFLAGS=$compile_rustflags"
            cargo nextest run
            --cargo-profile dev-optim
            -P test_gammaloop
        )
    else
        cmd=(
            cargo nextest run
            --cargo-profile dev-optim
            -P test_gammaloop
        )
    fi
    if [ "$run_ignored" -eq 1 ]; then
        # Explicit slow/failing selections must bypass the profile filter that
        # excludes those classes from the curated base suite.
        cmd+=(--run-ignored all --ignore-default-filter)
    fi
    for package in "${gammaloop_packages[@]}"; do
        cmd+=(-p "$package")
    done
    # Keep the known ARM/SymJIT backend mismatch out of the curated suite until
    # Symbolica updates its pinned SymJIT dependency. The test remains available
    # for direct nextest invocations.
    known_broken_filter='not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)'
    if [ ${#filter_terms[@]} -gt 0 ]; then
        filterset="${filter_terms[0]}"
        for term in "${filter_terms[@]:1}"; do
            filterset="${filterset} or ${term}"
        done
        filterset="(${filterset}) and ${known_broken_filter}"
    else
        filterset="${known_broken_filter}"
    fi
    cmd+=(-E "$filterset")
    if [ ${#nextest_args[@]} -gt 0 ]; then
        cmd+=("${nextest_args[@]}")
    fi

    printf 'Running:'
    printf ' %q' "${cmd[@]}"
    printf '\n'
    "${cmd[@]}"

test-all:
    cargo nextest run --workspace --cargo-profile release -P local_test_all

# Run workspace nextest via Nix (same as CI)
test-nix-all:
    nix build --impure .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-nextest

# Run workspace nextest via Nix
test-nix:
    just test-nix-all

# Generate coverage report via Nix
coverage-nix:
    nix build .#packages.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-llvm-coverage

# Run all CI checks locally (same as CI)
ci-checks: clippy-nix fmt-check-nix audit-nix deny-nix doc-nix doctest-nix test-nix

# Run tests in release mode (faster execution)
test-release TEST_NAME="":
    #!/usr/bin/env bash
    if [ -n "{{ TEST_NAME }}" ]; then
        cargo nextest run {{ TEST_NAME }} --release --profile ci
    else
        cargo nextest run --release --profile ci
    fi

# Run tests in release mode (faster execution)
test-ci TEST_NAME="":
    just --unstable _test-ci "{{ TEST_NAME }}" ""

# Run tests in release mode on GitHub macOS runners.
test-ci-mac TEST_NAME="":
    just --unstable _test-ci "{{ TEST_NAME }}" 'not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)'

_test-ci TEST_NAME="" NEXTEST_FILTERSET="":
    #!/usr/bin/env bash
    set -euo pipefail
    gammaloop_packages=(
        gammaloop-api
        gammalooprs
        idenso
        linnest
        linnet
        spenso
        spenso-hep-lib
        spenso-macros
        vakint
        gammaloop-integration-tests
    )
    cmd=(
        cargo nextest run
        --cargo-profile {{ ci_cargo_profile }}
        --profile ci_gammaloop
        --locked
        --no-fail-fast
    )
    for package in "${gammaloop_packages[@]}"; do
        cmd+=(-p "$package")
    done
    if [ -n "{{ NEXTEST_FILTERSET }}" ]; then
        cmd+=(-E '{{ NEXTEST_FILTERSET }}')
    fi
    if [ -n "{{ TEST_NAME }}" ]; then
        cmd+=({{ TEST_NAME }})
    fi
    "${cmd[@]}"

# Run the Python API integration tests explicitly after preparing the Python env.
test-python-api:
    cargo nextest run -p gammaloop-integration-tests --features python-api-tests --test test_python_api --profile ci --locked --no-fail-fast

# Run tests in release mode with full parallelism
test-fast TEST_NAME="":
    #!/usr/bin/env bash
    if [ -n "{{ TEST_NAME }}" ]; then
        cargo nextest r {{ TEST_NAME }} --cargo-profile dev-optim
    else
        cargo nextest  r --cargo-profile dev-optim
    fi

# Build cargo dependencies via Nix (useful for caching) - uses release mode
build-deps-nix:
    nix build .#packages.$(nix eval --impure --raw --expr 'builtins.currentSystem').cargoArtifacts

# Build everything in release mode for maximum performance
build-release-all:
    cargo build --workspace --release

# Quick development cycle: build deps, then build in release
dev-release: build-deps-nix build-release-all

# List all tests including unit tests (verbose output)
list-all-tests:
    cargo nextest list

# Run integration test by name (searches across all test files)

# Example: just test evaluate_1l_scalar_vacuum
test_one TEST_NAME:
    cargo nextest run {{ TEST_NAME }}

# Run integration tests with nocapture (shows output during test execution)
test-verbose TEST_NAME:
    cargo nextest run {{ TEST_NAME }} --nocapture

# Run gammaloop
run *ARGS:
    cargo run -p gammaloop-api --bin gammaloop --features ufo_support --profile dev-optim -- {{ ARGS }}
