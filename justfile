# Gammaloop build and development commands

# Build gammaloop Python CLI with UFO support and dev-optim profile
build-cli:
    cargo build -p gammaloop-api --bin gammaloop --features ufo_support --profile dev-optim

# Build gammaloop CLI without Python (no pyo3)
build-cli-no-pyo3:
    cargo build -p gammaloop-api --bin gammaloop --no-default-features --features cli,no_pyo3 --profile dev-optim

# Build gammaloop Python CLI with UFO support and stable ABI (dev-optim profile)
build-cli-abi:
	cargo build -p gammaloop-api --bin gammaloop --features ufo_support,python_abi --profile dev-optim

# Build gammaloop Python CLI in release mode
build-cli-release:
	cargo build -p gammaloop-api --bin gammaloop --features ufo_support --release

# Build gammaloop Python CLI in release mode with stable ABI
build-cli-release-abi:
	cargo build -p gammaloop-api --bin gammaloop --features ufo_support,python_abi --release

# Build gammaloop Python API with UFO support and dev-optim profile
build-api:
	maturin develop -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_api --profile=dev-optim

# Build gammaloop Python API with UFO support and stable ABI (dev-optim profile)
build-api-abi:
    maturin develop -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_abi --profile=dev-optim

# Build gammaloop Python API with UFO support and dev-optim profile
build-api-wheel:
	maturin build -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_api --profile=dev-optim

# Build gammaloop Python API with UFO support and stable ABI (dev-optim profile)
build-api-abi-weel:
    maturin build -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_abi --profile=dev-optim

# Build gammaloop Python API with UFO support and release profile
build-api-release:
	maturin develop -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_api --profile=release

# Build gammaloop Python API with UFO support and stable ABI (release profile)
build-api-abi-release:
    maturin develop -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_abi --profile=release

# Build gammaloop Python API wheel with UFO support and release profile
build-api-release-wheel:
	maturin build -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_api --profile=release

# Build gammaloop Python API weel with UFO support and stable ABI (release profile)
build-api-abi-release-wheel:
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
    cargo doc --workspace --no-deps --locked

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
        cargo clippy --workspace --all-targets --locked {{ lint_args }}
    else
        cargo clippy --workspace --all-targets --locked
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
        three-dimensional-reps
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

    existing_rustflags="${RUSTFLAGS-}"

    if [ "$enforce_warnings_as_errors" -eq 1 ]; then
        compile_rustflags="${existing_rustflags:+$existing_rustflags }-Dwarnings"
        cmd=(
            env
            "RUSTFLAGS=$compile_rustflags"
            cargo nextest run
            --cargo-profile dev-optim
            -P test_gammaloop
            --run-ignored all
        )
    else
        cmd=(
            cargo nextest run
            --cargo-profile dev-optim
            -P test_gammaloop
            --run-ignored all
        )
    fi
    for package in "${gammaloop_packages[@]}"; do
        cmd+=(-p "$package")
    done
    if [ ${#filter_terms[@]} -gt 0 ]; then
        cmd+=(--ignore-default-filter)
        filterset="${filter_terms[0]}"
        for term in "${filter_terms[@]:1}"; do
            filterset="${filterset} or ${term}"
        done
        cmd+=(-E "$filterset")
    fi
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
        --cargo-profile dev-optim
        --profile ci_gammaloop
        --locked
        --no-fail-fast
        --run-ignored all
    )
    for package in "${gammaloop_packages[@]}"; do
        cmd+=(-p "$package")
    done
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
