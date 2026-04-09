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

# Run clippy via Nix (same as CI)
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

# Run all nextest partitions locally (same as CI)
test-nix-all:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-nextest-partition-1
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-nextest-partition-2
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-nextest-partition-3
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-nextest-partition-4
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-nextest-partition-5
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-nextest-partition-6

# Run specific nextest partition via Nix
test-nix-partition PARTITION:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-nextest-partition-{{ PARTITION }}

# Run nextest (all tests) via Nix
test-nix:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-nextest

# Generate coverage report via Nix
coverage-nix:
    nix build .#packages.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-llvm-coverage

# Run all CI checks locally (release mode)
ci-checks: clippy-nix fmt-check-nix audit-nix deny-nix doc-nix test-nix

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
    if [ -n "{{ TEST_NAME }}" ]; then
       cargo nextest run {{ TEST_NAME }} --workspace --profile ci --locked --no-fail-fast {{ TEST_NAME }}
    else
        cargo nextest run --workspace --profile ci --locked --no-fail-fast
    fi

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
test TEST_NAME:
    cargo nextest run {{ TEST_NAME }}

# Run integration tests with nocapture (shows output during test execution)
test-verbose TEST_NAME:
    cargo nextest run {{ TEST_NAME }} --nocapture

# Run gammaloop
run *ARGS:
    cargo run -p gammaloop-api --bin gammaloop --features ufo_support --profile dev-optim -- {{ ARGS }}
