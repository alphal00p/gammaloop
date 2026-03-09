# Aliases for common commands

alias b := build
alias t := test
alias ta := test-all
alias c := check
alias ca := check-all
alias f := fmt
alias cl := clippy
alias w := wasm

default: wasm

# Build WASM target for linnest
wasm:
    cargo build -p linnest --release --target wasm32-unknown-unknown --features custom

prepare-templates:
    just wasm
    cp target/wasm32-unknown-unknown/release/linnest.wasm clinnet/templates/linnest.wasm

# Build all packages
build:
    just prepare-templates
    cargo build --all

build-clinnet:
    just prepare-templates
    cargo build -p clinnet

# Build release version
build-release:
    cargo build --release --all --features serde

# Run tests for linnet package (default)
test:
    cargo test -p linnet --features serde

# Run tests for linnest package only
test-linnest:
    cargo test -p linnest

# Run all tests with features
test-all:
    cargo test --all --features serde

# Quick check compilation for linnet (default)
check:
    cargo check -p linnet --features serde

# Check linnest package
check-linnest:
    cargo check -p linnest

# Check all packages with features
check-all:
    cargo check --all --features serde

# Run binary with arguments
run *args:
    cargo run --bin linnet-cli -- {{ args }}

# Format all code
fmt:
    cargo fmt --all

# Run clippy linter with strict warnings
clippy:
    cargo clippy --all --features serde -- -D warnings

# Clean build artifacts
clean:
    cargo clean

# Build and open documentation
docs:
    cargo doc --all --features serde --open

# Check documentation builds without warnings
check-docs:
    cargo doc --all --features serde --no-deps

# Run benchmarks (default)
bench:
    cargo bench

# Run benchmarks with custom getrandom
bench-custom:
    cargo bench --features linnest/custom

# Run specific benchmark
bench-layout:
    cargo bench layout

# Development workflow: format, clippy, test
dev: fmt clippy test-all

# CI workflow: check, test, docs
ci: check-all test-all check-docs

# Insta snapshot commands
insta-review:
    cargo insta review

insta-accept:
    cargo insta accept

insta-test:
    cargo insta test

# Full release preparation
release-prep: fmt clippy test-all check-docs wasm

# =============================================================================
# Agent Notes
# =============================================================================

# Open agent notes for editing
notes:
    ${EDITOR:-vim} .agent/notes.md

# Open agent knowledge base
knowledge:
    ${EDITOR:-vim} .agent/README.md

# Quick note - add timestamped entry to notes
note text:
    echo "\n### $(date '+%Y-%m-%d %H:%M')" >> .agent/notes.md
    echo "{{ text }}" >> .agent/notes.md
    echo "" >> .agent/notes.md
