# Gammaloop build commands

# Build gammaloop API with UFO support and dev-optim profile
build-api:
    cargo build -p gammaloop-api --bin gammaloop --features ufo_support --profile dev-optim

# Build gammaloop API in release mode
build-api-release:
    cargo build -p gammaloop-api --bin gammaloop --features ufo_support --release

# Build all workspace packages
build-all:
    cargo build --workspace

# Clean build artifacts
clean:
    cargo clean

# Run tests
test:
    cargo test

# Run tests for gammaloop-api specifically
test-api:
    cargo test -p gammaloop-api

# Check code without building
check:
    cargo check

# Format code
fmt:
    cargo fmt

# Run clippy linter
clippy:
    cargo clippy

# Run gammaloop
run *ARGS:
    cargo run -p gammaloop-api --bin gammaloop --features ufo_support --profile dev-optim -- {{ARGS}}
