# Gammaloop build commands

# Build gammaloop Python CLI with UFO support and dev-optim profile
build-cli:
    cargo build -p gammaloop-api --bin gammaloop --features ufo_support --profile dev-optim

# Build gammaloop Python CLI in release mode
build-cli-release:
    cargo build -p gammaloop-api --bin gammaloop --features ufo_support --release

# Build gammaloop Python API with UFO support and dev-optim profile
build-api:
	maturin develop -m gammaloop-api/Cargo.toml --features=ufo_support,python_api --profile=dev-optim

# Build gammaloop Python API in release mode
build-api-release:
	maturin develop -m gammaloop-api/Cargo.toml --features=ufo_support,python_api --release

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
