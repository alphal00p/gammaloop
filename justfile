# Gammaloop build and development commands

# Build gammaloop Python CLI with UFO support and dev-optim profile
build-cli:
    cargo build -p gammaloop-api --bin gammaloop --features ufo_support --profile dev-optim

# Build gammaloop Python CLI in release mode
build-cli-release:
    cargo build -p gammaloop-api --bin gammaloop --features ufo_support --release

# Build gammaloop Python API with UFO support and dev-optim profile
build-api:
	maturin develop -m gammaloop-api/Cargo.toml --features=ufo_support,python_api --profile=dev-optim

# Build all workspace packages
build-all:
    cargo build --workspace

# Clean build artifacts
clean:
    cargo clean

# Check code without building
check:
    cargo check

# Format code
fmt:
    cargo fmt

# Run clippy linter
clippy:
    cargo clippy


# List all tests including unit tests (verbose output)
list-all-tests:
    cargo nextest list

# Run integration test by name (searches across all test files)
# Example: just test-integration evaluate_1l_scalar_vacuum
test TEST_NAME:
    cargo nextest run {{TEST_NAME}}

# Run integration tests with nocapture (shows output during test execution)
test-verbose TEST_NAME:
    cargo nextest run {{TEST_NAME}} --nocapture

# Run gammaloop
run *ARGS:
    cargo run -p gammaloop-api --bin gammaloop --features ufo_support --profile dev-optim -- {{ARGS}}
