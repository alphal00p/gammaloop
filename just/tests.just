# Ordinary local and CI recipes share the strict fast policy in nextest.toml.
[positional-arguments]
test *args:
    just _test-ci "$@"

# Classes narrow the fast lane; selecting `slow` does not enable ignored tests
# or the long cases reserved for main. With no class, select the whole fast lane.
[positional-arguments]
test_gammaloop *args:
    #!/usr/bin/env bash
    set -euo pipefail
    raw_args=("$@")
    cargo_profile_args=(--cargo-profile {{ ci_cargo_profile }})
    filter_terms=()
    nextest_args=()
    passthrough_mode=0
    for arg in "${raw_args[@]}"; do
        if [ "$passthrough_mode" -eq 1 ]; then
            nextest_args+=("$arg")
            continue
        fi
        case "$arg" in
            --)
                passthrough_mode=1
                nextest_args+=(--)
                ;;
            base)
                filter_terms+=('not test(/(^|::)(slow|failing)::/)')
                ;;
            important|slow|failing)
                filter_terms+=("test(/(^|::)${arg}::/)")
                ;;
            --release)
                cargo_profile_args=(--release)
                ;;
            --test=*|--fail-fast|--ff|--no-fail-fast|--nff)
                nextest_args+=("$arg")
                ;;
            *)
                echo "unknown test_gammaloop selection: $arg" >&2
                echo 'classes: base important slow failing; test filters follow --' >&2
                exit 1
                ;;
        esac
    done
    cmd=(cargo nextest run --workspace "${cargo_profile_args[@]}"
         --profile ci_gammaloop --locked --no-fail-fast)
    if [ ${#filter_terms[@]} -gt 0 ]; then
        filterset="${filter_terms[0]}"
        for term in "${filter_terms[@]:1}"; do
            filterset="$filterset or $term"
        done
        cmd+=(-E "$filterset")
    fi
    cmd+=("${nextest_args[@]}")
    "${cmd[@]}"

# Run the scalar LU cases admitted to the strict fast lane.
test_LU_scalar_xs:
    #!/usr/bin/env bash
    set -euo pipefail
    export CARGO_BUILD_JOBS=1 NEXTEST_TEST_THREADS=1 NEXTEST_RETRIES=0
    export NEXTEST_STATUS_LEVEL="${NEXTEST_STATUS_LEVEL:-pass}"
    export RUST_MIN_STACK=67108864 INSTA_UPDATE=no
    # Keep the matrix on its normal assertion, cleanup and evaluator paths.
    unset GL_ALL_LOG_FILTER INSTA_FORCE_PASS GAMMALOOP_TESTS_NO_CLEAN_STATE
    unset SPENSO_NETWORK_PROFILE SPENSO_NETWORK_PROFILE_VERBOSE SPENSO_NETWORK_PROFILE_ATOM_BYTES
    unset SPENSO_NETWORK_LAZY_TENSOR_SUMS SPENSO_NETWORK_MAX_LAZY_DISTRIBUTED_TERMS
    export GL_DISPLAY_FILTER=off GL_LOGFILE_FILTER=off
    mkdir -p target/test_LU_scalar_xs
    run_dir="$(mktemp -d target/test_LU_scalar_xs/run.XXXXXX)"
    exec python3 bin/ram_watchdog.py --limit-gb 30 \
        --log "$run_dir/watchdog.jsonl" -- \
        just test_gammaloop --test=test_runs --fail-fast -- \
        scalar_3l_cross_section_inspects::

test-all:
    just _test-ci

# Explicit release compilation retains the same runtime selection and cap.
[positional-arguments]
test-release *args:
    cargo nextest run --workspace --release --profile ci_gammaloop --locked --no-fail-fast "$@"

[positional-arguments]
test-ci *args:
    just _test-ci "$@"

# The SymJIT backend mismatch is platform-specific; Linux runs this comparison.
[positional-arguments]
test-ci-mac *args:
    just _test-ci -E 'not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)' "$@"

[positional-arguments]
_test-ci *args:
    cargo nextest run --workspace --cargo-profile {{ ci_cargo_profile }} --profile ci_gammaloop --locked --no-fail-fast "$@"

# Run the Python API integration tests explicitly after preparing the Python env.
test-python-api:
    cargo nextest run -p gammaloop-integration-tests --cargo-profile {{ ci_cargo_profile }} --features python-api-tests --test test_python_api --profile ci_gammaloop --locked --no-fail-fast

[positional-arguments]
test-fast *args:
    just _test-ci "$@"

# Inventory includes skipped and ignored tests; listing does not run their bodies.
list-all-tests:
    cargo nextest list --workspace --cargo-profile {{ ci_cargo_profile }} --locked --profile ci_full

# Example: just test_one evaluate_1l_scalar_vacuum
[positional-arguments]
test_one TEST_NAME:
    just _test-ci "$@"

# Nocapture serializes the selected fast tests and retains their timeout.
[positional-arguments]
test-verbose TEST_NAME:
    just _test-ci "$@" --nocapture
