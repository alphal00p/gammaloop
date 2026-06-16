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
        three-dimensional-reps
        vakint
        gammaloop-integration-tests
    )
    gammaloop_module_classes=(important slow failing)
    gammaloop_base_excluded_classes=(slow failing)
    raw_args=({{ args }})
    cargo_profile_args=(--cargo-profile dev-optim)
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
            --release)
                cargo_profile_args=(--release)
                ;;
            --test=*|--fail-fast|--ff|--no-fail-fast|--nff)
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
            "${cargo_profile_args[@]}"
            -P test_gammaloop
        )
    else
        cmd=(
            cargo nextest run
            "${cargo_profile_args[@]}"
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
    # Slow/failing selections bypass the profile filter and include ignored
    # tests. Physical PySecDec validation remains an explicit manual command.
    filterset="(${filterset}) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))"
    cmd+=(-E "$filterset")
    if [ ${#nextest_args[@]} -gt 0 ]; then
        cmd+=("${nextest_args[@]}")
    fi

    printf 'Running:'
    printf ' %q' "${cmd[@]}"
    printf '\n'
    "${cmd[@]}"

# Run all 166 scalar LU cross-section cases, including slow cases, in release mode.
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
        just test_gammaloop base slow --release --test=test_runs --fail-fast -- \
        scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects:: \
        scalar_3l_cross_section_inspects::slow::

test-all:
    cargo nextest run --workspace --cargo-profile release -P local_test_all

# Run tests in release mode (faster execution)
test-release TEST_NAME="":
    cargo nextest run {{ TEST_NAME }} --release --profile ci

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
        three-dimensional-reps
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
    cargo nextest r {{ TEST_NAME }} --cargo-profile dev-optim

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
