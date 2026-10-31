#!/usr/bin/env bash
# Build a Hestia row, then execute licensed runners even when their build is cached.
set -euo pipefail

result_dir=target/ci-results
mkdir -p "$result_dir"
rm -f "$result_dir/metrics.csv"
rm -rf "$result_dir/reports"
artifact_name=${CHECK_ATTR//[^a-zA-Z0-9_.-]/_}
echo "artifact-name=nix-metrics-$artifact_name" >> "$GITHUB_OUTPUT"
printf 'check=%s\ncommit=%s\nrun=%s/%s/actions/runs/%s\nattempt=%s\n' \
    "$CHECK_ATTR" "$GITHUB_SHA" "$GITHUB_SERVER_URL" \
    "$GITHUB_REPOSITORY" "$GITHUB_RUN_ID" "$GITHUB_RUN_ATTEMPT" \
    > "$result_dir/provenance.txt"

# Derivations contain only redistributable inputs and launchers. Test output
# is produced below as the invoking user, outside the Nix store and its hooks.
exec 3> >(tee "$result_dir/nix-build.log" >&2)
log_pid=$!
set +e
env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_PATH \
    -u SYMBOLICA_LICENSE_SIGNED -u SYMBOLICA_LICENSE_FILE nix build \
    --impure \
    --keep-going \
    --out-link "$RUNNER_TEMP/nix-check-result" \
    --print-build-logs \
    --print-out-paths \
    ".#$CHECK_ATTR" \
    > "$result_dir/output-paths.txt" 2>&3
build_status=$?
exec 3>&-
wait "$log_pid"
log_status=$?
set -e
printf 'build-exit-code=%s\n' "$build_status" >> "$result_dir/provenance.txt"

run_status=0
run_log_status=0
check_name=${CHECK_ATTR##*.}
case "$check_name" in
    gammaloop-nextest-binaries*) ;;
    alphal00p-docs|gammaloop-doctest|gammaloop-nextest|gammaloop-nextest-*)
        if [ "$build_status" -eq 0 ]; then
            mkdir -p "$result_dir/reports"
            set +e
            "$RUNNER_TEMP/nix-check-result/bin/nix-ci-check-$check_name" \
                "$PWD/$result_dir/reports" 2>&1 | tee "$result_dir/test-run.log"
            run_statuses=("${PIPESTATUS[@]}")
            set -e
            run_status=${run_statuses[0]}
            run_log_status=${run_statuses[1]}
            printf 'run-exit-code=%s\n' "$run_status" >> "$result_dir/provenance.txt"
        fi
        ;;
esac

reports=()
if [ -d "$result_dir/reports" ]; then
    while IFS= read -r -d '' report; do reports+=("$report"); done \
        < <(find "$result_dir/reports" -type f -name '*.xml' -print0)
fi
export_status=0
if [ "${#reports[@]}" -gt 0 ]; then
    # Only metadata CSV is uploaded; raw JUnit and licensed output stay local.
    python3 .github/scripts/export-test-timings.py \
        --check "$CHECK_ATTR" \
        --output "$result_dir/metrics.csv" \
        "${reports[@]}" || export_status=$?
fi
for status in "$build_status" "$run_status" "$log_status" "$run_log_status" "$export_status"; do
    if [ "$status" -ne 0 ]; then
        exit "$status"
    fi
done
