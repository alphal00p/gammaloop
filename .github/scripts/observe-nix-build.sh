#!/usr/bin/env bash
# Keep output paths on stdout. Put build diagnostics on the existing stderr log.
set -uo pipefail

report_headroom() {
    local phase=$1 path
    TZ=UTC printf '[%(%FT%TZ)T] Nix scratch headroom (%s, bytes):\n' -1 "$phase" >&2
    for path in /nix/store "${RUNNER_TEMP:-${TMPDIR:-/tmp}}" "$PWD" "$build_dir"; do
        if [ -n "$path" ] && [ -d "$path" ]; then
            df -B1 --output=avail,target -- "$path" >&2 || true
        fi
    done
}

max_jobs=$(nix config show max-jobs 2>/dev/null) || max_jobs=unknown
cores=$(nix config show cores 2>/dev/null) || cores=unknown
build_dir=$(nix config show build-dir 2>/dev/null) || build_dir=
TZ=UTC printf '[%(%FT%TZ)T] Nix scheduler max-jobs=%s cores=%s (0 cores means all available CPUs)\n' \
    -1 "$max_jobs" "$cores" >&2
report_headroom before

exec 3> >(
    while IFS= read -r line || [ -n "$line" ]; do
        TZ=UTC printf '[%(%FT%TZ)T] %s\n' -1 "$line" >&2
    done
)
logger=$!
nix "$@" 2>&3 3>&-
status=$?
exec 3>&-
wait "$logger" || true
TZ=UTC printf '[%(%FT%TZ)T] Nix build exit=%d\n' -1 "$status" >&2
report_headroom after
exit "$status"
