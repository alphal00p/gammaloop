#!/usr/bin/env bash
# Build the benchmark_suite harness with its build identity embedded and freeze a copy.
#
#   build-benchmark-binary.sh OUT_DIR [--profile release|ci-optim|...]
#
# Embeds VAKINT_BENCH_BUILD_COMMIT (git HEAD), VAKINT_BENCH_BUILD_DIRTY (modified or untracked paths
# outside crates/vakint/docs, which never enter the binary) and VAKINT_BENCH_BUILD_PROFILE, then copies
# the binary to OUT_DIR/benchmark_suite-<profile>-<commit12>[-dirty<N>] (+ .sha256). The release profile
# (fat LTO) is the one documented for timings; make-benchmark.sh --install refuses a dirty or unidentified build.
set -euo pipefail
if (( $# < 1 )); then awk 'NR > 1 && /^#/ { print; next } NR > 1 { exit }' "$0"; exit 2; fi
OUT=$1; shift
PROFILE=release
while (( $# )); do
  case $1 in
    --profile) PROFILE=$2; shift 2;;
    *) echo "unknown option $1" >&2; exit 2;;
  esac
done
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel)
COMMIT=$(git -C "$REPO" rev-parse HEAD)
DIRTY=$(git -C "$REPO" status --porcelain --untracked-files=normal -- . ':(exclude)crates/vakint/docs' | wc -l)
export VAKINT_BENCH_BUILD_COMMIT=$COMMIT VAKINT_BENCH_BUILD_DIRTY=$DIRTY VAKINT_BENCH_BUILD_PROFILE=$PROFILE
echo "building benchmark_suite at $COMMIT (dirty paths outside docs: $DIRTY), profile $PROFILE"
(cd "$REPO" && cargo build --locked --profile "$PROFILE" -p vakint --features feynkit-ingress,experimental-rustred \
  --example benchmark_suite)
DIR=$PROFILE; [[ $PROFILE == dev ]] && DIR=debug
TARGET=${CARGO_TARGET_DIR:-$REPO/target}
SRC=$TARGET/$DIR/examples/benchmark_suite
mkdir -p "$OUT"
NAME=benchmark_suite-$PROFILE-${COMMIT:0:12}$( (( DIRTY > 0 )) && echo "-dirty$DIRTY" || true)
cp "$SRC" "$OUT/$NAME.tmp" && mv "$OUT/$NAME.tmp" "$OUT/$NAME"
(cd "$OUT" && sha256sum "$NAME" > "$NAME.sha256")
"$OUT/$NAME" identity --data-root / --out "$OUT/$NAME.identity.json" 2>/dev/null || true
echo "$OUT/$NAME"
