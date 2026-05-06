#!/usr/bin/env bash
set -euo pipefail

bench_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$bench_dir/../../../.." && pwd)"
out_dir="${OUT_DIR:-/tmp/kurvst-coil-bench}"
count="${COUNT:-250}"
repeats="${REPEATS:-5}"
turns="${TURNS:-18}"
samples="${SAMPLES_PER_PERIOD:-16}"

mkdir -p "$out_dir"

cetz_cmd=(
  typst compile
  --root "$repo_root/crates/kurvst/typst"
  --input "count=$count"
  --input "turns=$turns"
  --input "samples-per-period=$samples"
  "$bench_dir/cetz-coil.typ"
  "$out_dir/cetz-coil.pdf"
)

kurvst_cmd=(
  typst compile
  --root "$repo_root/crates/kurvst/typst"
  --input "count=$count"
  --input "turns=$turns"
  --input "samples-per-period=$samples"
  "$bench_dir/kurvst-coil.typ"
  "$out_dir/kurvst-coil.pdf"
)

printf 'coil benchmark: count=%s turns=%s samples-per-period=%s repeats=%s\n' "$count" "$turns" "$samples" "$repeats"
printf 'outputs: %s\n' "$out_dir"

if command -v hyperfine >/dev/null 2>&1; then
  hyperfine --warmup 1 --runs "$repeats" \
    --command-name "cetz coil" "${cetz_cmd[*]}" \
    --command-name "kurvst coil" "${kurvst_cmd[*]}"
else
  for name in "cetz coil" "kurvst coil"; do
    printf '\n%s\n' "$name"
    for run in $(seq 1 "$repeats"); do
      if [ "$name" = "cetz coil" ]; then
        /usr/bin/time -p "${cetz_cmd[@]}"
      else
        /usr/bin/time -p "${kurvst_cmd[@]}"
      fi 2>&1 | sed "s/^/run $run: /"
    done
  done
fi
