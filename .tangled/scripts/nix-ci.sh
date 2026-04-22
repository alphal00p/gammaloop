#!/usr/bin/env bash
set -euo pipefail

run_nix() {
  NIX_CONFIG=$'experimental-features = nix-command flakes' nix "$@"
}

current_system() {
  run_nix eval --impure --raw --expr builtins.currentSystem
}

announce() {
  printf '\n==> %s (%s)\n' "$1" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
}

run_step() {
  local name="$1"
  shift
  announce "$name"
  "$@"
}

full_ci() {
  local system
  local partition_count
  local partition

  system="$(current_system)"
  partition_count="$(run_nix eval --raw ".#ci.${system}.partitionCount")"

  run_step "Build cargo artifacts" \
    run_nix build .#cargoArtifacts --print-build-logs

  run_step "Build nextest archive" \
    run_nix build --impure ".#packages.${system}.gammaloop-nextest-archive" --print-build-logs

  run_step "Check formatting" \
    run_nix build ".#checks.${system}.gammaloop-fmt" --print-build-logs

  run_step "Run clippy" \
    run_nix build ".#checks.${system}.gammaloop-clippy" --print-build-logs

  for partition in $(seq 1 "${partition_count}"); do
    run_step "Run nextest partition ${partition}/${partition_count}" \
      run_nix build --impure ".#checks.${system}.gammaloop-nextest-partition-${partition}" --print-build-logs
  done

  run_step "Run doctests" \
    run_nix develop --command cargo test --doc --verbose

  run_step "Build documentation" \
    run_nix build ".#checks.${system}.gammaloop-doc" --print-build-logs

  run_step "Generate coverage report" \
    run_nix build ".#packages.${system}.gammaloop-llvm-coverage" --print-build-logs
}

build_package() {
  run_step "Build main package" \
    run_nix build .#default --print-build-logs
}

main() {
  case "${1:-}" in
    full-ci)
      full_ci
      ;;
    package)
      build_package
      ;;
    *)
      echo "usage: $0 {full-ci|package}" >&2
      exit 1
      ;;
  esac
}

main "$@"
