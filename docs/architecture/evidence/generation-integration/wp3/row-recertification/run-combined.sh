#!/usr/bin/env bash
set -euo pipefail
# Execute the exact nextestCheckFor CPU-budget calculation; keep its maximum16
# worker setting, then let the checked-in resource override limit owned parents.
nextest_threads=$(nproc)
if (( NIX_BUILD_CORES > 0 && NIX_BUILD_CORES < nextest_threads )); then
  nextest_threads=$NIX_BUILD_CORES
fi
export NEXTEST_TEST_THREADS=$((nextest_threads < 16 ? nextest_threads : 16))
export FEYNKIT_ORACLE_THREADS=1,8,32
printf 'NIX_BUILD_CORES=%s nproc=%s NEXTEST_TEST_THREADS=%s\n' "$NIX_BUILD_CORES" "$(nproc)" "$NEXTEST_TEST_THREADS" > target/wp3-units19-21-recertification/runtime-budget.log
selection='test(outputs_do_not_depend_on_index_registration_order_) or test(=index_registration_child_failure_propagates) or test(=ym_ggg_zero_color_L3) or test(=xs_qcdqed_addx_canonical_scalar_lr_L3)'
args=(--locked --cargo-profile ci-optim --profile ci_gammaloop -p feynkit-generator --test generation_oracle -E "$selection")
cargo nextest list "${args[@]}" --message-format json > target/wp3-units19-21-recertification/combined-inventory.json 2> target/wp3-units19-21-recertification/combined-inventory.log
cargo nextest show-config test-groups "${args[@]}" > target/wp3-units19-21-recertification/effective-test-groups.log 2>&1
cargo nextest run "${args[@]}" > target/wp3-units19-21-recertification/combined-normal.log 2>&1
