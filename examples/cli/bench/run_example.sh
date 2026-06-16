#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

./bench.py \
  --profile release \
  --graphs G1 G1_no_num H1 H1_no_num \
  --compare-eval \
  --profiles cff optimized_cff \
  --output ./outputs/run_example_output.json \
  "$@"

./bench.py \
  --profile release \
  --graphs G1 G1_no_num H1 H1_no_num \
  --compare-eval \
  --profiles cff optimized_cff ltd optimized_ltd \
  --3drep \
  --output ./outputs/run_example_output_3drep.json \
  "$@"
