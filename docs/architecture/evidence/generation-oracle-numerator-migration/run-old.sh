#!/usr/bin/env bash
set -uo pipefail
pool="$1"
dump="$2"
unset FEYNKIT_ORACLE_RECORD FEYNKIT_ORACLE_DIAGNOSTIC UPDATE_EXPECT EXPECT FORCE_UPDATE INSTA_UPDATE
export FEYNKIT_ORACLE_THREADS="$pool"
export FEYNKIT_ORACLE_DUMP="$dump"
if [[ -z "${SYMBOLICA_LICENSE:-}" ]]; then
  SYMBOLICA_LICENSE="$(<"${XDG_CONFIG_HOME:-$HOME/.config}/symbolica/license")"
  export SYMBOLICA_LICENSE
fi
printf 'source=%s\n' '16bfd2c84ea1fdeba0c8087cc6cdc8f1082f11a8'
printf 'pool=%s\n' "$FEYNKIT_ORACLE_THREADS"
printf 'dump=%s\n' "$FEYNKIT_ORACLE_DUMP"
printf 'record=unset diagnostic=unset update=unset\n'
printf 'elf=%s\n' 'target/oracle-numerator-migration/old/generation_oracle'
printf 'argv=%s\n' '--test-threads 1 --nocapture xs_qed_eedd_L1 xs_qed_eemm_L1 xs_qed_eeee_1pi_L1'
mkdir -p "$dump"
target/oracle-numerator-migration/old/generation_oracle \
  --test-threads 1 --nocapture \
  xs_qed_eedd_L1 xs_qed_eemm_L1 xs_qed_eeee_1pi_L1
status=$?
printf 'exit=%s\n' "$status"
exit "$status"
