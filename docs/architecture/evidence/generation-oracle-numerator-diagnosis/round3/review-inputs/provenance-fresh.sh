#!/usr/bin/env bash
# usage: fresh.sh <route-name> <workspace-dir> <expected-elf> ; runs dump + no-dump pristine assertions
set -u
B=/shared/lcnbr/dev/gammaloop
R=$B/generation-oracle-numerator-diagnosis/.local/review-oracle-qed-provenance-opus-round2-20261029184414-1303495
name=$1; w=$2; elf=$3; O=$R/runs/$name; mkdir -p $O/dump
if [ -z "${SYMBOLICA_LICENSE:-}" ]; then export SYMBOLICA_LICENSE="$(<"${XDG_CONFIG_HOME:-$HOME/.config}/symbolica/license")"; fi
unset FEYNKIT_ORACLE_RECORD FEYNKIT_ORACLE_DUMP FEYNKIT_ORACLE_THREADS FEYNKIT_ORACLE_DIAGNOSTIC CARGO_TARGET_DIR
cd $B/$w
{ echo "start $(date -u +%FT%TZ)"; echo "pwd $PWD"; echo "elf-before $(sha256sum $elf)"; env | grep -E '^(FEYNKIT_|CARGO_|RUST)' | sort; echo "license-set ${SYMBOLICA_LICENSE:+yes}"; } > $O/meta.txt
( FEYNKIT_ORACLE_THREADS=1 FEYNKIT_ORACLE_DUMP=$O/dump nix develop -c sh -c 'env | grep -E "^(FEYNKIT_|CARGO_TARGET|RUSTFLAGS)" | sort; exec cargo test --locked --profile dev-optim -p feynkit-generator --test generation_oracle -- --test-threads 1 --nocapture xs_qed_eedd_L1 xs_qed_eemm_L1 xs_qed_eeee_1pi_L1' ) > $O/dump.log 2>&1
echo "dump-exit $?" >> $O/meta.txt
( FEYNKIT_ORACLE_THREADS=1 nix develop -c sh -c 'env | grep -E "^(FEYNKIT_|CARGO_TARGET|RUSTFLAGS)" | sort; exec cargo test --locked --profile dev-optim -p feynkit-generator --test generation_oracle -- --test-threads 1 --nocapture xs_qed_eedd_L1 xs_qed_eemm_L1 xs_qed_eeee_1pi_L1' ) > $O/nodump.log 2>&1
echo "nodump-exit $?" >> $O/meta.txt
nix develop -c rustc -vV > $O/rustc.txt 2>&1
{ echo "elf-after $(sha256sum $elf)"; grep -o 'Running [^ ]* ([^)]*)' $O/dump.log $O/nodump.log; echo "end $(date -u +%FT%TZ)"; } >> $O/meta.txt
