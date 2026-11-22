#!/usr/bin/env bash
# MATRIX-R45 native cell: run_cell.sh PROFILE VARIANT MODE  -> $M/runs/PROFILE/VARIANT/MODE
# (inputs from gen_matrix.py under $M/inputs, or the review's native inputs as astra-<name>; classify with classify.py)
# M, BIN_PR19, BIN_U7B, PRODENV, ASTRA override the defaults below (the 2026-11-21 run).
#   PROFILE pr19 | u7bk3c (U7b, uncontracted check:root,child,ir-hard + strict, the rv-a-L3-check profile) | u7bk3i (U7b, check:ir-hard + strict: K3 only) | u7bk3on (U7b, on:ir-hard)
#   MODE    generic (GAMMALOOP_IR_HARD=off) | child | strict (child + GAMMALOOP_IR_HARD_DEPTH_STRICT=1) | guard (child + GAMMALOOP_IR_HARD_CHILD_GUARD=1)
# every runtime variable of r6-nf1.env (supervisor keys excluded), no ambient GAMMALOOP_/VAKINT_/IDENSO_/BETA_, one core, memgate 2 GiB.
set -u
M=${M:-/common/build-cache/ben/consol/irh-depth/matrix}
P=$1; V=$2; MODE=$3
case $P in
  pr19) BIN=${BIN_PR19:-/common/build-cache/ben/consol/irh-depth/round6/bin/irhr6-pr19-33a9355352-beta-fdb34991};;
  u7bk3c|u7bk3i|u7bk3on) BIN=${BIN_U7B:-/common/build-cache/ben/consol/irh-depth/round6/bin/irhr6-u7b-981f92ac93-beta-fdb34991};;
  *) echo "profile $P"; exit 2;;
esac
ENV=${PRODENV:-/common/ben/beta/results/five-loop-r5/env/r6-nf1.env}
O=$M/runs/$P/$V/$MODE; [ -e $O ] && { echo "exists: $O"; exit 2; }; mkdir -p $O
INPUT=$M/inputs/$V; [ -d $INPUT ] || INPUT=${ASTRA:-/shared/vhirschi/t3code/astra_5L/review-20261121/IRH-DEPTH/native/inputs}/${V#astra-}
[ -d $INPUT/L3 ] || { echo "no input $V"; exit 2; }
for v in $(env | grep -oE '^(GAMMALOOP|VAKINT|IDENSO|BETA)_[A-Z0-9_]*' | grep -v '^GAMMALOOP_DYNDRV_STORE_CONFIG$'); do unset $v; done
while IFS= read -r kv; do export "$kv"; done < <(grep -E '^[A-Z_]+=' $ENV | grep -vE '^(BIN|ROOT|LIST|WORKERS|JOBS|CHUNK|UNIT_MEMORY_MAX|WAIT_FOR|BETA_PROCESS_RSS_LIMIT_GIB)=')
case $P in
  u7bk3c) export GAMMALOOP_UNCONTRACTED_NUMERATOR=check:root,child,ir-hard GAMMALOOP_UNCONTRACTED_STRICT=1;;
  u7bk3i) export GAMMALOOP_UNCONTRACTED_NUMERATOR=check:ir-hard GAMMALOOP_UNCONTRACTED_STRICT=1;;
  u7bk3on) export GAMMALOOP_UNCONTRACTED_NUMERATOR=on:ir-hard;;
esac
case $MODE in
  generic) export GAMMALOOP_IR_HARD=off;;
  child) export GAMMALOOP_IR_HARD=child;;
  strict) export GAMMALOOP_IR_HARD=child GAMMALOOP_IR_HARD_DEPTH_STRICT=1;;
  guard) export GAMMALOOP_IR_HARD=child GAMMALOOP_IR_HARD_CHILD_GUARD=1;;
  *) echo "mode $MODE"; exit 2;;
esac
export BETA_PROCESS_RSS_LIMIT_GIB=2 RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
: "${SYMBOLICA_LICENSE:?SYMBOLICA_LICENSE must be in the environment}"
export SYMBOLICA_LICENSE
export BETA_MEMORY_LOG=$O/memory.jsonl
env | grep -E '^(BETA|VAKINT|MIMALLOC|RUST_MIN|GLIBC|GAMMALOOP|RAYON|IDENSO)_' | grep -v SYMBOLICA | sort > $O/env.txt
CMD=($BIN compute --from 3 --through 3 --jobs 1 --cores 1 --tensor-parallelism serial --method rstar
  --insertions extract --grouping off --delta-cache on --input-root $INPUT
  --master-cache /tmp/gammaloop-rstar-rustred-data/.vakint-external-data/v0.2.0
  --output $O/out --projection metric --colour-mode invariants --diagrams=0)
echo "${CMD[*]}" > $O/cmd.txt
echo "start $(date -u +%FT%TZ) bin_sha=$(sha256sum $BIN | cut -c1-64) input_sha=$(sha256sum $INPUT/L3/*.dot | cut -c1-64)" > $O/meta.txt
/common/ben/beta/scripts/memgate/memgate run irhm-$P-$V-$MODE 2 -- nice -n 5 timeout 900 "${CMD[@]}" > $O/log 2>&1
rc=$?
echo "end $(date -u +%FT%TZ) rc=$rc" >> $O/meta.txt
exit $rc
