#!/usr/bin/env bash
# Regenerate crates/vakint/docs/benchmark.pdf from an external-data root and a benchmark binary.
#
#   make-benchmark.sh DATA_ROOT BENCHMARK_BINARY OUT_DIR [--install [--accept-pending FILE]] [--draft]
#                     [run_benchmarks.py options...]
#
# DATA_ROOT         vakint-external-data checkout (1L-5L packs, with master values)
# BENCHMARK_BINARY  the benchmark_suite example of this crate, built (with its build identity) by
#                   docs/benchmark/build-benchmark-binary.sh BIN_DIR   (release profile, fat LTO)
# OUT_DIR           measurements (resumable), results.json, graphs.json, benchmark.pdf, steps.json
# --install         only if tools/receipt.py gate passes (not a draft, no pending cell unless FILE is a
#                   user-approved exception, clean identified release build of this tree, same tools, cases
#                   and data): copy results/graphs to docs/benchmark/data/{suite,graphs}.json, the PDF to
#                   docs/benchmark.pdf, write docs/benchmark/data/receipt.json. publish.sh then runs
#                   tools/receipt.py check --data-root <published data tree>.
# --draft           DRAFT watermark (never installable)
# Other options go to tools/run_benchmarks.py (--cpus, --skip-heavy, --heavy-only, --only, --cells, --repetitions,
# --max-sibling-busy, --max-run-delay, --max-foreign-cpu, --contention-retries, --allow-shared-l3,
# --heavy-n-cores, --heavy-cpus, --form, --assemble-only). Needs python3, typst (0.15) and, for the
# FORM backends, form on PATH (or --form).
set -euo pipefail
if (( $# < 3 )); then awk 'NR > 1 && /^#/ { print; next } NR > 1 { exit }' "$0"; exit 2; fi
ROOT=$(realpath "$1"); BIN=$(realpath "$2"); OUT=$3; shift 3
mkdir -p "$OUT"; OUT=$(realpath "$OUT")
HERE=$(cd "$(dirname "$0")" && pwd); DOCS=$(dirname "$HERE")
INSTALL=0; DRAFT=false; PASS=(); ACCEPT=()
while (( $# )); do
  case $1 in
    --install) INSTALL=1;;
    --accept-pending) ACCEPT=(--accept-pending "$(realpath "$2")"); shift;;
    --draft) DRAFT=true; PASS+=(--draft);;
    *) PASS+=("$1");;
  esac
  shift
done
if (( INSTALL )) && [[ $DRAFT == true ]]; then echo "--install refuses a --draft run" >&2; exit 2; fi
STEPS=$OUT/steps.json
step() { # name start-epoch
  python3 -I - "$STEPS" "$1" "$2" "$EPOCHREALTIME" <<'EOF'
import json, sys, datetime
p, name, sec = sys.argv[1], sys.argv[2], float(sys.argv[4]) - float(sys.argv[3])
try: d = json.load(open(p))
except Exception: d = {"steps": []}
d["steps"].append({"step": name, "seconds": sec, "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")})
json.dump(d, open(p, "w"), indent=1)
EOF
}
t=$EPOCHREALTIME
python3 -I "$HERE/tools/run_benchmarks.py" --binary "$BIN" --data-root "$ROOT" --out-dir "$OUT" "${PASS[@]}"
step measure "$t"
t=$EPOCHREALTIME
python3 -I "$HERE/tools/reduction_graphs.py" "$ROOT" "$OUT/graphs.json"
step graphs "$t"
t=$EPOCHREALTIME
typst compile --root / --input results="$OUT/results.json" --input graphs="$OUT/graphs.json" \
  --input draft=$DRAFT "$DOCS/benchmark.typ" "$OUT/benchmark.pdf"
step typst "$t"
echo "wrote $OUT/benchmark.pdf"
if (( INSTALL )); then
  python3 -I "$HERE/tools/receipt.py" gate --out "$OUT" --data-root "$ROOT" "${ACCEPT[@]}"
  cp "$OUT/results.json" "$HERE/data/suite.json"
  cp "$OUT/graphs.json" "$HERE/data/graphs.json"
  # Rebuild from the installed copies so the committed PDF is reproducible from the repository.
  (cd "$DOCS" && typst compile --input draft=$DRAFT benchmark.typ benchmark.pdf)
  python3 -I "$HERE/tools/receipt.py" write --out "$OUT" --data-root "$ROOT" --pdf "$DOCS/benchmark.pdf" \
    --receipt "$HERE/data/receipt.json" "${ACCEPT[@]}"
  echo "installed docs/benchmark.pdf and docs/benchmark/data/{suite,graphs,receipt}.json"
fi
