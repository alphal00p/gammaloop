#!/usr/bin/env bash
# Regenerate crates/vakint/docs/benchmark.pdf from an external-data root and a benchmark binary.
#
#   make-benchmark.sh DATA_ROOT BENCHMARK_BINARY OUT_DIR [--install [--accept-pending FILE]] [--draft]
#                     [run_benchmarks.py options...]
#   make-benchmark.sh DATA_ROOT - OUT_DIR --refresh-tables [--install]
#
# DATA_ROOT         vakint-external-data checkout (1L-5L packs, with master values)
# BENCHMARK_BINARY  the benchmark_suite example of this crate, built (with its build identity) by
#                   docs/benchmark/build-benchmark-binary.sh BIN_DIR   (release profile, fat LTO)
# OUT_DIR           measurements (resumable), results.json, graphs.json, tables.json, PRECISION-TABLES.md,
#                   benchmark.pdf, steps.json
# --install         only if tools/receipt.py gate passes (not a draft, no pending cell unless FILE is a
#                   user-approved exception, clean identified release build of this tree, same tools, cases
#                   and data, precision tables regenerated from this data root): copy results/graphs/tables to
#                   docs/benchmark/data/{suite,graphs,tables}.json and PRECISION-TABLES.md, the PDF to
#                   docs/benchmark.pdf, write docs/benchmark/data/receipt.json. publish.sh then runs
#                   tools/receipt.py check --data-root <published data tree>.
# Precision tables (tools/precision_tables.py, data root only): status 'final' with --install (needs a clean git
# checkout of the data tree), else 'provisional'; PRECISION_TABLES_LABEL=<text> names the tree in their header.
# --refresh-tables  a data release after the timing run (publish.sh pin-vakint, tools/bench_refresh.sh): no
#                   measurement (BENCHMARK_BINARY may be '-'); the timings are the committed
#                   docs/benchmark/data/suite.json and graphs.json; regenerates the reduction graphs (must equal
#                   the committed ones) and the precision tables from DATA_ROOT, rebuilds the PDF. With --install:
#                   tools/receipt.py refresh (refuses unless DATA_ROOT's reduction payloads equal the measured
#                   tree's, receipt rule 2) installs tables.json, PRECISION-TABLES.md and docs/benchmark.pdf and
#                   records the tables' data tree in the receipt (the timing provenance is kept).
# --draft           DRAFT watermark (never installable)
# Other options go to tools/run_benchmarks.py (--cpus, --skip-heavy, --heavy-only, --only, --cells, --repetitions,
# --max-sibling-busy, --max-run-delay, --max-foreign-cpu, --contention-retries, --allow-shared-l3,
# --heavy-n-cores, --heavy-cpus, --form, --assemble-only). Needs python3, typst (0.15; on PATH or TYPST=<path>)
# and, for the FORM backends, form on PATH (or --form).
set -euo pipefail
if (( $# < 3 )); then awk 'NR > 1 && /^#/ { print; next } NR > 1 { exit }' "$0"; exit 2; fi
ROOT=$(realpath "$1"); BIN=$2; [[ $BIN == - ]] || BIN=$(realpath "$BIN"); OUT=$3; shift 3
mkdir -p "$OUT"; OUT=$(realpath "$OUT")
HERE=$(cd "$(dirname "$0")" && pwd); DOCS=$(dirname "$HERE")
TYPST=${TYPST:-typst}
command -v "$TYPST" > /dev/null || { echo "typst not found (put typst 0.15 on PATH or set TYPST=<path>)" >&2; exit 2; }
INSTALL=0; DRAFT=false; REFRESH=0; PASS=(); ACCEPT=()
while (( $# )); do
  case $1 in
    --install) INSTALL=1;;
    --accept-pending) ACCEPT=(--accept-pending "$(realpath "$2")"); shift;;
    --draft) DRAFT=true; PASS+=(--draft);;
    --refresh-tables) REFRESH=1;;
    *) PASS+=("$1");;
  esac
  shift
done
if (( INSTALL )) && [[ $DRAFT == true ]]; then echo "--install refuses a --draft run" >&2; exit 2; fi
if (( REFRESH )); then
  if [[ $DRAFT == true ]] || (( ${#PASS[@]} )) || (( ${#ACCEPT[@]} )); then
    echo "--refresh-tables takes only --install (no measurement options)" >&2; exit 2; fi
  for f in suite.json graphs.json receipt.json; do
    [[ -f $HERE/data/$f ]] || { echo "--refresh-tables: no committed docs/benchmark/data/$f (measure first)" >&2; exit 2; }
  done
  cp "$HERE/data/suite.json" "$OUT/results.json"
  python3 -I "$HERE/tools/reduction_graphs.py" "$ROOT" "$OUT/graphs.json"
  TSTATUS=provisional; (( INSTALL )) && TSTATUS=final
  python3 -I "$HERE/tools/precision_tables.py" --data-root "$ROOT" --json "$OUT/tables.json" --md "$OUT/PRECISION-TABLES.md" \
    --status $TSTATUS ${PRECISION_TABLES_LABEL:+--label "$PRECISION_TABLES_LABEL"}
  "$TYPST" compile --root / --input results="$OUT/results.json" --input graphs="$HERE/data/graphs.json" \
    --input tables="$OUT/tables.json" --input draft=false "$DOCS/benchmark.typ" "$OUT/benchmark.pdf"
  echo "wrote $OUT/benchmark.pdf (timings: committed suite.json; tables: $ROOT)"
  if (( INSTALL )); then
    # the gate part of 'refresh' runs before any install copy (a refusal leaves the docs untouched)
    python3 -I "$HERE/tools/receipt.py" refresh --out "$OUT" --data-root "$ROOT" --dry-run
    cp "$OUT/tables.json" "$HERE/data/tables.json"
    cp "$OUT/PRECISION-TABLES.md" "$HERE/data/PRECISION-TABLES.md"
    (cd "$DOCS" && "$TYPST" compile --input draft=false benchmark.typ benchmark.pdf)
    python3 -I "$HERE/tools/receipt.py" refresh --out "$OUT" --data-root "$ROOT" --pdf "$DOCS/benchmark.pdf" \
      --receipt "$HERE/data/receipt.json"
    echo "installed the precision tables of $ROOT and docs/benchmark.pdf (timings and their receipt kept)"
  fi
  exit 0
fi
[[ $BIN != - ]] || { echo "BENCHMARK_BINARY '-' only with --refresh-tables" >&2; exit 2; }
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
TSTATUS=provisional; (( INSTALL )) && TSTATUS=final
python3 -I "$HERE/tools/precision_tables.py" --data-root "$ROOT" --json "$OUT/tables.json" --md "$OUT/PRECISION-TABLES.md" \
  --status $TSTATUS ${PRECISION_TABLES_LABEL:+--label "$PRECISION_TABLES_LABEL"}
step tables "$t"
t=$EPOCHREALTIME
"$TYPST" compile --root / --input results="$OUT/results.json" --input graphs="$OUT/graphs.json" \
  --input tables="$OUT/tables.json" --input draft=$DRAFT "$DOCS/benchmark.typ" "$OUT/benchmark.pdf"
step typst "$t"
echo "wrote $OUT/benchmark.pdf"
if (( INSTALL )); then
  python3 -I "$HERE/tools/receipt.py" gate --out "$OUT" --data-root "$ROOT" "${ACCEPT[@]}"
  cp "$OUT/results.json" "$HERE/data/suite.json"
  cp "$OUT/graphs.json" "$HERE/data/graphs.json"
  cp "$OUT/tables.json" "$HERE/data/tables.json"
  cp "$OUT/PRECISION-TABLES.md" "$HERE/data/PRECISION-TABLES.md"
  # Rebuild from the installed copies so the committed PDF is reproducible from the repository.
  (cd "$DOCS" && "$TYPST" compile --input draft=$DRAFT benchmark.typ benchmark.pdf)
  python3 -I "$HERE/tools/receipt.py" write --out "$OUT" --data-root "$ROOT" --pdf "$DOCS/benchmark.pdf" \
    --receipt "$HERE/data/receipt.json" "${ACCEPT[@]}"
  echo "installed docs/benchmark.pdf and docs/benchmark/data/{suite,graphs,tables,receipt}.json, PRECISION-TABLES.md"
fi
