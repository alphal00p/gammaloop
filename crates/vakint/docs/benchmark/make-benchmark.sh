#!/usr/bin/env bash
# Regenerate crates/vakint/docs/benchmark.pdf from an external-data root and a benchmark binary.
#
#   make-benchmark.sh DATA_ROOT BENCHMARK_BINARY OUT_DIR [--install] [--draft] [run_benchmarks.py options...]
#
# DATA_ROOT         vakint-external-data checkout (1L-5L packs, with master values)
# BENCHMARK_BINARY  the benchmark_suite example of this crate, built with
#                   cargo build --release -p vakint --features feynkit-ingress,experimental-rustred --example benchmark_suite
# OUT_DIR           measurements (resumable), results.json, graphs.json, benchmark.pdf, steps.json
# --install         copy results/graphs to docs/benchmark/data/{suite,graphs}.json, the PDF to docs/benchmark.pdf
#                   and write docs/benchmark/data/receipt.json (identities of data, binary, source and PDF)
# --draft           DRAFT watermark
# Other options go to tools/run_benchmarks.py (--cpus, --skip-heavy, --heavy-only, --only, --repetitions,
# --heavy-n-cores, --heavy-cpus, --form, --assemble-only). Needs python3, typst (0.15) and, for the
# FORM backends, form on PATH (or --form).
set -euo pipefail
if (( $# < 3 )); then sed -n '2,19p' "$0"; exit 2; fi
ROOT=$(realpath "$1"); BIN=$(realpath "$2"); OUT=$3; shift 3
mkdir -p "$OUT"; OUT=$(realpath "$OUT")
HERE=$(cd "$(dirname "$0")" && pwd); DOCS=$(dirname "$HERE")
INSTALL=0; DRAFT=false; PASS=()
for a in "$@"; do
  case $a in
    --install) INSTALL=1;;
    --draft) DRAFT=true; PASS+=(--draft);;
    *) PASS+=("$a");;
  esac
done
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
  cp "$OUT/results.json" "$HERE/data/suite.json"
  cp "$OUT/graphs.json" "$HERE/data/graphs.json"
  # Rebuild from the installed copies so the committed PDF is reproducible from the repository.
  (cd "$DOCS" && typst compile --input draft=$DRAFT benchmark.typ benchmark.pdf)
  python3 -I - "$OUT" "$HERE/data/receipt.json" "$DOCS/benchmark.pdf" "$ROOT" <<'EOF'
import json, sys, hashlib, datetime
out, receipt, pdf, root = sys.argv[1:5]
ident = json.load(open(out + "/identity.json"))
h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
json.dump({"schema": "vakint.benchmark-receipt.v1",
           "installed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
           "data_root": root, "root_manifest_sha256": ident.get("root_manifest_sha256"),
           "loop_manifest_sha256": ident.get("loop_manifest_sha256"),
           "executable_sha256": ident.get("executable_sha256"), "harness_source_sha256": ident.get("harness_source_sha256"),
           "vakint_version": ident.get("vakint_version"),
           "suite_sha256": h(out + "/results.json"), "graphs_sha256": h(out + "/graphs.json"), "pdf_sha256": h(pdf)},
          open(receipt, "w"), indent=1)
EOF
  echo "installed docs/benchmark.pdf and docs/benchmark/data/{suite,graphs,receipt}.json"
fi
