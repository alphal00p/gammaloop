set -uo pipefail
source target/mark-geometry-corrected/environment.sh
python3 target/mark-geometry-corrected/check-docs.py > target/mark-geometry-corrected/registry.log 2>&1
echo "$?" > target/mark-geometry-corrected/registry.exit
for note in generation-mark-geometry generation-drawing-ports generation-drawing-notices kurvst-arrowhead-engine nix-crane-cache-reuse; do
 typst compile --features html --root . --creation-timestamp 0 "docs/architecture/$note.typ" "target/mark-geometry-corrected/$note.pdf" > "target/mark-geometry-corrected/docs-$note.log" 2>&1
 echo "$?" > "target/mark-geometry-corrected/docs-$note.exit"
done
just check-symbolica-feature-isolation > target/mark-geometry-corrected/kernel.log 2>&1
echo "$?" > target/mark-geometry-corrected/kernel.exit
rustc -Vv > target/mark-geometry-corrected/rustc-version.log
cargo -V > target/mark-geometry-corrected/cargo-version.log
typst --version > target/mark-geometry-corrected/typst-version.log
