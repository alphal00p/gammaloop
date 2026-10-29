set -uo pipefail
source target/mark-geometry/environment.sh
python3 target/mark-geometry/check-docs.py > target/mark-geometry/registry.log 2>&1
echo "$?" > target/mark-geometry/registry.exit
for note in generation-mark-geometry generation-drawing-ports generation-drawing-notices kurvst-arrowhead-engine nix-crane-cache-reuse; do
 typst compile --features html --root . --creation-timestamp 0 "docs/architecture/$note.typ" "target/mark-geometry/$note.pdf" > "target/mark-geometry/docs-$note.log" 2>&1
 echo "$?" > "target/mark-geometry/docs-$note.exit"
 typst compile --features html,bundle --format bundle --root . --creation-timestamp 0 "docs/architecture/$note.typ" "target/mark-geometry/$note-html" >> "target/mark-geometry/docs-$note.log" 2>&1
 echo "$?" > "target/mark-geometry/html-$note.exit"
done
just check-symbolica-feature-isolation > target/mark-geometry/kernel.log 2>&1
echo "$?" > target/mark-geometry/kernel.exit
rustc -Vv > target/mark-geometry/rustc-version.log
cargo -V > target/mark-geometry/cargo-version.log
typst --version > target/mark-geometry/typst-version.log
