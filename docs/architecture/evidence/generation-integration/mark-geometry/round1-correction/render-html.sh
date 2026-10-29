set -uo pipefail
source target/mark-geometry-corrected/environment.sh
for note in generation-mark-geometry generation-drawing-ports generation-drawing-notices kurvst-arrowhead-engine nix-crane-cache-reuse; do
 printf '#document("index.html")[#include "/docs/architecture/%s.typ"]\n' "$note" > "target/mark-geometry-corrected/$note-wrapper.typ"
 typst compile --features html,bundle --format bundle --root . --creation-timestamp 0 "target/mark-geometry-corrected/$note-wrapper.typ" "target/mark-geometry-corrected/$note-html" > "target/mark-geometry-corrected/html-$note.log" 2>&1
 echo "$?" > "target/mark-geometry-corrected/html-$note.exit"
done
typst compile --features html,bundle --format bundle --root . --creation-timestamp 0 target/mark-geometry-corrected/generation-mark-geometry-wrapper.typ target/mark-geometry-corrected/note-html-reproduced > target/mark-geometry-corrected/note-html-reproduction.log 2>&1 && diff -r target/mark-geometry-corrected/generation-mark-geometry-html target/mark-geometry-corrected/note-html-reproduced >> target/mark-geometry-corrected/note-html-reproduction.log 2>&1
echo "$?" > target/mark-geometry-corrected/note-html-reproduction.exit
