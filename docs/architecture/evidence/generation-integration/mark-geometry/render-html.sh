set -uo pipefail
source target/mark-geometry/environment.sh
for note in generation-mark-geometry generation-drawing-ports generation-drawing-notices kurvst-arrowhead-engine nix-crane-cache-reuse; do
 printf '#document("index.html")[#include "/docs/architecture/%s.typ"]\n' "$note" > "target/mark-geometry/$note-wrapper.typ"
 typst compile --features html,bundle --format bundle --root . --creation-timestamp 0 "target/mark-geometry/$note-wrapper.typ" "target/mark-geometry/$note-html" > "target/mark-geometry/html-$note.log" 2>&1
 echo "$?" > "target/mark-geometry/html-$note.exit"
done
typst compile --features html,bundle --format bundle --root . --creation-timestamp 0 target/mark-geometry/generation-mark-geometry-wrapper.typ target/mark-geometry/note-html-reproduced > target/mark-geometry/note-html-reproduction.log 2>&1 && diff -r target/mark-geometry/generation-mark-geometry-html target/mark-geometry/note-html-reproduced >> target/mark-geometry/note-html-reproduction.log 2>&1
echo "$?" > target/mark-geometry/note-html-reproduction.exit
