set -u
source target/mark-geometry-attribution/environment.sh
out=target/mark-geometry-attribution
printf '#document("index.html")[#include "/docs/architecture/generation-mark-geometry.typ"]\n' > "$out/generation-mark-geometry-wrapper.typ"
for pass in first reproduced; do
 typst compile --features html --root . --creation-timestamp 0 docs/architecture/generation-mark-geometry.typ "$out/note-$pass.pdf" > "$out/note-$pass.log" 2>&1
 echo $? > "$out/note-$pass.exit"
 typst compile --features html,bundle --format bundle --root . --creation-timestamp 0 "$out/generation-mark-geometry-wrapper.typ" "$out/note-$pass-html" > "$out/note-$pass-html.log" 2>&1
 echo $? > "$out/note-$pass-html.exit"
done
cmp "$out/note-first.pdf" "$out/note-reproduced.pdf"
echo $? > "$out/note-pdf-reproduction.exit"
diff -r "$out/note-first-html" "$out/note-reproduced-html"
echo $? > "$out/note-html-reproduction.exit"
