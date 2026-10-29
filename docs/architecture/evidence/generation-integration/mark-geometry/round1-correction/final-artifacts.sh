set -uo pipefail
source target/mark-geometry-corrected/environment.sh
cargo nextest list --locked -p kurvst -p linnest --features linnest/svg --cargo-profile dev-optim --list-type binaries-only --message-format json > target/mark-geometry-corrected/native-test-binaries.json 2> target/mark-geometry-corrected/native-test-binaries.log
echo "$?" > target/mark-geometry-corrected/native-test-binaries.exit
for theme in light dark; do
 for format in svg png pdf; do
  typst compile --root . --creation-timestamp 0 --input "theme=$theme" --format "$format" crates/kurvst/typst/examples/marks.typ "target/mark-geometry-corrected/marks-$theme-reproduced.$format" > "target/mark-geometry-corrected/gallery-$theme-$format-reproduction.log" 2>&1 && cmp "target/mark-geometry-corrected/marks-$theme.$format" "target/mark-geometry-corrected/marks-$theme-reproduced.$format"
  echo "$?" > "target/mark-geometry-corrected/gallery-$theme-$format-reproduction.exit"
 done
done
typst compile --features html --root . --creation-timestamp 0 docs/architecture/generation-mark-geometry.typ target/mark-geometry-corrected/note-reproduced.pdf > target/mark-geometry-corrected/note-reproduction.log 2>&1 && cmp target/mark-geometry-corrected/generation-mark-geometry.pdf target/mark-geometry-corrected/note-reproduced.pdf
echo "$?" > target/mark-geometry-corrected/note-reproduction.exit
typst compile --features html,bundle --format bundle --root . --creation-timestamp 0 target/mark-geometry-corrected/generation-mark-geometry-wrapper.typ target/mark-geometry-corrected/note-html-reproduced >> target/mark-geometry-corrected/note-reproduction.log 2>&1 && diff -r target/mark-geometry-corrected/generation-mark-geometry-html target/mark-geometry-corrected/note-html-reproduced >> target/mark-geometry-corrected/note-reproduction.log 2>&1
echo "$?" > target/mark-geometry-corrected/note-html-reproduction.exit
