set -uo pipefail
source target/mark-geometry/environment.sh
cargo nextest list --locked -p kurvst -p linnest --features linnest/svg --cargo-profile dev-optim --list-type binaries-only --message-format json > target/mark-geometry/native-test-binaries.json 2> target/mark-geometry/native-test-binaries.log
echo "$?" > target/mark-geometry/native-test-binaries.exit
for theme in light dark; do
 for format in svg png pdf; do
  typst compile --root . --creation-timestamp 0 --input "theme=$theme" --format "$format" crates/kurvst/typst/examples/marks.typ "target/mark-geometry/marks-$theme-reproduced.$format" > "target/mark-geometry/gallery-$theme-$format-reproduction.log" 2>&1 && cmp "target/mark-geometry/marks-$theme.$format" "target/mark-geometry/marks-$theme-reproduced.$format"
  echo "$?" > "target/mark-geometry/gallery-$theme-$format-reproduction.exit"
 done
done
typst compile --features html --root . --creation-timestamp 0 docs/architecture/generation-mark-geometry.typ target/mark-geometry/note-reproduced.pdf > target/mark-geometry/note-reproduction.log 2>&1 && cmp target/mark-geometry/generation-mark-geometry.pdf target/mark-geometry/note-reproduced.pdf
echo "$?" > target/mark-geometry/note-reproduction.exit
typst compile --features html,bundle --format bundle --root . --creation-timestamp 0 docs/architecture/generation-mark-geometry.typ target/mark-geometry/note-html-reproduced >> target/mark-geometry/note-reproduction.log 2>&1 && diff -r target/mark-geometry/generation-mark-geometry-html target/mark-geometry/note-html-reproduced >> target/mark-geometry/note-reproduction.log 2>&1
echo "$?" > target/mark-geometry/note-html-reproduction.exit
