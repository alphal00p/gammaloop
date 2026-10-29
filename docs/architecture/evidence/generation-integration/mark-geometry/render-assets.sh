set -uo pipefail
source target/mark-geometry/environment.sh
just docs-svg-assets > target/mark-geometry/website-assets.log 2>&1
echo "$?" > target/mark-geometry/website-assets.exit
scripts/render-docs-svg-assets.sh target/mark-geometry/reproduced > target/mark-geometry/website-reproduction.log 2>&1
status=$?
if [ "$status" -eq 0 ]; then while IFS= read -r -d '' asset; do cmp "$asset" "${asset#target/mark-geometry/reproduced/}" || status=1; done < <(find target/mark-geometry/reproduced -name '*.svg' -print0); fi
echo "$status" > target/mark-geometry/website-reproduction.exit
for example in soft-triangle thrown-out xbox-sum; do
 typst compile --features html --root . --creation-timestamp 0 "crates/linnest/typst/examples/$example.typ" "crates/linnest/typst/examples/$example.pdf" > "target/mark-geometry/example-$example.log" 2>&1
 echo "$?" > "target/mark-geometry/example-$example.exit"
 typst compile --features html --root . --creation-timestamp 0 "crates/linnest/typst/examples/$example.typ" "target/mark-geometry/$example-reproduced.pdf" >> "target/mark-geometry/example-$example.log" 2>&1 && cmp "crates/linnest/typst/examples/$example.pdf" "target/mark-geometry/$example-reproduced.pdf"
 echo "$?" > "target/mark-geometry/example-$example-reproduction.exit"
done
typst compile --features html --root . --creation-timestamp 0 crates/kurvst/typst/docs/manual.typ crates/kurvst/typst/docs/manual.pdf > target/mark-geometry/manual.log 2>&1
echo "$?" > target/mark-geometry/manual.exit
typst compile --features html --root . --creation-timestamp 0 crates/kurvst/typst/docs/manual.typ target/mark-geometry/manual-reproduced.pdf >> target/mark-geometry/manual.log 2>&1 && cmp crates/kurvst/typst/docs/manual.pdf target/mark-geometry/manual-reproduced.pdf
echo "$?" > target/mark-geometry/manual-reproduction.exit
for theme in light dark; do
 for format in svg png pdf; do
  typst compile --root . --creation-timestamp 0 --input "theme=$theme" --format "$format" crates/kurvst/typst/examples/marks.typ "target/mark-geometry/marks-$theme.$format" > "target/mark-geometry/gallery-$theme-$format.log" 2>&1
  echo "$?" > "target/mark-geometry/gallery-$theme-$format.exit"
 done
 typst compile --root . --creation-timestamp 0 --input "theme=$theme" --format png docs/assets/typst/portal-graphs/graphs/portal-graph-epem-ttbar-cut.typ "target/mark-geometry/diagram-$theme.png" > "target/mark-geometry/diagram-$theme.log" 2>&1
 echo "$?" > "target/mark-geometry/diagram-$theme.exit"
done
