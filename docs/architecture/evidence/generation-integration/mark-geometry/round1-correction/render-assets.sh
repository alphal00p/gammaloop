set -uo pipefail
source target/mark-geometry-corrected/environment.sh
just docs-svg-assets > target/mark-geometry-corrected/website-assets.log 2>&1
echo "$?" > target/mark-geometry-corrected/website-assets.exit
scripts/render-docs-svg-assets.sh target/mark-geometry-corrected/reproduced > target/mark-geometry-corrected/website-reproduction.log 2>&1
status=$?
if [ "$status" -eq 0 ]; then while IFS= read -r -d '' asset; do cmp "$asset" "${asset#target/mark-geometry-corrected/reproduced/}" || status=1; done < <(find target/mark-geometry-corrected/reproduced -name '*.svg' -print0); fi
echo "$status" > target/mark-geometry-corrected/website-reproduction.exit
for example in soft-triangle thrown-out xbox-sum; do
 typst compile --features html --root . --creation-timestamp 0 "crates/linnest/typst/examples/$example.typ" "crates/linnest/typst/examples/$example.pdf" > "target/mark-geometry-corrected/example-$example.log" 2>&1
 echo "$?" > "target/mark-geometry-corrected/example-$example.exit"
 typst compile --features html --root . --creation-timestamp 0 "crates/linnest/typst/examples/$example.typ" "target/mark-geometry-corrected/$example-reproduced.pdf" >> "target/mark-geometry-corrected/example-$example.log" 2>&1 && cmp "crates/linnest/typst/examples/$example.pdf" "target/mark-geometry-corrected/$example-reproduced.pdf"
 echo "$?" > "target/mark-geometry-corrected/example-$example-reproduction.exit"
done
typst compile --features html --root . --creation-timestamp 0 crates/kurvst/typst/docs/manual.typ crates/kurvst/typst/docs/manual.pdf > target/mark-geometry-corrected/manual.log 2>&1
echo "$?" > target/mark-geometry-corrected/manual.exit
typst compile --features html --root . --creation-timestamp 0 crates/kurvst/typst/docs/manual.typ target/mark-geometry-corrected/manual-reproduced.pdf >> target/mark-geometry-corrected/manual.log 2>&1 && cmp crates/kurvst/typst/docs/manual.pdf target/mark-geometry-corrected/manual-reproduced.pdf
echo "$?" > target/mark-geometry-corrected/manual-reproduction.exit
for theme in light dark; do
 for format in svg png pdf; do
  typst compile --root . --creation-timestamp 0 --input "theme=$theme" --format "$format" crates/kurvst/typst/examples/marks.typ "target/mark-geometry-corrected/marks-$theme.$format" > "target/mark-geometry-corrected/gallery-$theme-$format.log" 2>&1
  echo "$?" > "target/mark-geometry-corrected/gallery-$theme-$format.exit"
 done
 typst compile --root . --creation-timestamp 0 --input "theme=$theme" --format png docs/assets/typst/portal-graphs/graphs/portal-graph-epem-ttbar-cut.typ "target/mark-geometry-corrected/diagram-$theme.png" > "target/mark-geometry-corrected/diagram-$theme.log" 2>&1
 echo "$?" > "target/mark-geometry-corrected/diagram-$theme.exit"
done
