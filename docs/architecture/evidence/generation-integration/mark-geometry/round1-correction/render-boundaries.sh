set -uo pipefail
source target/mark-geometry-corrected/environment.sh
for theme in light dark; do
 for format in svg png pdf; do
  typst compile --root . --creation-timestamp 0 --input "theme=$theme" --format "$format" crates/kurvst/typst/examples/mark-boundaries.typ "target/mark-geometry-corrected/boundaries-$theme.$format" > "target/mark-geometry-corrected/boundaries-$theme-$format.log" 2>&1
  echo "$?" > "target/mark-geometry-corrected/boundaries-$theme-$format.exit"
  typst compile --root . --creation-timestamp 0 --input "theme=$theme" --format "$format" crates/kurvst/typst/examples/mark-boundaries.typ "target/mark-geometry-corrected/boundaries-$theme-reproduced.$format" >> "target/mark-geometry-corrected/boundaries-$theme-$format.log" 2>&1 && cmp "target/mark-geometry-corrected/boundaries-$theme.$format" "target/mark-geometry-corrected/boundaries-$theme-reproduced.$format"
  echo "$?" > "target/mark-geometry-corrected/boundaries-$theme-$format-reproduction.exit"
 done
done
