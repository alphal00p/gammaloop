set -u
bash target/mark-geometry-round3/render-assets.sh
bash target/mark-geometry-round3/render-boundaries.sh
source target/mark-geometry-round3/environment.sh
out=target/mark-geometry-round3
for theme in light dark; do
 for format in svg png pdf; do
 typst compile --root . --creation-timestamp 0 --input "theme=$theme" --format "$format" crates/kurvst/typst/examples/marks.typ "$out/marks-$theme-reproduced.$format" > "$out/gallery-$theme-$format-reproduction.log" 2>&1 && cmp "$out/marks-$theme.$format" "$out/marks-$theme-reproduced.$format"
 echo $? > "$out/gallery-$theme-$format-reproduction.exit"
 done
done
