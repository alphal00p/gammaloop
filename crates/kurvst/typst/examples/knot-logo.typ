#import "../src/lib.typ" as kurvst
#import "@preview/cetz:0.5.1" as cetz

#set page(width: auto, height: auto, margin: 10mm, fill: none)

// The GammaLoop logo (page 1) and its construction (page 2). The website marks
// in docs/assets/typst/marks import `logo` and `construction` from this file.

// Letter stroke width and the gap cut on each side of an over-crossing stroke,
// in canvas units.
#let width = 0.6
#let gap = 0.15

// Stadium-shaped "o" whose right edge starts at `x`.
#let o(x) = kurvst.path(
  kurvst.arc((x - 1, 0.5), 1, 0deg, 180deg),
  kurvst.line((x - 2, 0.5), (x - 2, -0.5)),
  kurvst.arc((x - 1, -0.5), 1, 180deg, 360deg),
  kurvst.line((x, -0.5), (x, 0.5)),
  kurvst.close(),
)

// The gamma rises out of its loop into the bar, which passes over everything.
#let over = kurvst.path(
  kurvst.cubic((-2.6, -3.5), (-2.6, -1.7), (-0.86, 0), (0, 0)),
  kurvst.line((0, 0), (9, 0)),
)

// The gamma, bar, and P form one continuous stroke.
#let gamma-p = kurvst.path(
  kurvst.cubic((-3.5, 0), (-0.9, -0.9), (-0.6, -2.8), (-0.6, -3.5)),
  kurvst.arc((-1.6, -3.5), 1, 0deg, -180deg),
  over,
  kurvst.arc((9, 0.75), 0.75, -90deg, 90deg),
  kurvst.line((9, 1.5), (7.8, 1.5)),
  kurvst.line((7.8, 1.5), (7.8, -1.6)),
)

#let l = kurvst.path(
  kurvst.line((0, 1.5), (0, -1.5)),
  kurvst.line((0, -1.5), (1.4, -1.5)),
)

// Draw an outline as one compound path with a closed subpath per contour, so
// CeTZ boolean operations see holes instead of bridged contours.
#let solid(outline, ..style) = {
  let contours = ()
  for element in kurvst.elements(outline) {
    if element.kind == "move" { contours.push(()) }
    contours.last().push(element)
  }
  cetz.draw.compound-path(
    contours.map(contour => kurvst.to-cetz(kurvst.from-elements(contour), close: true)).join(),
    ..style,
  )
}

#let union(..shapes) = shapes.pos().reduce((a, b) => cetz.draw.boolean(a, b, op: "union"))

// Each stroke's centerline and filled outline, in drawing order.
#let pieces = (
  (centerline: gamma-p, outline: kurvst.outline(gamma-p, width: width)),
  (centerline: l, outline: kurvst.outline(l, width: width, cap: "square")),
  ..(4, 7).map(x => (centerline: o(x), outline: kurvst.outline(o(x), width: width))),
)

#let letters = union(..pieces.map(piece => solid(piece.outline)))

// Strips hugging both sides of the over stroke; their butt ends leave the
// strokes it continues into attached.
#let strips(..style) = cetz.draw.boolean(
  solid(kurvst.outline(over, width: width + 2 * gap)),
  solid(kurvst.outline(over, width: width)),
  op: "difference",
  ..style,
)

// The logo: every letter outline, minus the strips around the over stroke.
#let logo(fill) = cetz.canvas(length: 15mm, cetz.draw.boolean(
  letters,
  strips(),
  op: "difference",
  fill: fill,
  stroke: none,
))

// Every operand of the boolean operations above. Letter outlines are
// translucent (by `tint`) so overlaps show, the cut strips are `cut`, and
// centerlines are dashed `line` strokes.
#let construction(
  colors: (blue, orange, green, purple),
  cut: red,
  line: black,
  tint: 70%,
) = cetz.canvas(length: 15mm, {
  for (piece, color) in pieces.zip(colors) {
    solid(piece.outline, fill: color.transparentize(tint), stroke: color + 0.6pt)
  }
  strips(fill: cut.transparentize(40%), stroke: cut + 0.6pt)
  for piece in pieces {
    kurvst.to-cetz(piece.centerline, stroke: (paint: line, thickness: 0.4pt, dash: "dashed"))
  }
})

#logo(black)
#pagebreak()
#construction()
