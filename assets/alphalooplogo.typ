// Compile from the repository root:
// typst compile --root . assets/alphalooplogo.typ assets/alphalooplogo-light.svg
// typst compile --root . --input theme=dark assets/alphalooplogo.typ assets/alphalooplogo-dark.svg
#import "../crates/kurvst/typst/src/lib.typ" as kurvst
#import "@preview/cetz:0.5.2" as cetz
#import "../crates/kurvst/typst/examples/knot-logo.typ": (
  gap, l, o, solid, union, width,
)
#import "../docs/assets/typst/theme.typ": palette

#set document(title: "AlphaLoop logo")
#set page(width: auto, height: auto, margin: 10mm, fill: none)

// The compact alpha shares LOOP's baseline and sits close to the L.
// Its upper arm flows into the woven bar; its lower tail stays free.
#let over = kurvst.path(
  kurvst.cubic((-2.1, -1.5), (-1.1, -1.5), (-1.3, 0), (-0.5, 0)),
  kurvst.line((-0.5, 0), (9, 0)),
)

// One uninterrupted stroke makes α, the bar, and the P bowl.
#let alpha-p = kurvst.path(
  kurvst.cubic((-0.6, -1.5), (-1.35, -1.5), (-0.95, 0.5), (-2.1, 0.5)),
  kurvst.arc((-2.1, -0.5), 1, 90deg, 270deg),
  over,
  kurvst.arc((9, 0.75), 0.75, -90deg, 90deg),
  kurvst.line((9, 1.5), (7.8, 1.5)),
  kurvst.line((7.8, 1.5), (7.8, -1.6)),
)

#let letters = union(
  solid(kurvst.outline(alpha-p, width: width)),
  solid(kurvst.outline(l, width: width, cap: "square")),
  ..(4, 7).map(x => solid(kurvst.outline(o(x), width: width))),
)

// Cut real holes around the overpass, rather than painting a background
// colour: the exported mark stays transparent on any surface.
#let strips = cetz.draw.boolean(
  solid(kurvst.outline(over, width: width + 2 * gap)),
  solid(kurvst.outline(over, width: width)),
  op: "difference",
)

#cetz.canvas(length: 15mm, cetz.draw.boolean(
  letters,
  strips,
  op: "difference",
  fill: palette.ink,
  stroke: none,
))
