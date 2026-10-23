#import "../src/lib.typ" as kurvst
#import "@preview/cetz:0.5.1" as cetz

#set page(width: auto, height: auto, margin: 10mm)

// Stadium-shaped "o" whose right edge starts at `x`.
#let o(x) = kurvst.path(
  kurvst.arc((x - 1, 0.5), 1, 0deg, 180deg),
  kurvst.line((x - 2, 0.5), (x - 2, -0.5)),
  kurvst.arc((x - 1, -0.5), 1, 180deg, 360deg),
  kurvst.line((x, -0.5), (x, 0.5)),
  kurvst.close(),
)

// The gamma, top bar, and P bowl form one continuous stroke.
#let gamma-p = kurvst.path(
  kurvst.cubic((-3.5, 0), (-0.9, -0.9), (-0.6, -2.8), (-0.6, -3.5)),
  kurvst.arc((-1.6, -3.5), 1, 0deg, -180deg),
  kurvst.cubic((-2.6, -3.5), (-2.6, -1.7), (-0.86, 0), (0, 0)),
  kurvst.line((0, 0), (9, 0)),
  kurvst.arc((9, 0.75), 0.75, -90deg, 90deg),
  kurvst.line((9, 1.5), (7.8, 1.5)),
  kurvst.line((7.8, 1.5), (7.8, -1.6)),
)

#let l = kurvst.path(
  kurvst.line((0, 1.5), (0, -1.5)),
  kurvst.line((0, -1.5), (1.4, -1.5)),
)

// Stroke width in canvas units (the canvas uses 15mm per unit).
#let width = 25pt / 15mm

#let solid-logo(fill: black) = {
  kurvst.to-cetz(kurvst.outline(gamma-p, width: width), fill: none, stroke: fill)
  kurvst.to-cetz(kurvst.outline(l, width: width, cap: "square"), fill: none, stroke: fill)
  for x in (4, 7) {
    kurvst.to-cetz(kurvst.outline(o(x), width: width), fill: none, stroke: fill)
  }
}

#cetz.canvas(length: 15mm, solid-logo())

