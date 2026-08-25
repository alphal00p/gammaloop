#import "@preview/cetz:0.5.1"
#import "../theme.typ": palette

#set page(width: 3.7cm, height: 3.7cm, margin: 1mm, fill: none)

// Local Unitarity collaboration mark, ported from the original Typst artwork.
#let logo(origin, cap: "square", weight: 7pt) = {
  let line-stroke = (thickness: weight, paint: palette.ink)
  let integral-stroke = (cap: cap, thickness: weight, paint: palette.ink)
  let square-stroke = (cap: "square", thickness: weight, paint: palette.ink)
  let half-gap = float(weight / 1cm)
  let at(x, y) = (rel: (x, y), to: origin)

  // Leave real gaps where the integral crosses the three upright strokes.
  // Unlike a surface-coloured underlay, these remain transparent everywhere.
  let left-crossing = 0.9865245
  let middle-crossing = 1.0
  let right-crossing = 1.0134755

  cetz.draw.line(
    at(0.8, 0),
    at(0, 0),
    at(0, left-crossing - half-gap),
    stroke: line-stroke,
  )
  cetz.draw.line(
    at(0, left-crossing + half-gap),
    at(0, 1.5),
    stroke: line-stroke,
  )
  cetz.draw.arc(
    (rel: (0, 0)),
    start: 180deg,
    stop: 0deg,
    stroke: square-stroke,
    radius: 0.5,
  )
  cetz.draw.line(
    at(1, 1.5),
    at(1, middle-crossing + half-gap),
    stroke: line-stroke,
  )
  cetz.draw.line(
    at(1, middle-crossing - half-gap),
    at(1, 0.5),
    stroke: line-stroke,
  )
  cetz.draw.arc(
    (rel: (0, 0)),
    start: -180deg,
    stop: 0deg,
    stroke: square-stroke,
    radius: 0.5,
  )
  cetz.draw.line(
    at(2, 0.5),
    at(2, right-crossing - half-gap),
    stroke: line-stroke,
  )
  cetz.draw.line(
    at(2, right-crossing + half-gap),
    at(2, 2),
    at(1.2, 2),
    stroke: line-stroke,
  )

  cetz.draw.bezier(
    at(-0.5, 1.5),
    at(2.5, 0.5),
    at(-0.5, 0.25),
    at(2.5, 1.75),
    stroke: integral-stroke,
  )
}

#v(1fr)
#h(1fr)
#cetz.canvas(
  baseline: (0, 0),
  length: 1cm,
  { logo((0, 0)) },
)
#h(1fr)
#v(1fr)
