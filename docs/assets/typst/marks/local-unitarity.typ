#import "@preview/cetz:0.5.1"
#import "../theme.typ": palette

#set page(width: 3.7cm, height: 3.7cm, margin: 1mm, fill: palette.surface)

// Local Unitarity collaboration mark, ported from the original Typst artwork.
#let logo(origin, cap: "square", weight: 7pt) = {
  let line-stroke = (thickness: weight, paint: palette.ink)
  let integral-stroke = (cap: cap, thickness: weight, paint: palette.ink)
  let square-stroke = (cap: "square", thickness: weight, paint: palette.ink)
  let crossing-mask = (thickness: weight * 2, paint: palette.surface)

  cetz.draw.line(
    (rel: (0.8, 0), to: origin),
    (rel: (-0.8, 0)),
    (rel: (0, 1.5)),
    stroke: line-stroke,
  )
  cetz.draw.arc(
    (rel: (0, 0)),
    start: 180deg,
    stop: 0deg,
    stroke: square-stroke,
    radius: 0.5,
  )
  cetz.draw.line((rel: (0, 0)), (rel: (0, -1)), stroke: line-stroke)
  cetz.draw.arc(
    (rel: (0, 0)),
    start: -180deg,
    stop: 0deg,
    stroke: square-stroke,
    radius: 0.5,
  )
  cetz.draw.line(
    (rel: (0, 0)),
    (rel: (0, 1.5)),
    (rel: (-0.8, 0)),
    stroke: line-stroke,
  )

  // Mask the crossing before drawing the integral curve on top.
  cetz.draw.bezier(
    (rel: (-0.5, 1.5), to: origin),
    (rel: (2.5, 0.5), to: origin),
    (rel: (-0.5, 0.25), to: origin),
    (rel: (2.5, 1.75), to: origin),
    stroke: crossing-mask,
  )
  cetz.draw.bezier(
    (rel: (-0.5, 1.5), to: origin),
    (rel: (2.5, 0.5), to: origin),
    (rel: (-0.5, 0.25), to: origin),
    (rel: (2.5, 1.75), to: origin),
    stroke: integral-stroke,
  )
}

#v(1fr)
#h(1fr)
#cetz.canvas(
  baseline: (0, 0),
  background: palette.surface,
  length: 1cm,
  { logo((0, 0)) },
)
#h(1fr)
#v(1fr)
