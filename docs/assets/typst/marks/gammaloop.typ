#import "@preview/cetz:0.5.1"
#import "../theme.typ": palette

#set page(height: auto, width: auto, margin: 10mm, fill: none)

// GammaLoop mark, ported from the original Typst/CeTZ artwork. Keeping the
// geometry here makes the logo editable without touching generated SVG paths.
#cetz.canvas(length: 15mm, {
  import cetz.draw: *

  let weight = 25pt
  let line-stroke = (thickness: weight, paint: palette.ink)
  let square-stroke = (cap: "square", thickness: weight, paint: palette.ink)
  let crossing-mask = (thickness: 48pt, paint: palette.surface)
  let origin = (0, 0)

  line((rel: (0, 0), to: origin), (rel: (9, 0)), stroke: line-stroke)
  arc(
    (rel: (0.2, 0)),
    start: -90deg,
    stop: 90deg,
    stroke: square-stroke,
    radius: 0.75,
  )
  line((rel: (-0.2, 0)), (rel: (-1, 0)), stroke: square-stroke)
  line((), (rel: (0, -1)), stroke: line-stroke)
  line((rel: (0, -1)), (rel: (0, -1.3)), stroke: line-stroke)

  bezier(
    (rel: (-0.6, -3.5), to: origin),
    (rel: (-3.5, 0), to: origin),
    (rel: (-0.6, -2.8), to: origin),
    (rel: (-0.9, -0.9), to: origin),
    stroke: line-stroke,
  )
  bezier(
    origin,
    (rel: (-2.6, -3.5), to: origin),
    (rel: (-0.86, 0), to: origin),
    (rel: (-2.6, -1.7), to: origin),
    stroke: crossing-mask,
  )
  bezier(
    origin,
    (rel: (-2.6, -3.5), to: origin),
    (rel: (-0.86, 0), to: origin),
    (rel: (-2.6, -1.7), to: origin),
    stroke: (cap: "round", thickness: weight, paint: palette.ink),
  )
  arc((rel: (0, 0)), start: -180deg, stop: 0deg, stroke: line-stroke)

  line((rel: (0, 0.8), to: origin), (rel: (0, 0.7)), stroke: square-stroke)
  line(
    (rel: (0, -0.8), to: origin),
    (rel: (0, -0.7)),
    (rel: (1.4, 0)),
    stroke: square-stroke,
  )

  arc(
    (rel: (4, 0.5), to: origin),
    start: 0deg,
    stop: 180deg,
    stroke: line-stroke,
  )
  arc((rel: (0, -1)), start: -180deg, stop: 0deg, stroke: line-stroke)
  arc(
    (rel: (7, 0.5), to: origin),
    start: 0deg,
    stop: 180deg,
    stroke: line-stroke,
  )
  arc((rel: (0, -1)), start: -180deg, stop: 0deg, stroke: line-stroke)
})
