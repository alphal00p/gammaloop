#import "../src/lib.typ" as kurvst
#import "@preview/cetz:0.3.4" as cetz

#set page(width: 120mm, height: 80mm, margin: 10mm)

#let segment = (
  start: (x: 0, y: 0),
  ctrl-a: (x: 1.1, y: 0.6),
  ctrl-b: (x: 2.2, y: -0.6),
  end: (x: 3.3, y: 0),
)

#let custom = (
  kind: "points",
  name: "hook",
  interpolation: "linear",
  points: (
    (at: 0, x: 0, y: 0),
    (at: 0.25, x: 0.15, y: 1),
    (at: 0.5, x: 0, y: 0),
    (at: 0.75, x: -0.15, y: -1),
    (at: 1, x: 0, y: 0),
  ),
)

#let paths = (
  kurvst.pattern-segment(segment, pattern: kurvst.wave(), amplitude: 0.14, wavelength: 0.7),
  kurvst.pattern-segment(segment, pattern: kurvst.zigzag(), amplitude: 0.14, wavelength: 0.7),
  kurvst.pattern-segment(segment, pattern: kurvst.coil(longitudinal-scale: 1.6), amplitude: 0.14, wavelength: 0.7),
  kurvst.pattern-segment(segment, pattern: custom, amplitude: 0.14, wavelength: 0.7),
)

#cetz.canvas({
  for (index, path) in paths.enumerate() {
    cetz.draw.line((0, -index * 0.65), (3.3, -index * 0.65), stroke: gray + 0.25pt)
    let shifted = (
      points: path.points.map(point => (x: point.x, y: point.y - index * 0.65)),
      curves: path.curves.map(curve => (
        start: (x: curve.start.x, y: curve.start.y - index * 0.65),
        ctrl-a: (x: curve.ctrl-a.x, y: curve.ctrl-a.y - index * 0.65),
        ctrl-b: (x: curve.ctrl-b.x, y: curve.ctrl-b.y - index * 0.65),
        end: (x: curve.end.x, y: curve.end.y - index * 0.65),
      )),
      segments: path.segments,
    )
    kurvst.cetz-pattern(shifted, stroke: black + 0.6pt)
  }
})
