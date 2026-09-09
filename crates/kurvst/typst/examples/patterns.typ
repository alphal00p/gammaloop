#import "../src/lib.typ" as kurvst
#import "@preview/cetz:0.5.1" as cetz

#set page(width: 120mm, height: 80mm, margin: 10mm)

#let segment = (
  start: (0, 0),
  control-start: (1.1, 0.6),
  control-end: (2.2, -0.6),
  end: (3.3, 0),
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

#let path = kurvst.from-cubic(segment)
#let paths = (
  kurvst.pattern(path, pattern: kurvst.wave(), amplitude: 0.14, wavelength: 0.7),
  kurvst.pattern(path, pattern: kurvst.zigzag(), amplitude: 0.14, wavelength: 0.7),
  kurvst.pattern(path, pattern: kurvst.coil(longitudinal-scale: 1.6), amplitude: 0.14, wavelength: 0.7),
  kurvst.pattern(path, pattern: custom, amplitude: 0.14, wavelength: 0.7),
)

#cetz.canvas({
  for (index, path) in paths.enumerate() {
    cetz.draw.line((0, -index * 0.65), (3.3, -index * 0.65), stroke: gray + 0.25pt)
    cetz.draw.group({
      cetz.draw.translate(y: -index * 0.65)
      kurvst.to-cetz(path, stroke: black + 0.6pt)
    })
  }
})
