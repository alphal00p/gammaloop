#import "../src/lib.typ" as kurvst
#import "@preview/cetz:0.5.1" as cetz

#set page(width: 120mm, height: 70mm, margin: 10mm)

#let segment = (
  start: (0, 0),
  control-start: (1.0, 0.7),
  control-end: (2.2, -0.7),
  end: (3.3, 0),
)

#let path = kurvst.from-cubic(segment)
#let left = kurvst.parallel(path, distance: 0.18)
#let right = kurvst.parallel(path, distance: -0.18)

#cetz.canvas({
  kurvst.to-cetz(path, stroke: gray + 0.35pt)
  kurvst.to-cetz(left, stroke: red + 0.7pt)
  kurvst.to-cetz(right, stroke: blue + 0.7pt)
})
