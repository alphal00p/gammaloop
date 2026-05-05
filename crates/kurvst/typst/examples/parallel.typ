#import "../src/lib.typ" as kurvst
#import "@preview/cetz:0.3.4" as cetz

#set page(width: 120mm, height: 70mm, margin: 10mm)

#let segment = (
  start: (x: 0, y: 0),
  ctrl-a: (x: 1.0, y: 0.7),
  ctrl-b: (x: 2.2, y: -0.7),
  end: (x: 3.3, y: 0),
)

#let left = kurvst.parallel-segment(segment, distance: 0.18)
#let right = kurvst.parallel-segment(segment, distance: -0.18)

#cetz.canvas({
  kurvst.cetz-bezier(segment, stroke: gray + 0.35pt)
  kurvst.cetz-path(left, stroke: red + 0.7pt)
  kurvst.cetz-path(right, stroke: blue + 0.7pt)
})
