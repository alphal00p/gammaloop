#import "../src/lib.typ": draw, graph, layout
#import graph: build, dot, edge, edges, node, nodes, parse, sink, source

#set page(width: 120mm, height: 75mm, margin: 10mm)

#let g = build({
  node(<a>, pos: graph.pos(x: -2.0, y: 0, mode: "pin"))
  node(<c>, pos: graph.pos(ref: <a>, dx: 4.0, dy: 0, mode: "pin"))
  edge(<e0>, source(<a>), sink(<c>), pos: graph.pos(x: 0, y: 1.25, mode: "pin"))
  edge(<e1>, source(<a>), sink(<c>), pos: graph.pos(x: 0, y: 1.25, mode: "pin"))
  edge(<e2>, source(<a>), sink(<c>), pos: graph.pos(x: 0, y: 1.25, mode: "pin"))
},
  name: "parallel demo",
)

#let offsets = (-0.16, 0, 0.16)
#let colors = (rgb("#d72638"), gray, rgb("#355c9a"))
#let parallel-limits(edge) = if edge.eid == 1 {
  (:)
} else {
  (length: 1.6, ratio: 0.5)
}
#let source-style(edge) = (
  stroke: (paint: colors.at(edge.eid), thickness: 0.75pt, cap: "round"),
  offset: offsets.at(edge.eid),
) + parallel-limits(edge)
#let sink-style(edge) = (
  stroke: (paint: colors.at(edge.eid), thickness: 0.75pt, cap: "round"),
  offset: offsets.at(edge.eid),
) + parallel-limits(edge)

#draw(
  layout(g, steps: 1, epochs: 1, label-steps: 0),
  unit: 1.35,
  node-radius: 0.25,
  source-style: source-style,
  sink-style: sink-style,
)
