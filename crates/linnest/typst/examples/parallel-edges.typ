#import "../src/lib.typ": draw, graph, layout

#set page(width: 120mm, height: 75mm, margin: 10mm)

#let a = 0
#let c = 1
#let g = graph.build(
  name: "parallel demo",
  nodes: (
    graph.node(name: "a", pos: graph.pos(x: -2.0, y: 0, mode: "pin")),
    graph.node(name: "c", pos: graph.pos(ref: a, dx: 4.0, dy: 0, mode: "pin")),
  ),
  edges: (
    graph.edge(source: (node: a), sink: (node: c), pos: graph.pos(x: 0, y: 1.25, mode: "pin")),
    graph.edge(source: (node: a), sink: (node: c), pos: graph.pos(x: 0, y: 1.25, mode: "pin")),
    graph.edge(source: (node: a), sink: (node: c), pos: graph.pos(x: 0, y: 1.25, mode: "pin")),
  ),
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
