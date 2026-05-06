#import "../src/lib.typ": draw, graph, layout

#set page(width: 120mm, height: 75mm, margin: 10mm)

#let b = graph.builder(name: "parallel demo")
#let (node: a, builder: b) = graph.node(b, name: "a", pin: graph.pin(x: -2.0, y: 0))
#let (node: c, builder: b) = graph.node(b, name: "c", pin: graph.pin(x: 2.0, y: 0))
#let b = graph.edge(b, source: (node: a), sink: (node: c), pin: graph.pin(x: 0, y: 1.25))
#let b = graph.edge(b, source: (node: a), sink: (node: c), pin: graph.pin(x: 0, y: 1.25))
#let b = graph.edge(b, source: (node: a), sink: (node: c), pin: graph.pin(x: 0, y: 1.25))

#let offsets = (-0.16, 0, 0.16)
#let colors = (rgb("#d72638"), gray, rgb("#355c9a"))
#let parallel-limits(edge) = if edge.eid == 1 {
  (:)
} else {
  (parallel-length: 1.6, parallel-ratio: 0.5)
}
#let source-style(edge) = (
  stroke: (paint: colors.at(edge.eid), thickness: 0.75pt, cap: "round"),
  parallel-offset: offsets.at(edge.eid),
) + parallel-limits(edge)
#let sink-style(edge) = (
  stroke: (paint: colors.at(edge.eid), thickness: 0.75pt, cap: "round"),
  parallel-offset: offsets.at(edge.eid),
) + parallel-limits(edge)

#draw(
  layout(graph.finish(b), steps: 1, epochs: 1, label-steps: 0),
  unit: 1.35,
  node-radius: 0.25,
  source-style: source-style,
  sink-style: sink-style,
)
