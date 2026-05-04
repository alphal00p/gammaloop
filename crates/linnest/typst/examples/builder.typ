#import "../src/lib.typ": draw, graph, layout, subgraph

#let b = graph.builder(name: "example")
#let (node: a, builder: b) = graph.node(b, name: "a")
#let (node: c, builder: b) = graph.node(b, name: "c")
#let b = graph.edge(
  b,
  source: (node: a, compass: "e"),
  sink: (node: c, compass: "w"),
  statements: (
    color: "0055ff",
    source-color: "d72638",
    sink-color: "1b7f4c",
    label: "a-c",
  ),
)
#let g = layout(graph.finish(b))
#let east = subgraph.compass(g, "e")
#let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.label]
#let source-style(edge) = (stroke: rgb("#" + edge.source-color) + 0.5pt)
#let sink-style(edge) = (stroke: rgb("#" + edge.sink-color) + 0.5pt)

= Builder Example

#draw(g, subgraph: east, edge-label: edge-label, source-style: source-style, sink-style: sink-style)

#table(
  columns: (auto, 1fr),
  [nodes], [#graph.nodes(g).len()],
  [edges], [#graph.edges(g).len()],
  [east hedges], [#subgraph.hedges(east).join(", ")],
  [dot chars], [#graph.dot(g).len()],
)
