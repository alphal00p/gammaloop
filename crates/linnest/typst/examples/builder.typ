#import "../src/lib.typ": graph, layout, subgraph

#let b = graph.builder(name: "example")
#let (node: a, builder: b) = graph.node(b, name: "a")
#let (node: c, builder: b) = graph.node(b, name: "c")
#let b = graph.edge(b, source: (node: a, compass: "e"), sink: (node: c, compass: "w"))
#let g = layout(graph.finish(b), seed: "2", steps: "5")
#let east = subgraph.compass(g, "e")

= Builder Example

#table(
  columns: (auto, 1fr),
  [nodes], [#graph.nodes(g).len()],
  [edges], [#graph.edges(g).len()],
  [east hedges], [#subgraph.hedges(east).join(", ")],
  [dot chars], [#graph.dot(g).len()],
)
