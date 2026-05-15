#import "../src/lib.typ": draw, graph, layout, subgraph

#let a = 0
#let c = 1
#let g = graph.build(
  name: "example",
  nodes: (graph.node(name: "a"), graph.node(name: "c")),
  edges: (
    graph.edge(
      source: (node: a, compass: "e"),
      sink: (node: c, compass: "w"),
      statements: (
        color: "0055ff",
        source-color: "d72638",
        sink-color: "1b7f4c",
        label: "a-c",
      ),
    ),
  ),
)
#let g = layout(g)
#let east = subgraph.compass(g, "e")
#let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.label]
#let source-style(edge) = (stroke: rgb("#" + edge.source-color) + 0.5pt)
#let sink-style(edge) = (stroke: rgb("#" + edge.sink-color) + 0.5pt)

= Build Example

#draw(g, subgraph: east, edge-label: edge-label, source-style: source-style, sink-style: sink-style)

#table(
  columns: (auto, 1fr),
  [nodes], [#graph.nodes(g).len()],
  [edges], [#graph.edges(g).len()],
  [east hedges], [#subgraph.hedges(east).join(", ")],
  [dot chars], [#graph.dot(g).len()],
)
