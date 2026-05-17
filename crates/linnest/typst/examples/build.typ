#import "../src/lib.typ": draw, edge, graph, layout, node, sink, source, subgraph

#let g = graph.build({
  node(<a>, label: "A")
  node(<c>, label: "C")
  edge(
    source(<a>, name: <a-c-source>, id: 0, compass: "e"),
    <a-c>,
    sink(<c>, name: <a-c-sink>, id: 1, compass: "w"),
    label: "a-c",
    statements: (
      color: "0055ff",
      source-color: "d72638",
      sink-color: "1b7f4c",
    ),
  )
},
  name: "example",
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
