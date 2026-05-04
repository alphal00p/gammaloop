#import "../src/lib.typ": draw, graph, layout, subgraph

#let b = graph.builder(
  name: "example",
  edge-statements: (
    eval_label: "(text(fill: rgb(\"#{color}\"))[{label}])",
    eval_source: "(stroke: rgb(\"#{source-color}\") + 0.5pt)",
    eval_sink: "(stroke: rgb(\"#{sink-color}\") + 0.5pt)",
  ),
)
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

= Builder Example

#draw(g)

#table(
  columns: (auto, 1fr),
  [nodes], [#graph.nodes(g).len()],
  [edges], [#graph.edges(g).len()],
  [east hedges], [#subgraph.hedges(east).join(", ")],
  [dot chars], [#graph.dot(g).len()],
)
