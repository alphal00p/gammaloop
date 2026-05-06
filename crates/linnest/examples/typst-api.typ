#import "../../clinnet/templates/linnest.typ": draw, graph, layout, subgraph

#let b = graph.builder(
  name: "constructed",
  statements: (full_num: "x + y"),
)
#let (node: a, builder: b) = graph.node(b, name: "a", statements: (fill-color: "cfe8ff"))
#let (node: c, builder: b) = graph.node(b, name: "b", statements: (fill-color: "d6f5d6"))
#let b = graph.edge(
  b,
  source: (node: a, compass: "e", statement: "out"),
  sink: (node: c, compass: "w", statement: "in"),
  statements: (
    color: "8a2be2",
    source-color: "d72638",
    sink-color: "1b7f4c",
    label: "a-to-b",
  ),
)
#let b = graph.edge(
  b,
  source: none,
  sink: (node: a, compass: "n"),
  statements: (
    color: "666666",
    source-color: "666666",
    sink-color: "355c9a",
    label: "incoming",
  ),
)
#let b = graph.edge(
  b,
  source: (node: c, compass: "s"),
  sink: none,
  statements: (
    color: "666666",
    source-color: "8f5d2a",
    sink-color: "666666",
    label: "outgoing",
  ),
)
#let raw-graph = graph.finish(b)

#let g = layout(raw-graph)
#let info = graph.info(g)
#let nodes = graph.nodes(g)
#let edges = graph.edges(g)
#let dot = graph.dot(g)
#let north = subgraph.compass(g, "n")
#let internal-edge = subgraph.bits(g, (true, true, false, false))
#let north-label = subgraph.to-label(north)
#let north-hedges = subgraph.hedges(north)
#let internal-label = subgraph.to-label(internal-edge)
#let north-edges = graph.edges(g, subgraph: north)
#let internal-nodes = graph.nodes(g, subgraph: internal-edge)
#let internal-has-hedge-zero = subgraph.contains(internal-edge, 0)

= Linnest Typst API Example

This example imports the local `linnest.typ` wrapper, builds a graph object through the builder object API, lays it out, and queries it through subgraph objects.

#let node-style(node) = (fill: rgb("#" + node.fill-color))
#let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.label]
#let source-style(edge) = (stroke: rgb("#" + edge.source-color) + 0.5pt)
#let sink-style(edge) = (stroke: rgb("#" + edge.sink-color) + 0.5pt)

#draw(g, node-style: node-style, edge-label: edge-label, source-style: source-style, sink-style: sink-style)

#table(
  columns: (auto, 1fr),
  inset: 6pt,
  stroke: 0.5pt,
[graph name], [#info.name],
[nodes], [#nodes.len()],
[edges], [#edges.len()],
[DOT characters], [#dot.len()],
[north subgraph], [#north-label],
  [north hedges], [#north-hedges.join(", ")],
  [north edges], [#north-edges.len()],
  [internal edge subgraph], [#internal-label],
  [internal edge nodes], [#internal-nodes.len()],
  [internal has hedge 0], [#internal-has-hedge-zero],
)

== First Nodes

#table(
  columns: (auto, auto, auto, auto),
  inset: 4pt,
  stroke: 0.5pt,
  [node], [name], [x], [y],
  ..nodes.slice(0, calc.min(5, nodes.len())).map(node => (
    [#node.node],
    [#node.name],
    [#node.pos.x],
    [#node.pos.y],
  )).flatten(),
)

== First Edges

#table(
  columns: (auto, auto, auto, auto),
  inset: 4pt,
  stroke: 0.5pt,
  [edge], [orientation], [source], [sink],
  ..edges.slice(0, calc.min(5, edges.len())).map(edge => (
    [#edge.edge],
    [#edge.orientation],
    [#if edge.source == none { "none" } else { edge.source.node }],
    [#if edge.sink == none { "none" } else { edge.sink.node }],
  )).flatten(),
)
