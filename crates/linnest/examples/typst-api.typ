#import "../../clinnet/templates/linnest.typ": (
  default-layout-config,
  graph-compass-subgraph,
  graph-edges,
  graph-info,
  graph-nodes,
  layout-graph,
  parse-graphs,
)

#let dot = read("gl303.dot")
#let parsed = parse-graphs(dot)
#let graph = layout-graph(parsed.first(), config: default-layout-config + (seed: "2", steps: "5"))
#let info = graph-info(graph)
#let nodes = graph-nodes(graph)
#let edges = graph-edges(graph)
#let north = graph-compass-subgraph(graph, "n")
#let north-edges = graph-edges(graph, subgraph: north)

= Linnest Typst API Example

This example imports the local `linnest.typ` wrapper, reads `gl303.dot`, lays out the first graph, and queries the archived graph bytes through the WASM plugin.

#table(
  columns: (auto, 1fr),
  inset: 6pt,
  stroke: 0.5pt,
  [graph name], [#info.name],
  [graphs parsed], [#parsed.len()],
  [nodes], [#nodes.len()],
  [edges], [#edges.len()],
  [north subgraph], [#north],
  [north edges], [#north-edges.len()],
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
