#import "../../clinnet/templates/linnest.typ": (
  builder,
  builder-add-edge,
  builder-add-node,
  default-layout-config,
  finish-builder,
  graph-edges,
  graph-info,
  graph-nodes,
  layout-graph,
  subgraph-contains-hedge,
  subgraph-from-bits,
  subgraph-from-compass,
  subgraph-hedges,
  subgraph-label,
)

#let builder = builder(spec: (
  name: "constructed",
  statements: (full_num: "x + y"),
))
#let a = builder-add-node(builder, name: "a", statements: (eval: "(fill: blue.lighten(70%))"))
#let builder = a.builder
#let b = builder-add-node(builder, name: "b", statements: (eval: "(fill: green.lighten(70%))"))
#let builder = b.builder
#let builder = builder-add-edge(
  builder,
  source: (node: a.node, compass: "e", statement: "out"),
  sink: (node: b.node, compass: "w", statement: "in"),
  statements: (label: "a-to-b"),
)
#let builder = builder-add-edge(
  builder,
  source: none,
  sink: (node: a.node, compass: "n"),
  statements: (label: "incoming"),
)
#let builder = builder-add-edge(
  builder,
  source: (node: b.node, compass: "s"),
  sink: none,
  statements: (label: "outgoing"),
)
#let raw-graph = finish-builder(builder)

#let graph = layout-graph(raw-graph, config: default-layout-config + (seed: "2", steps: "5"))
#let info = graph-info(graph)
#let nodes = graph-nodes(graph)
#let edges = graph-edges(graph)
#let north = subgraph-from-compass(graph, "n")
#let internal-edge = subgraph-from-bits(graph, (true, true, false, false))
#let north-label = subgraph-label(north)
#let north-hedges = subgraph-hedges(north)
#let internal-label = subgraph-label(internal-edge)
#let north-edges = graph-edges(graph, subgraph: north)
#let internal-nodes = graph-nodes(graph, subgraph: internal-edge)
#let internal-has-hedge-zero = subgraph-contains-hedge(internal-edge, 0)

= Linnest Typst API Example

This example imports the local `linnest.typ` wrapper, builds an archived graph through the archived builder API, lays it out, and queries it through archived subgraph objects.

#table(
  columns: (auto, 1fr),
  inset: 6pt,
  stroke: 0.5pt,
  [graph name], [#info.name],
  [nodes], [#nodes.len()],
  [edges], [#edges.len()],
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
