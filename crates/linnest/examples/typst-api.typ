#import "../typst/src/lib.typ": draw, graph, layout, subgraph
#import graph: build, dot, edge, edges, node, nodes, parse, sink, source

#let raw-graph = build({
  node(<a>, statements: (fill-color: "cfe8ff"))
  node(<b>, statements: (fill-color: "d6f5d6"))
  edge(
    source(<a>, compass: "e", statement: "out"),
    <a-to-b>,
    sink(<b>, compass: "w", statement: "in"),
    label: [a-to-b],
    statements: (
      color: "8a2be2",
      source-color: "d72638",
      sink-color: "1b7f4c",
    ),
  )
  edge(
    <incoming>,
    sink(<a>, compass: "n"),
    label: [incoming],
    statements: (
      color: "666666",
      source-color: "666666",
      sink-color: "355c9a",
    ),
  )
  edge(
    source(<b>, compass: "s"),
    <outgoing>,
    label: [outgoing],
    statements: (
      color: "666666",
      source-color: "8f5d2a",
      sink-color: "666666",
    ),
  )
},
  name: "constructed",
  statements: (full_num: "x + y"),
)

#let g = layout(raw-graph)
#let info = graph.info(g)
#let node-records = nodes(g)
#let edge-records = edges(g)
#let dot-text = dot(g)
#let north = subgraph.compass(g, "n")
#let internal-edge = subgraph.bits(g, (true, true, false, false))
#let north-label = subgraph.to-label(north)
#let north-hedges = subgraph.hedges(north)
#let internal-label = subgraph.to-label(internal-edge)
#let north-edges = edges(g, subgraph: north)
#let internal-nodes = nodes(g, subgraph: internal-edge)
#let internal-has-hedge-zero = subgraph.contains(internal-edge, 0)

= Linnest Typst API Example

This example imports the canonical Linnest Typst package source, builds a graph object from node and edge specs, lays it out, and queries it through subgraph objects.

#let node-style(node) = (fill: rgb("#" + node.fill-color))
#let edge-label(edge) = text(fill: rgb("#" + edge.color))[#edge.label]
#let source-style(edge) = (stroke: rgb("#" + edge.source-color) + 0.5pt)
#let sink-style(edge) = (stroke: rgb("#" + edge.sink-color) + 0.5pt)

#draw(g, node-style: node-style, edge-label: edge-label, source-style: source-style, sink-style: sink-style)

== Positioned Graph Without Layout

#let positioned = build({
  node(<left>, label: [left], pos: graph.pos(x: 0, y: 0))
  node(<right>, label: [right], pos: graph.pos(ref: <left>, dx: 2.5, dy: 0))
  edge(source(<left>), sink(<right>))
  edge(source(<right>), <right-out>, pos: graph.pos(ref: <right>, dx: 0.9, dy: 0.7))
},
  name: "positioned",
)

#draw(positioned, source-style: (stroke: black + 0.7pt), sink-style: (stroke: black + 0.7pt))

#table(
  columns: (auto, 1fr),
  inset: 6pt,
  stroke: 0.5pt,
[graph name], [#info.name],
[nodes], [#node-records.len()],
[edges], [#edge-records.len()],
[DOT characters], [#dot-text.len()],
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
  ..node-records.slice(0, calc.min(5, node-records.len())).map(node => (
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
  ..edge-records.slice(0, calc.min(5, edge-records.len())).map(edge => (
    [#edge.edge],
    [#edge.orientation],
    [#if edge.source == none { "none" } else { edge.source.node }],
    [#if edge.sink == none { "none" } else { edge.sink.node }],
  )).flatten(),
)
