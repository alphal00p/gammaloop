#let _plugin = plugin("../linnest.wasm")

/// Construct an archived subgraph from a base62 label.
///
/// ```example
/// #let g = graph.build(
///   nodes: ((name: "a"), (name: "b")),
///   edges: ((source: (node: 0, compass: "e"), sink: (node: 1)),),
/// )
/// #let east = subgraph.compass(g, "e")
/// #let same = subgraph.label(g, subgraph.to-label(east))
/// #subgraph.hedges(same).len()
/// ```
/// -> bytes
#let label(graph, label) = _plugin.graph_archived_subgraph(bytes(graph), cbor.encode(label))

/// Construct an archived subgraph from a boolean hedge array.
///
/// ```example
/// #let g = graph.build(
///   nodes: ((name: "a"), (name: "b")),
///   edges: ((source: (node: 0), sink: (node: 1)),),
/// )
/// #let first = subgraph.bits(g, (true, false))
/// #subgraph.contains(first, 0)
/// ```
/// -> bytes
#let bits(graph, bits) = _plugin.graph_archived_subgraph(bytes(graph), cbor.encode(bits))

/// Construct an archived subgraph from a DOT compass point such as `"n"` or
/// `"s"`.
///
/// ```example
/// #let g = graph.build(
///   nodes: ((name: "a"), (name: "b")),
///   edges: ((source: (node: 0, compass: "e"), sink: (node: 1)),),
/// )
/// #subgraph.hedges(subgraph.compass(g, "e")).len()
/// ```
/// -> bytes
#let compass(graph, compass) = {
  _plugin.graph_archived_compass_subgraph(bytes(graph), cbor.encode(compass))
}

/// Convert an archived subgraph to its base62 label.
///
/// ```example
/// #let g = graph.build(
///   nodes: ((name: "a"), (name: "b")),
///   edges: ((source: (node: 0), sink: (node: 1)),),
/// )
/// #subgraph.to-label(subgraph.bits(g, (true, false)))
/// ```
/// -> string
#let to-label(subgraph) = cbor(_plugin.subgraph_label(bytes(subgraph)))

/// Return the hedge indices included in an archived subgraph.
///
/// ```example
/// #let g = graph.build(
///   nodes: ((name: "a"), (name: "b")),
///   edges: ((source: (node: 0), sink: (node: 1)),),
/// )
/// #subgraph.hedges(subgraph.bits(g, (true, false)))
/// ```
/// -> array
#let hedges(subgraph) = cbor(_plugin.subgraph_hedges(bytes(subgraph)))

/// Test whether an archived subgraph includes a hedge.
///
/// ```example
/// #let g = graph.build(
///   nodes: ((name: "a"), (name: "b")),
///   edges: ((source: (node: 0), sink: (node: 1)),),
/// )
/// #subgraph.contains(subgraph.bits(g, (true, false)), 0)
/// ```
/// -> bool
#let contains(subgraph, hedge) = {
  cbor(_plugin.subgraph_contains_hedge(bytes(subgraph), cbor.encode(hedge)))
}
