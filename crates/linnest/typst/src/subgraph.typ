#let _plugin = plugin("../linnest.wasm")

/// Construct an subgraph object from a base62 label.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>, compass: "e"), graph.sink(<b>))
/// })
/// #let east = subgraph.compass(g, "e")
/// #let same = subgraph.label(g, subgraph.to-label(east))
/// #subgraph.hedges(same).len()
/// ```
/// -> bytes
#let label(
  /// Graph object whose half-edge set the label refers to. -> bytes
  graph,
  /// Base62 subgraph label returned by @to-label or produced by Linnest. -> string
  label,
) = _plugin.graph_archived_subgraph(bytes(graph), cbor.encode(label))

/// Construct an subgraph object from a boolean hedge array.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>), graph.sink(<b>))
/// })
/// #let first = subgraph.bits(g, (true, false))
/// #subgraph.contains(first, 0)
/// ```
/// -> bytes
#let bits(
  /// Graph object whose half-edge order defines the bit array. -> bytes
  graph,
  /// Boolean array selecting half edges by graph half-edge index. -> array
  bits,
) = _plugin.graph_archived_subgraph(bytes(graph), cbor.encode(bits))

/// Construct an subgraph object from a DOT compass point such as `"n"` or
/// `"s"`.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>, compass: "e"), graph.sink(<b>))
/// })
/// #subgraph.hedges(subgraph.compass(g, "e")).len()
/// ```
/// -> bytes
#let compass(
  /// Graph object to filter by half-edge compass statement. -> bytes
  graph,
  /// DOT compass point such as `"n"`, `"s"`, `"e"`, or `"w"`. -> string
  compass,
) = {
  _plugin.graph_archived_compass_subgraph(bytes(graph), cbor.encode(compass))
}

/// Convert an subgraph object to its base62 label.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>), graph.sink(<b>))
/// })
/// #subgraph.to-label(subgraph.bits(g, (true, false)))
/// ```
/// -> string
#let to-label(
  /// Subgraph object returned by this module or by `graph.cycles`/`graph.forests`. -> bytes
  subgraph,
) = cbor(_plugin.subgraph_label(bytes(subgraph)))

/// Return the hedge indices included in an subgraph object.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>), graph.sink(<b>))
/// })
/// #subgraph.hedges(subgraph.bits(g, (true, false)))
/// ```
/// -> array
#let hedges(
  /// Subgraph object to inspect. -> bytes
  subgraph,
) = cbor(_plugin.subgraph_hedges(bytes(subgraph)))

/// Test whether an subgraph object includes a hedge.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>), graph.sink(<b>))
/// })
/// #subgraph.contains(subgraph.bits(g, (true, false)), 0)
/// ```
/// -> bool
#let contains(
  /// Subgraph object to inspect. -> bytes
  subgraph,
  /// Half-edge index to test. -> int
  hedge,
) = {
  cbor(_plugin.subgraph_contains_hedge(bytes(subgraph), cbor.encode(hedge)))
}
