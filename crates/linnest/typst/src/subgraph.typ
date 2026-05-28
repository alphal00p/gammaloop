#import "graph.typ" as graph-module

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
  /// Graph object whose half-edge set the label refers to. -> dictionary
  graph,
  /// Base62 subgraph label returned by @to-label or produced by Linnest. -> string
  label,
) = _plugin.graph_archived_subgraph(graph-module.graph-bytes(graph), cbor.encode(label))

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
  /// Graph object whose half-edge order defines the bit array. -> dictionary
  graph,
  /// Boolean array selecting half edges by graph half-edge index. -> array
  bits,
) = _plugin.graph_archived_subgraph(graph-module.graph-bytes(graph), cbor.encode(bits))

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
  /// Graph object to filter by half-edge compass statement. -> dictionary
  graph,
  /// DOT compass point such as `"n"`, `"s"`, `"e"`, or `"w"`. -> string
  compass,
) = {
  _plugin.graph_archived_compass_subgraph(graph-module.graph-bytes(graph), cbor.encode(compass))
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

#let _edge-half(edge, key) = {
  let half = edge.at(key + "-half-edge", default: none)
  if half == none {
    edge.at(key, default: none)
  } else {
    half
  }
}

#let _half-edge-count(graph) = {
  let count = 0
  for edge in graph-module.edges(graph) {
    for endpoint in (_edge-half(edge, "source"), _edge-half(edge, "sink")) {
      if endpoint != none {
        count = calc.max(count, endpoint.hedge + 1)
      }
    }
  }
  count
}

/// Return the complement of a subgraph's selected half edges.
/// -> bytes
#let complement(
  /// Graph object whose half-edge order defines the result. -> dictionary
  graph,
  /// Subgraph object to invert. -> bytes
  selected,
) = {
  let selected-hedges = hedges(selected)
  bits(graph, range(_half-edge-count(graph)).map(i => not selected-hedges.contains(i)))
}

/// Test whether either half-edge of an edge record is selected.
/// -> bool
#let contains-edge(
  /// Subgraph object to inspect. -> bytes
  subgraph,
  /// Edge record from `graph.edges(...)` or draw callback data. -> dictionary
  edge,
) = {
  let source = _edge-half(edge, "source")
  let sink = _edge-half(edge, "sink")
  let source-selected = source != none and contains(subgraph, source.hedge)
  let sink-selected = sink != none and contains(subgraph, sink.hedge)
  source-selected or sink-selected
}

/// Compute breadth-first node depths inside a selected edge set.
/// -> dictionary
#let node-depths(
  /// Graph object to inspect. -> dictionary
  graph,
  /// Subgraph whose selected half edges define the traversal edges. -> bytes
  selected,
  /// Root node index. -> int
  root: 0,
) = {
  let adjacency = (:)
  for edge in graph-module.edges(graph) {
    let source = _edge-half(edge, "source")
    let sink = _edge-half(edge, "sink")
    if source != none and sink != none and contains-edge(selected, edge) {
      let source-key = str(source.node)
      let sink-key = str(sink.node)
      adjacency.insert(source-key, adjacency.at(source-key, default: ()) + (sink.node,))
      adjacency.insert(sink-key, adjacency.at(sink-key, default: ()) + (source.node,))
    }
  }

  let depths = (:)
  let frontier = (root,)
  depths.insert(str(root), 0)
  while frontier.len() > 0 {
    let next = ()
    for node in frontier {
      for child in adjacency.at(str(node), default: ()) {
        if not (str(child) in depths) {
          depths.insert(str(child), depths.at(str(node)) + 1)
          next.push(child)
        }
      }
    }
    frontier = next
  }
  depths
}
