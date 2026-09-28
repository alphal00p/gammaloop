// Subgraphs pair native selection bytes with Typst-only data and hedge-data
// (keyed by str(hedge ID)). Topology signatures bind selections to incidence
// and names, not graph identity or drawing state; raw bytes are not public inputs.

#import "graph.typ" as graph-module
#import "impl/subgraph.typ" as _impl

#let _plugin = plugin("../linnest.wasm")

/// Construct a subgraph object from a base62 label.
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
/// -> dictionary
#let label(
  /// Graph object whose half-edge set the label refers to. -> dictionary
  graph,
  /// Base62 subgraph label returned by @to-label or produced by Linnest. -> string
  label,
) = _impl.wrap(
  _plugin.graph_archived_subgraph(
    graph-module.graph-bytes(graph),
    cbor.encode(label),
  ),
  topology: _impl.topology(graph),
)

/// Construct a subgraph object from a boolean hedge array.
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
/// -> dictionary
#let bits(
  /// Graph object whose half-edge order defines the bit array. -> dictionary
  graph,
  /// Boolean array selecting half edges by graph half-edge index. -> array
  bits,
) = _impl.wrap(
  _plugin.graph_archived_subgraph(
    graph-module.graph-bytes(graph),
    cbor.encode(bits),
  ),
  topology: _impl.topology(graph),
)

/// Construct a subgraph object from a DOT compass point such as `"n"` or
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
/// -> dictionary
#let compass(
  /// Graph object to filter by half-edge compass statement. -> dictionary
  graph,
  /// DOT compass point such as `"n"`, `"s"`, `"e"`, or `"w"`. -> string
  compass,
) = _impl.wrap(
  _plugin.graph_archived_compass_subgraph(
    graph-module.graph-bytes(graph),
    cbor.encode(compass),
  ),
  topology: _impl.topology(graph),
)

/// Convert a subgraph object to its base62 label.
///
/// The label encodes only membership, not topology or annotations.
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
  /// Subgraph object returned by this module or by `graph.cycles`/`graph.forests`. -> dictionary
  subgraph,
) = cbor(_plugin.subgraph_label(_impl.subgraph-bytes(subgraph)))

/// Return the hedge indices included in a subgraph object.
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
  /// Subgraph object to inspect. -> dictionary
  subgraph,
) = cbor(_plugin.subgraph_hedges(_impl.subgraph-bytes(subgraph)))

/// Test whether a subgraph object includes a hedge.
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
  /// Subgraph object to inspect. -> dictionary
  subgraph,
  /// Half-edge index to test. -> int
  hedge,
) = {
  cbor(_plugin.subgraph_contains_hedge(
    _impl.subgraph-bytes(subgraph),
    cbor.encode(hedge),
  ))
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

/// Select the union of node incidences, whole edges, edge halves, and exact hedges.
///
/// Node and edge references are Typst label names or numeric IDs. `source` and
/// `sink` refer to the structural halves, independently of drawing orientation.
/// Missing references or absent requested halves are errors. An isolated node
/// contributes no hedges: subgraphs represent half-edge sets, not node sets.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.edge(<D1>, graph.source(<a>))
///   graph.edge(<D6>, graph.sink(<a>))
/// })
/// #let left = subgraph.select(g, source: (<D1>,), sink: (<D6>,))
/// #subgraph.hedges(left)
/// ```
/// -> dictionary
#let select(
  /// Graph whose records define names and IDs. -> dictionary
  g,
  /// Select all half edges incident to these nodes. -> array
  nodes: (),
  /// Select both available halves of these edges. -> array
  edges: (),
  /// Select only the source half of these edges. -> array
  source: (),
  /// Select only the sink half of these edges. -> array
  sink: (),
  /// Select these exact numeric half-edge IDs. -> array
  hedges: (),
) = {
  let node-records = graph-module.nodes(g)
  let edge-records = graph-module.edges(g)
  let resolved = (:)
  for (kind, refs) in (nodes: nodes, edges: edges, source: source, sink: sink) {
    assert(
      type(refs) == array,
      message: "subgraph.select: " + kind + " must be an array",
    )
    let records = if kind == "nodes" { node-records } else { edge-records }
    let field = if kind == "nodes" { "node" } else { "edge" }
    resolved.insert(kind, refs.map(value => {
      assert(
        type(value) == int or type(value) == type(<name>),
        message: "subgraph.select: "
          + field
          + " reference must be an integer or Typst label",
      )
      let record = records.find(record => if type(value) == int {
        record.at(field) == value
      } else {
        record.at("name", default: none) == value
      })
      assert(
        record != none,
        message: "subgraph.select: unknown " + field + " " + repr(value),
      )
      record.at(field)
    }))
  }

  let selected = (:)
  let all-hedges = ()
  for edge in edge-records {
    for flow in ("source", "sink") {
      let half = _edge-half(edge, flow)
      let requested = edge.edge in resolved.at(flow)
      assert(
        not requested or half != none,
        message: "subgraph.select: edge "
          + str(edge.edge)
          + " has no "
          + flow
          + " half-edge",
      )
      if half != none {
        all-hedges.push(half.hedge)
        if (
          requested
            or edge.edge in resolved.edges
            or half.node in resolved.nodes
        ) {
          selected.insert(str(half.hedge), true)
        }
      }
    }
  }
  assert(
    type(hedges) == array,
    message: "subgraph.select: hedges must be an array",
  )
  for hedge in hedges {
    assert(
      type(hedge) == int and hedge in all-hedges,
      message: "subgraph.select: unknown half-edge ID " + repr(hedge),
    )
    selected.insert(str(hedge), true)
  }
  bits(g, range(all-hedges.len()).map(i => str(i) in selected))
}

/// Attach arbitrary Typst data without serializing it through Wasm.
///
/// `data: auto` preserves subgraph-wide data; `hedge: none` preserves annotations.
/// A hedge callback receives the original endpoint record plus numeric `edge`,
/// label-valued `edge-name` (or `none`), and `flow` (`"source"` or `"sink"`).
/// Its `data` is the existing subgraph annotation; `graph-data` holds the original
/// half-edge data. The return value replaces the annotation, including `none`.
/// Other non-function hedge values replace every selected hedge's annotation.
/// Only selected hedges retain annotations; callback values are never evaluated
/// unless supplied as the `hedge` updater itself.
/// -> dictionary
#let with-data(
  /// Compatible graph supplying endpoint records and graph data. -> dictionary
  g,
  /// Subgraph to annotate. -> dictionary
  s,
  /// Replacement subgraph-wide data, or `auto` to preserve it. -> any
  data: auto,
  /// Half-edge updater, constant annotation, or `none` to preserve annotations. -> any
  hedge: none,
) = {
  let s = _impl.validate(g, s)
  let annotations = s.hedge-data
  if hedge != none {
    annotations = (:)
    let selected = hedges(s)
    for edge in graph-module.edges(g) {
      for flow in ("source", "sink") {
        let half = _edge-half(edge, flow)
        if half != none and half.hedge in selected {
          let key = str(half.hedge)
          let value = if type(hedge) == function {
            hedge(
              half
                + (
                  edge: edge.edge,
                  edge-name: edge.at("name", default: none),
                  flow: flow,
                  data: s.hedge-data.at(key, default: none),
                  graph-data: half.at("data", default: none),
                ),
            )
          } else { hedge }
          annotations.insert(key, value)
        }
      }
    }
  }
  _impl.wrap(
    s.bytes,
    data: if data == auto { s.data } else { data },
    hedge-data: annotations,
    topology: s.topology,
  )
}

/// Return a selected hedge's annotation, or `none` when it has no annotation.
/// Unselected and out-of-range half-edge IDs are errors.
/// -> any
#let hedge-data(
  /// Subgraph whose annotation is read. -> dictionary
  s,
  /// Selected numeric half-edge ID. -> int
  h,
) = {
  assert(
    type(h) == int and h in hedges(s),
    message: "subgraph.hedge-data: half-edge ID is not selected: " + repr(h),
  )
  s.hedge-data.at(str(h), default: none)
}

/// Return the complement of a subgraph's selected half edges.
///
/// Subgraph-wide data is preserved; newly selected hedges have no annotations.
/// -> dictionary
#let complement(
  /// Graph object whose half-edge order defines the result. -> dictionary
  graph,
  /// Subgraph object to invert. -> dictionary
  selected,
) = {
  let selected = _impl.validate(graph, selected)
  let selected-hedges = hedges(selected)
  let result = bits(
    graph,
    range(_half-edge-count(graph)).map(i => not selected-hedges.contains(i)),
  )
  result.data = selected.data
  result
}

/// Test whether either half-edge of an edge record is selected.
/// -> bool
#let contains-edge(
  /// Subgraph object to inspect. -> dictionary
  subgraph,
  /// Edge record from `graph.edges(...)` or draw callback data. -> dictionary
  edge,
) = {
  let selected = hedges(subgraph)
  let source = _edge-half(edge, "source")
  let sink = _edge-half(edge, "sink")
  (
    (source != none and source.hedge in selected)
      or (sink != none and sink.hedge in selected)
  )
}

/// Compute breadth-first node depths inside a selected edge set.
/// -> dictionary
#let node-depths(
  /// Graph object to inspect. -> dictionary
  graph,
  /// Subgraph whose selected half edges define the traversal edges. -> dictionary
  selected,
  /// Root node index. -> int
  root: 0,
) = {
  let selected = _impl.validate(graph, selected)
  let adjacency = (:)
  for edge in graph-module.edges(graph) {
    let source = _edge-half(edge, "source")
    let sink = _edge-half(edge, "sink")
    if source != none and sink != none and contains-edge(selected, edge) {
      let source-key = str(source.node)
      let sink-key = str(sink.node)
      adjacency.insert(
        source-key,
        adjacency.at(source-key, default: ()) + (sink.node,),
      )
      adjacency.insert(
        sink-key,
        adjacency.at(sink-key, default: ()) + (source.node,),
      )
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
