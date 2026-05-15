#let _plugin = plugin("../linnest.wasm")

#let _edge-render-statements(
  statements,
  source-style-eval: none,
  sink-style-eval: none,
  label-eval: none,
) = {
  let result = statements
  if source-style-eval != none {
    result = result + (source-style-eval: source-style-eval)
  }
  if sink-style-eval != none {
    result = result + (sink-style-eval: sink-style-eval)
  }
  if label-eval != none {
    result = result + (label-eval: label-eval)
  }
  result
}

#let _statement-value(value) = if type(value) == str { value } else { str(value) }

#let _point-statement(point) = {
  if type(point) == str {
    point
  } else if type(point) == dictionary and point.keys().contains("x") and point.keys().contains("y") {
    str(point.x) + "," + str(point.y)
  } else if type(point) == array and point.len() == 2 {
    str(point.at(0)) + "," + str(point.at(1))
  } else {
    panic("graph position values must be strings, (x:, y:) dictionaries, or two-item arrays")
  }
}

#let _statements-with-point(statements, key, point) = {
  let result = statements
  if point != none {
    result.insert(key, _point-statement(point))
  }
  result
}

#let _statements-with-value(statements, key, value) = {
  let result = statements
  if value != none {
    result.insert(key, _statement-value(value))
  }
  result
}

#let _node-spec(node) = {
  let statements = _statements-with-point(node.at("statements", default: (:)), "shift", node.at("shift", default: none))
  let clean = node
  for key in ("shift",) {
    if clean.keys().contains(key) {
      let _ = clean.remove(key)
    }
  }
  clean + (statements: statements)
}

#let _node-args(name: none, index: none, pos: none, shift: none, statements: (:)) = _node-spec((
  name: name,
  index: index,
  pos: pos,
  shift: shift,
  statements: statements,
))

#let _edge-statements(edge) = {
  let statements = _statements-with-value(
    _statements-with-value(
      _statements-with-point(
        _statements-with-point(
          edge.at("statements", default: (:)),
          "shift",
          edge.at("shift", default: none),
        ),
        "label-pos",
        edge.at("label-pos", default: none),
      ),
      "label-angle",
      edge.at("label-angle", default: none),
    ),
    "bend",
    edge.at("bend", default: none),
  )
  _edge-render-statements(
    statements,
    source-style-eval: edge.at("source-style-eval", default: none),
    sink-style-eval: edge.at("sink-style-eval", default: none),
    label-eval: edge.at("label-eval", default: none),
  )
}

#let _edge-spec(edge) = {
  let statements = _edge-statements(edge)
  let clean = edge
  for key in (
    "shift",
    "label-pos",
    "label-angle",
    "bend",
    "source-style-eval",
    "sink-style-eval",
    "label-eval",
  ) {
    if clean.keys().contains(key) {
      let _ = clean.remove(key)
    }
  }
  clean + (statements: statements)
}

#let _edge-args(
  source: none,
  sink: none,
  orientation: "default",
  flow: none,
  id: none,
  pos: none,
  shift: none,
  label-pos: none,
  label-angle: none,
  bend: none,
  statements: (:),
  source-style-eval: none,
  sink-style-eval: none,
  label-eval: none,
) = _edge-spec((
  source: source,
  sink: sink,
  orientation: orientation,
  flow: flow,
  id: id,
  pos: pos,
  shift: shift,
  label-pos: label-pos,
  label-angle: label-angle,
  bend: bend,
  statements: statements,
  source-style-eval: source-style-eval,
  sink-style-eval: sink-style-eval,
  label-eval: label-eval,
))

/// Create a grouped placement coordinate.
///
/// `side: "+"` keeps the solved coordinate non-negative and `side: "-"`
/// keeps it non-positive. Groups are layout constraints and therefore require
/// pin placement, which is the `graph.pos` default.
///
/// ```example
/// #graph.group("right", side: "+")
/// ```
/// -> dictionary
#let group(name, side: none) = {
  if side != none and not (side == "+" or side == "-" or side == "positive" or side == "negative") {
    panic("graph.group: side must be none, \"+\", \"-\", \"positive\", or \"negative\"")
  }
  (kind: "group", name: _statement-value(name), side: side)
}

/// Create a first-class graph placement.
///
/// The default `mode: "pin"` turns numeric and grouped coordinates into layout
/// constraints and also makes the coordinates immediately drawable without a
/// layout pass. Use `mode: "start"` when the coordinate should only seed the
/// layout.
/// `ref` may reference a previously created node index and combines with `dx`
/// and `dy`.
///
/// ```example
/// #graph.pos(x: 0, y: graph.group("row"), mode: "pin")
/// ```
/// -> dictionary
#let pos(x: none, y: none, ref: none, dx: none, dy: none, mode: "pin") = {
  if mode != "start" and mode != "pin" {
    panic("graph.pos: mode must be \"start\" or \"pin\"")
  }
  let result = (mode: mode)
  if x != none {
    result.insert("x", x)
  }
  if y != none {
    result.insert("y", y)
  }
  if ref != none {
    result.insert("ref", ref)
  }
  if dx != none {
    result.insert("dx", dx)
  }
  if dy != none {
    result.insert("dy", dy)
  }
  result
}

/// Parse one or more DOT digraphs into graph objects.
///
/// ```example
/// #let graphs = graph.parse("digraph first { a -> b }")
/// #graphs.len()
/// ```
/// -> array
#let parse(input) = cbor(_plugin.parse_graph(bytes(input)))

/// Build one graph object from node and edge arrays.
///
/// Use @node and @edge to create reusable spec dictionaries, then pass those
/// arrays to @build.
///
/// ```example
/// #let g = graph.build(
///   name: "demo",
///   nodes: (graph.node(name: "a"), graph.node(name: "b")),
///   edges: (graph.edge(source: (node: 0), sink: (node: 1)),),
/// )
/// #graph.info(g).name
/// ```
/// -> bytes
#let build(
  /// Graph name. -> none | string
  name: none,

  /// Graph-level DOT statements. -> dictionary
  statements: (:),

  /// Default DOT statements for edges. String values may use `{name}`
  /// placeholders that expand from each edge's merged statement dictionary.
  /// -> dictionary
  edge-statements: (:),

  /// Default always-evaluated Typst style for source half-edges.
  /// This is merged into `edge-statements` as `source-style-eval`.
  /// -> none | string
  source-style-eval: none,

  /// Default always-evaluated Typst style for sink half-edges.
  /// This is merged into `edge-statements` as `sink-style-eval`.
  /// -> none | string
  sink-style-eval: none,

  /// Default always-evaluated Typst edge label template.
  /// This is merged into `edge-statements` as `label-eval`.
  /// -> none | string
  label-eval: none,

  /// Default DOT statements for nodes. -> dictionary
  node-statements: (:),

  /// Node specifications. Each node may define `name`, `index`, `pos`, `shift`,
  /// and `statements`. This matches the per-node parameters accepted by @node.
  /// -> array
  nodes: (),

  /// Edge specifications. Each edge may define `source`, `sink`, `orientation`,
  /// `flow`, `id`, `pos`, `shift`, `label-pos`, `label-angle`, `bend`,
  /// `statements`, `source-style-eval`, `sink-style-eval`, and `label-eval`.
  /// The `source` and `sink` fields use the same half-edge dictionaries accepted
  /// by @edge. -> array
  edges: (),
) = {
  let edge-statements = _edge-render-statements(
    edge-statements,
    source-style-eval: source-style-eval,
    sink-style-eval: sink-style-eval,
    label-eval: label-eval,
  )
  _plugin.graph_from_spec(cbor.encode((
    name: name,
    statements: statements,
    edge-statements: edge-statements,
    node-statements: node-statements,
    nodes: nodes.map(_node-spec),
    edges: edges.map(_edge-spec),
  )))
}

/// Create a node spec dictionary for @build.
///
/// `shift` is normalized into `statements.shift`.
///
/// ```example
/// #let spec = graph.node(name: "a", shift: (1, 0))
/// #spec.statements.shift
/// ```
/// -> dictionary
#let node(name: none, index: none, pos: none, shift: none, statements: (:)) = _node-args(
  name: name,
  index: index,
  pos: pos,
  shift: shift,
  statements: statements,
)

/// Create a paired or external edge spec dictionary for @build.
///
/// Set `source` to `none` or `sink` to `none` for an external edge.
/// Half-edge dictionaries accept `node`, `statement`, `id`, `port-label`,
/// `compass`, and `in-subgraph`. Drawing helpers like `shift`, `label-pos`,
/// `bend`, and `label-eval` are normalized into `statements`.
///
/// ```example
/// #let spec = graph.edge(
///   source: (node: 0, compass: "e"),
///   sink: (node: 1, compass: "w"),
///   statements: (label: "a-c"),
/// )
/// #spec.statements.label
/// ```
/// -> dictionary
#let edge(
  source: none,
  sink: none,
  orientation: "default",
  flow: none,
  id: none,
  pos: none,
  shift: none,
  label-pos: none,
  label-angle: none,
  bend: none,
  statements: (:),
  source-style-eval: none,
  sink-style-eval: none,
  label-eval: none,
) = _edge-args(
  source: source,
  sink: sink,
  orientation: orientation,
  flow: flow,
  id: id,
  pos: pos,
  shift: shift,
  label-pos: label-pos,
  label-angle: label-angle,
  bend: bend,
  statements: statements,
  source-style-eval: source-style-eval,
  sink-style-eval: sink-style-eval,
  label-eval: label-eval,
)

/// Return graph metadata.
///
/// The result has `name`, `global-statements`, `edge-statements`, and
/// `node-statements`.
///
/// ```example
/// #let g = graph.build(name: "demo", nodes: ((name: "a"),))
/// #graph.info(g).name
/// ```
/// -> dictionary
#let info(graph) = cbor(_plugin.graph_info(graph))

/// Serialize a graph object to DOT.
///
/// ```example
/// #let g = graph.build(
///   name: "demo",
///   nodes: ((name: "a"), (name: "b")),
///   edges: ((source: (node: 0), sink: (node: 1)),),
/// )
/// #graph.dot(g).contains("digraph demo")
/// ```
/// -> string
#let dot(graph) = cbor(_plugin.graph_dot(bytes(graph)))

/// Return node records, optionally filtered by an subgraph object.
///
/// ```example
/// #let g = graph.build(nodes: ((name: "a"), (name: "b")))
/// #graph.nodes(g).map(node => node.name).join(", ")
/// ```
/// -> array
#let nodes(graph, subgraph: none) = {
  if subgraph == none {
    cbor(_plugin.graph_nodes(graph))
  } else {
    cbor(_plugin.graph_nodes_of_subgraph(bytes(graph), bytes(subgraph)))
  }
}

/// Return edge records, optionally filtered by an subgraph object.
///
/// ```example
/// #let g = graph.build(
///   nodes: ((name: "a"), (name: "b")),
///   edges: ((source: (node: 0, compass: "e"), sink: (node: 1)),),
/// )
/// #let east = subgraph.compass(g, "e")
/// #graph.edges(g, subgraph: east).len()
/// ```
/// -> array
#let edges(graph, subgraph: none) = {
  if subgraph == none {
    cbor(_plugin.graph_edges(graph))
  } else {
    cbor(_plugin.graph_edges_of_subgraph(bytes(graph), bytes(subgraph)))
  }
}

/// Join two graphs by matching dangling half-edge data on `key`.
///
/// Supported key values are `"statement"`, `"port-label"`, `"compass"`, and
/// `"id"`.
///
/// ```example
/// #let left = graph.build(
///   nodes: ((name: "a"),),
///   edges: ((source: none, sink: (node: 0, statement: "j")),),
/// )
/// #let right = graph.build(
///   nodes: ((name: "b"),),
///   edges: ((source: (node: 0, statement: "j"), sink: none),),
/// )
/// #graph.edges(graph.join(left, right, key: "statement")).len()
/// ```
/// -> bytes
#let join(left, right, key: "statement") = {
  _plugin.graph_join_by_hedge_key(bytes(left), bytes(right), cbor.encode((key: key)))
}

/// Return subgraph objects for the graph's cycle basis.
///
/// ```example
/// #let g = graph.build(
///   nodes: ((name: "a"), (name: "b")),
///   edges: ((source: (node: 0), sink: (node: 1)),),
/// )
/// #graph.cycles(g).len()
/// ```
/// -> array
#let cycles(graph) = cbor(_plugin.graph_cycle_basis(bytes(graph)))

/// Return subgraph objects for the graph's spanning forests.
///
/// ```example
/// #let g = graph.build(
///   nodes: ((name: "a"), (name: "b")),
///   edges: ((source: (node: 0), sink: (node: 1)),),
/// )
/// #graph.forests(g).len()
/// ```
/// -> array
#let forests(graph) = cbor(_plugin.graph_spanning_forests(bytes(graph)))
