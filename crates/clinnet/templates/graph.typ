#let _plugin = plugin("./linnest.wasm")

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

#let _statements-with-pin(statements, pin) = {
  if pin == none {
    statements
  } else {
    statements + (pin: _statement-value(pin))
  }
}

#let _node-spec(node) = {
  let statements = _statements-with-pin(
    node.at("statements", default: (:)),
    node.at("pin", default: none),
  )
  let clean = node
  if clean.keys().contains("pin") {
    let _ = clean.remove("pin")
  }
  clean + (statements: statements)
}

#let _edge-spec(edge) = {
  let statements = _edge-render-statements(
    _statements-with-pin(
      edge.at("statements", default: (:)),
      edge.at("pin", default: none),
    ),
    source-style-eval: edge.at("source-style-eval", default: none),
    sink-style-eval: edge.at("sink-style-eval", default: none),
    label-eval: edge.at("label-eval", default: none),
  )
  let clean = edge
  for key in ("pin", "source-style-eval", "sink-style-eval", "label-eval") {
    if clean.keys().contains(key) {
      let _ = clean.remove(key)
    }
  }
  clean + (statements: statements)
}

/// Format a layout pin statement.
///
/// Pins constrain node positions or edge control-point positions during
/// layout. Numeric `x` and `y` values fix coordinates. Use @pin-group for
/// grouped constraints such as shared external-edge columns.
///
/// ```example
/// #graph.pin(x: -2, y: 0)
/// ```
/// -> string
#let pin(x: none, y: none, group: none) = {
  if group != none {
    if x != none or y != none {
      panic("graph.pin: group cannot be combined with x or y")
    }
    _statement-value(group)
  } else {
    let parts = ()
    if x != none {
      parts.push("x:" + _statement-value(x))
    }
    if y != none {
      parts.push("y:" + _statement-value(y))
    }
    if parts.len() == 0 {
      panic("graph.pin requires x, y, or group")
    }
    parts.join(",")
  }
}

/// Format a grouped layout-pin reference.
///
/// `side: "+"` keeps the final coordinate non-negative and `side: "-"`
/// keeps it non-positive. Omit `side` to link the coordinate without a sign
/// constraint.
///
/// ```example
/// #graph.pin(x: graph.pin-group("right", side: "+"), y: graph.pin-group("row-1"))
/// ```
/// -> string
#let pin-group(name, side: none) = {
  let prefix = if side == none {
    "@"
  } else if side == "+" or side == "positive" {
    "@+"
  } else if side == "-" or side == "negative" {
    "@-"
  } else {
    panic("graph.pin-group: side must be none, \"+\", \"-\", \"positive\", or \"negative\"")
  }
  prefix + _statement-value(name)
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
/// This is convenience sugar over @builder, @node, @edge, and @finish.
///
/// ```example
/// #let g = graph.build(
///   name: "demo",
///   nodes: ((name: "a"), (name: "b")),
///   edges: ((source: (node: 0), sink: (node: 1)),),
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

  /// Node specifications. Each node may define `name`, `index`, `pin`, and
  /// `statements`. This matches the per-node parameters accepted by @node.
  /// -> array
  nodes: (),

  /// Edge specifications. Each edge may define `source`, `sink`,
  /// `orientation`, `flow`, `id`, `pin`, `statements`, `source-style-eval`,
  /// `sink-style-eval`, and `label-eval`. The `source` and `sink` fields use
  /// the same endpoint dictionaries accepted by @edge. -> array
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

/// Create a builder object.
///
/// `statements`, `node-statements`, and `edge-statements` set graph-level DOT
/// attributes that are carried into the finished graph. String values may use
/// `{name}` placeholders; edge defaults expand after per-edge statements are
/// merged, so an `edge-statements`, `source-style-eval`, or `sink-style-eval` value can
/// refer to a `label` supplied by @edge.
///
/// ```example
/// #let b = graph.builder(
///   name: "demo",
///   statements: (full_num: "x + y"),
///   edge-statements: (display-label: "{label}"),
///   source-style-eval: "(stroke: red + 0.5pt)",
///   sink-style-eval: "(stroke: blue + 0.5pt)",
/// )
/// #let (node: a, builder: b) = graph.node(b, name: "a")
/// #a
/// ```
/// -> bytes
#let builder(
  name: none,
  statements: (:),
  node-statements: (:),
  edge-statements: (:),
  source-style-eval: none,
  sink-style-eval: none,
  label-eval: none,
) = {
  let edge-statements = _edge-render-statements(
    edge-statements,
    source-style-eval: source-style-eval,
    sink-style-eval: sink-style-eval,
    label-eval: label-eval,
  )
  _plugin.graph_builder(cbor.encode((
    name: name,
    statements: statements,
    node-statements: node-statements,
    edge-statements: edge-statements,
  )))
}

/// Add a node to a builder object.
///
/// Returns a dictionary with fields `node` and `builder`, so callers can write
/// `#let (node: a, builder: b) = graph.node(b, name: "a")`. The returned
/// builder is the value to pass to the next @node, @edge, or @finish call.
///
/// ```example
/// #let b = graph.builder(name: "demo")
/// #let (node: a, builder: b) = graph.node(b, name: "a")
/// #let (node: c, builder: b) = graph.node(b, name: "c")
/// #repr((a, c))
/// ```
/// -> dictionary
#let node(builder, name: none, index: none, pin: none, statements: (:)) = {
  cbor(_plugin.graph_builder_add_node(
    bytes(builder),
    cbor.encode((
      name: name,
      index: index,
      statements: _statements-with-pin(statements, pin),
    )),
  ))
}

/// Add a paired or external edge to a builder object.
///
/// Set `source` to `none` or `sink` to `none` for an external edge.
/// Endpoint dictionaries accept `node`, `statement`, `id`, `port-label`,
/// `compass`, and `in-subgraph`.
///
/// ```example
/// #let b = graph.builder(name: "demo")
/// #let (node: a, builder: b) = graph.node(b, name: "a")
/// #let (node: c, builder: b) = graph.node(b, name: "c")
/// #let b = graph.edge(
///   b,
///   source: (node: a, compass: "e"),
///   sink: (node: c, compass: "w"),
///   statements: (label: "a-c"),
/// )
/// #graph.edges(graph.finish(b)).len()
/// ```
/// -> bytes
#let edge(
  builder,
  source: none,
  sink: none,
  orientation: "default",
  flow: none,
  id: none,
  pin: none,
  statements: (:),
  source-style-eval: none,
  sink-style-eval: none,
  label-eval: none,
) = {
  let statements = _edge-render-statements(
    _statements-with-pin(statements, pin),
    source-style-eval: source-style-eval,
    sink-style-eval: sink-style-eval,
    label-eval: label-eval,
  )
  _plugin.graph_builder_add_edge(
    bytes(builder),
    cbor.encode((
      source: source,
      sink: sink,
      orientation: orientation,
      flow: flow,
      id: id,
      statements: statements,
    )),
  )
}

/// Finish a builder object into a graph object.
///
/// ```example
/// #let b = graph.builder(name: "demo")
/// #let (node: a, builder: b) = graph.node(b, name: "a")
/// #let g = graph.finish(b)
/// #graph.info(g).name
/// ```
/// -> bytes
#let finish(builder) = _plugin.graph_builder_finish(bytes(builder))

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
