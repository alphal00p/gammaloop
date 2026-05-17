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

#let _kind(item) = if type(item) == dictionary { item.at("linnest-kind", default: none) } else { none }

#let _is-label(value) = type(value) == label

#let _expect-no-extra-args(context_, args) = {
  let named = args.named()
  if named.len() > 0 {
    panic(context_ + ": unknown named argument " + repr(named.keys().first()))
  }
}

#let _label-key(value, context_) = {
  if type(value) != label {
    panic(context_ + ": expected a Typst label")
  }
  str(value)
}

#let _check-name(value, context_) = {
  if value != none and type(value) != label {
    panic(context_ + ": name must be a Typst label")
  }
}

#let _check-id(value, context_) = {
  if value != none and type(value) != int {
    panic(context_ + ": id must be an integer")
  }
}

#let _collect-named-id(name, id, used, names, context_) = {
  _check-name(name, context_)
  _check-id(id, context_)
  if id != none {
    if used.contains(id) {
      panic(context_ + ": duplicate id " + repr(id))
    }
    used.push(id)
  }
  if name != none {
    let key = str(name)
    for entry in names {
      if entry.name == key {
        panic(context_ + ": duplicate name " + repr(name))
      }
    }
    names.push((name: key, id: id))
  }
  (used: used, names: names)
}

#let _allocate-named-ids(used, names) = {
  let resolved = (:)
  let next = 0
  for entry in names {
    let id = entry.id
    if id == none {
      while used.contains(next) {
        next += 1
      }
      id = next
      used.push(next)
    }
    resolved.insert(entry.name, id)
  }
  resolved
}

#let _resolve-named-id(name, id, resolved, context_) = {
  _check-name(name, context_)
  _check-id(id, context_)
  if name == none {
    id
  } else {
    let key = str(name)
    if key not in resolved.keys() {
      panic(context_ + ": unresolved name " + repr(name))
    }
    let resolved-id = resolved.at(key)
    if id != none and id != resolved-id {
      panic(context_ + ": name/id mismatch for " + repr(name))
    }
    resolved-id
  }
}

#let _flatten-items(values) = {
  let result = ()
  for value in values {
    if value == none {
      continue
    } else if type(value) == array {
      result += _flatten-items(value)
    } else {
      result.push(value)
    }
  }
  result
}

#let _node-name-map(nodes) = {
  let map = (:)
  for (index, item) in nodes.enumerate() {
    let name = item.at("name", default: none)
    if type(name) == label {
      let key = _label-key(name, "graph.node")
      if map.keys().contains(key) {
        panic("graph.node: duplicate node name <" + key + ">")
      }
      map.insert(key, index)
    }
  }
  map
}

#let _resolve-node-ref(value, node-keys, context_) = {
  if type(value) == int {
    value
  } else if type(value) == label {
    let key = str(value)
    if key not in node-keys.keys() {
      panic(context_ + ": unknown node name <" + key + ">")
    }
    node-keys.at(key)
  } else {
    panic(context_ + ": node reference must be an integer or Typst label name")
  }
}

#let _resolve-pos(pos, node-keys) = {
  if type(pos) == dictionary and pos.at("ref", default: none) != none {
    let resolved = pos
    resolved.ref = _resolve-node-ref(pos.ref, node-keys, "graph.pos")
    resolved
  } else {
    pos
  }
}

#let _collect-node-ids(nodes) = {
  let used = ()
  let names = ()
  for item in nodes {
    let name = if type(item.at("name", default: none)) == label { item.name } else { none }
    let state = _collect-named-id(name, item.at("id", default: none), used, names, "graph.node")
    used = state.used
    names = state.names
  }
  _allocate-named-ids(used, names)
}

#let _collect-edge-ids(edges) = {
  let used = ()
  let names = ()
  for item in edges {
    let state = _collect-named-id(item.at("name", default: none), item.at("id", default: none), used, names, "graph.edge")
    used = state.used
    names = state.names
  }
  _allocate-named-ids(used, names)
}

#let _collect-half-ids(edges) = {
  let used = ()
  let names = ()
  for item in edges {
    for side in ("source", "sink") {
      let half = item.at(side, default: none)
      if half != none {
        let state = _collect-named-id(half.at("name", default: none), half.at("id", default: none), used, names, "graph." + side)
        used = state.used
        names = state.names
      }
    }
  }
  _allocate-named-ids(used, names)
}

#let _node-spec(node) = {
  let statements = _statements-with-point(node.at("statements", default: (:)), "shift", node.at("shift", default: none))
  let clean = node
  let name = clean.at("name", default: none)
  for key in ("linnest-kind", "key", "id", "shift", "name") {
    if clean.keys().contains(key) {
      let _ = clean.remove(key)
    }
  }
  if name != none {
    if type(name) == label {
      clean.insert("name", str(name))
    } else {
      clean.insert("name", name)
    }
  }
  clean + (statements: statements)
}

#let _resolved-node-spec(node, node-ids, node-keys) = {
  let name = if type(node.at("name", default: none)) == label { node.name } else { none }
  let index = _resolve-named-id(name, node.at("id", default: none), node-ids, "graph.node")
  let item = node + (pos: _resolve-pos(node.at("pos", default: none), node-keys))
  if index != none {
    item.insert("index", index)
  }
  _node-spec(item)
}

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
    "linnest-kind",
    "shift",
    "label-pos",
    "label-angle",
    "bend",
    "source-style-eval",
    "sink-style-eval",
    "label-eval",
    "name",
  ) {
    if clean.keys().contains(key) {
      let _ = clean.remove(key)
    }
  }
  clean + (statements: statements)
}

#let _resolved-half-spec(half, node-keys, half-ids, context_) = {
  if half == none {
    return none
  }
  let result = (
    node: _resolve-node-ref(half.node, node-keys, context_),
  )
  let id = _resolve-named-id(half.at("name", default: none), half.at("id", default: none), half-ids, context_)
  if id != none {
    result.insert("id", id)
  }
  if half.at("statement", default: none) != none {
    result.insert("statement", _statement-value(half.statement))
  }
  if half.at("compass", default: none) != none {
    result.insert("compass", _statement-value(half.compass))
  }
  result
}

#let _resolved-edge-spec(edge, node-keys, edge-ids, half-ids) = {
  let id = _resolve-named-id(edge.at("name", default: none), edge.at("id", default: none), edge-ids, "graph.edge")
  let item = edge + (
    source: _resolved-half-spec(edge.at("source", default: none), node-keys, half-ids, "graph.source"),
    sink: _resolved-half-spec(edge.at("sink", default: none), node-keys, half-ids, "graph.sink"),
    pos: _resolve-pos(edge.at("pos", default: none), node-keys),
  )
  if id != none {
    item.id = id
  }
  _edge-spec(item)
}

#let _split-items(items, nodes, edges) = {
  let resolved-nodes = _flatten-items(nodes)
  let resolved-edges = _flatten-items(edges)
  for item in _flatten-items(items) {
    let kind = _kind(item)
    if kind == "node" {
      resolved-nodes.push(item)
    } else if kind == "edge" {
      resolved-edges.push(item)
    } else {
      panic("graph.build: expected graph.node(..) or graph.edge(..), got " + repr(item))
    }
  }
  (nodes: resolved-nodes, edges: resolved-edges)
}

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

/// Build one graph object from a stream of node and edge items.
///
/// Use @node, @source, @sink, and @edge to create graph items. Positional items
/// may be passed as comma-separated arguments or yielded from a Typst code
/// block.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>), <e>, graph.sink(<b>))
/// }, name: "demo")
/// #graph.info(g).name
/// ```
/// -> bytes
#let build(
  /// Node and edge items returned by @node and @edge. -> array
  ..items,

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

  /// Additional node items or raw node specs. -> array
  /// -> array
  nodes: (),

  /// Additional edge items or raw edge specs. -> array
  edges: (),
) = {
  let split = _split-items(items.pos(), nodes, edges)
  let node-keys = _node-name-map(split.nodes)
  let node-ids = _collect-node-ids(split.nodes)
  let edge-ids = _collect-edge-ids(split.edges)
  let half-ids = _collect-half-ids(split.edges)
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
    nodes: split.nodes.map(node => _resolved-node-spec(node, node-ids, node-keys)),
    edges: split.edges.map(edge => _resolved-edge-spec(edge, node-keys, edge-ids, half-ids)),
  )))
}

/// Create a graph node item for @build.
///
/// A Typst label is the node name used by @source, @sink, and @pos. The
/// optional numeric `id` chooses the node order/index. `label` is display
/// metadata stored in the DOT `label` statement.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>, id: 0, label: "A")
/// })
/// #graph.nodes(g).first().name
/// ```
/// -> array
#let node(
  ..args,
  /// Typst node name for references. -> none | label
  name: none,
  /// Numeric node order/index override. -> none | int
  id: none,
  /// Visible/DOT node label, stored in `statements.label`. -> none | string
  label: none,
  /// Node placement. -> none | dictionary
  pos: none,
  /// Drawing shift stored as a statement. -> none | string | array | dictionary
  shift: none,
  /// Additional DOT node statements. -> dictionary
  statements: (:),
) = {
  _expect-no-extra-args("graph.node", args)
  let pos-args = args.pos()
  if pos-args.len() > 1 {
    panic("graph.node: expected at most one positional name")
  }
  let resolved-name = name
  if pos-args.len() == 1 {
    let arg = pos-args.first()
    if _is-label(arg) {
      if resolved-name != none {
        panic("graph.node: name specified twice")
      }
      resolved-name = arg
    } else {
      panic("graph.node: positional argument must be a Typst label name")
    }
  }
  _check-name(resolved-name, "graph.node")
  _check-id(id, "graph.node")
  let resolved-statements = statements
  if label != none {
    if resolved-statements.keys().contains("label") {
      panic("graph.node: label specified twice")
    }
    resolved-statements.insert("label", _statement-value(label))
  }
  ((
    linnest-kind: "node",
    id: id,
    name: resolved-name,
    pos: pos,
    shift: shift,
    statements: resolved-statements,
  ),)
}

/// Create a source half-edge endpoint.
///
/// `node` may be a node name like `<a>` or a numeric node index. `name` gives
/// the half-edge a Typst name; `id` is a numeric half-edge order/index
/// override.
///
/// ```example
/// #graph.source(<a>, name: <h1>, id: 0, compass: "e")
/// ```
/// -> dictionary
#let source(node, name: none, id: none, statement: none, compass: none) = {
  _check-name(name, "graph.source")
  _check-id(id, "graph.source")
  (linnest-kind: "source", node: node, name: name, id: id, statement: statement, compass: compass)
}

/// Create a sink half-edge endpoint.
///
/// `node` may be a node name like `<a>` or a numeric node index. `name` gives
/// the half-edge a Typst name; `id` is a numeric half-edge order/index
/// override.
///
/// ```example
/// #graph.sink(<b>, name: <h2>, id: 1, compass: "w")
/// ```
/// -> dictionary
#let sink(node, name: none, id: none, statement: none, compass: none) = {
  _check-name(name, "graph.sink")
  _check-id(id, "graph.sink")
  (linnest-kind: "sink", node: node, name: name, id: id, statement: statement, compass: compass)
}

/// Create a graph edge item for @build.
///
/// Positional arguments may contain one @source, one @sink, and optionally one
/// Typst label used as the edge name. The numeric `id` chooses the edge order
/// and the visible/DOT label is `label`.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>), <e>, graph.sink(<b>), label: "p")
/// })
/// ```
/// -> array
#let edge(
  /// Source/sink half-edges and optional edge name. -> any
  ..args,
  /// Typst edge name. -> none | label
  name: none,
  /// Numeric edge order/index override. -> none | int
  id: none,
  /// Edge orientation: `"default"`, `"reversed"`, or `"undirected"`. -> string
  orientation: "default",
  /// Visible/DOT edge label, stored in `statements.label`. -> none | string
  label: none,
  /// Edge placement. -> none | dictionary
  pos: none,
  /// Drawing shift stored as a statement. -> none | string | array | dictionary
  shift: none,
  /// Edge label position stored as a statement. -> none | string | array | dictionary
  label-pos: none,
  /// Edge label angle stored as a statement. -> none | int | float | string
  label-angle: none,
  /// Edge bend stored as a statement. -> none | int | float | string
  bend: none,
  /// Additional DOT edge statements. -> dictionary
  statements: (:),
  /// Always-evaluated Typst style for the source half-edge. -> none | string
  source-style-eval: none,
  /// Always-evaluated Typst style for the sink half-edge. -> none | string
  sink-style-eval: none,
  /// Always-evaluated Typst edge label template. -> none | string
  label-eval: none,
) = {
  _expect-no-extra-args("graph.edge", args)
  let resolved-id = id
  let resolved-name = name
  let resolved-source = none
  let resolved-sink = none
  for arg in args.pos() {
    let kind = _kind(arg)
    if kind == "source" {
      if resolved-source != none {
        panic("graph.edge: source specified twice")
      }
      resolved-source = arg
    } else if kind == "sink" {
      if resolved-sink != none {
        panic("graph.edge: sink specified twice")
      }
      resolved-sink = arg
    } else if _is-label(arg) {
      if resolved-name != none {
        panic("graph.edge: name specified twice")
      }
      resolved-name = arg
    } else {
      panic("graph.edge: positional arguments must be graph.source(..), graph.sink(..), or an edge name")
    }
  }
  _check-name(resolved-name, "graph.edge")
  _check-id(resolved-id, "graph.edge")
  if resolved-source == none and resolved-sink == none {
    panic("graph.edge: expected a source or sink half-edge")
  }
  let flow = if resolved-source != none and resolved-sink == none {
    "source"
  } else if resolved-source == none and resolved-sink != none {
    "sink"
  } else {
    none
  }
  let resolved-statements = statements
  if label != none {
    if resolved-statements.keys().contains("label") {
      panic("graph.edge: label specified twice")
    }
    resolved-statements.insert("label", _statement-value(label))
  }
  ((
    linnest-kind: "edge",
    source: resolved-source,
    sink: resolved-sink,
    orientation: orientation,
    flow: flow,
    name: resolved-name,
    id: resolved-id,
    pos: pos,
    shift: shift,
    label-pos: label-pos,
    label-angle: label-angle,
    bend: bend,
    statements: resolved-statements,
    source-style-eval: source-style-eval,
    sink-style-eval: sink-style-eval,
    label-eval: label-eval,
  ),)
}

/// Return graph metadata.
///
/// The result has `name`, `global-statements`, `edge-statements`, and
/// `node-statements`.
///
/// ```example
/// #let g = graph.build({ graph.node(<a>) }, name: "demo")
/// #graph.info(g).name
/// ```
/// -> dictionary
#let info(graph) = cbor(_plugin.graph_info(graph))

/// Serialize a graph object to DOT.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>), graph.sink(<b>))
/// }, name: "demo")
/// #graph.dot(g).contains("digraph demo")
/// ```
/// -> string
#let dot(graph) = cbor(_plugin.graph_dot(bytes(graph)))

/// Return node records, optionally filtered by an subgraph object.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
/// })
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
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>, compass: "e"), graph.sink(<b>))
/// })
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
/// Supported key values are `"statement"`, `"compass"`, and `"id"`.
///
/// ```example
/// #let left = graph.build({
///   graph.node(<a>)
///   graph.edge(graph.sink(<a>, statement: "j"))
/// })
/// #let right = graph.build({
///   graph.node(<b>)
///   graph.edge(graph.source(<b>, statement: "j"))
/// })
/// #graph.edges(graph.join(left, right, key: "statement")).len()
/// ```
/// -> bytes
#let join(left, right, key: "statement") = {
  _plugin.graph_join_by_hedge_key(bytes(left), bytes(right), cbor.encode((key: key)))
}

/// Return subgraph objects for the graph's cycle basis.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>), graph.sink(<b>))
/// })
/// #graph.cycles(g).len()
/// ```
/// -> array
#let cycles(graph) = cbor(_plugin.graph_cycle_basis(bytes(graph)))

/// Return subgraph objects for the graph's spanning forests.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>), graph.sink(<b>))
/// })
/// #graph.forests(g).len()
/// ```
/// -> array
#let forests(graph) = cbor(_plugin.graph_spanning_forests(bytes(graph)))
