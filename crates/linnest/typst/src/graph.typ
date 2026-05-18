// Public graph API.
//
// Keep user-facing constructors, queries, and modifiers here. The
// serialization-heavy implementation lives in `graph-impl.typ`.

#import "graph-impl.typ" as _impl

/// Build one graph object from a stream of node and edge items.
///
/// Use @node, @source, @sink, and @edge to create graph items. Positional items
/// may be passed as comma-separated arguments or yielded from a Typst code
/// block.
/// -> bytes
#let build(
  /// Node and edge items returned by @node and @edge.
  /// The best way to use these is to use a code scope, so that the edges and nodes append each other:
  ///
  /// ```example
  /// #let g = build({
  ///   node(<a>, label: [a])
  ///   node(<b>, label: [b])
  ///   node(<c>, label: [c])
  ///   edge(source(<a>), sink(<b>))
  ///   edge(source(<b>), sink(<c>))
  ///   edge(source(<a>), sink(<c>))
  /// })
  /// >>>#align(center+horizon, draw(layout(g)))
  /// ```
  /// -> array
  ..items,
  /// Graph name.
  /// ```example
  /// #let g = build(name: "My Graph", {
  /// })
  /// #info(g).name
  /// ```
  ///  -> none | string
  name: none,
  /// Opaque graph payload. Any CBOR-encodable Typst value.
  ///
  /// ```example
  /// #let g = build(payload: (a:(b:(1, ))), {
  /// })
  /// #info(g).payload
  /// ```
  /// -> any
  payload: none,
  /// Flat graph statements, that get turned into a string to string dictionary in rust. Used by DOT parsing; values cannot nest.
  /// ```example
  /// #let g = build(statements: (a:1), {
  /// })
  /// #info(g).global-statements
  /// ```
  ///  -> dictionary
  statements: (:),
  /// Flat default node statements. Used by DOT; values cannot nest. -> dictionary
  default-node-statements: (:),
  /// Flat default edge statements. Used by DOT parsing; values cannot nest. These are applied/delegated to all edges.
  /// ```example
  /// #let g = build(default-edge-statements: (a:1), {
  ///    node(<a>)
  ///    edge(source(<a>))
  ///    edge(sink(<a>))
  /// })
  /// #edges(g).map(e=>e.statements)
  /// ```
  ///-> dictionary
  default-edge-statements: (:),
  /// Default payload merged into every node payload. Captured node payload fields override it.
  /// ```example
  /// #let g = build(default-node-payload: (a:1), {
  ///    node(<a>)
  ///    node(<b>,a:2)
  ///    node(<c>,b:(a:1))
  /// })
  /// #nodes(g).map(n=>n.payload)
  /// ```
  ///  -> any
  default-node-payload: none,
  /// Default payload merged into every edge payload. Captured edge payload fields override it.
  /// ```example
  /// #let g = build(default-edge-payload: (a:1), {
  ///    node(<a>)
  ///    edge(source(<a>))
  ///    edge(sink(<a>))
  /// })
  /// #edges(g).map(e=>e.payload)
  /// ```
  /// -> any
  default-edge-payload: none,
  /// Default payload merged into every source half-edge payload. Captured source payload fields override it.
  /// ```example
  /// #let g = build(default-source-payload: (a:1), {
  ///    node(<a>)
  ///    edge(source(<a>))
  ///    edge(sink(<a>))
  /// })
  /// #edges(g).map(e=>if e.source != none {e.source.payload} else {none})
  /// ```
  /// -> any
  default-source-payload: none,
  /// Default payload merged into every sink half-edge payload. Captured sink payload fields override it.
  /// ```example
  /// #let g = build(default-sink-payload: (a:1), {
  ///    node(<a>)
  ///    edge(source(<a>))
  ///    edge(sink(<a>))
  /// })
  /// #edges(g).map(e=>if e.sink != none {e.sink.payload} else {none})
  /// ```
  /// -> any
  default-sink-payload: none,
) = _impl.build(
  ..items,
  name: name,
  payload: payload,
  statements: statements,
  default-edge-statements: default-edge-statements,
  default-node-payload: default-node-payload,
  default-edge-payload: default-edge-payload,
  default-source-payload: default-source-payload,
  default-sink-payload: default-sink-payload,
  default-node-statements: default-node-statements,
  nodes: (),
  edges: (),
)

/// Parse one or more DOT digraphs into graph objects.
///
/// Default payloads are applied before `eval-*` fields are evaluated, so
/// default payload strings can refer to parsed record fields such as `#name`.
/// Parsed fields take precedence over default payload fields, so DOT
/// `label="..."` overrides `default-node-payload: (label: ...)`.
/// ````example
/// #let a = ```dot
/// digraph {
/// ext [style=invis]
/// a [id=0 pos="0,0!"]
/// b [id=1 pos="ref(node:0)+4,0!"]
/// a -> b [id=0 pos="ref(node:1)+0,1!"]
/// b -> c [id=1 pos="ref(edge:0)+1,0!"]
/// c [id=2 pos="x:2!" label="$c_nu$"]
/// d [id=3 pos="x:2!,y:0.1"]
/// ext -> a [id=2 pos="x:@-left!,y:@edge0!"]
/// }
/// ```
/// #let g = parse(a.text,default-node-payload:(label:"#name"),eval-node-fields:"label")
/// >>>#align(center+horizon, draw(layout(g.at(0))))
///
/// ````
/// -> array
#let parse(
  /// DOT source text containing one or more `digraph` definitions. -> string | bytes
  input,
  /// Default payload merged into every node payload. Captured node payload fields override it.
  ///
  /// ````example
  /// #let a = ```dot
  /// digraph {
  /// node [particle="g"]
  /// a -> b -> c -> d -> a
  /// a -> a
  /// b -> d
  /// c [particle="q"]
  /// }
  /// ```
  /// #let g = parse(a.text,default-node-payload:(particle:"g"),eval-node-fields:"particle").at(0)
  /// #nodes(g).map(n=>n)
  /// ````
  ///   -> any
  default-node-payload: none,
  /// Node statement fields to evaluate into `graph.nodes(g).at(i).payload`.
  /// ````example
  /// #let a = ```dot
  /// digraph {
  /// a -> b -> c -> d -> a
  /// a -> a
  /// b -> d
  /// }
  /// ```
  /// #let g = parse(a.text,default-node-payload:(label:"#name"),eval-node-fields:"label")
  /// >>>#align(center+horizon, draw(layout(g.at(0))))
  ///
  /// ````
  /// -> string | array
  eval-node-fields: (),
  /// Default payload merged into every edge payload. Captured edge payload fields override it. -> any
  default-edge-payload: none,
  /// Default payload merged into every source half-edge payload. Captured source payload fields override it. -> any
  default-source-payload: none,
  /// Default payload merged into every sink half-edge payload. Captured sink payload fields override it. -> any
  default-sink-payload: none,
  /// Graph statement fields to evaluate into `graph.info(g).payload`. -> string | array
  eval-graph-fields: (),
  /// Edge statement fields to evaluate into `graph.edges(g).at(i).payload`. -> string | array
  eval-edge-fields: (),
  /// Source half-edge fields to evaluate into `edge.source.payload`. -> string | array
  eval-source-fields: (),
  /// Sink half-edge fields to evaluate into `edge.sink.payload`. -> string | array
  eval-sink-fields: (),
  /// Typst `eval` mode used for string field values. -> string
  eval-mode: "markup",
  /// Additional Typst names available while evaluating field values. -> dictionary
  scope: (:),
) = _impl.parse(
  input,
  default-node-payload: default-node-payload,
  default-edge-payload: default-edge-payload,
  default-source-payload: default-source-payload,
  default-sink-payload: default-sink-payload,
  eval-graph-fields: eval-graph-fields,
  eval-node-fields: eval-node-fields,
  eval-edge-fields: eval-edge-fields,
  eval-source-fields: eval-source-fields,
  eval-sink-fields: eval-sink-fields,
  eval-mode: eval-mode,
  scope: scope,
)

/// Create a graph node item for @build.
///
/// A Typst label is the node name used by @source, @sink, and @pos. The
/// optional numeric `id` fixes the resulting graph node index. Extra named
/// arguments are captured as node payload fields, so
/// `node(<a>, label: [A], color: red)` stores `(label: [A], color: red)`.
/// The default draw style uses `payload.label` as the visible node label when
/// present.
/// -> array
#let node(
  /// Optional positional node name. Must be a Typst label when provided; extra named arguments become payload fields. -> label
  ..args,
  /// Typst node name for references. -> none | label
  name: none,
  /// Numeric graph node index. Must be unique and in bounds when provided. -> none | int
  id: none,
  /// Node placement. -> none | dictionary
  pos: none,
  /// Drawing shift stored as a statement. -> none | string | array | dictionary
  shift: none,
  /// Additional flat node statements. Used by DOT; values cannot nest. -> dictionary
  statements: (:),
) = _impl.node(
  ..args,
  name: name,
  id: id,
  pos: pos,
  shift: shift,
  statements: statements,
)

/// Create a source half-edge endpoint.
///
/// `node` may be a node name like `<a>` or a numeric node index. `name` gives
/// the half-edge a Typst name; `id` is a numeric half-edge order/index override.
/// Extra named arguments are captured as source payload fields.
/// -> dictionary
#let source(
  /// Referenced node, either by Typst label name or numeric node index. -> label | int
  node,
  /// Extra named arguments become source payload fields. -> any
  ..args,
  /// Optional half-edge name used for later references. -> none | label
  name: none,
  /// Optional numeric half-edge index/order override. -> none | int
  id: none,
  /// DOT-ish flat statement used for matching/joining dangling half edges. -> none | string
  statement: none,
  /// DOT compass point such as `"n"`, `"s"`, `"e"`, or `"w"`. -> none | string
  compass: none,
) = {
  _impl.source(node, ..args, name: name, id: id, statement: statement, compass: compass)
}

/// Create a sink half-edge endpoint.
///
/// `node` may be a node name like `<a>` or a numeric node index. `name` gives
/// the half-edge a Typst name; `id` is a numeric half-edge order/index override.
/// Extra named arguments are captured as sink payload fields.
/// -> dictionary
#let sink(
  /// Referenced node, either by Typst label name or numeric node index. -> label | int
  node,
  /// Extra named arguments become sink payload fields. -> any
  ..args,
  /// Optional half-edge name used for later references. -> none | label
  name: none,
  /// Optional numeric half-edge index/order override. -> none | int
  id: none,
  /// DOT-ish flat statement used for matching/joining dangling half edges. -> none | string
  statement: none,
  /// DOT compass point such as `"n"`, `"s"`, `"e"`, or `"w"`. -> none | string
  compass: none,
) = {
  _impl.sink(node, ..args, name: name, id: id, statement: statement, compass: compass)
}

/// Create a graph edge item for @build.
///
/// Positional arguments may contain one @source, one @sink, and optionally one
/// Typst label used as the edge name. The numeric `id` chooses the edge order.
/// Extra named arguments are captured as edge payload fields, so
/// `edge(source(<a>), sink(<b>), particle: "g")` stores
/// `(particle: "g")`. The default draw style uses `payload.label` as the
/// visible edge label when present.
/// -> array
#let edge(
  /// Source/sink half-edges and optional edge name; extra named arguments become payload fields. -> any
  ..args,
  /// Typst edge name. -> none | label
  name: none,
  /// Numeric edge order/index override. -> none | int
  id: none,
  /// Edge orientation: `"default"`, `"reversed"`, or `"undirected"`. -> string
  orientation: "default",
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
  /// Additional flat edge statements. Used by DOT; values cannot nest. -> dictionary
  statements: (:),
) = _impl.edge(
  ..args,
  name: name,
  id: id,
  orientation: orientation,
  pos: pos,
  shift: shift,
  label-pos: label-pos,
  label-angle: label-angle,
  bend: bend,
  statements: statements,
)

/// Create a grouped placement coordinate.
/// -> dictionary
#let group(
  /// Group identifier shared by positions constrained to the same coordinate. -> string | int | bool
  name,
  /// Optional sign constraint: `"+"`, `"-"`, `"positive"`, or `"negative"`. -> none | string
  side: none,
) = _impl.group(name, side: side)

/// Create a first-class graph placement.
/// -> dictionary
#let pos(
  /// Absolute or grouped x coordinate. -> none | int | float | dictionary
  x: none,
  /// Absolute or grouped y coordinate. -> none | int | float | dictionary
  y: none,
  /// Node reference for relative placement, by name or numeric index. -> none | label | int
  ref: none,
  /// Relative x offset from `ref`. -> none | int | float
  dx: none,
  /// Relative y offset from `ref`. -> none | int | float
  dy: none,
  /// Placement mode: `"pin"` constrains layout, `"start"` only seeds it. -> string
  mode: "pin",
) = {
  _impl.pos(x: x, y: y, ref: ref, dx: dx, dy: dy, mode: mode)
}

/// Map graph metadata to new opaque payloads.
/// -> bytes
#let map(
  /// Graph object to transform. -> bytes
  graph_,
  /// Callback for graph metadata records. -> none | function
  graph: none,
  /// Callback for node records. -> none | function
  node: none,
  /// Callback for edge records. -> none | function
  edge: none,
  /// Callback for source half-edge records. -> none | function
  source: none,
  /// Callback for sink half-edge records. -> none | function
  sink: none,
) = {
  _impl.map(graph_, graph: graph, node: node, edge: edge, source: source, sink: sink)
}

/// Evaluate selected fields into opaque payload entries.
/// -> bytes
#let eval-fields(
  /// Graph object whose selected statement fields should be evaluated. -> bytes
  graph_,
  /// Graph statement fields to evaluate into `graph.info(g).payload`. -> string | array
  eval-graph-fields: (),
  /// Node statement fields to evaluate into node payloads. -> string | array
  eval-node-fields: (),
  /// Edge statement fields to evaluate into edge payloads. -> string | array
  eval-edge-fields: (),
  /// Source half-edge fields to evaluate into source payloads. -> string | array
  eval-source-fields: (),
  /// Sink half-edge fields to evaluate into sink payloads. -> string | array
  eval-sink-fields: (),
  /// Typst `eval` mode used for string field values. -> string
  eval-mode: "markup",
  /// Additional Typst names available while evaluating field values. -> dictionary
  scope: (:),
) = _impl.eval-fields(
  graph_,
  eval-graph-fields: eval-graph-fields,
  eval-node-fields: eval-node-fields,
  eval-edge-fields: eval-edge-fields,
  eval-source-fields: eval-source-fields,
  eval-sink-fields: eval-sink-fields,
  eval-mode: eval-mode,
  scope: scope,
)

/// Return graph metadata.
/// -> dictionary
#let info(
  /// Graph object returned by @build, @parse, #api-link("layout-", "layout"), or another graph API. -> bytes
  graph,
) = _impl.info(graph)

/// Serialize a graph object to DOT.
/// -> string
#let dot(
  /// Graph object to serialize. -> bytes
  graph,
) = _impl.dot(graph)

/// Return node records, optionally filtered by a subgraph object.
/// -> array
#let nodes(
  /// Graph object to inspect. -> bytes
  graph,
  /// Optional subgraph filter; only nodes incident to selected half edges are returned. -> none | bytes
  subgraph: none,
) = _impl.nodes(graph, subgraph: subgraph)

/// Return edge records, optionally filtered by a subgraph object.
/// -> array
#let edges(
  /// Graph object to inspect. -> bytes
  graph,
  /// Optional subgraph filter; only selected edges/half-edges are returned. -> none | bytes
  subgraph: none,
) = _impl.edges(graph, subgraph: subgraph)

/// Return one named node's opaque payload.
/// -> any
#let node-data(
  /// Graph object to inspect. -> bytes
  graph_,
  /// Node name as a Typst label or its string form. -> label | string
  name,
) = _impl.node-data(graph_, name)

/// Return one named edge's opaque payload.
/// -> any
#let edge-data(
  /// Graph object to inspect. -> bytes
  graph_,
  /// Edge name as a Typst label or its string form. -> label | string
  name,
) = _impl.edge-data(graph_, name)

/// Update one named node's opaque payload.
/// -> bytes
#let update-node-data(
  /// Graph object to update. -> bytes
  graph_,
  /// Node name as a Typst label or its string form. -> label | string
  name,
  /// Replacement payload or `(data, node) => new-data` callback. -> any | function
  update,
) = _impl.update-node-data(graph_, name, update)

/// Update one named edge's opaque payload.
/// -> bytes
#let update-edge-data(
  /// Graph object to update. -> bytes
  graph_,
  /// Edge name as a Typst label or its string form. -> label | string
  name,
  /// Replacement payload or `(data, edge) => new-data` callback. -> any | function
  update,
) = _impl.update-edge-data(graph_, name, update)

/// Join two graphs by matching dangling half-edge data on `key`.
/// -> bytes
#let join(
  /// Left graph object. -> bytes
  left,
  /// Right graph object. -> bytes
  right,
  /// Dangling half-edge match key: `"statement"`, `"compass"`, or `"id"`. -> string
  key: "statement",
) = _impl.join(left, right, key: key)

/// Return subgraph objects for the graph's cycle basis.
/// -> array
#let cycles(
  /// Graph object to analyze. -> bytes
  graph,
) = _impl.cycles(graph)

/// Return subgraph objects for the graph's spanning forests.
/// -> array
#let forests(
  /// Graph object to analyze. -> bytes
  graph,
) = _impl.forests(graph)

/// Decode an opaque payload returned by low-level graph APIs.
/// -> any
#let decode-payload(
  /// Payload bytes or byte array produced by the graph API. -> bytes | array | none
  value,
) = _impl.decode-payload(value)
