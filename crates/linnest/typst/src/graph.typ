// Public graph API.
//
// Keep user-facing constructors, queries, and modifiers here. The
// serialization-heavy implementation lives in `impl/graph.typ`.

#import "impl/graph.typ" as _impl

#let graph-bytes(graph) = _impl.graph-bytes(graph)
#let with-bytes(graph, graph-bytes) = _impl.with-bytes(graph, graph-bytes)

/// Build one graph object from a stream of node and edge items.
///
/// Use @node, @source, @sink, and @edge to create graph items. Positional items
/// may be passed as comma-separated arguments or yielded from a Typst code
/// block.
/// -> dictionary
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
  /// Native Typst graph data.
  ///
  /// ```example
  /// #let g = build(data: (a:(b:(1, ))), {
  /// })
  /// #info(g).data
  /// ```
  /// -> any
  data: none,
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
  /// Default data merged into every node data. Captured node data fields override it.
  /// ```example
  /// #let g = build(default-node-data: (a:1), {
  ///    node(<a>)
  ///    node(<b>,a:2)
  ///    node(<c>,b:(a:1))
  /// })
  /// #nodes(g).map(n=>n.data)
  /// ```
  ///  -> any
  default-node-data: none,
  /// Default data merged into every edge data. Captured edge data fields override it.
  /// ```example
  /// #let g = build(default-edge-data: (a:1,b:[$m_mu$]), {
  ///    node(<a>,id:0)
  ///    edge(source(<a>),id:1)
  ///    edge(sink(<a>),<e>,id:0,b:[#set text(font:"Reforma")
  ///    This is a test: ])
  /// })
  /// #edges(g).map(e=>e.data.b).join()
  /// ```
  /// -> any
  default-edge-data: none,
  /// Default data merged into every source half-edge data. Captured source data fields override it.
  /// ```example
  /// #let g = build(default-source-data: (a:1), {
  ///    node(<a>)
  ///    edge(source(<a>))
  ///    edge(sink(<a>))
  /// })
  /// #edges(g).map(e=>if e.source != none {e.source.data} else {none})
  /// ```
  /// -> any
  default-source-data: none,
  /// Default data merged into every sink half-edge data. Captured sink data fields override it.
  /// ```example
  /// #let g = build(default-sink-data: (a:1), {
  ///    node(<a>)
  ///    edge(source(<a>))
  ///    edge(sink(<a>))
  /// })
  /// #edges(g).map(e=>if e.sink != none {e.sink.data} else {none})
  /// ```
  /// -> any
  default-sink-data: none,
) = _impl.build(
  (
    name: name,
    data: data,
    statements: statements,
    default-edge-statements: default-edge-statements,
    default-node-data: default-node-data,
    default-edge-data: default-edge-data,
    default-source-data: default-source-data,
    default-sink-data: default-sink-data,
    default-node-statements: default-node-statements,
    nodes: (),
    edges: (),
  ),
  ..items,
)

/// Parse one or more DOT graphs into graph objects.
///
/// Default data are applied before `eval-*` fields are evaluated, so
/// default data strings can refer to parsed record fields such as `#str(name)`.
/// Parsed fields take precedence over default data fields, so DOT
/// `label="..."` overrides `default-node-data: (label: ...)`.
///
/// Linnest uses `dot-parser` for DOT syntax and then gives special meaning to a
/// small set of attributes. Every other attribute is preserved as a flat string
/// statement. Preserved statements are visible through @info, @nodes, @edges,
/// @map, and @eval-fields, but they do not become native Typst `data` unless a
/// matching `eval-*` argument is passed.
///
/// ````example
/// #let src = ```dot
/// digraph { a [label="A"]; }
/// ```
/// #let n = nodes(parse(src.text).first()).first()
/// #n.statements.label
/// #n.data
/// ````
///
/// Graph-level handling:
///
/// - The DOT graph name becomes `graph.info(g).name`. A graph-level `name`
///   attribute can name an anonymous graph; a named `graph`/`digraph` header
///   takes precedence.
///
/// ````example
/// #let src = ```dot
/// digraph header_wins { graph [name="ignored"]; a }
/// ```
/// #info(parse(src.text).first()).name
/// ````
///
/// - Top-level `key=value` statements and `graph [key=value]` attributes become
///   `graph.info(g).global-statements`, except for the consumed `name`. They do
///   not configure layout, even when a key looks like a layout option; pass
///   layout options directly to #api-link("layout-", "layout").
///
/// ````example
/// #let src = ```dot
/// digraph { steps=1; graph [subtitle="nested"]; a }
/// ```
/// #let statements = info(parse(src.text).first()).global-statements
/// #statements.steps
/// #statements.subtitle
/// ````
///
/// - `node [key=value]` and `edge [key=value]` become
///   `default-node-statements` and `default-edge-statements`. They are merged
///   into each parsed node or edge before local attributes, so local attributes
///   override defaults.
///
/// ````example
/// #let src = ```dot
/// digraph { node [color=gray]; edge [particle=g]; a [color=red]; a -> b }
/// ```
/// #let g = parse(src.text).first()
/// #nodes(g).first().statements.color
/// #edges(g).first().statements.particle
/// ````
///
/// Node handling:
///
/// - The DOT node id becomes the node `name` returned by @nodes. If the DOT node
///   id is numeric, it is consumed as the node index instead and the public
///   `name` is `none`.
///
/// ````example
/// #let src = ```dot
/// digraph { alpha; 1; }
/// ```
/// #nodes(parse(src.text).first()).map(n => n.name)
/// ````
///
/// - `id=<n>` is also consumed as an explicit node index. Explicit node indexes
///   must be unique and in bounds after parsing.
///
/// ````example
/// #let src = ```dot
/// digraph { a [id=1]; b [id=0]; }
/// ```
/// #nodes(parse(src.text).first()).map(n => n.name)
/// ````
///
/// - `style=invis` marks the DOT node as a dangling external endpoint. It is
///   not returned by @nodes. Its remaining attributes are copied onto the
///   external edge that touches it, except `shape` and `label`; `style` and
///   `id` are consumed.
///
/// ````example
/// #let src = ```dot
/// digraph { ext [style=invis, column=left, label=skip]; ext -> a [id=0]; }
/// ```
/// #let g = parse(src.text).first()
/// #nodes(g).map(n => n.name)
/// #edges(g).first().statements.column
/// #edges(g).first().statements.at("label", default: none)
/// ````
///
/// - Node `style` is consumed for every node; only `style=invis` has special
///   behavior. Use another attribute name for drawing style metadata that must
///   survive parsing.
///
/// ````example
/// #let src = ```dot
/// digraph { a [style=filled, "draw-style"=filled]; }
/// ```
/// #nodes(parse(src.text).first()).first().statements
/// ````
///
/// - `pos` is parsed as node placement. Simple numeric values are exposed as
///   `node.pos`; extended placement values are used by
///   #api-link("layout-", "layout"). `pin` is an explicit layout constraint and
///   is not exposed as a public statement.
///
/// ````example
/// #let src = ```dot
/// digraph { a [id=0, pos="0,0!"]; b [id=1, pos="ref(node:0)+2,0!"]; }
/// ```
/// #let parsed-nodes = nodes(parse(src.text).first())
/// #parsed-nodes.map(n => n.pos)
/// #parsed-nodes.map(n => n.statements.at("pin", default: none))
/// ````
///
/// - `shift` is parsed as a drawing shift and also remains available as a
///   statement. `eval` is retained for legacy round-tripping. Other node
///   attributes are preserved as node statements.
///
/// ````example
/// #let src = ```dot
/// digraph { a [shift="0.1,0", eval=legacy, label=A]; }
/// ```
/// #let n = nodes(parse(src.text).first()).first()
/// #n.shift
/// #n.statements
/// ````
///
/// Edge handling:
///
/// - `id=<n>` is consumed as the edge index. It is not preserved as a statement;
///   drawing callbacks expose the resulting stable edge index as `eid`.
///   Explicit edge indexes must be unique and in bounds. Without explicit ids,
///   edge order is parser-internal, not DOT input order.
///
/// ````example
/// #let src = ```dot
/// digraph { a -> b [id=1, label=later]; b -> c [id=0, label=first]; }
/// ```
/// #let parsed-edges = edges(parse(src.text).first())
/// #parsed-edges.map(e => e.edge)
/// #parsed-edges.map(e => e.statements.label)
/// #parsed-edges.map(e => e.statements.at("id", default: none))
/// ````
///
/// - `dir=forward`, `dir=back`, and `dir=none` become edge orientations
///   `"default"`, `"reversed"`, and `"undirected"`. If omitted, directed DOT
///   edges are `"default"` and undirected DOT edges are `"undirected"`.
///
/// ````example
/// #let src = ```dot
/// digraph { a -> b [id=0]; b -> c [id=1, dir=back]; c -> d [id=2, dir=none]; }
/// ```
/// #edges(parse(src.text).first()).map(e => e.orientation)
/// ````
///
/// - `source="..."` and `sink="..."` are consumed as endpoint `statement`
///   values and removed from edge statements. On an external edge, use the
///   attribute matching the real endpoint side in the DOT edge: `ext -> a`
///   uses `sink=...`, while `a -> ext` uses `source=...`.
///
/// ````example
/// #let src = ```dot
/// digraph { ext [style=invis]; ext -> a [id=0, sink=in]; a -> ext [id=1, source=out]; }
/// ```
/// #let endpoint-statement(endpoint) = if endpoint == none { none } else { endpoint.at("statement", default: none) }
/// #let parsed-edges = edges(parse(src.text).first())
/// #parsed-edges.map(e => endpoint-statement(e.source))
/// #parsed-edges.map(e => endpoint-statement(e.sink))
/// #parsed-edges.map(e => e.statements.at("source", default: none))
/// ````
///
/// - `pos` is parsed as the edge control-point placement. `pin` is an explicit
///   edge layout constraint. `shift`, `label-pos`, `label-angle`, and `bend`
///   are parsed into edge geometry fields; except for `pin`, they also remain
///   available as statements. Other edge attributes are preserved as edge
///   statements.
///
/// ````example
/// #let src = ```dot
/// digraph { a -> b [id=0, pos="0,1!", pin="x:@edge", shift="0.1,0", "label-pos"="0,1.2", "label-angle"="0.3rad", bend="0.4rad", particle=g]; }
/// ```
/// #let e = edges(parse(src.text).first()).first()
/// #e.pos
/// #e.shift
/// #e.label-pos
/// #e.label-angle
/// #e.bend
/// #e.statements.particle
/// #e.statements.at("pin", default: none)
/// ````
///
/// - Attributes copied from an invisible endpoint node are merged into the
///   external edge statements unless their keys are `shape` or `label`.
///
/// ````example
/// #let src = ```dot
/// digraph { ext [style=invis, column=left, shape=none, label=skip]; ext -> a [id=0]; }
/// ```
/// #edges(parse(src.text).first()).first().statements
/// ````
///
/// Port and half-edge handling:
///
/// - DOT ports of the form `node:port` and `node:port:compass` are preserved on
///   source/sink endpoint records. Numeric ports also assign explicit half-edge
///   ids. Non-numeric ports are exposed only as `port-label`.
///
/// ````example
/// #let src = ```dot
/// digraph { a:left:e -> b:0:w [id=0]; }
/// ```
/// #let e = edges(parse(src.text).first()).first()
/// #e.source.port-label
/// #e.sink.hedge
/// ````
///
/// - Compass values are exposed as `compass`; valid DOT compass names include
///   `n`, `ne`, `e`, `se`, `s`, `sw`, `w`, `nw`, `c`, and `_`.
///
/// ````example
/// #let src = ```dot
/// digraph { a:n -> b:sw [id=0]; }
/// ```
/// #let e = edges(parse(src.text).first()).first()
/// #e.source.compass
/// #e.sink.compass
/// ````
///
/// - Explicit half-edge ids from numeric ports must be unique and in bounds.
///   They are retained for half-edge ordering and for @join with `key: "id"`;
///   query and eval records expose the resulting half-edge index as `hedge`.
///
/// ````example
/// #let src = ```dot
/// digraph { a:1 -> b:0 [id=0]; }
/// ```
/// #let e = edges(parse(src.text).first()).first()
/// #e.source.hedge
/// #e.sink.hedge
/// ````
///
/// Placement fields:
///
/// - `pos="x,y"` gives a starting coordinate; `pos="x,y!"` pins both axes.
///
/// ````example
/// #let src = ```dot
/// digraph { a [id=0, pos="0,0"]; b [id=1, pos="2,0!"]; }
/// ```
/// #let parsed-nodes = nodes(parse(src.text).first())
/// #parsed-nodes.map(n => n.pos)
/// #parsed-nodes.map(n => n.statements.at("pos-mode", default: none))
/// ````
///
/// - `pos="ref(node:<id>)+dx,dy!"` and `pos="ref(edge:<id>)+dx,dy!"` reference
///   explicit DOT node or edge ids, not names and not implicit parse order.
///
/// ````example
/// #let src = ```dot
/// digraph { a [id=0, pos="1,1!"]; b [id=1, pos="ref(node:0)+2,0!"]; a -> b [id=0, pos="ref(node:1)+0,1!"]; }
/// ```
/// #let g = parse(src.text).first()
/// #nodes(g).at(1).pos
/// #edges(g).first().pos
/// ````
///
/// - Axis form accepts `x:<coord>` and `y:<coord>` entries. `!` applies to the
///   individual axis, for example `pos="x:2!,y:0"`.
///
/// ````example
/// #let src = ```dot
/// digraph { a [id=0, pos="x:2!,y:0"]; }
/// ```
/// #nodes(parse(src.text).first()).first().pos
/// ````
///
/// - Group coordinates use `@name`, `@+name`, or `@-name` in axis form and must
///   be pinned with `!`, for example `pos="x:@-left!,y:@row!"`.
///
/// ````example
/// #let src = ```dot
/// digraph { ext [style=invis]; ext -> a [id=0, pos="x:@-left!,y:@row!"]; }
/// ```
/// #edges(layout(parse(src.text).first(), layout-algo: "tree")).first().pos
/// ````
///
/// - `pin` accepts numeric point constraints, `x:<coord>`, `y:<coord>`, and
///   grouped constraints such as `x:@left`, `x:@+right`, or `@row`.
///
/// ````example
/// #let src = ```dot
/// digraph { a -> b [id=0, pin="x:@+right"]; }
/// ```
/// #edges(layout(parse(src.text).first(), layout-algo: "tree")).first().statements.at("pin", default: none)
/// ````
/// -> array
#let parse(
  /// DOT source text containing one or more `graph` or `digraph` definitions. -> string | bytes
  input,
  /// Default data merged into every node data. Captured node data fields override it, but bare dot statements don't set data fields. To turn statements into data fields, use `eval-node-fields`.
  ///
  /// ````example
  /// #let a = ```dot
  /// digraph {
  /// a -> b -> c -> d -> a
  /// a -> a
  /// b -> d
  /// c [particle="q"]
  /// }
  /// ```
  /// #let g = parse(a.text,default-node-data:(particle:"g")).at(0)
  /// #nodes(g).map(n=>n.data.particle)
  /// ````
  ///   -> any
  default-node-data: none,
  /// Node fields to evaluate into node data. The eval scope includes node
  /// statements plus `node`, `name`, and `pos`.
  ///
  /// ````example
  /// #let src = ```dot
  /// digraph { a [label="#str(name)"]; }
  /// ```
  /// #nodes(parse(src.text, eval-node-fields: "label").first()).first().data.label
  /// ````
  /// -> string | array
  eval-node-fields: (),
  /// Default data merged into every edge data. Captured edge data fields override it, but bare dot statements don't set data fields. To turn statements into data fields, use `eval-edge-fields`.
  /// 
  /// -> any
  default-edge-data: none,
  /// Edge statement fields to evaluate into `graph.edges(g).at(i).data`. The eval
  /// scope includes edge statements plus `edge`, `orientation`, endpoint records,
  /// placement fields, and existing edge data. Drawing callbacks later expose the
  /// edge index as `eid`.
  ///
  /// ````example
  /// #let src = ```dot
  /// digraph { a -> b [id=0, label="#orientation"]; }
  /// ```
  /// #edges(parse(src.text, eval-edge-fields: "label").first()).first().data.label
  /// ````
  /// -> string | array
  eval-edge-fields: (),
  /// Default data merged into every source half-edge data. Captured source data fields override it. -> any
  default-source-data: none,
  /// Source half-edge fields to evaluate into `edge.source.data`. The eval scope
  /// includes source endpoint fields `statement`, `port-label`, `compass`, `node`,
  /// and `hedge`, plus the surrounding edge fields.
  ///
  /// ````example
  /// #let src = ```dot
  /// digraph { a:left:e -> b [id=0, source="#port-label"]; }
  /// ```
  /// #edges(parse(src.text, eval-source-fields: "statement").first()).first().source.data.statement
  /// ````
  /// -> string | array
  eval-source-fields: (),
  /// Default data merged into every sink half-edge data. Captured sink data fields override it. -> any
  default-sink-data: none,
  /// Sink half-edge fields to evaluate into `edge.sink.data`. The eval scope
  /// includes sink endpoint fields `statement`, `port-label`, `compass`, `node`,
  /// and `hedge`, plus the surrounding edge fields.
  ///
  /// ````example
  /// #let src = ```dot
  /// digraph { a -> b:0:w [id=0, sink="#str(hedge)"]; }
  /// ```
  /// #edges(parse(src.text, eval-sink-fields: "statement").first()).first().sink.data.statement
  /// ````
  /// -> string | array
  eval-sink-fields: (),
  /// Graph statement fields to evaluate into `graph.info(g).data`. The eval scope
  /// includes graph statements and direct graph record fields.
  ///
  /// ````example
  /// #let src = ```dot
  /// digraph {
  ///   title = "Strong graph";
  /// }
  /// ```
  /// #info(parse(src.text, eval-graph-fields: "title").first()).data.title
  /// ````
  /// -> string | array
  eval-graph-fields: (),
  /// Typst `eval` mode used for string field values. -> string
  eval-mode: "markup",
  /// Additional Typst names available while evaluating field values. -> dictionary
  scope: (:),
) = _impl.parse(input, (
  default-node-data: default-node-data,
  default-edge-data: default-edge-data,
  default-source-data: default-source-data,
  default-sink-data: default-sink-data,
  eval-graph-fields: eval-graph-fields,
  eval-node-fields: eval-node-fields,
  eval-edge-fields: eval-edge-fields,
  eval-source-fields: eval-source-fields,
  eval-sink-fields: eval-sink-fields,
  eval-mode: eval-mode,
  scope: scope,
))

/// Create a graph node item for @build.
///
/// A Typst label is the node name used by @source, @sink, and @pos. The
/// optional numeric `id` fixes the resulting graph node index. Extra named
/// arguments are captured as node data fields, so
/// `node(<a>, label: [A], color: red)` stores `(label: [A], color: red)`.
/// The default draw style uses `data.label` as the visible node label when
/// present.
/// -> array
#let node(
  /// Optional positional node name. Must be a Typst label when provided; extra named arguments become data fields. -> label
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
  (
    name: name,
    id: id,
    pos: pos,
    shift: shift,
    statements: statements,
  ),
  ..args,
)

/// Create a source half-edge endpoint.
///
/// `node` may be a node name like `<a>` or a numeric node index. `name` gives
/// the half-edge a Typst name; `id` is a numeric half-edge order/index override.
/// Extra named arguments are captured as source data fields.
/// -> dictionary
#let source(
  /// Referenced node, either by Typst label name or numeric node index. -> label | int
  node,
  /// Extra named arguments become source data fields. -> any
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
  _impl.source(node, (name: name, id: id, statement: statement, compass: compass), ..args)
}

/// Create a sink half-edge endpoint.
///
/// `node` may be a node name like `<a>` or a numeric node index. `name` gives
/// the half-edge a Typst name; `id` is a numeric half-edge order/index override.
/// Extra named arguments are captured as sink data fields.
/// -> dictionary
#let sink(
  /// Referenced node, either by Typst label name or numeric node index. -> label | int
  node,
  /// Extra named arguments become sink data fields. -> any
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
  _impl.sink(node, (name: name, id: id, statement: statement, compass: compass), ..args)
}

/// Create a graph edge item for @build.
///
/// Positional arguments may contain one @source, one @sink, and optionally one
/// Typst label used as the edge name. The numeric `id` chooses the edge order.
/// Extra named arguments are captured as edge data fields, so
/// `edge(source(<a>), sink(<b>), particle: "g")` stores
/// `(particle: "g")`. The default draw style uses `data.label` as the
/// visible edge label when present.
/// -> array
#let edge(
  /// Source/sink half-edges and optional edge name; extra named arguments become data fields. -> any
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
  (
    name: name,
    id: id,
    orientation: orientation,
    pos: pos,
    shift: shift,
    label-pos: label-pos,
    label-angle: label-angle,
    bend: bend,
    statements: statements,
  ),
  ..args,
)

/// Create a grouped placement coordinate.
///
/// `side: "+"` keeps the solved coordinate non-negative and `side: "-"`
/// keeps it non-positive. Groups are layout constraints and therefore require
/// pin placement, which is the @pos default.
///
/// ```example
/// #group("right", side: "+")
/// ```
/// -> dictionary
#let group(
  /// Group identifier shared by positions constrained to the same coordinate. -> string | int | bool
  name,
  /// Optional sign constraint: `"+"`, `"-"`, `"positive"`, or `"negative"`. -> none | string
  side: none,
) = _impl.group(name, side)

/// Create a first-class graph placement.
///
/// The default `mode: "pin"` turns numeric and grouped coordinates into layout
/// constraints and also makes the coordinates immediately drawable without a
/// layout pass. Use `mode: "start"` when the coordinate should only seed the
/// layout.
///
/// ```example
/// #pos(x: 0, y: group("row"), mode: "pin")
/// ```
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
  _impl.pos((x: x, y: y, ref: ref, dx: dx, dy: dy, mode: mode))
}

/// Map graph metadata to new native data.
///
/// The callbacks receive decoded records plus a `fields` dictionary containing
/// merged statements and direct record fields. A callback returns `none` to
/// leave the record unchanged, or `(data: value)` to set new native data.
///
/// ```example
/// #let g = build({ node(<a>) })
/// #let g = map(g, node: node => (data: (label: [A])))
/// #nodes(g).first().data.label
/// ```
/// -> dictionary
#let map(
  /// Graph object to transform. -> dictionary
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
  _impl.map(graph_, (graph: graph, node: node, edge: edge, source: source, sink: sink))
}

/// Evaluate selected fields into native data entries.
///
/// Each selected field is read from the record's merged `fields` dictionary,
/// evaluated in a scope containing those fields, and written to `data.<field>`.
///
/// ```example
/// #let g = build({
///   node(<a>)
///   node(<b>)
///   edge(source(<a>), sink(<b>), statements: (mom: "p"))
/// }, default-edge-statements: (display-label: "$#mom$"))
/// #let g = eval-fields(g, eval-edge-fields: ("display-label",))
/// #edges(g).first().data.at("display-label")
/// ```
/// -> dictionary
#let eval-fields(
  /// Graph object whose selected statement fields should be evaluated. -> dictionary
  graph_,
  /// Graph statement fields to evaluate into `graph.info(g).data`. -> string | array
  eval-graph-fields: (),
  /// Node statement fields to evaluate into node data. -> string | array
  eval-node-fields: (),
  /// Edge statement fields to evaluate into edge data. -> string | array
  eval-edge-fields: (),
  /// Source half-edge fields to evaluate into source data. -> string | array
  eval-source-fields: (),
  /// Sink half-edge fields to evaluate into sink data. -> string | array
  eval-sink-fields: (),
  /// Typst `eval` mode used for string field values. -> string
  eval-mode: "markup",
  /// Additional Typst names available while evaluating field values. -> dictionary
  scope: (:),
) = _impl.eval-fields(graph_, (
  eval-graph-fields: eval-graph-fields,
  eval-node-fields: eval-node-fields,
  eval-edge-fields: eval-edge-fields,
  eval-source-fields: eval-source-fields,
  eval-sink-fields: eval-sink-fields,
  eval-mode: eval-mode,
  scope: scope,
))

/// Return graph metadata.
///
/// The result has `name`, `global-statements`, `default-edge-statements`, and
/// `default-node-statements`.
///
/// ```example
/// #let g = build({ node(<a>) }, name: "demo")
/// #info(g).name
/// ```
/// -> dictionary
#let info(
  /// Graph object returned by @build, @parse, #api-link("layout-", "layout"), or another graph API. -> dictionary
  graph,
) = _impl.info(graph)

/// Serialize a graph object to DOT.
///
/// ```example
/// #let g = build({
///   node(<a>)
///   node(<b>)
///   edge(source(<a>), sink(<b>))
/// }, name: "demo")
/// #dot(g).contains("digraph demo")
/// ```
/// -> string
#let dot(
  /// Graph object to serialize. -> dictionary
  graph,
) = _impl.dot(graph)

/// Return node records, optionally filtered by a subgraph object.
///
/// Node `name` values are Typst labels when present.
///
/// ```example
/// #let g = build({
///   node(<a>)
///   node(<b>)
/// })
/// #nodes(g).map(node => str(node.name)).join(", ")
/// ```
/// -> array
#let nodes(
  /// Graph object to inspect. -> dictionary
  graph,
  /// Optional subgraph filter; only nodes incident to selected half edges are returned. -> none | bytes
  subgraph: none,
) = _impl.nodes(graph, subgraph)

/// Return edge records, optionally filtered by a subgraph object.
///
/// Edge `name` values are Typst labels when present.
///
/// ```example
/// #let g = build({
///   node(<a>)
///   node(<b>)
///   edge(source(<a>, compass: "e"), sink(<b>))
/// })
/// #edges(g).len()
/// ```
/// -> array
#let edges(
  /// Graph object to inspect. -> dictionary
  graph,
  /// Optional subgraph filter; only selected edges/half-edges are returned. -> none | bytes
  subgraph: none,
) = _impl.edges(graph, subgraph)

/// Return one named node's native data.
///
/// `name` is a Typst label such as `<a>` or the corresponding string name.
///
/// ```example
/// #let g = build({ node(<a>, label: [A]) })
/// #node-data(g, <a>).label
/// ```
/// -> any
#let node-data(
  /// Graph object to inspect. -> dictionary
  graph_,
  /// Node name as a Typst label or its string form. -> label | string
  name,
) = _impl.node-data(graph_, name)

/// Return one named edge's native data.
///
/// `name` is a Typst label such as `<e>` or the corresponding string name.
///
/// ```example
/// #let g = build({
///   node(<a>)
///   node(<b>)
///   edge(source(<a>), <e>, sink(<b>), label: [$p$])
/// })
/// #edge-data(g, <e>).label
/// ```
/// -> any
#let edge-data(
  /// Graph object to inspect. -> dictionary
  graph_,
  /// Edge name as a Typst label or its string form. -> label | string
  name,
) = _impl.edge-data(graph_, name)

/// Update one named node's native data.
///
/// `update` may be a replacement data value or a function
/// `(data, node) => new-data`.
///
/// ```example
/// #let g = build({ node(<a>) })
/// #let g = update-node-data(g, <a>, (label: [A]))
/// #nodes(g).first().data.label
/// ```
/// -> dictionary
#let update-node-data(
  /// Graph object to update. -> dictionary
  graph_,
  /// Node name as a Typst label or its string form. -> label | string
  name,
  /// Replacement data or `(data, node) => new-data` callback. -> any | function
  update,
) = _impl.update-node-data(graph_, name, update)

/// Update one named edge's native data.
///
/// `update` may be a replacement data value or a function
/// `(data, edge) => new-data`.
///
/// ```example
/// #let g = build({
///   node(<a>)
///   node(<b>)
///   edge(source(<a>), <e>, sink(<b>))
/// })
/// #let g = update-edge-data(g, <e>, (label: [$p$]))
/// #edges(g).first().data.label
/// ```
/// -> dictionary
#let update-edge-data(
  /// Graph object to update. -> dictionary
  graph_,
  /// Edge name as a Typst label or its string form. -> label | string
  name,
  /// Replacement data or `(data, edge) => new-data` callback. -> any | function
  update,
) = _impl.update-edge-data(graph_, name, update)

/// Join two graphs by matching dangling half-edge statements or ids on `key`.
///
/// Supported key values are `"statement"`, `"compass"`, and `"id"`.
///
/// ```example
/// #let left = build({
///   node(<a>)
///   edge(sink(<a>, statement: "j"))
/// })
/// #let right = build({
///   node(<b>)
///   edge(source(<b>, statement: "j"))
/// })
/// #edges(join(left, right, key: "statement")).len()
/// ```
/// -> dictionary
#let join(
  /// Left graph object. -> dictionary
  left,
  /// Right graph object. -> dictionary
  right,
  /// Dangling half-edge match key: `"statement"`, `"compass"`, or `"id"`. -> string
  key: "statement",
) = _impl.join(left, right, key)

/// Return subgraph objects for the graph's cycle basis.
///
/// ```example
/// #let g = build({
///   node(<a>)
///   node(<b>)
///   edge(source(<a>), sink(<b>))
/// })
/// #cycles(g).len()
/// ```
/// -> array
#let cycles(
  /// Graph object to analyze. -> dictionary
  graph,
) = _impl.cycles(graph)

/// Return subgraph objects for the graph's spanning forests.
///
/// ```example
/// #let g = build({
///   node(<a>)
///   node(<b>)
///   edge(source(<a>), sink(<b>))
/// })
/// #forests(g).len()
/// ```
/// -> array
#let forests(
  /// Graph object to analyze. -> dictionary
  graph,
) = _impl.forests(graph)
