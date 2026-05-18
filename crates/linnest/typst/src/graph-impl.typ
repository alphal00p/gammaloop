// Internal graph implementation. Public users should import `graph.typ`.

#let _plugin = plugin("../linnest.wasm")

#let _encode-data-value(value) = {
  if type(value) == content {
    (linnest-data-kind: "content", value: cbor.encode(value))
  } else if type(value) == array {
    value.map(_encode-data-value)
  } else if type(value) == dictionary {
    let result = (:)
    for key in value.keys() {
      result.insert(key, _encode-data-value(value.at(key)))
    }
    result
  } else {
    value
  }
}

#let _encode-data(value) = if value == none { none } else { cbor.encode(_encode-data-value(value)) }

#let _data-bytes(value) = if type(value) == array { bytes(value) } else { value }

#let _math-text-source(text_) = "\"" + text_.replace("\\", "\\\\").replace("\"", "\\\"") + "\""

#let _sequence-content(children, decode-content) = {
  let result = []
  for child in children {
    result += decode-content(child)
  }
  result
}

#let _math-source(value) = {
  let atom = value => {
    let source = _math-source(value)
    if source == none {
      none
    } else if type(value) == dictionary and value.at("func", default: none) in ("symbol", "text") {
      source
    } else {
      "(" + source + ")"
    }
  }
  if type(value) != dictionary {
    none
  } else {
    let func = value.at("func", default: none)
    if func == "symbol" {
      value.at("text", default: "")
    } else if func == "text" {
      _math-text-source(value.at("text", default: ""))
    } else if func == "space" {
      " "
    } else if func == "sequence" {
      value.children.map(_math-source).join("")
    } else if func == "attach" {
      let base = atom(value.base)
      if base == none {
        return none
      }
      let result = base
      let bottom = value.at("b", default: none)
      if bottom != none {
        let source = atom(bottom)
        if source == none {
          return none
        }
        result += "_" + source
      }
      let top = value.at("t", default: none)
      if top != none {
        let source = atom(top)
        if source == none {
          return none
        }
        result += "^" + source
      }
      result
    } else {
      none
    }
  }
}

#let _decode-content(value, decode) = {
  let func = value.at("func", default: none)
  if func == "text" {
    text(value.text)
  } else if func == "space" {
    [ ]
  } else if func == "sequence" {
    _sequence-content(value.children, child => _decode-content(child, decode))
  } else if func == "strong" {
    strong(_decode-content(value.body, decode))
  } else if func == "emph" {
    emph(_decode-content(value.body, decode))
  } else if func == "equation" {
    let source = _math-source(value.body)
    if source == none {
      text(repr(value))
    } else {
      eval("$" + source + "$", mode: "markup")
    }
  } else {
    text(repr(value))
  }
}

#let _decode-data-value(value) = {
  if type(value) == array {
    value.map(_decode-data-value)
  } else if type(value) == dictionary {
    let kind = value.at("linnest-data-kind", default: none)
    if kind == "content" {
      _decode-content(cbor(_data-bytes(value.value)), _decode-data-value)
    } else {
      let result = (:)
      for key in value.keys() {
        result.insert(key, _decode-data-value(value.at(key)))
      }
      result
    }
  } else {
    value
  }
}

#let decode-data(value) = if value == none { none } else { _decode-data-value(cbor(_data-bytes(value))) }

#let _merge-data(default, data) = {
  if default == none {
    data
  } else if data == none {
    default
  } else if type(default) == dictionary and type(data) == dictionary {
    default + data
  } else {
    data
  }
}

#let _decode-data-field(record) = {
  let result = record
  let data = result.at("data", default: none)
  if data != none {
    result.data = decode-data(data)
  }
  result
}

#let _decode-edge-data(edge) = {
  let result = _decode-data-field(edge)
  for side in ("source", "sink") {
    let endpoint = result.at(side, default: none)
    if endpoint != none {
      result.insert(side, _decode-data-field(endpoint))
    }
  }
  result
}

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

#let _statements-with-value(statements, key, value, context_) = {
  let result = statements
  if value != none {
    result.insert(key, _statement-value(value, context_ + "." + key))
  }
  result
}

#let _statement-value(value, context_) = {
  if type(value) == str {
    value
  } else if type(value) == int or type(value) == float or type(value) == bool {
    str(value)
  } else {
    panic(context_ + ": statement values must be flat scalars; use data for nested data or Typst content")
  }
}

#let _flat-statements(statements, context_) = {
  if type(statements) != dictionary {
    panic(context_ + ": statements must be a flat dictionary")
  }
  let result = (:)
  for key in statements.keys() {
    let value = statements.at(key)
    if value != none {
      result.insert(key, _statement-value(value, context_ + "." + key))
    }
  }
  result
}

#let _kind(item) = if type(item) == dictionary { item.at("linnest-kind", default: none) } else { none }

#let _is-label(value) = type(value) == label

#let _data-from-args(context_, args) = {
  let named = args.named()
  if named.keys().contains("data") {
    panic(context_ + ": use direct named data fields instead of data: (...)")
  }
  if named.len() == 0 {
    none
  } else {
    named
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
  } else if value != none and value < 0 {
    panic(context_ + ": id must be non-negative")
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

#let _check-edge-names(edges) = {
  let seen = ()
  for item in edges {
    let name = item.at("name", default: none)
    if type(name) == label {
      let key = _label-key(name, "graph.edge")
      if seen.contains(key) {
        panic("graph.edge: duplicate edge name <" + key + ">")
      }
      seen.push(key)
    }
  }
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

#let _node-spec(node, default-data: none) = {
  let statements = _statements-with-point(_flat-statements(node.at("statements", default: (:)), "graph.node statements"), "shift", node.at("shift", default: none))
  let data = _encode-data(_merge-data(default-data, node.at("data", default: none)))
  let clean = node
  let name = clean.at("name", default: none)
  for key in ("linnest-kind", "key", "id", "shift", "name", "data") {
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
  if data != none {
    clean.insert("data", data)
  }
  clean + (statements: statements)
}

#let _resolved-node-spec(node, node-keys, default-data: none) = {
  let index = node.at("id", default: none)
  let item = node + (pos: _resolve-pos(node.at("pos", default: none), node-keys))
  if index != none {
    item.insert("index", index)
  }
  _node-spec(item, default-data: default-data)
}

#let _edge-statements(edge) = {
  let statements = _statements-with-value(
    _statements-with-value(
      _statements-with-point(
        _statements-with-point(
          _flat-statements(edge.at("statements", default: (:)), "graph.edge statements"),
          "shift",
          edge.at("shift", default: none),
        ),
        "label-pos",
        edge.at("label-pos", default: none),
      ),
      "label-angle",
      edge.at("label-angle", default: none),
      "graph.edge statements",
    ),
    "bend",
    edge.at("bend", default: none),
    "graph.edge statements",
  )
  statements
}

#let _edge-spec(edge, default-data: none) = {
  let statements = _edge-statements(edge)
  let name = edge.at("name", default: none)
  if name != none {
    statements.insert("__linnest-edge-name", _label-key(name, "graph.edge"))
  }
  let data = _encode-data(_merge-data(default-data, edge.at("data", default: none)))
  let clean = edge
  for key in (
    "linnest-kind",
    "shift",
    "label-pos",
    "label-angle",
    "bend",
    "data",
    "name",
  ) {
    if clean.keys().contains(key) {
      let _ = clean.remove(key)
    }
  }
  if data != none {
    clean.insert("data", data)
  }
  clean + (statements: statements)
}

#let _resolved-half-spec(half, node-keys, context_, default-data: none) = {
  if half == none {
    return none
  }
  let result = (
    node: _resolve-node-ref(half.node, node-keys, context_),
  )
  let id = half.at("id", default: none)
  if id != none {
    result.insert("id", id)
  }
  let data = _encode-data(_merge-data(default-data, half.at("data", default: none)))
  if data != none {
    result.insert("data", data)
  }
  if half.at("statement", default: none) != none {
    result.insert("statement", _statement-value(half.statement, context_ + ".statement"))
  }
  if half.at("compass", default: none) != none {
    result.insert("compass", _statement-value(half.compass, context_ + ".compass"))
  }
  result
}

#let _resolved-edge-spec(
  edge,
  node-keys,
  default-edge-data: none,
  default-source-data: none,
  default-sink-data: none,
) = {
  let item = edge + (
    source: _resolved-half-spec(
      edge.at("source", default: none),
      node-keys,
      "graph.source",
      default-data: default-source-data,
    ),
    sink: _resolved-half-spec(
      edge.at("sink", default: none),
      node-keys,
      "graph.sink",
      default-data: default-sink-data,
    ),
    pos: _resolve-pos(edge.at("pos", default: none), node-keys),
  )
  _edge-spec(item, default-data: default-edge-data)
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
  (kind: "group", name: _statement-value(name, "graph.group name"), side: side)
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

#let _eval-field-list(fields, context_) = {
  if fields == none {
    ()
  } else if type(fields) == str {
    (fields,)
  } else if type(fields) == array {
    for field in fields {
      if type(field) != str {
        panic(context_ + ": eval field names must be strings")
      }
    }
    fields
  } else {
    panic(context_ + ": eval fields must be a string or an array of strings")
  }
}

#let _eval-dot-value(value, eval-mode, scope) = {
  if type(value) == str {
    eval(value, mode: eval-mode, scope: scope)
  } else {
    value
  }
}

#let _record-fields(record) = {
  let fields = (:)
  let data = record.at("data", default: none)
  if type(data) == dictionary {
    fields += data
  }
  let statements = record.at("statements", default: (:))
  if type(statements) == dictionary {
    fields += statements
  }
  for key in record.keys() {
    if key != "statements" and key != "fields" {
      let value = record.at(key)
      if value != none {
        fields.insert(key, value)
      }
    }
  }
  fields
}

#let _record-with-fields(record, base-fields: (:), extra: (:)) = {
  let fields = base-fields + _record-fields(record)
  record + extra + (fields: fields)
}

#let _data-field-value(record, field) = {
  let fields = record.at("fields", default: _record-fields(record))
  if type(fields) == dictionary and field in fields.keys() {
    fields.at(field)
  } else {
    record.at(field, default: none)
  }
}

#let _record-eval-scope(record, scope) = {
  let fields = record.at("fields", default: _record-fields(record))
  let result = scope + fields
  for key in record.keys() {
    if key != "statements" and key != "fields" {
      let value = record.at(key)
      if value != none {
        result.insert(key, value)
      }
    }
  }
  result + (fields: fields, record: record)
}

#let _eval-data-value(record, fields, eval-mode, scope) = {
  let existing = record.at("data", default: none)
  let data = if type(existing) == dictionary { existing } else { (:) }
  let changed = false
  let local-scope = _record-eval-scope(record, scope)
  for field in fields {
    let value = _data-field-value(record, field)
    if value != none {
      data.insert(field, _eval-dot-value(value, eval-mode, local-scope))
      changed = true
    }
  }
  if not changed {
    none
  } else {
    data
  }
}

#let _mapped-data(callback, record, context_) = {
  if callback == none {
    return none
  }
  let result = callback(record)
  if result == none {
    none
  } else if type(result) != dictionary {
    panic(context_ + ": callback must return none or a dictionary")
  } else {
    result.at("data", default: none)
  }
}

#let _data-patch(callback, record, context_) = {
  let data = _mapped-data(callback, record, context_)
  if data == none {
    none
  } else {
    _encode-data(data)
  }
}

#let _default-data-patch(default, record) = {
  let data = _merge-data(default, record.at("data", default: none))
  if data == none {
    none
  } else {
    (data: data)
  }
}

/// Map graph metadata to new opaque data.
///
/// The callbacks receive decoded records plus a `fields` dictionary containing
/// merged statements and direct record fields. A callback returns `none` to
/// leave the record unchanged, or `(data: value)` to set a new opaque data.
///
/// ```example
/// #let g = graph.build({ graph.node(<a>) })
/// #let g = graph.map(g, node: node => (data: (label: [A])))
/// #graph.nodes(g).first().data.label
/// ```
/// -> bytes
#let map(
  graph_,
  /// Graph metadata callback. -> none | function
  graph: none,
  /// Node metadata callback. -> none | function
  node: none,
  /// Edge metadata callback. -> none | function
  edge: none,
  /// Source half-edge metadata callback. -> none | function
  source: none,
  /// Sink half-edge metadata callback. -> none | function
  sink: none,
) = {
  let patch = (:)
  let info = _decode-data-field(cbor(_plugin.graph_info(bytes(graph_))))
  let graph-record = _record-with-fields(info + (statements: info.at("global-statements", default: (:))))
  let graph-data = _data-patch(graph, graph-record, "graph.map graph")
  if graph-data != none {
    patch.insert("data", graph-data)
  }

  let node-patches = ()
  if node != none {
    for node-record in cbor(_plugin.graph_nodes(bytes(graph_))).map(_decode-data-field) {
      let data = _data-patch(node, _record-with-fields(node-record), "graph.map node")
      if data != none {
        node-patches.push((index: node-record.node, data: data))
      }
    }
  }
  if node-patches.len() > 0 {
    patch.insert("nodes", node-patches)
  }

  let edge-patches = ()
  if edge != none or source != none or sink != none {
    for edge-source in cbor(_plugin.graph_edges(bytes(graph_))).map(_decode-edge-data) {
      let edge-patch = (index: edge-source.edge)
      let edge-record = _record-with-fields(edge-source)
      let edge-fields = edge-record.fields
      let data = _data-patch(edge, edge-record, "graph.map edge")
      if data != none {
        edge-patch.insert("data", data)
      }

      let source-record = edge-source.at("source", default: none)
      if source-record != none {
        let source-data = _data-patch(
          source,
          _record-with-fields(source-record, base-fields: edge-fields, extra: (edge: edge-record)),
          "graph.map source",
        )
        if source-data != none {
          edge-patch.insert("source", source-data)
        }
      }

      let sink-record = edge-source.at("sink", default: none)
      if sink-record != none {
        let sink-data = _data-patch(
          sink,
          _record-with-fields(sink-record, base-fields: edge-fields, extra: (edge: edge-record)),
          "graph.map sink",
        )
        if sink-data != none {
          edge-patch.insert("sink", sink-data)
        }
      }

      if edge-patch.len() > 1 {
        edge-patches.push(edge-patch)
      }
    }
  }
  if edge-patches.len() > 0 {
    patch.insert("edges", edge-patches)
  }

  if patch.len() == 0 {
    graph_
  } else {
    _plugin.graph_with_data(bytes(graph_), cbor.encode(patch))
  }
}

#let _apply-default-data(
  graph_,
  default-node-data: none,
  default-edge-data: none,
  default-source-data: none,
  default-sink-data: none,
) = {
  if (
    default-node-data == none
      and default-edge-data == none
      and default-source-data == none
      and default-sink-data == none
  ) {
    return graph_
  }
  map(
    graph_,
    node: record => _default-data-patch(default-node-data, record),
    edge: record => _default-data-patch(default-edge-data, record),
    source: record => _default-data-patch(default-source-data, record),
    sink: record => _default-data-patch(default-sink-data, record),
  )
}

/// Evaluate selected fields into opaque data entries.
///
/// Each selected field is read from the record's merged `fields` dictionary,
/// evaluated in a scope containing those fields, and written to
/// `data.<field>`.
///
/// ```example
/// #let g = graph.build(
///   {
///     graph.node(<a>)
///     graph.node(<b>)
///     graph.edge(graph.source(<a>), graph.sink(<b>), statements: (mom: "p"))
///   },
///   default-edge-statements: (display-label: "$#mom$"),
/// )
/// #let g = graph.eval-fields(g, eval-edge-fields: ("display-label",))
/// #graph.edges(g).first().data.at("display-label")
/// ```
/// -> bytes
#let eval-fields(
  graph_,

  /// Graph statement fields to evaluate into `graph.info(g).data`. -> string | array
  eval-graph-fields: (),

  /// Node statement fields to evaluate into `graph.nodes(g).at(i).data`. -> string | array
  eval-node-fields: (),

  /// Edge statement fields to evaluate into `graph.edges(g).at(i).data`. -> string | array
  eval-edge-fields: (),

  /// Source half-edge fields to evaluate into `edge.source.data`. -> string | array
  eval-source-fields: (),

  /// Sink half-edge fields to evaluate into `edge.sink.data`. -> string | array
  eval-sink-fields: (),

  /// Typst `eval` mode used for string field values. -> string
  eval-mode: "markup",

  /// Additional Typst names available while evaluating field values. -> dictionary
  scope: (:),
) = {
  let graph-fields = _eval-field-list(eval-graph-fields, "graph.eval-fields")
  let node-fields = _eval-field-list(eval-node-fields, "graph.eval-fields")
  let edge-fields = _eval-field-list(eval-edge-fields, "graph.eval-fields")
  let source-fields = _eval-field-list(eval-source-fields, "graph.eval-fields")
  let sink-fields = _eval-field-list(eval-sink-fields, "graph.eval-fields")
  if graph-fields.len() == 0 and node-fields.len() == 0 and edge-fields.len() == 0 and source-fields.len() == 0 and sink-fields.len() == 0 {
    return graph_
  }
  map(
    graph_,
    graph: record => {
      let data = _eval-data-value(record, graph-fields, eval-mode, scope)
      if data == none { none } else { (data: data) }
    },
    node: record => {
      let data = _eval-data-value(record, node-fields, eval-mode, scope)
      if data == none { none } else { (data: data) }
    },
    edge: record => {
      let data = _eval-data-value(record, edge-fields, eval-mode, scope)
      if data == none { none } else { (data: data) }
    },
    source: record => {
      let data = _eval-data-value(record, source-fields, eval-mode, scope)
      if data == none { none } else { (data: data) }
    },
    sink: record => {
      let data = _eval-data-value(record, sink-fields, eval-mode, scope)
      if data == none { none } else { (data: data) }
    },
  )
}

/// Parse one or more DOT digraphs into graph objects.
///
/// Default data are applied before `eval-*` fields are evaluated.
/// Parsed fields take precedence over default data fields.
///
/// ```example
/// #let graphs = graph.parse("digraph first { a -> b }")
/// #graphs.len()
/// ```
/// -> array
#let parse(
  input,

  /// Default data merged into every node data. Captured node data fields override it. -> any
  default-node-data: none,

  /// Default data merged into every edge data. Captured edge data fields override it. -> any
  default-edge-data: none,

  /// Default data merged into every source half-edge data. Captured source data fields override it. -> any
  default-source-data: none,

  /// Default data merged into every sink half-edge data. Captured sink data fields override it. -> any
  default-sink-data: none,

  /// Graph statement fields to evaluate into `graph.info(g).data`. -> string | array
  eval-graph-fields: (),

  /// Node statement fields to evaluate into `graph.nodes(g).at(i).data`. -> string | array
  eval-node-fields: (),

  /// Edge statement fields to evaluate into `graph.edges(g).at(i).data`. -> string | array
  eval-edge-fields: (),

  /// Source half-edge fields to evaluate into `edge.source.data`. -> string | array
  eval-source-fields: (),

  /// Sink half-edge fields to evaluate into `edge.sink.data`. -> string | array
  eval-sink-fields: (),

  /// Typst `eval` mode used for string field values. -> string
  eval-mode: "markup",

  /// Additional Typst names available while evaluating field values. -> dictionary
  scope: (:),
) = {
  let graph-fields = _eval-field-list(eval-graph-fields, "graph.parse")
  let node-fields = _eval-field-list(eval-node-fields, "graph.parse")
  let edge-fields = _eval-field-list(eval-edge-fields, "graph.parse")
  let source-fields = _eval-field-list(eval-source-fields, "graph.parse")
  let sink-fields = _eval-field-list(eval-sink-fields, "graph.parse")
  let graphs = cbor(_plugin.parse_graph(bytes(input)))
  let needs-eval = graph-fields.len() != 0 or node-fields.len() != 0 or edge-fields.len() != 0 or source-fields.len() != 0 or sink-fields.len() != 0
  let needs-defaults = (
    default-node-data != none
      or default-edge-data != none
      or default-source-data != none
      or default-sink-data != none
  )
  if not needs-eval and not needs-defaults {
    return graphs
  }
  graphs.map(graph => {
    let result = graph
    if needs-defaults {
      result = _apply-default-data(
        result,
        default-node-data: default-node-data,
        default-edge-data: default-edge-data,
        default-source-data: default-source-data,
        default-sink-data: default-sink-data,
      )
    }
    if needs-eval {
      result = eval-fields(
        result,
        eval-graph-fields: graph-fields,
        eval-node-fields: node-fields,
        eval-edge-fields: edge-fields,
        eval-source-fields: source-fields,
        eval-sink-fields: sink-fields,
        eval-mode: eval-mode,
        scope: scope,
      )
    }
    result
  })
}

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

  /// Opaque graph data. Any CBOR-encodable Typst value. -> any
  data: none,

  /// Flat graph statements. Used by DOT; values cannot nest. -> dictionary
  statements: (:),

  /// Flat default edge statements. Used by DOT; values cannot nest.
  /// -> dictionary
  default-edge-statements: (:),

  /// Default data merged into every node data. Captured node data fields override it. -> any
  default-node-data: none,

  /// Default data merged into every edge data. Captured edge data fields override it. -> any
  default-edge-data: none,

  /// Default data merged into every source half-edge data. Captured source data fields override it. -> any
  default-source-data: none,

  /// Default data merged into every sink half-edge data. Captured sink data fields override it. -> any
  default-sink-data: none,

  /// Flat default node statements. Used by DOT; values cannot nest. -> dictionary
  default-node-statements: (:),

  /// Additional node items or raw node specs. -> array
  /// -> array
  nodes: (),

  /// Additional edge items or raw edge specs. -> array
  edges: (),
) = {
  let split = _split-items(items.pos(), nodes, edges)
  let node-keys = _node-name-map(split.nodes)
  _check-edge-names(split.edges)
  _plugin.graph_from_spec(cbor.encode((
    name: name,
    data: _encode-data(data),
    statements: _flat-statements(statements, "graph.build statements"),
    default-edge-statements: _flat-statements(default-edge-statements, "graph.build default-edge-statements"),
    default-node-statements: _flat-statements(default-node-statements, "graph.build default-node-statements"),
    nodes: split.nodes.map(node => _resolved-node-spec(
      node,
      node-keys,
      default-data: default-node-data,
    )),
    edges: split.edges.map(edge => _resolved-edge-spec(
      edge,
      node-keys,
      default-edge-data: default-edge-data,
      default-source-data: default-source-data,
      default-sink-data: default-sink-data,
    )),
  )))
}

/// Create a graph node item for @build.
///
/// A Typst label is the node name used by @source, @sink, and @pos. The
/// optional numeric `id` fixes the resulting graph node index. Extra named
/// arguments are captured as node data fields. The default draw style uses
/// `data.label` as the visible node label when present.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>, id: 0, label: [A])
/// })
/// #graph.nodes(g).first().name
/// ```
/// -> array
#let node(
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
) = {
  let resolved-data = _data-from-args("graph.node", args)
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
  ((
    linnest-kind: "node",
    id: id,
    name: resolved-name,
    data: resolved-data,
    pos: pos,
    shift: shift,
    statements: statements,
  ),)
}

/// Create a source half-edge endpoint.
///
/// `node` may be a node name like `<a>` or a numeric node index. `name` gives
/// the half-edge a Typst name; `id` is a numeric half-edge order/index
/// override. Extra named arguments are captured as source data fields.
///
/// ```example
/// #graph.source(<a>, name: <h1>, id: 0, kind: "out", compass: "e")
/// ```
/// -> dictionary
#let source(
  node,
  ..args,
  name: none,
  id: none,
  statement: none,
  compass: none,
) = {
  if args.pos().len() > 0 {
    panic("graph.source: expected only one positional node reference")
  }
  let data = _data-from-args("graph.source", args)
  _check-name(name, "graph.source")
  _check-id(id, "graph.source")
  (linnest-kind: "source", node: node, name: name, id: id, data: data, statement: statement, compass: compass)
}

/// Create a sink half-edge endpoint.
///
/// `node` may be a node name like `<a>` or a numeric node index. `name` gives
/// the half-edge a Typst name; `id` is a numeric half-edge order/index
/// override. Extra named arguments are captured as sink data fields.
///
/// ```example
/// #graph.sink(<b>, name: <h2>, id: 1, kind: "in", compass: "w")
/// ```
/// -> dictionary
#let sink(
  node,
  ..args,
  name: none,
  id: none,
  statement: none,
  compass: none,
) = {
  if args.pos().len() > 0 {
    panic("graph.sink: expected only one positional node reference")
  }
  let data = _data-from-args("graph.sink", args)
  _check-name(name, "graph.sink")
  _check-id(id, "graph.sink")
  (linnest-kind: "sink", node: node, name: name, id: id, data: data, statement: statement, compass: compass)
}

/// Create a graph edge item for @build.
///
/// Positional arguments may contain one @source, one @sink, and optionally one
/// Typst label used as the edge name. The numeric `id` chooses the edge order
/// and extra named arguments are captured as edge data fields. The default
/// draw style uses `data.label` as the visible edge label when present.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>), <e>, graph.sink(<b>), label: [$p$], particle: "g")
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
) = {
  let resolved-data = _data-from-args("graph.edge", args)
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
  ((
    linnest-kind: "edge",
    source: resolved-source,
    sink: resolved-sink,
    orientation: orientation,
    flow: flow,
    name: resolved-name,
    id: resolved-id,
    data: resolved-data,
    pos: pos,
    shift: shift,
    label-pos: label-pos,
    label-angle: label-angle,
    bend: bend,
    statements: statements,
  ),)
}

/// Return graph metadata.
///
/// The result has `name`, `global-statements`, `default-edge-statements`, and
/// `default-node-statements`.
///
/// ```example
/// #let g = graph.build({ graph.node(<a>) }, name: "demo")
/// #graph.info(g).name
/// ```
/// -> dictionary
#let info(graph) = _decode-data-field(cbor(_plugin.graph_info(bytes(graph))))

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
  let records = if subgraph == none {
    cbor(_plugin.graph_nodes(bytes(graph)))
  } else {
    cbor(_plugin.graph_nodes_of_subgraph(bytes(graph), bytes(subgraph)))
  }
  records.map(_decode-data-field)
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
  let records = if subgraph == none {
    cbor(_plugin.graph_edges(bytes(graph)))
  } else {
    cbor(_plugin.graph_edges_of_subgraph(bytes(graph), bytes(subgraph)))
  }
  records.map(_decode-edge-data)
}

#let _name-key(value, context_) = {
  if type(value) == label {
    _label-key(value, context_)
  } else if type(value) == str {
    value
  } else {
    panic(context_ + ": expected a Typst label or string name")
  }
}

#let _find-record-by-name(records, key) = {
  let found = false
  for record in records {
    if record.at("name", default: none) == key {
      found = true
    }
  }
  found
}

#let _data-update(update, record, context_) = {
  let data = record.at("data", default: none)
  let result = if type(update) == function {
    update(data, record)
  } else {
    update
  }
  if result == none {
    none
  } else {
    (data: result)
  }
}

#let _named-data-callback(key, update, context_) = record => {
  if record.at("name", default: none) == key {
    _data-update(update, record, context_)
  } else {
    none
  }
}

/// Return one named node's opaque data.
///
/// `name` is a Typst label such as `<a>` or the corresponding string name.
///
/// ```example
/// #let g = graph.build({ graph.node(<a>, label: [A]) })
/// #graph.node-data(g, <a>).label
/// ```
/// -> any
#let node-data(graph_, name) = {
  let key = _name-key(name, "graph.node-data")
  decode-data(cbor(_plugin.graph_node_data_by_name(bytes(graph_), cbor.encode(key))))
}

/// Return one named edge's opaque data.
///
/// `name` is a Typst label such as `<e>` or the corresponding string name.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>), <e>, graph.sink(<b>), label: [$p$])
/// })
/// #graph.edge-data(g, <e>).label
/// ```
/// -> any
#let edge-data(graph_, name) = {
  let key = _name-key(name, "graph.edge-data")
  decode-data(cbor(_plugin.graph_edge_data_by_name(bytes(graph_), cbor.encode(key))))
}

/// Update one named node's opaque data.
///
/// `name` is a Typst label such as `<a>` or the corresponding string name.
/// `update` may be a replacement data value or a function
/// `(data, node) => new-data`.
///
/// ```example
/// #let g = graph.build({ graph.node(<a>) })
/// #let g = graph.update-node-data(g, <a>, (label: [A]))
/// #graph.nodes(g).first().data.label
/// ```
/// -> bytes
#let update-node-data(graph_, name, update) = {
  let key = _name-key(name, "graph.update-node-data")
  if update == none {
    return graph_
  }
  if type(update) != function {
    return _plugin.graph_set_node_data_by_name(
      bytes(graph_),
      cbor.encode((name: key, data: _encode-data(update))),
    )
  }
  if not _find-record-by-name(nodes(graph_), key) {
    panic("graph.update-node-data: no node named " + repr(name))
  }
  map(graph_, node: _named-data-callback(key, update, "graph.update-node-data"))
}

/// Update one named edge's opaque data.
///
/// `name` is a Typst label such as `<e>` or the corresponding string name.
/// `update` may be a replacement data value or a function
/// `(data, edge) => new-data`.
///
/// ```example
/// #let g = graph.build({
///   graph.node(<a>)
///   graph.node(<b>)
///   graph.edge(graph.source(<a>), <e>, graph.sink(<b>))
/// })
/// #let g = graph.update-edge-data(g, <e>, (label: [$p$]))
/// #graph.edges(g).first().data.label
/// ```
/// -> bytes
#let update-edge-data(graph_, name, update) = {
  let key = _name-key(name, "graph.update-edge-data")
  if update == none {
    return graph_
  }
  if type(update) != function {
    return _plugin.graph_set_edge_data_by_name(
      bytes(graph_),
      cbor.encode((name: key, data: _encode-data(update))),
    )
  }
  if not _find-record-by-name(edges(graph_), key) {
    panic("graph.update-edge-data: no edge named " + repr(name))
  }
  map(graph_, edge: _named-data-callback(key, update, "graph.update-edge-data"))
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
