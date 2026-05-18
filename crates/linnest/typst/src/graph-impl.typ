// Internal graph implementation. Public users should import `graph.typ`.

#let _plugin = plugin("../linnest.wasm")

#let _encode-payload-value(value) = {
  if type(value) == content {
    (linnest-payload-kind: "content", value: cbor.encode(value))
  } else if type(value) == array {
    value.map(_encode-payload-value)
  } else if type(value) == dictionary {
    let result = (:)
    for key in value.keys() {
      result.insert(key, _encode-payload-value(value.at(key)))
    }
    result
  } else {
    value
  }
}

#let _encode-payload(value) = if value == none { none } else { cbor.encode(_encode-payload-value(value)) }

#let _payload-bytes(value) = if type(value) == array { bytes(value) } else { value }

#let _sequence-content(children, decode-content) = {
  let result = []
  for child in children {
    result += decode-content(child)
  }
  result
}

#let _math-source(value) = {
  if type(value) != dictionary {
    none
  } else {
    let func = value.at("func", default: none)
    if func == "symbol" or func == "text" {
      value.at("text", default: "")
    } else if func == "space" {
      " "
    } else if func == "sequence" {
      value.children.map(_math-source).join("")
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

#let _decode-payload-value(value) = {
  if type(value) == array {
    value.map(_decode-payload-value)
  } else if type(value) == dictionary {
    let kind = value.at("linnest-payload-kind", default: none)
    if kind == "content" {
      _decode-content(cbor(_payload-bytes(value.value)), _decode-payload-value)
    } else {
      let result = (:)
      for key in value.keys() {
        result.insert(key, _decode-payload-value(value.at(key)))
      }
      result
    }
  } else {
    value
  }
}

#let decode-payload(value) = if value == none { none } else { _decode-payload-value(cbor(_payload-bytes(value))) }

#let _payload-with-label(payload, label, context_) = {
  if label == none {
    payload
  } else if payload == none {
    (label: label)
  } else if type(payload) == dictionary {
    if payload.keys().contains("label") {
      panic(context_ + ": label specified twice")
    }
    payload + (label: label)
  } else {
    panic(context_ + ": label requires payload to be a dictionary")
  }
}

#let _merge-payload(default, payload) = {
  if default == none {
    payload
  } else if payload == none {
    default
  } else if type(default) == dictionary and type(payload) == dictionary {
    default + payload
  } else {
    payload
  }
}

#let _decode-payload-field(record) = {
  let result = record
  let payload = result.at("payload", default: none)
  if payload != none {
    result.payload = decode-payload(payload)
  }
  result
}

#let _decode-edge-payloads(edge) = {
  let result = _decode-payload-field(edge)
  for side in ("source", "sink") {
    let endpoint = result.at(side, default: none)
    if endpoint != none {
      result.insert(side, _decode-payload-field(endpoint))
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
    panic(context_ + ": statement values must be flat scalars; use payload for nested data or Typst content")
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

#let _payload-from-args(context_, args) = {
  let named = args.named()
  if named.keys().contains("payload") {
    panic(context_ + ": use direct named payload fields instead of payload: (...)")
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

#let _node-spec(node, default-payload: none) = {
  let statements = _statements-with-point(_flat-statements(node.at("statements", default: (:)), "graph.node statements"), "shift", node.at("shift", default: none))
  let payload = _encode-payload(_merge-payload(
    default-payload,
    _payload-with-label(node.at("payload", default: none), node.at("label", default: none), "graph.node"),
  ))
  let clean = node
  let name = clean.at("name", default: none)
  for key in ("linnest-kind", "key", "id", "shift", "name", "label", "payload") {
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
  if payload != none {
    clean.insert("payload", payload)
  }
  clean + (statements: statements)
}

#let _resolved-node-spec(node, node-ids, node-keys, default-payload: none) = {
  let name = if type(node.at("name", default: none)) == label { node.name } else { none }
  let index = _resolve-named-id(name, node.at("id", default: none), node-ids, "graph.node")
  let item = node + (pos: _resolve-pos(node.at("pos", default: none), node-keys))
  if index != none {
    item.insert("index", index)
  }
  _node-spec(item, default-payload: default-payload)
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

#let _edge-spec(edge, default-payload: none) = {
  let statements = _edge-statements(edge)
  let name = edge.at("name", default: none)
  if name != none {
    statements.insert("__linnest-edge-name", _label-key(name, "graph.edge"))
  }
  let payload = _encode-payload(_merge-payload(
    default-payload,
    _payload-with-label(edge.at("payload", default: none), edge.at("label", default: none), "graph.edge"),
  ))
  let clean = edge
  for key in (
    "linnest-kind",
    "shift",
    "label-pos",
    "label-angle",
    "bend",
    "label",
    "payload",
    "name",
  ) {
    if clean.keys().contains(key) {
      let _ = clean.remove(key)
    }
  }
  if payload != none {
    clean.insert("payload", payload)
  }
  clean + (statements: statements)
}

#let _resolved-half-spec(half, node-keys, half-ids, context_, default-payload: none) = {
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
  let payload = _encode-payload(_merge-payload(default-payload, half.at("payload", default: none)))
  if payload != none {
    result.insert("payload", payload)
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
  edge-ids,
  half-ids,
  default-edge-payload: none,
  default-source-payload: none,
  default-sink-payload: none,
) = {
  let id = _resolve-named-id(edge.at("name", default: none), edge.at("id", default: none), edge-ids, "graph.edge")
  let item = edge + (
    source: _resolved-half-spec(
      edge.at("source", default: none),
      node-keys,
      half-ids,
      "graph.source",
      default-payload: default-source-payload,
    ),
    sink: _resolved-half-spec(
      edge.at("sink", default: none),
      node-keys,
      half-ids,
      "graph.sink",
      default-payload: default-sink-payload,
    ),
    pos: _resolve-pos(edge.at("pos", default: none), node-keys),
  )
  if id != none {
    item.id = id
  }
  _edge-spec(item, default-payload: default-edge-payload)
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
  let payload = record.at("payload", default: none)
  if type(payload) == dictionary {
    fields += payload
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

#let _payload-field-value(record, field) = {
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

#let _eval-payload-value(record, fields, eval-mode, scope) = {
  let existing = record.at("payload", default: none)
  let payload = if type(existing) == dictionary { existing } else { (:) }
  let changed = false
  let local-scope = _record-eval-scope(record, scope)
  for field in fields {
    let value = _payload-field-value(record, field)
    if value != none {
      payload.insert(field, _eval-dot-value(value, eval-mode, local-scope))
      changed = true
    }
  }
  if not changed {
    none
  } else {
    payload
  }
}

#let _mapped-payload(callback, record, context_) = {
  if callback == none {
    return none
  }
  let result = callback(record)
  if result == none {
    none
  } else if type(result) != dictionary {
    panic(context_ + ": callback must return none or a dictionary")
  } else {
    result.at("payload", default: none)
  }
}

#let _payload-patch(callback, record, context_) = {
  let payload = _mapped-payload(callback, record, context_)
  if payload == none {
    none
  } else {
    _encode-payload(payload)
  }
}

#let _default-payload-patch(default, record) = {
  let payload = _merge-payload(default, record.at("payload", default: none))
  if payload == none {
    none
  } else {
    (payload: payload)
  }
}

/// Map graph metadata to new opaque payloads.
///
/// The callbacks receive decoded records plus a `fields` dictionary containing
/// merged statements and direct record fields. A callback returns `none` to
/// leave the record unchanged, or `(payload: value)` to set a new opaque payload.
///
/// ```example
/// #let g = graph.build({ graph.node(<a>) })
/// #let g = graph.map(g, node: node => (payload: (label: [A])))
/// #graph.nodes(g).first().payload.label
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
  let info = _decode-payload-field(cbor(_plugin.graph_info(bytes(graph_))))
  let graph-record = _record-with-fields(info + (statements: info.at("global-statements", default: (:))))
  let graph-payload = _payload-patch(graph, graph-record, "graph.map graph")
  if graph-payload != none {
    patch.insert("payload", graph-payload)
  }

  let node-patches = ()
  if node != none {
    for node-record in cbor(_plugin.graph_nodes(bytes(graph_))).map(_decode-payload-field) {
      let payload = _payload-patch(node, _record-with-fields(node-record), "graph.map node")
      if payload != none {
        node-patches.push((index: node-record.node, payload: payload))
      }
    }
  }
  if node-patches.len() > 0 {
    patch.insert("nodes", node-patches)
  }

  let edge-patches = ()
  if edge != none or source != none or sink != none {
    for edge-source in cbor(_plugin.graph_edges(bytes(graph_))).map(_decode-edge-payloads) {
      let edge-patch = (index: edge-source.edge)
      let edge-record = _record-with-fields(edge-source)
      let edge-fields = edge-record.fields
      let payload = _payload-patch(edge, edge-record, "graph.map edge")
      if payload != none {
        edge-patch.insert("payload", payload)
      }

      let source-record = edge-source.at("source", default: none)
      if source-record != none {
        let source-payload = _payload-patch(
          source,
          _record-with-fields(source-record, base-fields: edge-fields, extra: (edge: edge-record)),
          "graph.map source",
        )
        if source-payload != none {
          edge-patch.insert("source", source-payload)
        }
      }

      let sink-record = edge-source.at("sink", default: none)
      if sink-record != none {
        let sink-payload = _payload-patch(
          sink,
          _record-with-fields(sink-record, base-fields: edge-fields, extra: (edge: edge-record)),
          "graph.map sink",
        )
        if sink-payload != none {
          edge-patch.insert("sink", sink-payload)
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
    _plugin.graph_with_payloads(bytes(graph_), cbor.encode(patch))
  }
}

#let _apply-default-payloads(
  graph_,
  default-node-payload: none,
  default-edge-payload: none,
  default-source-payload: none,
  default-sink-payload: none,
) = {
  if (
    default-node-payload == none
      and default-edge-payload == none
      and default-source-payload == none
      and default-sink-payload == none
  ) {
    return graph_
  }
  map(
    graph_,
    node: record => _default-payload-patch(default-node-payload, record),
    edge: record => _default-payload-patch(default-edge-payload, record),
    source: record => _default-payload-patch(default-source-payload, record),
    sink: record => _default-payload-patch(default-sink-payload, record),
  )
}

/// Evaluate selected fields into opaque payload entries.
///
/// Each selected field is read from the record's merged `fields` dictionary,
/// evaluated in a scope containing those fields, and written to
/// `payload.<field>`.
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
/// #graph.edges(g).first().payload.at("display-label")
/// ```
/// -> bytes
#let eval-fields(
  graph_,

  /// Graph statement fields to evaluate into `graph.info(g).payload`. -> string | array
  eval-graph-fields: (),

  /// Node statement fields to evaluate into `graph.nodes(g).at(i).payload`. -> string | array
  eval-node-fields: (),

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
      let payload = _eval-payload-value(record, graph-fields, eval-mode, scope)
      if payload == none { none } else { (payload: payload) }
    },
    node: record => {
      let payload = _eval-payload-value(record, node-fields, eval-mode, scope)
      if payload == none { none } else { (payload: payload) }
    },
    edge: record => {
      let payload = _eval-payload-value(record, edge-fields, eval-mode, scope)
      if payload == none { none } else { (payload: payload) }
    },
    source: record => {
      let payload = _eval-payload-value(record, source-fields, eval-mode, scope)
      if payload == none { none } else { (payload: payload) }
    },
    sink: record => {
      let payload = _eval-payload-value(record, sink-fields, eval-mode, scope)
      if payload == none { none } else { (payload: payload) }
    },
  )
}

/// Parse one or more DOT digraphs into graph objects.
///
/// Default payloads are applied before `eval-*` fields are evaluated.
/// Parsed fields take precedence over default payload fields.
///
/// ```example
/// #let graphs = graph.parse("digraph first { a -> b }")
/// #graphs.len()
/// ```
/// -> array
#let parse(
  input,

  /// Default payload merged into every node payload. Captured node payload fields override it. -> any
  default-node-payload: none,

  /// Default payload merged into every edge payload. Captured edge payload fields override it. -> any
  default-edge-payload: none,

  /// Default payload merged into every source half-edge payload. Captured source payload fields override it. -> any
  default-source-payload: none,

  /// Default payload merged into every sink half-edge payload. Captured sink payload fields override it. -> any
  default-sink-payload: none,

  /// Graph statement fields to evaluate into `graph.info(g).payload`. -> string | array
  eval-graph-fields: (),

  /// Node statement fields to evaluate into `graph.nodes(g).at(i).payload`. -> string | array
  eval-node-fields: (),

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
) = {
  let graph-fields = _eval-field-list(eval-graph-fields, "graph.parse")
  let node-fields = _eval-field-list(eval-node-fields, "graph.parse")
  let edge-fields = _eval-field-list(eval-edge-fields, "graph.parse")
  let source-fields = _eval-field-list(eval-source-fields, "graph.parse")
  let sink-fields = _eval-field-list(eval-sink-fields, "graph.parse")
  let graphs = cbor(_plugin.parse_graph(bytes(input)))
  let needs-eval = graph-fields.len() != 0 or node-fields.len() != 0 or edge-fields.len() != 0 or source-fields.len() != 0 or sink-fields.len() != 0
  let needs-defaults = (
    default-node-payload != none
      or default-edge-payload != none
      or default-source-payload != none
      or default-sink-payload != none
  )
  if not needs-eval and not needs-defaults {
    return graphs
  }
  graphs.map(graph => {
    let result = graph
    if needs-defaults {
      result = _apply-default-payloads(
        result,
        default-node-payload: default-node-payload,
        default-edge-payload: default-edge-payload,
        default-source-payload: default-source-payload,
        default-sink-payload: default-sink-payload,
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

  /// Opaque graph payload. Any CBOR-encodable Typst value. -> any
  payload: none,

  /// Flat graph statements. Used by DOT; values cannot nest. -> dictionary
  statements: (:),

  /// Flat default edge statements. Used by DOT; values cannot nest.
  /// -> dictionary
  default-edge-statements: (:),

  /// Default payload merged into every node payload. Captured node payload fields override it. -> any
  default-node-payload: none,

  /// Default payload merged into every edge payload. Captured edge payload fields override it. -> any
  default-edge-payload: none,

  /// Default payload merged into every source half-edge payload. Captured source payload fields override it. -> any
  default-source-payload: none,

  /// Default payload merged into every sink half-edge payload. Captured sink payload fields override it. -> any
  default-sink-payload: none,

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
  let node-ids = _collect-node-ids(split.nodes)
  let edge-ids = _collect-edge-ids(split.edges)
  let half-ids = _collect-half-ids(split.edges)
  _plugin.graph_from_spec(cbor.encode((
    name: name,
    payload: _encode-payload(payload),
    statements: _flat-statements(statements, "graph.build statements"),
    default-edge-statements: _flat-statements(default-edge-statements, "graph.build default-edge-statements"),
    default-node-statements: _flat-statements(default-node-statements, "graph.build default-node-statements"),
    nodes: split.nodes.map(node => _resolved-node-spec(
      node,
      node-ids,
      node-keys,
      default-payload: default-node-payload,
    )),
    edges: split.edges.map(edge => _resolved-edge-spec(
      edge,
      node-keys,
      edge-ids,
      half-ids,
      default-edge-payload: default-edge-payload,
      default-source-payload: default-source-payload,
      default-sink-payload: default-sink-payload,
    )),
  )))
}

/// Create a graph node item for @build.
///
/// A Typst label is the node name used by @source, @sink, and @pos. The
/// optional numeric `id` chooses the node order/index. `label` is stored as
/// `payload.label`. Extra named arguments are captured as node payload fields.
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
  /// Numeric node order/index override. -> none | int
  id: none,
  /// Visible node label, stored as `payload.label`. -> any
  label: none,
  /// Node placement. -> none | dictionary
  pos: none,
  /// Drawing shift stored as a statement. -> none | string | array | dictionary
  shift: none,
  /// Additional flat node statements. Used by DOT; values cannot nest. -> dictionary
  statements: (:),
) = {
  let resolved-payload = _payload-from-args("graph.node", args)
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
    payload: resolved-payload,
    label: label,
    pos: pos,
    shift: shift,
    statements: statements,
  ),)
}

/// Create a source half-edge endpoint.
///
/// `node` may be a node name like `<a>` or a numeric node index. `name` gives
/// the half-edge a Typst name; `id` is a numeric half-edge order/index
/// override. Extra named arguments are captured as source payload fields.
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
  let payload = _payload-from-args("graph.source", args)
  _check-name(name, "graph.source")
  _check-id(id, "graph.source")
  (linnest-kind: "source", node: node, name: name, id: id, payload: payload, statement: statement, compass: compass)
}

/// Create a sink half-edge endpoint.
///
/// `node` may be a node name like `<a>` or a numeric node index. `name` gives
/// the half-edge a Typst name; `id` is a numeric half-edge order/index
/// override. Extra named arguments are captured as sink payload fields.
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
  let payload = _payload-from-args("graph.sink", args)
  _check-name(name, "graph.sink")
  _check-id(id, "graph.sink")
  (linnest-kind: "sink", node: node, name: name, id: id, payload: payload, statement: statement, compass: compass)
}

/// Create a graph edge item for @build.
///
/// Positional arguments may contain one @source, one @sink, and optionally one
/// Typst label used as the edge name. The numeric `id` chooses the edge order
/// and the visible label is `label`, stored as `payload.label`. Extra named
/// arguments are captured as edge payload fields.
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
  /// Visible edge label, stored as `payload.label`. -> any
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
  /// Additional flat edge statements. Used by DOT; values cannot nest. -> dictionary
  statements: (:),
) = {
  let resolved-payload = _payload-from-args("graph.edge", args)
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
    payload: resolved-payload,
    label: label,
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
#let info(graph) = _decode-payload-field(cbor(_plugin.graph_info(bytes(graph))))

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
  records.map(_decode-payload-field)
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
  records.map(_decode-edge-payloads)
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
  let data = record.at("payload", default: none)
  let result = if type(update) == function {
    update(data, record)
  } else {
    update
  }
  if result == none {
    none
  } else {
    (payload: result)
  }
}

#let _named-data-callback(key, update, context_) = record => {
  if record.at("name", default: none) == key {
    _data-update(update, record, context_)
  } else {
    none
  }
}

/// Return one named node's opaque payload.
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
  decode-payload(cbor(_plugin.graph_node_payload_by_name(bytes(graph_), cbor.encode(key))))
}

/// Return one named edge's opaque payload.
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
  decode-payload(cbor(_plugin.graph_edge_payload_by_name(bytes(graph_), cbor.encode(key))))
}

/// Update one named node's opaque payload.
///
/// `name` is a Typst label such as `<a>` or the corresponding string name.
/// `update` may be a replacement data value or a function
/// `(data, node) => new-data`.
///
/// ```example
/// #let g = graph.build({ graph.node(<a>) })
/// #let g = graph.update-node-data(g, <a>, (label: [A]))
/// #graph.nodes(g).first().payload.label
/// ```
/// -> bytes
#let update-node-data(graph_, name, update) = {
  let key = _name-key(name, "graph.update-node-data")
  if update == none {
    return graph_
  }
  if type(update) != function {
    return _plugin.graph_set_node_payload_by_name(
      bytes(graph_),
      cbor.encode((name: key, payload: _encode-payload(update))),
    )
  }
  if not _find-record-by-name(nodes(graph_), key) {
    panic("graph.update-node-data: no node named " + repr(name))
  }
  map(graph_, node: _named-data-callback(key, update, "graph.update-node-data"))
}

/// Update one named edge's opaque payload.
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
/// #graph.edges(g).first().payload.label
/// ```
/// -> bytes
#let update-edge-data(graph_, name, update) = {
  let key = _name-key(name, "graph.update-edge-data")
  if update == none {
    return graph_
  }
  if type(update) != function {
    return _plugin.graph_set_edge_payload_by_name(
      bytes(graph_),
      cbor.encode((name: key, payload: _encode-payload(update))),
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
