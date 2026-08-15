// Internal graph implementation. Public users should import `graph.typ`.

#let _plugin = plugin("../../linnest.wasm")

#let _graph-kind = "linnest-graph"
#let _style-key = "linnest-style"

#let _empty-native-data() = (
  graph: none,
  nodes: (),
  edges: (),
  hedges: (),
)

#let _is-graph(value) = type(value) == dictionary and value.at("linnest-kind", default: none) == _graph-kind

#let _assert-graph(graph) = {
  if not _is-graph(graph) {
    panic("graph: expected a Linnest graph object returned by graph.build or graph.parse")
  }
  graph
}

#let graph-bytes(graph) = _assert-graph(graph).bytes

#let _native-data(graph) = _assert-graph(graph).native-data

#let _as-bytes(value) = if type(value) == array { bytes(value) } else { value }

#let _graph-object(graph-bytes, native-data) = (
  linnest-kind: _graph-kind,
  bytes: _as-bytes(graph-bytes),
  native-data: native-data,
)

#let with-bytes(graph, graph-bytes) = _graph-object(graph-bytes, _native-data(graph))

#let _array-at(values, index) = if index == none {
  none
} else if type(values) == array {
  values.at(index, default: none)
} else {
  none
}

#let _array-set(values, index, value) = {
  let result = ()
  let len = calc.max(values.len(), index + 1)
  for i in range(len) {
    result.push(if i == index { value } else { values.at(i, default: none) })
  }
  result
}

#let _with-native-data(graph, native-data) = _graph-object(_assert-graph(graph).bytes, native-data)

#let _payload(value) = if value == none { none } else { cbor.encode(value) }

#let _decode-payload(value) = if value == none { none } else { cbor(_as-bytes(value)) }

#let _payload-field(payload, key) = if type(payload) == dictionary {
  payload.at(key, default: none)
} else {
  none
}

#let _get-native-data(native-data, kind, index) = {
  if kind == "graph" {
    native-data.graph
  } else if kind == "node" {
    _array-at(native-data.nodes, index)
  } else if kind == "edge" {
    _array-at(native-data.edges, index)
  } else if kind == "hedge" {
    _array-at(native-data.hedges, index)
  } else {
    none
  }
}

#let _clean-statements(statements) = statements

#let _with-build-payload(item, data-key) = {
  let result = item
  let payload = (data-key: data-key)
  let name = result.at("name", default: none)
  if result.at("linnest-kind", default: none) == "edge" and name != none {
    if type(name) != label {
      panic("graph.edge: name must be a Typst label")
    }
    payload.name = str(name)
  }
  result.payload = _payload(payload)
  result
}

#let _with-build-payloads(items) = {
  let result = ()
  for (index, item) in items.enumerate() {
    result.push(_with-build-payload(item, index))
  }
  result
}

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

#let _name-label(value, context_) = {
  if value == none {
    none
  } else if type(value) == label {
    value
  } else if type(value) == str {
    label(value)
  } else {
    panic(context_ + ": name must be a Typst label or string")
  }
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

#let _node-spec(node) = {
  let statements = _statements-with-point(_flat-statements(node.at("statements", default: (:)), "graph.node statements"), "shift", node.at("shift", default: none))
  let clean = node
  let name = clean.at("name", default: none)
  let payload = clean.at("payload", default: none)
  for key in ("linnest-kind", "key", "id", "shift", "name", "data", "payload") {
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
    clean.insert("data", payload)
  }
  clean + (statements: statements)
}

#let _resolved-node-spec(node, node-keys) = {
  let index = node.at("id", default: none)
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

#let _edge-spec(edge) = {
  let statements = _edge-statements(edge)
  let payload = edge.at("payload", default: none)
  let clean = edge
  for key in (
    "linnest-kind",
    "shift",
    "label-pos",
    "label-angle",
    "bend",
    "data",
    "name",
    "payload",
  ) {
    if clean.keys().contains(key) {
      let _ = clean.remove(key)
    }
  }
  if payload != none {
    clean.insert("data", payload)
  }
  clean + (statements: statements)
}

#let _resolved-half-spec(half, node-keys, context_) = {
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
) = {
  let item = edge + (
    source: _resolved-half-spec(
      edge.at("source", default: none),
      node-keys,
      "graph.source",
    ),
    sink: _resolved-half-spec(
      edge.at("sink", default: none),
      node-keys,
      "graph.sink",
    ),
    pos: _resolve-pos(edge.at("pos", default: none), node-keys),
  )
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
#let group(name, side) = {
  if side != none and not (side == "+" or side == "-" or side == "positive" or side == "negative") {
    panic("graph.group: side must be none, \"+\", \"-\", \"positive\", or \"negative\"")
  }
  (kind: "group", name: _statement-value(name, "graph.group name"), side: side)
}
#let _axis-placement-kind = "axis-placement"
#let _axis-placement(mode, value, context_) = {
  if value == none {
    panic(context_ + ": expected an axis coordinate value")
  }
  (kind: _axis-placement-kind, mode: mode, value: value)
}
#let pin(value) = _axis-placement("pin", value, "graph.pin")
#let start(value) = _axis-placement("start", value, "graph.start")
#let _axis-value(value, context_) = {
  if type(value) == dictionary and value.at("kind", default: none) == _axis-placement-kind {
    let mode = value.at("mode", default: none)
    if mode != "start" and mode != "pin" {
      panic(context_ + ": axis mode must be \"start\" or \"pin\"")
    }
    (value: value.at("value", default: none), mode: mode)
  } else {
    (value: value, mode: none)
  }
}
#let pos(options) = {
  let x = _axis-value(options.x, "graph.pos x")
  let y = _axis-value(options.y, "graph.pos y")
  let ref = options.ref
  let dx = options.dx
  let dy = options.dy
  let mode = options.mode
  let x-mode = options.at("x-mode", default: x.mode)
  let y-mode = options.at("y-mode", default: y.mode)
  if mode != "start" and mode != "pin" {
    panic("graph.pos: mode must be \"start\" or \"pin\"")
  }
  if x-mode != none and x-mode != "start" and x-mode != "pin" {
    panic("graph.pos: x-mode must be none, \"start\", or \"pin\"")
  }
  if y-mode != none and y-mode != "start" and y-mode != "pin" {
    panic("graph.pos: y-mode must be none, \"start\", or \"pin\"")
  }
  let result = (mode: mode)
  if x-mode != none {
    result.insert("x-mode", x-mode)
  }
  if y-mode != none {
    result.insert("y-mode", y-mode)
  }
  if x.value != none {
    result.insert("x", x.value)
  }
  if y.value != none {
    result.insert("y", y.value)
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

#let _record-with-fields(record, base-fields, extra) = {
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

#let _canvas-length(unit) = {
  if type(unit) in (int, float) {
    unit * 1em
  } else {
    unit
  }
}

#let _call(value, data) = {
  if type(value) == function {
    value(data)
  } else {
    value
  }
}

#let _style(value, data) = {
  let value = _call(value, data)
  if value == none { (:) } else { value }
}

#let _as-content(value) = if value == none {
  none
} else if type(value) == str or type(value) == label {
  [#str(value)]
} else {
  value
}

#let _content(value, data, default) = {
  if value == auto {
    default
  } else {
    _as-content(_call(value, data))
  }
}

#let _data-label(record) = {
  let data = record.at("data", default: none)
  if type(data) == dictionary {
    data.at("label", default: none)
  } else {
    none
  }
}

#let _data-fields(record) = {
  let data = record.at("data", default: none)
  if type(data) == dictionary { data } else { (:) }
}

#let _padding-side(padding, side) = {
  if type(padding) in (int, float) {
    0
  } else if type(padding) == dictionary {
    padding.at(side, default: 0)
  } else {
    0
  }
}

#let _numeric(value, default: 0) = {
  if type(value) in (int, float) {
    value
  } else {
    default
  }
}

#let _node-style-data(record, style-record) = {
  let record-label = _data-label(record)
  let statement-label = record.statements.at("label", default: none)
  let data-label = if record-label == none { statement-label } else { record-label }
  (
    style-record.at("scope", default: (:))
      + record.statements
      + _data-fields(record)
      + (
        vid: record.node,
        node: record,
        name: record.name,
        data: record.at("data", default: none),
        label: data-label,
      )
  )
}

#let _node-layout-statements(record, style-record) = {
  let node-data = _node-style-data(record, style-record)
  let default-label-value = node-data.at("label", default: record.name)
  let default-label = _as-content(default-label-value)
  let label = _content(style-record.at("node-label", default: auto), node-data, default-label)
  let content-style = _style(style-record.at("node-label-style", default: (:)), node-data)
  let draw-style = _style(style-record.at("node-style", default: (:)), node-data)
  let unit = style-record.at("unit", default: 1)
  let length = _canvas-length(unit)

  let width = 0.0
  let height = 0.0
  if label != none {
    let size = measure(text(top-edge: "cap-height", bottom-edge: "bounds", label))
    width = calc.abs(size.width.to-absolute() / length.to-absolute())
    height = calc.abs(size.height.to-absolute() / length.to-absolute())
  }

  let padding = content-style.at("padding", default: 0)
  width += _numeric(_padding-side(padding, "left")) + _numeric(_padding-side(padding, "right"))
  height += _numeric(_padding-side(padding, "top")) + _numeric(_padding-side(padding, "bottom"))
  if type(padding) in (int, float) {
    width += 2 * padding
    height += 2 * padding
  }

  let radius = draw-style.at("radius", default: none)
  if type(radius) in (int, float) {
    width = calc.max(width, 2 * radius)
    height = calc.max(height, 2 * radius)
  } else if type(radius) == array and radius.len() > 0 {
    let radius = calc.max(..radius)
    width = calc.max(width, 2 * radius)
    height = calc.max(height, 2 * radius)
  }

  ("layout-width": width, "layout-height": height)
}

#let _edge-style-data(record, style-record) = {
  let source-half-edge = record.at("source", default: none)
  let sink-half-edge = record.at("sink", default: none)
  let source-statement = if source-half-edge == none { none } else { source-half-edge.statement }
  let sink-statement = if sink-half-edge == none { none } else { sink-half-edge.statement }
  let ext = source-half-edge == none or sink-half-edge == none
  let record-label = _data-label(record)
  let statement-label = record.statements.at("label", default: none)
  let data-label = if record-label == none { statement-label } else { record-label }
  (
    style-record.at("scope", default: (:))
      + record.statements
      + _data-fields(record)
      + (
        eid: record.edge,
        edge: record,
        source-statement: source-statement,
        sink-statement: sink-statement,
        source-half-edge: source-half-edge,
        sink-half-edge: sink-half-edge,
        data: record.at("data", default: none),
        label: data-label,
        orientation: record.orientation,
        ext: ext,
      )
  )
}

#let _edge-label-layout-statements(record, style-record) = {
  let edge-data = _edge-style-data(record, style-record)
  let label = _content(
    style-record.at("edge-label", default: none),
    edge-data,
    _as-content(edge-data.at("label", default: none)),
  )
  if label == none {
    return none
  }

  let content-style = _style(style-record.at("edge-label-style", default: (:)), edge-data)
  let length = _canvas-length(style-record.at("unit", default: 1))
  let size = measure(text(top-edge: "cap-height", bottom-edge: "bounds", label))
  let padding = content-style.at("padding", default: 0)
  let width = calc.abs(size.width.to-absolute() / length.to-absolute())
  let height = calc.abs(size.height.to-absolute() / length.to-absolute())
  width += _numeric(_padding-side(padding, "left")) + _numeric(_padding-side(padding, "right"))
  height += _numeric(_padding-side(padding, "top")) + _numeric(_padding-side(padding, "bottom"))
  if type(padding) in (int, float) {
    width += 2 * padding
    height += 2 * padding
  }
  ("label-width": width, "label-height": height)
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

#let _structural-keys(kind) = {
  if kind == "node" {
    ("pos", "shift", "statements")
  } else if kind == "edge" {
    ("pos", "shift", "label-pos", "label-angle", "bend", "statements")
  } else {
    ()
  }
}

#let _mapped-patch(callback, record, kind, context_) = {
  if callback == none {
    return (data: none, structural: none)
  }
  let result = callback(record)
  if result == none {
    return (data: none, structural: none)
  }
  if type(result) != dictionary {
    panic(context_ + ": callback must return none or a dictionary")
  }

  let structural-keys = _structural-keys(kind)
  let structural = (:)
  for key in structural-keys {
    if result.keys().contains(key) {
      if key == "statements" {
        structural.insert(key, _flat-statements(result.at(key), context_ + " statements"))
      } else {
        structural.insert(key, result.at(key))
      }
    }
  }

  let data = none
  if result.keys().contains("data") {
    data = result.data
  } else {
    let data-patch = (:)
    for key in result.keys() {
      if not structural-keys.contains(key) {
        data-patch.insert(key, result.at(key))
      }
    }
    if data-patch.len() != 0 {
      let existing = record.at("data", default: none)
      data = if type(existing) == dictionary {
        existing + data-patch
      } else {
        data-patch
      }
    }
  }

  (
    data: data,
    structural: if structural.len() == 0 { none } else { structural },
  )
}

#let _default-data-patch(default, record) = {
  let data = _merge-data(default, record.at("data", default: none))
  if data == none {
    none
  } else {
    (data: data)
  }
}

#let _native-data-from-build(graph-bytes_, graph-data, nodes, edges, default-node-data, default-edge-data, default-source-data, default-sink-data) = {
  let native-data = _empty-native-data()
  native-data.graph = graph-data

  for node-record in cbor(_plugin.graph_nodes(graph-bytes_)) {
    let payload = _decode-payload(node-record.at("data", default: none))
    let key = _payload-field(payload, "data-key")
    if key != none {
      let node = nodes.at(key)
      let data = _merge-data(default-node-data, node.at("data", default: none))
      if data != none {
        native-data.nodes = _array-set(native-data.nodes, node-record.node, data)
      }
    }
  }

  for edge-record in cbor(_plugin.graph_edges(graph-bytes_)) {
    let payload = _decode-payload(edge-record.at("data", default: none))
    let key = _payload-field(payload, "data-key")
    if key != none {
      let edge = edges.at(key)
      let data = _merge-data(default-edge-data, edge.at("data", default: none))
      if data != none {
        native-data.edges = _array-set(native-data.edges, edge-record.edge, data)
      }

      let source-record = edge-record.at("source", default: none)
      let source-item = edge.at("source", default: none)
      if source-record != none and source-item != none {
        let source-data = _merge-data(default-source-data, source-item.at("data", default: none))
        if source-data != none {
          native-data.hedges = _array-set(native-data.hedges, source-record.hedge, source-data)
        }
      }

      let sink-record = edge-record.at("sink", default: none)
      let sink-item = edge.at("sink", default: none)
      if sink-record != none and sink-item != none {
        let sink-data = _merge-data(default-sink-data, sink-item.at("data", default: none))
        if sink-data != none {
          native-data.hedges = _array-set(native-data.hedges, sink-record.hedge, sink-data)
        }
      }
    }
  }

  native-data
}

#let _info-record(graph) = {
  let result = cbor(_plugin.graph_info(graph-bytes(graph)))
  result.data = _native-data(graph).graph
  result
}

#let _native-record-data(record, native-data, kind, index) = {
  let result = record
  if result.at("statements", default: none) != none {
    result.statements = _clean-statements(result.statements)
  }
  if kind == "node" and result.at("name", default: none) != none {
    result.name = _name-label(result.name, "graph.nodes")
  }
  result.data = _get-native-data(native-data, kind, index)
  result
}

#let _native-edge-record(record, native-data) = {
  let payload = _decode-payload(record.at("data", default: none))
  let result = _native-record-data(record, native-data, "edge", record.edge)
  let name = _payload-field(payload, "name")
  if name != none {
    result.name = _name-label(name, "graph.edges")
  } else if result.at("name", default: none) != none {
    result.name = _name-label(result.name, "graph.edges")
  }
  for side in ("source", "sink") {
    let endpoint = result.at(side, default: none)
    if endpoint != none {
      result.insert(side, _native-record-data(endpoint, native-data, "hedge", endpoint.hedge))
    }
  }
  result
}

#let _node-records(graph, subgraph) = {
  let native-data = _native-data(graph)
  let records = if subgraph == none {
    cbor(_plugin.graph_nodes(graph-bytes(graph)))
  } else {
    cbor(_plugin.graph_nodes_of_subgraph(graph-bytes(graph), bytes(subgraph)))
  }
  records.map(record => _native-record-data(record, native-data, "node", record.node))
}

#let _edge-records(graph, subgraph) = {
  let native-data = _native-data(graph)
  let records = if subgraph == none {
    cbor(_plugin.graph_edges(graph-bytes(graph)))
  } else {
    cbor(_plugin.graph_edges_of_subgraph(graph-bytes(graph), bytes(subgraph)))
  }
  records.map(record => _native-edge-record(record, native-data))
}
#let map(graph_, callbacks) = {
  let graph = callbacks.graph
  let node = callbacks.node
  let edge = callbacks.edge
  let source = callbacks.source
  let sink = callbacks.sink
  let changed = false
  let structural-changed = false
  let structural-patches = (nodes: (), edges: ())
  let native-data = _native-data(graph_)
  let info = _info-record(graph_)
  let graph-record = _record-with-fields(info + (statements: info.at("global-statements", default: (:))), (:), (:))
  let graph-patch = _mapped-patch(graph, graph-record, "graph", "graph.map graph")
  if graph-patch.data != none {
    native-data.graph = graph-patch.data
    changed = true
  }

  if node != none {
    for node-record in _node-records(graph_, none) {
      let node-record = _record-with-fields(node-record, (:), (:))
      let patch = _mapped-patch(node, node-record, "node", "graph.map node")
      if patch.data != none {
        native-data.nodes = _array-set(native-data.nodes, node-record.node, patch.data)
        changed = true
      }
      if patch.structural != none {
        let structural = patch.structural + (index: node-record.node)
        structural-patches.nodes.push(structural)
        structural-changed = true
      }
    }
  }

  if edge != none or source != none or sink != none {
    for edge-source in _edge-records(graph_, none) {
      let edge-record = _record-with-fields(edge-source, (:), (:))
      let edge-fields = edge-record.fields
      let patch = _mapped-patch(edge, edge-record, "edge", "graph.map edge")
      if patch.data != none {
        native-data.edges = _array-set(native-data.edges, edge-source.edge, patch.data)
        changed = true
      }
      if patch.structural != none {
        let structural = patch.structural + (index: edge-source.edge)
        structural-patches.edges.push(structural)
        structural-changed = true
      }

      let source-record = edge-source.at("source", default: none)
      if source-record != none {
        let source-patch = _mapped-patch(
          source,
          _record-with-fields(source-record, edge-fields, (edge: edge-record)),
          "hedge",
          "graph.map source",
        )
        if source-patch.data != none {
          native-data.hedges = _array-set(native-data.hedges, source-record.hedge, source-patch.data)
          changed = true
        }
      }

      let sink-record = edge-source.at("sink", default: none)
      if sink-record != none {
        let sink-patch = _mapped-patch(
          sink,
          _record-with-fields(sink-record, edge-fields, (edge: edge-record)),
          "hedge",
          "graph.map sink",
        )
        if sink-patch.data != none {
          native-data.hedges = _array-set(native-data.hedges, sink-record.hedge, sink-patch.data)
          changed = true
        }
      }
    }
  }

  if structural-changed {
    graph_ = _graph-object(
      _plugin.graph_apply_structural_patches(graph-bytes(graph_), cbor.encode(structural-patches)),
      native-data,
    )
  }

  if not changed and not structural-changed {
    graph_
  } else if not structural-changed {
    _with-native-data(graph_, native-data)
  } else {
    graph_
  }
}

#let style(graph_, options) = {
  let style-record = (
    node-label: options.node-label,
    node-label-style: options.node-label-style,
    node-style: options.node-style,
    edge-label: options.edge-label,
    edge-label-style: options.edge-label-style,
    scope: options.scope,
    unit: options.unit,
  )
  map(graph_, (
    graph: graph-record => {
      let data = graph-record.at("data", default: none)
      data = if type(data) == dictionary {
        data
      } else if data == none {
        (:)
      } else {
        (value: data)
      }
      data.insert(_style-key, style-record)
      (data: data)
    },
    node: node-record => (
      statements: _node-layout-statements(node-record, style-record),
    ),
    edge: edge-record => {
      let statements = _edge-label-layout-statements(edge-record, style-record)
      if statements == none {
        none
      } else {
        (statements: statements)
      }
    },
    source: none,
    sink: none,
  ))
}


#let _apply-default-data(graph_, default-node-data, default-edge-data, default-source-data, default-sink-data) = {
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
    (
      graph: none,
      node: record => _default-data-patch(default-node-data, record),
      edge: record => _default-data-patch(default-edge-data, record),
      source: record => _default-data-patch(default-source-data, record),
      sink: record => _default-data-patch(default-sink-data, record),
    ),
  )
}
#let eval-fields(graph_, options) = {
  let eval-graph-fields = options.eval-graph-fields
  let eval-node-fields = options.eval-node-fields
  let eval-edge-fields = options.eval-edge-fields
  let eval-source-fields = options.eval-source-fields
  let eval-sink-fields = options.eval-sink-fields
  let eval-mode = options.eval-mode
  let scope = options.scope
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
    (
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
    ),
  )
}
#let parse(input, options) = {
  let default-node-data = options.default-node-data
  let default-edge-data = options.default-edge-data
  let default-source-data = options.default-source-data
  let default-sink-data = options.default-sink-data
  let eval-graph-fields = options.eval-graph-fields
  let eval-node-fields = options.eval-node-fields
  let eval-edge-fields = options.eval-edge-fields
  let eval-source-fields = options.eval-source-fields
  let eval-sink-fields = options.eval-sink-fields
  let eval-mode = options.eval-mode
  let scope = options.scope
  let graph-fields = _eval-field-list(eval-graph-fields, "graph.parse")
  let node-fields = _eval-field-list(eval-node-fields, "graph.parse")
  let edge-fields = _eval-field-list(eval-edge-fields, "graph.parse")
  let source-fields = _eval-field-list(eval-source-fields, "graph.parse")
  let sink-fields = _eval-field-list(eval-sink-fields, "graph.parse")
  let graphs = cbor(_plugin.parse_graph(bytes(input))).map(graph => _graph-object(graph, _empty-native-data()))
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
        default-node-data,
        default-edge-data,
        default-source-data,
        default-sink-data,
      )
    }
    if needs-eval {
      result = eval-fields(
        result,
        (
          eval-graph-fields: graph-fields,
          eval-node-fields: node-fields,
          eval-edge-fields: edge-fields,
          eval-source-fields: source-fields,
          eval-sink-fields: sink-fields,
          eval-mode: eval-mode,
          scope: scope,
        ),
      )
    }
    result
  })
}
#let build(options, ..items) = {
  let name = options.name
  let data = options.data
  let statements = options.statements
  let default-edge-statements = options.default-edge-statements
  let default-node-data = options.default-node-data
  let default-edge-data = options.default-edge-data
  let default-source-data = options.default-source-data
  let default-sink-data = options.default-sink-data
  let default-node-statements = options.default-node-statements
  let nodes = options.nodes
  let edges = options.edges
  let split = _split-items(items.pos(), nodes, edges)
  let keyed-nodes = _with-build-payloads(split.nodes)
  let keyed-edges = _with-build-payloads(split.edges)
  let node-keys = _node-name-map(keyed-nodes)
  _check-edge-names(keyed-edges)
  let graph-bytes_ = _plugin.graph_from_spec(cbor.encode((
    name: name,
    data: none,
    statements: _flat-statements(statements, "graph.build statements"),
    default-edge-statements: _flat-statements(default-edge-statements, "graph.build default-edge-statements"),
    default-node-statements: _flat-statements(default-node-statements, "graph.build default-node-statements"),
    nodes: keyed-nodes.map(node => _resolved-node-spec(
      node,
      node-keys,
    )),
    edges: keyed-edges.map(edge => _resolved-edge-spec(
      edge,
      node-keys,
    )),
  )))
  let native-data = _native-data-from-build(
    graph-bytes_,
    data,
    keyed-nodes,
    keyed-edges,
    default-node-data,
    default-edge-data,
    default-source-data,
    default-sink-data,
  )
  _graph-object(graph-bytes_, native-data)
}
#let node(options, ..args) = {
  let name = options.name
  let id = options.id
  let pos = options.pos
  let shift = options.shift
  let statements = options.statements
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
#let source(node, options, ..args) = {
  let name = options.name
  let id = options.id
  let statement = options.statement
  let compass = options.compass
  if args.pos().len() > 0 {
    panic("graph.source: expected only one positional node reference")
  }
  let data = _data-from-args("graph.source", args)
  _check-name(name, "graph.source")
  _check-id(id, "graph.source")
  (linnest-kind: "source", node: node, name: name, id: id, data: data, statement: statement, compass: compass)
}
#let sink(node, options, ..args) = {
  let name = options.name
  let id = options.id
  let statement = options.statement
  let compass = options.compass
  if args.pos().len() > 0 {
    panic("graph.sink: expected only one positional node reference")
  }
  let data = _data-from-args("graph.sink", args)
  _check-name(name, "graph.sink")
  _check-id(id, "graph.sink")
  (linnest-kind: "sink", node: node, name: name, id: id, data: data, statement: statement, compass: compass)
}
#let edge(options, ..args) = {
  let name = options.name
  let id = options.id
  let orientation = options.orientation
  let pos = options.pos
  let shift = options.shift
  let label-pos = options.label-pos
  let label-angle = options.label-angle
  let bend = options.bend
  let statements = options.statements
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
#let info(graph) = _info-record(graph)
#let dot(graph) = cbor(_plugin.graph_dot(graph-bytes(graph)))
#let nodes(graph, subgraph) = _node-records(graph, subgraph)
#let edges(graph, subgraph) = _edge-records(graph, subgraph)

#let _name-key(value, context_) = {
  if type(value) == label {
    _label-key(value, context_)
  } else if type(value) == str {
    value
  } else {
    panic(context_ + ": expected a Typst label or string name")
  }
}

#let _record-by-name(records, key, context_) = {
  for record in records {
    let name = record.at("name", default: none)
    if name != none and _name-key(name, context_) == key {
      return record
    }
  }
  panic(context_ + ": no record named " + repr(key))
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
  let name = record.at("name", default: none)
  if name != none and _name-key(name, context_) == key {
    _data-update(update, record, context_)
  } else {
    none
  }
}
#let node-data(graph_, name) = {
  let key = _name-key(name, "graph.node-data")
  _record-by-name(nodes(graph_, none), key, "graph.node-data").data
}
#let edge-data(graph_, name) = {
  let key = _name-key(name, "graph.edge-data")
  _record-by-name(edges(graph_, none), key, "graph.edge-data").data
}
#let update-node-data(graph_, name, update) = {
  let key = _name-key(name, "graph.update-node-data")
  if update == none {
    return graph_
  }
  let _ = _record-by-name(nodes(graph_, none), key, "graph.update-node-data")
  map(graph_, (
    graph: none,
    node: _named-data-callback(key, update, "graph.update-node-data"),
    edge: none,
    source: none,
    sink: none,
  ))
}
#let update-edge-data(graph_, name, update) = {
  let key = _name-key(name, "graph.update-edge-data")
  if update == none {
    return graph_
  }
  let _ = _record-by-name(edges(graph_, none), key, "graph.update-edge-data")
  map(graph_, (
    graph: none,
    node: none,
    edge: _named-data-callback(key, update, "graph.update-edge-data"),
    source: none,
    sink: none,
  ))
}
#let join(left, right, key) = {
  _graph-object(_plugin.graph_join_by_hedge_key(graph-bytes(left), graph-bytes(right), cbor.encode((key: key))), _empty-native-data())
}
#let cycles(graph) = cbor(_plugin.graph_cycle_basis(graph-bytes(graph)))
#let forests(graph) = cbor(_plugin.graph_spanning_forests(graph-bytes(graph)))
