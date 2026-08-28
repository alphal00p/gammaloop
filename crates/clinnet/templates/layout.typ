#import "crates/linnest/typst/src/lib.typ": (
  draw, graph, layout as apply-layout, subgraph,
)

#let _dictionary(value, context_) = if value == none {
  (:)
} else if type(value) == dictionary {
  value
} else {
  panic(context_ + " must be a dictionary")
}

#let _indexed(values, index) = if type(values) == array {
  values.at(index, default: none)
} else if type(values) == dictionary {
  values.at(str(index), default: none)
} else {
  none
}

#let _drawing-patch(value, structural-keys) = {
  if value == none {
    return none
  }
  let drawing = _dictionary(value, "config.elements entry")
  let patch = (data: drawing)
  for key in structural-keys {
    if drawing.keys().contains(key) {
      patch.insert(key, drawing.at(key))
    }
  }
  patch
}

// Native drawing values are correlated by the Rust graph indices before any
// generated style or layout pass runs. Structural positions are patched into
// the graph while the complete dictionary remains available as record data.
#let _attach-elements(g, elements) = {
  let elements = _dictionary(elements, "config.elements")
  let graph-data = elements.at("graph", default: none)
  let nodes = elements.at("nodes", default: ())
  let edges = elements.at("edges", default: ())
  let hedges = elements.at("hedges", default: ())
  graph.map(
    g,
    graph: if graph-data == none { none } else { _ => (data: graph-data) },
    node: node => _drawing-patch(
      _indexed(nodes, node.node),
      ("pos", "shift", "statements"),
    ),
    edge: edge => _drawing-patch(
      _indexed(edges, edge.edge),
      ("pos", "shift", "label-pos", "label-angle", "bend", "statements"),
    ),
    source: half-edge => _drawing-patch(
      _indexed(hedges, half-edge.hedge),
      ("statement", "port-label", "compass"),
    ),
    sink: half-edge => _drawing-patch(
      _indexed(hedges, half-edge.hedge),
      ("statement", "port-label", "compass"),
    ),
  )
}

#let _call(value, record) = if type(value) == function { value(record) } else { value }
#let _style(value, record) = {
  let value = _call(value, record)
  if value == none { (:) } else { value }
}
#let _record-style(record, key) = _style(record.at(key, default: (:)), record)
#let _overlay-style(base, patch) = {
  if patch == none or (type(patch) == dictionary and patch.len() == 0) {
    base
  } else if base == none {
    patch
  } else if type(patch) == array {
    let base = if type(base) == array {
      if base.len() > 0 { base } else { ((:),) }
    } else {
      (base,)
    }
    patch.enumerate().map(((index, layer)) => (
      _dictionary(base.at(index, default: base.last()), "style layer")
        + _dictionary(layer, "style layer")
    ))
  } else if type(base) == array {
    base.map(layer => _dictionary(layer, "style layer") + patch)
  } else if type(base) == dictionary {
    base + patch
  } else {
    patch
  }
}
#let _record-drawing-style(record, key) = _record-style(record, key)
#let _record-data(record) = {
  let data = record.at("data", default: none)
  if type(data) == dictionary { data } else { (:) }
}

#let _base-layers(style, patch, record) = {
  let patch = _call(patch, record)
  let patches = if type(patch) == array { patch } else { (patch,) }
  let layers = if type(style) == array { style } else { (style,) }
  let base = layers.at(0, default: (:))
  let decorated = patches.map(patch => _overlay-style(
    base,
    _dictionary(patch, "EdgeDrawing.decoration layer"),
  ))
  decorated + layers.slice(calc.min(1, layers.len()))
}

// Per-edge decoration patches the base layer. Multiple explicit decorations
// create multiple decorated layers while retaining any later custom layers.
#let _element-edge-style(style, record) = {
  let data = _record-data(record)
  if data.keys().contains("decoration") {
    let decoration = _call(data.at("decoration"), record)
    let patch = if decoration == none {
      (pattern: none)
    } else if type(decoration) in (dictionary, array) {
      decoration
    } else {
      (pattern: decoration)
    }
    style = _base-layers(style, patch, record)
  }
  style
}
#let _merged-style(value, record, key) = {
  _overlay-style(_style(value, record), _record-drawing-style(record, key))
}
#let _label(value, record, fallback) = {
  let data = record.at("data", default: none)
  if type(data) == dictionary and data.keys().contains("label") {
    data.at("label")
  } else if value == auto {
    fallback
  } else {
    _call(value, record)
  }
}

#let _constraint-number(value, context_, integer: false) = {
  if value == none {
    return none
  }
  if (integer and type(value) != int) or (not integer and type(value) not in (int, float)) {
    panic(context_ + if integer { " must be an integer" } else { " must be a number" })
  }
  if value < 0 {
    panic(context_ + " must be non-negative")
  }
  value
}

#let _statement-number(record, key, default: 0) = {
  let value = record.at("statements", default: (:)).at(key, default: none)
  if value == none {
    default
  } else if type(value) in (int, float) {
    value
  } else {
    float(str(value))
  }
}

// Convert first-class per-element layout constraints to the statement names
// consumed by Linnest's layered layout. Size bounds clamp graph.style's final
// measured dimensions, so explicit labels and node styles still affect layout.
#let _apply-element-layout-constraints(g) = graph.map(
  g,
  node: node => {
    let data = _record-data(node)
    let statements = (:)
    if data.keys().contains("rank") and data.at("rank") != none {
      statements.insert("layout-rank", _constraint-number(
        data.at("rank"),
        "NodeDrawing.rank",
        integer: true,
      ))
    }
    let minimum = _constraint-number(
      data.at("minimum-size", default: none),
      "NodeDrawing.minimum_size",
    )
    let maximum = _constraint-number(
      data.at("maximum-size", default: none),
      "NodeDrawing.maximum_size",
    )
    if minimum != none and maximum != none and minimum > maximum {
      panic("NodeDrawing.minimum_size must not exceed maximum_size")
    }
    if minimum != none or maximum != none {
      for key in ("layout-width", "layout-height") {
        let value = _statement-number(node, key)
        if minimum != none {
          value = calc.max(value, minimum)
        }
        if maximum != none {
          value = calc.min(value, maximum)
        }
        statements.insert(key, value)
      }
    }
    if statements.len() == 0 { none } else { (statements: statements) }
  },
  edge: edge => {
    let data = _record-data(edge)
    if not data.keys().contains("minimum-length") or data.at("minimum-length") == none {
      none
    } else {
      (statements: (minlen: _constraint-number(
        data.at("minimum-length"),
        "EdgeDrawing.minimum_length",
        integer: true,
      )))
    }
  },
  source: none,
  sink: none,
)

#let _edge-rank-same(g) = {
  let groups = ()
  for edge in graph.edges(g) {
    let data = _record-data(edge)
    if data.keys().contains("same-rank") {
      let value = data.at("same-rank")
      if value != none and type(value) != bool {
        panic("EdgeDrawing.same_rank must be a boolean")
      }
      if value == true and edge.source != none and edge.sink != none {
        groups.push((edge.source.node, edge.sink.node))
      }
    }
  }
  groups
}

#let _point-x(point) = if type(point) == array { point.at(0) } else { point.x }
#let _point-y(point) = if type(point) == array { point.at(1) } else { point.y }
#let _point-difference(left, right) = (
  x: _point-x(left) - _point-x(right),
  y: _point-y(left) - _point-y(right),
)
#let _point-length(point) = calc.sqrt(
  _point-x(point) * _point-x(point) + _point-y(point) * _point-y(point),
)

#let _edge-label-offset-point(edge, nodes, offset) = {
  let base = edge.at("label-pos", default: none)
  if base == none {
    base = edge.pos
  }
  let radial = _point-difference(base, edge.pos)
  let length = _point-length(radial)
  let direction = radial
  if length <= 1e-9 {
    let tangent = if edge.source != none and edge.sink != none {
      _point-difference(nodes.at(edge.sink.node).pos, nodes.at(edge.source.node).pos)
    } else if edge.source != none {
      _point-difference(edge.pos, nodes.at(edge.source.node).pos)
    } else if edge.sink != none {
      _point-difference(nodes.at(edge.sink.node).pos, edge.pos)
    } else {
      (x: 1, y: 0)
    }
    length = _point-length(tangent)
    direction = if length <= 1e-9 {
      (x: 0, y: 1)
    } else {
      (x: -tangent.y / length, y: tangent.x / length)
    }
  } else {
    direction = (x: radial.x / length, y: radial.y / length)
  }
  (
    _point-x(base) + offset * direction.x,
    _point-y(base) + offset * direction.y,
  )
}

// Label relaxation selects the side first. A per-edge offset then moves the
// result farther along that radial direction; coincident labels use the left
// normal of the oriented edge as a deterministic fallback.
#let _apply-label-offsets(g) = {
  let nodes = graph.nodes(g)
  graph.map(
    g,
    node: none,
    edge: edge => {
      let data = _record-data(edge)
      if not data.keys().contains("label-offset") or data.at("label-offset") == none {
        none
      } else {
        let offset = data.at("label-offset")
        if type(offset) not in (int, float) {
          panic("EdgeDrawing.label_offset must be a number")
        }
        (label-pos: _edge-label-offset-point(edge, nodes, offset))
      }
    },
    source: none,
    sink: none,
  )
}

#let _bool-array(value) = type(value) == array and value.all(item => type(item) == bool)

#let _subgraph-bits(g, value, context_) = {
  if type(value) == bytes {
    return value
  }
  if type(value) == function {
    return value(g)
  }
  if not _bool-array(value) {
    panic(context_ + " must be a boolean half-edge array or a module function")
  }
  subgraph.bits(g, value)
}

#let _rank-subgraph(g, value) = {
  if type(value) == bytes {
    return value
  }
  if type(value) == function {
    return value(g)
  }
  if type(value) != array or not value.all(item => type(item) == int) {
    panic("config.layouts rank-same entries must be node-index arrays or module functions")
  }
  let nodes = graph.nodes(g)
  let count = 0
  for edge in graph.edges(g) {
    for endpoint in (edge.source, edge.sink) {
      if endpoint != none {
        count = calc.max(count, endpoint.hedge + 1)
      }
    }
  }
  let bits = range(count).map(_ => false)
  for index in value {
    if index < 0 or index >= nodes.len() {
      panic("config.layouts rank-same node index is out of bounds")
    }
    for edge in graph.edges(g) {
      for endpoint in (edge.source, edge.sink) {
        if endpoint != none and endpoint.node == index {
          bits.at(endpoint.hedge) = true
        }
      }
    }
  }
  subgraph.bits(g, bits)
}

#let _layout-pass(g, value) = {
  let pass = _dictionary(value, "config.layouts entry")
  if pass.keys().contains("subgraph") and pass.subgraph != none {
    pass.subgraph = _subgraph-bits(g, pass.subgraph, "config.layouts subgraph")
  }
  let groups = pass.at("rank-same", default: ())
  if type(groups) != array {
    panic("config.layouts rank-same must be an array")
  }
  groups += _edge-rank-same(g)
  if groups.len() > 0 {
    pass.insert("rank-same", groups.map(group => _rank-subgraph(g, group)))
  }
  pass
}

#let _draw-subgraph(g, value) = {
  if value == none or type(value) == bytes or type(value) == function {
    if type(value) == function { value(g) } else { value }
  } else if _bool-array(value) {
    if value.len() == 0 { () } else { _subgraph-bits(g, value, "config.draw subgraph") }
  } else if type(value) == array {
    value.map(item => {
      if type(item) == dictionary and item.keys().contains("subgraph") {
        item + (subgraph: _subgraph-bits(g, item.subgraph, "config.draw subgraph entry"),)
      } else {
        _subgraph-bits(g, item, "config.draw subgraph entry")
      }
    })
  } else {
    panic("config.draw subgraph must be a boolean half-edge array or an array of them")
  }
}

#let _half-edge-style(record, side) = {
  let half-edge = record.at(side + "-half-edge", default: none)
  if half-edge == none or type(half-edge.at("data", default: none)) != dictionary {
    (:)
  } else {
    let data = half-edge.data
    let routing = (:)
    if data.keys().contains("anchor") {
      routing.insert(side + "-anchor", data.at("anchor"))
    }
    if data.keys().contains("routing") {
      routing.insert("route", data.at("routing"))
    }
    let explicit = _style(
      data.at(side + "-style", default: data.at("style", default: (:))),
      record,
    )
    _overlay-style(routing, explicit)
  }
}

#let _endpoint-style(
  configured,
  has-configured,
  record,
  side,
) = {
  let style = if has-configured { _call(configured, record) } else { (:) }
  let data = record.at("data", default: none)
  if type(data) == dictionary and data.keys().contains("routing") {
    style = _overlay-style(style, (route: data.at("routing")))
  }
  style = _overlay-style(style, _record-drawing-style(record, "edge-style"))
  style = _overlay-style(style, _record-style(record, side + "-style"))
  style = _overlay-style(style, _half-edge-style(record, side))
  _element-edge-style(style, record)
}

#let _effective-options(defaults, style-options, draw-options) = {
  let result = defaults + style-options
  for key in (
    "scope",
    "unit",
    "node-label",
    "node-label-style",
    "node-style",
    "edge-label",
    "edge-label-style",
  ) {
    if draw-options.keys().contains(key) {
      result.insert(key, draw-options.at(key))
    }
  }
  result
}

#let _callbacks(options) = {
  let node-label = options.at("node-label", default: auto)
  let node-label-style = options.at("node-label-style", default: (:))
  let node-style = options.at("node-style", default: (:))
  let edge-label = options.at("edge-label", default: none)
  let edge-label-style = options.at("edge-label-style", default: (:))
  (
    node-label: node => _label(node-label, node, node.at("name", default: none)),
    node-label-style: node => (
      _style(node-label-style, node) + _record-style(node, "node-label-style")
    ),
    node-style: node => _merged-style(node-style, node, "node-style"),
    edge-label: edge => _label(edge-label, edge, none),
    edge-label-style: edge => (
      _style(edge-label-style, edge) + _record-style(edge, "edge-label-style")
    ),
  )
}

#let _draw-options(config, g) = {
  let options = _dictionary(config.at("draw", default: (:)), "config.draw")
  if options.keys().contains("subgraph") {
    options.subgraph = _draw-subgraph(g, options.subgraph)
  }
  let has-source-style = options.keys().contains("source-style")
  let has-sink-style = options.keys().contains("sink-style")
  let source-style = options.at("source-style", default: (:))
  let sink-style = options.at("sink-style", default: (:))
  options.source-style = edge => _endpoint-style(
    source-style,
    has-source-style,
    edge,
    "source",
  )
  options.sink-style = edge => _endpoint-style(
    sink-style,
    has-sink-style,
    edge,
    "sink",
  )
  if not options.keys().contains("title") {
    options.title = config.at("title", default: auto)
  }
  // graph.style stores the final callbacks. `draw` uses these sentinel values
  // to retrieve them without applying the local layer a second time.
  options.node-label = auto
  options.node-label-style = (:)
  options.node-style = (:)
  options.edge-label = none
  options.edge-label-style = (:)
  options
}

// Apply final style before layout so native labels, custom node shapes and
// label padding contribute to the measured layout dimensions.
#let layout(config) = {
  let path = config.at("data-path", default: none)
  if path == none {
    panic("render config requires data-path")
  }
  let style-options = _dictionary(config.at("style", default: (:)), "config.style")
  let draw-options = _dictionary(config.at("draw", default: (:)), "config.draw")
  let passes = if config.keys().contains("layouts") {
    config.at("layouts")
  } else {
    ((:),)
  }
  if passes == none {
    passes = ()
  }
  if type(passes) != array {
    panic("config.layouts must be an array")
  }

  for parsed in graph.parse(read(path)) {
    let g = _attach-elements(parsed, config.at("elements", default: (:)))
    let defaults = (
      scope: (:),
      unit: 1.5,
      node-label: auto,
      node-label-style: (padding: 0.08),
      node-style: (:),
      edge-label: none,
      edge-label-style: (:),
    )
    let effective = _effective-options(defaults, style-options, draw-options)
    let callbacks = _callbacks(effective)
    g = graph.style(g, ..(effective + callbacks))
    g = _apply-element-layout-constraints(g)
    for pass in passes {
      g = apply-layout(g, .._layout-pass(g, pass))
    }
    g = _apply-label-offsets(g)
    draw(g, .._draw-options(config, g))
  }
}
