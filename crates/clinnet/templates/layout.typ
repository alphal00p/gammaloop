#import "crates/linnest/typst/src/lib.typ": (
  draw, graph, layout as apply-layout, physics, subgraph,
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
#let _record-drawing-style(record, key) = {
  let selected = _style(record.at("selector-style", default: (:)), record)
  _overlay-style(selected, _record-style(record, key))
}
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

// Per-edge decoration patches the generated particle layer. Multiple explicit
// layers create multiple particle decorations while retaining a generated
// momentum-arrow layer. Momentum style only patches that final arrow layer.
#let _element-edge-style(style, record, momentum-arrows: false) = {
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
  if not momentum-arrows or not data.keys().contains("momentum-style") {
    return style
  }
  let layers = if type(style) == array { style } else { (style,) }
  if layers.len() < 2 {
    return style
  }
  let momentum = _call(data.at("momentum-style"), record)
  let base = layers.slice(0, layers.len() - 1)
  if momentum == none {
    return if base.len() == 1 { base.first() } else { base }
  }
  let patches = if type(momentum) == array { momentum } else { (momentum,) }
  base + patches.map(patch => _overlay-style(
    layers.last(),
    _dictionary(patch, "EdgeDrawing.momentum_style layer"),
  ))
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

#let _node-index-label(node) = [$n_#node.vid$]

#let _field(record, name, default: none) = {
  let data = record.at("data", default: none)
  if type(data) == dictionary {
    if data.keys().contains(name) {
      return data.at(name)
    }
    let alias = if name == "cut-id" { "is_cut" } else if name == "is_cut" { "cut-id" } else { none }
    if alias != none and data.keys().contains(alias) {
      return data.at(alias)
    }
  }
  let fields = record.at("fields", default: record.at("statements", default: (:)))
  if fields.keys().contains(name) {
    fields.at(name)
  } else if name == "cut-id" {
    fields.at("is_cut", default: default)
  } else if name == "is_cut" {
    fields.at("cut-id", default: default)
  } else {
    default
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

#let _rank(values, value) = {
  for (index, candidate) in values.enumerate() {
    if candidate == value {
      return index
    }
  }
  none
}

#let _physics-options(config) = {
  let options = _dictionary(config.at("physics", default: (:)), "config.physics")
  for key in ("show-node-index", "typst-fields") {
    if options.keys().contains(key) {
      let _ = options.remove(key)
    }
  }
  // Extend the model-neutral aliases rather than replacing them wholesale;
  // explicit user entries are merged last and therefore have final precedence.
  options.map = physics.default-map + options.at("map", default: (:))
  options
}

#let _physics-callbacks(options) = physics.style(
  // Native V1 values are data. Executable Typst enters only through imported
  // module references, never through field evaluation.
  typst-fields: "plain",
  ..options,
)

#let _physics-edge-label(edge, options) = (_physics-callbacks(options).edge-label)(edge)

#let _particle-content(edge, options) = {
  let entry = physics.edge-entry(
    edge,
    map: options.at("map", default: physics.default-map),
    default: options.at("default", default: physics.default-edge),
  )
  physics.label-content(
    entry.at("label", default: none),
    edge,
    mode: "plain",
    map: options.at("map", default: physics.default-map),
    scope: options.at("scope", default: (:)),
  )
}

#let _join-label(left, right, separator) = if left == none {
  right
} else if right == none {
  left
} else {
  left + separator + right
}

#let _momentum-index-label(edge, options) = {
  let index = edge.at("momentum-index", default: physics.edge-index(edge))
  if index == none {
    none
  } else if options.keys().contains("edge-index-prefix") {
    let prefix = options.at("edge-index-prefix")
    if prefix == none { [#index] } else { prefix + [#index] }
  } else {
    [$p_#index$]
  }
}

#let _mode-edge-label(edge, options) = {
  let show-particle = options.at("show-particle", default: true)
  let show-edge-index = options.at("show-edge-index", default: true)
  let label = if show-particle { _particle-content(edge, options) } else { none }
  let index = if show-edge-index { _momentum-index-label(edge, options) } else { none }
  let result = _join-label(
    label,
    index,
    options.at("label-separator", default: [, ]),
  )
  let show-extra = options.at("show-momentum", default: false) or options.at(
    "show-half-edge-index",
    default: false,
  )
  if not show-extra {
    return result
  }
  let extra-options = options + (
    show-edge-index: false,
    show-particle: false,
  )
  let extra = _physics-edge-label(edge, extra-options)
  _join-label(extra, result, options.at("label-separator", default: [, ]))
}

#let _momentum-edge-label(edge, options) = {
  let arrows = options.at("momentum-arrows", default: false)
  let show-index = arrows and options.at("show-edge-index", default: true)
  let label-options = if show-index {
    options + (show-edge-index: false,)
  } else {
    options
  }
  let label = _physics-edge-label(edge, label-options)
  if show-index {
    _join-label(label, _momentum-index-label(edge, options), options.at(
      "label-separator",
      default: [, ],
    ))
  } else {
    label
  }
}

#let _generated-edge-label(edge, mode, options) = {
  let data = edge.at("data", default: none)
  if mode != "generic" and type(data) == dictionary and data.at(
    "mode-side",
    default: none,
  ) != none {
    _mode-edge-label(edge, options)
  } else {
    _momentum-edge-label(edge, options)
  }
}

#let _mode-edge-label-style(edge) = {
  let data = edge.at("data", default: none)
  let anchor = if type(data) == dictionary {
    data.at("label-anchor", default: none)
  } else {
    none
  }
  if type(anchor) == str { (anchor: anchor) } else { (:) }
}

// External-edge order fixes amplitude labels. Cross-section partners instead
// share the incoming edge's rank by matching their typed cut ID.
#let _autogen-mode-edges(g, mode) = {
  if mode == "generic" {
    return g
  }
  let match-cut = mode == "cross-section"
  let left = ()
  let right = ()
  for edge in graph.edges(g) {
    let id = if match-cut { _field(edge, "cut-id") } else { edge.edge }
    if id != none and edge.source == none and edge.sink != none {
      left.push(id)
    } else if id != none and edge.source != none and edge.sink == none {
      right.push(id)
    }
  }
  graph.map(g, edge: edge => {
    let side = if edge.source == none and edge.sink != none {
      "left"
    } else if edge.source != none and edge.sink == none {
      "right"
    } else {
      none
    }
    if side == none {
      return none
    }
    let id = if match-cut { _field(edge, "cut-id") } else { edge.edge }
    let ids = if match-cut or side == "left" { left } else { right }
    let rank = _rank(ids, id)
    if rank == none {
      return none
    }
    let generated = (
      "mode-side": side,
      "label-anchor": if side == "left" { "east" } else { "west" },
      "momentum-index": if match-cut { rank } else { edge.edge },
    )
    if mode == "amplitude" and not (
      edge.at("pos-x-set", default: false)
        or edge.at("pos-y-set", default: false)
    ) {
      let x = if side == "left" {
        graph.group("left", side: "-")
      } else {
        graph.group("right", side: "+")
      }
      generated.insert(
        "pos",
        graph.pos(
          x: x,
          y: graph.start(((ids.len() - 1) / 2 - rank) * 10),
          dx: if side == "left" { -10 } else { 10 },
        ),
      )
    }
    // Explicit per-edge drawing values stay above generated mode metadata.
    let explicit = edge.at("data", default: none)
    let patch = (:)
    for key in generated.keys() {
      if type(explicit) != dictionary or not explicit.keys().contains(key) {
        patch.insert(key, generated.at(key))
      }
    }
    patch
  })
}

#let _resolved-mode(g, requested) = {
  if requested != "auto" {
    return requested
  }
  let external = false
  let cut = false
  for edge in graph.edges(g) {
    if edge.source == none or edge.sink == none {
      external = true
      if _field(edge, "cut-id") != none {
        cut = true
      }
    }
  }
  if cut { "cross-section" } else if external { "amplitude" } else { "generic" }
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
    let selected = _style(data.at("selector-style", default: (:)), record)
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
    _overlay-style(_overlay-style(selected, routing), explicit)
  }
}

#let _endpoint-style(
  generated,
  configured,
  has-configured,
  record,
  side,
  momentum-arrows: false,
) = {
  let style = _style(generated, record)
  if has-configured {
    let explicit = _call(configured, record)
    style = if explicit == none { none } else { _overlay-style(style, explicit) }
  }
  let data = record.at("data", default: none)
  if type(data) == dictionary and data.keys().contains("routing") {
    style = _overlay-style(style, (route: data.at("routing")))
  }
  style = _overlay-style(style, _record-drawing-style(record, "edge-style"))
  style = _overlay-style(style, _record-style(record, side + "-style"))
  style = _element-edge-style(style, record, momentum-arrows: momentum-arrows)
  _overlay-style(style, _half-edge-style(record, side))
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

#let _draw-options(config, physics-callbacks, g) = {
  let options = _dictionary(config.at("draw", default: (:)), "config.draw")
  let momentum-arrows = _dictionary(
    config.at("physics", default: (:)),
    "config.physics",
  ).at("momentum-arrows", default: false)
  if options.keys().contains("subgraph") {
    options.subgraph = _draw-subgraph(g, options.subgraph)
  }
  let has-source-style = options.keys().contains("source-style")
  let has-sink-style = options.keys().contains("sink-style")
  let source-style = options.at("source-style", default: (:))
  let sink-style = options.at("sink-style", default: (:))
  options.source-style = edge => _endpoint-style(
    physics-callbacks.source-style,
    source-style,
    has-source-style,
    edge,
    "source",
    momentum-arrows: momentum-arrows,
  )
  options.sink-style = edge => _endpoint-style(
    physics-callbacks.sink-style,
    sink-style,
    has-sink-style,
    edge,
    "sink",
    momentum-arrows: momentum-arrows,
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
  let requested-mode = config.at("mode", default: "generic")
  if not ("auto", "generic", "amplitude", "cross-section").contains(requested-mode) {
    panic("config.mode must be auto, generic, amplitude, or cross-section")
  }
  let physics-value = config.at("physics", default: none)
  let physics-enabled = physics-value != none
  let physics-config = _dictionary(physics-value, "config.physics")
  let physics-options = if physics-enabled { _physics-options(config) } else { (:) }
  let physics-callbacks = if physics-enabled {
    _physics-callbacks(physics-options)
  } else {
    (source-style: (:), sink-style: (:))
  }
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
    let mode = _resolved-mode(g, requested-mode)
    g = _autogen-mode-edges(g, mode)
    let defaults = (
      scope: (:),
      unit: 1.5,
      node-label: if physics-enabled and physics-config.at(
        "show-node-index",
        default: false,
      ) { _node-index-label } else { auto },
      node-label-style: (padding: 0.08),
      node-style: (:),
      edge-label: if physics-enabled {
        edge => _generated-edge-label(edge, mode, physics-options)
      } else {
        none
      },
      edge-label-style: if physics-enabled { _mode-edge-label-style } else { (:) },
    )
    let effective = _effective-options(defaults, style-options, draw-options)
    let callbacks = _callbacks(effective)
    g = graph.style(g, ..(effective + callbacks))
    g = _apply-element-layout-constraints(g)
    let mode-layout = if mode == "amplitude" {
      (
        k-spring: 4.5,
        eps: 1e-7,
        step: 0.6,
        gamma-dangling: 2.3,
        label-length-scale: 1.2,
        label-steps: 100,
        directional-force: 4.5,
        label-layout: "dangling-tangent",
      )
    } else if mode == "cross-section" {
      (
        length-scale: 0.4,
        z-spring-growth: 1.01,
        label-length-scale: 1.2,
        label-steps: 100,
        label-layout: "dangling-tangent",
      )
    } else {
      (:)
    }
    for pass in passes {
      g = apply-layout(g, ..(mode-layout + _layout-pass(g, pass)))
    }
    g = _apply-label-offsets(g)
    draw(g, .._draw-options(config, physics-callbacks, g))
  }
}
