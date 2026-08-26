// Shared GammaLoop graph layout. Location-specific Linnest and edge-style
// dependencies are injected by the save-dot and documentation adapters.

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

#let _attach-elements(g, elements, graph: none) = {
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
#let _merge-style(base, patch) = {
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
      _dictionary(base.at(index, default: base.last()), "edge style layer")
        + _dictionary(layer, "edge style layer")
    ))
  } else if type(base) == array {
    base.map(layer => _dictionary(layer, "edge style layer") + patch)
  } else if type(base) == dictionary {
    base + patch
  } else {
    patch
  }
}
#let _record-drawing-style(record, key) = {
  let selected = _style(record.at("selector-style", default: (:)), record)
  _merge-style(selected, _record-style(record, key))
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
  let decorated = patches.map(patch => _merge-style(
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
  base + patches.map(patch => _merge-style(
    layers.last(),
    _dictionary(patch, "EdgeDrawing.momentum_style layer"),
  ))
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
    _merge-style(_merge-style(selected, routing), explicit)
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
    style = if explicit == none { none } else { _merge-style(style, explicit) }
  }
  let data = record.at("data", default: none)
  if type(data) == dictionary and data.keys().contains("routing") {
    style = _merge-style(style, (route: data.at("routing")))
  }
  style = _merge-style(style, _record-drawing-style(record, "edge-style"))
  style = _merge-style(style, _record-style(record, side + "-style"))
  style = _element-edge-style(style, record, momentum-arrows: momentum-arrows)
  _merge-style(style, _half-edge-style(record, side))
}

#let edge-label-style(edge) = {
  let raw-edge = edge.at("edge", default: (:))
  let statements = raw-edge.at("statements", default: (:))
  let anchor = statements.at("label-anchor", default: edge.at(
    "label-anchor",
    default: none,
  ))
  if type(anchor) == str {
    (anchor: anchor)
  } else {
    (:)
  }
}

#let _node-index-label(node) = {
  let index = node.vid
  [$n_#index$]
}

#let _node-label(value, node) = {
  let native = node.at("data", default: none)
  if type(native) == dictionary and native.keys().contains("label") {
    native.at("label")
  } else if value == auto {
    node.at("label", default: node.at("name", default: none))
  } else {
    _call(value, node)
  }
}

#let _edge-label(value, edge) = {
  let native = edge.at("data", default: none)
  if type(native) == dictionary and native.keys().contains("label") {
    native.at("label")
  } else {
    _call(value, edge)
  }
}

#let _external-label-options(options) = (
  include-particle: options.at("show-particle", default: true),
  include-index: options.at("show-edge-index", default: true),
)

#let _particle-label(
  edge,
  index,
  typst-fields,
  physics: none,
  particle-map: (:),
  include-particle: true,
  include-index: true,
) = {
  let native = edge.at("data", default: none)
  let scope = edge.fields
  if type(native) == dictionary {
    scope += native
  }
  scope += (edge: edge.edge)
  let label = if include-particle {
    let entry = physics.edge-entry(scope, map: particle-map, default: (
      label: none,
    ))
    physics.label-content(
      entry.at("label", default: none),
      scope,
      mode: typst-fields,
      map: particle-map,
      scope: scope,
    )
  } else {
    none
  }
  let index-label = if include-index { [$(p_#index)$] } else { none }
  if label == none {
    index-label
  } else if index-label == none {
    label
  } else {
    label + index-label
  }
}

#let _momentum-edge-label(
  edge,
  typst-fields,
  edge-style-options,
  physics: none,
  edge-style: none,
) = {
  let native = edge.at("data", default: none)
  if type(native) == dictionary and native.keys().contains("label") {
    return native.at("label")
  }
  if type(native) == dictionary and native.keys().contains("mode-label") {
    return native.at("mode-label")
  }
  let momentum-arrows = edge-style-options.at("momentum-arrows", default: false)
  let show-momentum-index = (
    momentum-arrows
      and edge-style-options.at(
        "show-edge-index",
        default: true,
      )
  )
  let label-options = if show-momentum-index {
    (
      edge-style-options
        + (
          show-edge-index: false,
        )
    )
  } else {
    edge-style-options
  }
  let label = edge-style.edge-label(
    edge,
    typst-fields: typst-fields,
    ..label-options,
  )
  if not show-momentum-index {
    return label
  }
  let index = edge.at("momentum-index", default: physics.edge-index(edge))
  if index == none {
    label
  } else if label == none {
    [$p_#index$]
  } else {
    label + [$(p_#index)$]
  }
}

#let _rank(ids, id) = {
  for (rank, value) in ids.enumerate() {
    if value == id {
      return rank
    }
  }
  none
}

#let _field(record, name, default: none) = {
  let native = record.at("data", default: none)
  if type(native) == dictionary {
    if native.keys().contains(name) {
      return native.at(name)
    }
    if name == "is_cut" and native.keys().contains("cut-id") {
      return native.at("cut-id")
    }
  }
  record
    .at("fields", default: record.at("statements", default: (:)))
    .at(name, default: default)
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
#let _apply-element-layout-constraints(g, graph: none) = graph.map(
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

#let _edge-rank-same(g, graph: none) = {
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
#let _apply-label-offsets(g, graph: none) = {
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

#let _subgraph-bits(g, value, context_, subgraph: none) = {
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

#let _rank-subgraph(g, value, graph: none, subgraph: none) = {
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

#let _layout-pass(g, value, graph: none, subgraph: none) = {
  let pass = _dictionary(value, "config.layouts entry")
  if pass.keys().contains("subgraph") and pass.subgraph != none {
    pass.subgraph = _subgraph-bits(
      g,
      pass.subgraph,
      "config.layouts subgraph",
      subgraph: subgraph,
    )
  }
  let groups = pass.at("rank-same", default: ())
  if type(groups) != array {
    panic("config.layouts rank-same must be an array")
  }
  groups += _edge-rank-same(g, graph: graph)
  if groups.len() > 0 {
    pass.insert("rank-same", groups.map(group => _rank-subgraph(
      g,
      group,
      graph: graph,
      subgraph: subgraph,
    )))
  }
  pass
}

#let _draw-subgraph(g, value, subgraph: none) = {
  if value == none or type(value) == bytes or type(value) == function {
    if type(value) == function { value(g) } else { value }
  } else if _bool-array(value) {
    if value.len() == 0 {
      ()
    } else {
      _subgraph-bits(g, value, "config.draw subgraph", subgraph: subgraph)
    }
  } else if type(value) == array {
    value.map(item => {
      if type(item) == dictionary and item.keys().contains("subgraph") {
        item + (subgraph: _subgraph-bits(
          g,
          item.subgraph,
          "config.draw subgraph entry",
          subgraph: subgraph,
        ),)
      } else {
        _subgraph-bits(g, item, "config.draw subgraph entry", subgraph: subgraph)
      }
    })
  } else {
    panic("config.draw subgraph must be a boolean half-edge array or an array of them")
  }
}

#let _resolved-mode(
  g,
  graph: none,
  amplitude-mode: false,
  cross-section-mode: false,
  auto-mode: false,
) = {
  if not auto-mode {
    return (
      amplitude: amplitude-mode,
      cross-section: cross-section-mode,
    )
  }
  let external = false
  let cut = false
  for edge in graph.edges(g) {
    if edge.source == none or edge.sink == none {
      external = true
      if _field(edge, "is_cut") != none {
        cut = true
      }
    }
  }
  (
    amplitude: external and not cut,
    cross-section: cut,
  )
}

// Match GammaLoop external-edge conventions with outward-facing particle
// labels. A match field pairs the two sides of a cross section without exposing
// its internal sewing tag as the momentum index.
#let autogen-external-edge-fields(
  g,
  graph: none,
  physics: none,
  particle-map: (:),
  typst-fields: "plain",
  match-field: none,
  place: true,
  y-scale: 10,
  include-particle: true,
  include-index: true,
) = {
  let left = ()
  let right = ()
  for edge in graph.edges(g) {
    let id = if match-field == none {
      edge.edge
    } else {
      _field(edge, match-field)
    }
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

    let id = if match-field == none {
      edge.edge
    } else {
      _field(edge, match-field)
    }
    if id == none {
      return none
    }
    // Incoming order defines p_i; the sewing tag only maps its partner on the
    // right back to the same rank.
    let ids = if match-field != none {
      left
    } else if side == "left" {
      left
    } else {
      right
    }
    let rank = _rank(ids, id)
    if rank == none {
      return none
    }
    let index = if match-field == none { edge.edge } else { rank }
    let generated = (
      "label-anchor": if side == "left" { "east" } else { "west" },
      "momentum-index": index,
    )
    if place and not (
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
          y: graph.start(((ids.len() - 1) / 2 - rank) * y-scale),
          dx: if side == "left" { -y-scale } else { y-scale },
        ),
      )
    }
    let explicit-label = edge.fields.at(
      "display-label",
      default: edge.fields.at("label", default: none),
    )
    if explicit-label == none {
      let label = _particle-label(
        edge,
        index,
        typst-fields,
        physics: physics,
        particle-map: particle-map,
        include-particle: include-particle,
        include-index: include-index,
      )
      // An explicit false show flag must suppress the mode-generated portion,
      // so preserve `none` as a deliberate label rather than falling through
      // to the particle-map default.
      generated.insert("mode-label", label)
    }
    // Native per-element drawing data has final precedence over mode-generated
    // side labels, momentum indices and placements.
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

#let layout(
  input,
  draw: none,
  graph: none,
  subgraph: none,
  apply-layout: none,
  physics: none,
  edge-style: none,
  diagram-options: (:),
  split-edge: true,
  scope: (:),
  columns: 1fr,
  unit: 1,
  typst-fields: "plain",
  edge-style-options: (:),
  show-node-index: false,
  amplitude-mode: false,
  cross-section-mode: false,
  additional-data: (:),
  elements: (:),
  style-options: (:),
  layout-passes: ((:),),
  auto-mode: false,
) = {
  let graphs = graph.parse(input)
  let diags = ()
  for graph-bytes in graphs {
    graph-bytes = _attach-elements(graph-bytes, elements, graph: graph)
    let mode = _resolved-mode(
      graph-bytes,
      graph: graph,
      amplitude-mode: amplitude-mode,
      cross-section-mode: cross-section-mode,
      auto-mode: auto-mode,
    )
    let is-amplitude = mode.amplitude
    let is-cross-section = mode.cross-section
    // Explicit PhysicsOptions.map entries override the generated model map,
    // which itself extends the model-neutral Linnest aliases.
    let particle-map = (
      physics.default-map
        + edge-style.map
        + edge-style-options.at("map", default: (:))
    )
    let external-label-options = _external-label-options(edge-style-options)
    if is-amplitude {
      graph-bytes = autogen-external-edge-fields(
        graph-bytes,
        graph: graph,
        physics: physics,
        particle-map: particle-map,
        typst-fields: typst-fields,
        ..external-label-options,
      )
    } else if is-cross-section {
      graph-bytes = autogen-external-edge-fields(
        graph-bytes,
        graph: graph,
        physics: physics,
        particle-map: particle-map,
        typst-fields: typst-fields,
        match-field: "is_cut",
        place: false,
        ..external-label-options,
      )
    }
    let mode-layout-options = if is-amplitude {
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
    } else if is-cross-section {
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
    let default-node-label = if show-node-index { _node-index-label } else { auto }
    let default-edge-label = edge => _momentum-edge-label(
      edge,
      typst-fields,
      edge-style-options,
      physics: physics,
      edge-style: edge-style,
    )
    let requested-style = (
      scope: scope,
      unit: unit,
      node-label: default-node-label,
      node-label-style: (padding: 0.08),
      node-style: (:),
      edge-label: default-edge-label,
      edge-label-style: edge-label-style,
    ) + style-options
    for key in (
      "scope",
      "unit",
      "node-label",
      "node-label-style",
      "node-style",
      "edge-label",
      "edge-label-style",
    ) {
      if diagram-options.keys().contains(key) {
        requested-style.insert(key, diagram-options.at(key))
      }
    }
    let node-label = node => _node-label(requested-style.node-label, node)
    let node-label-style = node => (
      _style(requested-style.node-label-style, node)
        + _record-style(node, "node-label-style")
    )
    let node-style = node => (
      _style(requested-style.node-style, node)
        + _record-drawing-style(node, "node-style")
    )
    let edge-label = edge => _edge-label(requested-style.edge-label, edge)
    let edge-label-style-value = edge => (
      _style(requested-style.edge-label-style, edge)
        + _record-style(edge, "edge-label-style")
    )
    graph-bytes = graph.style(
      graph-bytes,
      scope: requested-style.scope,
      unit: requested-style.unit,
      node-label: node-label,
      node-label-style: node-label-style,
      node-style: node-style,
      edge-label: edge-label,
      edge-label-style: edge-label-style-value,
    )
    graph-bytes = _apply-element-layout-constraints(graph-bytes, graph: graph)
    let passes = if layout-passes == none { () } else { layout-passes }
    if type(passes) != array {
      panic("config.layouts must be an array")
    }
    for pass in passes {
      graph-bytes = apply-layout(
        graph-bytes,
        ..(
          mode-layout-options
            + additional-data
            + _layout-pass(
              graph-bytes,
              pass,
              graph: graph,
              subgraph: subgraph,
            )
        ),
      )
    }
    graph-bytes = _apply-label-offsets(graph-bytes, graph: graph)
    let generated-source-style = edge => edge-style.source-style(
      edge,
      typst-fields: typst-fields,
      ..edge-style-options,
    )
    let generated-sink-style = edge => edge-style.sink-style(
      edge,
      typst-fields: typst-fields,
      ..edge-style-options,
    )
    let has-source-style = diagram-options.keys().contains("source-style")
    let has-sink-style = diagram-options.keys().contains("sink-style")
    let source-style = diagram-options.at("source-style", default: (:))
    let sink-style = diagram-options.at("sink-style", default: (:))
    let momentum-arrows = edge-style-options.at("momentum-arrows", default: false)
    let draw-subgraph = diagram-options.at("subgraph", default: none)
    let options = (
      (
        scope: requested-style.scope,
        unit: requested-style.unit,
        title: auto,
        node-label: node-label,
        node-label-style: node-label-style,
        node-style: node-style,
        source-style: edge => _endpoint-style(
          generated-source-style,
          source-style,
          has-source-style,
          edge,
          "source",
          momentum-arrows: momentum-arrows,
        ),
        sink-style: edge => _endpoint-style(
          generated-sink-style,
          sink-style,
          has-sink-style,
          edge,
          "sink",
          momentum-arrows: momentum-arrows,
        ),
        edge-label: edge-label,
        edge-label-style: edge-label-style-value,
      )
        + diagram-options
    )
    // Keep per-element drawing fields above global draw callbacks even when a
    // caller supplied source-style or sink-style in diagram-options.
    options.source-style = edge => _endpoint-style(
      generated-source-style,
      source-style,
      has-source-style,
      edge,
      "source",
      momentum-arrows: momentum-arrows,
    )
    options.sink-style = edge => _endpoint-style(
      generated-sink-style,
      sink-style,
      has-sink-style,
      edge,
      "sink",
      momentum-arrows: momentum-arrows,
    )
    if draw-subgraph != none {
      options.subgraph = _draw-subgraph(
        graph-bytes,
        draw-subgraph,
        subgraph: subgraph,
      )
    }
    // graph.style stores the final callbacks. `draw` uses these sentinel values
    // to retrieve them without applying the local layer a second time.
    options.node-label = auto
    options.node-label-style = (:)
    options.node-style = (:)
    options.edge-label = none
    options.edge-label-style = (:)
    diags.push(draw(graph-bytes, ..options))
  }
  for d in diags {
    d
  }
}

// V1 template adapter. The caller supplies native Typst dictionaries and
// arrays; this layer never parses strings or evaluates source fragments.
#let render-layout(
  config,
  draw: none,
  graph: none,
  subgraph: none,
  apply-layout: none,
  physics: none,
  edge-style: none,
  diagram-options: (:),
) = {
  let path = config.at("data-path", default: none)
  if path == none {
    panic("render config requires data-path")
  }
  let mode = config.at("mode", default: "auto")
  if not ("auto", "generic", "amplitude", "cross-section").contains(mode) {
    panic("config.mode must be auto, generic, amplitude, or cross-section")
  }
  let style-options = _dictionary(config.at("style", default: (:)), "config.style")
  let draw-options = diagram-options + _dictionary(
    config.at("draw", default: (:)),
    "config.draw",
  )
  if not draw-options.keys().contains("title") {
    draw-options.title = config.at("title", default: auto)
  }
  let physics-options = _dictionary(
    config.at("physics", default: (:)),
    "config.physics",
  )
  let edge-style-options = physics-options
  for key in ("show-node-index", "typst-fields") {
    if edge-style-options.keys().contains(key) {
      let _ = edge-style-options.remove(key)
    }
  }
  let unit = draw-options.at(
    "unit",
    default: style-options.at("unit", default: 1.5),
  )
  layout(
    read(path),
    draw: draw,
    graph: graph,
    subgraph: subgraph,
    apply-layout: apply-layout,
    physics: physics,
    edge-style: edge-style,
    diagram-options: draw-options,
    scope: draw-options.at("scope", default: style-options.at("scope", default: (:))),
    unit: unit,
    // The V1 boundary accepts native values only; executable field evaluation
    // remains available solely through the standalone Linnest layout API.
    typst-fields: "plain",
    edge-style-options: edge-style-options,
    show-node-index: physics-options.at("show-node-index", default: false),
    amplitude-mode: mode == "amplitude",
    cross-section-mode: mode == "cross-section",
    elements: config.at("elements", default: (:)),
    style-options: style-options,
    layout-passes: if config.keys().contains("layouts") {
      config.at("layouts")
    } else {
      ((:),)
    },
    auto-mode: mode == "auto",
  )
}

// Bind one location's Linnest package and particle-style adapter once. The
// returned function preserves the public GammaLoop layout signature.
#let bind-layout(
  draw: none,
  graph: none,
  subgraph: none,
  apply-layout: none,
  physics: none,
  edge-style: none,
  diagram-options: (:),
) = layout.with(
  draw: draw,
  graph: graph,
  subgraph: subgraph,
  apply-layout: apply-layout,
  physics: physics,
  edge-style: edge-style,
  diagram-options: diagram-options,
)

// Bind the same location-specific dependencies to the V1 render(config)
// contract used by Clinnet-generated entrypoints.
#let bind-render(
  draw: none,
  graph: none,
  subgraph: none,
  apply-layout: none,
  physics: none,
  edge-style: none,
  diagram-options: (:),
) = render-layout.with(
  draw: draw,
  graph: graph,
  subgraph: subgraph,
  apply-layout: apply-layout,
  physics: physics,
  edge-style: edge-style,
  diagram-options: diagram-options,
)
