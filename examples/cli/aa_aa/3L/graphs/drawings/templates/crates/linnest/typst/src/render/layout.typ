#import "../draw.typ": _overlay-style, draw
#import "../graph.typ" as graph
#import "../layout.typ": _rank-subgraph, layout as apply-layout
#import "../subgraph.typ" as subgraph

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

// Native drawing values are correlated by the Rust graph indices before any
// generated style or layout pass runs. graph.map normalizes structural fields
// and merges the remaining native values into each record's existing data.
/// Attach indexed graph, node, edge and half-edge drawing values to a graph.
/// Call this before domain preparation and layout-graph. Entries are sparse
/// graph.map patches; omitted entries preserve the graph's native data.
#let attach-elements(g, elements) = {
  let elements = _dictionary(elements, "config.elements")
  let graph-data = elements.at("graph", default: none)
  let nodes = elements.at("nodes", default: ())
  let edges = elements.at("edges", default: ())
  let hedges = elements.at("hedges", default: ())
  graph.map(
    g,
    graph: if graph-data == none { none } else { _ => graph-data },
    node: node => _indexed(nodes, node.node),
    edge: edge => _indexed(edges, edge.edge),
    source: half-edge => _indexed(hedges, half-edge.hedge),
    sink: half-edge => _indexed(hedges, half-edge.hedge),
  )
}

#let _call(value, record) = if type(value) == function { value(record) } else {
  value
}
#let _style(value, record) = {
  let value = _call(value, record)
  if value == none { (:) } else { value }
}
#let _record-style(record, key) = _style(record.at(key, default: (:)), record)
#let _record-data(record) = {
  let data = record.at("data", default: none)
  if type(data) == dictionary { data } else { (:) }
}

#let _merged-style(value, record, key) = {
  _overlay-style(_style(value, record), _record-style(record, key))
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
  if (
    (integer and type(value) != int)
      or (not integer and type(value) not in (int, float))
  ) {
    panic(
      context_
        + if integer { " must be an integer" } else { " must be a number" },
    )
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
    if (
      not data.keys().contains("minimum-length")
        or data.at("minimum-length") == none
    ) {
      none
    } else {
      (
        statements: (
          minlen: _constraint-number(
            data.at("minimum-length"),
            "EdgeDrawing.minimum_length",
            integer: true,
          ),
        ),
      )
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
      _point-difference(
        nodes.at(edge.sink.node).pos,
        nodes.at(edge.source.node).pos,
      )
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
      if (
        not data.keys().contains("label-offset")
          or data.at("label-offset") == none
      ) {
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

#let _bool-array(value) = (
  type(value) == array and value.all(item => type(item) == bool)
)

#let _subgraph-bits(g, value, context_) = {
  let value = if type(value) == function { value(g) } else { value }
  if type(value) == dictionary {
    return subgraph._impl.validate(g, value)
  }
  if not _bool-array(value) {
    panic(
      context_
        + " must be a subgraph object, boolean half-edge array, or module function",
    )
  }
  subgraph.bits(g, value)
}

#let _layout-pass(g, value) = {
  let pass = _dictionary(value, "config.layouts entry")
  if pass.keys().contains("subgraph") and pass.subgraph != none {
    pass.subgraph = _subgraph-bits(g, pass.subgraph, "config.layouts subgraph")
  }
  let constraints = _dictionary(
    pass.at("constraints", default: (:)),
    "config.layouts constraints",
  )
  let groups = if constraints.keys().contains("same-rank") {
    constraints.at("same-rank")
  } else {
    pass.at("rank-same", default: ())
  }
  if type(groups) != array {
    panic("config.layouts rank-same/constraints.same-rank must be an array")
  }
  groups += _edge-rank-same(g)
  if groups.len() > 0 {
    pass.insert("rank-same", groups.map(group => _rank-subgraph(g, group)))
  }
  if constraints.keys().contains("same-rank") {
    let _ = constraints.remove("same-rank")
    pass.insert("constraints", constraints)
  }
  pass
}

#let _layout-passes(value) = if value == none {
  ()
} else if type(value) == dictionary {
  (value,)
} else if type(value) == array {
  value
} else {
  panic(
    "config.layouts must be one layout dictionary, an array of passes, or none",
  )
}

#let _draw-subgraph(g, value) = {
  let value = if type(value) == function { value(g) } else { value }
  if value == none {
    return none
  }
  if _bool-array(value) and value.len() == 0 {
    return ()
  }
  let multiple = type(value) == array and not _bool-array(value)
  let entries = if multiple { value } else { (value,) }
  let result = entries.map(item => {
    let item = if type(item) == function { item(g) } else { item }
    if (
      type(item) == dictionary
        and item.at("linnest-kind", default: none) != "linnest-subgraph"
    ) {
      (
        item
          + (
            subgraph: _subgraph-bits(
              g,
              item.at("subgraph", default: none),
              "config.draw subgraph entry",
            ),
          )
      )
    } else {
      _subgraph-bits(g, item, "config.draw subgraph entry")
    }
  })
  if multiple { result } else { result.first() }
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
    node-label: node => _label(node-label, node, node.at(
      "name",
      default: none,
    )),
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
  if not options.keys().contains("title") {
    options.title = config.at("title", default: auto)
  }
  // graph.style stores the final callbacks. `draw` retrieves them through these
  // sentinel values without applying the element-local layer a second time.
  options.node-label = auto
  options.node-label-style = (:)
  options.node-style = (:)
  options.edge-label = none
  options.edge-label-style = (:)
  options
}

// Apply final style before layout so native labels, custom node shapes and
// label padding contribute to the measured layout dimensions.
// Element attachment and domain preparation are explicit preceding operations.
#let layout-graph(config, g) = {
  let style-options = _dictionary(
    config.at("style", default: (:)),
    "config.style",
  )
  let draw-options = _dictionary(config.at("draw", default: (:)), "config.draw")
  let passes = _layout-passes(if config.keys().contains("layouts") {
    config.at("layouts")
  } else {
    ((:),)
  })

  let defaults = (
    scope: (:),
    unit: 1.5,
    node-label: node => [$n_(#node.vid)$],
    node-label-style: (padding: 0.08),
    node-style: (:),
    edge-label: edge => [$e_(#edge.eid)$],
    edge-label-style: (:),
  )
  let prepared-style = _record-data(graph.info(g)).at(
    "linnest-style",
    default: (:),
  )
  let effective = _effective-options(
    defaults + prepared-style,
    style-options,
    draw-options,
  )
  let callbacks = _callbacks(effective)
  g = graph.style(g, ..(effective + callbacks))
  g = _apply-element-layout-constraints(g)
  let layout-defaults = _dictionary(
    config.at("layout-defaults", default: (:)),
    "config.layout-defaults",
  )
  for pass in passes {
    g = apply-layout(g, .._layout-pass(
      g,
      layout-defaults + _dictionary(pass, "config.layouts entry"),
    ))
  }
  g = _apply-label-offsets(g)
  draw(g, .._draw-options(config, g))
}

#let layout-spec(config) = {
  let path = config.at("graph-spec-path", default: none)
  if path == none {
    panic("render config requires graph-spec-path")
  }
  let g = graph.from-spec(read(path, encoding: none))
  let g = attach-elements(g, config.at("elements", default: (:)))
  layout-graph(config, g)
}
