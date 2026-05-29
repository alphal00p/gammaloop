// Internal drawing implementation. Public users should import `draw.typ`.

#import "@preview/cetz:0.5.1" as cetz
#import "../curve.typ" as curve-api
#import "../graph.typ" as graph-api
#import "../subgraph.typ" as subgraph-api

#let _style-key = "linnest-style"

#let _point-x(p) = if type(p) == array { p.at(0) } else { p.x }
#let _point-y(p) = if type(p) == array { p.at(1) } else { p.y }
#let _point(p) = (_point-x(p), _point-y(p))
#let _point-add(a, b) = (_point-x(a) + _point-x(b), _point-y(a) + _point-y(b))
#let _point-sub(a, b) = (_point-x(a) - _point-x(b), _point-y(a) - _point-y(b))
#let _point-scale(p, factor) = (_point-x(p) * factor, _point-y(p) * factor)
#let _point-lerp(a, b, t) = _point-add(a, _point-scale(_point-sub(b, a), t))
#let _route-points(endpoint) = {
  if endpoint == none {
    ()
  } else {
    endpoint.at("route-points", default: ()).map(_point)
  }
}

#let _statement-number(record, key, default: none) = {
  let value = record.statements.at(key, default: none)
  if value == none {
    default
  } else if type(value) in (int, float) {
    value
  } else {
    float(str(value))
  }
}

#let _canvas-length(unit) = {
  if type(unit) in (int, float) {
    unit * 1em
  } else {
    unit
  }
}

#let _graph-style(info) = {
  let data = info.at("data", default: none)
  if type(data) == dictionary {
    data.at(_style-key, default: (:))
  } else {
    (:)
  }
}

#let _radius-outset(radius) = {
  if type(radius) == array {
    calc.max(..radius)
  } else {
    radius
  }
}

#let _edge-geometry-defaults = (
  offset: 0,
  length: none,
  ratio: none,
  resolve-length: "min",
  accuracy: 0.001,
  optimize: true,
  offset-side: none,
)

#let _edge-routing-defaults = (
  source-anchor: auto,
  sink-anchor: auto,
  anchor-control-distance: auto,
  route: "edge-pos",
  route-points: "ignore",
)

#let _mark-defaults = (
  mark-position: "end",
  mark-orientation: "path",
  mark-direction: "forward",
)

#let _pattern-defaults = (
  pattern: none,
  pattern-amplitude: 0.1,
  pattern-wavelength: 1.0,
  pattern-phase: 0,
  pattern-samples-per-period: 16,
  pattern-coil-longitudinal-scale: 1.25,
  pattern-accuracy: 0.001,
)

#let _without-keys(style, keys) = {
  let clean = style
  for key in keys {
    if clean.keys().contains(key) {
      let _ = clean.remove(key)
    }
  }
  clean
}

#let _without-pattern-style(style) = _without-keys(style, _pattern-defaults.keys())

#let _draw-style(style) = {
  let clean = _without-pattern-style(style)
  clean = _without-keys(clean, _edge-geometry-defaults.keys())
  clean = _without-keys(clean, _edge-routing-defaults.keys())
  _without-keys(clean, _mark-defaults.keys())
}

#let _without-mark-style(style) = {
  let clean = style
  if clean.keys().contains("mark") {
    let _ = clean.remove("mark")
  }
  clean
}

#let _start-mark-entry(entry) = if type(entry) == dictionary {
  entry
} else {
  (symbol: entry)
}

#let _backward-mark(mark) = if mark == none {
  none
} else if type(mark) == dictionary {
  let clean = mark
  let entry = none
  if clean.keys().contains("end") {
    entry = clean.end
    let _ = clean.remove("end")
  } else if clean.keys().contains("symbol") {
    entry = clean.symbol
    let _ = clean.remove("symbol")
  } else if clean.keys().contains("start") {
    entry = clean.start
    let _ = clean.remove("start")
  }
  if entry == none {
    mark
  } else {
    clean + (start: _start-mark-entry(entry))
  }
} else {
  (start: _start-mark-entry(mark))
}

#let _style-value(style, key) = {
  if _edge-geometry-defaults.keys().contains(key) {
    style.at(key, default: _edge-geometry-defaults.at(key))
  } else if _edge-routing-defaults.keys().contains(key) {
    style.at(key, default: _edge-routing-defaults.at(key))
  } else if _pattern-defaults.keys().contains(key) {
    style.at(key, default: _pattern-defaults.at(key))
  } else if _mark-defaults.keys().contains(key) {
    style.at(key, default: _mark-defaults.at(key))
  } else {
    style.at(key, default: none)
  }
}

#let _edge-geometry(style) = _edge-geometry-defaults + style

#let _edge-routing(style) = _edge-routing-defaults + style

#let _pattern-style(style) = _pattern-defaults + style

#let _style-offset(style) = {
  _style-value(style, "offset")
}

#let _resolve-length-method(style) = {
  _style-value(style, "resolve-length")
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

#let _style-layers(value, data) = {
  let value = _call(value, data)
  if value == none {
    ((:),)
  } else if type(value) == array {
    value
  } else {
    (value,)
  }
}

#let _style-layer(layers, index) = {
  if index < layers.len() {
    layers.at(index)
  } else {
    none
  }
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

#let _pattern-name(style) = _style-value(style, "pattern")

#let _has-pattern(style) = {
  let pattern = _pattern-name(style)
  pattern != none and pattern != "normal" and pattern != "curve"
}

#let _has-mark(style) = style.at("mark", default: none) != none

#let _mark-position(style) = _style-value(style, "mark-position")

#let _mark-orientation(style) = _style-value(style, "mark-orientation")

#let _mark-direction(style) = _style-value(style, "mark-direction")

#let _dangling-mark-style(style) = {
  if _mark-position(style) == "center-if-dangling" {
    style + (mark-position: "center")
  } else {
    style
  }
}

#let _edge-has-source(edge) = edge.at("source-half-edge", default: none) != none

#let _edge-has-sink(edge) = edge.at("sink-half-edge", default: none) != none

#let _edge-oriented-mark-half(edge) = {
  let has-source = _edge-has-source(edge)
  let has-sink = _edge-has-sink(edge)
  if has-source and has-sink {
    if edge.at("orientation", default: "default") == "reversed" { "sink" } else { "source" }
  } else if has-source {
    "source"
  } else if has-sink {
    "sink"
  } else {
    none
  }
}

#let _orient-mark-style(style, edge, half) = {
  if style == none or not _has-mark(style) or _mark-orientation(style) != "edge" {
    style
  } else {
    let orientation = edge.at("orientation", default: "default")
    let oriented-half = _edge-oriented-mark-half(edge)
    if orientation == "undirected" or half != oriented-half {
      _without-mark-style(style)
    } else if orientation == "reversed" {
      style + (
        mark: _backward-mark(style.mark),
        mark-direction: "backward",
      )
    } else {
      style + (mark-direction: "forward")
    }
  }
}

#let _same-layer-geometry(source-style, sink-style) = {
  let same = _style-offset(source-style) == _style-offset(sink-style)
  same = (
    same and _style-value(source-style, "length") == _style-value(sink-style, "length")
  )
  same = same and _style-value(source-style, "ratio") == _style-value(sink-style, "ratio")
  same = same and _resolve-length-method(source-style) == _resolve-length-method(sink-style)
  same = (
    same and _style-value(source-style, "accuracy") == _style-value(sink-style, "accuracy")
  )
  same = (
    same and _style-value(source-style, "optimize") == _style-value(sink-style, "optimize")
  )
  same = same and _style-value(source-style, "offset-side") == _style-value(sink-style, "offset-side")
  same
}

#let _same-pattern-geometry(source-style, sink-style) = {
  let same = _has-pattern(source-style)
  same = same and _pattern-name(source-style) == _pattern-name(sink-style)
  same = same and _same-layer-geometry(source-style, sink-style)
  same = (
    same and _style-value(source-style, "pattern-amplitude") == _style-value(sink-style, "pattern-amplitude")
  )
  same = (
    same
      and _style-value(source-style, "pattern-wavelength") == _style-value(sink-style, "pattern-wavelength")
  )
  same = same and _style-value(source-style, "pattern-phase") == _style-value(sink-style, "pattern-phase")
  same = (
    same
      and _style-value(source-style, "pattern-samples-per-period")
        == _style-value(sink-style, "pattern-samples-per-period")
  )
  same = (
    same
      and _style-value(source-style, "pattern-coil-longitudinal-scale")
        == _style-value(sink-style, "pattern-coil-longitudinal-scale")
  )
  same = (
    same and _style-value(source-style, "pattern-accuracy") == _style-value(sink-style, "pattern-accuracy")
  )
  same
}

#let _same-draw-path-style(source-style, sink-style) = {
  let same-path-geometry = if _has-pattern(source-style) or _has-pattern(sink-style) {
    _same-pattern-geometry(source-style, sink-style)
  } else {
    _same-layer-geometry(source-style, sink-style)
  }
  same-path-geometry and _draw-style(source-style) == _draw-style(sink-style)
}

#let _center-outset(base-length, style, source-outset, sink-outset) = {
  let geometry = _edge-geometry(style)
  curve-api.center-outset(
    base-length,
    length: geometry.length,
    ratio: geometry.ratio,
    resolve-length: geometry.resolve-length,
    start-outset: source-outset,
    end-outset: sink-outset,
  )
}

#let _clamp(value, low, high) = calc.min(high, calc.max(low, value))

#let _node-anchor-point(box, anchor) = {
  let anchor = if anchor == auto or anchor == none { "center" } else { anchor }
  let center = box.center
  let half-width = box.width / 2
  let half-height = box.height / 2
  if anchor == "center" {
    center
  } else if anchor == "north" {
    (center.at(0), center.at(1) + half-height)
  } else if anchor == "south" {
    (center.at(0), center.at(1) - half-height)
  } else if anchor == "east" {
    (center.at(0) + half-width, center.at(1))
  } else if anchor == "west" {
    (center.at(0) - half-width, center.at(1))
  } else if anchor == "north-east" {
    (center.at(0) + half-width, center.at(1) + half-height)
  } else if anchor == "north-west" {
    (center.at(0) - half-width, center.at(1) + half-height)
  } else if anchor == "south-east" {
    (center.at(0) + half-width, center.at(1) - half-height)
  } else if anchor == "south-west" {
    (center.at(0) - half-width, center.at(1) - half-height)
  } else {
    center
  }
}

#let _anchor-control-point(anchor, point, amount) = {
  if anchor == "north" {
    (point.at(0), point.at(1) + amount)
  } else if anchor == "south" {
    (point.at(0), point.at(1) - amount)
  } else if anchor == "east" {
    (point.at(0) + amount, point.at(1))
  } else if anchor == "west" {
    (point.at(0) - amount, point.at(1))
  } else if anchor == "north-east" {
    (point.at(0) + amount, point.at(1) + amount)
  } else if anchor == "north-west" {
    (point.at(0) - amount, point.at(1) + amount)
  } else if anchor == "south-east" {
    (point.at(0) + amount, point.at(1) - amount)
  } else if anchor == "south-west" {
    (point.at(0) - amount, point.at(1) - amount)
  } else {
    point
  }
}

#let _segment-length(segment, accuracy) = {
  curve-api.length(curve-api.from-cubic(segment), accuracy: accuracy)
}

#let _visible-half-outsets(
  base-length,
  half-start,
  half-length,
  style,
  source-outset,
  sink-outset,
) = {
  let center-outset = _center-outset(
    base-length,
    style,
    source-outset,
    sink-outset,
  )
  let visible-start = source-outset + center-outset
  let visible-end = base-length - sink-outset - center-outset
  let local-start = _clamp(visible-start - half-start, 0, half-length)
  let local-end = _clamp(visible-end - half-start, 0, half-length)
  if local-end <= local-start {
    none
  } else {
    (
      start: local-start,
      end: half-length - local-end,
    )
  }
}

#let _layer(segment, style, start-outset, end-outset, label-pos, center-outset) = {
  let geometry = _edge-geometry(style)
  let side-point = if geometry.offset-side == "label" { label-pos } else { none }
  let length = geometry.length
  let ratio = geometry.ratio
  if center-outset != auto {
    start-outset = start-outset + center-outset
    end-outset = end-outset + center-outset
    length = none
    ratio = none
  }
  curve-api.layer(
    curve-api.from-cubic(segment),
    offset: geometry.offset,
    length: length,
    ratio: ratio,
    resolve-length: geometry.resolve-length,
    start-outset: start-outset,
    end-outset: end-outset,
    side-point: if side-point == none { none } else { _point(side-point) },
    accuracy: geometry.accuracy,
    optimize: geometry.optimize,
  )
}

#let _geometry-segments(
  segment,
  style,
  start-outset,
  end-outset,
  label-pos,
  center-outset,
) = {
  curve-api.segments(_layer(
    segment,
    style,
    start-outset,
    end-outset,
    label-pos,
    center-outset,
  ))
}

#let _geometry-path-segments(
  path,
  style,
  start-outset,
  end-outset,
  label-pos,
  center-outset,
) = {
  let geometry = _edge-geometry(style)
  let side-point = if geometry.offset-side == "label" { label-pos } else { none }
  let length = geometry.length
  let ratio = geometry.ratio
  if center-outset != auto {
    start-outset = start-outset + center-outset
    end-outset = end-outset + center-outset
    length = none
    ratio = none
  }
  curve-api.segments(curve-api.layer(
    path,
    offset: geometry.offset,
    length: length,
    ratio: ratio,
    resolve-length: geometry.resolve-length,
    start-outset: start-outset,
    end-outset: end-outset,
    side-point: if side-point == none { none } else { _point(side-point) },
    accuracy: geometry.accuracy,
    optimize: geometry.optimize,
  ))
}

#let _half-geometry(segment, style, outsets, label-pos) = {
  if outsets == none {
    ()
  } else {
    _geometry-segments(
      segment,
      style,
      outsets.start,
      outsets.end,
      label-pos,
      0,
    )
  }
}

#let _half-path-geometry(path, style, outsets, label-pos) = {
  if outsets == none {
    ()
  } else {
    _geometry-path-segments(
      path,
      style,
      outsets.start,
      outsets.end,
      label-pos,
      0,
    )
  }
}

#let _auto-anchor-control-distance(start, route, end) = {
  let dx = calc.abs(_point-x(end) - _point-x(start))
  let dy = calc.abs(_point-y(end) - _point-y(start))
  calc.max(0.45, calc.min(4.0, 0.18 * dx + 0.3 * dy))
}

#let _anchor-control-distance(source-style, sink-style, start, route, end) = {
  let value = _style-value(source-style, "anchor-control-distance")
  if value == auto {
    value = _style-value(sink-style, "anchor-control-distance")
  }
  if value == auto {
    _auto-anchor-control-distance(start, route, end)
  } else {
    value
  }
}

#let _anchor-points(start, anchor, route, amount, reverse: false) = {
  if anchor == auto or anchor == none or anchor == "center" {
    (start, route)
  } else if reverse {
    (route, _anchor-control-point(anchor, start, amount), start)
  } else {
    (start, _anchor-control-point(anchor, start, amount), route)
  }
}

#let _anchor-control-guide(anchor, point, fallback, amount) = {
  if anchor == auto or anchor == none or anchor == "center" {
    fallback
  } else {
    _anchor-control-point(anchor, point, amount)
  }
}

#let _anchored-cubic-route-split(start, source-anchor, route, sink-anchor, end, amount) = {
  let source-guide = _anchor-control-guide(source-anchor, start, _point-lerp(start, end, 1 / 3), amount)
  let sink-guide = _anchor-control-guide(sink-anchor, end, _point-lerp(end, start, 1 / 3), amount)
  let control-sum = _point-scale(
    _point-sub(_point-scale(route, 8), _point-add(start, end)),
    1 / 3,
  )
  let delta = _point-scale(
    _point-sub(control-sum, _point-add(source-guide, sink-guide)),
    1 / 2,
  )
  let control-start = _point-add(source-guide, delta)
  let control-end = _point-add(sink-guide, delta)
  let start-control = _point-lerp(start, control-start, 1 / 2)
  let control-control = _point-lerp(control-start, control-end, 1 / 2)
  let control-end-mid = _point-lerp(control-end, end, 1 / 2)
  let source-control-end = _point-lerp(start-control, control-control, 1 / 2)
  let sink-control-start = _point-lerp(control-control, control-end-mid, 1 / 2)
  (
    source: curve-api.cubic(start, start-control, source-control-end, route),
    sink: curve-api.cubic(route, sink-control-start, control-end-mid, end),
    curve: curve-api.cubic(start, control-start, control-end, end),
  )
}

#let _straight-through-route-split(start, route, end) = (
  source: curve-api.line(start, route),
  sink: curve-api.line(route, end),
  curve: curve-api.path(curve-api.line(start, route), curve-api.line(route, end)),
)

#let _split-edge-geometry(
  source-path,
  sink-path,
  curve,
  source-style,
  sink-style,
  source-outset,
  sink-outset,
  label-pos,
  accuracy,
) = {
  let source-length = curve-api.length(source-path, accuracy: accuracy)
  let sink-length = curve-api.length(sink-path, accuracy: accuracy)
  let base-length = source-length + sink-length
  let source-outsets = _visible-half-outsets(
    base-length,
    0,
    source-length,
    source-style,
    source-outset,
    sink-outset,
  )
  let source-geometry = _half-path-geometry(source-path, source-style, source-outsets, label-pos)
  let sink-geometry = if _same-layer-geometry(source-style, sink-style) {
    let sink-outsets = _visible-half-outsets(
      base-length,
      source-length,
      sink-length,
      source-style,
      source-outset,
      sink-outset,
    )
    _half-path-geometry(sink-path, source-style, sink-outsets, label-pos)
  } else {
    let sink-outsets = _visible-half-outsets(
      base-length,
      source-length,
      sink-length,
      sink-style,
      source-outset,
      sink-outset,
    )
    _half-path-geometry(sink-path, sink-style, sink-outsets, label-pos)
  }
  let whole-geometry = if _same-layer-geometry(source-style, sink-style) {
    let center-outset = _center-outset(base-length, source-style, source-outset, sink-outset)
    _geometry-path-segments(curve, source-style, source-outset, sink-outset, label-pos, center-outset)
  } else {
    none
  }
  (
    source: source-geometry,
    sink: sink-geometry,
    curve: curve,
    whole: whole-geometry,
  )
}

#let _anchored-edge-geometry-halves(
  edge,
  node-boxes,
  source-style,
  sink-style,
  options,
) = {
  let omega = options.omega
  let source-outset = options.source-outset
  let sink-outset = options.sink-outset
  let accuracy = options.accuracy
  let label-pos = options.label-pos
  let source-anchor = _style-value(source-style, "source-anchor")
  let sink-anchor = _style-value(sink-style, "sink-anchor")
  let source-box = node-boxes.at(edge.source.node)
  let sink-box = node-boxes.at(edge.sink.node)
  let start = _node-anchor-point(source-box, source-anchor)
  let end = _node-anchor-point(sink-box, sink-anchor)
  let route-mode = _style-value(source-style, "route")
  let route = _point(edge.pos)
  let source-route = _route-points(edge.source)
  let sink-route = _route-points(edge.sink)
  let source-start-outset = if source-anchor == auto { source-outset } else { 0 }
  let sink-end-outset = if sink-anchor == auto { sink-outset } else { 0 }
  let route-points-mode = _style-value(source-style, "route-points")
  if route-points-mode != "through" {
    route-points-mode = _style-value(sink-style, "route-points")
  }
  if route-mode == "straight-through" {
    let split = _straight-through-route-split(start, route, end)
    return _split-edge-geometry(
      split.source,
      split.sink,
      split.curve,
      source-style,
      sink-style,
      source-start-outset,
      sink-end-outset,
      label-pos,
      accuracy,
    )
  }
  if route-mode != "direct" and route-points-mode == "through" and (source-route.len() > 0 or sink-route.len() > 0) {
    let source-points = (start, ..source-route, route)
    let sink-points = (..sink-route.rev(), end)
    let split = curve-api.split-through(
      (..source-points, ..sink-points),
      omega: omega,
      accuracy: accuracy,
    )
    let source-span-count = source-points.len() - 1
    let source-path = curve-api.path(..split.parts.slice(0, source-span-count))
    let sink-path = curve-api.path(..split.parts.slice(source-span-count))
    return _split-edge-geometry(
      source-path,
      sink-path,
      split.curve,
      source-style,
      sink-style,
      source-start-outset,
      sink-end-outset,
      label-pos,
      accuracy,
    )
  }
  if route-mode in ("direct", "edge-pos") {
    let amount = _anchor-control-distance(source-style, sink-style, start, route, end)
    let split = _anchored-cubic-route-split(start, source-anchor, route, sink-anchor, end, amount)
    return _split-edge-geometry(
      split.source,
      split.sink,
      split.curve,
      source-style,
      sink-style,
      source-start-outset,
      sink-end-outset,
      label-pos,
      accuracy,
    )
  } else if route-mode == "hobby-through" {
    let amount = _anchor-control-distance(source-style, sink-style, start, route, end)
    let source-guide = _anchor-control-guide(source-anchor, start, _point-lerp(start, route, 1 / 3), amount)
    let sink-guide = _anchor-control-guide(sink-anchor, end, _point-lerp(end, route, 1 / 3), amount)
    let split = curve-api.split-through(
      (start, source-guide, route, sink-guide, end),
      omega: omega,
      accuracy: accuracy,
    )
    let source-path = curve-api.path(split.parts.at(0), split.parts.at(1))
    let sink-path = curve-api.path(split.parts.at(2), split.parts.at(3))
    return _split-edge-geometry(
      source-path,
      sink-path,
      split.curve,
      source-style,
      sink-style,
      source-start-outset,
      sink-end-outset,
      label-pos,
      accuracy,
    )
  }
  let amount = _anchor-control-distance(source-style, sink-style, start, route, end)
  let source-path = curve-api.hobby-spline(
    _anchor-points(start, source-anchor, route, amount),
    omega: omega,
    accuracy: accuracy,
  )
  let sink-path = curve-api.hobby-spline(
    _anchor-points(end, sink-anchor, route, amount, reverse: true),
    omega: omega,
    accuracy: accuracy,
  )
  (
    source: _geometry-path-segments(source-path, source-style, source-start-outset, 0, label-pos, 0),
    sink: _geometry-path-segments(sink-path, sink-style, 0, sink-end-outset, label-pos, 0),
    curve: curve-api.hobby-spline((start, route, end), omega: omega, accuracy: accuracy),
    whole: none,
  )
}

#let edge-halves(edge, nodes, options) = {
  let omega = options.omega
  let source-outset = options.source-outset
  let sink-outset = options.sink-outset
  let accuracy = options.accuracy
  if edge.source == none or edge.sink == none {
    panic("edge-halves currently requires a paired edge with source and sink")
  }

  let points = (
    _point(nodes.at(edge.source.node).pos),
    _point(edge.pos),
    _point(nodes.at(edge.sink.node).pos),
  )
  let split = curve-api.split-through(
    points,
    omega: omega,
    start-outset: source-outset,
    end-outset: sink-outset,
    accuracy: accuracy,
  )

  (
    source: split.parts.at(0),
    sink: split.parts.at(1),
    curve: split.curve,
  )
}
#let to-cetz-edge-halves(
  edge,
  nodes,
  options,
) = {
  let halves = edge-halves(
    edge,
    nodes,
    (
      omega: options.omega,
      source-outset: options.source-outset,
      sink-outset: options.sink-outset,
      accuracy: options.accuracy,
    ),
  )
  curve-api.to-cetz(halves.source, unit: options.unit, ..options.source-style)
  curve-api.to-cetz(halves.sink, unit: options.unit, ..options.sink-style)
}

#let _edge-geometry-halves(
  edge,
  nodes,
  node-boxes,
  source-style,
  sink-style,
  options,
) = {
  if (
    _style-value(source-style, "source-anchor") != auto
      or _style-value(sink-style, "sink-anchor") != auto
  ) {
    return _anchored-edge-geometry-halves(edge, node-boxes, source-style, sink-style, options)
  }

  let omega = options.omega
  let source-outset = options.source-outset
  let sink-outset = options.sink-outset
  let accuracy = options.accuracy
  let label-pos = options.label-pos
  let route-mode = _style-value(source-style, "route")
  if route-mode == "straight-through" {
    let start = _point(nodes.at(edge.source.node).pos)
    let route = _point(edge.pos)
    let end = _point(nodes.at(edge.sink.node).pos)
    let split = _straight-through-route-split(start, route, end)
    return _split-edge-geometry(
      split.source,
      split.sink,
      split.curve,
      source-style,
      sink-style,
      source-outset,
      sink-outset,
      label-pos,
      accuracy,
    )
  }
  let base = edge-halves(edge, nodes, (
    omega: omega,
    source-outset: 0,
    sink-outset: 0,
    accuracy: accuracy,
  ))
  let source-path = base.source
  let sink-path = base.sink
  let source-length = curve-api.length(source-path, accuracy: accuracy)
  let sink-length = curve-api.length(sink-path, accuracy: accuracy)
  let base-length = source-length + sink-length
  let source-outsets = _visible-half-outsets(
    base-length,
    0,
    source-length,
    source-style,
    source-outset,
    sink-outset,
  )
  let source-geometry = _half-path-geometry(source-path, source-style, source-outsets, label-pos)
  let sink-geometry = if _same-layer-geometry(source-style, sink-style) {
    let sink-outsets = _visible-half-outsets(
      base-length,
      source-length,
      sink-length,
      source-style,
      source-outset,
      sink-outset,
    )
    _half-path-geometry(sink-path, source-style, sink-outsets, label-pos)
  } else {
    let sink-outsets = _visible-half-outsets(
      base-length,
      source-length,
      sink-length,
      sink-style,
      source-outset,
      sink-outset,
    )
    _half-path-geometry(sink-path, sink-style, sink-outsets, label-pos)
  }
  let whole-geometry = if _same-layer-geometry(source-style, sink-style) {
    let center-outset = _center-outset(base-length, source-style, source-outset, sink-outset)
    _geometry-path-segments(base.curve, source-style, source-outset, sink-outset, label-pos, center-outset)
  } else {
    none
  }
  (
    source: source-geometry,
    sink: sink-geometry,
    curve: base.curve,
    whole: whole-geometry,
  )
}

#let _pattern(segment, style, phase, anchor-start, anchor-end) = {
  let pattern-style = _pattern-style(style)
  curve-api.pattern(
    curve-api.from-cubic(segment),
    pattern: pattern-style.pattern,
    amplitude: pattern-style.pattern-amplitude,
    wavelength: pattern-style.pattern-wavelength,
    phase: if phase == auto { pattern-style.pattern-phase } else { phase },
    samples-per-period: pattern-style.pattern-samples-per-period,
    coil-longitudinal-scale: pattern-style.pattern-coil-longitudinal-scale,
    anchor-start: anchor-start,
    anchor-end: anchor-end,
    accuracy: pattern-style.pattern-accuracy,
  )
}

#let _bezier-element(segment, style) = {
  cetz.draw.bezier(
    _point(segment.start),
    _point(segment.end),
    _point(segment.control-start),
    _point(segment.control-end),
    .._draw-style(style),
  )
}

#let _segments-elements(segments, style, phase, anchor-start, anchor-end) = {
  let elements = ()
  let length = 0
  if _has-pattern(style) {
    let current-phase = if phase == auto { _style-value(style, "pattern-phase") } else { phase }
    let wavelength = _style-value(style, "pattern-wavelength")
    for (index, segment) in segments.enumerate() {
      let piece = _pattern(
        segment,
        style,
        current-phase,
        anchor-start and index == 0,
        anchor-end and index == segments.len() - 1,
      )
      elements.push(curve-api.to-cetz(piece, .._draw-style(style)))
      let piece-length = _segment-length(segment, _style-value(style, "pattern-accuracy"))
      length = length + piece-length
      current-phase = current-phase + 2 * calc.pi * piece-length / wavelength
    }
  } else {
    let mark-index = if _mark-direction(style) == "backward" { 0 } else { segments.len() - 1 }
    for (index, segment) in segments.enumerate() {
      if _has-mark(style) and index == mark-index {
        elements.push(_bezier-element(segment, style))
      } else if _has-mark(style) {
        elements.push(_bezier-element(segment, _without-mark-style(style)))
      } else {
        elements.push(curve-api.to-cetz(curve-api.from-cubic(segment), .._draw-style(style)))
      }
    }
  }
  (elements: elements, length: length)
}

#let _segment-elements(segment, style, phase, anchor-start, anchor-end, label-pos) = {
  let path = _layer(segment, style, 0, 0, label-pos, auto)
  if _has-mark(style) and _mark-position(style) == "center" and not _has-pattern(style) {
    let length = curve-api.length(path, accuracy: _style-value(style, "accuracy"))
    let first = curve-api.trim(path, end-outset: length / 2, accuracy: _style-value(style, "accuracy"))
    let second = curve-api.trim(path, start-outset: length / 2, accuracy: _style-value(style, "accuracy"))
    let backward = _mark-direction(style) == "backward"
    let first-style = if backward { _without-mark-style(style) } else { style }
    let second-style = if backward { style } else { _without-mark-style(style) }
    let first-elements = _segments-elements(
      curve-api.segments(first),
      first-style,
      phase,
      anchor-start,
      false,
    ).elements
    let second-elements = _segments-elements(
      curve-api.segments(second),
      second-style,
      phase,
      false,
      anchor-end,
    ).elements
    let pieces = ()
    if backward {
      for element in first-elements {
        pieces.push(element)
      }
      for element in second-elements {
        pieces.push(element)
      }
    } else {
      for element in second-elements {
        pieces.push(element)
      }
      for element in first-elements {
        pieces.push(element)
      }
    }
    (elements: pieces, length: length)
  } else {
    _segments-elements(
      curve-api.segments(path),
      style,
      phase,
      anchor-start,
      anchor-end,
    )
  }
}

#let _pattern-edge-halves(halves, source-style, sink-style) = {
  let elements = ()
  let whole = halves.at("whole", default: none)
  if whole != none and _same-draw-path-style(source-style, sink-style) {
    return _segments-elements(whole, source-style, auto, true, true).elements
  }
  if _same-pattern-geometry(source-style, sink-style) {
    let source = _segments-elements(halves.source, source-style, auto, true, false)
    let wavelength = _style-value(source-style, "pattern-wavelength")
    let phase = _style-value(source-style, "pattern-phase")
    let sink-phase = phase + 2 * calc.pi * source.length / wavelength
    let sink-elements = _segments-elements(halves.sink, sink-style, sink-phase, false, true).elements
    if _has-mark(source-style) and not _has-mark(sink-style) {
      for element in sink-elements {
        elements.push(element)
      }
      for element in source.elements {
        elements.push(element)
      }
    } else {
      for element in source.elements {
        elements.push(element)
      }
      for element in sink-elements {
        elements.push(element)
      }
    }
  } else {
    let source-elements = _segments-elements(halves.source, source-style, auto, true, true).elements
    let sink-elements = _segments-elements(halves.sink, sink-style, auto, true, true).elements
    if _has-mark(source-style) and not _has-mark(sink-style) {
      for element in sink-elements {
        elements.push(element)
      }
      for element in source-elements {
        elements.push(element)
      }
    } else {
      for element in source-elements {
        elements.push(element)
      }
      for element in sink-elements {
        elements.push(element)
      }
    }
  }
  elements
}

#let _pattern-line(start, end, style, label-pos) = {
  _segment-elements(curve-api.line-segment(_point(start), _point(end)), style, auto, true, true, label-pos).elements
}

#let _bent-line-segment(start, end, bend) = {
  let start = _point(start)
  let end = _point(end)
  if bend == none or bend == 0 {
    return curve-api.line-segment(start, end)
  }

  let dx = _point-x(end) - _point-x(start)
  let dy = _point-y(end) - _point-y(start)
  let length = calc.sqrt(dx * dx + dy * dy)
  if length == 0 {
    return curve-api.line-segment(start, end)
  }

  let amount = calc.sin(bend) * length * 0.5
  let control = (
    (_point-x(start) + _point-x(end)) / 2 - dy / length * amount,
    (_point-y(start) + _point-y(end)) / 2 + dx / length * amount,
  )
  curve-api.segments(curve-api.quad(start, control, end)).first()
}

#let _pattern-dangling(start, end, bend, style, label-pos) = {
  _segment-elements(_bent-line-segment(start, end, bend), style, auto, true, true, label-pos).elements
}

#let _node-outset(style, node-outset) = {
  if node-outset == auto {
    _radius-outset(style.at("radius", default: 0))
  } else {
    node-outset
  }
}

#let _fit-node-radius(ctx, label, minimum, padding) = {
  if label == none {
    minimum
  } else {
    let size = measure(text(top-edge: "cap-height", bottom-edge: "bounds", label))
    let width = calc.abs(size.width / ctx.length) + 2 * padding
    let height = calc.abs(size.height / ctx.length) + 2 * padding
    calc.max(minimum, calc.sqrt(width * width + height * height) / 2)
  }
}

#let _node-radius(ctx, label, style, minimum, padding) = {
  let radius = style.at("radius", default: auto)
  if radius == auto {
    _fit-node-radius(ctx, label, minimum, padding)
  } else {
    radius
  }
}

#let _debug-level(debug) = {
  if type(debug) == bool {
    if debug { 1 } else { 0 }
  } else {
    debug
  }
}

#let _in-subgraph(hedges, half-edge) = {
  hedges != none and half-edge != none and hedges.contains(half-edge.hedge)
}

#let _node-pos(node) = {
  if node.pos == none {
    panic("draw: graph node has no position; call layout(...) or pass pos: to graph.node/build")
  }
  node.pos
}

#let _midpoint(a, b) = ((_point-x(a) + _point-x(b)) / 2, (_point-y(a) + _point-y(b)) / 2)

#let _edge-pos(edge, nodes) = {
  if edge.pos != none {
    edge.pos
  } else if edge.source != none and edge.sink != none {
    _midpoint(_node-pos(nodes.at(edge.source.node)), _node-pos(nodes.at(edge.sink.node)))
  } else {
    panic("draw: dangling graph edge has no position; call layout(...) or pass pos: to graph.edge/build")
  }
}
#let draw(graph, options) = {
  let scope = options.scope
  let unit = options.unit
  let title = options.title
  let subgraph = options.subgraph
  let debug = options.debug
  let node-radius = options.node-radius
  let node-min-radius = options.node-min-radius
  let node-label-padding = options.node-label-padding
  let node-fill = options.node-fill
  let node-stroke = options.node-stroke
  let node-outset = options.node-outset
  let node-label-style = options.node-label-style
  let node-style = options.node-style
  let node-label = options.node-label
  let draw-node = options.draw-node
  let edge-stroke = options.edge-stroke
  let edge-offset = options.edge-offset
  let edge-length = options.edge-length
  let edge-ratio = options.edge-ratio
  let edge-resolve-length = options.edge-resolve-length
  let edge-accuracy = options.edge-accuracy
  let edge-optimize = options.edge-optimize
  let source-style = options.source-style
  let sink-style = options.sink-style
  let edge-label = options.edge-label
  let edge-label-style = options.edge-label-style
  let edge-omega = options.edge-omega
  let edge-trim-accuracy = options.edge-trim-accuracy
  let padding = options.padding
  let debug-edge-radius = options.debug-edge-radius
  let debug-edge-fill = options.debug-edge-fill
  let debug-edge-stroke = options.debug-edge-stroke
  let debug-edge-label-fill = options.debug-edge-label-fill
  let subgraph-edge-style = options.subgraph-edge-style
  let subgraph-edge-underlay = options.subgraph-edge-underlay
  let info = graph-api.info(graph)
  let graph-style = _graph-style(info)
  unit = if unit == auto { graph-style.at("unit", default: 1) } else { unit }
  scope = graph-style.at("scope", default: (:)) + scope
  node-label = if node-label == auto {
    graph-style.at("node-label", default: auto)
  } else {
    node-label
  }
  let graph-node-label-style = graph-style.at("node-label-style", default: (:))
  let graph-node-style = graph-style.at("node-style", default: (:))
  let graph-edge-label = graph-style.at("edge-label", default: none)
  edge-label = if edge-label == none { graph-edge-label } else { edge-label }
  let graph-edge-label-style = graph-style.at("edge-label-style", default: (:))
  let debug-level = _debug-level(debug)
  let diagram-content = cetz.canvas(
    length: _canvas-length(unit),
    debug: debug-level >= 1,
    padding: padding,
    {
      cetz.draw.get-ctx(ctx => {
        let nodes = graph-api.nodes(graph)
        let edges = graph-api.edges(graph)
        let elements = ()
        let node-elements = ()
        let node-outsets = ()
        let node-boxes = ()
        let debug-level = _debug-level(debug)
        let subgraph-hedges = if subgraph == none { none } else { subgraph-api.hedges(subgraph) }

        for (i, v) in nodes.enumerate() {
          let pos = _node-pos(v)
          let node = v + (pos: pos)
          let record-label = _data-label(v)
          let statement-label = v.statements.at("label", default: none)
          let data-label = if record-label == none { statement-label } else { record-label }
          let node-data = (
            scope
              + v.statements
              + _data-fields(v)
              + (
                vid: i,
                node: node,
                name: node.name,
                data: v.at("data", default: none),
                label: data-label,
              )
          )
          let default-label-value = node-data.at("label", default: node.name)
          let default-label = _as-content(default-label-value)
          let label = _content(node-label, node-data, default-label)
          let node-style-value = (
            (
              radius: node-radius,
              fill: node-fill,
              stroke: node-stroke,
            )
              + _style(graph-node-style, node-data)
              + _style(node-style, node-data)
          )
          let node-style = (
            node-style-value
              + (
                radius: _node-radius(ctx, label, node-style-value, node-min-radius, node-label-padding),
              )
          )
          let node-label-draw-style = _style(graph-node-label-style, node-data) + node-label-style
          let layout-width = _statement-number(v, "layout-width", default: none)
          let layout-height = _statement-number(v, "layout-height", default: none)
          let node-width = if layout-width == none { 2 * node-style.radius } else { layout-width }
          let node-height = if layout-height == none { 2 * node-style.radius } else { layout-height }
          let box = (
            name: "n" + str(i),
            center: _point(pos),
            pos: pos,
            width: node-width,
            height: node-height,
            unit: unit,
            label: label,
            label-style: node-label-draw-style,
            style: node-style,
            radius: node-style.radius,
            node: node,
          )
          node-boxes.push(box)
          node-outsets.push(if draw-node == auto {
            _node-outset(node-style, node-outset)
          } else if node-outset == auto {
            calc.max(node-width, node-height) / 2
          } else {
            node-outset
          })

          if draw-node == auto {
            node-elements.push(cetz.draw.circle(_point(pos), name: box.name, ..node-style))
          } else {
            node-elements.push(draw-node(node-data, box))
          }
          if label != none {
            if draw-node == auto {
              node-elements.push(cetz.draw.content(
                _point(pos),
                label,
                padding: 0,
                ..node-label-draw-style,
              ))
            }
          }
        }

        for (i, e) in edges.enumerate() {
          let edge-pos = _edge-pos(e, nodes)
          let edge = e + (pos: edge-pos)
          let source-half-edge = edge.source
          let sink-half-edge = edge.sink
          let source-statement = if source-half-edge == none { none } else { source-half-edge.statement }
          let sink-statement = if sink-half-edge == none { none } else { sink-half-edge.statement }
          let source-node = if source-half-edge == none { none } else { node-boxes.at(source-half-edge.node).node }
          let sink-node = if sink-half-edge == none { none } else { node-boxes.at(sink-half-edge.node).node }
          let ext = source-half-edge == none or sink-half-edge == none
          let data = edge.statements
          let record-label = _data-label(edge)
          let statement-label = data.at("label", default: none)
          let data-label = if record-label == none { statement-label } else { record-label }
          let orientation = edge.orientation
          let label-pos = if edge.label-pos == none { edge.pos } else { edge.label-pos }
          let edge-data = (
            scope
              + data
              + _data-fields(edge)
              + (
                eid: edge.edge,
                edge: edge,
                source-statement: source-statement,
                sink-statement: sink-statement,
                source-half-edge: source-half-edge,
                sink-half-edge: sink-half-edge,
                source-node: source-node,
                sink-node: sink-node,
                data: edge.at("data", default: none),
                label: data-label,
                label-pos: label-pos,
                orientation: orientation,
                ext: ext,
              )
          )
          let source-in-subgraph = _in-subgraph(subgraph-hedges, source-half-edge)
          let sink-in-subgraph = _in-subgraph(subgraph-hedges, sink-half-edge)

          let geometry-style = _edge-geometry-defaults + (
            offset: edge-offset,
            length: edge-length,
            ratio: edge-ratio,
            resolve-length: edge-resolve-length,
            accuracy: edge-accuracy,
            optimize: edge-optimize,
          )
          let source-style-layers = _style-layers(source-style, edge-data)
          let sink-style-layers = _style-layers(sink-style, edge-data)
          let layer-count = calc.max(source-style-layers.len(), sink-style-layers.len())
          let ev-label = _content(edge-label, edge-data, _as-content(data-label))
          let edge-label-draw-style = _style(graph-edge-label-style, edge-data) + _style(edge-label-style, edge-data)

          for layer-index in range(0, layer-count) {
            let source-layer = _style-layer(source-style-layers, layer-index)
            let sink-layer = _style-layer(sink-style-layers, layer-index)
            let source-style-value = if source-layer == none {
              none
            } else {
              (stroke: edge-stroke) + geometry-style + source-layer
            }
            let sink-style-value = if sink-layer == none {
              none
            } else {
              (stroke: edge-stroke) + geometry-style + sink-layer
            }
            let source-draw-style = if source-style-value == none {
              none
            } else if source-in-subgraph and not subgraph-edge-underlay {
              source-style-value + subgraph-edge-style
            } else {
              source-style-value
            }
            source-draw-style = _orient-mark-style(source-draw-style, edge-data, "source")
            let sink-draw-style = if sink-style-value == none {
              none
            } else if sink-in-subgraph and not subgraph-edge-underlay {
              sink-style-value + subgraph-edge-style
            } else {
              sink-style-value
            }
            sink-draw-style = _orient-mark-style(sink-draw-style, edge-data, "sink")

            if source-style-value == none and sink-style-value == none {
              ()
            } else if source-half-edge != none and sink-half-edge != none {
              let source-geometry-style = if source-style-value == none { sink-style-value } else { source-style-value }
              let sink-geometry-style = if sink-style-value == none { source-style-value } else { sink-style-value }
              let halves = _edge-geometry-halves(
                edge,
                nodes,
                node-boxes,
                source-geometry-style,
                sink-geometry-style,
                (
                  omega: edge-omega,
                  source-outset: node-outsets.at(source-half-edge.node),
                  sink-outset: node-outsets.at(sink-half-edge.node),
                  accuracy: edge-trim-accuracy,
                  label-pos: label-pos,
                ),
              )
              if source-style-value != none and source-in-subgraph and subgraph-edge-underlay {
                for element in (
                  _segments-elements(
                    halves.source,
                    _without-mark-style(_without-pattern-style(source-style-value)) + subgraph-edge-style,
                    auto,
                    true,
                    true,
                  ).elements
                ) {
                  elements.push(element)
                }
              }
              if sink-style-value != none and sink-in-subgraph and subgraph-edge-underlay {
                for element in (
                  _segments-elements(
                    halves.sink,
                    _without-mark-style(_without-pattern-style(sink-style-value)) + subgraph-edge-style,
                    auto,
                    true,
                    true,
                  ).elements
                ) {
                  elements.push(element)
                }
              }
              if source-style-value != none and sink-style-value != none {
                for element in _pattern-edge-halves(halves, source-draw-style, sink-draw-style) {
                  elements.push(element)
                }
              } else if source-style-value != none {
                for element in _segments-elements(halves.source, source-draw-style, auto, true, true).elements {
                  elements.push(element)
                }
              } else if sink-style-value != none {
                for element in _segments-elements(halves.sink, sink-draw-style, auto, true, true).elements {
                  elements.push(element)
                }
              }
	            } else if source-half-edge != none {
              let source-geometry-style = if source-style-value == none { sink-style-value } else { source-style-value }
              let source-anchor = _style-value(source-geometry-style, "source-anchor")
              let line-start = if source-anchor == auto {
                curve-api.outset-point(
                  _point(_node-pos(nodes.at(source-half-edge.node))),
                  _point(edge.pos),
                  distance: node-outsets.at(source-half-edge.node),
                )
              } else {
                _node-anchor-point(node-boxes.at(source-half-edge.node), source-anchor)
              }
              let bend = edge.at("bend", default: none)
              if source-style-value != none and source-in-subgraph and subgraph-edge-underlay {
                for element in _pattern-dangling(
                  line-start,
                  edge.pos,
                  bend,
                  _without-mark-style(_without-pattern-style(_dangling-mark-style(source-style-value))) + subgraph-edge-style,
                  label-pos,
                ) {
                  elements.push(element)
                }
              }
              if source-style-value != none {
                for element in _pattern-dangling(line-start, edge.pos, bend, _dangling-mark-style(source-draw-style), label-pos) {
                  elements.push(element)
                }
              }
	            } else if sink-half-edge != none {
              let sink-geometry-style = if sink-style-value == none { source-style-value } else { sink-style-value }
              let sink-anchor = _style-value(sink-geometry-style, "sink-anchor")
              let line-end = if sink-anchor == auto {
                curve-api.outset-point(
                  _point(_node-pos(nodes.at(sink-half-edge.node))),
                  _point(edge.pos),
                  distance: node-outsets.at(sink-half-edge.node),
                )
              } else {
                _node-anchor-point(node-boxes.at(sink-half-edge.node), sink-anchor)
              }
              let bend = edge.at("bend", default: none)
              if sink-style-value != none and sink-in-subgraph and subgraph-edge-underlay {
                for element in _pattern-dangling(
                  edge.pos,
                  line-end,
                  bend,
                  _without-mark-style(_without-pattern-style(_dangling-mark-style(sink-style-value))) + subgraph-edge-style,
                  label-pos,
                ) {
                  elements.push(element)
                }
              }
              if sink-style-value != none {
                for element in _pattern-dangling(edge.pos, line-end, bend, _dangling-mark-style(sink-draw-style), label-pos) {
                  elements.push(element)
                }
              }
            }
          }

          if ev-label != none {
            elements.push(cetz.draw.content(_point(label-pos), ev-label, padding: 0, ..edge-label-draw-style))
          }

          if debug-level >= 2 {
            elements.push(cetz.draw.circle(
              _point(edge.pos),
              radius: debug-edge-radius,
              fill: debug-edge-fill,
              stroke: debug-edge-stroke,
            ))
            elements.push(cetz.draw.content(
              _point(edge.pos),
              text(size: 0.65em, fill: debug-edge-label-fill)[e#i],
              padding: 0,
            ))
          }
        }

        for element in node-elements {
          elements.push(element)
        }

        for element in elements {
          element
        }
      })
    },
  )

  if title == none {
    diagram-content
  } else {
    let shown-title = if title == auto { info.name } else { title }
    grid(
      align: center + top,
      gutter: 1em,
      [#shown-title],
      diagram-content,
    )
  }
}
