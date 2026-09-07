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
#let _point-length(p) = calc.sqrt(
  _point-x(p) * _point-x(p) + _point-y(p) * _point-y(p),
)
#let _point-distance(a, b) = _point-length(_point-sub(b, a))
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
  shift: 0,
  accuracy: 0.001,
  optimize: true,
  offset-side: none,
  split-gap: 0,
)

#let _edge-crossing-defaults = (
  crossing-under: none,
  crossing-gap: 0.55,
)

#let _edge-label-defaults = (
  label: none,
  label-style: (:),
  label-gap: 0.15,
  label-side: auto,
)

#let _edge-routing-defaults = (
  source-anchor: auto,
  sink-anchor: auto,
  anchor-control-distance: auto,
  route: "edge-pos",
  route-points: "ignore",
  dangling-tangent: auto,
)

#let _mark-defaults = (
  mark-position: "end",
  mark-shift: 0,
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

#let _without-pattern-style(style) = _without-keys(
  style,
  _pattern-defaults.keys(),
)

#let _draw-style(style) = {
  let clean = _without-pattern-style(style)
  clean = _without-keys(clean, _edge-geometry-defaults.keys())
  clean = _without-keys(clean, _edge-crossing-defaults.keys())
  clean = _without-keys(clean, _edge-routing-defaults.keys())
  clean = _without-keys(clean, _edge-label-defaults.keys())
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
  } else if _edge-crossing-defaults.keys().contains(key) {
    style.at(key, default: _edge-crossing-defaults.at(key))
  } else if _edge-routing-defaults.keys().contains(key) {
    style.at(key, default: _edge-routing-defaults.at(key))
  } else if _edge-label-defaults.keys().contains(key) {
    style.at(key, default: _edge-label-defaults.at(key))
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

#let _style-dictionary(value) = if type(value) == dictionary {
  value
} else {
  panic("style layer must be a dictionary")
}

#let _overlay-style(base, patch) = {
  if (
    patch == none
      or patch == auto
      or (type(patch) == dictionary and patch.len() == 0)
  ) {
    return base
  }
  if base == none {
    return patch
  }
  if type(patch) == array {
    let bases = if type(base) == array {
      if base.len() == 0 { ((:),) } else { base }
    } else {
      (base,)
    }
    return patch
      .enumerate()
      .map(((index, layer)) => (
        _style-dictionary(bases.at(index, default: bases.last()))
          + _style-dictionary(layer)
      ))
  }
  if type(base) == array {
    return base.map(layer => (
      _style-dictionary(layer) + _style-dictionary(patch)
    ))
  }
  _style-dictionary(base) + _style-dictionary(patch)
}

#let _edge-style-layers(default, endpoint, data) = {
  let style = _call(default, data)
  if style == none {
    return ()
  }
  if style == auto {
    style = (:)
  }
  let local = _call(data.at("edge-style", default: auto), data)
  if local == none {
    return ()
  }
  style = _overlay-style(style, local)
  style = _overlay-style(style, _call(endpoint, data))
  if type(style) == array { style } else { (style,) }
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

#let _mark-ratio(style) = {
  let position = _mark-position(style)
  if position == "center" {
    0.5
  } else if type(position) in (int, float) {
    if position < 0 or position > 1 {
      panic("draw: numeric mark-position must be between 0 and 1")
    }
    position
  } else {
    none
  }
}

#let _positioned-mark-style(style) = if (
  _mark-position(style) == "center-if-dangling"
) {
  (
    style
      + (
        mark-position: if _style-value(style, "mark-direction") == "backward" {
          0
        } else { 1 },
      )
  )
} else {
  style
}

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
    if edge.at("orientation", default: "default") == "reversed" {
      "sink"
    } else { "source" }
  } else if has-source {
    "source"
  } else if has-sink {
    "sink"
  } else {
    none
  }
}

#let _orient-mark-style(style, edge, half) = {
  if (
    style == none or not _has-mark(style) or _mark-orientation(style) != "edge"
  ) {
    style
  } else {
    let orientation = edge.at("orientation", default: "default")
    let oriented-half = _edge-oriented-mark-half(edge)
    if orientation == "undirected" or half != oriented-half {
      _without-mark-style(style)
    } else if orientation == "reversed" {
      (
        style
          + (
            mark: _backward-mark(style.mark),
            mark-direction: "backward",
          )
      )
    } else {
      style + (mark-direction: "forward")
    }
  }
}

#let _same-layer-geometry(source-style, sink-style) = {
  let same = _style-offset(source-style) == _style-offset(sink-style)
  same = (
    same
      and _style-value(source-style, "length")
        == _style-value(sink-style, "length")
  )
  same = (
    same
      and _style-value(source-style, "ratio")
        == _style-value(sink-style, "ratio")
  )
  same = (
    same
      and _resolve-length-method(source-style)
        == _resolve-length-method(sink-style)
  )
  same = (
    same
      and _style-value(source-style, "shift")
        == _style-value(sink-style, "shift")
  )
  same = (
    same
      and _style-value(source-style, "accuracy")
        == _style-value(sink-style, "accuracy")
  )
  same = (
    same
      and _style-value(source-style, "optimize")
        == _style-value(sink-style, "optimize")
  )
  same = (
    same
      and _style-value(source-style, "offset-side")
        == _style-value(sink-style, "offset-side")
  )
  same = (
    same
      and _style-value(source-style, "split-gap")
        == _style-value(sink-style, "split-gap")
  )
  same
}

#let _same-pattern-geometry(source-style, sink-style) = {
  let same = _has-pattern(source-style)
  same = same and _pattern-name(source-style) == _pattern-name(sink-style)
  same = same and _same-layer-geometry(source-style, sink-style)
  same = (
    same
      and _style-value(source-style, "pattern-amplitude")
        == _style-value(sink-style, "pattern-amplitude")
  )
  same = (
    same
      and _style-value(source-style, "pattern-wavelength")
        == _style-value(sink-style, "pattern-wavelength")
  )
  same = (
    same
      and _style-value(source-style, "pattern-phase")
        == _style-value(sink-style, "pattern-phase")
  )
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
    same
      and _style-value(source-style, "pattern-accuracy")
        == _style-value(sink-style, "pattern-accuracy")
  )
  same
}

#let _same-draw-path-style(source-style, sink-style) = {
  let same-path-geometry = if (
    _has-pattern(source-style) or _has-pattern(sink-style)
  ) {
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
  let shift = _clamp(
    _style-value(style, "shift"),
    -center-outset,
    center-outset,
  )
  let visible-start = source-outset + center-outset + shift
  let visible-end = base-length - sink-outset - center-outset + shift
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

#let _path-layer(
  path,
  style,
  start-outset,
  end-outset,
  label-pos,
  center-outset,
) = {
  let geometry = _edge-geometry(style)
  let side-point = if geometry.offset-side == "label" { label-pos } else {
    none
  }
  let length = geometry.length
  let ratio = geometry.ratio
  if center-outset != auto {
    let shift = _clamp(geometry.shift, -center-outset, center-outset)
    start-outset = start-outset + center-outset + shift
    end-outset = end-outset + center-outset - shift
    length = none
    ratio = none
    geometry.shift = 0
  }
  curve-api.layer(
    path,
    offset: geometry.offset,
    length: length,
    ratio: ratio,
    resolve-length: geometry.resolve-length,
    shift: geometry.shift,
    start-outset: start-outset,
    end-outset: end-outset,
    side-point: if side-point == none { none } else { _point(side-point) },
    accuracy: geometry.accuracy,
    optimize: geometry.optimize,
  )
}

#let _layer(
  segment,
  style,
  start-outset,
  end-outset,
  label-pos,
  center-outset,
) = {
  _path-layer(
    curve-api.from-cubic(segment),
    style,
    start-outset,
    end-outset,
    label-pos,
    center-outset,
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
  curve-api.segments(_path-layer(
    path,
    style,
    start-outset,
    end-outset,
    label-pos,
    center-outset,
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

#let _merged-half-outsets(outsets, half-length, start: 0, end: 0) = {
  if outsets == none {
    return none
  }
  let start = calc.max(outsets.start, start)
  let end = calc.max(outsets.end, end)
  if start + end >= half-length { none } else { (start: start, end: end) }
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

#let _route-aware-anchor-amount(point, amount, route-points, fallback) = {
  let target = if route-points.len() > 0 { route-points.first() } else {
    fallback
  }
  let distance = _point-distance(point, target)
  if distance == 0 {
    amount
  } else {
    calc.min(amount, 0.55 * distance)
  }
}

#let _anchored-cubic-route-split(
  start,
  source-anchor,
  route,
  sink-anchor,
  end,
  amount,
) = {
  let source-guide = _anchor-control-guide(
    source-anchor,
    start,
    _point-lerp(start, route, 1 / 3),
    amount,
  )
  let sink-guide = _anchor-control-guide(
    sink-anchor,
    end,
    _point-lerp(end, route, 1 / 3),
    amount,
  )
  let route-direction = _point-sub(end, start)
  let route-direction-length = _point-length(route-direction)
  let route-handle = if route-direction-length == 0 {
    (0, 0)
  } else {
    let max-handle = (
      calc.min(_point-distance(start, route), _point-distance(route, end)) / 3
    )
    _point-scale(
      route-direction,
      calc.min(amount, max-handle) / route-direction-length,
    )
  }
  let source-route-guide = _point-sub(route, route-handle)
  let sink-route-guide = _point-add(route, route-handle)
  let source = curve-api.cubic(start, source-guide, source-route-guide, route)
  let sink = curve-api.cubic(route, sink-route-guide, sink-guide, end)
  (
    source: source,
    sink: sink,
    curve: curve-api.path(source, sink),
  )
}

#let _straight-through-route-split(start, route, end) = (
  source: curve-api.line(start, route),
  sink: curve-api.line(route, end),
  curve: curve-api.path(curve-api.line(start, route), curve-api.line(
    route,
    end,
  )),
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
  let source-split-outset = calc.min(
    source-length,
    calc.max(0, _style-value(source-style, "split-gap")) / 2,
  )
  let sink-split-outset = calc.min(
    sink-length,
    calc.max(0, _style-value(sink-style, "split-gap")) / 2,
  )
  let source-outsets = _merged-half-outsets(
    _visible-half-outsets(
      base-length,
      0,
      source-length,
      source-style,
      source-outset,
      sink-outset,
    ),
    source-length,
    end: source-split-outset,
  )
  let source-geometry = _half-path-geometry(
    source-path,
    source-style,
    source-outsets,
    label-pos,
  )
  let sink-geometry = if _same-layer-geometry(source-style, sink-style) {
    let sink-outsets = _merged-half-outsets(
      _visible-half-outsets(
        base-length,
        source-length,
        sink-length,
        source-style,
        source-outset,
        sink-outset,
      ),
      sink-length,
      start: sink-split-outset,
    )
    _half-path-geometry(sink-path, source-style, sink-outsets, label-pos)
  } else {
    let sink-outsets = _merged-half-outsets(
      _visible-half-outsets(
        base-length,
        source-length,
        sink-length,
        sink-style,
        source-outset,
        sink-outset,
      ),
      sink-length,
      start: sink-split-outset,
    )
    _half-path-geometry(sink-path, sink-style, sink-outsets, label-pos)
  }
  let split-gap = source-split-outset + sink-split-outset
  let whole-geometry = if (
    _same-layer-geometry(source-style, sink-style) and split-gap == 0
  ) {
    let center-outset = _center-outset(
      base-length,
      source-style,
      source-outset,
      sink-outset,
    )
    _geometry-path-segments(
      curve,
      source-style,
      source-outset,
      sink-outset,
      label-pos,
      center-outset,
    )
  } else {
    none
  }
  (
    source: source-geometry,
    sink: sink-geometry,
    curve: curve,
    whole: whole-geometry,
    split-gap: split-gap,
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
  let source-start-outset = if source-anchor == auto { source-outset } else {
    0
  }
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
  if (
    route-mode != "direct"
      and route-points-mode == "through"
      and (source-route.len() > 0 or sink-route.len() > 0)
  ) {
    let amount = _anchor-control-distance(
      source-style,
      sink-style,
      start,
      route,
      end,
    )
    let source-amount = _route-aware-anchor-amount(
      start,
      amount,
      source-route,
      route,
    )
    let sink-amount = _route-aware-anchor-amount(end, amount, sink-route, route)
    let source-guide = _anchor-control-guide(
      source-anchor,
      start,
      _point-lerp(start, route, 1 / 3),
      source-amount,
    )
    let sink-guide = _anchor-control-guide(
      sink-anchor,
      end,
      _point-lerp(end, route, 1 / 3),
      sink-amount,
    )
    let source-points = (start, source-guide, ..source-route, route)
    let sink-points = (..sink-route.rev(), sink-guide, end)
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
    let amount = _anchor-control-distance(
      source-style,
      sink-style,
      start,
      route,
      end,
    )
    let split = _anchored-cubic-route-split(
      start,
      source-anchor,
      route,
      sink-anchor,
      end,
      amount,
    )
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
    let amount = _anchor-control-distance(
      source-style,
      sink-style,
      start,
      route,
      end,
    )
    let source-guide = _anchor-control-guide(
      source-anchor,
      start,
      _point-lerp(start, route, 1 / 3),
      amount,
    )
    let sink-guide = _anchor-control-guide(
      sink-anchor,
      end,
      _point-lerp(end, route, 1 / 3),
      amount,
    )
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
  let amount = _anchor-control-distance(
    source-style,
    sink-style,
    start,
    route,
    end,
  )
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
  let halves = _split-edge-geometry(
    source-path,
    sink-path,
    curve-api.hobby-spline(
      (start, route, end),
      omega: omega,
      accuracy: accuracy,
    ),
    source-style,
    sink-style,
    source-start-outset,
    sink-end-outset,
    label-pos,
    accuracy,
  )
  // The fallback halves contain anchor guides that are absent from `curve`.
  halves.whole = none
  halves
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
  let split-outset = calc.max(0, options.at("split-gap", default: 0)) / 2
  let source = split.parts.at(0)
  let sink = split.parts.at(1)
  let source-split-outset = 0
  let sink-split-outset = 0
  if split-outset > 0 {
    let source-length = curve-api.length(source, accuracy: accuracy)
    let sink-length = curve-api.length(sink, accuracy: accuracy)
    source-split-outset = calc.min(split-outset, source-length)
    sink-split-outset = calc.min(split-outset, sink-length)
    source = if source-split-outset == source-length {
      curve-api.path()
    } else {
      curve-api.trim(
        source,
        end-outset: source-split-outset,
        accuracy: accuracy,
      )
    }
    sink = if sink-split-outset == sink-length {
      curve-api.path()
    } else {
      curve-api.trim(sink, start-outset: sink-split-outset, accuracy: accuracy)
    }
  }

  (
    source: source,
    sink: sink,
    curve: split.curve,
    split-gap: source-split-outset + sink-split-outset,
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
      split-gap: options.at("split-gap", default: 0),
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
    return _anchored-edge-geometry-halves(
      edge,
      node-boxes,
      source-style,
      sink-style,
      options,
    )
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
    split-gap: 0,
    accuracy: accuracy,
  ))
  _split-edge-geometry(
    base.source,
    base.sink,
    base.curve,
    source-style,
    sink-style,
    source-outset,
    sink-outset,
    label-pos,
    accuracy,
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
    let current-phase = if phase == auto {
      _style-value(style, "pattern-phase")
    } else { phase }
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
      let piece-length = _segment-length(segment, _style-value(
        style,
        "pattern-accuracy",
      ))
      length = length + piece-length
      current-phase = current-phase + 2 * calc.pi * piece-length / wavelength
    }
  } else {
    let mark-index = if _mark-direction(style) == "backward" { 0 } else {
      segments.len() - 1
    }
    for (index, segment) in segments.enumerate() {
      if _has-mark(style) and index == mark-index {
        elements.push(_bezier-element(segment, style))
      } else if _has-mark(style) {
        elements.push(_bezier-element(segment, _without-mark-style(style)))
      } else {
        elements.push(curve-api.to-cetz(
          curve-api.from-cubic(segment),
          .._draw-style(style),
        ))
      }
    }
  }
  (elements: elements, length: length)
}

#let _segments-path(segments) = {
  curve-api.path(..segments.map(curve-api.from-cubic))
}

#let _positioned-mark-elements(path, style, phase, anchor-start, anchor-end) = {
  let mark-ratio = _mark-ratio(style)
  let length = curve-api.length(path, accuracy: _style-value(style, "accuracy"))
  let accuracy = _style-value(style, "accuracy")
  let inset = calc.min(length / 2, calc.max(length * 1e-9, accuracy))
  let mark-at = calc.max(inset, calc.min(
    length - inset,
    mark-ratio * length + _style-value(style, "mark-shift"),
  ))
  let first = curve-api.trim(
    path,
    end-outset: length - mark-at,
    accuracy: accuracy,
  )
  let second = curve-api.trim(path, start-outset: mark-at, accuracy: accuracy)
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
  let pieces = if backward {
    first-elements + second-elements
  } else {
    second-elements + first-elements
  }
  (elements: pieces, length: length)
}

#let _mark-carrier-elements(path, style) = {
  if style == none or not _has-mark(style) {
    return ()
  }
  let carrier = _without-pattern-style(style)
  carrier = _without-keys(carrier, _edge-geometry-defaults.keys())
  carrier = _without-keys(carrier, _edge-crossing-defaults.keys())
  carrier = _without-keys(carrier, _edge-routing-defaults.keys())
  carrier = _without-keys(carrier, _edge-label-defaults.keys())
  carrier = carrier + (stroke: none)
  carrier = _positioned-mark-style(carrier)
  let mark-ratio = _mark-ratio(carrier)
  if mark-ratio == none {
    _segments-elements(
      curve-api.segments(path),
      carrier,
      auto,
      true,
      true,
    ).elements
  } else {
    _positioned-mark-elements(path, carrier, auto, true, true).elements
  }
}

#let _derived-path-elements(path, style, phase, anchor-start, anchor-end) = {
  let segments = curve-api.segments(path)
  if segments.len() == 0 {
    return (elements: (), length: 0)
  }
  let positioned-style = _positioned-mark-style(style)
  let mark-ratio = _mark-ratio(positioned-style)
  if _has-mark(style) and mark-ratio != none {
    if _has-pattern(style) {
      let painted = _segments-elements(
        segments,
        _without-mark-style(style),
        phase,
        anchor-start,
        anchor-end,
      )
      (
        elements: painted.elements
          + _mark-carrier-elements(path, positioned-style),
        length: painted.length,
      )
    } else {
      _positioned-mark-elements(
        path,
        positioned-style,
        phase,
        anchor-start,
        anchor-end,
      )
    }
  } else {
    _segments-elements(segments, style, phase, anchor-start, anchor-end)
  }
}

#let _derived-segments-elements(
  segments,
  style,
  phase,
  anchor-start,
  anchor-end,
) = {
  if segments.len() == 0 {
    (elements: (), length: 0)
  } else {
    _derived-path-elements(
      _segments-path(segments),
      style,
      phase,
      anchor-start,
      anchor-end,
    )
  }
}

#let _path-elements(path, style, phase, anchor-start, anchor-end, label-pos) = {
  _derived-path-elements(
    _path-layer(path, style, 0, 0, label-pos, auto),
    style,
    phase,
    anchor-start,
    anchor-end,
  )
}

#let _center-mark-style(source-style, sink-style) = {
  for style in (source-style, sink-style) {
    if (
      style != none
        and _has-mark(style)
        and (
          _mark-ratio(style) != none
            or _mark-position(style) == "center-if-dangling"
        )
    ) {
      return style
    }
  }
  none
}

#let _paired-carrier-style(path, halves, style) = {
  if style == none or _mark-position(style) != "center-if-dangling" {
    return style
  }
  let accuracy = _style-value(style, "accuracy")
  let total = curve-api.length(path, accuracy: accuracy)
  let source-length = if halves.source.len() == 0 {
    0
  } else {
    curve-api.length(_segments-path(halves.source), accuracy: accuracy)
  }
  (
    style
      + (
        mark-position: if total <= accuracy { 0.5 } else {
          source-length / total
        },
      )
  )
}

#let _path-mid-frame(path, accuracy) = {
  let segments = curve-api.segments(path)
  if segments.len() == 0 {
    return none
  }
  let total = curve-api.length(path, accuracy: accuracy)
  if total <= accuracy {
    let segment = segments.first()
    return (
      point: segment.start,
      tangent: curve-api.cubic-tangent(segment, 0),
    )
  }
  let prefix = curve-api.trim(path, end-outset: total / 2, accuracy: accuracy)
  let prefix-segments = curve-api.segments(prefix)
  let segment = if prefix-segments.len() == 0 {
    segments.first()
  } else {
    prefix-segments.last()
  }
  (
    point: segment.end,
    tangent: curve-api.cubic-tangent(segment, 1),
  )
}

#let _label-side(style, frame, label-pos, label-origin) = {
  let side = _style-value(style, "label-side")
  if side == auto {
    if label-pos == none or label-origin == none {
      1
    } else {
      let toward = _point-sub(_point(label-pos), _point(label-origin))
      let cross = (
        _point-x(frame.tangent) * _point-y(toward)
          - _point-y(frame.tangent) * _point-x(toward)
      )
      if cross < 0 { -1 } else { 1 }
    }
  } else if side in ("left", "+") {
    1
  } else if side in ("right", "-") {
    -1
  } else if type(side) in (int, float) {
    if side < 0 { -1 } else { 1 }
  } else {
    panic("draw: label-side must be auto, left, right, +, -, or a number")
  }
}

#let _layer-label(style, data) = if style == none {
  none
} else {
  _as-content(_call(_style-value(style, "label"), data))
}

#let _has-layer-label(layers, data) = layers.any(style => (
  _layer-label(style, data) != none
))

#let _layer-label-element(ctx, path, style, label-pos, data) = {
  let label = _layer-label(style, data)
  if label == none {
    return none
  }
  let edge = data.at("edge", default: none)
  let label-origin = if edge == none { none } else {
    edge.at("pos", default: none)
  }
  let frame = _path-mid-frame(path, _style-value(style, "accuracy"))
  if frame == none {
    return none
  }
  let tangent-length = _point-length(frame.tangent)
  let normal = if tangent-length <= 1e-9 {
    (0, 1)
  } else {
    (
      -_point-y(frame.tangent) / tangent-length,
      _point-x(frame.tangent) / tangent-length,
    )
  }
  let size = measure(label)
  let width = calc.abs(size.width / ctx.length)
  let height = calc.abs(size.height / ctx.length)
  let clearance = (
    calc.abs(_point-x(normal)) * width / 2
      + calc.abs(_point-y(normal)) * height / 2
  )
  let side = _label-side(style, frame, label-pos, label-origin)
  let position = _point-add(
    frame.point,
    _point-scale(
      normal,
      side * (clearance + calc.max(0, _style-value(style, "label-gap"))),
    ),
  )
  let label-style = _style(_style-value(style, "label-style"), data)
  cetz.draw.content(_point(position), label, padding: 0, ..label-style)
}

#let _paired-layer-label-element(
  ctx,
  halves,
  source-style,
  sink-style,
  label-pos,
  data,
) = {
  let source-label = _layer-label(source-style, data)
  let sink-label = _layer-label(sink-style, data)
  let style = if source-label != none { source-style } else { sink-style }
  if style == none {
    return none
  }
  let segments = if source-label != none and sink-label != none {
    if halves.whole == none { halves.source + halves.sink } else {
      halves.whole
    }
  } else if source-label != none {
    halves.source
  } else {
    halves.sink
  }
  if segments.len() == 0 {
    none
  } else {
    _layer-label-element(ctx, _segments-path(segments), style, label-pos, data)
  }
}

#let _segment-elements(
  segment,
  style,
  phase,
  anchor-start,
  anchor-end,
  label-pos,
) = {
  _path-elements(
    curve-api.from-cubic(segment),
    style,
    phase,
    anchor-start,
    anchor-end,
    label-pos,
  )
}

#let _pattern-edge-halves(halves, source-style, sink-style) = {
  let elements = ()
  let whole = halves.at("whole", default: none)
  if whole != none and _same-draw-path-style(source-style, sink-style) {
    return _derived-segments-elements(
      whole,
      source-style,
      auto,
      true,
      true,
    ).elements
  }
  if _same-pattern-geometry(source-style, sink-style) {
    let source = _derived-segments-elements(
      halves.source,
      source-style,
      auto,
      true,
      false,
    )
    let wavelength = _style-value(source-style, "pattern-wavelength")
    let phase = _style-value(source-style, "pattern-phase")
    let hidden-length = halves.at("split-gap", default: 0)
    let sink-phase = (
      phase + 2 * calc.pi * (source.length + hidden-length) / wavelength
    )
    let sink-elements = _derived-segments-elements(
      halves.sink,
      sink-style,
      sink-phase,
      false,
      true,
    ).elements
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
    let source-elements = _derived-segments-elements(
      halves.source,
      source-style,
      auto,
      true,
      true,
    ).elements
    let sink-elements = _derived-segments-elements(
      halves.sink,
      sink-style,
      auto,
      true,
      true,
    ).elements
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
  _segment-elements(
    curve-api.line-segment(_point(start), _point(end)),
    style,
    auto,
    true,
    true,
    label-pos,
  ).elements
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

#let _dangling-path(start, end, bend, style, dangling-at-start: false) = {
  let segment = _bent-line-segment(start, end, bend)
  let tangent = _style-value(style, "dangling-tangent")
  if tangent == auto {
    return curve-api.from-cubic(segment)
  }
  if tangent != "horizontal" and tangent != "vertical" {
    panic(
      "draw: dangling-tangent must be auto, \"horizontal\", or \"vertical\"",
    )
  }

  let distance = _point-distance(start, end)
  if distance == 0 {
    return curve-api.from-cubic(segment)
  }
  let control-distance = _style-value(style, "anchor-control-distance")
  if control-distance == auto or control-distance <= 0 {
    control-distance = distance / 3
  } else {
    control-distance = calc.min(control-distance, 0.9 * distance)
  }
  let delta = if tangent == "horizontal" {
    _point-x(end) - _point-x(start)
  } else {
    _point-y(end) - _point-y(start)
  }
  let direction = if delta < 0 { -1 } else { 1 }
  if tangent == "horizontal" and dangling-at-start {
    segment.control-start = (
      _point-x(start) + direction * control-distance,
      _point-y(start),
    )
  } else if tangent == "horizontal" {
    segment.control-end = (
      _point-x(end) - direction * control-distance,
      _point-y(end),
    )
  } else if dangling-at-start {
    segment.control-start = (
      _point-x(start),
      _point-y(start) + direction * control-distance,
    )
  } else {
    segment.control-end = (
      _point-x(end),
      _point-y(end) - direction * control-distance,
    )
  }
  curve-api.from-cubic(segment)
}

#let _pattern-dangling(path, style, label-pos) = {
  _path-elements(path, style, auto, true, true, label-pos).elements
}

#let _crossing-target(crossing-paths, under, eid) = {
  let targets = if type(under) == int {
    crossing-paths.ids
  } else if type(under) == label {
    crossing-paths.names
  } else {
    panic(
      "draw: crossing-under on edge "
        + str(eid)
        + " must be an integer edge id or Typst edge name",
    )
  }
  let key = str(under)
  if not targets.keys().contains(key) {
    let shown = if type(under) == label { "<" + key + ">" } else { key }
    panic(
      "draw: crossing-under on edge "
        + str(eid)
        + " refers to unknown or invisible edge "
        + shown,
    )
  }
  let target = targets.at(key)
  if target.eid == eid {
    panic(
      "draw: crossing-under on edge " + str(eid) + " cannot refer to itself",
    )
  }
  target
}

#let _crossing-arclengths(path, style, crossing-paths, eid) = {
  let under = _style-value(style, "crossing-under")
  if under == none {
    return ()
  }
  let target = _crossing-target(crossing-paths, under, eid)
  let accuracy = _style-value(style, "accuracy")
  curve-api
    .intersections(path, target.path, accuracy: accuracy)
    .map(hit => hit.at("distance-a"))
}

#let _cut-path-elements(path, style, cuts, mark-style: none) = {
  let paint-style = _without-mark-style(style)
  let gap = calc.max(0, _style-value(style, "crossing-gap"))
  if cuts.len() == 0 or gap == 0 {
    let elements = _segments-elements(
      curve-api.segments(path),
      paint-style,
      auto,
      true,
      true,
    ).elements
    if mark-style != none {
      elements += _mark-carrier-elements(path, mark-style)
    }
    return elements
  }
  let accuracy = _style-value(style, "accuracy")
  let total = curve-api.length(path, accuracy: accuracy)
  let cursor = 0
  let elements = ()
  for center in cuts {
    let cut-start = calc.max(cursor, calc.min(total, center - gap / 2))
    if cut-start > cursor {
      let piece = curve-api.trim(
        path,
        start-outset: cursor,
        end-outset: total - cut-start,
        accuracy: accuracy,
      )
      let phase = if _has-pattern(style) {
        (
          float(_style-value(style, "pattern-phase"))
            + 2 * calc.pi * cursor / _style-value(style, "pattern-wavelength")
        )
      } else {
        auto
      }
      elements += _segments-elements(
        curve-api.segments(piece),
        paint-style,
        phase,
        cursor == 0,
        false,
      ).elements
    }
    cursor = calc.max(cursor, calc.min(total, center + gap / 2))
  }
  if cursor < total {
    let piece = curve-api.trim(path, start-outset: cursor, accuracy: accuracy)
    let phase = if _has-pattern(style) {
      (
        float(_style-value(style, "pattern-phase"))
          + 2 * calc.pi * cursor / _style-value(style, "pattern-wavelength")
      )
    } else {
      auto
    }
    elements += _segments-elements(
      curve-api.segments(piece),
      paint-style,
      phase,
      false,
      true,
    ).elements
  }
  if mark-style != none {
    elements += _mark-carrier-elements(path, mark-style)
  }
  elements
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
    let size = measure(text(
      top-edge: "cap-height",
      bottom-edge: "bounds",
      label,
    ))
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

#let _subgraph-record(entry, default-style) = {
  if type(entry) == dictionary {
    let graph = entry.at("subgraph", default: entry.at("graph", default: none))
    if graph == none {
      panic("draw: subgraph dictionary entries need a `subgraph` field")
    }
    (
      hedges: subgraph-api.hedges(graph),
      edge-style: entry.at(
        "edge-style",
        default: entry.at(
          "style",
          default: entry.at("subgraph-edge-style", default: default-style),
        ),
      ),
    )
  } else {
    (
      hedges: subgraph-api.hedges(entry),
      edge-style: default-style,
    )
  }
}

#let _subgraph-records(subgraph, default-style) = {
  if subgraph == none {
    ()
  } else if type(subgraph) == array {
    subgraph
      .filter(entry => entry != none)
      .map(entry => _subgraph-record(entry, default-style))
  } else {
    (_subgraph-record(subgraph, default-style),)
  }
}

#let _subgraph-edge-styles(records, half-edge, edge-data) = {
  if half-edge == none {
    ()
  } else {
    let styles = ()
    for record in records {
      if record.hedges.contains(half-edge.hedge) {
        styles.push(_style(record.edge-style, edge-data))
      }
    }
    styles
  }
}

#let _last-style(styles) = {
  if styles.len() == 0 { none } else { styles.last() }
}

#let _node-pos(node) = {
  if node.pos == none {
    panic(
      "draw: graph node has no position; call layout(...) or pass pos: to graph.node/build",
    )
  }
  node.pos
}

#let _midpoint(a, b) = (
  (_point-x(a) + _point-x(b)) / 2,
  (_point-y(a) + _point-y(b)) / 2,
)

#let _edge-pos(edge, nodes) = {
  if edge.pos != none {
    edge.pos
  } else if edge.source != none and edge.sink != none {
    _midpoint(_node-pos(nodes.at(edge.source.node)), _node-pos(nodes.at(
      edge.sink.node,
    )))
  } else {
    panic(
      "draw: dangling graph edge has no position; call layout(...) or pass pos: to graph.edge/build",
    )
  }
}

#let _half-edge-label(half-edge, show-ids) = {
  let data = half-edge.at("data", default: (:))
  if type(data) == dictionary and data.keys().contains("label") {
    _as-content(data.label)
  } else if show-ids {
    [$h_(#half-edge.hedge)$]
  } else {
    none
  }
}

#let _half-edge-label-pos(segments, t, side, fallback) = {
  if segments == none or segments.len() == 0 {
    return fallback
  }
  let target = t * segments.len()
  let segment = segments.last()
  let segment-t = 1
  for (index, candidate) in segments.enumerate() {
    if target >= index and target < index + 1 {
      segment = candidate
      segment-t = target - index
    }
  }
  let point = curve-api.cubic-point(segment, segment-t)
  let tangent = curve-api.cubic-tangent(segment, segment-t)
  let length = _point-length(tangent)
  let normal = if length <= 1e-9 {
    (0, 0.12 * side)
  } else {
    _point-scale((-_point-y(tangent), _point-x(tangent)), 0.12 * side / length)
  }
  _point-add(point, normal)
}

#let _drawing-edge-record(e, nodes, node-boxes, scope) = {
  let edge = e + (pos: _edge-pos(e, nodes))
  let source-half-edge = edge.source
  let sink-half-edge = edge.sink
  let source-statement = if source-half-edge == none { none } else {
    source-half-edge.statement
  }
  let sink-statement = if sink-half-edge == none { none } else {
    sink-half-edge.statement
  }
  let source-node = if source-half-edge == none { none } else {
    node-boxes.at(source-half-edge.node).node
  }
  let sink-node = if sink-half-edge == none { none } else {
    node-boxes.at(sink-half-edge.node).node
  }
  let data = edge.statements
  let record-label = _data-label(edge)
  let statement-label = data.at("label", default: none)
  let data-label = if record-label == none { statement-label } else {
    record-label
  }
  let label-pos = if edge.label-pos == none { edge.pos } else { edge.label-pos }
  let aliases = (
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
    orientation: edge.orientation,
    ext: source-half-edge == none or sink-half-edge == none,
  )
  let fields = data + _data-fields(edge) + aliases
  (
    edge: edge,
    source-half-edge: source-half-edge,
    sink-half-edge: sink-half-edge,
    data-label: data-label,
    label-pos: label-pos,
    edge-data: (
      scope + fields + (fields: fields)
    ),
  )
}

#let _has-crossing-style(style) = (
  style != none and style.at("crossing-under", default: none) != none
)

#let _has-crossing-layer(layers) = layers.any(_has-crossing-style)

#let _has-visible-stroke(style) = if style == none {
  false
} else {
  let stroke = style.at("stroke", default: none)
  (
    stroke != none
      and (
        type(stroke) != dictionary or stroke.at("paint", default: black) != none
      )
  )
}

#let _primary-edge-path(
  record,
  nodes,
  node-boxes,
  node-outsets,
  geometry-style,
  edge-stroke,
  edge-omega,
  edge-trim-accuracy,
) = {
  let source-layer = _style-layer(record.source-style-layers, 0)
  let sink-layer = _style-layer(record.sink-style-layers, 0)
  let source-style = if source-layer == none {
    none
  } else {
    (stroke: edge-stroke) + geometry-style + source-layer
  }
  let sink-style = if sink-layer == none {
    none
  } else {
    (stroke: edge-stroke) + geometry-style + sink-layer
  }
  if (
    not _has-visible-stroke(source-style)
      and not _has-visible-stroke(sink-style)
  ) {
    return none
  }

  let edge = record.edge
  let source-half-edge = record.source-half-edge
  let sink-half-edge = record.sink-half-edge
  if source-half-edge != none and sink-half-edge != none {
    let source-geometry-style = if source-style == none { sink-style } else {
      source-style
    }
    let sink-geometry-style = if sink-style == none { source-style } else {
      sink-style
    }
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
        label-pos: record.label-pos,
      ),
    )
    let segments = if halves.whole == none {
      halves.source + halves.sink
    } else {
      halves.whole
    }
    if segments.len() == 0 { none } else { _segments-path(segments) }
  } else if source-half-edge != none {
    let source-geometry-style = if source-style == none { sink-style } else {
      source-style
    }
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
    let path = _dangling-path(
      line-start,
      edge.pos,
      edge.at("bend", default: none),
      source-geometry-style,
    )
    _path-layer(path, source-geometry-style, 0, 0, record.label-pos, auto)
  } else {
    let sink-geometry-style = if sink-style == none { source-style } else {
      sink-style
    }
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
    let path = _dangling-path(
      edge.pos,
      line-end,
      edge.at("bend", default: none),
      sink-geometry-style,
      dangling-at-start: true,
    )
    _path-layer(path, sink-geometry-style, 0, 0, record.label-pos, auto)
  }
}

#let _crossing-path-index(
  records,
  nodes,
  node-boxes,
  node-outsets,
  geometry-style,
  edge-stroke,
  edge-omega,
  edge-trim-accuracy,
) = {
  let ids = (:)
  let names = (:)
  for record in records {
    let path = _primary-edge-path(
      record,
      nodes,
      node-boxes,
      node-outsets,
      geometry-style,
      edge-stroke,
      edge-omega,
      edge-trim-accuracy,
    )
    if path != none {
      let entry = (eid: record.edge.edge, path: path)
      ids.insert(str(entry.eid), entry)
      let name = record.edge.at("name", default: none)
      if name != none {
        names.insert(str(name), entry)
      }
    }
  }
  (ids: ids, names: names)
}

#let draw(graph, options) = {
  let scope = options.scope
  let unit = options.unit
  let title = options.title
  let subgraph = options.subgraph
  let debug = options.debug
  let show-half-edge-ids = options.show-half-edge-ids
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
  let edge-style = options.edge-style
  let edge-offset = options.edge-offset
  let edge-length = options.edge-length
  let edge-ratio = options.edge-ratio
  let edge-resolve-length = options.edge-resolve-length
  let edge-accuracy = options.edge-accuracy
  let edge-optimize = options.edge-optimize
  let edge-split-gap = options.edge-split-gap
  let edge-dangling-tangent = options.edge-dangling-tangent
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
        let subgraph-records = _subgraph-records(subgraph, subgraph-edge-style)

        for (i, v) in nodes.enumerate() {
          let pos = _node-pos(v)
          let node = v + (pos: pos)
          let record-label = _data-label(v)
          let statement-label = v.statements.at("label", default: none)
          let data-label = if record-label == none { statement-label } else {
            record-label
          }
          let aliases = (
            vid: i,
            node: node,
            name: node.name,
            data: v.at("data", default: none),
            label: data-label,
          )
          let fields = v.statements + _data-fields(v) + aliases
          let node-data = scope + fields + (fields: fields)
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
                radius: _node-radius(
                  ctx,
                  label,
                  node-style-value,
                  node-min-radius,
                  node-label-padding,
                ),
              )
          )
          let node-label-draw-style = (
            _style(graph-node-label-style, node-data) + node-label-style
          )
          let layout-width = _statement-number(v, "layout-width", default: none)
          let layout-height = _statement-number(
            v,
            "layout-height",
            default: none,
          )
          let node-width = if layout-width == none {
            2 * node-style.radius
          } else { layout-width }
          let node-height = if layout-height == none {
            2 * node-style.radius
          } else { layout-height }
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
            node-elements.push(cetz.draw.circle(
              _point(pos),
              name: box.name,
              ..node-style,
            ))
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

        let geometry-style = (
          _edge-geometry-defaults
            + (
              offset: edge-offset,
              length: edge-length,
              ratio: edge-ratio,
              resolve-length: edge-resolve-length,
              accuracy: edge-accuracy,
              optimize: edge-optimize,
              split-gap: edge-split-gap,
              dangling-tangent: edge-dangling-tangent,
            )
        )
        let edge-records = ()
        for e in edges {
          let record = _drawing-edge-record(e, nodes, node-boxes, scope)
          let source-subgraph-styles = _subgraph-edge-styles(
            subgraph-records,
            record.source-half-edge,
            record.edge-data,
          )
          let sink-subgraph-styles = _subgraph-edge-styles(
            subgraph-records,
            record.sink-half-edge,
            record.edge-data,
          )
          let source-style-layers = _edge-style-layers(
            edge-style,
            source-style,
            record.edge-data,
          )
          let sink-style-layers = _edge-style-layers(
            edge-style,
            sink-style,
            record.edge-data,
          )
          let hidden = (
            source-style-layers.len() == 0 and sink-style-layers.len() == 0
          )
          let has-attached-label = (
            _has-layer-label(source-style-layers, record.edge-data)
              or _has-layer-label(sink-style-layers, record.edge-data)
          )
          edge-records.push(
            record
              + (
                source-subgraph-styles: source-subgraph-styles,
                sink-subgraph-styles: sink-subgraph-styles,
                source-in-subgraph: source-subgraph-styles.len() > 0,
                sink-in-subgraph: sink-subgraph-styles.len() > 0,
                source-style-layers: source-style-layers,
                sink-style-layers: sink-style-layers,
                layer-count: calc.max(
                  source-style-layers.len(),
                  sink-style-layers.len(),
                ),
                ev-label: if hidden or has-attached-label {
                  none
                } else {
                  _content(edge-label, record.edge-data, _as-content(
                    record.data-label,
                  ))
                },
                edge-label-draw-style: (
                  _style(graph-edge-label-style, record.edge-data)
                    + _style(edge-label-style, record.edge-data)
                ),
              ),
          )
        }
        let has-crossings = edge-records.any(record => (
          _has-crossing-layer(record.source-style-layers)
            or _has-crossing-layer(record.sink-style-layers)
            or (
              not subgraph-edge-underlay
                and (
                  _has-crossing-style(_last-style(
                    record.source-subgraph-styles,
                  ))
                    or _has-crossing-style(_last-style(
                      record.sink-subgraph-styles,
                    ))
                )
            )
        ))
        let crossing-paths = if has-crossings {
          _crossing-path-index(
            edge-records,
            nodes,
            node-boxes,
            node-outsets,
            geometry-style,
            edge-stroke,
            edge-omega,
            edge-trim-accuracy,
          )
        } else {
          (ids: (:), names: (:))
        }

        for (i, record) in edge-records.enumerate() {
          let edge = record.edge
          let source-half-edge = record.source-half-edge
          let sink-half-edge = record.sink-half-edge
          let label-pos = record.label-pos
          let edge-data = record.edge-data
          let source-subgraph-styles = record.source-subgraph-styles
          let sink-subgraph-styles = record.sink-subgraph-styles
          let source-in-subgraph = record.source-in-subgraph
          let sink-in-subgraph = record.sink-in-subgraph
          let source-style-layers = record.source-style-layers
          let sink-style-layers = record.sink-style-layers
          let layer-count = record.layer-count
          let ev-label = record.ev-label
          let edge-label-draw-style = record.edge-label-draw-style
          let source-label-segments = none
          let sink-label-segments = none
          let self-loop = (
            source-half-edge != none
              and sink-half-edge != none
              and source-half-edge.node == sink-half-edge.node
          )

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
              source-style-value + _last-style(source-subgraph-styles)
            } else {
              source-style-value
            }
            source-draw-style = _orient-mark-style(
              source-draw-style,
              edge-data,
              "source",
            )
            let sink-draw-style = if sink-style-value == none {
              none
            } else if sink-in-subgraph and not subgraph-edge-underlay {
              sink-style-value + _last-style(sink-subgraph-styles)
            } else {
              sink-style-value
            }
            sink-draw-style = _orient-mark-style(
              sink-draw-style,
              edge-data,
              "sink",
            )
            let source-crossing-under = if source-draw-style == none {
              none
            } else {
              _style-value(source-draw-style, "crossing-under")
            }
            let sink-crossing-under = if sink-draw-style == none {
              none
            } else {
              _style-value(sink-draw-style, "crossing-under")
            }
            let source-crossing-gap = if source-draw-style == none {
              none
            } else {
              _style-value(source-draw-style, "crossing-gap")
            }
            let sink-crossing-gap = if sink-draw-style == none {
              none
            } else {
              _style-value(sink-draw-style, "crossing-gap")
            }
            let crossing-under = if source-crossing-under != none {
              source-crossing-under
            } else {
              sink-crossing-under
            }
            if (
              crossing-under != none
                and subgraph-edge-underlay
                and (source-in-subgraph or sink-in-subgraph)
            ) {
              panic(
                "draw: crossing-under is incompatible with a subgraph edge underlay",
              )
            }

            if source-style-value == none and sink-style-value == none {
              ()
            } else if source-half-edge != none and sink-half-edge != none {
              let source-geometry-style = if source-style-value == none {
                sink-style-value
              } else { source-style-value }
              let sink-geometry-style = if sink-style-value == none {
                source-style-value
              } else { sink-style-value }
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
              if source-label-segments == none and halves.source.len() > 0 {
                source-label-segments = halves.source
              }
              if sink-label-segments == none and halves.sink.len() > 0 {
                sink-label-segments = halves.sink
              }
              if (
                source-style-value != none
                  and source-in-subgraph
                  and subgraph-edge-underlay
              ) {
                for subgraph-style in source-subgraph-styles {
                  for element in (
                    _segments-elements(
                      halves.source,
                      _without-mark-style(
                        _without-pattern-style(source-style-value),
                      )
                        + subgraph-style,
                      auto,
                      true,
                      true,
                    ).elements
                  ) {
                    elements.push(element)
                  }
                }
              }
              if (
                sink-style-value != none
                  and sink-in-subgraph
                  and subgraph-edge-underlay
              ) {
                for subgraph-style in sink-subgraph-styles {
                  for element in (
                    _segments-elements(
                      halves.sink,
                      _without-mark-style(
                        _without-pattern-style(sink-style-value),
                      )
                        + subgraph-style,
                      auto,
                      true,
                      true,
                    ).elements
                  ) {
                    elements.push(element)
                  }
                }
              }
              if crossing-under != none {
                let source-paint-style = if source-draw-style == none {
                  none
                } else {
                  _without-mark-style(source-draw-style)
                }
                let sink-paint-style = if sink-draw-style == none {
                  none
                } else {
                  _without-mark-style(sink-draw-style)
                }
                if (
                  source-draw-style == none
                    or sink-draw-style == none
                    or source-crossing-under != sink-crossing-under
                    or source-crossing-gap != sink-crossing-gap
                    or halves.whole == none
                    or not _same-draw-path-style(
                      source-paint-style,
                      sink-paint-style,
                    )
                ) {
                  panic(
                    "draw: crossing-under on paired edge "
                      + str(edge-data.eid)
                      + " requires one continuous source/sink style layer",
                  )
                }
                let path = _segments-path(halves.whole)
                let mark-style = if _has-mark(source-draw-style) {
                  source-draw-style
                } else if _has-mark(sink-draw-style) {
                  sink-draw-style
                } else {
                  none
                }
                mark-style = _paired-carrier-style(path, halves, mark-style)
                let cuts = _crossing-arclengths(
                  path,
                  source-paint-style,
                  crossing-paths,
                  edge-data.eid,
                )
                for element in _cut-path-elements(
                  path,
                  source-paint-style,
                  cuts,
                  mark-style: mark-style,
                ) {
                  elements.push(element)
                }
              } else if (
                source-style-value != none and sink-style-value != none
              ) {
                let source-paint-style = _without-mark-style(source-draw-style)
                let sink-paint-style = _without-mark-style(sink-draw-style)
                let center-mark-style = _center-mark-style(
                  source-draw-style,
                  sink-draw-style,
                )
                if (
                  center-mark-style != none
                    and halves.whole != none
                    and _same-draw-path-style(
                      source-paint-style,
                      sink-paint-style,
                    )
                ) {
                  let path = _segments-path(halves.whole)
                  for element in (
                    _segments-elements(
                      halves.whole,
                      source-paint-style,
                      auto,
                      true,
                      true,
                    ).elements
                  ) {
                    elements.push(element)
                  }
                  let carrier-style = _paired-carrier-style(
                    path,
                    halves,
                    center-mark-style,
                  )
                  for element in _mark-carrier-elements(path, carrier-style) {
                    elements.push(element)
                  }
                } else {
                  for element in _pattern-edge-halves(
                    halves,
                    source-draw-style,
                    sink-draw-style,
                  ) {
                    elements.push(element)
                  }
                }
              } else if source-style-value != none {
                for element in (
                  _derived-segments-elements(
                    halves.source,
                    source-draw-style,
                    auto,
                    true,
                    true,
                  ).elements
                ) {
                  elements.push(element)
                }
              } else if sink-style-value != none {
                for element in (
                  _derived-segments-elements(
                    halves.sink,
                    sink-draw-style,
                    auto,
                    true,
                    true,
                  ).elements
                ) {
                  elements.push(element)
                }
              }
              let attached-label = _paired-layer-label-element(
                ctx,
                halves,
                source-draw-style,
                sink-draw-style,
                label-pos,
                edge-data,
              )
              if attached-label != none {
                elements.push(attached-label)
              }
            } else if source-half-edge != none {
              let source-geometry-style = if source-style-value == none {
                sink-style-value
              } else { source-style-value }
              let source-anchor = _style-value(
                source-geometry-style,
                "source-anchor",
              )
              let line-start = if source-anchor == auto {
                curve-api.outset-point(
                  _point(_node-pos(nodes.at(source-half-edge.node))),
                  _point(edge.pos),
                  distance: node-outsets.at(source-half-edge.node),
                )
              } else {
                _node-anchor-point(
                  node-boxes.at(source-half-edge.node),
                  source-anchor,
                )
              }
              let bend = edge.at("bend", default: none)
              let dangling-path = _dangling-path(
                line-start,
                edge.pos,
                bend,
                source-geometry-style,
              )
              let visible-path = _path-layer(
                dangling-path,
                source-geometry-style,
                0,
                0,
                label-pos,
                auto,
              )
              if (
                source-label-segments == none and source-geometry-style != none
              ) {
                source-label-segments = _geometry-path-segments(
                  dangling-path,
                  source-geometry-style,
                  0,
                  0,
                  label-pos,
                  auto,
                )
              }
              if (
                source-style-value != none
                  and source-in-subgraph
                  and subgraph-edge-underlay
              ) {
                for subgraph-style in source-subgraph-styles {
                  for element in _pattern-dangling(
                    dangling-path,
                    _without-mark-style(_without-pattern-style(
                      _dangling-mark-style(source-style-value),
                    ))
                      + subgraph-style,
                    label-pos,
                  ) {
                    elements.push(element)
                  }
                }
              }
              if source-style-value != none {
                let draw-style = _dangling-mark-style(source-draw-style)
                let cuts = _crossing-arclengths(
                  visible-path,
                  draw-style,
                  crossing-paths,
                  edge-data.eid,
                )
                let drawn = if (
                  _style-value(draw-style, "crossing-under") == none
                ) {
                  _pattern-dangling(dangling-path, draw-style, label-pos)
                } else {
                  _cut-path-elements(
                    visible-path,
                    _without-mark-style(draw-style),
                    cuts,
                    mark-style: if _has-mark(draw-style) { draw-style } else {
                      none
                    },
                  )
                }
                for element in drawn {
                  elements.push(element)
                }
                let attached-label = _layer-label-element(
                  ctx,
                  visible-path,
                  draw-style,
                  label-pos,
                  edge-data,
                )
                if attached-label != none {
                  elements.push(attached-label)
                }
              }
            } else if sink-half-edge != none {
              let sink-geometry-style = if sink-style-value == none {
                source-style-value
              } else { sink-style-value }
              let sink-anchor = _style-value(sink-geometry-style, "sink-anchor")
              let line-end = if sink-anchor == auto {
                curve-api.outset-point(
                  _point(_node-pos(nodes.at(sink-half-edge.node))),
                  _point(edge.pos),
                  distance: node-outsets.at(sink-half-edge.node),
                )
              } else {
                _node-anchor-point(
                  node-boxes.at(sink-half-edge.node),
                  sink-anchor,
                )
              }
              let bend = edge.at("bend", default: none)
              let dangling-path = _dangling-path(
                edge.pos,
                line-end,
                bend,
                sink-geometry-style,
                dangling-at-start: true,
              )
              let visible-path = _path-layer(
                dangling-path,
                sink-geometry-style,
                0,
                0,
                label-pos,
                auto,
              )
              if sink-label-segments == none and sink-geometry-style != none {
                sink-label-segments = _geometry-path-segments(
                  dangling-path,
                  sink-geometry-style,
                  0,
                  0,
                  label-pos,
                  auto,
                )
              }
              if (
                sink-style-value != none
                  and sink-in-subgraph
                  and subgraph-edge-underlay
              ) {
                for subgraph-style in sink-subgraph-styles {
                  for element in _pattern-dangling(
                    dangling-path,
                    _without-mark-style(_without-pattern-style(
                      _dangling-mark-style(sink-style-value),
                    ))
                      + subgraph-style,
                    label-pos,
                  ) {
                    elements.push(element)
                  }
                }
              }
              if sink-style-value != none {
                let draw-style = _dangling-mark-style(sink-draw-style)
                let cuts = _crossing-arclengths(
                  visible-path,
                  draw-style,
                  crossing-paths,
                  edge-data.eid,
                )
                let drawn = if (
                  _style-value(draw-style, "crossing-under") == none
                ) {
                  _pattern-dangling(dangling-path, draw-style, label-pos)
                } else {
                  _cut-path-elements(
                    visible-path,
                    _without-mark-style(draw-style),
                    cuts,
                    mark-style: if _has-mark(draw-style) { draw-style } else {
                      none
                    },
                  )
                }
                for element in drawn {
                  elements.push(element)
                }
                let attached-label = _layer-label-element(
                  ctx,
                  visible-path,
                  draw-style,
                  label-pos,
                  edge-data,
                )
                if attached-label != none {
                  elements.push(attached-label)
                }
              }
            }
          }

          if ev-label != none {
            elements.push(cetz.draw.content(
              _point(label-pos),
              ev-label,
              padding: 0,
              ..edge-label-draw-style,
            ))
          }

          for half-edge in (source-half-edge, sink-half-edge) {
            if half-edge != none {
              let label = _half-edge-label(half-edge, show-half-edge-ids)
              if label != none {
                let node-pos = _node-pos(nodes.at(half-edge.node))
                let fallback-direction = _point-sub(edge.pos, node-pos)
                let fallback-length = _point-length(fallback-direction)
                let fallback-normal = if fallback-length <= 1e-9 {
                  (0, 0.12)
                } else {
                  _point-scale(
                    (
                      -_point-y(fallback-direction),
                      _point-x(fallback-direction),
                    ),
                    0.12 / fallback-length,
                  )
                }
                let fallback = _point-add(
                  _point-lerp(node-pos, edge.pos, 0.38),
                  fallback-normal,
                )
                let source = (
                  source-half-edge != none
                    and half-edge.hedge == source-half-edge.hedge
                )
                let segments = if source {
                  source-label-segments
                } else {
                  sink-label-segments
                }
                let t = if source or self-loop { 0.38 } else { 0.62 }
                let side = if source { 1 } else { -1 }
                elements.push(cetz.draw.content(
                  _half-edge-label-pos(segments, t, side, fallback),
                  text(size: 0.65em)[#label],
                  padding: 0,
                ))
              }
            }
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
      if options.draw-after != none {
        _call(options.draw-after, graph)
      }
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
