// Internal drawing implementation. Public users should import `draw.typ`.

#import "@preview/cetz:0.5.1" as cetz
#import "../curve.typ" as curve-api
#import "../graph.typ" as graph-api
#import "../subgraph.typ" as subgraph-api

#let _plugin = plugin("../../linnest.wasm")

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

// Preserve the caller's span indexes while removing solver-collapsed knots
// from the Hobby solve. Empty spans keep the half-edge split at its anchor.
#let _routed-split-through(points, omega: 1.0, accuracy: 0.001) = {
  let knots = ()
  let indexes = ()
  for point in points.map(_point) {
    if knots.len() == 0 or _point-distance(knots.last(), point) > 2.220446049250313e-16 {
      knots.push(point)
    }
    indexes.push(knots.len() - 1)
  }
  let split = if knots.len() < 2 {
    (curve: curve-api.path(), parts: ())
  } else {
    curve-api.split-through(knots, omega: omega, accuracy: accuracy)
  }
  (
    curve: split.curve,
    parts: range(calc.max(0, points.len() - 1)).map(i => {
      if indexes.at(i) == indexes.at(i + 1) { curve-api.path() }
      else { split.parts.at(indexes.at(i)) }
    }),
  )
}

#let _trim-routed-path(path, start-outset: 0, end-outset: 0, accuracy: 0.001) = {
  if curve-api.elements(path).len() == 0 or (start-outset == 0 and end-outset == 0) {
    return path
  }
  if start-outset + end-outset >= curve-api.length(path, accuracy: accuracy) {
    return curve-api.path()
  }
  curve-api.trim(path, start-outset: start-outset, end-outset: end-outset, accuracy: accuracy)
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

// Typst preserves links as transparent SVG rectangles. Floating content keeps
// these identity targets out of the graph's bounds and visible drawing styles.
#let _identity-href(kind, id, details) = {
  // SVG link attributes escape quotes but need JSON escapes for XML metacharacters.
  let details = json.encode(details).replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
  "#linnet-" + kind + "-" + str(id) + "?" + details
}

#let _identity-target(pos, href, width: 8pt, height: 8pt) = {
  ((position: _point(pos), body: link(href, box(width: width, height: height))),)
}

#let _draw-identity-targets(ctx, targets) = {
  if "content-many" not in cetz.draw {
    let result = cetz.process.many(ctx, targets.map(target => {
      cetz.draw.floating(cetz.draw.content(target.position, target.body, padding: 0))
    }).flatten(), compute-bounds: false)
    return (ctx: result.ctx, drawables: result.drawables)
  }
  cetz.draw.content-many(ctx, targets, padding: 0, tags: (cetz.drawable.TAG.no-bounds,))
}

#let _edge-identity-targets(ctx, parts, hrefs) = {
  let targets = ()
  // Preserve arc-length quarters and the control-polygon sampling bound while
  // batching path construction, trimming and de Casteljau evaluation.
  for (region, points) in curve-api.region-samples(
    parts, regions: 4, unit: ctx.length / 1pt, step: 4,
  ) {
    let body = link(hrefs.at(region), box(width: 8pt, height: 8pt))
    targets += points.map(point => (position: point, body: body))
  }
  targets
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
  label-only: false,
  label: none,
  label-style: (:),
  label-shift: 0,
  label-slide: true,
  label-path: none,
  label-gap: 0.15,
  label-side: auto,
)

#let _edge-routing-defaults = (
  source-anchor: auto,
  sink-anchor: auto,
  anchor-control-distance: auto,
  route: "edge-pos",
  route-points: "through",
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
  pattern-fit: false,
  pattern-phase: 0,
  pattern-samples-per-period: 16,
  pattern-coil-longitudinal-scale: 1.25,
  pattern-endpoint-slope: 0,
  pattern-natural-endpoints: false,
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

#let _half-edge-style(record, side) = {
  let half-edge = record.at(side + "-half-edge", default: none)
  if half-edge == none { return (:) }
  let data = half-edge.at("data", default: none)
  if type(data) != dictionary { return (:) }
  let routing = (:)
  if data.keys().contains("anchor") {
    routing.insert(side + "-anchor", data.anchor)
  }
  if data.keys().contains("routing") {
    routing.insert("route", data.routing)
  }
  _overlay-style(routing, _call(
    data.at(side + "-style", default: data.at("style", default: (:))), record,
  ))
}

#let _edge-style-layers(defaults, options, data) = {
  // Resolve logical-edge callbacks once, before composing the two endpoints.
  let base = _call(defaults.at("edge-style", default: (:)), data)
  let configured = _call(options.edge-style, data)
  let local = _call(data.at("edge-style", default: auto), data)
  if base == none or configured == none or local == none {
    return (source: (), sink: ())
  }
  if base == auto { base = (:) }
  let routing = if data.keys().contains("routing") {
    (route: data.routing)
  } else { (:) }
  let decoration = if data.keys().contains("decoration") {
    let value = _call(data.decoration, data)
    if value == none { (pattern: none) }
    else if type(value) in (dictionary, array) { value }
    else { (pattern: value) }
  } else { auto }
  let result = (:)
  for side in ("source", "sink") {
    let style = _overlay-style(base, _call(
      defaults.at(side + "-style", default: (:)), data,
    ))
    style = _overlay-style(style, configured)
    style = _overlay-style(style, routing)
    style = _overlay-style(style, local)
    style = _overlay-style(style, _call(options.at(side + "-style"), data))
    style = _overlay-style(style, _call(data.at(side + "-style", default: (:)), data))
    style = _overlay-style(style, _half-edge-style(data, side))
    let layers = if type(style) == array { style } else { (style,) }
    // Per-edge decoration patches the base layer. Multiple explicit decorations
    // create multiple decorated layers while retaining any later custom layers.
    if decoration != auto {
      let patches = if type(decoration) == array { decoration } else { (decoration,) }
      let first = layers.at(0, default: (:))
      layers = (patches.map(patch => _overlay-style(first, patch))
        + layers.slice(calc.min(1, layers.len())))
    }
    result.insert(side, layers)
  }
  result
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
  (not _style-value(style, "label-only")
    and pattern != none and pattern != "normal" and pattern != "curve")
}

#let _has-mark(style) = (
  not _style-value(style, "label-only") and style.at("mark", default: none) != none
)

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
      and _style-value(source-style, "pattern-fit")
        == _style-value(sink-style, "pattern-fit")
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
      and _style-value(source-style, "pattern-endpoint-slope")
        == _style-value(sink-style, "pattern-endpoint-slope")
  )
  same = (
    same
      and _style-value(source-style, "pattern-natural-endpoints")
        == _style-value(sink-style, "pattern-natural-endpoints")
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
  (same-path-geometry
    and _style-value(source-style, "label-only") == _style-value(sink-style, "label-only")
    and _draw-style(source-style) == _draw-style(sink-style))
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
  if curve-api.elements(path).len() == 0 { return path }
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

#let _anchor-control-distance(style, start, route, end) = {
  let value = _style-value(style, "anchor-control-distance")
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
  source-amount,
  sink-amount,
) = {
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
  let route-direction = _point-sub(end, start)
  let route-direction-length = _point-length(route-direction)
  // Share a middle handle to keep the two halves tangent-continuous.
  let route-handle = if route-direction-length == 0 {
    (0, 0)
  } else {
    let max-handle = (
      calc.min(_point-distance(start, route), _point-distance(route, end)) / 3
    )
    _point-scale(
      route-direction,
      calc.min(source-amount, sink-amount, max-handle) / route-direction-length,
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
  pattern-curve: none,
) = {
  let shared = _same-layer-geometry(source-style, sink-style)
  let whole-path = if shared {
    _path-layer(
      if pattern-curve == none { curve } else { pattern-curve },
      source-style,
      source-outset,
      sink-outset,
      label-pos,
      auto,
    )
  } else { none }
  if (
    shared
      and _style-value(source-style, "offset-side") == "label"
      and curve-api.elements(whole-path).len() > 0
  ) {
    // The label chooses one side of the full carrier, including both paint halves.
    let resolved = (offset: whole-path.offset, offset-side: none)
    source-style += resolved
    sink-style += resolved
  }
  let geometries = ()
  let split-gap = 0
  for (index, style) in (source-style, sink-style).enumerate() {
    // Window distances belong to the offset geometry, including split layers.
    let offset-style = style + (length: none, ratio: none, shift: 0)
    let paths = (source-path, sink-path).map(path => _path-layer(
      path,
      offset-style,
      0,
      0,
      label-pos,
      auto,
    ))
    let lengths = ()
    for path in paths {
      lengths.push(curve-api.length(path, accuracy: accuracy))
    }
    let half-length = lengths.at(index)
    let split-outset = calc.min(
      half-length,
      calc.max(0, _style-value(style, "split-gap")) / 2,
    )
    split-gap += split-outset
    let outsets = _merged-half-outsets(
      _visible-half-outsets(
        lengths.sum(),
        if index == 0 { 0 } else { lengths.first() },
        half-length,
        style,
        source-outset,
        sink-outset,
      ),
      half-length,
      start: if index == 1 { split-outset } else { 0 },
      end: if index == 0 { split-outset } else { 0 },
    )
    geometries.push(_half-path-geometry(
      paths.at(index),
      style + (offset: 0),
      outsets,
      label-pos,
    ))
  }
  let whole-geometry = if shared { curve-api.segments(whole-path) } else {
    none
  }
  (
    source: geometries.first(),
    sink: geometries.last(),
    curve: curve,
    whole: if split-gap == 0 { whole-geometry } else { none },
    pattern-whole: whole-geometry,
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
  let source-amount = _anchor-control-distance(source-style, start, route, end)
  let sink-amount = _anchor-control-distance(sink-style, start, route, end)
  if (
    route-mode != "direct"
      and route-points-mode == "through"
      and (source-route.len() > 0 or sink-route.len() > 0)
  ) {
    let source-amount = _route-aware-anchor-amount(
      start,
      source-amount,
      source-route,
      route,
    )
    let sink-amount = _route-aware-anchor-amount(
      end,
      sink-amount,
      sink-route,
      route,
    )
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
    let split = _routed-split-through(
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
    let split = _anchored-cubic-route-split(
      start,
      source-anchor,
      route,
      sink-anchor,
      end,
      source-amount,
      sink-amount,
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
  let source-path = curve-api.hobby-spline(
    _anchor-points(start, source-anchor, route, source-amount),
    omega: omega,
    accuracy: accuracy,
  )
  let sink-path = curve-api.hobby-spline(
    _anchor-points(end, sink-anchor, route, sink-amount, reverse: true),
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
    pattern-curve: curve-api.path(source-path, sink-path),
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

  let carrier = _data-fields(edge).at("layout-carrier", default: none)
  let source-route = _route-points(edge.source)
  let sink-route = _route-points(edge.sink)
  let points = (
    _point(nodes.at(edge.source.node).pos),
    ..source-route,
    _point(edge.pos),
    ..sink-route.rev(),
    _point(nodes.at(edge.sink.node).pos),
  )
  let routed = carrier != none or source-route.len() > 0 or sink-route.len() > 0
  let split = if carrier != none {
    _routed-split-through(carrier, omega: omega, accuracy: accuracy)
  } else if routed {
    _routed-split-through(points, omega: omega, accuracy: accuracy)
  } else {
    curve-api.split-through(points,
      omega: omega, start-outset: source-outset,
      end-outset: sink-outset, accuracy: accuracy,
    )
  }
  let split-outset = calc.max(0, options.at("split-gap", default: 0)) / 2
  let boundary = source-route.len() + 1
  let half-length = if carrier == none { 0 } else {
    curve-api.length(split.curve, accuracy: accuracy) / 2
  }
  let source = if carrier != none {
    curve-api.trim(split.curve, end-outset: half-length, accuracy: accuracy)
  } else if source-route.len() == 0 { split.parts.at(0) }
    else { curve-api.path(..split.parts.slice(0, boundary)) }
  let sink = if carrier != none {
    curve-api.trim(split.curve, start-outset: half-length, accuracy: accuracy)
  } else if sink-route.len() == 0 { split.parts.at(boundary) }
    else { curve-api.path(..split.parts.slice(boundary)) }
  if routed {
    source = _trim-routed-path(source, start-outset: source-outset, accuracy: accuracy)
    sink = _trim-routed-path(sink, end-outset: sink-outset, accuracy: accuracy)
  }
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
  let routed-edge = if (
    _style-value(source-style, "route-points") == "ignore"
      and _style-value(sink-style, "route-points") == "ignore"
  ) {
    edge + (
      data: _data-fields(edge) + (layout-carrier: none),
      source: edge.source + (route-points: ()),
      sink: edge.sink + (route-points: ()),
    )
  } else { edge }
  let base = edge-halves(routed-edge, nodes, (
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


#let _bezier-element(segment, style) = {
  if _has-mark(style) {
    return _mark-carrier-elements(
      curve-api.from-cubic(segment),
      style,
      paint: true,
    ).flatten()
  }
  cetz.draw.bezier(
    _point(segment.start),
    _point(segment.end),
    _point(segment.control-start),
    _point(segment.control-end),
    .._draw-style(style),
  )
}

#let _segments-path(segments) = {
  curve-api.path(..segments.map(curve-api.from-cubic))
}

#let _pattern-path(
  path,
  style,
  phase: auto,
  anchor-start: true,
  anchor-end: true,
  split-at: (),
) = {
  let pattern-style = _pattern-style(style)
  let length = curve-api.length(path, accuracy: pattern-style.pattern-accuracy)
  let pattern = pattern-style.pattern
  let wavelength = pattern-style.pattern-wavelength
  let samples = pattern-style.pattern-samples-per-period
  let phase = if phase == auto { pattern-style.pattern-phase } else { phase }
  assert(wavelength > 0, message: "pattern-wavelength must be positive")
  // Fit the physical endpoints before dividing the decorated path for painting.
  let natural = (
    pattern-style.pattern-natural-endpoints and anchor-start and anchor-end
  )
  if natural and length > 0 {
    assert(
      type(pattern) == str and pattern.trim() in ("coil", "helix", "spring"),
      message: "pattern-natural-endpoints requires a built-in coil",
    )
    pattern = curve-api.coil(
      fit-length: length,
      amplitude: pattern-style.pattern-amplitude,
      wavelength: wavelength,
      samples-per-period: samples,
      longitudinal-scale: pattern-style.pattern-coil-longitudinal-scale,
    )
    wavelength = length
    samples = pattern.points.len() - 1
    phase = 0
  } else if (
    pattern-style.pattern-fit and anchor-start and anchor-end and length > 0
  ) {
    wavelength = length / calc.max(1, calc.round(length / wavelength))
  }
  let patterned = curve-api.pattern(
    path,
    pattern: pattern,
    amplitude: pattern-style.pattern-amplitude,
    wavelength: wavelength,
    phase: phase,
    samples-per-period: samples,
    coil-longitudinal-scale: pattern-style.pattern-coil-longitudinal-scale,
    anchor-start: anchor-start,
    anchor-end: anchor-end,
    endpoint-slope: pattern-style.pattern-endpoint-slope,
    accuracy: pattern-style.pattern-accuracy,
    split-at: split-at,
  )
  patterned
}

#let _segments-elements(segments, style, phase, anchor-start, anchor-end) = {
  if _style-value(style, "label-only") { return (elements: (), length: 0) }
  if _has-mark(style) {
    return _derived-path-elements(
      _segments-path(segments),
      style,
      phase,
      anchor-start,
      anchor-end,
    )
  }
  let elements = ()
  let length = 0
  if _has-pattern(style) {
    let pattern-style = _pattern-style(style)
    let path = _segments-path(segments)
    length = curve-api.length(path, accuracy: pattern-style.pattern-accuracy)
    let pattern = pattern-style.pattern
    let wavelength = pattern-style.pattern-wavelength
    let samples = pattern-style.pattern-samples-per-period
    let phase = if phase == auto { pattern-style.pattern-phase } else { phase }
    assert(wavelength > 0, message: "pattern-wavelength must be positive")
    // Fit only complete edges; split styles and crossing gaps retain phase continuity.
    let natural = pattern-style.pattern-natural-endpoints and anchor-start and anchor-end
    if natural and length > 0 {
      assert(
        type(pattern) == str and pattern.trim() in ("coil", "helix", "spring"),
        message: "pattern-natural-endpoints requires a built-in coil",
      )
      pattern = curve-api.coil(
        fit-length: length,
        amplitude: pattern-style.pattern-amplitude,
        wavelength: wavelength,
        samples-per-period: samples,
        longitudinal-scale: pattern-style.pattern-coil-longitudinal-scale,
      )
      wavelength = length
      samples = pattern.points.len() - 1
      phase = 0
    } else if pattern-style.pattern-fit and anchor-start and anchor-end and length > 0 {
      wavelength = length / calc.max(1, calc.round(length / wavelength))
    }
    let draw-style = _draw-style(style)
    let unit = draw-style.remove("unit", default: 1)
    let patterned = curve-api.pattern-to-cetz(
      path,
      pattern: pattern,
      amplitude: pattern-style.pattern-amplitude,
      wavelength: wavelength,
      phase: phase,
      samples-per-period: samples,
      coil-longitudinal-scale: pattern-style.pattern-coil-longitudinal-scale,
      anchor-start: anchor-start,
      anchor-end: anchor-end,
      endpoint-slope: pattern-style.pattern-endpoint-slope,
      accuracy: pattern-style.pattern-accuracy,
      unit: unit, style: draw-style,
    )
    elements.push(patterned)
  } else {
    for segment in segments {
      elements.push(curve-api.to-cetz(
        curve-api.from-cubic(segment),
        .._draw-style(style),
      ))
    }
  }
  (elements: elements, length: length)
}

// Windows are distances on the carrier, not on the longer decorated curve.
#let _path-window-elements(path, style, windows) = {
  let windows = windows.filter(window => window.end > window.start)
  if windows.len() == 0 { return () }
  let parts = if _has-pattern(style) {
    let boundaries = windows.map(window => (window.start, window.end)).flatten()
    _pattern-path(path, style, split-at: boundaries).parts
  } else { () }
  let total = curve-api.length(path, accuracy: _style-value(style, "accuracy"))
  let elements = ()
  for (index, window) in windows.enumerate() {
    if not _style-value(window.style, "label-only") {
      let piece = if _has-pattern(style) { parts.at(2 * index + 1) } else {
        curve-api.trim(
          path,
          start-outset: window.start,
          end-outset: total - window.end,
          accuracy: _style-value(style, "accuracy"),
        )
      }
      elements.push(curve-api.to-cetz(piece, .._draw-style(_without-mark-style(
        window.style,
      ))))
    }
  }
  elements
}

#let _half-paint-windows(halves, source-style, sink-style, total) = {
  if halves.at("split-gap", default: 0) == 0 and _same-draw-path-style(
    _without-mark-style(source-style), _without-mark-style(sink-style),
  ) {
    return ((start: 0, end: total, style: source-style),)
  }
  let accuracy = _style-value(source-style, "accuracy")
  let source-end = calc.min(total, curve-api.length(
    _segments-path(halves.source),
    accuracy: accuracy,
  ))
  let sink-start = if halves.at("split-gap", default: 0) == 0 {
    source-end
  } else {
    calc.max(
      source-end,
      total - curve-api.length(_segments-path(halves.sink), accuracy: accuracy),
    )
  }
  (
    (start: 0, end: source-end, style: source-style),
    (start: sink-start, end: total, style: sink-style),
  )
}

#let _mark-carrier-elements(path, style, paint: false) = {
  if style == none or not _has-mark(style) { return () }
  let segments = curve-api.segments(path)
  if segments.len() == 0 { return () }
  let style = _positioned-mark-style(style)
  let plain-style = _draw-style(style) + (mark: none)
  let element = if segments.len() == 1 {
    _bezier-element(segments.first(), plain-style)
  } else { curve-api.to-cetz(path, ..plain-style) }
  (
    (
      ctx => {
        let resolved = cetz.styles.resolve(
          ctx.style,
          merge: _draw-style(style),
          root: "bezier",
        )
        if ctx.at("_perspective-projection", default: false) {
          resolved.mark.transform-shape = true
        }
        let transform = ctx.transform
        let plain = element.first()(ctx + (transform: cetz.matrix.ident(4)))
        let carrier = plain.drawables.first()
        // Match CeTZ's coordinate space so physical mark sizes survive canvas transforms.
        let flat = not resolved.mark.transform-shape
        if flat {
          carrier = cetz
            .drawable
            .apply-transform(
              cetz.matrix.mul-mat(
                cetz.matrix.transform-scale((1, 1, 0)),
                transform,
              ),
              carrier,
            )
            .first()
        }
        let total = cetz.path-util.length(carrier.segments)
        let ratio = _mark-ratio(style)
        let mark-ctx = ctx + (resolve-coordinate: ())
        let marks = ()
        let cuts = (none, none)
        if total > 0 {
          for (side, root) in ("start", "end").enumerate() {
            let at = if ratio == none { 0 } else {
              if side == 0 {
                ratio * total + _style-value(style, "mark-shift")
              } else {
                (1 - ratio) * total - _style-value(style, "mark-shift")
              }
            }
            for (index, entry) in cetz
              .mark
              .process-style(mark-ctx, resolved.mark, root, total)
              .enumerate() {
              if entry.symbol == none { continue }
              if entry.pos != none { at = entry.pos }
              at += entry.offset
              let (shape, defaults) = cetz.mark-shapes.get-mark(
                mark-ctx,
                entry.symbol,
              )
              for flag in ("reverse", "flip", "harpoon") {
                entry.at(flag) = (
                  entry.at(flag) != defaults.at(flag, default: false)
                )
              }
              entry.mark = none
              let mark = cetz.mark._eval-mark-shape-and-anchors(
                mark-ctx,
                shape(entry),
                entry,
              )
              let symbol = ctx.marks.mnemonics.at(
                entry.symbol,
                default: entry.symbol,
              )
              let builtin = symbol not in ctx.marks.marks
              if builtin {
                symbol = cetz
                  .mark-shapes
                  .mnemonics
                  .at(symbol, default: (symbol, (:)))
                  .first()
              }
              // Triangle/straight stroke anchors are not their geometric tip/back;
              // in particular, straight's declared base is beside its tip.
              if builtin and symbol in ("triangle", "straight") {
                mark.tip = (0, 0, 0)
                mark.base = (entry.length, 0, 0)
                mark.center = cetz.vector.lerp(mark.tip, mark.base, 0.5)
                mark.reverse-tip = mark.base
                mark.reverse-base = mark.tip
                mark.reverse-center = mark.center
              }
              // Let CeTZ resolve reversal, slant and anchors on a straight reference,
              // then map both contacts to a chord of the full, unpatterned carrier.
              let reference = cetz.drawable.line-strip((mark.tip, mark.base))
              // Keep the contact pair even for a zero-length mark.
              reference.segments = ((mark.tip, false, (("l", mark.base),)),)
              mark.drawables.push(reference)
              mark = cetz.mark.transform-mark(
                entry,
                mark,
                (0, 0, 0),
                (1, 0, 0),
                reverse: entry.reverse,
                slant: entry.slant,
                flip: entry.flip,
                harpoon: entry.harpoon,
              )
              let reference = mark.drawables.pop().segments.first()
              let tip = reference.first()
              let back = reference.last().last().last()
              let axis = cetz.vector.sub(back, tip)
              let length = cetz.vector.len(axis)
              let compression = if length == 0 { 1 } else {
                calc.min(1, total / length)
              }
              let direction = if axis.first() < 0 { -1 } else { 1 }
              let stations = (tip, back).map(p => (
                at
                  + if length == 0 { 0 } else {
                    cetz.vector.dot(p, axis) / length * direction * compression
                  }
              ))
              let low = calc.min(..stations)
              let span = calc.abs(stations.last() - stations.first())
              let inward = calc.clamp(low, 0, calc.max(0, total - span)) - low
              stations = stations.map(s => calc.clamp(s + inward, 0, total))
              // Repeated controls need no special direction: only the two sampled points matter.
              let contacts = stations.map(s => cetz.path-util.point-at(
                carrier.segments,
                if side == 1 { total - s } else { s },
              ))
              let chord = cetz.vector.sub(
                contacts.last().point,
                contacts.first().point,
              )
              let angle = if cetz.vector.len(chord) == 0 { 0deg } else {
                calc.atan2(..chord.slice(0, 2))
              }
              let reference-angle = if length == 0 { 0deg } else {
                calc.atan2(..axis.slice(0, 2))
              }
              let alignment = cetz.matrix.mul-mat(
                cetz.matrix.transform-translate(..contacts.first().point),
                cetz.matrix.transform-rotate-z(angle),
                cetz.matrix.transform-scale((
                  if length == 0 { 1 } else { cetz.vector.len(chord) / length },
                  1,
                  1,
                )),
                cetz.matrix.transform-rotate-z(-reference-angle),
                cetz.matrix.transform-translate(..cetz.vector.scale(tip, -1)),
              )
              marks += cetz.drawable.apply-transform(alignment, mark.drawables)
              // An interior mark overlays the propagator; shortening would leave a gap.
              if (
                ratio == none
                  and entry.shorten-to != none
                  and (entry.shorten-to == auto or index <= entry.shorten-to)
              ) {
                let contact = if builtin and symbol == "straight" { 0 } else {
                  if stations.first() > stations.last() { 0 } else { 1 }
                }
                cuts.at(side) = contacts.at(contact)
              }
              at += length + entry.sep
            }
          }
        }
        let painted = carrier
        // Split at the very same sampled parameters, keeping the shaft on its curve.
        // CeTZ 0.5.1 shorten-to uses different cubic samples; snap-to cannot fix
        // both contacts without changing the subcurve. The carrier is one joined path.
        if paint and cuts != (none, none) {
          let (origin, closed, commands) = carrier.segments.first()
          let bezier = cetz.path-util.bezier
          let bounds = cuts
            .enumerate()
            .map(((side, cut)) => {
              if cut == none { return side * commands.len() }
              let (kind, ..args) = commands.at(cut.segment-index)
              let t = if kind == "l" {
                cut.distance / cetz.vector.dist(cut.previous-point, args.last())
              } else {
                bezier.cubic-t-for-distance(
                  cut.previous-point,
                  args.last(),
                  ..args.slice(0, 2),
                  cut.distance,
                  samples: cetz.path-util.number-of-samples(auto),
                )
              }
              cut.segment-index + calc.clamp(t, 0, 1)
            })
          painted.segments = ()
          if bounds.first() < bounds.last() {
            let first = calc.floor(bounds.first())
            let last = calc.ceil(bounds.last()) - 1
            if first > 0 { origin = commands.at(first - 1).last() }
            commands = commands.slice(first, last + 1)
            let end = bounds.last() - last
            // Both ends may be trimmed within this cubic; rescale the second split.
            let start = (
              (bounds.first() - first) / if first == last { end } else { 1 }
            )
            for (side, t) in (end, start).enumerate() {
              let index = if side == 0 { commands.len() - 1 } else { 0 }
              let previous = if index == 0 { origin } else {
                commands.at(index - 1).last()
              }
              let (kind, ..args) = commands.at(index)
              if kind == "c" {
                let (s, e, c1, c2) = bezier
                  .split(previous, args.last(), ..args.slice(0, 2), t)
                  .at(side)
                if side == 1 { origin = s }
                args = (c1, c2, e)
              } else if side == 0 {
                args = (cetz.vector.lerp(previous, args.last(), t),)
              } else { origin = cetz.vector.lerp(previous, args.last(), t) }
              commands.at(index) = (kind, ..args)
            }
            painted.segments = ((origin, closed, commands),)
          }
        }
        if not paint {
          painted.stroke = none
          painted.fill = none
        }
        let drawables = (
          (painted,) + cetz.drawable.apply-tags(marks, cetz.drawable.TAG.mark)
        )
        if not flat {
          drawables = cetz.drawable.apply-transform(transform, drawables)
        }
        (
          ..plain,
          ctx: plain.ctx + (transform: transform),
          anchors: plain.anchors.with(transform: transform),
          drawables: drawables,
        )
      },
    ),
  )
}

#let _derived-path-elements(path, style, phase, anchor-start, anchor-end) = {
  let segments = curve-api.segments(path)
  if segments.len() == 0 {
    return (elements: (), length: 0)
  }
  if _has-mark(style) {
    if _has-pattern(style) {
      let painted = _segments-elements(
        segments,
        _without-mark-style(style),
        phase,
        anchor-start,
        anchor-end,
      )
      (
        elements: painted.elements + _mark-carrier-elements(path, style),
        length: painted.length,
      )
    } else {
      (
        elements: _mark-carrier-elements(path, style, paint: true),
        length: curve-api.length(path, accuracy: _style-value(
          style,
          "accuracy",
        )),
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

#let _path-mid-frame(path, accuracy, shift: 0) = {
  let segments = curve-api.segments(path)
  if segments.len() == 0 {
    return none
  }
  let total = curve-api.length(path, accuracy: accuracy)
  let at = calc.clamp(total / 2 + shift, 0, total)
  let (segment, t) = if at == 0 {
    (segments.first(), 0)
  } else if at == total {
    (segments.last(), 1)
  } else if total <= accuracy {
    (segments.first(), 0)
  } else {
    let prefix = curve-api.segments(curve-api.trim(
      path, end-outset: total - at, accuracy: accuracy,
    ))
    if prefix.len() == 0 { (segments.first(), 0) } else { (prefix.last(), 1) }
  }
  (
    point: if t == 0 { segment.start } else { segment.end },
    tangent: curve-api.cubic-tangent(segment, t),
  )
}

#let _label-side(style, frame, label-pos, label-origin) = {
  let side = _style-value(style, "label-side")
  if side == auto {
    if label-pos == none or label-origin == none {
      if _style-value(style, "offset") < 0 { -1 } else { 1 }
    } else {
      let toward = _point-sub(_point(label-pos), _point(label-origin))
      if _point-length(toward) <= 1e-9 {
        return if _style-value(style, "offset") < 0 { -1 } else { 1 }
      }
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

// Collision footprints use actual rendered shafts and marks. Keep separate
// pieces: the empty space between an arrow and its text is not occupied.
#let _annotation-bounds(ctx, elements) = {
  let painted = cetz.process.many(ctx, elements.flatten(), compute-bounds: false)
  let arrow-bounds = ()
  for drawable in painted.drawables {
    if cetz.drawable.TAG.hidden in drawable.tags { continue }
    if drawable.type == "path" and drawable.at("stroke", default: none) == none and drawable.at("fill", default: none) == none { continue }
    let stroke = cetz.util.resolve-stroke(drawable.at("stroke", default: none))
    let radius = cetz.util.resolve-number(ctx, stroke.thickness) / 2
    let pieces = ()
    if drawable.type == "path" {
      if cetz.drawable.TAG.mark in drawable.tags {
        pieces.push(cetz.path-util.bounds(drawable.segments))
      } else {
        for (start, _, commands) in drawable.segments {
          for command in commands {
            pieces.push(cetz.path-util.bounds(((start, false, (command,)),)))
            start = command.last()
          }
        }
      }
    } else if drawable.type == "content" {
      let (x, y, _) = drawable.pos
      pieces.push(((x - drawable.width / 2, y - drawable.height / 2),
        (x + drawable.width / 2, y + drawable.height / 2)))
    }
    for points in pieces {
      if points.len() == 0 { continue }
      arrow-bounds.push((
        left: calc.min(..points.map(p => p.at(0))) - radius,
        right: calc.max(..points.map(p => p.at(0))) + radius,
        bottom: calc.min(..points.map(p => p.at(1))) - radius,
        top: calc.max(..points.map(p => p.at(1))) + radius,
      ))
    }
  }
  arrow-bounds
}

// Flatten only the finite attachment geometry. Affine transformation happens
// before the flatness test, so rotated text and scaled canvases use one metric.
#let _label-path-lines(ctx, path, accuracy) = {
  let lines = ()
  if path == none { return lines }
  let pending = curve-api.segments(path).rev().map(segment => (
    points: (segment.start, segment.control-start, segment.control-end, segment.end).map(
      point => _point(cetz.matrix.mul4x4-vec3(ctx.transform, (..point, 0))),
    ),
    depth: 0,
  ))
  while pending.len() > 0 {
    let current = pending.pop()
    let (a, b, c, d) = current.points
    let chord = _point-sub(d, a)
    let squared = chord.at(0) * chord.at(0) + chord.at(1) * chord.at(1)
    let distance = point => {
      let relative = _point-sub(point, a)
      let at = if squared <= 1e-18 { 0 } else {
        calc.clamp((relative.at(0) * chord.at(0) + relative.at(1) * chord.at(1)) / squared, 0, 1)
      }
      _point-distance(point, _point-lerp(a, d, at))
    }
    let flatness = calc.max(distance(b), distance(c))
    if flatness <= accuracy or current.depth >= 10 {
      lines.push((a, d))
    } else {
      let ab = _point-lerp(a, b, 0.5)
      let bc = _point-lerp(b, c, 0.5)
      let cd = _point-lerp(c, d, 0.5)
      let abc = _point-lerp(ab, bc, 0.5)
      let bcd = _point-lerp(bc, cd, 0.5)
      let mid = _point-lerp(abc, bcd, 0.5)
      pending.push((points: (mid, bcd, cd, d), depth: current.depth + 1))
      pending.push((points: (a, ab, abc, mid), depth: current.depth + 1))
    }
  }
  lines
}

// Find actual clearance from a measured quad to the finite attachment path.
// The arc frame is unchanged; only its outward normal offset is adjusted.
#let _label-attachment-offset(lines, corners, frame, outward, gap, initial) = {
  if lines.len() == 0 { return initial }
  let dot = (a, b) => _point-x(a) * _point-x(b) + _point-y(a) * _point-y(b)
  let cross = (a, b) => _point-x(a) * _point-y(b) - _point-y(a) * _point-x(b)
  let direction = _point(outward)
  let scale = _point-length(direction)
  if scale <= 1e-9 { return initial }
  let normal = _point-scale(direction, 1 / scale)
  let tangent = (-normal.at(1), normal.at(0))
  let project = point => (dot(point, tangent), dot(point, normal))
  let quad = (corners.at(0), corners.at(1), corners.at(3), corners.at(2)).map(project)
  let orientation = cross(_point-sub(quad.at(1), quad.at(0)), _point-sub(quad.at(2), quad.at(1)))
  let left = calc.min(..quad.map(point => point.at(0)))
  let right = calc.max(..quad.map(point => point.at(0)))
  let radius = gap * scale
  // Intersect the vertical normal ray with a finite capsule. Its rectangle and
  // endpoint discs give exact clearance to each straight boundary segment.
  let capsule = (a, b) => {
    if calc.min(a.at(0), b.at(0)) > radius or calc.max(a.at(0), b.at(0)) < -radius { return none }
    let low = calc.inf
    let high = -calc.inf
    for point in (a, b) {
      if calc.abs(point.at(0)) <= radius {
        let reach = calc.sqrt(calc.max(0, radius * radius - point.at(0) * point.at(0)))
        low = calc.min(low, point.at(1) - reach)
        high = calc.max(high, point.at(1) + reach)
      }
    }
    let delta = _point-sub(b, a)
    let length = _point-length(delta)
    if length > 1e-12 {
      let along = _point-scale(delta, 1 / length)
      let across = (-along.at(1), along.at(0))
      let bottom = -calc.inf
      let top = calc.inf
      for (axis, minimum, maximum) in ((along, 0, length), (across, -radius, radius)) {
        let intercept = dot(a, axis)
        let slope = axis.at(1)
        if calc.abs(slope) <= 1e-12 {
          if -intercept < minimum or -intercept > maximum {
            bottom = calc.inf
            top = -calc.inf
            break
          }
        } else {
          let first = (minimum + intercept) / slope
          let last = (maximum + intercept) / slope
          bottom = calc.max(bottom, calc.min(first, last))
          top = calc.min(top, calc.max(first, last))
        }
      }
      if bottom <= top {
        low = calc.min(low, bottom)
        high = calc.max(high, top)
      }
    }
    if low <= high { (low, high) } else { none }
  }
  let intervals = ()
  for (start, end) in lines {
    let start = project(_point-sub(start, frame))
    let end = project(_point-sub(end, frame))
    if calc.min(start.at(0), end.at(0)) > right + radius or calc.max(start.at(0), end.at(0)) < left - radius { continue }
    // Sweep the measured quad along this finite segment. Each quad edge is
    // exposed at one endpoint; side transitions add the corner trajectories.
    let first = quad.map(corner => _point-sub(start, corner))
    let last = quad.map(corner => _point-sub(end, corner))
    let delta = _point-sub(end, start)
    let at-end = range(4).map(i => cross(
      _point-sub(first.at(calc.rem(i + 1, 4)), first.at(i)), delta,
    ) * orientation < 0)
    let boundary = ()
    for i in range(4) {
      let next = calc.rem(i + 1, 4)
      let points = if at-end.at(i) { last } else { first }
      boundary.push((points.at(i), points.at(next)))
      if calc.abs(orientation) <= 1e-18 { boundary.push((last.at(i), last.at(next))) }
      if at-end.at(i) != at-end.at(calc.rem(i + 3, 4)) or calc.abs(orientation) <= 1e-18 {
        boundary.push((first.at(i), last.at(i)))
      }
    }
    let low = calc.inf
    let high = -calc.inf
    for (a, b) in boundary {
      let interval = capsule(a, b)
      if interval != none {
        low = calc.min(low, interval.at(0))
        high = calc.max(high, interval.at(1))
      }
    }
    if low <= high { intervals.push((low, high)) }
  }
  // Curves can give disjoint collision intervals. Leave the component touching
  // the attachment frame, rather than jumping over a clear gap to a remote bend.
  let reach = 0
  let attached = false
  for (low, high) in intervals.sorted(key: interval => interval.at(0)) {
    if high < 0 { continue }
    if low > reach + 1e-10 { break }
    attached = true
    reach = calc.max(reach, high)
  }
  // Explicit independent arrow/label shifts may separate the whole finite
  // arrow from this normal ray; retain that deliberate tangential placement.
  if attached { reach / scale } else { initial }
}

#let _layer-label-element(ctx, path, style, label-pos, data, href: none, fixed-position: none, placement: false) = {
  let label = _layer-label(style, data)
  if label == none { return none }
  if href != none { label = link(href, label) }
  let edge = data.at("edge", default: none)
  let label-origin = if edge == none { none } else { edge.at("pos", default: none) }
  let accuracy = _style-value(style, "accuracy")
  let shift = _style-value(style, "label-shift")
  let frame = if fixed-position != none {
    (point: _point(fixed-position), tangent: (1, 0))
  } else { _path-mid-frame(path, accuracy, shift: shift) }
  if frame == none { return none }
  let attached = _style-value(style, "label-path")
  let side = _label-side(
    if attached == none { style } else { style + (offset: _style-value(attached, "offset")) },
    frame,
    if attached != none and _style-value(attached, "offset-side") != "label" { none } else { label-pos },
    label-origin,
  )
  let label-style = _style(_style-value(style, "label-style"), data)
  let gap = calc.max(0, _style-value(style, "label-gap"))
  // Explicit anchors are deliberate placement, not a request for box clearance.
  let clear-box = label-style.at("anchor", default: auto) in (auto, "auto")
  if clear-box { label-style.anchor = "center" }
  // CeTZ owns text bounds, wrapping, padding and rotation; measure its actual box.
  let element = cetz.draw.content((0, 0), label, padding: 0, ..label-style)
  let measured = element.first()(ctx)
  let origin = cetz.matrix.mul4x4-vec3(ctx.transform, (0, 0, 0))
  let corners = ()
  for anchor in ("north-west", "north-east", "south-west", "south-east") {
    corners.push(cetz.vector.sub((measured.anchors)(anchor), origin))
  }
  let total = if path == none { 0 } else { curve-api.length(path, accuracy: accuracy) }
  let preferred = calc.clamp(total / 2 + shift, 0, total)
  let movable = placement and fixed-position == none and clear-box and _style-value(style, "label-slide") and total > accuracy
  let sides = (side,)
  if movable and _style-value(style, "label-side") == auto { sides.push(-side) }
  // Momentum annotations share one full carrier. The short arrow is painted
  // only after choosing the label position, with exactly the same arc shift.
  let paths = sides.map(candidate-side => if attached == none { path } else {
    _path-layer(path, attached + (
      offset: calc.abs(_style-value(attached, "offset")) * candidate-side,
      offset-side: none, length: none, ratio: none, resolve-length: "none", shift: 0,
    ), 0, 0, none, auto)
  })
  let carrier-lines = if clear-box and fixed-position == none {
    _label-path-lines(ctx, path, accuracy)
  } else { () }
  let path-options = paths.enumerate().map(((path-index, candidate-path)) => {
    let candidate-side = sides.at(path-index)
    let total = if candidate-path == none { 0 } else { curve-api.length(candidate-path, accuracy: accuracy) }
    let preferred = calc.clamp(total / 2 + shift, 0, total)
    let low = 0
    let high = total
    if movable and attached != none {
      let centered = _path-layer(candidate-path, attached + (offset: 0, offset-side: none, shift: 0), 0, 0, none, auto)
      let half-length = curve-api.length(centered, accuracy: accuracy) / 2
      let relative-shift = _style-value(attached, "shift") - shift
      // Bound the shared displacement so neither end of the arrow is clamped
      // independently of its label. Deliberate pinned placements still clamp.
      low = calc.clamp(half-length - relative-shift, 0, total)
      high = calc.clamp(total - half-length - relative-shift, low, total)
      preferred = calc.clamp(preferred, low, high)
    }
    (
      total: total, preferred: preferred, low: low, high: high,
      lines: carrier-lines,
    )
  })
  let candidate-at = (path-index, at, sampled-frame: none, base: none) => {
    let candidate-path = paths.at(path-index)
    let candidate-side = sides.at(path-index)
    let options = path-options.at(path-index)
    let total = options.total
    let preferred = options.preferred
    let at = calc.clamp(at, options.low, options.high)
    let frame = if sampled-frame != none { sampled-frame } else if fixed-position != none { frame } else {
      _path-mid-frame(candidate-path, accuracy, shift: at - total / 2)
    }
    let tangent-length = _point-length(frame.tangent)
    let normal = if tangent-length <= 1e-9 { (0, candidate-side) } else {
      _point-scale((-_point-y(frame.tangent), _point-x(frame.tangent)), candidate-side / tangent-length)
    }
    let outward = cetz.vector.sub(
      cetz.matrix.mul4x4-vec3(ctx.transform, (..normal, 0)), origin,
    )
    let normal-squared = cetz.vector.dot(outward, outward)
    let nearest = if clear-box and normal-squared > 1e-18 {
      calc.min(..corners.map(corner => cetz.vector.dot(corner, outward))) / normal-squared
    } else { 0 }
    let path-shift = at - total / 2 + if attached == none { 0 } else { _style-value(attached, "shift") - shift }
    let arrow-bounds = ()
    let attachment-lines = options.lines
    if attached != none {
      let arrow-style = attached + (offset: 0, offset-side: none, shift: path-shift)
      let arrow-path = _path-layer(candidate-path, arrow-style, 0, 0, none, auto)
      if clear-box and fixed-position == none {
        // The finite arrow and its physical carrier are one attachment:
        // a wide label must not cross the carrier beyond the short arrow.
        attachment-lines = _label-path-lines(ctx, arrow-path, accuracy) + carrier-lines
      }
      if placement {
        arrow-bounds = _annotation-bounds(ctx, _derived-path-elements(
          arrow-path, arrow-style, auto, true, true,
        ).elements)
      }
    }
    let normal-offset = if clear-box and fixed-position == none {
      _label-attachment-offset(
        attachment-lines, corners,
        _point(cetz.matrix.mul4x4-vec3(ctx.transform, (..frame.point, 0))),
        outward, gap, gap - nearest,
      )
    } else { gap - nearest }
    let position = if base != none and (fixed-position != none or normal-offset == gap - nearest) {
      base.position
    } else if fixed-position != none { frame.point } else {
      _point-add(frame.point, _point-scale(normal, normal-offset))
    }
    let center = cetz.matrix.mul4x4-vec3(ctx.transform, (..position, 0))
    let bounds = corners.map(corner => cetz.vector.add(center, corner))
    (
      position: position,
      at: at,
      side: candidate-side,
      path-index: path-index,
      path-shift: path-shift,
      arrow-bounds: arrow-bounds,
      corners: bounds,
      bounds: if base != none and position == base.position { base.bounds } else {
        (
          left: calc.min(..bounds.map(p => p.at(0))),
          right: calc.max(..bounds.map(p => p.at(0))),
          bottom: calc.min(..bounds.map(p => p.at(1))),
          top: calc.max(..bounds.map(p => p.at(1))),
        )
      },
      // Prefer the original side when both placements are otherwise equivalent.
      cost: if base != none { base.cost } else { (
        (if total <= accuracy { 0 } else { 0.002 * calc.pow((at - preferred) / total, 2) })
          + if candidate-side == side { 0 } else { 0.00002 }
      ) },
    )
  }
  let candidates = ()
  for (path-index, options) in path-options.enumerate() {
    let positions = (options.preferred,)
    if movable {
      // Stay inside the carrier; endpoints belong to the vertex labels.
      positions += range(3, 30).map(i => calc.clamp(options.total * i / 32, options.low, options.high))
        .dedup().filter(at => calc.abs(at - options.preferred) > accuracy)
    }
    let candidate-path = paths.at(path-index)
    let candidate-side = sides.at(path-index)
    // Batch prefix trims and coordinate transforms before attaching the finite
    // arrow footprint and enforcing its shared clearance with the text.
    let frames = if fixed-position != none { (frame,) } else {
      curve-api.frames(candidate-path, positions.map(at => options.total / 2 + (at - options.total / 2)), accuracy: accuracy)
    }
    let bases = cbor(_plugin.label_candidates(cbor.encode((
      frames: frames.map(frame => (point: frame.point.map(float), tangent: frame.tangent.map(float))),
      positions: positions.map(float), transform: ctx.transform.map(row => row.map(float)),
      origin: origin.map(float), corners: corners.map(corner => corner.map(float)),
      total: float(options.total), preferred: float(options.preferred),
      accuracy: float(accuracy), side: float(candidate-side), preferred-side: float(side),
      path-index: path-index,
      path-shift: float(if attached == none { 0 } else { _style-value(attached, "shift") - shift }),
      gap: float(gap), clear-box: clear-box, fixed: fixed-position != none,
    ))))
    for (at, sampled-frame, base) in positions.zip(frames, bases) {
      candidates.push(candidate-at(path-index, at, sampled-frame: sampled-frame, base: base))
    }
  }
  if attached == none and sides.len() == 2 {
    let count = calc.quo(candidates.len(), 2)
    candidates = range(count).map(i => (candidates.at(i), candidates.at(i + count))).flatten()
  }
  if placement {
    (label: label, style: label-style, candidates: candidates, edge: data.at("eid", default: none), paths: paths, path-style: attached,
      candidate-at: candidate-at, ranges: path-options.map(options => (options.low, options.high)))
  } else {
    cetz.draw.content(candidates.first().position, label, padding: 0, ..label-style)
  }
}

// Base pair padding is shared between the two boxes. Extra label padding
// enlarges only label boxes; the collision overlay uses the same clearances.
#let _label-collision-padding = (labels: 0.35, obstacles: 0.08, extra: 0.6)

// Use the painted paths, including waves, coils and crossing gaps. Collision
// scoring must not turn their empty surrounding space into an edge obstacle.
// Native dash styles use the continuous stroke envelope, as for a solid edge.
#let _edge-collision-lines(ctx, elements) = {
  let painted = cetz.process.many(ctx, elements.flatten(), compute-bounds: false)
  let lines = ()
  let canvas = ctx + (transform: cetz.matrix.ident(4))
  for drawable in painted.drawables {
    if drawable.type != "path" or cetz.drawable.TAG.hidden in drawable.tags { continue }
    let stroke = cetz.util.resolve-stroke(drawable.at("stroke", default: none))
    if drawable.at("stroke", default: none) == none { continue }
    let radius = cetz.util.resolve-number(ctx, stroke.thickness) / 2
    for (origin, closed, commands) in drawable.segments {
      let start = origin
      for command in commands {
        let pieces = ()
        if command.first() == "l" {
          for end in command.slice(1) {
            pieces.push((_point(start), _point(end)))
            start = end
          }
        } else if command.first() == "c" {
          let segment = (start: _point(start), control-start: _point(command.at(1)),
            control-end: _point(command.at(2)), end: _point(command.at(3)))
          pieces = _label-path-lines(canvas, _segments-path((segment,)), 0.005)
          start = command.last()
        }
        lines += pieces.map(((a, b)) => (start: a, end: b, radius: radius))
      }
      if closed and start != origin {
        lines.push((start: _point(start), end: _point(origin), radius: radius))
      }
    }
  }
  lines
}

// Clip physical edge segments to the measured, possibly rotated text quad.
// Only the stroke's radius expands this box (conservatively at corners):
// no label collision padding and
// no force at all outside an actual intersection. Integrating covered length
// makes the penalty independent of how finely the edge was subdivided.
#let _label-edge-intersection(candidate, lines) = {
  if lines.len() == 0 { return 0 }
  let corners = candidate.corners
  let origin = _point(corners.at(0))
  let x = _point-sub(corners.at(1), origin)
  let y = _point-sub(corners.at(2), origin)
  let determinant = x.at(0) * y.at(1) - x.at(1) * y.at(0)
  if calc.abs(determinant) < 1e-12 { return 0 }
  let width = _point-length(x)
  let height = _point-length(y)
  let local = point => {
    let delta = _point-sub(point, origin)
    ((delta.at(0) * y.at(1) - delta.at(1) * y.at(0)) / determinant,
      (x.at(0) * delta.at(1) - x.at(1) * delta.at(0)) / determinant)
  }
  let length = 0
  for line in lines {
    let box = candidate.bounds
    if (calc.max(line.start.at(0), line.end.at(0)) + line.radius < box.left
      or calc.min(line.start.at(0), line.end.at(0)) - line.radius > box.right
      or calc.max(line.start.at(1), line.end.at(1)) + line.radius < box.bottom
      or calc.min(line.start.at(1), line.end.at(1)) - line.radius > box.top) { continue }
    let a = local(line.start)
    let b = local(line.end)
    let low = 0
    let high = 1
    for axis in range(2) {
      // The inverse transform's row norm also handles sheared canvases.
      let padding = line.radius * (if axis == 0 { height } else { width }) / calc.abs(determinant)
      let delta = b.at(axis) - a.at(axis)
      if calc.abs(delta) < 1e-12 {
        if a.at(axis) < -padding or a.at(axis) > 1 + padding { high = -1; break }
      } else {
        let from = (-padding - a.at(axis)) / delta
        let to = (1 + padding - a.at(axis)) / delta
        low = calc.max(low, calc.min(from, to))
        high = calc.min(high, calc.max(from, to))
      }
    }
    length += calc.max(0, high - low) * _point-distance(line.start, line.end)
  }
  6 * length / calc.max(1e-9, calc.min(width, height))
}

// Optimize arc length and automatic side choices with fixed normal clearance.
// Joint proposals let annotations pass one another without requiring a single
// temporarily worse move; a short joint descent settles remaining collisions.
#let _relax-label-placements(placements, obstacles, label-padding: _label-collision-padding.extra, fixed-arrows: (), edge-lines: ()) = {
  assert(label-padding >= 0, message: "label-collision-padding must be non-negative")
  if placements.all(label => label.candidates.len() == 1) {
    return placements.map(label => label.candidates.first())
  }
  let coordinated = placements.any(label => label.candidates.any(candidate => (
    candidate.keys().contains("arrow-bounds")
  )))
  // Normal clearance handles ordinary carriers. Only actual text intersections,
  // vertices and self-loop bends contribute static annotation penalties.
  let obstacles = if coordinated { obstacles.filter(obstacle => (
    obstacle.at("edge", default: none) == none or obstacle.at("self-loop", default: false)
  )) } else { obstacles }
  let geometry = cbor.encode((
    placements: placements.map(label => (
      edge: label.at("edge", default: none),
      candidates: label.candidates.map(candidate => (
        bounds: candidate.bounds,
        cost: candidate.cost + if edge-lines.len() > 0 and candidate.keys().contains("corners") {
          _label-edge-intersection(candidate, edge-lines)
        } else { 0 },
        arrow-bounds: candidate.at("arrow-bounds", default: ()),
      )),
    )),
    obstacles: obstacles, label-padding: label-padding,
    fixed-arrows: fixed-arrows, coordinated: coordinated,
    obstacle-padding: _label-collision-padding.obstacles,
  ))
  let math = (
    temperatures: range(84).map(sweep => 0.15 * calc.pow(0.9, sweep)),
    exponentials: (), pair-padding: _label-collision-padding.labels,
  )
  let result = cbor(_plugin.label_search(geometry, cbor.encode(math)))
  // Native exp is only a prediction. Its value only affects comparisons with
  // RNG thresholds, so Typst can certify all decisions through their interval.
  // A mismatch replays with exact host values, extending the verified prefix.
  while result.missing.len() > 0 {
    let verified = true
    for check in result.missing {
      let probability = calc.exp(check.argument)
      math.exponentials.push((check.argument, probability))
      // Negate the original comparison to certify every rejection.
      verified = (verified
        and (check.lower == none or check.lower < probability)
        and (check.upper == none or not (check.upper < probability)))
    }
    if verified { break }
    result = cbor(_plugin.label_search(geometry, cbor.encode(math)))
  }
  placements.enumerate().map(((i, label)) => {
    label.candidates.at(result.choices.at(i)) + (cost: result.costs.at(i))
  })
}

#let _paired-layer-label-element(
  ctx,
  halves,
  source-style,
  sink-style,
  label-pos,
  data,
  href: none,
  placement: false,
) = {
  let source-label = _layer-label(source-style, data)
  let sink-label = _layer-label(sink-style, data)
  let style = if source-label != none { source-style } else { sink-style }
  if style == none {
    return none
  }
  if source-label != none and sink-label != none {
    let source-path = _style-value(source-style, "label-path")
    let sink-path = _style-value(sink-style, "label-path")
    if source-path != none and sink-path != none and not _has-mark(source-path) and _has-mark(sink-path) {
      style.label-path = source-path + (mark: sink-path.mark)
    }
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
    _layer-label-element(ctx, _segments-path(segments), style, label-pos, data, href: href, placement: placement)
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
  let pattern-whole = halves.at("pattern-whole", default: whole)
  if (
    pattern-whole != none and _same-pattern-geometry(source-style, sink-style)
  ) {
    let path = _segments-path(pattern-whole)
    let total = curve-api.length(path, accuracy: _style-value(
      source-style,
      "accuracy",
    ))
    elements = _path-window-elements(path, source-style, _half-paint-windows(
      halves,
      source-style,
      sink-style,
      total,
    ))
    // Marks retain their original carrier and placement, above both painted halves.
    for (segments, style) in (
      (halves.source, source-style),
      (halves.sink, sink-style),
    ) {
      elements += _mark-carrier-elements(_segments-path(segments), style)
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

#let _dangling-path(
  start, end, bend, style, dangling-at-start: false,
  route-points: (), node-point: none, node-outset: 0, omega: 1.0, carrier: none,
) = {
  if (carrier != none or route-points.len() > 0) and _style-value(style, "route-points") != "ignore" {
    let start = if not dangling-at-start and node-point != none { _point(node-point) } else { _point(start) }
    let end = if dangling-at-start and node-point != none { _point(node-point) } else { _point(end) }
    let route = if dangling-at-start { route-points.rev() } else { route-points }
    // Explicit box anchors override only the carrier endpoint; interior knots
    // still follow the accepted route. Automatic anchors resolve to node centers.
    let points = if carrier == none { (start, ..route.map(_point), end) } else {
      (start, ..carrier.slice(1, carrier.len() - 1), end)
    }
    let path = _routed-split-through(points, omega: omega).curve
    let segments = curve-api.segments(path)
    if segments.len() > 0 and _style-value(style, "dangling-tangent") != auto {
      let index = if dangling-at-start { 0 } else { segments.len() - 1 }
      let segment = segments.at(index)
      // Keep the interior Hobby tangent while imposing the existing external
      // endpoint tangent on its adjacent span only.
      let endpoint = curve-api.segments(_dangling-path(
        segment.start, segment.end, none, style, dangling-at-start: dangling-at-start,
      )).first()
      if dangling-at-start { segment.control-start = endpoint.control-start }
      else { segment.control-end = endpoint.control-end }
      segments.at(index) = segment
      path = curve-api.path(..segments.map(curve-api.from-cubic))
    }
    return _trim-routed-path(path,
      start-outset: if dangling-at-start { 0 } else { node-outset },
      end-outset: if dangling-at-start { node-outset } else { 0 },
    )
  }

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

#let _cut-path-elements(
  path,
  style,
  cuts,
  mark-style: none,
  paint-windows: none,
) = {
  let accuracy = _style-value(style, "accuracy")
  let total = curve-api.length(path, accuracy: accuracy)
  let gap = calc.max(0, _style-value(style, "crossing-gap"))
  let windows = if paint-windows == none {
    ((start: 0, end: total, style: style),)
  } else { paint-windows }
  let visible = ()
  let cursor = 0
  if gap > 0 {
    for center in cuts.sorted() {
      let cut-start = _clamp(center - gap / 2, cursor, total)
      if cut-start > cursor { visible.push((start: cursor, end: cut-start)) }
      cursor = _clamp(center + gap / 2, cursor, total)
    }
  }
  if cursor < total { visible.push((start: cursor, end: total)) }
  let clipped = ()
  for window in windows {
    for interval in visible {
      let start = calc.max(window.start, interval.start)
      let end = calc.min(window.end, interval.end)
      if end > start { clipped.push(window + (start: start, end: end)) }
    }
  }
  let elements = _path-window-elements(path, style, clipped)
  if mark-style != none { elements += _mark-carrier-elements(path, mark-style) }
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

#let _subgraph-record(graph, entry, default-style) = {
  let selected = entry
  let style = default-style
  if (
    type(entry) == dictionary
      and entry.at("linnest-kind", default: none) != "linnest-subgraph"
  ) {
    selected = entry.at("subgraph", default: entry.at("graph", default: none))
    if selected == none {
      panic("draw: subgraph dictionary entries need a `subgraph` field")
    }
    style = entry.at(
      "edge-style",
      default: entry.at(
        "style",
        default: entry.at("subgraph-edge-style", default: default-style),
      ),
    )
  }
  (
    hedges: subgraph-api.hedges(subgraph-api._impl.validate(graph, selected)),
    edge-style: style,
  )
}

#let _subgraph-records(graph, subgraph, default-style) = {
  if subgraph == none {
    ()
  } else if type(subgraph) == array {
    subgraph
      .filter(entry => entry != none)
      .map(entry => _subgraph-record(graph, entry, default-style))
  } else {
    (_subgraph-record(graph, subgraph, default-style),)
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
  style != none and not _style-value(style, "label-only")
    and style.at("crossing-under", default: none) != none
)

#let _has-crossing-layer(layers) = layers.any(_has-crossing-style)

#let _has-visible-stroke(style) = if style == none or _style-value(style, "label-only") {
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
      route-points: _route-points(source-half-edge),
      carrier: _data-fields(edge).at("layout-carrier", default: none),
      node-point: if source-anchor == auto { _node-pos(nodes.at(source-half-edge.node)) } else { none },
      node-outset: if source-anchor == auto { node-outsets.at(source-half-edge.node) } else { 0 },
      omega: edge-omega,
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
      route-points: _route-points(sink-half-edge),
      carrier: _data-fields(edge).at("layout-carrier", default: none),
      node-point: if sink-anchor == auto { _node-pos(nodes.at(sink-half-edge.node)) } else { none },
      node-outset: if sink-anchor == auto { node-outsets.at(sink-half-edge.node) } else { 0 },
      omega: edge-omega,
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
  let edge-offset = options.edge-offset
  let edge-length = options.edge-length
  let edge-ratio = options.edge-ratio
  let edge-resolve-length = options.edge-resolve-length
  let edge-accuracy = options.edge-accuracy
  let edge-optimize = options.edge-optimize
  let edge-split-gap = options.edge-split-gap
  let edge-dangling-tangent = options.edge-dangling-tangent
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
    (
      ctx => {
        let nodes = graph-api.nodes(graph)
        let edges = graph-api.edges(graph)
        let elements = ()
        let node-elements = ()
        let node-outsets = ()
        let node-boxes = ()
        let node-targets = ()
        let edge-targets = ()
        let label-placements = ()
        let fixed-arrows = ()
        let edge-lines = ()
        let label-obstacles = ()
        let debug-level = _debug-level(debug)
        let subgraph-records = _subgraph-records(
          graph,
          subgraph,
          subgraph-edge-style,
        )

        for (i, v) in nodes.enumerate() {
          let boundary = v.at("boundary", default: none) != none
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
          let label = if boundary { none } else {
            _content(node-label, node-data, default-label)
          }
          let node-style-value = if boundary {
            (radius: 0, fill: none, stroke: none)
          } else {
            (
              (
                radius: node-radius,
                fill: node-fill,
                stroke: node-stroke,
              )
                + _style(graph-node-style, node-data)
                + _style(node-style, node-data)
            )
          }
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
          let node-label-draw-style = if boundary { (:) } else {
            _style(graph-node-label-style, node-data) + node-label-style
          }
          let layout-width = _statement-number(v, "layout-width", default: none)
          let layout-height = _statement-number(
            v,
            "layout-height",
            default: none,
          )
          let node-width = if boundary { 0 } else if layout-width == none {
            2 * node-style.radius
          } else { layout-width }
          let node-height = if boundary { 0 } else if layout-height == none {
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
          if not boundary {
            let corners = ()
            for sign in ((-1, -1), (-1, 1), (1, -1), (1, 1)) {
              corners.push(cetz.matrix.mul4x4-vec3(ctx.transform, (
                _point-x(pos) + sign.at(0) * node-width / 2,
                _point-y(pos) + sign.at(1) * node-height / 2,
                0,
              )))
            }
            label-obstacles.push((
              left: calc.min(..corners.map(p => p.at(0))),
              right: calc.max(..corners.map(p => p.at(0))),
              bottom: calc.min(..corners.map(p => p.at(1))),
              top: calc.max(..corners.map(p => p.at(1))),
            ))
          }
          node-outsets.push(if boundary { 0 } else if draw-node == auto {
            _node-outset(node-style, node-outset)
          } else if node-outset == auto {
            calc.max(node-width, node-height) / 2
          } else {
            node-outset
          })

          if boundary {
            // Cut boundary nodes are endpoint positions, not interactions.
            node-elements.push(cetz.draw.anchor(box.name, box.center))
          } else if draw-node == auto {
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
          if not boundary {
            let details = (
              edges: edges.filter(e => (
                (e.source != none and e.source.node == i)
                  or (e.sink != none and e.sink.node == i)
              )).map(e => e.edge),
            )
            // Display transformations may retain the source graph's inspection identities.
            details += node-data.at("inspection", default: (:))
            let name = node-data.at("feynkit-name", default: node.name)
            if name != none {
              details.insert("name", str(name))
            }
            let label-size = if label == none { (width: 0pt, height: 0pt) } else {
              measure(label)
            }
            node-targets += _identity-target(
              pos,
              _identity-href("node", details.at("node", default: i), details),
              width: calc.max(10pt, node-width * ctx.length, label-size.width),
              height: calc.max(10pt, node-height * ctx.length, label-size.height),
            )
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
          let layers = _edge-style-layers(graph-style, options, record.edge-data)
          let hidden = layers.source.len() == 0 and layers.sink.len() == 0
          let ordinary-label = if hidden { none } else {
            _content(edge-label, record.edge-data, _as-content(record.data-label))
          }
          let edge-label-draw-style = (
            _style(graph-edge-label-style, record.edge-data)
              + _style(edge-label-style, record.edge-data)
          )
          for side in ("source", "sink") {
            layers.at(side) = layers.at(side).map(layer => {
              // With no label direction, the signed offset sets the preferred
              // label side. Keep automatic labels free to flip during placement.
              if (
                _style-value(layer, "offset-side") == "label"
                  and _point-distance(record.label-pos, record.edge.pos) <= 1e-9
              ) {
                layer.offset-side = none
                layer.offset = layer.at("offset", default: geometry-style.offset)
              }
              let label = _call(_style-value(layer, "label"), record.edge-data)
              layer = if label == auto {
                layer + (
                  label: ordinary-label,
                  label-style: edge-label-draw-style
                    + _style(_style-value(layer, "label-style"), record.edge-data),
                )
              } else { layer + (label: _as-content(label)) }
              // A momentum arrow remains visible when its text is disabled.
              let attached = _style-value(layer, "label-path")
              if layer.label == none and attached != none {
                layer + attached + (label: none, label-path: none, label-only: false)
              } else { layer }
            })
          }
          let source-style-layers = layers.source
          let sink-style-layers = layers.sink
          let has-attached-label = (
            (record.source-half-edge != none
              and _has-layer-label(source-style-layers, record.edge-data))
              or (record.sink-half-edge != none
                and _has-layer-label(sink-style-layers, record.edge-data))
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
                  ordinary-label
                },
                edge-label-draw-style: edge-label-draw-style,
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
        let crossing-paths = if has-crossings or edge-records.any(record => (
          record.ev-label != none and _statement-number(record.edge, "layout-label-gap") != none
        )) {
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
          let inspection = edge-data.at("inspection", default: (:))
          let details = (
            edge: edge.edge,
            source: if source-half-edge == none { none } else { source-half-edge.node },
            sink: if sink-half-edge == none { none } else { sink-half-edge.node },
            orientation: edge.orientation,
          ) + inspection
          if edge.name != none {
            details.insert("name", str(edge.name))
          }
          for key in ("particle", "pdg", "external-state", "external-index", "external-name", "momentum") {
            let value = edge-data.at(key, default: none)
            if type(value) in (str, int, float, bool) {
              details.insert(key, value)
            }
          }
          let href = _identity-href("edge", details.edge, details)
          let half-hrefs = ()
          let edge-hrefs = ()
          for (side, half-edge, pair) in (
            ("source", source-half-edge, sink-half-edge),
            ("sink", sink-half-edge, source-half-edge),
          ) {
            let owner = if half-edge == none { pair } else { half-edge }
            let flow = if half-edge != none { side } else if side == "source" { "sink" } else { "source" }
            let other = if flow == "source" { "sink" } else { "source" }
            let hedge = inspection.at(flow + "-hedge", default: owner.hedge)
            let half-details = details + ("half-edge": hedge)
            edge-hrefs.push(_identity-href("edge", details.edge, half-details))
            half-hrefs.push(_identity-href("halfedge", hedge, half-details + (
              node: owner.node,
              flow: flow,
              pair: inspection.at(other + "-hedge", default: if half-edge == none or pair == none { none } else { pair.hedge }),
            )))
          }
          let hrefs = (half-hrefs.at(0), edge-hrefs.at(0), edge-hrefs.at(1), half-hrefs.at(1))

          for layer-index in range(0, layer-count) {
            let layer-start = elements.len()
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
            } else if (
              source-in-subgraph and not subgraph-edge-underlay
                and not _style-value(source-style-value, "label-only")
            ) {
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
            } else if (
              sink-in-subgraph and not subgraph-edge-underlay
                and not _style-value(sink-style-value, "label-only")
            ) {
              sink-style-value + _last-style(sink-subgraph-styles)
            } else {
              sink-style-value
            }
            sink-draw-style = _orient-mark-style(
              sink-draw-style,
              edge-data,
              "sink",
            )
            let source-crossing-under = if not _has-crossing-style(source-draw-style) {
              none
            } else {
              _style-value(source-draw-style, "crossing-under")
            }
            let sink-crossing-under = if not _has-crossing-style(sink-draw-style) {
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
              edge-targets += _edge-identity-targets(
                ctx,
                ((halves.source, source-draw-style), (halves.sink, sink-draw-style)).map(
                  ((segments, style)) => (
                    segments: segments,
                    visible: style != none and (_has-visible-stroke(style) or _has-mark(style)),
                  ),
                ),
                hrefs,
              )
              if (
                source-style-value != none
                  and not _style-value(source-style-value, "label-only")
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
                  and not _style-value(sink-style-value, "label-only")
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
                let carrier = halves.at("pattern-whole", default: halves.whole)
                let same-geometry = if source-paint-style == none or sink-paint-style == none {
                  false
                } else if _has-pattern(source-paint-style) or _has-pattern(sink-paint-style) {
                  _same-pattern-geometry(source-paint-style, sink-paint-style)
                } else { _same-layer-geometry(source-paint-style, sink-paint-style) }
                if (
                  source-draw-style == none
                    or sink-draw-style == none
                    or source-crossing-under != sink-crossing-under
                    or source-crossing-gap != sink-crossing-gap
                    or carrier == none
                    or not same-geometry
                ) {
                  panic(
                    "draw: crossing-under on paired edge "
                      + str(edge-data.eid)
                      + " requires continuous source/sink geometry with the same pattern",
                  )
                }
                let path = _segments-path(carrier)
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
                  paint-windows: _half-paint-windows(halves,
                    source-paint-style, sink-paint-style,
                    curve-api.length(path, accuracy: _style-value(source-paint-style, "accuracy")),
                  ),
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
                href: href,
                placement: true,
              )
              if attached-label != none {
                label-placements.push(attached-label + (hrefs: hrefs))
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
                route-points: _route-points(source-half-edge),
                carrier: _data-fields(edge).at("layout-carrier", default: none),
                node-point: if source-anchor == auto { _node-pos(nodes.at(source-half-edge.node)) } else { none },
                node-outset: if source-anchor == auto { node-outsets.at(source-half-edge.node) } else { 0 },
                omega: edge-omega,
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
                  and not _style-value(source-style-value, "label-only")
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
                  not _has-crossing-style(draw-style)
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
                  href: href,
                  placement: true,
                )
                if attached-label != none {
                  label-placements.push(attached-label + (hrefs: hrefs))
                }
                if _has-visible-stroke(draw-style) or _has-mark(draw-style) {
                  edge-targets += _edge-identity-targets(
                    ctx, ((segments: curve-api.segments(visible-path), visible: true),), hrefs,
                  )
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
                route-points: _route-points(sink-half-edge),
                carrier: _data-fields(edge).at("layout-carrier", default: none),
                node-point: if sink-anchor == auto { _node-pos(nodes.at(sink-half-edge.node)) } else { none },
                node-outset: if sink-anchor == auto { node-outsets.at(sink-half-edge.node) } else { 0 },
                omega: edge-omega,
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
                  and not _style-value(sink-style-value, "label-only")
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
                  not _has-crossing-style(draw-style)
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
                  href: href,
                  placement: true,
                )
                if attached-label != none {
                  label-placements.push(attached-label + (hrefs: hrefs))
                }
                if _has-visible-stroke(draw-style) or _has-mark(draw-style) {
                  edge-targets += _edge-identity-targets(
                    ctx, ((segments: curve-api.segments(visible-path), visible: true),), hrefs,
                  )
                }
              }
            }
            // An unlabelled bounded offset layer is a fixed decoration, such
            // as an external momentum arrow whose text stays at the endpoint.
            // Other painted paths participate only when they intersect text.
            if (source-draw-style, sink-draw-style).any(style => (
              style != none and not _style-value(style, "label-only")
                and _layer-label(style, edge-data) == none
                and _style-value(style, "offset") != 0
                and (_style-value(style, "length") != none or _style-value(style, "ratio") != none)
            )) {
              fixed-arrows += _annotation-bounds(ctx, elements.slice(layer-start))
            } else {
              edge-lines += _edge-collision-lines(ctx, elements.slice(layer-start))
            }
          }

          // Only self-loop bends repel annotations: a remote part of the loop
          // can approach from another side. Ordinary carriers are handled by the
          // fixed normal gap. Waves and coils occupy a band, not just a line.
          let radius = 0.06 + calc.max(0, ..(source-style-layers + sink-style-layers).map(style => (
            if _style-value(style, "pattern") == none { 0 } else {
              calc.abs(_style-value(style, "pattern-amplitude"))
            }
          )))
          let pad-x = radius * (calc.abs(ctx.transform.at(0).at(0)) + calc.abs(ctx.transform.at(0).at(1)))
          let pad-y = radius * (calc.abs(ctx.transform.at(1).at(0)) + calc.abs(ctx.transform.at(1).at(1)))
          if self-loop {
            let obstacle-segments = (source-label-segments, sink-label-segments)
              .filter(segments => segments != none).flatten()
            label-obstacles += cbor(_plugin.label_obstacle_boxes(cbor.encode((
              segments: obstacle-segments, transform: ctx.transform.map(row => row.map(float)),
              pad-x: pad-x, pad-y: pad-y,
            )))).map(box => (edge: edge.edge, self-loop: self-loop, ..box))
          }

          if ev-label != none {
            let gap = _statement-number(edge, "layout-label-gap")
            let carrier = crossing-paths.ids.at(str(edge.edge), default: none)
            let follows-path = gap != none and carrier != none
            label-placements.push(_layer-label-element(
              ctx,
              if follows-path { carrier.path } else { none },
              (label: ev-label, label-style: edge-label-draw-style, label-gap: if gap == none { 0 } else { gap }),
              label-pos,
              edge-data,
              href: href,
              placement: true,
              fixed-position: if follows-path { none } else { label-pos },
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

        label-placements = label-placements.filter(label => label != none)
        let label-padding = options.label-collision-padding
        let relaxed = _relax-label-placements(label-placements, label-obstacles, label-padding: label-padding, fixed-arrows: fixed-arrows, edge-lines: edge-lines)
        for (label, candidate) in label-placements.zip(relaxed) {
          if label.path-style != none {
            let style = label.path-style + (offset: 0, offset-side: none, shift: candidate.path-shift)
            let path = _path-layer(label.paths.at(candidate.path-index), style, 0, 0, none, auto)
            elements += _derived-path-elements(path, style, auto, true, true).elements
            edge-targets += _edge-identity-targets(
              ctx, ((segments: curve-api.segments(path), visible: true),), label.hrefs,
            )
          }
          elements.push(cetz.draw.content(
            candidate.position, label.label, padding: 0, ..label.style,
          ))
        }

        for element in node-elements {
          elements.push(element)
        }

        if options.debug-label-collisions {
          // Bounds are already in canvas coordinates. Reset the transform only
          // inside this floating overlay, so it cannot change graph bounds.
          elements.push(cetz.draw.floating(cetz.draw.scope({
            cetz.draw.set-transform(none)
            for (boxes, padding, color, dashed) in (
              (label-obstacles, _label-collision-padding.obstacles / 2, rgb("#f59e0b"), false),
              (relaxed.map(candidate => candidate.bounds), _label-collision-padding.labels / 2 + label-padding, rgb("#a855f7"), true),
              (relaxed.map(candidate => candidate.bounds), _label-collision-padding.obstacles / 2 + label-padding, rgb("#06b6d4"), false),
              (relaxed.map(candidate => candidate.arrow-bounds).flatten() + fixed-arrows, _label-collision-padding.labels / 2, rgb("#22c55e"), true),
            ) {
              for box in boxes {
                cetz.draw.rect(
                  (box.left - padding, box.bottom - padding),
                  (box.right + padding, box.top + padding),
                  fill: color.transparentize(94%),
                  stroke: (paint: color, thickness: 0.4pt, dash: if dashed { "dashed" } else { "solid" }),
                )
              }
            }
          })))
        }

        // Measure the graph before overlays and reuse its processed drawables,
        // so custom drawing callbacks run once and named anchors remain available.
        let drawable-elements = ()
        for element in elements.flatten() {
          if element != none { drawable-elements.push(element) }
        }
        let rendered = cetz.process.many(
          ctx, drawable-elements,
          compute-bounds: type(options.draw-after) == function,
        )
        let overlay = options.draw-after
        if type(overlay) == function {
          let low = if rendered.bounds == none { (0, 0, 0) } else { rendered.bounds.low }
          let high = if rendered.bounds == none { (0, 0, 0) } else { rendered.bounds.high }
          overlay = overlay(graph, (
            left: low.at(0), right: high.at(0),
            bottom: low.at(1), top: high.at(1),
            width: high.at(0) - low.at(0), height: high.at(1) - low.at(1),
          ))
        }
        let after = cetz.process.many(
          rendered.ctx,
          if overlay == none { () } else { overlay.flatten().filter(element => element != none) },
          compute-bounds: false,
        )
        // Nodes come last so incident edge targets cannot intercept node hits.
        let targets = _draw-identity-targets(after.ctx, edge-targets + node-targets)
        (ctx: targets.ctx, drawables: rendered.drawables + after.drawables + targets.drawables)
      },
    ),
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
