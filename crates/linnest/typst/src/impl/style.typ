// Edge style resolution shared by `graph.style`, which measures what layout
// needs to know about drawing, and `draw`, which paints with it.

#let _style-key = "linnest-style"

#let _canvas-length(unit) = {
  if type(unit) in (int, float) {
    unit * 1em
  } else {
    unit
  }
}

#let _edge-geometry-defaults = (
  offset: 0,
  length: none,
  ratio: none,
  resolve-length: "min",
  shift: 0,
  accuracy: 0.001,
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

#let _style-defaults = (_mark-defaults + _pattern-defaults + _edge-label-defaults
  + _edge-routing-defaults + _edge-crossing-defaults + _edge-geometry-defaults)

#let _without-keys(style, keys) = {
  let clean = style
  for key in keys {
    if key in clean {
      let _ = clean.remove(key)
    }
  }
  clean
}

#let _without-pattern-style(style) = _without-keys(
  style,
  _pattern-defaults.keys(),
)

#let _draw-style(style) = _without-keys(style, _style-defaults.keys())

#let _without-mark-style(style) = {
  let clean = style
  if clean.keys().contains("mark") {
    let _ = clean.remove("mark")
  }
  clean
}

#let _style-value(style, key) = style.at(key, default: _style-defaults.at(key, default: none))

#let _edge-geometry(style) = _edge-geometry-defaults + style

#let _edge-routing(style) = _edge-routing-defaults + style

#let _pattern-style(style) = _pattern-defaults + style

#let _style-offset(style) = {
  _style-value(style, "offset")
}

#let _resolve-length-method(style) = {
  _style-value(style, "resolve-length")
}

// Waves and coils occupy a band, not just a line: the half width of an edge's
// painted band beside its carrier is their amplitude plus 0.06 for the stroke.
#let _painted-band(layers) = 0.06 + calc.max(0, ..layers.map(style => (
  if _style-value(style, "pattern") == none { 0 } else {
    calc.abs(_style-value(style, "pattern-amplitude"))
  }
)))

// How far an edge's paint reaches beside its carrier, in units of
// `unit-length`: half the stroke of each layer painted along the carrier, plus
// the amplitude of a coil, wave or zigzag. A layer without a stroke style
// takes the default edge stroke, 0.1em.
#let _painted-reach(layers, unit-length) = calc.max(0, ..layers.map(style => {
  if style == none or _style-value(style, "label-only") or _style-value(style, "offset") != 0 {
    return 0
  }
  let stroke = style.at("stroke", default: auto)
  let thickness = if stroke == none {
    return 0
  } else if type(stroke) == length {
    stroke
  } else if type(stroke) == dictionary {
    if stroke.at("paint", default: auto) == none { return 0 }
    stroke.at("thickness", default: auto)
  } else if type(stroke) == std.stroke {
    stroke.thickness
  } else { auto }
  if thickness == auto { thickness = 0.1em }
  let half = if type(thickness) == length {
    thickness.to-absolute() / unit-length.to-absolute() / 2
  } else { thickness / 2 }
  half + if _style-value(style, "pattern") == none { 0 } else {
    calc.abs(_style-value(style, "pattern-amplitude"))
  }
}))

// Carrier-following labels clear the painted band, and plain lines by at
// least their label gap.
#let _carrier-label-gap(layers, gap) = calc.max(
  gap, _painted-band(layers), ..layers.map(style => _style-value(style, "label-gap")),
)

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


#let _layer-label(style, data) = if style == none {
  none
} else {
  _as-content(_call(_style-value(style, "label"), data))
}

// An unlabelled bounded offset layer is a fixed decoration, such as an
// external momentum arrow whose text stays at the endpoint.
#let _fixed-arrow-layer(style, data) = (
  style != none and not _style-value(style, "label-only")
    and _layer-label(style, data) == none
    and _style-value(style, "offset") != 0
    and (_style-value(style, "length") != none or _style-value(style, "ratio") != none)
)
