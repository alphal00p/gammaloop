// Internal physics edge-style implementation. Public users should import
// `../physics-edge-style.typ`.

#let _record-data(record) = {
  if type(record) != dictionary {
    (:)
  } else {
    let data = record.at("data", default: none)
    if type(data) == dictionary { data } else { (:) }
  }
}

#let _data-field(record, field, default) = {
  _record-data(record).at(field, default: default)
}

#let _edge-data-field(edge, field, default) = {
  edge.at(field, default: _data-field(edge, field, default))
}

#let _half-record(edge, half) = {
  let half-edge = edge.at(half + "-half-edge", default: none)
  if half-edge != none {
    half-edge
  } else {
    let value = edge.at(half, default: none)
    if type(value) == dictionary { value } else { none }
  }
}

#let _half-data-field(edge, half, field, default) = {
  _data-field(_half-record(edge, half), field, default)
}

#let _base-half-style(edge, half, options) = {
  let api = options.api
  let map = options.map
  let default = options.default
  let orientation-split = options.orientation-split
  let entry = (api.edge-entry)(edge, map: map, default: default)
  let source = entry.at("source", default: (:))
  let sink = entry.at("sink", default: source)
  let base = if half == "source" {
    source
  } else if orientation-split {
    sink
  } else {
    source
  }
  let has-source = edge.at("source-half-edge", default: none) != none
  let has-sink = edge.at("sink-half-edge", default: none) != none
  let needs-arrow = entry.at("fermion-arrow", default: false)
  let arrow-half = (
    needs-arrow
      and (
        (half == "source" and has-source) or (half == "sink" and has-sink)
      )
  )
  if arrow-half {
    (
      base
        + (
          mark: entry.at("fermion-arrow-mark", default: api.fermion-arrow-mark),
          mark-position: "center-if-dangling",
          mark-orientation: "edge",
        )
    )
  } else {
    base
  }
}

#let _single-end-mark(mark, arrow-stroke) = {
  if type(mark) == str {
    mark = mark.trim("\"")
    if mark == "none" { mark = none } else if mark == "auto" { mark = auto }
  }
  let base = (end: "straight", stroke: arrow-stroke, scale: 1.1)
  if mark == auto { base } else if mark == none { none } else if (
    type(mark) == dictionary
  ) {
    let clean = mark
    if clean.keys().contains("start") { let _ = clean.remove("start") }
    base + clean
  } else { base + (end: mark) }
}

#let _momentum-arrow-layer(edge, options) = {
  let side = options.momentum-arrow-side
  let side = if side == auto { "auto" } else { str(side).trim("\"") }
  assert(
    side in ("auto", "left", "right"),
    message: "momentum-arrow-side must be auto, left, or right",
  )
  let stroke = options.momentum-arrow-stroke
  if stroke == none { stroke = options.api.momentum-arrow-defaults.stroke }
  let offset = options.momentum-arrow-offset
  if side == "auto" {
    // Use the edge's bend relative to its chord, independently of label content
    // or relaxation. Straight, dangling and closed edges use the signed offset.
    let direction = if offset < 0 { -1 } else { 1 }
    let source = edge.at("source-node", default: none)
    let sink = edge.at("sink-node", default: none)
    let record = edge.at("edge", default: none)
    let point = if type(record) == dictionary {
      record.at("pos", default: none)
    } else { none }
    if type(point) == array { point = (x: point.at(0), y: point.at(1)) }
    if source != none and sink != none and point != none {
      let a = source.pos
      let b = sink.pos
      let cross = (b.x - a.x) * (point.y - a.y) - (b.y - a.y) * (point.x - a.x)
      if calc.abs(cross) > 1e-9 {
        direction *= if cross < 0 { -1 } else { 1 }
      }
    }
    side = if direction < 0 { "right" } else { "left" }
  }
  let anchor = options.momentum-label-anchor
  // Resolved sides are relative to source -> sink and shared with the label carrier.
  // The gap is in graph units: to the box for auto, or to an explicit anchor.
  let geometry = (
    offset: calc.abs(offset) * if side == "left" { 1 } else { -1 },
    offset-side: none,
    label-side: side,
    label-gap: options.momentum-label-gap,
    label-style: (
      anchor: if type(anchor) == str { anchor.trim("\"") } else { anchor },
    ),
  )
  let shift = options.momentum-arrow-shift
  let label-shift = options.momentum-label-shift
  // Auto clears the complete label box; explicit anchors use only the gap.
  // Labels follow the full invisible path, so endpoint clamps never depend on
  // arrow length. Linnest resolves label:auto and owns all generic overlays.
  (
    geometry
      + (
        length: options.momentum-arrow-length,
        ratio: options.momentum-arrow-ratio,
        resolve-length: if options.momentum-arrow-ratio == none {
          "length"
        } else { "min" },
        shift: shift,
        stroke: stroke,
        pattern: none,
        mark: if options.show-mark {
          _single-end-mark(options.momentum-arrow-mark, stroke)
        } else { none },
        // A numeric end position uses Linnest's continuous paired-edge carrier,
        // keeping one arrowhead when a shift crosses the source/sink split.
        mark-position: 1,
        mark-orientation: "path",
      ),
    geometry
      + (
        length: none,
        ratio: none,
        resolve-length: "none",
        shift: 0,
        label-only: true,
        label: auto,
        label-shift: if label-shift == auto { shift } else { label-shift },
      ),
  )
}

#let _half-style(edge, half, options) = {
  let style = _base-half-style(edge, half, options)
  let enabled = _edge-data-field(
    edge,
    "momentum-arrows",
    options.momentum-arrows,
  )
  if enabled not in (true, "true", "\"true\"") { return style }
  // Domain controls create particle, arrow, and label layers directly; generic
  // styles and decorations are composed once by Linnest after these defaults.
  for key in (
    "momentum-arrow-offset",
    "momentum-arrow-length",
    "momentum-arrow-ratio",
    "momentum-arrow-stroke",
    "momentum-arrow-mark",
    "momentum-arrow-side",
    "momentum-arrow-shift",
    "momentum-label-gap",
    "momentum-label-shift",
    "momentum-label-anchor",
  ) {
    let value = _half-data-field(edge, half, key, _edge-data-field(
      edge,
      key,
      options.at(key),
    ))
    if (
      key
        in (
          "momentum-arrow-offset",
          "momentum-arrow-length",
          "momentum-arrow-ratio",
          "momentum-arrow-shift",
          "momentum-label-gap",
          "momentum-label-shift",
        )
        and type(value) == str
    ) {
      let value-text = value.trim("\"")
      value = if value-text == "none" { none } else if value-text == "auto" {
        auto
      } else { float(value-text) }
    }
    options.insert(key, value)
  }
  let arrow-half = if edge.at("sink-half-edge", default: none) != none {
    "sink"
  } else { "source" }
  let layers = _momentum-arrow-layer(edge, options + (show-mark: half == arrow-half))
  let anchor = _edge-data-field(edge, "label-anchor", none)
  let dangling = (
    edge.at("source-half-edge", default: none) == none
      or edge.at("sink-half-edge", default: none) == none
  )
  // Ordered external labels stay outside their endpoint at Linnest's solved
  // label position. Explicit momentum placement still selects the path carrier.
  if (
    dangling and type(anchor) == str and anchor.trim("\"") in ("east", "west")
      and options.momentum-arrow-shift == 0
      and options.momentum-label-shift == auto
      and options.momentum-label-anchor == auto
  ) {
    layers = layers.slice(0, 1)
  }
  (style, ..layers)
}

#let source-style(edge, options) = _half-style(edge, "source", options)

#let sink-style(edge, options) = _half-style(edge, "sink", options)

#let _prefixed-content(prefix, value, api) = {
  if value == none {
    none
  } else if prefix == none {
    [#((api.text-value)(value))]
  } else {
    prefix + [#((api.text-value)(value))]
  }
}

#let _join-content(pieces, separator) = {
  let out = none
  for piece in pieces {
    if piece != none {
      if out == none {
        out = piece
      } else {
        out = out + separator + piece
      }
    }
  }
  out
}

#let _selected-label(edge, options) = {
  let api = options.api
  let pieces = ()
  if options.show-particle != false {
    let particle = (api.edge-entry)(
      edge,
      map: options.map,
      default: options.default,
    ).at("label", default: none)
    if particle != none {
      particle = (api.label-content)(
        particle,
        edge,
        map: options.map,
        scope: options.scope,
      )
      pieces.push(if options.particle-prefix == none { particle } else {
        options.particle-prefix + particle
      })
    }
  }
  if options.show-momentum {
    pieces.push([$q_(#edge.eid)$])
  }
  if options.show-edge-index {
    pieces.push(_prefixed-content(
      options.edge-index-prefix,
      (api.edge-index)(edge, fields: options.edge-index-fields),
      api,
    ))
  }
  if options.show-half-edge-index {
    pieces.push(_prefixed-content(
      options.half-edge-index-prefix,
      (api.dangling-half-edge-index)(edge),
      api,
    ))
  }
  let joined = _join-content(pieces, options.label-separator)
  if joined == none { none } else {
    text(size: options.label-size, fill: options.label-fill)[#joined]
  }
}

#let edge-label(edge, options) = {
  let api = options.api
  // Explicit labels take precedence over generated particle/momentum portions.
  // Explicit false show flags suppress their portion, including all labels,
  // rather than falling through to the particle-map default.
  for (index, record) in (_record-data(edge), edge).enumerate() {
    for key in ("display-label", "label") {
      if (
        record.keys().contains(key) and (index == 0 or record.at(key) != none)
      ) {
        return (api.label-content)(
          record.at(key),
          edge,
          mode: options.typst-fields,
          map: options.map,
          scope: options.scope,
        )
      }
    }
  }
  _selected-label(edge, options)
}
