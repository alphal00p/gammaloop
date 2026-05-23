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
  _data-field(edge, field, default)
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
  let arrow-half = needs-arrow and (
    (half == "source" and has-source)
      or (half == "sink" and has-sink)
  )
  if arrow-half {
    base + (
      mark: entry.at("fermion-arrow-mark", default: api.fermion-arrow-mark),
      mark-position: "center-if-dangling",
      mark-orientation: "edge",
    )
  } else {
    base
  }
}

#let _single-end-mark(mark, arrow-stroke) = {
  let base = (
    end: "straight",
    stroke: arrow-stroke,
    scale: .75,
  )
  if mark == auto {
    base
  } else if mark == none {
    none
  } else if type(mark) == dictionary {
    let clean = mark
    if clean.keys().contains("start") {
      let _ = clean.remove("start")
    }
    if clean.keys().contains("end") {
      base + clean
    } else {
      base + clean + (end: "barbed")
    }
  } else {
    base + (end: mark)
  }
}

#let _momentum-arrow-layer(options) = {
  let offset = options.offset
  let length = options.length
  let ratio = options.ratio
  let stroke = options.stroke
  let mark = options.mark
  let show-mark = options.show-mark
  let default-stroke = options.default-stroke
  let arrow-stroke = if stroke == none {
    default-stroke
  } else {
    stroke
  }
  let arrow-mark = _single-end-mark(mark, arrow-stroke)
  let layer = (
    offset: offset,
    length: length,
    ratio: ratio,
    resolve-length: "min",
    offset-side: "label",
    stroke: arrow-stroke,
  )
  if show-mark {
    if arrow-mark == none {
      layer
    } else {
      layer + (mark: arrow-mark)
    }
  } else {
    layer
  }
}

#let _has-source(edge) = {
  edge.at("source-half-edge", default: none) != none
}

#let _has-sink(edge) = {
  edge.at("sink-half-edge", default: none) != none
}

#let _momentum-arrow-half(edge) = {
  if _has-sink(edge) {
    "sink"
  } else if _has-source(edge) {
    "source"
  } else {
    none
  }
}

#let _half-style(edge, half, options) = {
  let api = options.api
  let map = options.map
  let default = options.default
  let typst-fields = options.typst-fields
  let scope = options.scope
  let orientation-split = options.orientation-split
  let momentum-arrows = options.momentum-arrows
  let momentum-arrow-offset = options.momentum-arrow-offset
  let momentum-arrow-length = options.momentum-arrow-length
  let momentum-arrow-ratio = options.momentum-arrow-ratio
  let momentum-arrow-stroke = options.momentum-arrow-stroke
  let momentum-arrow-mark = options.momentum-arrow-mark
  let style = _base-half-style(edge, half, (
    map: map,
    default: default,
    orientation-split: orientation-split,
    api: api,
  ))
  style = style + (api.style-dict)(edge.at(half + "-style", default: none), edge, mode: typst-fields, map: map, scope: scope)
  style = style + (api.style-dict)(_half-data-field(edge, half, "style", none), edge, mode: typst-fields, map: map, scope: scope)
  if momentum-arrows {
    let show-mark = _momentum-arrow-half(edge) == half
    (
      style,
      _momentum-arrow-layer((
        offset: momentum-arrow-offset,
        length: momentum-arrow-length,
        ratio: momentum-arrow-ratio,
        stroke: momentum-arrow-stroke,
        mark: momentum-arrow-mark,
        show-mark: show-mark,
        default-stroke: api.momentum-arrow-defaults.stroke,
      )),
    )
  } else {
    style
  }
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
  let show-momentum = options.show-momentum
  let show-edge-index = options.show-edge-index
  let show-half-edge-index = options.show-half-edge-index
  let show-particle = options.show-particle
  let momentum-fields = options.momentum-fields
  let edge-index-fields = options.edge-index-fields
  let momentum-prefix = options.momentum-prefix
  let edge-index-prefix = options.edge-index-prefix
  let half-edge-index-prefix = options.half-edge-index-prefix
  let particle-prefix = options.particle-prefix
  let label-separator = options.label-separator
  let label-size = options.label-size
  let label-fill = options.label-fill
  let pieces = ()
  if show-momentum {
    pieces.push(_prefixed-content(momentum-prefix, (api.momentum-value)(edge, fields: momentum-fields), api))
  }
  if show-edge-index {
    pieces.push(_prefixed-content(edge-index-prefix, (api.edge-index)(edge, fields: edge-index-fields), api))
  }
  if show-half-edge-index {
    pieces.push(_prefixed-content(half-edge-index-prefix, (api.dangling-half-edge-index)(edge), api))
  }
  if show-particle {
    pieces.push(_prefixed-content(particle-prefix, (api.particle-name)(edge), api))
  }
  let joined = _join-content(pieces, label-separator)
  if joined == none {
    none
  } else {
    text(size: label-size, fill: label-fill)[#joined]
  }
}

#let _default-label(edge, options) = {
  let api = options.api
  let map = options.map
  let default = options.default
  let typst-fields = options.typst-fields
  let scope = options.scope
  let data-label = _edge-data-field(
    edge,
    "display-label",
    _edge-data-field(edge, "label", none),
  )
  if data-label != none {
    (api.label-content)(data-label, edge, mode: typst-fields, map: map, scope: scope)
  } else {
    let label-value = edge.at(
      "display-label",
      default: edge.at("label", default: none),
    )
    if label-value != none {
      (api.label-content)(label-value, edge, mode: typst-fields, map: map, scope: scope)
    } else {
      (api.label-content)(
        (api.edge-entry)(edge, map: map, default: default).at("label", default: none),
        edge,
        map: map,
        scope: scope,
      )
    }
  }
}

#let edge-label(edge, options) = {
  let selected = _selected-label(
    edge,
    (
      show-momentum: options.show-momentum,
      show-edge-index: options.show-edge-index,
      show-half-edge-index: options.show-half-edge-index,
      show-particle: options.show-particle,
      momentum-fields: options.momentum-fields,
      edge-index-fields: options.edge-index-fields,
      momentum-prefix: options.momentum-prefix,
      edge-index-prefix: options.edge-index-prefix,
      half-edge-index-prefix: options.half-edge-index-prefix,
      particle-prefix: options.particle-prefix,
      label-separator: options.label-separator,
      label-size: options.label-size,
      label-fill: options.label-fill,
      api: options.api,
    ),
  )
  if selected != none {
    selected
  } else {
    _default-label(edge, (
      map: options.map,
      default: options.default,
      typst-fields: options.typst-fields,
      scope: options.scope,
      api: options.api,
    ))
  }
}
