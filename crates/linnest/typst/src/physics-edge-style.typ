#import "@preview/mitex:0.2.6" as mitex

#let mi = mitex.mi

/// Conventional massive-particle stroke thickness.
/// -> length
#let massive = 1.0pt

/// Conventional massless-particle stroke thickness.
/// -> length
#let massless = 0.55pt

/// Dash pattern used for scalar-style lines.
/// -> array
#let dashed = (0.1em, 0.45em)

/// Dotted dash style used for ghost-style lines.
/// -> string
#let dotted = "dotted"

/// Build a rounded CeTZ stroke dictionary accepted by `linnest.draw`.
/// -> dictionary
#let stroke-style(c: black, thickness: massless, dash: none) = {
  let stroke = (paint: c, thickness: thickness, cap: "round")
  if dash == none {
    (stroke: stroke)
  } else {
    (stroke: stroke + (dash: dash))
  }
}

/// Source-half stroke helper.
/// -> dictionary
#let source-stroke(c: black, thickness: massless, dash: none) = {
  stroke-style(c: c, thickness: thickness, dash: dash)
}

/// Sink-half stroke helper. The default lightening makes the two halves encode
/// the graph's source/sink split without requiring arrowheads.
/// -> dictionary
#let sink-stroke(c: black, thickness: massless, dash: none, lighten: 45%) = {
  stroke-style(c: c.lighten(lighten), thickness: thickness, dash: dash)
}

/// Photon-style wave pattern.
/// -> dictionary
#let wave = (
  pattern: "wave",
  pattern-amplitude: 0.14,
  pattern-wavelength: 0.55,
)

/// Gluon-style coil pattern.
/// -> dictionary
#let coil = (
  pattern: "coil",
  pattern-amplitude: 0.14,
  pattern-wavelength: 0.55,
  pattern-coil-longitudinal-scale: 1.6,
)

/// Weak-boson-style zigzag pattern.
/// -> dictionary
#let zigzag = (
  pattern: "zigzag",
  pattern-amplitude: 0.14,
  pattern-wavelength: 0.55,
)

/// Fallback style entry used when an edge has no known particle type.
/// -> dictionary
#let default-edge = (
  source: source-stroke(),
  sink: sink-stroke(),
  label: none,
)

/// Small built-in particle map for standalone physics diagrams. GammaLoop
/// generated `edge-style.typ` files pass their model-specific map to `style`.
/// -> dictionary
#let default-map = (
  "a": (
    source: source-stroke(c: black, thickness: massless) + wave,
    sink: sink-stroke(c: black, thickness: massless) + wave,
    label: mi(`{\gamma}`),
  ),
  "photon": (
    source: source-stroke(c: black, thickness: massless) + wave,
    sink: sink-stroke(c: black, thickness: massless) + wave,
    label: mi(`{\gamma}`),
  ),
  "g": (
    source: source-stroke(c: black, thickness: massless) + coil,
    sink: sink-stroke(c: black, thickness: massless) + coil,
    label: mi(`{g}`),
  ),
  "gluon": (
    source: source-stroke(c: black, thickness: massless) + coil,
    sink: sink-stroke(c: black, thickness: massless) + coil,
    label: mi(`{g}`),
  ),
  "fermion": (
    source: source-stroke(c: blue, thickness: massless),
    sink: sink-stroke(c: blue, thickness: massless),
    label: mi(`{f}`),
  ),
  "scalar": (
    source: source-stroke(c: black, thickness: massive, dash: dashed),
    sink: sink-stroke(c: black, thickness: massive, dash: dashed),
    label: mi(`{\phi}`),
  ),
)

/// Return an edge's particle name, stripping the quotes often present in DOT
/// statement values.
/// -> none | string
#let particle-name(edge) = {
  let particle = edge.at("particle", default: none)
  if particle == none {
    none
  } else {
    str(particle).trim("\"")
  }
}

/// Return the style-map entry for an edge.
/// -> dictionary
#let edge-entry(edge, map: default-map, default: default-edge) = {
  let particle = particle-name(edge)
  if particle == none {
    default
  } else {
    map.at(particle, default: default)
  }
}

/// Convert DOT-ish values to plain text.
/// -> string
#let text-value(value) = str(value).trim("\"")

/// Replace `{field}` placeholders in a string-valued style template.
/// `{{` and `}}` produce literal braces.
/// -> none | string
#let interpolate-template(template, edge) = {
  if template == none {
    none
  } else {
    let text = text-value(template)
    text = text.replace("{{", "\u{e000}").replace("}}", "\u{e001}")
    for key in edge.keys() {
      let value = edge.at(key)
      if value != none {
        text = text.replace("{" + str(key) + "}", text-value(value))
      }
    }
    text.replace("\u{e000}", "{").replace("\u{e001}", "}")
  }
}

#let eval-scope(map: default-map) = (
  mi: mi,
  massive: massive,
  massless: massless,
  dashed: dashed,
  dotted: dotted,
  stroke-style: stroke-style,
  source-stroke: source-stroke,
  sink-stroke: sink-stroke,
  wave: wave,
  coil: coil,
  zigzag: zigzag,
  default-edge: default-edge,
  default-map: default-map,
  map: map,
)

/// Evaluate a string-valued style template after placeholder interpolation.
/// -> any
#let eval-template(template, edge, map: default-map, scope: (:)) = {
  eval(interpolate-template(template, edge), scope: eval-scope(map: map) + scope + edge)
}

#let _assert-typst-fields-mode(mode) = {
  if mode != "plain" and mode != "eval" {
    panic("typst-fields must be \"plain\" or \"eval\"")
  }
}

/// Convert a style label value to content. In `mode: "eval"`, string values are
/// evaluated after placeholder interpolation.
/// -> none | content
#let label-content(value, edge, mode: "plain", map: default-map, scope: (:)) = {
  _assert-typst-fields-mode(mode)
  if value == none {
    none
  } else if type(value) == content {
    value
  } else {
    let text = interpolate-template(value, edge)
    if mode == "eval" {
      eval(text, scope: eval-scope(map: map) + scope + edge)
    } else {
      [#text]
    }
  }
}

/// Convert a style value to a dictionary. In `mode: "eval"`, string values are
/// evaluated after placeholder interpolation.
/// -> dictionary
#let style-dict(value, edge, mode: "plain", map: default-map, scope: (:)) = {
  _assert-typst-fields-mode(mode)
  if value == none {
    (:)
  } else if type(value) == dictionary {
    value
  } else if mode == "eval" {
    eval-template(value, edge, map: map, scope: scope)
  } else {
    (:)
  }
}

#let _base-half-style(edge, half, map: default-map, default: default-edge, orientation-split: true) = {
  let entry = edge-entry(edge, map: map, default: default)
  let source = entry.at("source", default: (:))
  let sink = entry.at("sink", default: source)
  if half == "source" {
    source
  } else if orientation-split {
    sink
  } else {
    source
  }
}

#let _momentum-arrow-layer(
  offset: 0.16,
  length: 1.0,
  ratio: 0.5,
  stroke: none,
  mark: none,
  show-mark: false,
) = {
  let arrow-stroke = if stroke == none {
    (paint: black, thickness: 0.55pt, cap: "round")
  } else {
    stroke
  }
  let arrow-mark = if mark == none {
    (
      end: "barbed",
      stroke: arrow-stroke,
      fill: black,
      scale: 1.35,
    )
  } else if type(mark) == dictionary {
    (
      (
        stroke: arrow-stroke,
        fill: black,
      )
        + mark
    )
  } else {
    (
      end: mark,
      stroke: arrow-stroke,
      fill: black,
    )
  }
  let layer = (
    offset: offset,
    length: length,
    ratio: ratio,
    resolve-length: "min",
    offset-side: "label",
    stroke: arrow-stroke,
  )
  if show-mark {
    layer + (mark: arrow-mark)
  } else {
    layer
  }
}

#let _half-style(
  edge,
  half,
  map: default-map,
  default: default-edge,
  typst-fields: "plain",
  scope: (:),
  orientation-split: true,
  momentum-arrows: false,
  momentum-arrow-offset: 0.16,
  momentum-arrow-length: 1.0,
  momentum-arrow-ratio: 0.5,
  momentum-arrow-stroke: none,
  momentum-arrow-mark: none,
) = {
  let style = _base-half-style(edge, half, map: map, default: default, orientation-split: orientation-split)
  style = style + style-dict(edge.at(half + "-style", default: none), edge, mode: typst-fields, map: map, scope: scope)
  style = style + style-dict(edge.at(half + "-style-eval", default: none), edge, mode: "eval", map: map, scope: scope)
  if momentum-arrows {
    let full-edge = edge.at("source", default: none) != none and edge.at("sink", default: none) != none
    let show-mark = if full-edge { half == "sink" } else { true }
    (
      style,
      _momentum-arrow-layer(
        offset: momentum-arrow-offset,
        length: momentum-arrow-length,
        ratio: momentum-arrow-ratio,
        stroke: momentum-arrow-stroke,
        mark: momentum-arrow-mark,
        show-mark: show-mark,
      ),
    )
  } else {
    style
  }
}

/// Style callback for the source half edge.
///
/// `momentum-arrows: true` adds one centered black CeTZ-mark decoration while
/// the main edge remains drawn with its normal particle style.
/// -> dictionary | array
#let source-style(
  edge,
  map: default-map,
  default: default-edge,
  typst-fields: "plain",
  scope: (:),
  orientation-split: true,
  momentum-arrows: false,
  momentum-arrow-offset: 0.16,
  momentum-arrow-length: 1.0,
  momentum-arrow-ratio: 0.5,
  momentum-arrow-stroke: none,
  momentum-arrow-mark: none,
) = _half-style(
  edge,
  "source",
  map: map,
  default: default,
  typst-fields: typst-fields,
  scope: scope,
  orientation-split: orientation-split,
  momentum-arrows: momentum-arrows,
  momentum-arrow-offset: momentum-arrow-offset,
  momentum-arrow-length: momentum-arrow-length,
  momentum-arrow-ratio: momentum-arrow-ratio,
  momentum-arrow-stroke: momentum-arrow-stroke,
  momentum-arrow-mark: momentum-arrow-mark,
)

/// Style callback for the sink half edge.
///
/// `momentum-arrows: true` adds one centered black CeTZ-mark decoration toward
/// the sink node, so momentum arrows always flow from source to sink
/// independently of `edge.orientation`.
/// -> dictionary | array
#let sink-style(
  edge,
  map: default-map,
  default: default-edge,
  typst-fields: "plain",
  scope: (:),
  orientation-split: true,
  momentum-arrows: false,
  momentum-arrow-offset: 0.16,
  momentum-arrow-length: 1.0,
  momentum-arrow-ratio: 0.5,
  momentum-arrow-stroke: none,
  momentum-arrow-mark: none,
) = _half-style(
  edge,
  "sink",
  map: map,
  default: default,
  typst-fields: typst-fields,
  scope: scope,
  orientation-split: orientation-split,
  momentum-arrows: momentum-arrows,
  momentum-arrow-offset: momentum-arrow-offset,
  momentum-arrow-length: momentum-arrow-length,
  momentum-arrow-ratio: momentum-arrow-ratio,
  momentum-arrow-stroke: momentum-arrow-stroke,
  momentum-arrow-mark: momentum-arrow-mark,
)

#let _field-value(edge, fields) = {
  if fields == none {
    none
  } else if type(fields) == array {
    let value = none
    for field in fields {
      if value == none {
        value = edge.at(field, default: none)
      }
    }
    value
  } else {
    edge.at(fields, default: none)
  }
}

/// Return the momentum field used by optional edge labels.
/// -> none | any
#let momentum-value(edge, fields: ("momentum", "mom", "q")) = _field-value(edge, fields)

/// Return the edge index used by optional edge labels. The DOT `id` statement
/// wins over the renderer-local `eid`.
/// -> none | any
#let edge-index(edge, fields: ("id", "eid")) = _field-value(edge, fields)

/// Return the exposed half-edge index for dangling edges only.
/// -> none | any
#let dangling-half-edge-index(edge) = {
  if not edge.at("ext", default: false) {
    none
  } else {
    let endpoint = if edge.at("source", default: none) == none {
      edge.at("sink-endpoint", default: none)
    } else if edge.at("sink", default: none) == none {
      edge.at("source-endpoint", default: none)
    } else {
      none
    }
    if endpoint == none {
      none
    } else {
      endpoint.at("statement", default: endpoint.at("id", default: endpoint.at("hedge", default: none)))
    }
  }
}

#let _prefixed-content(prefix, value) = {
  if value == none {
    none
  } else if prefix == none {
    [#text-value(value)]
  } else {
    prefix + [#text-value(value)]
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

#let _selected-label(
  edge,
  show-momentum: false,
  show-edge-index: false,
  show-half-edge-index: false,
  show-particle: false,
  momentum-fields: ("momentum", "mom", "q"),
  edge-index-fields: ("id", "eid"),
  momentum-prefix: none,
  edge-index-prefix: [e],
  half-edge-index-prefix: [h],
  particle-prefix: none,
  label-separator: [, ],
  label-size: 7pt,
  label-fill: red,
) = {
  let pieces = ()
  if show-momentum {
    pieces.push(_prefixed-content(momentum-prefix, momentum-value(edge, fields: momentum-fields)))
  }
  if show-edge-index {
    pieces.push(_prefixed-content(edge-index-prefix, edge-index(edge, fields: edge-index-fields)))
  }
  if show-half-edge-index {
    pieces.push(_prefixed-content(half-edge-index-prefix, dangling-half-edge-index(edge)))
  }
  if show-particle {
    pieces.push(_prefixed-content(particle-prefix, particle-name(edge)))
  }
  let joined = _join-content(pieces, label-separator)
  if joined == none {
    none
  } else {
    text(size: label-size, fill: label-fill)[#joined]
  }
}

#let _default-label(edge, map: default-map, default: default-edge, typst-fields: "plain", scope: (:)) = {
  let label-eval = edge.at("label-eval", default: none)
  if label-eval != none {
    label-content(label-eval, edge, mode: "eval", map: map, scope: scope)
  } else {
    let label-template = edge.at(
      "display-label",
      default: edge.at("label-template", default: edge.at("label", default: none)),
    )
    if label-template != none {
      label-content(label-template, edge, mode: typst-fields, map: map, scope: scope)
    } else {
      label-content(
        edge-entry(edge, map: map, default: default).at("label", default: none),
        edge,
        map: map,
        scope: scope,
      )
    }
  }
}

/// Edge-label callback.
///
/// By default this preserves explicit `label-eval`, `display-label`,
/// `label-template`, `label`, and particle-map labels. Set any `show-*` option
/// to build a label from selected metadata instead:
///
/// ```example
/// #let callbacks = physics.style(
///   momentum-arrows: true,
///   show-edge-index: true,
///   show-particle: true,
/// )
/// ```
/// -> none | content
#let edge-label(
  edge,
  map: default-map,
  default: default-edge,
  typst-fields: "plain",
  scope: (:),
  show-momentum: false,
  show-edge-index: false,
  show-half-edge-index: false,
  show-particle: false,
  momentum-fields: ("momentum", "mom", "q"),
  edge-index-fields: ("id", "eid"),
  momentum-prefix: none,
  edge-index-prefix: [e],
  half-edge-index-prefix: [h],
  particle-prefix: none,
  label-separator: [, ],
  label-size: 7pt,
  label-fill: red,
) = {
  let selected = _selected-label(
    edge,
    show-momentum: show-momentum,
    show-edge-index: show-edge-index,
    show-half-edge-index: show-half-edge-index,
    show-particle: show-particle,
    momentum-fields: momentum-fields,
    edge-index-fields: edge-index-fields,
    momentum-prefix: momentum-prefix,
    edge-index-prefix: edge-index-prefix,
    half-edge-index-prefix: half-edge-index-prefix,
    particle-prefix: particle-prefix,
    label-separator: label-separator,
    label-size: label-size,
    label-fill: label-fill,
  )
  if selected != none {
    selected
  } else {
    _default-label(edge, map: map, default: default, typst-fields: typst-fields, scope: scope)
  }
}

/// Bundle source-style, sink-style, and edge-label callbacks with shared
/// options for `linnest.draw`.
///
/// ```example
/// #let b = graph.builder(name: "physics")
/// #let (node: a, builder: b) = graph.node(b, name: "a")
/// #let (node: b-node, builder: b) = graph.node(b, name: "b")
/// #let b = graph.edge(
///   b,
///   source: (node: a),
///   sink: (node: b-node),
///   statements: (particle: "g", id: "7"),
/// )
/// #let callbacks = physics.style(momentum-arrows: true, show-edge-index: true)
/// #draw(
///   layout(graph.finish(b)),
///   source-style: callbacks.source-style,
///   sink-style: callbacks.sink-style,
///   edge-label: callbacks.edge-label,
/// )
/// ```
/// -> dictionary
#let style(
  map: default-map,
  default: default-edge,
  typst-fields: "plain",
  scope: (:),
  orientation-split: true,
  momentum-arrows: false,
  momentum-arrow-offset: 0.16,
  momentum-arrow-length: 1.0,
  momentum-arrow-ratio: 0.5,
  momentum-arrow-stroke: none,
  momentum-arrow-mark: none,
  show-momentum: false,
  show-edge-index: false,
  show-half-edge-index: false,
  show-particle: false,
  momentum-fields: ("momentum", "mom", "q"),
  edge-index-fields: ("id", "eid"),
  momentum-prefix: none,
  edge-index-prefix: [e],
  half-edge-index-prefix: [h],
  particle-prefix: none,
  label-separator: [, ],
  label-size: 7pt,
  label-fill: red,
) = (
  source-style: edge => source-style(
    edge,
    map: map,
    default: default,
    typst-fields: typst-fields,
    scope: scope,
    orientation-split: orientation-split,
    momentum-arrows: momentum-arrows,
    momentum-arrow-offset: momentum-arrow-offset,
    momentum-arrow-length: momentum-arrow-length,
    momentum-arrow-ratio: momentum-arrow-ratio,
    momentum-arrow-stroke: momentum-arrow-stroke,
    momentum-arrow-mark: momentum-arrow-mark,
  ),
  sink-style: edge => sink-style(
    edge,
    map: map,
    default: default,
    typst-fields: typst-fields,
    scope: scope,
    orientation-split: orientation-split,
    momentum-arrows: momentum-arrows,
    momentum-arrow-offset: momentum-arrow-offset,
    momentum-arrow-length: momentum-arrow-length,
    momentum-arrow-ratio: momentum-arrow-ratio,
    momentum-arrow-stroke: momentum-arrow-stroke,
    momentum-arrow-mark: momentum-arrow-mark,
  ),
  edge-label: edge => edge-label(
    edge,
    map: map,
    default: default,
    typst-fields: typst-fields,
    scope: scope,
    show-momentum: show-momentum,
    show-edge-index: show-edge-index,
    show-half-edge-index: show-half-edge-index,
    show-particle: show-particle,
    momentum-fields: momentum-fields,
    edge-index-fields: edge-index-fields,
    momentum-prefix: momentum-prefix,
    edge-index-prefix: edge-index-prefix,
    half-edge-index-prefix: half-edge-index-prefix,
    particle-prefix: particle-prefix,
    label-separator: label-separator,
    label-size: label-size,
    label-fill: label-fill,
  ),
)
