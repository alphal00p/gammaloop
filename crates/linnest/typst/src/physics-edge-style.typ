// Public physics edge-style API.
//
// Keep documented defaults here. The implementation lives in
// `impl/physics-edge-style.typ` and takes explicit values/options.

#import "@preview/mitex:0.2.6" as mitex
#import "impl/physics-edge-style.typ" as _impl

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

/// CeTZ marker used for fermion particle-flow arrows on the main edge.
/// -> dictionary
#let fermion-arrow-mark = (
  end: (
    symbol: ">",
    fill: black,
    stroke: black + 0.2pt,
    anchor: "center",
    shorten-to: auto,
  ),
  scale: 0.75,
)

/// Mark an edge-map entry as a fermion so the main edge receives one
/// particle-flow arrow that follows `edge.orientation`. This is separate from
/// optional momentum arrows.
/// -> dictionary
#let fermion-flow = (
  fermion-arrow: true,
)

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
    fermion-arrow: true,
    label: mi(`{f}`),
  ),
  "scalar": (
    source: source-stroke(c: black, thickness: massive, dash: dashed),
    sink: sink-stroke(c: black, thickness: massive, dash: dashed),
    label: mi(`{\phi}`),
  ),
)

#let _momentum-arrow-stroke = (paint: black, thickness: 0.55pt, cap: "round")

/// Default style for source-to-sink momentum arrow layers.
/// -> dictionary
#let momentum-arrow-defaults = (
  offset: 0.46,
  length: 5.0,
  ratio: 0.5,
  stroke: _momentum-arrow-stroke,
  mark: auto,
)

/// Return an edge's particle name, stripping the quotes often present in DOT
/// statement values.
/// -> none | string
#let particle-name(edge) = {
  let particle = edge.at("particle", default: none)
  if particle == none {
    let data = edge.at("data", default: none)
    if type(data) == dictionary {
      particle = data.at("particle", default: none)
    }
  }
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

#let eval-scope(map: default-map) = (
  mi: mi,
  massive: massive,
  massless: massless,
  dashed: dashed,
  dotted: dotted,
  stroke-style: stroke-style,
  source-stroke: source-stroke,
  sink-stroke: sink-stroke,
  fermion-arrow-mark: fermion-arrow-mark,
  fermion-flow: fermion-flow,
  wave: wave,
  coil: coil,
  zigzag: zigzag,
  default-edge: default-edge,
  default-map: default-map,
  map: map,
)

/// Evaluate a string-valued style expression in the edge scope.
/// -> any
#let eval-expression(value, edge, map: default-map, scope: (:)) = {
  eval(text-value(value), scope: eval-scope(map: map) + scope + edge)
}

#let _assert-typst-fields-mode(mode) = {
  if mode != "plain" and mode != "eval" {
    panic("typst-fields must be \"plain\" or \"eval\"")
  }
}

/// Convert a style label value to content. In `mode: "eval"`, string values are
/// evaluated in the edge scope.
/// -> none | content
#let label-content(value, edge, mode: "plain", map: default-map, scope: (:)) = {
  _assert-typst-fields-mode(mode)
  if value == none {
    none
  } else if type(value) == content {
    value
  } else if mode == "eval" {
    eval-expression(value, edge, map: map, scope: scope)
  } else {
    [#text-value(value)]
  }
}

/// Convert a style value to a dictionary. In `mode: "eval"`, string values are
/// evaluated in the edge scope.
/// -> dictionary
#let style-dict(value, edge, mode: "plain", map: default-map, scope: (:)) = {
  _assert-typst-fields-mode(mode)
  if value == none {
    (:)
  } else if type(value) == dictionary {
    value
  } else if mode == "eval" {
    eval-expression(value, edge, map: map, scope: scope)
  } else {
    (:)
  }
}

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

#let _has-source(edge) = {
  edge.at("source-half-edge", default: none) != none
}

#let _has-sink(edge) = {
  edge.at("sink-half-edge", default: none) != none
}

/// Return the exposed half-edge index for dangling edges only.
/// -> none | any
#let dangling-half-edge-index(edge) = {
  if not edge.at("ext", default: false) {
    none
  } else {
    let half-edge = if not _has-source(edge) {
      edge.at("sink-half-edge", default: none)
    } else if not _has-sink(edge) {
      edge.at("source-half-edge", default: none)
    } else {
      none
    }
    if half-edge == none {
      none
    } else {
      half-edge.at("statement", default: half-edge.at("id", default: half-edge.at("hedge", default: none)))
    }
  }
}

#let _api() = eval-scope() + (
  momentum-arrow-defaults: momentum-arrow-defaults,
  particle-name: particle-name,
  edge-entry: edge-entry,
  text-value: text-value,
  label-content: label-content,
  style-dict: style-dict,
  momentum-value: momentum-value,
  edge-index: edge-index,
  dangling-half-edge-index: dangling-half-edge-index,
)

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
  momentum-arrow-offset: momentum-arrow-defaults.offset,
  momentum-arrow-length: momentum-arrow-defaults.length,
  momentum-arrow-ratio: momentum-arrow-defaults.ratio,
  momentum-arrow-stroke: none,
  momentum-arrow-mark: momentum-arrow-defaults.mark,
) = _impl.source-style(edge, (
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
  api: _api(),
))

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
  momentum-arrow-offset: momentum-arrow-defaults.offset,
  momentum-arrow-length: momentum-arrow-defaults.length,
  momentum-arrow-ratio: momentum-arrow-defaults.ratio,
  momentum-arrow-stroke: none,
  momentum-arrow-mark: momentum-arrow-defaults.mark,
) = _impl.sink-style(edge, (
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
  api: _api(),
))

/// Edge-label callback.
///
/// By default this preserves data `display-label` or `label`, then explicit
/// `display-label`, `label`, and particle-map labels. Set any
/// `show-*` option to build a label from selected metadata instead:
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
) = _impl.edge-label(edge, (
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
  api: _api(),
))

/// Bundle source-style, sink-style, and edge-label callbacks with shared
/// options for `linnest.draw`.
///
/// ````example
/// #let g = build({
///   node(<a>,pos:pos(y:0))
///   node(<b>,pos:pos(y:0))
///   edge(source(<a>), <g-edge>, sink(<b>), particle: "g")
///   edge(source(<a>), <g-edge2>, sink(<b>), particle: "g")
///   edge(source(<a>), <g-edge3>, sink(<b>), particle: "g")
/// },
///   name: "physics",
/// )
/// #let a = ```dot
/// digraph {
/// 0 [dod="-100" num="1"];
///	1 [dod="-100" num="1"];
///	2 [dod="-100" num="1"];
///	3 [dod="-100" num="1"];
///	4 [dod="-100" num="1"];
///	5 [dod="-100" num="1"];
///	6 [dod="-100" num="1"];
///	exte0	 [style=invis];
///	exte0	-> 2:0	 [id=0  dod="-100" lmb_rep="P(0,a___)" num="1" particle="d" pin="x:@-left"];
///	exte1	 [style=invis];
///	exte1	-> 1:1	 [id=1 dir=back  dod="-100" lmb_rep="P(1,a___)" num="1" particle="d~" pin="x:@-left"];
///	exte2	 [style=invis];
///	exte2	-> 1:2	 [id=2 dir=none  dod="-100" lmb_rep="P(2,a___)" mass="0" num="1" particle="H" pin="x:@-left"];
///	exte3	 [style=invis];
///	exte3	-> 2:3	 [id=3 dir=none  dod="-100" is_dummy="true" lmb_rep="0" mass="0" num="1" particle="H" pin="x:@-left"];
///	exte4	 [style=invis];
///	6:4	-> exte4	 [id=4 dir=none  dod="-100" lmb_rep="P(3,a___)" mass="0" num="1" particle="H" pin="x:@+right"];
///	exte5	 [style=invis];
///	1:5	-> exte5	 [id=5 dir=none  dod="-100" is_dummy="true" lmb_rep="0" mass="0" num="1" particle="H" pin="x:@+right"];
///	exte6	 [style=invis];
///	5:6	-> exte6	 [id=6 dir=none  dod="-100" lmb_rep="P(4,a___)" num="1" particle="H" pin="x:@+right"];
///	exte7	 [style=invis];
///	4:7	-> exte7	 [id=7 dir=none  dod="-100" lmb_rep="P(5,a___)" num="1" particle="H" pin="x:@+right"];
///	exte8	 [style=invis];
///	3:8	-> exte8	 [id=8 dir=none  dod="-100" lmb_rep="-1*P(3,a___)+-1*P(4,a___)+-1*P(5,a___)+P(0,a___)+P(1,a___)+P(2,a___)" num="1" particle="H" pin="x:@+right"];
///	0:9	-> 1:10	 [id=9  dod="-100" lmb_id="1" lmb_rep="K(1,a___)" num="1" particle="t"];
///	0:11	-> 2:12	 [id=10 dir=none  dod="-100" lmb_rep="-1*P(0,a___)+K(0,a___)" num="1" particle="g"];
///	5:13	-> 0:14	 [id=11  dod="-100" lmb_rep="-1*P(0,a___)+K(0,a___)+K(1,a___)" num="1" particle="t"];
///	6:24	-> 1:18	 [id=12 dir=none  dod="-100" lmb_rep="-1*P(3,a___)+K(0,a___)" mass="((Q(2,spenso::cind(0)))^2+(Q(2,spenso::cind(1)))^2*-1+(Q(2,spenso::cind(2)))^2*-1+(Q(2,spenso::cind(3)))^2*-1)^(1/2)" num="1" particle="g"];
///	1:15	-> 4:16	 [id=13  dod="-100" lmb_rep="-1*P(3,a___)+K(0,a___)+K(1,a___)+P(1,a___)+P(2,a___)" num="1" particle="t"];
///	2:17	-> 6:23	 [id=14  dod="-100" lmb_id="0" lmb_rep="K(0,a___)" num="1" particle="d"];
///	4:19	-> 3:20	 [id=15  dod="-100" lmb_rep="-1*P(3,a___)+-1*P(5,a___)+K(0,a___)+K(1,a___)+P(1,a___)+P(2,a___)" num="1" particle="t"];
///	3:21	-> 5:22	 [id=16  dod="-100" lmb_rep="-1*P(0,a___)+K(0,a___)+K(1,a___)+P(4,a___)" num="1" particle="t"];
///}
///
/// ```
/// #let (source-style, sink-style, edge-label, ..callbacks) = physics.style(momentum-arrows: true, show-edge-index: true)
/// #draw(
///   layout(g,beta: 100),
///   source-style: source-style,
///   sink-style:   sink-style,
///   edge-label: edge-label,
/// )
/// #let g = parse(a.text).at(0) 
/// #draw(
///   layout(layout(g,beta: 100,epochs:300,steps:200),subgraph: subgraph.bits(g,(false,false,false,false,false,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true))),
///   source-style: source-style,
///   sink-style:   sink-style,
///   edge-label: edge-label,
/// )
/// ````
/// -> dictionary
#let style(
  map: default-map,
  default: default-edge,
  typst-fields: "plain",
  scope: (:),
  orientation-split: true,
  momentum-arrows: false,
  momentum-arrow-offset: momentum-arrow-defaults.offset,
  momentum-arrow-length: momentum-arrow-defaults.length,
  momentum-arrow-ratio: momentum-arrow-defaults.ratio,
  momentum-arrow-stroke: none,
  momentum-arrow-mark: momentum-arrow-defaults.mark,
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
