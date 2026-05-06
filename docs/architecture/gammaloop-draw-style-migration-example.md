# GammaLoop Draw Style Migration Example

This document shows the intended shape of the GammaLoop drawing migration from
Fletcher `eval_*` strings to pure Typst style callbacks consumed by the
Linnest/CeTZ `draw` API.

## Generated `edge-style.typ`

GammaLoop can keep generating one model-specific file, but the file should
export regular Typst callbacks instead of strings that are later evaluated.

```typ
#import "@preview/mitex:0.2.6": *

#let massive = 1.0pt
#let massless = 0.55pt
#let dashed = (0.1em, 0.45em)
#let dotted = "dotted"

#let stroke-style(c: black, thickness: massless, dash: none) = {
  let base = (paint: c, thickness: thickness, cap: "round")
  if dash == none {
    (stroke: base)
  } else {
    (stroke: base + (dash: dash))
  }
}

#let source-stroke(c: black, thickness: massless, dash: none) = {
  stroke-style(c: c, thickness: thickness, dash: dash)
}

#let sink-stroke(c: black, thickness: massless, dash: none) = {
  stroke-style(c: c.lighten(45%), thickness: thickness, dash: dash)
}

#let wave = (
  pattern: "wave",
  pattern-amplitude: 0.14,
  pattern-wavelength: 0.55,
)

#let coil = (
  pattern: "coil",
  pattern-amplitude: 0.14,
  pattern-wavelength: 0.55,
  pattern-coil-longitudinal-scale: 1.6,
)

#let zigzag = (
  pattern: "zigzag",
  pattern-amplitude: 0.14,
  pattern-wavelength: 0.55,
)

#let default-edge = (
  source: source-stroke(),
  sink: sink-stroke(),
  label: none,
)

// This map is still generated from the model in Rust.
#let map = (
  "a": (
    source: source-stroke(c: black) + wave,
    sink: sink-stroke(c: black) + wave,
    label: mi(`gamma`),
  ),
  "g": (
    source: source-stroke(c: black) + coil,
    sink: sink-stroke(c: black) + coil,
    label: mi(`g`),
  ),
  "W+": (
    source: source-stroke(c: blue, thickness: massive) + zigzag,
    sink: sink-stroke(c: blue, thickness: massive) + zigzag,
    label: mi(`W^+`),
  ),
  "phi": (
    source: source-stroke(c: black, dash: dashed),
    sink: sink-stroke(c: black, dash: dashed),
    label: mi(`phi`),
  ),
  "ghG": (
    source: source-stroke(c: black, dash: dotted),
    sink: sink-stroke(c: black, dash: dotted),
    label: mi(`c`),
  ),
)

#let particle-name(edge) = {
  let particle = edge.at("particle", default: none)
  if particle == none {
    none
  } else {
    str(particle).trim("\"")
  }
}

#let edge-entry(edge) = {
  let particle = particle-name(edge)
  if particle == none {
    default-edge
  } else {
    map.at(particle, default: default-edge)
  }
}

#let eval-scope = (
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
  map: map,
)

#let text-value(value) = str(value).trim("\"")

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

#let eval-template(template, edge) = eval(interpolate-template(template, edge), scope: eval-scope + edge)

#let _assert-typst-fields-mode(mode) = {
  if mode != "plain" and mode != "eval" {
    panic("typst-fields must be \"plain\" or \"eval\"")
  }
}

#let label-content(value, edge, mode: "plain") = {
  _assert-typst-fields-mode(mode)
  if value == none {
    none
  } else if type(value) == content {
    value
  } else {
    let text = interpolate-template(value, edge)
    if mode == "eval" {
      eval(text, scope: eval-scope + edge)
    } else {
      [#text]
    }
  }
}

#let style-dict(value, edge, mode: "plain") = {
  _assert-typst-fields-mode(mode)
  if value == none {
    (:)
  } else if type(value) == dictionary {
    value
  } else if mode == "eval" {
    eval-template(value, edge)
  } else {
    (:)
  }
}

#let source-style(edge, typst-fields: "plain") = {
  let style = edge-entry(edge).source
  style = style + style-dict(edge.at("source-style", default: none), edge, mode: typst-fields)
  style + style-dict(edge.at("source-style-eval", default: none), edge, mode: "eval")
}

#let sink-style(edge, typst-fields: "plain") = {
  let style = edge-entry(edge).sink
  style = style + style-dict(edge.at("sink-style", default: none), edge, mode: typst-fields)
  style + style-dict(edge.at("sink-style-eval", default: none), edge, mode: "eval")
}

#let edge-label(edge, typst-fields: "plain") = {
  let label-eval = edge.at("label-eval", default: none)
  if label-eval != none {
    label-content(label-eval, edge, mode: "eval")
  } else {
    let label-template = edge.at(
      "display-label",
      default: edge.at("label-template", default: edge.at("label", default: none)),
    )
    if label-template != none {
      label-content(label-template, edge, mode: typst-fields)
    } else {
      label-content(edge-entry(edge).label, edge)
    }
  }
}
```

The important part is that `source-style`, `sink-style`, and `edge-label` are
ordinary Typst functions. Linnest calls them with the merged edge data:
DOT statements such as `particle`, plus fields like `eid`, `orientation`,
`source-half-edge`, `sink-half-edge`, and `edge`.

Interpolation is still data-driven by default: DOT statements can provide
`display-label` or `label-template` with `{field}` placeholders, for example
`display-label="{particle} edge {eid}"`. The callback replaces placeholders
from the edge data without evaluating the string as Typst code. Literal braces
can be written as `{{` and `}}`.

The embedded `layout` template also accepts `typst-fields: "plain" | "eval"`.
In `"eval"` mode, known render fields are interpolated and then passed to
Typst's `eval`. This applies to `label`, `display-label`, `label-template`,
`source-style`, and `sink-style`. Explicit `label-eval`,
`source-style-eval`, and `sink-style-eval` fields are always evaluated because
their names are already an opt-in to executable Typst. Structural fields such as `particle`, `id`,
`source`, and `sink` remain raw data for lookup and layout.

The default embedded figure template forwards `sys.inputs.typst-fields`, so the
CLI can opt into executable render fields with `--input typst-fields=eval`.

## Embedded `layout.typ` Wiring

The embedded drawing template can import that generated file and pass the
callbacks directly to `draw`.

```typ
#import "linnest.typ": draw, graph, layout as apply-layout
#import "edge-style.typ" as edge-style

#let layout(
  input,
  scope: (:),
  columns: 1fr,
  unit: 1,
  typst-fields: "plain",
  additional-data: (:),
) = {
  let graphs = graph.parse(input)
  let diagrams = ()

  for graph-bytes in graphs {
    let graph-bytes = apply-layout(graph-bytes, ..additional-data)
    diagrams.push(draw(
      graph-bytes,
      scope: scope,
      unit: unit,
      title: auto,
      source-style: edge => edge-style.source-style(edge, typst-fields: typst-fields),
      sink-style: edge => edge-style.sink-style(edge, typst-fields: typst-fields),
      edge-label: edge => edge-style.edge-label(edge, typst-fields: typst-fields),
    ))
  }

  for diagram in diagrams {
    diagram
  }
}
```

With this shape, GammaLoop owns the particle-to-style policy in one generated
file, while Linnest owns all geometry: edge splitting, node outsets, label
placement, subgraph underlays, and Kurvst path patterns.
