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

#let label-content(value, edge) = {
  if value == none {
    none
  } else if type(value) == content {
    value
  } else {
    [#interpolate-template(value, edge)]
  }
}

#let source-style(edge) = edge-entry(edge).source
#let sink-style(edge) = edge-entry(edge).sink
#let edge-label(edge) = {
  let label-template = edge.at("display-label", default: edge.at("label-template", default: none))
  if label-template != none {
    label-content(label-template, edge)
  } else if edge.at("label", default: none) != none {
    label-content(edge.at("label"), edge)
  } else {
    label-content(edge-entry(edge).label, edge)
  }
}
```

The important part is that `source-style`, `sink-style`, and `edge-label` are
ordinary Typst functions. Linnest calls them with the merged edge data:
DOT statements such as `particle`, plus fields like `eid`, `orientation`,
`source-endpoint`, `sink-endpoint`, and `edge`.

Interpolation is still data-driven: DOT statements can provide
`display-label` or `label-template` with `{field}` placeholders, for example
`display-label="{particle} edge {id}"`. The callback replaces placeholders from
the edge data without evaluating the string as Typst code. Literal braces can
be written as `{{` and `}}`.

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
      source-style: edge-style.source-style,
      sink-style: edge-style.sink-style,
      edge-label: edge-style.edge-label,
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
