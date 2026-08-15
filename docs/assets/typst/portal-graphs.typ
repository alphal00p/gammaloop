#import "../../../crates/linnest/typst/src/lib.typ": draw, graph, physics
#import graph: build, edge, node, pos, sink, source

#let variant = sys.inputs.at("variant", default: "amplitude")
#let dark = sys.inputs.at("theme", default: "light") == "dark"
#let ink = if dark { rgb("#f8effa") } else { rgb("#3d2645") }
#let accent = if dark { rgb("#b893c7") } else { rgb("#745080") }
#let soft = if dark { rgb("#d8b9e3") } else { rgb("#745080") }

#let page-size = if variant == "field" { (240mm, 150mm) } else { (96mm, 54mm) }
#set page(
  width: page-size.at(0),
  height: page-size.at(1),
  margin: 0pt,
  fill: none,
)

#let line(c, thickness: 0.72pt, dash: none) = physics.stroke-style(
  c: c,
  thickness: thickness,
  dash: dash,
)
#let flow-mark = (
  end: (
    symbol: ">",
    fill: ink,
    stroke: ink + 0.18pt,
    anchor: "center",
    shorten-to: auto,
  ),
  scale: 0.68,
)
#let photon = line(accent, thickness: 0.62pt) + physics.wave
#let gluon = line(soft, thickness: 0.62pt) + physics.coil
#let fermion = line(ink)
#let plain = line(soft, thickness: 0.62pt)
#let cut = line(accent, thickness: 0.8pt, dash: (0.18em, 0.34em))
#let particles = (
  "a": (source: photon, sink: photon, label: none),
  "g": (source: gluon, sink: gluon, label: none),
  "fermion": (
    source: fermion,
    sink: fermion,
    fermion-arrow: true,
    fermion-arrow-mark: flow-mark,
    label: none,
  ),
  "line": (source: plain, sink: plain, label: none),
  "cut": (source: cut, sink: cut, label: none),
)
#let (source-style, sink-style, ..callbacks) = physics.style(map: particles)

#let graph-drawing(g, unit: 8mm, node-radius: 0.09, padding: 0.2) = draw(
  g,
  unit: unit,
  node-label: none,
  node-radius: node-radius,
  node-fill: accent,
  node-stroke: ink + 0.28pt,
  source-style: source-style,
  sink-style: sink-style,
  edge-label: none,
  edge-omega: 1.05,
  padding: padding,
)

#let amplitude() = {
  let g = build(
    {
      node(<left>, id: 0, pos: pos(x: -2.1, y: 0))
      node(<right>, id: 1, pos: pos(x: 2.1, y: 0))
      edge(
        <incoming>,
        sink(<left>),
        id: 0,
        pos: pos(x: -3.8, y: 0),
        particle: "a",
      )
      edge(
        source(<left>),
        <loop-0>,
        sink(<right>),
        id: 1,
        pos: pos(x: 0, y: -1.45),
        particle: "fermion",
      )
      edge(
        source(<right>),
        <loop-1>,
        sink(<left>),
        id: 2,
        pos: pos(x: 0, y: -0.48),
        particle: "fermion",
      )
      edge(
        source(<left>),
        <loop-2>,
        sink(<right>),
        id: 3,
        pos: pos(x: 0, y: 0.48),
        particle: "fermion",
      )
      edge(
        source(<right>),
        <loop-3>,
        sink(<left>),
        id: 4,
        pos: pos(x: 0, y: 1.45),
        particle: "fermion",
      )
      edge(
        source(<right>),
        <outgoing>,
        id: 5,
        pos: pos(x: 3.8, y: 0),
        particle: "a",
      )
    },
    name: "three-loop amplitude",
  )
  graph-drawing(g, unit: 10.5mm)
}

#let cross-section() = {
  let g = build(
    {
      node(<left-top>, id: 0, pos: pos(x: -2.1, y: -0.85))
      node(<left-bottom>, id: 1, pos: pos(x: -2.1, y: 0.85))
      node(<right-top>, id: 2, pos: pos(x: 2.1, y: -0.85))
      node(<right-bottom>, id: 3, pos: pos(x: 2.1, y: 0.85))

      edge(
        <in-0>,
        sink(<left-top>),
        id: 0,
        pos: pos(x: -3.9, y: -0.85),
        particle: "g",
      )
      edge(
        <in-1>,
        sink(<left-bottom>),
        id: 1,
        pos: pos(x: -3.9, y: 0.85),
        particle: "g",
      )
      edge(
        source(<right-top>),
        <out-0>,
        id: 2,
        pos: pos(x: 3.9, y: -0.85),
        particle: "g",
      )
      edge(
        source(<right-bottom>),
        <out-1>,
        id: 3,
        pos: pos(x: 3.9, y: 0.85),
        particle: "g",
      )

      edge(
        source(<left-top>),
        sink(<left-bottom>),
        id: 4,
        pos: pos(x: -3.0, y: 0),
        particle: "fermion",
      )
      edge(
        source(<left-bottom>),
        sink(<left-top>),
        id: 5,
        pos: pos(x: -2.1, y: 0),
        particle: "fermion",
      )
      edge(
        source(<left-top>),
        sink(<left-bottom>),
        id: 6,
        pos: pos(x: -1.2, y: 0),
        particle: "fermion",
      )
      edge(
        source(<right-bottom>),
        sink(<right-top>),
        id: 7,
        pos: pos(x: 1.2, y: 0),
        particle: "fermion",
      )
      edge(
        source(<right-top>),
        sink(<right-bottom>),
        id: 8,
        pos: pos(x: 2.1, y: 0),
        particle: "fermion",
      )
      edge(
        source(<right-bottom>),
        sink(<right-top>),
        id: 9,
        pos: pos(x: 3.0, y: 0),
        particle: "fermion",
      )

      edge(
        source(<left-top>),
        sink(<right-top>),
        id: 10,
        pos: pos(x: 0, y: -1.45),
        particle: "cut",
      )
      edge(
        source(<right-bottom>),
        sink(<left-bottom>),
        id: 11,
        pos: pos(x: 0, y: 1.45),
        particle: "cut",
      )
    },
    name: "five-loop cross section",
  )
  graph-drawing(g, unit: 10mm)
}

#let topology-field() = {
  let g = build(
    {
      // A field of small, disconnected loop topologies. Their placement is
      // first-class graph data, so the background is deterministic.
      node(id: 0, pos: pos(x: -10, y: -5))
      node(id: 1, pos: pos(x: -7.8, y: -5))
      edge(
        source(0),
        sink(1),
        id: 0,
        pos: pos(x: -8.9, y: -5.8),
        particle: "line",
      )
      edge(
        source(1),
        sink(0),
        id: 1,
        pos: pos(x: -8.9, y: -5),
        particle: "line",
      )
      edge(
        source(0),
        sink(1),
        id: 2,
        pos: pos(x: -8.9, y: -4.2),
        particle: "line",
      )

      node(id: 2, pos: pos(x: -5.6, y: -5))
      node(id: 3, pos: pos(x: -2.8, y: -5))
      edge(
        source(2),
        sink(3),
        id: 3,
        pos: pos(x: -4.2, y: -6.2),
        particle: "fermion",
      )
      edge(
        source(3),
        sink(2),
        id: 4,
        pos: pos(x: -4.2, y: -5.4),
        particle: "fermion",
      )
      edge(
        source(2),
        sink(3),
        id: 5,
        pos: pos(x: -4.2, y: -4.6),
        particle: "fermion",
      )
      edge(
        source(3),
        sink(2),
        id: 6,
        pos: pos(x: -4.2, y: -3.8),
        particle: "fermion",
      )

      node(id: 4, pos: pos(x: 0.2, y: -5.7))
      node(id: 5, pos: pos(x: 2.1, y: -4.1))
      node(id: 6, pos: pos(x: 4.0, y: -5.7))
      edge(
        source(4),
        sink(5),
        id: 7,
        pos: pos(x: 1.0, y: -5.0),
        particle: "line",
      )
      edge(
        source(5),
        sink(6),
        id: 8,
        pos: pos(x: 3.1, y: -5.0),
        particle: "line",
      )
      edge(
        source(6),
        sink(4),
        id: 9,
        pos: pos(x: 2.1, y: -6.6),
        particle: "line",
      )
      edge(
        source(4),
        sink(6),
        id: 10,
        pos: pos(x: 2.1, y: -5.5),
        particle: "cut",
      )

      node(id: 7, pos: pos(x: 7.0, y: -5.8))
      node(id: 8, pos: pos(x: 8.7, y: -4.2))
      node(id: 9, pos: pos(x: 10.4, y: -5.8))
      edge(source(7), sink(8), id: 11, pos: pos(x: 7.8, y: -5.0), particle: "g")
      edge(source(8), sink(9), id: 12, pos: pos(x: 9.6, y: -5.0), particle: "g")
      edge(source(9), sink(7), id: 13, pos: pos(x: 8.7, y: -6.7), particle: "g")

      node(id: 10, pos: pos(x: -10.2, y: -0.8))
      node(id: 11, pos: pos(x: -7.4, y: -0.8))
      edge(
        <field-in-0>,
        sink(10),
        id: 14,
        pos: pos(x: -11.5, y: -0.8),
        particle: "a",
      )
      edge(
        source(10),
        sink(11),
        id: 15,
        pos: pos(x: -8.8, y: -2.0),
        particle: "fermion",
      )
      edge(
        source(11),
        sink(10),
        id: 16,
        pos: pos(x: -8.8, y: -1.2),
        particle: "fermion",
      )
      edge(
        source(10),
        sink(11),
        id: 17,
        pos: pos(x: -8.8, y: -0.4),
        particle: "fermion",
      )
      edge(
        source(11),
        sink(10),
        id: 18,
        pos: pos(x: -8.8, y: 0.4),
        particle: "fermion",
      )
      edge(
        source(11),
        <field-out-0>,
        id: 19,
        pos: pos(x: -6.1, y: -0.8),
        particle: "a",
      )

      node(id: 12, pos: pos(x: -4.7, y: -1.8))
      node(id: 13, pos: pos(x: -2.5, y: -1.8))
      node(id: 14, pos: pos(x: -2.5, y: 0.4))
      node(id: 15, pos: pos(x: -4.7, y: 0.4))
      edge(
        source(12),
        sink(13),
        id: 20,
        pos: pos(x: -3.6, y: -2.25),
        particle: "line",
      )
      edge(
        source(13),
        sink(14),
        id: 21,
        pos: pos(x: -2.05, y: -0.7),
        particle: "line",
      )
      edge(
        source(14),
        sink(15),
        id: 22,
        pos: pos(x: -3.6, y: 0.85),
        particle: "line",
      )
      edge(
        source(15),
        sink(12),
        id: 23,
        pos: pos(x: -5.15, y: -0.7),
        particle: "line",
      )
      edge(
        source(12),
        sink(14),
        id: 24,
        pos: pos(x: -3.6, y: -0.7),
        particle: "cut",
      )

      node(id: 16, pos: pos(x: 0.5, y: -1.8))
      node(id: 17, pos: pos(x: 3.3, y: -1.8))
      edge(
        source(16),
        sink(17),
        id: 25,
        pos: pos(x: 1.9, y: -3.0),
        particle: "fermion",
      )
      edge(
        source(17),
        sink(16),
        id: 26,
        pos: pos(x: 1.9, y: -2.2),
        particle: "fermion",
      )
      edge(
        source(16),
        sink(17),
        id: 27,
        pos: pos(x: 1.9, y: -1.4),
        particle: "fermion",
      )
      edge(
        source(17),
        sink(16),
        id: 28,
        pos: pos(x: 1.9, y: -0.6),
        particle: "fermion",
      )

      node(id: 18, pos: pos(x: 6.0, y: -1.8))
      node(id: 19, pos: pos(x: 8.0, y: -0.2))
      node(id: 20, pos: pos(x: 10.0, y: -1.8))
      edge(
        source(18),
        sink(19),
        id: 29,
        pos: pos(x: 6.9, y: -1.0),
        particle: "line",
      )
      edge(
        source(19),
        sink(20),
        id: 30,
        pos: pos(x: 9.1, y: -1.0),
        particle: "line",
      )
      edge(
        source(20),
        sink(18),
        id: 31,
        pos: pos(x: 8.0, y: -2.7),
        particle: "line",
      )
      edge(
        source(18),
        sink(20),
        id: 32,
        pos: pos(x: 8.0, y: -1.6),
        particle: "line",
      )

      node(id: 21, pos: pos(x: -9.8, y: 3.5))
      node(id: 22, pos: pos(x: -7.0, y: 3.5))
      edge(
        source(21),
        sink(22),
        id: 33,
        pos: pos(x: -8.4, y: 2.3),
        particle: "g",
      )
      edge(
        source(22),
        sink(21),
        id: 34,
        pos: pos(x: -8.4, y: 3.1),
        particle: "g",
      )
      edge(
        source(21),
        sink(22),
        id: 35,
        pos: pos(x: -8.4, y: 3.9),
        particle: "g",
      )
      edge(
        source(22),
        sink(21),
        id: 36,
        pos: pos(x: -8.4, y: 4.7),
        particle: "g",
      )

      node(id: 23, pos: pos(x: -4.5, y: 2.7))
      node(id: 24, pos: pos(x: -2.8, y: 4.4))
      node(id: 25, pos: pos(x: -1.1, y: 2.7))
      edge(
        source(23),
        sink(24),
        id: 37,
        pos: pos(x: -3.7, y: 3.5),
        particle: "fermion",
      )
      edge(
        source(24),
        sink(25),
        id: 38,
        pos: pos(x: -1.9, y: 3.5),
        particle: "fermion",
      )
      edge(
        source(25),
        sink(23),
        id: 39,
        pos: pos(x: -2.8, y: 1.8),
        particle: "fermion",
      )
      edge(
        source(23),
        sink(25),
        id: 40,
        pos: pos(x: -2.8, y: 2.9),
        particle: "cut",
      )

      node(id: 26, pos: pos(x: 1.4, y: 3.5))
      node(id: 27, pos: pos(x: 4.2, y: 3.5))
      edge(
        <field-in-1>,
        sink(26),
        id: 41,
        pos: pos(x: 0.1, y: 3.5),
        particle: "a",
      )
      edge(
        source(26),
        sink(27),
        id: 42,
        pos: pos(x: 2.8, y: 2.3),
        particle: "line",
      )
      edge(
        source(27),
        sink(26),
        id: 43,
        pos: pos(x: 2.8, y: 3.1),
        particle: "line",
      )
      edge(
        source(26),
        sink(27),
        id: 44,
        pos: pos(x: 2.8, y: 3.9),
        particle: "line",
      )
      edge(
        source(27),
        sink(26),
        id: 45,
        pos: pos(x: 2.8, y: 4.7),
        particle: "line",
      )
      edge(
        source(27),
        <field-out-1>,
        id: 46,
        pos: pos(x: 5.5, y: 3.5),
        particle: "a",
      )

      node(id: 28, pos: pos(x: 7.2, y: 2.7))
      node(id: 29, pos: pos(x: 9.0, y: 4.3))
      node(id: 30, pos: pos(x: 10.8, y: 2.7))
      edge(
        source(28),
        sink(29),
        id: 47,
        pos: pos(x: 8.0, y: 3.5),
        particle: "line",
      )
      edge(
        source(29),
        sink(30),
        id: 48,
        pos: pos(x: 10.0, y: 3.5),
        particle: "line",
      )
      edge(
        source(30),
        sink(28),
        id: 49,
        pos: pos(x: 9.0, y: 1.8),
        particle: "line",
      )
      edge(
        source(28),
        sink(30),
        id: 50,
        pos: pos(x: 9.0, y: 2.9),
        particle: "cut",
      )
    },
    name: "loop topology field",
  )
  graph-drawing(g, unit: 8.4mm, node-radius: 0.055, padding: 0.08)
}

#let artwork = if variant == "amplitude" {
  amplitude()
} else if variant == "cross-section" {
  cross-section()
} else if variant == "field" {
  topology-field()
} else {
  panic("unknown portal graph variant: " + variant)
}

#align(center + horizon, artwork)
