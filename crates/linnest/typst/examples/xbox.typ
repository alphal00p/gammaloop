#set page(width: auto, height: auto, margin: 5mm)

#import "../src/lib.typ": draw, graph, layout, layouts
#import graph: *
#import "map-style.typ" as feynman

#let base-layout = layouts.options(
  spring: (strength: 4, length: 0.45),
  repulsion: (
    strength: 8,
    centering: 0.005,
    edge-node: 0.01,
    edge-edge: 0.1,
    dangling: 0.75,
    dangling-centroid: 0,
  ),
  constraints: (side-strength: 10),
  labels: (steps: 100, distance: 0.9, spring: 12, repulsion: 8),
  solver: (
    steps: 50,
    epochs: 50,
    step: 0.81,
    max-movement: 0.4,
    cooling: 0.85,
    z-spring: 2,
    z-spring-growth: 1,
  ),
)
#let opened-layout = layouts.options(
  base: base-layout,
  spring: (strength: 11, length: 0.51),
  repulsion: (
    strength: 4.8,
    centering: 0.015,
    edge-node: 0.07,
    edge-edge: 0.08,
    dangling: 1.8,
    dangling-centroid: 3.6,
  ),
  constraints: (side-strength: 5),
  solver: (step: 0.09, max-movement: 0.35, cooling: 0.92, z-spring: 1),
)
#let cut-layout = layouts.options(
  base: base-layout,
  repulsion: (strength: 10, edge-node: 0.5, dangling: 2.0),
)
#let diagram(g, options: base-layout) = draw(
  layout(graph.style(g, ..feynman.graph-style), ..options),
  ..feynman.draw-style,
)

#context {
  let xbox = graph.build(default-edge-data: (particle: "d"), {
    node(<a>, pos: pos(x: start(-2), y: pin(3)))
    node(<b>, pos: pos(x: start(-2), y: pin(-3)))
    node(<c>, pos: pos(x: start(2), y: pin(3)))
    node(<d>, pos: pos(x: start(2), y: pin(-3)))

    edge(source(<a>), pos: pos(x: group("in", side: "-"), y: pin(3)))
    edge(sink(<a>), source(<b>))
    edge(sink(<b>), pos: pos(x: group("in", side: "-"), y: pin(-3)))
    edge(sink(<c>), pos: pos(x: group("out", side: "+"), y: pin(3)))
    edge(source(<c>), sink(<d>))
    edge(source(<d>), pos: pos(x: group("out", side: "+"), y: pin(-3)))

    edge(
      source(<a>),
      sink(<d>),
      particle: "a",
      crossing-under: <bridge>,
      crossing-gap: 0.9,
      momentum-arrow-shift: -0.55,
      momentum-arrow-offset: 0.80,
    )
    edge(
      source(<b>),
      <bridge>,
      sink(<c>),
      particle: "g",
      momentum-arrow-shift: 0.55,
      momentum-arrow-offset: 0.80,
    )
  })

  let opened-in-x = group("in", side: "-", start: -7)
  let opened-out-x = group("out", side: "+", start: 7)
  let opened-cut-y = group("cut", start: 0)
  let xbox-opened = graph.build(default-edge-data: (particle: "d"), {
    node(<a>, pos: pos(x: start(0), y: pin(0)))
    node(<b>, pos: pos(x: start(-2), y: pin(-3)))
    node(<c>, pos: pos(x: start(-2), y: pin(3)))
    node(<d>, pos: pos(x: start(2), y: pin(0)))

    edge(source(<a>), sink(<c>), momentum-arrow-shift: -0.55)
    edge(
      sink(<a>),
      source(<b>),
      pos: pos(x: group("lower-inner", side: "+", start: 0.5)),
      momentum-arrow-shift: 0.35,
    )
    edge(sink(<b>), pos: pos(x: opened-in-x, y: pin(-4)), bend: 0.18)
    edge(sink(<c>), pos: pos(x: opened-in-x, y: pin(4)))
    edge(
      sink(<d>),
      pos: pos(x: opened-out-x, y: pin(4)),
      momentum-label-shift: -1.8,
    )
    edge(
      source(<d>),
      pos: pos(x: opened-out-x, y: pin(-4)),
      bend: -0.18,
      crossing-under: <bridge>,
      crossing-gap: 0.9,
      fermion-arrow-shift: 1.15,
    )

    edge(
      source(<a>),
      sink(<d>),
      particle: "a",
      momentum-arrow-shift: -0.45,
      momentum-arrow-offset: 0.80,
    )
    edge(
      source(<b>),
      <bridge>,
      particle: "g",
      pos: pos(x: opened-out-x, y: opened-cut-y),
      bend: -0.55,
      momentum-label-shift: -1.0,
    )
    edge(sink(<c>), particle: "g", pos: pos(x: opened-in-x, y: opened-cut-y))
  })

  let in-x = group("in", side: "-", start: -7)
  let out-x = group("out", side: "+", start: 7)
  let top-y = group("top", side: "+", start: 5)
  let bottom-y = group("bottom", side: "-", start: -5)
  let upper-middle-y = group("upper-middle", side: "+", start: 2)
  let lower-middle-y = group("lower-middle", side: "-", start: -2)
  let xbox-cut = graph.build(default-edge-data: (particle: "d"), {
    node(<a>, pos: pos(x: pin(-3), y: top-y))
    node(<b>, pos: pos(x: pin(-3), y: bottom-y))
    node(<c>, pos: pos(x: pin(3), y: top-y))
    node(<d>, pos: pos(x: pin(3), y: bottom-y))

    edge(source(<a>), pos: pos(x: in-x, y: top-y), momentum-label-shift: 1.0)
    edge(
      sink(<a>),
      source(<b>),
      crossing-under: <bridge>,
      crossing-gap: 1.5,
      fermion-arrow-shift: 1.15,
      momentum-label-shift: -1.5,
    )
    edge(sink(<b>), pos: pos(x: in-x, y: bottom-y))

    node(<h1>, hidden: true, pos: pos(x: in-x, y: lower-middle-y))
    node(<h2>, hidden: true, pos: pos(x: out-x, y: upper-middle-y))
    edge(
      source(<h1>),
      <bridge>,
      sink(<h2>),
      particle: "g",
      pos: pos(x: pin(0)),
      route: "straight-through",
      momentum-arrow-shift: 0.55,
      momentum-arrow-offset: 0.80,
    )

    edge(sink(<c>), pos: pos(x: out-x, y: top-y))
    edge(
      source(<c>),
      sink(<d>),
      crossing-under: <bridge>,
      crossing-gap: 1.5,
      fermion-arrow-shift: 1.15,
    )
    edge(source(<d>), pos: pos(x: out-x, y: bottom-y))

    edge(
      source(<a>),
      particle: "g",
      pos: pos(x: in-x, y: upper-middle-y),
      momentum-label-shift: 0.8,
    )
    edge(sink(<d>), particle: "g", pos: pos(x: out-x, y: lower-middle-y))
    edge(
      source(<b>),
      sink(<c>),
      particle: "a",
      crossing-under: <bridge>,
      crossing-gap: 1.5,
      momentum-arrow-shift: -0.45,
      momentum-label-shift: 0.8,
      momentum-arrow-offset: 0.80,
    )
  })

  align(center, stack(
    dir: ttb,
    spacing: 20mm,
    diagram(xbox),
    diagram(xbox-opened, options: opened-layout),
    diagram(xbox-cut, options: cut-layout),
  ))
}
