#set page(width: auto, height: auto, margin: 5mm)

#import "../src/lib.typ": draw, graph, layout
#import graph: *
#import "map-style.typ": *

#let layout-box(
  graph_,
  gamma-dangling: 0.75,
  gamma-dangling-centroid: 0,
  gamma-ev: 0.01,
  gamma-ee: 0.1,
  g-center: 0.005,
  k-spring: 4,
  beta: 8,
  directional-force: 10,
  length-scale: 0.45,
  step: 0.81,
  delta: 0.4,
  cool: 0.85,
  z-spring: 2,
  z-spring-growth: 1,
) = position-labels(layout(
  graph_,
  steps: 50,
  epochs: 50,
  step: step,
  delta: delta,
  cool: cool,
  label-steps: 100,
  label-length-scale: 0.9,
  label-spring: 12,
  label-charge: 8,
  gamma-dangling: gamma-dangling,
  gamma-dangling-centroid: gamma-dangling-centroid,
  gamma-ev: gamma-ev,
  k-spring: k-spring,
  beta: beta,
  gamma-ee: gamma-ee,
  g-center: g-center,
  directional-force: directional-force,
  length-scale: length-scale,
  z-spring: z-spring,
  z-spring-growth: z-spring-growth,
))
#let draw-box(graph_) = draw(
  graph_,
  padding: 1.5,
  edge-dangling-tangent: "horizontal",
  source-style: source-style,
  sink-style: sink-style,
)

#context {
  let xbox = style(graph.build({
    node(<a>, pos: pos(x: start(-2), y: pin(3)))
    node(<b>, pos: pos(x: start(-2), y: pin(-3)))
    node(<c>, pos: pos(x: start(2), y: pin(3)))
    node(<d>, pos: pos(x: start(2), y: pin(-3)))

    edge(
      source(<a>),
      particle: "d",
      pos: pos(x: group("in", side: "-"), y: pin(3)),
    )
    edge(sink(<a>), source(<b>), particle: "d")
    edge(
      sink(<b>),
      particle: "d",
      pos: pos(x: group("in", side: "-"), y: pin(-3)),
    )

    edge(
      sink(<c>),
      particle: "d",
      pos: pos(x: group("out", side: "+"), y: pin(3)),
    )
    edge(source(<c>), sink(<d>), particle: "d")
    edge(
      source(<d>),
      particle: "d",
      pos: pos(x: group("out", side: "+"), y: pin(-3)),
    )

    edge(
      source(<a>),
      sink(<d>),
      particle: "a",
      crossing-under: 7,
      crossing-gap: 0.9,
      momentum-arrow-shift: -0.55,
      momentum-arrow-offset: 0.80,
    )
    edge(
      source(<b>),
      sink(<c>),
      particle: "g",
      momentum-arrow-shift: 0.55,
      momentum-arrow-offset: 0.80,
    )
  }))
  xbox = layout-box(xbox)

  let opened-in-x = group("in", side: "-", start: -7)
  let opened-out-x = group("out", side: "+", start: 7)
  let opened-cut-y = group("cut", start: 0)
  let xbox-opened = style(graph.build({
    node(<a>, pos: pos(x: start(0), y: pin(0)))
    node(<b>, pos: pos(x: start(-2), y: pin(-3)))
    node(<c>, pos: pos(x: start(-2), y: pin(3)))
    node(<d>, pos: pos(x: start(2), y: pin(0)))

    edge(
      source(<a>),
      sink(<c>),
      particle: "d",
      momentum-arrow-shift: -0.55,
    )
    edge(
      sink(<a>),
      source(<b>),
      particle: "d",
      pos: pos(x: group("lower-inner", side: "+", start: 0.5)),
      momentum-arrow-shift: 0.35,
    )
    edge(
      sink(<b>),
      particle: "d",
      pos: pos(x: opened-in-x, y: pin(-4)),
      bend: 0.18,
    )
    edge(
      sink(<c>),
      particle: "d",
      pos: pos(x: opened-in-x, y: pin(4)),
    )
    edge(
      sink(<d>),
      particle: "d",
      pos: pos(x: opened-out-x, y: pin(4)),
      momentum-label-shift: -1.8,
    )
    edge(
      source(<d>),
      particle: "d",
      pos: pos(x: opened-out-x, y: pin(-4)),
      bend: -0.18,
      crossing-under: 7,
      crossing-gap: 0.9,
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
      particle: "g",
      pos: pos(x: opened-out-x, y: opened-cut-y),
      bend: -0.55,
      momentum-label-shift: -1.0,
    )
    edge(
      sink(<c>),
      particle: "g",
      pos: pos(x: opened-in-x, y: opened-cut-y),
    )
  }))
  xbox-opened = layout-box(
    xbox-opened,
    gamma-dangling: 1.8,
    gamma-dangling-centroid: 3.6,
    gamma-ev: 0.07,
    gamma-ee: 0.08,
    g-center: 0.015,
    k-spring: 11,
    beta: 4.8,
    directional-force: 5,
    length-scale: 0.51,
    step: 0.09,
    delta: 0.35,
    cool: 0.92,
    z-spring: 1,
    z-spring-growth: 1,
  )
  assert(
    graph.edges(xbox-opened).at(7).pos.y == graph.edges(xbox-opened).at(8).pos.y,
    message: "the two cut gluons must remain in one free y-coordinate group",
  )

  let in-x = group("in", side: "-", start: -7)
  let out-x = group("out", side: "+", start: 7)
  let top-y = group("top", side: "+", start: 5)
  let bottom-y = group("bottom", side: "-", start: -5)
  let upper-middle-y = group("upper-middle", side: "+", start: 2)
  let lower-middle-y = group("lower-middle", side: "-", start: -2)
  let xbox-cut = style(graph.build({
    node(<a>, pos: pos(x: pin(-3), y: top-y))
    node(<b>, pos: pos(x: pin(-3), y: bottom-y))
    node(<c>, pos: pos(x: pin(3), y: top-y))
    node(<d>, pos: pos(x: pin(3), y: bottom-y))

    edge(
      source(<a>),
      particle: "d",
      pos: pos(x: in-x, y: top-y),
      momentum-label-shift: 1.0,
    )
    edge(
      sink(<a>),
      source(<b>),
      particle: "d",
      crossing-under: 3,
      crossing-gap: 1.5,
    )
    edge(
      sink(<b>),
      particle: "d",
      pos: pos(x: in-x, y: bottom-y),
    )

    node(
      <h1>,
      hidden: true,
      pos: pos(x: in-x, y: lower-middle-y),
    )
    node(
      <h2>,
      hidden: true,
      pos: pos(x: out-x, y: upper-middle-y),
    )
    edge(
      source(<h1>),
      sink(<h2>),
      particle: "g",
      pos: pos(x: pin(0)),
      route: "straight-through",
      momentum-arrow-shift: 0.55,
      momentum-arrow-offset: 0.80,
    )

    edge(
      sink(<c>),
      particle: "d",
      pos: pos(x: out-x, y: top-y),
    )
    edge(
      source(<c>),
      sink(<d>),
      particle: "d",
      crossing-under: 3,
      crossing-gap: 1.5,
    )
    edge(
      source(<d>),
      particle: "d",
      pos: pos(x: out-x, y: bottom-y),
    )

    edge(
      source(<a>),
      particle: "g",
      pos: pos(x: in-x, y: upper-middle-y),
    )
    edge(
      sink(<d>),
      particle: "g",
      pos: pos(x: out-x, y: lower-middle-y),
    )
    edge(
      source(<b>),
      sink(<c>),
      particle: "a",
      crossing-under: 3,
      crossing-gap: 1.5,
      momentum-arrow-shift: -0.45,
      momentum-arrow-offset: 0.80,
    )
  }))
  xbox-cut = layout-box(
    xbox-cut,
    gamma-dangling: 2.0,
    gamma-ev: 0.5,
    beta: 10,
  )

  align(center, stack(
    dir: ttb,
    spacing: 20mm,
    draw-box(xbox),
    draw-box(xbox-opened),
    draw-box(xbox-cut),
  ))
}
