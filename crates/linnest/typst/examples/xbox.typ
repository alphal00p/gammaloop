#set page(width: auto, height: auto, margin: 5mm)

#import "../src/lib.typ": draw, graph, layout
#import graph: *
#import "map-style.typ": *

#let layout-box(
  graph_,
  gamma-dangling: 0.75,
  gamma-ev: 0.01,
  beta: 8,
  length-scale: 0.45,
) = position-labels(layout(
  graph_,
  steps: 50,
  epochs: 50,
  label-steps: 100,
  label-length-scale: 0.9,
  label-spring: 12,
  label-charge: 8,
  gamma-dangling: gamma-dangling,
  gamma-ev: gamma-ev,
  k-spring: 4,
  beta: beta,
  directional-force: 10,
  length-scale: length-scale,
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
      momentum-arrow-shift: 0.35,
    )
    edge(
      sink(<b>),
      particle: "d",
      pos: pos(x: group("in", side: "-"), y: pin(-4)),
      bend: 0.18,
    )
    edge(
      sink(<c>),
      particle: "d",
      pos: pos(x: group("in", side: "-"), y: pin(4)),
    )
    edge(
      sink(<d>),
      particle: "d",
      pos: pos(x: group("out", side: "+"), y: pin(4)),
      momentum-label-shift: -1.8,
    )
    edge(
      source(<d>),
      particle: "d",
      pos: pos(x: group("out", side: "+"), y: pin(-4)),
      bend: -0.18,
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
      pos: pos(x: group("out", side: "+"), y: pin(0)),
      bend: -0.55,
      momentum-label-shift: -1.0,
    )
    edge(
      sink(<c>),
      particle: "g",
      pos: pos(x: group("in", side: "-"), y: pin(0)),
    )
  }))
  xbox-opened = layout-box(
    xbox-opened,
    gamma-dangling: 3.0,
    gamma-ev: 0.05,
    beta: 10,
    length-scale: 0.4,
  )

  let xbox-cut = style(graph.build({
    node(<a>, pos: pos(x: pin(-3), y: pin(5)))
    node(<b>, pos: pos(x: pin(-3), y: pin(-5)))
    node(<c>, pos: pos(x: pin(3), y: pin(5)))
    node(<d>, pos: pos(x: pin(3), y: pin(-5)))

    edge(
      source(<a>),
      particle: "d",
      pos: pos(x: pin(-7), y: pin(5)),
      momentum-label-shift: 1.0,
    )
    edge(sink(<a>), source(<b>), particle: "d")
    edge(
      sink(<b>),
      particle: "d",
      pos: pos(x: pin(-7), y: pin(-5)),
    )

    node(
      <h1>,
      hidden: true,
      pos: pos(x: pin(-7), y: pin(-2)),
    )
    node(
      <h2>,
      hidden: true,
      pos: pos(x: pin(7), y: pin(2)),
    )
    edge(
      source(<h1>),
      sink(<h2>),
      particle: "g",
      cut: true,
      cut-gap: 0.70,
      pos: pos(x: pin(0), y: pin(0)),
      route: "straight-through",
      momentum-arrow-shift: 0.55,
      momentum-arrow-offset: 0.80,
    )

    edge(
      sink(<c>),
      particle: "d",
      pos: pos(x: pin(7), y: pin(5)),
    )
    edge(source(<c>), sink(<d>), particle: "d")
    edge(
      source(<d>),
      particle: "d",
      pos: pos(x: pin(7), y: pin(-5)),
    )

    edge(
      source(<a>),
      particle: "g",
      pos: pos(x: pin(-7), y: pin(2)),
    )
    edge(
      sink(<d>),
      particle: "g",
      pos: pos(x: pin(7), y: pin(-2)),
    )
    edge(
      source(<b>),
      sink(<c>),
      particle: "a",
      momentum-arrow-shift: -0.45,
      momentum-arrow-offset: 0.80,
    )
  }))
  xbox-cut = layout-box(
    xbox-cut,
    gamma-dangling: 2.0,
    gamma-ev: 0.05,
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
