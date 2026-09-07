#set page(width: auto, height: auto, margin: 5mm)

#import "../src/lib.typ": draw, graph, layout, layouts
#import graph: *
#import "map-style.typ" as feynman

#let base-layout = layouts.options(
  spring: (strength: 18, length: 0.25),
  repulsion: (
    strength: 12,
    centering: 0.0005,
    edge-node: .29,
    edge-edge: .1,
    dangling: 14.75,
    dangling-centroid: 50,
  ),
  constraints: (side-strength: 10),
  labels: (steps: 100, distance: 0.9, spring: 12, repulsion: 8),
  solver: (
    steps: 50,
    epochs: 50,
    step: 0.81,
    max-movement: 0.4,
    cooling: 0.85,
    z-spring: 10,
    z-spring-growth: 1.5,
  ),
)

#let diagram(g, options: base-layout) = draw(
  layout(graph.style(g, ..feynman.graph-style), ..options),
  ..feynman.draw-style,
)

#context {

  let d1 = group("d1", start: -5)
  let d2 = group("d2", start: 5)
  let xbox = graph.build(default-edge-data: (particle: "d"), {
    node(<a>, pos: pos(x: start(-2), y: d2))
    node(<b>, pos: pos(x: start(-2), y: d1))
    node(<c>, pos: pos(x: start(2), y: d2))
    node(<d>, pos: pos(x: start(2), y: d1))

    edge(sink(<a>), 
      momentum:[$p_1$],
      reverse: true,  orientation: "reversed", momentum-arrow-side: "left",
      pos: pos(x: group("in", side: "-"), y: d2)
    )
    edge(sink(<a>), source(<b>),momentum:[$p_2-k$],)
    edge(sink(<b>),momentum:[$p_2$], pos: pos(x: group("in", side: "-"), y: d1))
    edge(source(<c>), momentum:[$p_1$],orientation: "reversed", pos: pos(x: group("out", side: "+"), y: d2))
    edge(source(<c>), sink(<d>),momentum:[$k-p_1$],)
    edge(source(<d>),momentum:[$p_2$], momentum-arrow-side: "right",pos: pos(x: group("out", side: "+"), y: d1))

    edge(
      source(<a>),
      sink(<d>),
      momentum:[$p_1+p_2-k$],
      particle: "a",
      crossing-under: <bridge>,
      crossing-gap: 0.9,
      momentum-arrow-shift: 1.5, momentum-arrow-side: "right",
      // momentum-arrow-offset: 0.80,
    )
    edge(
      source(<b>),
      <bridge>,
      sink(<c>),momentum:[$k$], momentum-arrow-side: "left",
      particle: "g",
      momentum-arrow-shift: 1.5, 
    )
  })

  let opened-in-x = group("in", side: "-", start: -7)
  let opened-out-x = group("out", side: "+", start: 7)
  let opened-cut-y = group("cut", start: 0)
  let xbox-opened = graph.build(default-edge-data: (particle: "d"), {
    node(<a>, pos: pos(x: start(0)))
    node(<b>, pos: pos(x: start(-2)))
    node(<c>, pos: pos(x: start(-2)))
    node(<d>, pos: pos(x: start(2)))

    edge(sink(<a>), source(<c>), momentum:[$p_1$],orientation: "reversed")
    edge(
      sink(<a>),
      source(<b>),momentum:[$p_2-k$],
    )
    edge(sink(<b>), momentum:[$p_2$],pos: pos(x: opened-in-x, y: d1))
    edge(source(<c>),momentum:[$k-p_1$], pos: pos(x: opened-in-x, y: d2), momentum-arrow-side: "right",)
    edge(
      sink(<d>),momentum:[$k-p_1$],
      pos: pos(x: opened-out-x, y: d2),
      momentum-arrow-side: "left",
    )
    edge(
      source(<d>),momentum:[$p_2$],
      pos: pos(x: opened-out-x, y: d1),
      bend: -0.18,
      crossing-under: <bridge>,
      crossing-gap: 0.9,
      fermion-arrow-shift: 1.15,
    )

    edge(
      source(<a>),
      sink(<d>),momentum:[$p_1+p_2-k$],spring-length: 1.5,
      particle: "a",
      momentum-arrow-shift: -0.45,
      momentum-arrow-offset: 0.80,
    )
    edge(
      source(<b>),
      <bridge>,momentum:[$k$],
      particle: "g",
      pos: pos(x: opened-out-x, y: opened-cut-y),
      bend: -0.55,
      momentum-arrow-side: "right",
    )
    edge(sink(<c>), momentum:[$k$],particle: "g", pos: pos(x: opened-in-x, y: opened-cut-y), momentum-arrow-side: "left",)
  })

  let in-x = group("in", side: "-", start: -7)
  
  let in2-x = group("in2", side: "-", start: -3)
  let out-x = group("out", side: "+", start: 7)
  let out2-x = group("out2", side: "+", start: 3)
  let top-y = group("top", side: "+", start: 5)
  let bottom-y = group("bottom", side: "-", start: -5)
  let upper-middle-y = group("upper-middle", side: "+", start: 2)
  let lower-middle-y = group("lower-middle", side: "-", start: -2)
  let xbox-cut = graph.build(default-edge-data: (particle: "d",momentum:[]), {
    node(<a>, pos: pos(x:in2-x, y: top-y))
    node(<b>, pos: pos(x:in2-x, y: bottom-y))
    node(<c>, pos: pos(x: out2-x, y: top-y))
    node(<d>, pos: pos(x: out2-x, y: bottom-y))

    edge(source(<a>), pos: pos(x: in-x, y: top-y), momentum-label-shift: 1.0)
    edge(
      sink(<a>),
      source(<b>),momentum:[],  spring-length: 1.5,
      crossing-under: <bridge>,
      crossing-gap: 1.5,
      fermion-arrow-shift: 1.15,
      momentum-label-shift: -1.5,
    )
    edge(sink(<b>),momentum:[], pos: pos(x: in-x, y: bottom-y))

    node(<h1>, hidden: true, pos: pos(x: in-x, y: lower-middle-y))
    node(<h2>, hidden: true, pos: pos(x: out-x, y: upper-middle-y))
    edge(
      source(<h1>),
      <bridge>,
      sink(<h2>),momentum:[],
        spring-length: 3.5,
      particle: "g",
      pos: pos(x: pin(0)),
      route: "straight-through",
      momentum-arrow-shift: 0.55,
      momentum-arrow-offset: 0.80,
    )

    edge(sink(<c>),momentum:[], pos: pos(x: out-x, y: top-y))
    edge(
      source(<c>),
      sink(<d>),momentum:[],  spring-length: 1.5,
      crossing-under: <bridge>,
      crossing-gap: 1.5,
      fermion-arrow-shift: 1.15,
    )
    edge(source(<d>),momentum:[], pos: pos(x: out-x, y: bottom-y))

    edge(
      source(<a>),momentum:[],
      particle: "g",
      pos: pos(x: in-x, y: upper-middle-y),
      momentum-label-shift: 0.8,
    )
    edge(sink(<d>), particle: "g", pos: pos(x: out-x, y: lower-middle-y))
    edge(
      source(<b>),
      sink(<c>),momentum:[],
      particle: "a",
      crossing-under: <bridge>,
      crossing-gap: 1.5,
      momentum-arrow-shift: -0.45,
      momentum-label-shift: 0.8,
      momentum-arrow-offset: 0.80,
    )
  })

   grid(
   columns: 3,
    diagram(xbox),
    diagram(xbox-opened),
    diagram(xbox-cut),
  )
}
