#set page("a4", margin: 20mm)

#import "../src/lib.typ": draw, graph, layout, layouts
#import graph: *
#import "map-style.typ" as feynman
#import "@preview/cetz:0.5.1" as cetz

// Native font and stroke sizes are independent of the graph coordinate unit.
#let diagram-style = (
  font-size: 6pt,
  unit: 2.6mm,
  line-width: 0.5pt,
  cut-line-width: 0.6pt,
  padding: 0.3,
)
#set text(size: diagram-style.font-size)

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
    depth-scale: 1,
    flattening-end: 0.5,
  ),
)

// `cut-x` is in graph units; auto puts the cut just right of center.
#let diagram(g, options: base-layout, cut-x: auto) = draw(
  layout(
    graph.style(
      g,
      ..feynman.graph-style,
      unit: diagram-style.unit,
      node-style: node => {
        let style = feynman.node-style(node)
        if style.stroke != none {
          style.stroke += (thickness: diagram-style.line-width)
        }
        style
      },
    ),
    ..options,
  ),
  ..feynman.draw-style,
  padding: diagram-style.padding,
  draw-after: g => {
    // Exclude the invisible compactification spring and its anchor nodes.
    let nodes = graph.nodes(g).filter(n => n.name in (<a>, <b>, <c>, <d>))
    let edges = graph
      .edges(g)
      .filter(e => e.data.at("edge-style", default: auto) != none)
    let xs = nodes.map(n => n.pos.x)
    let ys = (nodes + edges).map(item => item.pos.y)
    let x = if cut-x == auto {
      calc.min(..xs) + 0.6 * (calc.max(..xs) - calc.min(..xs))
    } else { cut-x }
    cetz.draw.line(
      (x, calc.max(..ys) + 1.5),
      (x, calc.min(..ys) - 1.5),
      stroke: (
        paint: red,
        thickness: diagram-style.cut-line-width,
        dash: "dashed",
      ),
    )
  },
)

#context {

  let D1 =[$k-p_1$]
  let D2 = [$k-p_2$]
  let D3 = [$p_1$]
  let D4 = [$p_2$]
  let D5 = [$p_(12) - k$]
  let D6 = [$k$]


  let in-x = group("in", side: "-", start: -7)
  let out-x = group("out", side: "+", start: 7)
  let top = group("top", start: 5)
  let mid = group("mid", start: 0)
  let bot = group("bot", start: -5)



  let xbox = graph.build(default-edge-data: (particle: "d",momentum-arrow-offset: 0.4), {
    node(<a>, pos: pos(x: start(-2), y: top))
    node(<b>, pos: pos(x: start(-2), y: bot))
    node(<c>, pos: pos(x: start(2), y: top))
    node(<d>, pos: pos(x: start(2), y: bot))

    edge(
      source(<c>),
      momentum: D3,
      spring-length: .3,
      orientation: "reversed",
      momentum-label-gap: 0.4,
      pos: pos(x:out-x,y: top)
    )
    edge(
      sink(<a>),
      momentum: D3,
      reverse: true,
      orientation: "reversed",
      spring-length: .3,
      momentum-label-gap: 0.4,
      momentum-arrow-side: "left",
      pos: pos(x:in-x, y: top),
    )
    edge(source(<a>), sink(<b>),orientation: "reversed",  momentum: D2, momentum-label-gap: 0.25)
    edge(source(<c>), sink(<d>), momentum: D1,momentum-label-gap: 0.25)
    edge(
      source(<d>),
      spring-length: .3,
      momentum: D4,
      momentum-label-gap: 0.2,
      momentum-arrow-side: "right",
      pos: pos(x:out-x,y: bot )
    )
    edge(
      sink(<b>),
      momentum: D4,
      momentum-label-gap: 0.2,
      spring-length: .3,
      pos: pos(x: in-x, y: bot),
    )
    edge(
      source(<a>),
      sink(<d>),
      momentum: D5,
      particle: "a",
      crossing-under: <bridge>,
      crossing-gap: 0.9,
      momentum-arrow-shift: 1.5,
      momentum-label-gap: 0.05,
      momentum-arrow-side: "right",
      // momentum-arrow-offset: 0.80,
    )
    edge(
      source(<c>),
      <bridge>,
      sink(<b>),
      momentum: D6,
      momentum-arrow-side: "right",
      particle: "g",
      momentum-label-gap: 0.2,
      momentum-arrow-shift: -1.5,
    )
  })

  let xbox-opened = graph.build(default-edge-data: (particle: "d",momentum-arrow-offset: 0.4), {
    node(<a> )
    node(<b>)
    node(<c>)
    node(<d>)

    edge(source(<c>),sink(<a>), momentum-arrow-shift: 1.,  momentum:D3, orientation: "reversed")
    edge(
      sink(<a>),
      source(<b>),
      orientation: "reversed",  momentum: D2, momentum-label-gap: 0.2,
      spring-length: .3,
    )
    edge(
      source(<c>),
      momentum: D1,
      pos: pos(x: in-x, y: top,z:pin(0)),
      momentum-arrow-side: "right", momentum-label-gap: 0.15,  momentum-label-shift: -1.5,
    )
    edge(
      sink(<d>),
      momentum: D1,
      pos: pos(x: out-x, y: top, z: pin(0)),
      momentum-arrow-side: "left",
      momentum-arrow-length: 1.4, momentum-arrow-shift: -.4,
    )
  
    edge(
      source(<d>),
      momentum: D4,
      pos: pos(x: out-x, y: bot, z: pin(0)),
      bend: -0.18,
      crossing-under: <bridge>,
      crossing-gap: 0.9,
       momentum-arrow-length: 1.4, momentum-arrow-shift: .5,
        momentum-arrow-side: "left",
      // fermion-arrow-shift: 1.15,
    )
    edge(
      sink(<b>),
      momentum: D4,
      pos: pos(x: in-x, y: bot, z: pin(0)),
      spring-length: .3,
      momentum-arrow-shift: -1.3,
      momentum-arrow-side: "left",
    )

    edge(
      source(<a>),
      sink(<d>),
      momentum: D5,
      spring-length: 1.5,
      particle: "a",
      // momentum-arrow-offset: 0.80,
    )
    edge(
      source(<b>),
      <bridge>,
      momentum: D6,
      particle: "g",
      pos: pos(x: out-x, y: mid, z: pin(0)),
      bend: -0.55,
      momentum-arrow-side: "right",

    )
    edge(
      sink(<c>),
      momentum: D6,
      particle: "g",
      pos: pos(x: in-x, y:mid, z: pin(0)),
      momentum-arrow-side: "left",
      momentum-arrow-length: 1.4, momentum-arrow-shift: -0.4,
    )

    // Pull the external rows together without drawing another propagator.
    node(<compact-top>, hidden: true, pos: pos(x: pin(0), y: top))
    node(<compact-bottom>, hidden: true, pos: pos(x: pin(0), y: bot))
    edge(
      source(<compact-top>),
      sink(<compact-bottom>),
      momentum: [],
      style: none,
      spring-length: 0.25,
      pos: pos(x: pin(0)),
    )
  })

 
  let xbox-opened2 = graph.build(default-edge-data: (particle: "d",momentum-arrow-offset: 0.4), {
    node(<a> )
    node(<b>)
    node(<c>, pos: pos(y:start(-4)))
    node(<d>, pos: pos(y:start(0) ))
    
    edge(
      source(<c>),
      momentum: D3,
      orientation: "reversed",
      pos: pos(x: out-x, y: top, z: pin(0)),
      crossing-under: <bridge>,
      crossing-gap: 0.9,
    )
    edge(sink(<a>), momentum: D3, orientation: "reversed", pos: pos(x: in-x,y: top, z: pin(0)),momentum-arrow-length: 1.2,momentum-label-gap: 0.2, momentum-arrow-shift: -.6, momentum-arrow-side: "left")
    edge(
      sink(<a>),
      pos: pos(x: in-x, y: bot, z: pin(0)),
      momentum: D2,orientation: "reversed",
      momentum-arrow-side: "left",
      momentum-arrow-length: 1.4, momentum-arrow-shift: -1.4, momentum-label-gap: 0.2,
    )
    edge(
      source(<b>),
      pos: pos(x: out-x, y: bot, z: pin(0)),
      momentum: D2,orientation: "reversed", momentum-arrow-side: "right",
      spring-length: .3,momentum-arrow-length: .7, momentum-arrow-shift: 0.3, momentum-label-gap: 0.1,momentum-label-shift: -0.7,
    )

    edge(
      source(<c>),
      sink(<d>),
      momentum:D1,
      momentum-arrow-side: "left",momentum-arrow-length: 1.2, momentum-label-gap: 0.2,
    )

    edge(
      source(<d>),
      sink(<b>),
      momentum: D4,
      momentum-arrow-side: "right",momentum-arrow-shift: -0.4, 
      spring-length: .1,momentum-label-gap: 0.1,
    )

    edge(
      source(<a>),
      sink(<d>),
      momentum: D5,
      spring-length: .5,
      particle: "a",
      momentum-arrow-side: "right",momentum-label-gap: 0.05,
      // momentum-arrow-offset: 0.80,
    )
    edge(
      source(<b>),
      <bridge>,
      momentum: D6,
      particle: "g",
      pos: pos(x: out-x, y: mid, z: pin(0)),
      bend: -0.55, momentum-arrow-shift: .8,momentum-arrow-length: 1.5,
      momentum-arrow-side: "left",
    )
    edge(
      sink(<c>),
      momentum: D6,
      particle: "g", momentum-arrow-shift: .8, momentum-label-gap: 0.2,
      pos: pos(x: in-x, y: mid, z: pin(0)),
      momentum-arrow-side: "left",
      
    )

    // Pull the external rows together without drawing another propagator.
    node(<compact-top>, hidden: true, pos: pos(x: start(0), y: top))
    node(<compact-bottom>, hidden: true, pos: pos(x: start(0), y: bot))
    edge(
      source(<compact-top>),
      sink(<compact-bottom>),
      momentum: [],
      style: none,
      spring-length: 0.01,
      pos: pos(x: pin(0)),
    )
  })
  
  let mid2 = group("mid2", start: 2)
  let xbox-cut = graph.build(default-edge-data: (particle: "d", momentum-label-gap: 0.4,momentum-arrow-offset: 0.4), {
    node(<a> )
    node(<b>)
    node(<c>)
    node(<d>)
    
    edge(
      source(<c>),
      sink(<a>),
      momentum: D3,orientation: "reversed",
      spring-length: 1.5,
      crossing-under: <bridge>,
      crossing-gap: 0.8,
      fermion-arrow-shift: -0.5,
      momentum-arrow-shift: -.6,
    )

    edge(source(<a>), momentum:D2, pos: pos(x: in-x, y: bot,z: pin(0)), momentum-arrow-side: "right",momentum-arrow-shift: 0.7,momentum-label-gap: 0.1,orientation: "reversed")
    edge(sink(<b>), momentum: D2, pos: pos(x: out-x, y: bot,z: pin(0)), momentum-arrow-side: "left",momentum-label-gap: 0.1,momentum-label-shift: 0.5,orientation: "reversed", )

    node(<h1>, hidden: true, pos: pos(x: in-x, y: mid,z: pin(0)))
    node(<h2>, hidden: true, pos: pos(x: out-x, y: mid2,z: pin(0)))
    edge(
      sink(<h1>),
      <bridge>,
      source(<h2>),
      momentum: D6,
      spring-length: 3.5,
      particle: "g",
      momentum-arrow-shift: -4,
      momentum-arrow-offset: 0.40,
      momentum-label-offset: 0.10,
    )

    edge(sink(<d>), momentum: D1, momentum-arrow-side: "right",momentum-label-shift: .6,
    momentum-label-gap: 0.2,pos: pos(x: out-x, y: top,z: pin(0)))
    edge(
      source(<c>),
      momentum: D1,
      pos: pos(x: in-x, y: top,z: pin(0)),
      momentum-arrow-side: "right", momentum-label-shift: -.6,
      momentum-label-gap: 0.1,
    )

    edge(
      source(<d>),
      sink(<b>),
      momentum: D4,
      spring-length: 1.5,
      crossing-under: <bridge>,
      crossing-gap: .8,  momentum-arrow-shift: .4, momentum-label-gap: 0.2,
      fermion-arrow-shift: 1.15,
    )
    edge(
      source(<a>),
      sink(<d>),
      momentum: D5,
      particle: "a",
      crossing-under: <bridge>,
      crossing-gap: 1.5,
       momentum-label-gap: 0.01, momentum-arrow-side: "right",
      momentum-arrow-shift: -1.,
    )
    edge(
      sink(<b>),
      momentum: D6, momentum-arrow-side: "left",
      particle: "g",
      pos: pos(x: out-x, y: mid),
      momentum-arrow-shift: -0.5,
    )
    edge(sink(<c>),momentum:D6, particle: "g", pos: pos(x: in-x, y:mid2), momentum-arrow-side: "left",)
    // Pull the external rows together without drawing another propagator.
    node(<compact-top>, hidden: true, pos: pos(x: start(0), y: top))
    node(<compact-bottom>, hidden: true, pos: pos(x: start(0), y: bot))
    edge(
      source(<compact-top>),
      sink(<compact-bottom>),
      momentum: [],
      style: none,
      spring-length: 1.8,
      pos: pos(x: pin(0)),
    )
  })

  let diagrams = $

    #diagram(xbox, cut-x: -1)+
    #diagram(xbox-opened, cut-x: 0)+
    #diagram(xbox-opened2, cut-x: -1)+
    #diagram(xbox-cut) = op("disc")_(p_1^2)  op("disc")_(p_2^2) integral (dif ^4 k  )/(2 pi )^4 (cal(N)_times.square delta_+(k))/(  p_1^2 p_2^2 (p_2-k)^2 (k-p_1)^2)
  $
  diagrams
}
