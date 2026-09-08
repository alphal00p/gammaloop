#set page(height: auto, width: auto, margin: 5mm)

#import "../src/lib.typ": draw, graph, layout, layouts, subgraph
#import graph: *
#import "map-style.typ" as feynman
#import feynman: momentum as mom
#import "@preview/cetz:0.5.1" as cetz

// Native font and stroke sizes are independent of the graph coordinate unit.
#let diagram-style = (
  font-size: 6pt,
  unit: 2.6mm,
  line-width: 0.5pt,
  cut-line-width: 0.6pt,
  endpoint-box: (
    padding: (x: 0.5, y: 0.4),
    radius: 0.25,
    fills: (p1: blue.lighten(85%), p2: orange.lighten(85%)),
  ),
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

#let in-x = group("in", side: "-", start: -7)
#let out-x = group("out", side: "+", start: 7)
#let top = group("top", start: 5)
#let mid = group("mid", start: 0)
#let bot = group("bot", start: -5)
#let mid2 = group("mid2", start: 2)

// The physical K4 is shared; only cuts and drawing patches vary between views.
#let master = {
  node(<a>)
  node(<b>)
  node(<c>)
  node(<d>)
  edge(<D1>, source(<c>), sink(<d>), momentum: [$k-p_1$])
  edge(<D2>, source(<a>), sink(<b>), momentum: [$k-p_2$], orientation: "reversed")
  edge(<D3>, source(<c>), sink(<a>), momentum: [$p_1$], orientation: "reversed")
  edge(<D4>, source(<d>), sink(<b>), momentum: [$p_2$])
  edge(<D5>, source(<a>), sink(<d>), momentum: [$p_(12) - k$], particle: "a")
  edge(<D6>, source(<b>), sink(<c>), momentum: [$k$],particle:"g")
}

#let compact = {
  node(<compact-top>, hidden: true, pos: pos(x: start(0), y: top))
  node(<compact-bottom>, hidden: true, pos: pos(x: start(0), y: bot))
  edge(<compact>, source(<compact-top>), sink(<compact-bottom>),
    momentum: [], style: none, pos: pos(x: pin(0)))
}

// Each entry describes one crossing along the underlying source -> sink flow.
#let cut-side(g, x, crossings, ..selection) = subgraph.with-data(
  g, subgraph.select(g, ..selection), data: (x: x),
  hedge: h => {
    let points = crossings.at(str(h.edge-name))
    (winding: points.len(), crossings: points)
  },
)
#let boundary-position(item) = {
  let b = item.boundary
  (pos: pos(x: b.cut-data.x, y: b.data.crossings.at(b.crossing).y))
}
#let flatten-boundary = item => if item.boundary != none { (pos: pos(z: pin(0))) }

// `cut-x` is in graph units; auto puts the cut just right of center.
#let diagram(g, options: base-layout, cut-x: auto, cut-y: auto) = context {
  draw(
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
      // Group signed external momenta explicitly, including the through-gluon's
      // hidden endpoints. Shading follows the solved positions, not the layout seeds.
      let endpoints = graph.boundaries(g).map(b => b + (
        group: b.data.crossings.at(b.crossing).group,
      ))
      let box-style = diagram-style.endpoint-box
      cetz.draw.on-layer(-1, {
        for (momentum, fill) in box-style.fills {
          for side in ("left", "right") {
            let points = endpoints.filter(b => b.side == side and b.group == momentum).map(b => b.pos)
            if points.len() > 0 {
              let xs = points.map(p => p.x)
              let ys = points.map(p => p.y)
              cetz.draw.rect(
                (calc.min(..xs) - box-style.padding.x, calc.min(..ys) - box-style.padding.y),
                (calc.max(..xs) + box-style.padding.x, calc.max(..ys) + box-style.padding.y),
                radius: box-style.radius,
                fill: fill,
                stroke: none,
              )
            }
          }
        }
      })

      // Exclude the invisible compactification spring and its anchor nodes.
      let nodes = graph.nodes(g).filter(n => n.boundary == none and (
        n.data == none or not n.data.at("hidden", default: false)
      ))
      let edges = graph.edges(g).filter(e => e.data.at("edge-style", default: auto) != none)
      let xs = nodes.map(n => n.pos.x)
      let ys = (nodes + edges).map(item => item.pos.y)
      let x = if cut-x == auto {
        calc.min(..xs) + 0.6 * (calc.max(..xs) - calc.min(..xs))
      } else { cut-x }
      let y = if cut-y == auto { 1.5 } else { cut-y }
      cetz.draw.line(
        (x, calc.max(..ys) + y),
        (x, calc.min(..ys) - y),
        stroke: (paint: red, thickness: diagram-style.cut-line-width, dash: "dashed"),
      )
    },
  )
}

#{
  let xbox = {
    let g = graph.build(default-edge-data: (particle: "d") + mom(offset: .4), master)
    let crossings = (D3: ((group: "p1", y: top),), D4: ((group: "p2", y: bot),))
    let g = graph.cut(g,
      left: cut-side(g, in-x, crossings, sink: (<D3>, <D4>)),
      right: cut-side(g, out-x, crossings, source: (<D3>, <D4>)),
      boundary: boundary-position,
    )
    let nodes = (
      a: (pos: pos(x: start(-2), y: top)), b: (pos: pos(x: start(-2), y: bot)),
      c: (pos: pos(x: start(2), y: top)), d: (pos: pos(x: start(2), y: bot)),
    )
    let edges = (
      D1: mom(label: (gap: .25)), D2: mom(label: (gap: .25)),
      "D3.0": (statements: ("spring-length": .3)) + mom(label: (gap: .4)),
      "D3.1": (statements: ("spring-length": .3), reverse: true)
        + mom(side: "left", label: (gap: .4)),
      "D4.0": (statements: ("spring-length": .3)) + mom(side: "right", label: (gap: .2)),
      "D4.1": (statements: ("spring-length": .3)) + mom(label: (gap: .2)),
      D5: (crossing-under: <D6>, crossing-gap: 0.9) + mom(
        side: "right", shift: 1.5, label: (gap: .05),
        // offset: 0.80,
      ),
      D6: mom(side: "left", shift: 1.5, label: (gap: .2)),
    )
    graph.map(g, node: nodes, edge: edges)
  }

  // Pull the external rows together without drawing another propagator.
  let g = graph.build(default-edge-data: (particle: "d") + mom(offset: .4), master, compact)
  let xbox-opened = {
    let crossings = (
      D1: ((group: "p1", y: top),), D4: ((group: "p2", y: bot),), D6: ((group: "p1", y: mid),),
    )
    let g = graph.cut(g,
      left: cut-side(g, in-x, crossings, source: (<D1>,), sink: (<D6>, <D4>)),
      right: cut-side(g, out-x, crossings, sink: (<D1>,), source: (<D6>, <D4>)),
      boundary: boundary-position,
    )
    let nodes = (
      "compact-top": (pos: pos(x: pin(0), y: top)),
      "compact-bottom": (pos: pos(x: pin(0), y: bot)),
    )
    let edges = (
      compact: (statements: ("spring-length": 0.25)),
      "D1.0": mom(side: "right", label: (gap: .15, shift: -1.5)),
      "D1.1": mom(side: "left", length: 1.4, shift: -.4),
      D2: (statements: ("spring-length": .3)) + mom(label: (gap: .2)),
      D3: mom(shift: 1.),
      "D4.0": (bend: -0.18, crossing-under: <D6.0>, crossing-gap: 0.9,
        // fermion-arrow-shift: 1.15,
      ) + mom(side: "left", length: 1., shift: .5),
      "D4.1": (statements: ("spring-length": .3)) + mom(side: "left", shift: -1.3),
      D5: (statements: ("spring-length": 1.5),
        // ..mom(offset: 0.80),
      ),
      "D6.0": (bend: -0.55) + mom(side: "right"),
      "D6.1": mom(side: "left", length: 1.4, shift: -.4),
    )
    let g = graph.map(g, edge: flatten-boundary)
    graph.map(g, node: nodes, edge: edges)
  }

  let xbox-opened2 = {
    let crossings = (
      D2: ((group: "p2", y: bot),), D3: ((group: "p1", y: top),), D6: ((group: "p2", y: mid),),
    )
    let g = graph.cut(g,
      left: cut-side(g, in-x, crossings, source: (<D2>,), sink: (<D3>, <D6>)),
      right: cut-side(g, out-x, crossings, sink: (<D2>,), source: (<D3>, <D6>)),
      boundary: boundary-position,
    )
    let nodes = (c: (pos: pos(y: start(-4))), d: (pos: pos(y: start(0))))
    let edges = (
      compact: (statements: ("spring-length": 0.01)),
      D1: mom(side: "left", length: 1.2, label: (gap: .2)),
      "D2.0": mom(side: "right", length: 1., shift: .4, label: (gap: .2)),
      "D2.1": mom(side: "right", length: .7, shift: -.4,
        label: (gap: .2, shift: -.75, anchor: "south-west")),
      "D3.1": (crossing-under: <D6.1>, crossing-gap: 0.7)
        + mom(side: "left", length: 1.2, shift: -.6, label: (gap: .2)),
      D4: (statements: ("spring-length": .1)) + mom(side: "right", shift: -.4, label: (gap: .1)),
      D5: (statements: ("spring-length": .5)) + mom(side: "right", label: (gap: .05),
        // offset: 0.80,
      ),
      "D6.0": (bend: -0.55) + mom(side: "left", length: 1.5, shift: .8),
      "D6.1": mom(side: "left", shift: .8, label: (gap: .2)),
    )
    let g = graph.map(g, edge: flatten-boundary)
    graph.map(g, node: nodes, edge: edges)
  }

  let xbox-cut = {
    let crossings = (
      D1: ((group: "p1", y: mid2),),
      D2: ((group: "p2", y: bot),),
      // Along b -> c: right p2 stub, left p2 -> right p1 middle, left p1 stub.
      D6: ((group: "p2", y: mid), (group: "p1", y: top)),
    )
    let g = graph.cut(g,
      left: cut-side(g, in-x, crossings, source: (<D1>, <D2>), sink: (<D6>,)),
      right: cut-side(g, out-x, crossings, sink: (<D1>, <D2>), source: (<D6>,)),
      boundary: boundary-position,
    )
    let edges = (
      compact: (statements: ("spring-length": 1.8)),
      "D1.0": mom(side: "left", length: .8, label: (shift: 1, gap: .01, anchor: "north")),
      "D1.1": mom(side: "right", length: 1, shift: .5, label: (shift: 3.2, gap: .1)),
      "D2.0": mom(side: "right", shift: .4, label: (anchor: "south-east", gap: .2)),
      "D2.1": mom(side: "right", label: (gap: .3, shift: -.2, anchor: "west")),
      D3: (statements: ("spring-length": 1.5), crossing-under: <D6.1>, crossing-gap: 0.8,
        fermion-arrow-shift: -0.5) + mom(shift: -.6),
      D4: (statements: ("spring-length": .5)) + mom(label: (gap: .02)),
      D5: (crossing-under: <D6.1>, crossing-gap: 1.5) + mom(side: "right", shift: -1., label: (gap: .01)),
      "D6.0": mom(side: "left", shift: .5, label: (gap: .1)),
      "D6.1": (statements: ("spring-length": 3.5)) + mom(side: "left", shift: 3, label: (shift: 2.5, gap: .1)),
      "D6.2": mom(side: "right", shift: -.4, label: (gap: .1)),
    )

    let nodes = (
      a: (pos: pos(y:group("h"))),
      d: (pos: pos(y:group("h"))),
    )
    let g = graph.map(g,
      // Keep the D6 stubs free in depth, but flatten their middle endpoints.
      node: flatten-boundary,
      edge: e => (if e.origin.name != <D6> { flatten-boundary(e) }) + mom(label: (gap: .4)),
    )
    graph.map(g, node: nodes, edge: edges)
  }

  let diagrams = $
    #diagram(xbox, cut-x: -1)+
    #diagram(xbox-opened, cut-x: 0,cut-y: 1)+
    #diagram(xbox-opened2, cut-x: -2.1,cut-y: 1)+
    #diagram(xbox-cut,cut-y: 0.5) = op("disc")_(p_1^2)  op("disc")_(p_2^2) integral (dif ^d k  )/(2 pi )^d (N^frak(q q')_times.square delta^+_(q^2)(p_(12)-k) delta^+_0(k))/(  p_1^2 p_2^2 (k-p_2)^2 (k-p_1)^2)
  $
  diagrams
}
