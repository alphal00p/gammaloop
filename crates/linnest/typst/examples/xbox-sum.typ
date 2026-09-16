#import "xbox-common.typ": *
#set page(height: auto, width: 179mm, margin: 2mm)
#set text(size: diagram-style.font-size)
#show math.equation: set text(size: math-font-size)

// The physical K4 is shared; only cuts and drawing patches vary between views.
#let master = {
  vertices
  edge(<D1>, source(<c>), sink(<d>), momentum: [$k-p_1$])
  edge(
    <D2>,
    source(<a>),
    sink(<b>),
    momentum: [$k-p_2$],
    orientation: "reversed",
  )
  edge(<D3>, source(<c>), sink(<a>), momentum: [$p_1$], orientation: "reversed")
  edge(<D4>, source(<d>), sink(<b>), momentum: [$p_2$])
  edge(
    <D5>,
    source(<a>),
    sink(<d>),
    momentum: [$p_(12) - k$],
    particle: "a",
    show-momentum: false,
  )
  edge(<D6>, source(<b>), sink(<c>), momentum: [$k$], particle: "g")
}

// Each entry describes one crossing along the underlying source -> sink flow.
#let cut-side(g, x, crossings, ..selection) = subgraph.with-data(
  g,
  subgraph.select(g, ..selection),
  data: (x: x),
  hedge: h => {
    let points = crossings.at(str(h.edge-name))
    (winding: points.len(), crossings: points)
  },
)
#let boundary-position(item) = {
  let b = item.boundary
  (pos: pos(x: b.cut-data.x, y: b.data.crossings.at(b.crossing).y))
}
#let flatten-boundary = item => if item.boundary != none {
  (pos: pos(z: pin(0)))
}

#{
  let xbox = {
    let g = graph.build(default-edge-data: edge-data, master)
    let crossings = (D3: ((group: "p1", y: top),), D4: ((group: "p2", y: bot),))
    let g = graph.cut(
      g,
      left: cut-side(g, in-x, crossings, sink: (<D3>, <D4>)),
      right: cut-side(g, out-x, crossings, source: (<D3>, <D4>)),
      boundary: boundary-position,
    )
    let nodes = (
      a: (pos: pos(x: start(-2), y: top)),
      b: (pos: pos(x: start(-2), y: bot)),
      c: (pos: pos(x: start(2), y: top)),
      d: (pos: pos(x: start(2), y: bot)),
    )
    let edges = (
      D1: (spring-length: 1.3,show-momentum: false),
      D2: (spring-length: 1.3,show-momentum: false),
      "D3.0": (spring-length: .3) + mom(label: (gap: .4)),
      "D3.1": (spring-length: .3, reverse: true)
        + mom(side: "left", label: (gap: .4)),
      "D4.0": (spring-length: .3) + mom(side: "right", label: (gap: .2)),
      "D4.1": (spring-length: .3) + mom(label: (gap: .2)),
      D5: (
        spring-length: 1.3,
        crossing-under: <D6>,
        crossing-gap: 0.9,
        // momentum-arrow-offset: 0.80,
      ),
      D6: (spring-length: 1.3) + mom(side: "left", shift: 1.5, label: (gap: .2)),
    )
    graph.map(g, node: nodes, edge: edges)
  }

  // Pull the external rows together without drawing another propagator.
  let g = graph.build(default-edge-data: edge-data, master, compact)
  let xbox-opened = {
    let crossings = (
      D1: ((group: "p1", y: top),),
      D4: ((group: "p2", y: bot),),
      D6: ((group: "p1", y: mid),),
    )
    let g = graph.cut(
      g,
      left: cut-side(g, in-x, crossings, source: (<D1>,), sink: (<D6>, <D4>)),
      right: cut-side(g, out-x, crossings, sink: (<D1>,), source: (<D6>, <D4>)),
      boundary: boundary-position,
    )
    let nodes = (
      "compact-top": (pos: pos(x: pin(0), y: top)),
      "compact-bottom": (pos: pos(x: pin(0), y: bot)),
    )
    let edges = (
      compact: (spring-length: 0.3),
      "D1.0": (spring-length: 1.3)+mom(side: "right", label: (gap: .15, shift: -1.5)),
      "D1.1": (spring-length: 2) + mom(side: "right", label: (gap: .01, shift: -1.5, anchor: "south-east")),
      D2: (spring-length: .6, show-momentum: false),
      D3: (spring-length: .6, show-momentum: false),
      "D4.0": (
        spring-length: 2.,
        bend: -0.18,
        crossing-under: <D6.0>,
        crossing-gap: 0.9,
        fermion-arrow-shift: -.7,
      )
        + mom(side: "right", shift: .7, label: (gap: .015, shift: .9)),
      "D4.1": (spring-length: 0.1) + mom(side: "right", shift: -.5, label: (gap: .01)),
      D5: (spring-length: .5),
      "D6.0": (spring-length: 1.5, bend: -1.55) + mom(side: "right", shift:-.2,label: (gap: .1)),
      "D6.1": (spring-length: 1.4) + mom(side: "right", shift: -.2, label: (gap: .05, shift: .2)),
    )
    let g = graph.map(g, edge: flatten-boundary)
    graph.map(g, node: nodes, edge: edges)
  }

  let xbox-opened2 = {
    let crossings = (
      D2: ((group: "p2", y: bot),),
      D3: ((group: "p1", y: top),),
      D6: ((group: "p2", y: mid),),
    )
    let g = graph.cut(
      g,
      left: cut-side(g, in-x, crossings, source: (<D2>,), sink: (<D3>, <D6>)),
      right: cut-side(g, out-x, crossings, sink: (<D2>,), source: (<D3>, <D6>)),
      boundary: boundary-position,
    )
    let nodes = (c: (pos: pos(y: start(-4))), d: (pos: pos(y: start(0))))
    let edges = (
      compact: (spring-length: 1.7),
      D1: (spring-length: 1)+(show-momentum: false),
      "D2.0": (spring-length: 1.2)+mom(side: "left",  shift: .3, label: (gap: .1)),
      "D2.1": (spring-length: 1)+mom(side: "left",  shift: .3, label: (
        gap: 0.01,
      )),
      "D3.0": (spring-length: 1.2)
        + mom(shift: .6, ),
      "D3.1": (spring-length: 1.2,crossing-under: <D6.1>, crossing-gap: 0.7, fermion-arrow-shift: -.3)
        + mom(side: "left", shift: -.5, label: (gap: .2,shift:-.7)),
      D4: (spring-length: .1, show-momentum: false),
      D5: (spring-length: .1),
      "D6.0": (spring-length: 1)
        + mom(side: "left", label: (gap: .1)),
      "D6.1": (spring-length: 1.5)
        + mom(side: "left", shift: .8, label: (gap: .2)),
    )
    let g = graph.map(g, edge: flatten-boundary)
    graph.map(g, node: nodes, edge: edges)
  }
let g = graph.build(default-edge-data: edge-data, master)
  let xbox-cut = {
    let crossings = (
      D1: ((group: "p1", y: top),),
      D2: ((group: "p2", y: bot),),
      // Along b -> c: right p2 stub, left p2 -> right p1 middle, left p1 stub.
      D6: ((group: "p2", y: mid), (group: "p1", y: mid2)),
    )
    let g = graph.cut(
      g,
      left: cut-side(g, in-x, crossings, source: (<D1>, <D2>), sink: (<D6>,)),
      right: cut-side(g, out-x, crossings, sink: (<D1>, <D2>), source: (<D6>,)),
      boundary: boundary-position
    )
    let edges = (
      // compact: (spring-length: 10),
      "D1.0": (spring-length: 1.4)+mom(side: "right", length: .8, label: (gap: .01)),
      "D1.1": (spring-length: 1.6)+mom(side: "right", length: 1, shift: -.5, label: (gap: .1)),
      "D2.0": (spring-length: 2.3,  crossing-under: <D6.1>,
      crossing-gap: 0.8,)+mom(side: "left", shift: .8, label: (
        gap: .1,
        shift: .8,
      )),
      "D2.1": (spring-length: 2)+mom(side: "left", shift: -.3, label: (gap: .01, shift: -.2)),
      D3: (
        spring-length: 1.,
        show-momentum: false,
      ),
      D4: (spring-length: 1.2, show-momentum: false,  crossing-under: <D6.1>,
      crossing-gap: 0.8,  fermion-arrow-shift: 0.3),
      D5: (spring-length: .1),
      "D6.0": (spring-length: 2)+mom(side: "right", label: (shift: 0.2, gap: .15)),
      "D6.1": (
        spring-length: 6,
        shift: (0, -0.4),
        edge-style: (source-anchor: "east", sink-anchor: "west"),
        
      )
        + mom(side: "right",shift:0.2,label: (gap: .1)),
      "D6.2": mom(side: "right", label: (shift:1,gap: .1)),
    )
    let nodes = (a: (pos: pos(y: group("h"))), d: (pos: pos(y: group("h"))))
    let g = graph.map(
      g,
      // Keep the D6 stubs free in depth, but flatten their middle endpoints.
      node: flatten-boundary,
      edge: e => if e.origin.name != <D6> { flatten-boundary(e) },
    )
    graph.map(g, node: nodes, edge: edges)
  }

  $
    #diagram(xbox, cut-x: -1,cut-y: -.8, initial-cut: 0, draw-after: (g, bounds) => {
      // if draw-initials {
      //   let nodes = graph.nodes(g)
      //   // The two single-replacement cuts cross the other external leg;
      //   // the matching pair opens both, with each branch reaching the diagram boundary.
      //   for (name, side, style, both) in (
      //     (<c>, 1, 1, false),
      //     (<b>, -1, 2, false),
      //     (<c>, 1, 3, true),
      //     (<b>, -1, 3, true),
      //   ) {
      //     let p = nodes.find(n => n.name == name).pos
      //     let (near, far) = if side > 0 { (bounds.top, bounds.bottom) } else { (bounds.bottom, bounds.top) }
      //     let outer = if side > 0 { bounds.right } else { bounds.left }
      //     // Meet the top/bottom normally; the paired cut also exits the sides horizontally.
      //     let points = if both {
      //       ((p.x - side * .8, near), (outer, p.y - side * 1.7), (p.x - side * .8, p.y - side * .6), (p.x + side * .1, p.y - side * 1.7))
      //     } else {
      //       ((p.x - side * 1.25, near), (outer - side, far), (p.x - side * 1.25, p.y - side * 2), (outer - side, p.y - side * .8))
      //     }
      //     // Both branches cross D6 in the same direction, retaining its middle segment.
      //     if side < 0 { points = (points.at(1), points.at(0), points.at(3), points.at(2)) }
      //     cetz.draw.bezier(..points, stroke: initial-cut-styles.at(style))
      //   }
      // }
    })+
    #diagram(xbox-opened, cut-x: -.4, cut-y: -0.8, initial-cut: 1)+
    #diagram(xbox-opened2, cut-x: -1.5, cut-y: -0.8, initial-cut: 2)+
    #diagram(xbox-cut, cut-y: -0.8,cut-x:-0.3, initial-cut: 3) = op("disc")_(p_1^2) op("disc")_(p_2^2) integral (dif^d k)/(2 pi)^d (N^(q overline(q))_times.square delta^+_(q^2)(p_(12)-k) delta^+_0(k))/(p_1^2 p_2^2 (k-p_2)^2 (k-p_1)^2)
  $
}
