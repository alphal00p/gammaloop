#set page(height: auto, width: 179mm, margin: 5mm)

#import "../src/lib.typ": draw, graph, layout, layouts, subgraph
#import graph: *
#import "map-style.typ" as feynman
#import feynman: momentum as mom
#import "@preview/cetz:0.5.1" as cetz

// Native font and stroke sizes are independent of the graph coordinate unit.
#let graph-style = feynman.graph-style(unit: 2.6mm, line-width: 0.5pt)
#let diagram-style = (
  font-size: 6pt,
  cut-line-width: 0.8pt,

  endpoint-box: (
    padding: (x: 0.5, y: 0.4),
    radius: 0.25,
    fills: (p1: blue.lighten(85%), p2: orange.lighten(85%)),
  ),
  padding: 0.3,
)
#set text(size: diagram-style.font-size)

#let draw-initials = false;

// One stroke per opening: original, p1 replaced, p2 replaced, both replaced.
#let initial-cut-styles = (
  (paint: blue.darken(45%), dash: "dotted",thickness:1pt),
  (paint: blue.darken(15%), dash: "dashed",thickness:1pt),
  (paint: blue.lighten(10%), dash: (4pt, 1.5pt, 1pt, 1.5pt),thickness:1pt),
  (paint: blue.lighten(30%), dash: (4pt, 1.5pt, 1pt, 1.5pt, 1pt, 1.5pt),thickness:1pt),
).map(style => (thickness: diagram-style.cut-line-width) + style)

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
#let diagram(g, options: base-layout, cut-x: auto, cut-y: auto, initial-cut: none, draw-after: none) = context {
  draw(
    layout(
      graph.style(g, ..graph-style),
      ..options,
    ),
    ..feynman.draw-style,
    padding: diagram-style.padding,
    draw-after: (g, bounds) => {
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

      // Exclude the invisible compactification spring and its anchor nodes
      // when locating the cut; the renderer supplies the visible drawing bounds.
      let nodes = graph.nodes(g).filter(n => n.boundary == none and (
        n.data == none or not n.data.at("hidden", default: false)
      ))
      let xs = nodes.map(n => n.pos.x)
      let x = if cut-x == auto {
        calc.min(..xs) + 0.6 * (calc.max(..xs) - calc.min(..xs))
      } else { cut-x }
      let y = if cut-y == auto { 1.5 } else { cut-y }
      bounds = (
        left: bounds.left - box-style.padding.x,
        right: bounds.right + box-style.padding.x,
        top: bounds.top + y,
        bottom: bounds.bottom - y,
        width: bounds.width + 2 * box-style.padding.x,
        height: bounds.height + 2 * y,
      )
      cetz.draw.line(
        (x, bounds.top), (x, bounds.bottom),
        stroke: (paint: red, thickness: diagram-style.cut-line-width, dash: "dashed"),
      )
      if initial-cut != none and draw-initials {
        for side in ("left", "right") {
          let side-xs = endpoints.filter(b => b.side == side).map(b => b.pos.x)
          if side-xs.len() > 0 {
            let x = side-xs.first()
            cetz.draw.line((x, bounds.top), (x, bounds.bottom), stroke: initial-cut-styles.at(initial-cut))
          }
        }
      }
      if type(draw-after) == function { draw-after(g, bounds) } else { draw-after }
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
      "D3.0": (spring-length: .3) + mom(label: (gap: .4)),
      "D3.1": (spring-length: .3, reverse: true)
        + mom(side: "left", label: (gap: .4)),
      "D4.0": (spring-length: .3) + mom(side: "right", label: (gap: .2)),
      "D4.1": (spring-length: .3) + mom(label: (gap: .2)),
      D5: (spring-length: .2,crossing-under: <D6>, crossing-gap: 0.9) + mom(
        side: "right", shift: 1.5, label: (gap: .05),
        // offset: 0.80,
      ),
      D6:(spring-length: .2)+ mom(side: "left", shift: 1.5, label: (gap: .2)),
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
      compact: (spring-length: 0.25),
      "D1.0": mom(side: "right", label: (gap: .15, shift: -1.5)),
      "D1.1": (spring-length: .3)+mom(side: "left", length: 1.4, shift: -.4),
      D2: (spring-length: .3) + mom(label: (gap: .2)),
      D3: mom(shift: 1.),
      "D4.0": (spring-length: .3, bend: -0.18, crossing-under: <D6.0>, crossing-gap: 0.9,
        // fermion-arrow-shift: 1.15,
      ) + mom(side: "left", length: 1., shift: .5),
      "D4.1": (spring-length: .3) + mom(side: "left", shift: -1.3),
      D5: (spring-length: .5,),
      "D6.0": (spring-length: .5)+(bend: -1.55) + mom(side: "right"),
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
      compact: (spring-length: 0.01),
      D1: mom(side: "left", length: 1.2, label: (gap: .2)),
      "D2.0": mom(side: "right", length: 1., shift: .4, label: (gap: .2)),
      "D2.1": mom(side: "right", length: .7, shift: -.4,
        label: (gap: .2, shift: -.75, anchor: "south-west")),
      "D3.0": (spring-length: .2),
      "D3.1": (crossing-under: <D6.1>, crossing-gap: 0.7)
        + mom(side: "left", length: 1.2, shift: -.6, label: (gap: .2)),
      D4: (spring-length: .1) + mom(side: "right", shift: -.4, label: (gap: .1)),
      D5: (spring-length: .5) + mom(side: "right", label: (gap: .05),
        // offset: 0.80,
      ),
      "D6.0": (spring-length: .3)+(bend: -0.55) + mom(side: "left", length: 1.5, shift: .8),
      "D6.1": (spring-length: .3)+mom(side: "left", shift: .8, label: (gap: .2)),
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
      compact: (spring-length: 1.8),
      "D1.0": mom(side: "left", length: .8, label: (shift: 1, gap: .01, anchor: "north")),
      "D1.1": mom(side: "right", length: 1, shift: .5, label: (shift: 3.2, gap: .1)),
      "D2.0": mom(side: "right", shift: .4, label: (anchor: "south-east", gap: .2)),
      "D2.1": mom(side: "right", label: (gap: .3, shift: -.2, anchor: "west")),
      D3: (spring-length: 1.5, crossing-under: <D6.1>, crossing-gap: 0.8,
        fermion-arrow-shift: -0.5) + mom(shift: -.6),
      D4: (spring-length: .5) + mom(label: (gap: .02)),
      D5: (crossing-under: <D6.1>, crossing-gap: 1.5,spring-length: .2) + mom(side: "right", shift: -1., label: (gap: .01)),
      "D6.0": mom(side: "left", shift: .5, label: (gap: .1)),
      "D6.1": (spring-length: 2.) + mom(side: "left", shift: 3, label: (shift: 2.5, gap: .1)),
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
    #diagram(xbox, cut-x: -1, initial-cut: 0,
      draw-after: (g, bounds) => {
      if draw-initials{
      let nodes = graph.nodes(g)
      // The two single-replacement cuts cross the other external leg;
      // the matching pair opens both, with each branch reaching the diagram boundary.
      for (name, side, style, both) in (
        (<c>, 1, 1, false), (<b>, -1, 2, false),
        (<c>, 1, 3, true), (<b>, -1, 3, true),
      ) {
        let p = nodes.find(n => n.name == name).pos
        let (near, far) = if side > 0 { (bounds.top, bounds.bottom) } else { (bounds.bottom, bounds.top) }
        let outer = if side > 0 { bounds.right } else { bounds.left }
        // Meet the top/bottom normally; the paired cut also exits the sides horizontally.
        let points = if both {
          ((p.x - side * .8, near), (outer, p.y - side * 1.7),
            (p.x - side * .8, p.y - side * .6), (p.x + side * .1, p.y - side * 1.7))
        } else {
          ((p.x - side * 1.25, near), (outer - side, far),
            (p.x - side * 1.25, p.y - side * 2), (outer - side, p.y - side * .8))
        }
        // Both branches cross D6 in the same direction, retaining its middle segment.
        if side < 0 { points = (points.at(1), points.at(0), points.at(3), points.at(2)) }
        cetz.draw.bezier(
          ..points,
          stroke: initial-cut-styles.at(style),
        )
      }
    }
    }
    )+
    #diagram(xbox-opened, cut-x: 1.25,cut-y: 1, initial-cut: 1)#h(-2mm)+#h(-2mm)
    #diagram(xbox-opened2, cut-x: -2.5,cut-y: 1, initial-cut: 2)#h(-2mm)+
    #diagram(xbox-cut,cut-y: 0.5, initial-cut: 3) = op("disc")_(p_1^2)  op("disc")_(p_2^2) integral (dif ^d k  )/(2 pi )^d (N^(q overline(q))_times.square delta^+_(q^2)(p_(12)-k) delta^+_0(k))/(  p_1^2 p_2^2 (k-p_2)^2 (k-p_1)^2)
  $
  diagrams
}





#let better-layout = layouts.options(
  spring: (strength: 40, length: 0.2),
  repulsion: (
    strength: 10,
    centering: 0.005,
    edge-node: .39,
    edge-edge: .75,
    dangling: 5.5,
    dangling-centroid: 30,
  ),
  constraints: (side-strength: 10),
  labels: (steps: 100, distance: 0.9, spring: 12, repulsion: 8),
  solver: (
    steps: 50,
    epochs: 50,
    step: .81,
    max-movement: 0.2,
    cooling: 0.85,
    depth-scale: 1,
    flattening-end: 0.5,
  ),
)
#pagebreak()
#grid(columns:2,align:center+horizon, diagram(graph.build(default-edge-data: (particle: "d",spring-length:.3) + mom(offset: .4), {
  node(<a>)
  node(<b>)
  node(<c>)
  node(<d>)
  edge(<D1.1>, sink(<a>), momentum: [$p_1$],pos:pos(x:in-x,y:top),..mom(side: "right", label: ( gap: .01,)))
  edge(<D2>, source(<a>), sink(<b>), momentum: [$p_1+k$],..mom(side: "left",shift:-0.4,length: .8, label: (shift: -0.3, gap: .2,)))
  
  edge(<D3.1>, sink(<b>),orientation: "reversed", momentum: [$p_2-k$],pos:pos(x:in-x,y:bot),..mom(side: "left",shift: -0.4, label: ( shift: -0.7, gap: .01,)))
  edge(<D4>, source(<b>), sink(<c>), momentum: [$p_(12)$], particle: "a",spring-length:0.1,..mom(side: "right", label: ( gap: .1,)))
  edge(<D3.2>, source(<c>), orientation: "reversed", momentum: [$p_2-k$],pos:pos(x:out-x,y:bot),..mom(side: "left", label: ( gap: .2,shift:0.4)))
  edge(<D5>, source(<c>), sink(<d>), momentum: [$p_1+k$],..mom(side: "left",shift:-0.4,length: .8, label: (shift: -1, gap: .1,)))
  edge(<D1.2>, source(<d>), momentum: [$p_1$],pos:pos(x:out-x,y:top),..mom(side: "left", label: ( gap: .1,)))
  edge(<D6.1>, sink(<a>), momentum: [$k$],particle:"g",pos:pos(x:in-x,y:mid),..mom(side: "right", label: ( gap: .1,)))
  edge(<D6.2>, source(<d>), momentum: [$k$],particle:"g",pos:pos(x:out-x,y:mid),..mom(side: "left", label: ( gap: .2,)))
}

// ,compact
),options: better-layout,cut-y: 0.5,cut-x:0.3)


,diagram(graph.build(default-edge-data: (particle: "d") + mom(offset: .4), {
  node(<a>)
  node(<b>,pos:pos(y:mid))
  node(<c>,pos:pos(y:mid))
  node(<d>)
  edge(<D1.1>, sink(<a>), momentum: [$p_1$],pos:pos(x:in-x,y:top),..mom(side: "left", length: .8, label: (shift: 1, gap: .01, anchor: "south")))
  edge(<D2>, source(<a>), sink(<b>), momentum: [$p_1-k$],..mom(side: "right", length: .8, label: ( gap: .1)))
  
  edge(<D3.1>, sink(<b>),orientation: "reversed", momentum: [$p_2-k$],pos:pos(x:in-x,y:bot),spring-length:1.5,..mom(side: "right", length: .8, label: ( gap: .1)),crossing-under: <D6.2>, crossing-gap: .7)
  edge(<D4>, source(<b>), sink(<c>), momentum: [$p_(12) - 2k$], particle: "a",spring-length:1.5,..mom(side: "right", length: .8, label: (shift: 0.25, gap: .01, anchor: "north")))
  edge(<D3.2>, sink(<c>), momentum: [$p_1$],pos:pos(x:out-x,y:top),spring-length:1.5,crossing-under: <D6.1>, crossing-gap: .7)
  edge(<D5>, source(<c>), sink(<d>), momentum: [$p_2-2k$],..mom(side: "left", length: .8, label: (shift: -0.25, gap: .2)))
  edge(<D1.2>, source(<d>), momentum: [$p_2-k$],pos:pos(x:out-x,y:bot),..mom(side: "left",shift:1, length: .8, label: ( shift: 2, gap: .2,anchor:"south-west")))
  edge(<D6.1>, source(<a>), momentum: [$k$],particle:"g",pos:pos(x:out-x,y:mid), bend:.44,spring-length:2.5,..mom(side: "left",shift:2, length: .8, label: ( gap: .2)))
  edge(<D6.2>, sink(<d>), bend:-.54,momentum: [$k$],particle:"g",pos:pos(x:in-x,y:mid),spring-length:2.5,..mom(side: "right",shift:-2, length: .8, label: ( gap: .2)))
}


),options: better-layout,cut-x: 0.3,cut-y:0.5)
)


#diagram(graph.build(default-edge-data: (particle: "d",spring-length:.3) + mom(offset: .4), {
  node(<a>)
  node(<b>)
  node(<c>)
  node(<d>)
  edge(<D1.1>, sink(<a>),momentum:[],pos:pos(x:in-x,y:bot,z:pin(0)),..mom(side: "right", label: ( gap: .01,)))
  edge(<D2>, source(<a>), sink(<b>), momentum: [$p_1+k$],..mom(side: "left",shift:-0.4,length: .8, label: (shift: -0.3, gap: .2,)))
  
  edge(<D3.1>, sink(<b>),orientation: "reversed", momentum: [$p_2-k$],pos:pos(x:in-x,y:top,z:pin(0)),..mom(side: "left",shift: -0.4, label: ( shift: -0.7, gap: .01,)))
  edge(<D4>, source(<b>), sink(<c>), momentum: [$p_(12)$], particle: "a",spring-length:0.1,..mom(side: "right", label: ( gap: .1,)))
  edge(<D3.2>, source(<c>), orientation: "reversed", momentum: [$p_2-k$],pos:pos(x:out-x,y:bot,z:pin(0)),..mom(side: "left", label: ( gap: .2,shift:0.4)))
  edge(<D5>, source(<c>), sink(<d>), momentum: [$p_1+k$],..mom(side: "left",shift:-0.4,length: .8, label: (shift: -1, gap: .1,)))
  edge(<D1.2>, source(<d>), momentum: [$p_1$],pos:pos(x:out-x,y:top,z:pin(0)),..mom(side: "left", label: ( gap: .1,)))
  edge(<D6.1>, sink(<a>), momentum: [$k$],particle:"g",pos:pos(x:in-x,y:mid,z:pin(0)),..mom(side: "right", label: ( gap: .1,)))
  edge(<D6.2>, source(<d>), momentum: [$k$],particle:"g",pos:pos(x:out-x,y:mid,z:pin(0)),..mom(side: "left", label: ( gap: .2,)))
}

// ,compact
),options: better-layout,cut-y: 0.5,cut-x:0.3)