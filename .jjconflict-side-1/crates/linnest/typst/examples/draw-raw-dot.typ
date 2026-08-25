#import "../src/lib.typ": draw, graph, layout
#import graph: *


#let source-style(edge) = (
  stroke: (
    paint: if edge.ext { rgb("#596579") } else { rgb("#315f9f") },
    thickness: if edge.ext { 0.7pt } else { 0.9pt },
    cap: "round",
  ),
)
#let sink-style(edge) = (
  stroke: (
    paint: if edge.ext { rgb("#7a879b") } else { rgb("#79a0d0") },
    thickness: if edge.ext { 0.7pt } else { 0.9pt },
    cap: "round",
  ),
)
#let edge-label(edge) = [link #edge.eid]
#let edge-label-style(edge) = {
  let raw-edge = edge.at("edge", default: (:))
  let statements = raw-edge.at("statements", default: (:))
  let anchor = statements.at("label-anchor", default: edge.at("label-anchor", default: none))
  if type(anchor) == str {
    (anchor: anchor)
  } else {
    (:)
  }
}

#let _rank(ids, id) = {
  for (rank, value) in ids.enumerate() {
    if value == id {
      return rank
    }
  }
  none
}

#let _centered-y(ids, id, y-scale) = {
  let rank = _rank(ids, id)
  if rank == none {
    return 0
  }
  ((ids.len() - 1) / 2 - rank) * y-scale
}

#let autogen-boundary-edge-fields(g, y-scale: 10) = {
  let left = ()
  let right = ()
  for edge in graph.edges(g) {
    if edge.source == none and edge.sink != none {
      left.push(edge.edge)
    } else if edge.source != none and edge.sink == none {
      right.push(edge.edge)
    }
  }

  graph.map(g, edge: edge => {
    let side = if edge.source == none and edge.sink != none {
      "left"
    } else if edge.source != none and edge.sink == none {
      "right"
    } else {
      none
    }
    if side == none {
      return none
    }

    let ids = if side == "left" { left } else { right }
    let x = if side == "left" {
      graph.group("left", side: "-")
    } else {
      graph.group("right", side: "+")
    }
    let anchor = if side == "left" { "east" } else { "west" }
    let patch = (
      pos: graph.pos(x: x, y: graph.start(_centered-y(ids, edge.edge, y-scale))),
      "label-anchor": anchor,
    )
    patch
  })
}

#let draw_dot(doc) = {
  show raw: it => if it.at("lang") == "dot" {
    for g in parse(it.text) {
      let g = autogen-boundary-edge-fields(g)
      let g = graph.style(
        g,
        node-label: node => [item #node.vid],
        node-label-style: (padding: 0.12),
        node-style: (
          fill: rgb("#f3f6fa"),
          stroke: rgb("#34445c") + 0.8pt,
        ),
        edge-label: edge-label,
        edge-label-style: edge-label-style,
      )
      let laid-out = layout(
        g,
        layout-algo: "force",
        epochs: 30,
        steps: 40,
        k-spring: 4.5,
        eps: 1e-7,
        step: .6,
        gamma-dangling: 2.3,
        label-length-scale: 1.2,
        label-steps: 100,
        directional-force: 4.5,
        label-layout: "dangling-tangent",
      )


      (
        [#align(center, info(g).name)
          #edges(g).len() links

        ]
          + draw(
            laid-out,
            source-style: source-style,
            sink-style: sink-style,
          )
      )
    }
  }
  doc
}
