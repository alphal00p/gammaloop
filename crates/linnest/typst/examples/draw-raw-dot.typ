#import "../src/lib.typ": draw, graph, layout, physics, subgraph
#import graph: *
#import physics: *


#let fermion = (
  source: source-stroke(c: blue, thickness: massless),
  sink: sink-stroke(c: blue, thickness: massless),
  fermion-arrow: true,
)

#let massive-fermion = (
  source: source-stroke(c: blue, thickness: massive),
  sink: sink-stroke(c: blue, thickness: massive),
  fermion-arrow: true,
)

#let scalar = (
  source: source-stroke(c: black, thickness: massive, dash: dashed),
  sink: sink-stroke(c: black, thickness: massive, dash: dashed),
)

#let gluon = (
  source: source-stroke(c: black, thickness: massless) + coil,
  sink: sink-stroke(c: black, thickness: massless) + coil,
  label: mi(`{g}`),
)

#let photon = (
  source: source-stroke(c: black, thickness: massless) + wave,
  sink: sink-stroke(c: black, thickness: massless) + wave,
)


#let map = (
  "a": photon + (label: [$gamma$]),
  "mu+": fermion + (label: [$mu^+$]),
  "mu-": fermion + (label: [$mu^-$]),
  "e+": fermion + (label: [$e^+$]),
  "e-": fermion + (label: [$e^-$]),
  "q": fermion + (label: [q]),
  "d": fermion + (label: [d]),
  "t": massive-fermion + (label: [$t$]),
  "photon": photon,
  "g": gluon,
  "gluon": gluon,
  "fermion": fermion + (label: [$f$]),
  "scalar": scalar + (label: [$s$]),
  "ghG": scalar + (label: [$s$]),
)
#let (source-style, sink-style, edge-label, ..callbacks) = physics.style(
  map: map,
  momentum-arrows: false,
  show-particle: false,
  show-edge-index: false,
)
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

#let _particle-label(edge, edge-style-map) = {
  let scope = edge.fields + (edge: edge.edge)
  let label-value = physics.edge-entry(scope, map: edge-style-map, default: (label: none)).at("label", default: none)
  let particle-label = physics.label-content(label-value, scope, map: edge-style-map, scope: scope)
  if particle-label == none {
    return none
  }
  particle-label + [$(p_#edge.edge)$]
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

#let autogen-external-edge-fields(g, edge-style-map: map, y-scale: 10) = {
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
    let label = _particle-label(edge, edge-style-map)
    if label != none {
      patch.insert("label", label)
    }
    patch
  })
}

#let draw_dot(doc) = {
  show raw: it => if it.at("lang") == "dot" {
    for g in parse(it.text, eval-edge-fields: "label") {
      let g = autogen-external-edge-fields(g)
      let layed-out = layout(
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
          #edges

        ]
          + draw(
            layed-out,
            source-style: source-style,
            sink-style: sink-style,
            edge-label: edge-label,
            edge-label-style: edge-label-style,
            subgraph: subgraph.compass(layed-out, "s"),
          )
      )
    }
  }
  doc
}
