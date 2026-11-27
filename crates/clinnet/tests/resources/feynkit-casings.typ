#set page(width: auto, height: auto, margin: (x: 2mm, y: 2mm), fill: none)
#set text(size: 9pt)
#import "crates/linnest/typst/src/graph.typ" as graph
#import "crates/linnest/typst/src/subgraph.typ" as subgraph
#import "crates/linnest/typst/src/render/layout.typ" as renderer
#import "assets/embedded/drawing/templates/layout-core.typ" as physics-layout
#import "assets/embedded/drawing/templates/physics-edge-style.typ" as physics
#import physics: palette, massive, massless, dashed, dotted, source-stroke, sink-stroke, fermion-flow, wave, coil, zigzag
#import graph: build, edge, node, sink, source
#set text(fill: palette.ink)

#let particle-map = (
  "phi": (source:source-stroke(c: palette.ink, thickness: massive, dash: dashed), sink:sink-stroke(c: palette.ink, thickness: massive, dash: dashed), label:$ phi $),
)

#context {
  let raw = build({
    node(<v0>, id: 0, inspection: (node: 0, edges: (0,1,2,3,)), feynkit-name: "a", interaction: "V_SCALAR", numerator: "1")
    node(<v1>, id: 1, inspection: (node: 1, edges: (0,1,2,3,)), feynkit-name: "b", interaction: "V_SCALAR", numerator: "1")
    edge(source(<v0>), <e0>, sink(<v1>), id: 0, orientation: "undirected", particle: "phi", pdg: 9000001, directed: false, numerator: "1", inspection: (edge: 0, source: 0, sink: 1), momentum: "k_0", momentum-signature: (loops: (1,0,0,), external: ()))
    edge(source(<v0>), <e1>, sink(<v1>), id: 1, orientation: "undirected", particle: "phi", pdg: 9000001, directed: false, numerator: "1", inspection: (edge: 1, source: 0, sink: 1), momentum: "k_1", momentum-signature: (loops: (0,1,0,), external: ()))
    edge(source(<v0>), <e2>, sink(<v1>), id: 2, orientation: "undirected", particle: "phi", pdg: 9000001, directed: false, numerator: "1", inspection: (edge: 2, source: 0, sink: 1), momentum: "k_2", momentum-signature: (loops: (0,0,1,), external: ()))
    edge(source(<v0>), <e3>, sink(<v1>), id: 3, orientation: "undirected", particle: "phi", pdg: 9000001, directed: false, numerator: "1", inspection: (edge: 3, source: 0, sink: 1), momentum: "-k_0-k_1-k_2", momentum-signature: (loops: (-1,-1,-1,), external: ()))
  }, name: "basketball", data: (symmetry-factor: 1, overall-factor: "1", numerator: "1", loop-count: 3))
  raw = graph.map(raw, edge: e => (inspection: e.data.inspection + (
    source-hedge: if e.source == none { none } else { e.source.hedge },
    sink-hedge: if e.sink == none { none } else { e.sink.hedge },
  )))
  let outside-stroke = (paint: rgb("#777777").transparentize(55%), thickness: 0.6pt, dash: "dotted")
  let highlight-0 = subgraph.select(raw, edges: (0, 1, 2, 3, ))
  let style-0 = edge => (stroke: physics.source-style(edge, map: particle-map).stroke + (paint: rgb("#3b82f6"), thickness: 1.45pt,))
  let node-0 = (nodes: (0, 1, ), style: (:))
  let highlight-1 = subgraph.select(raw, edges: (0, 1, 2, ))
  let style-1 = edge => (stroke: physics.source-style(edge, map: particle-map).stroke + (paint: rgb("#8b5cf6"), thickness: 1.55pt,))
  let node-1 = (nodes: (0, 1, ), style: (:))
  let highlight-2 = subgraph.select(raw, edges: (0, 1, ))
  let style-2 = edge => (stroke: physics.source-style(edge, map: particle-map).stroke + (paint: rgb("#e4504f"), thickness: 1.65pt,))
  let node-2 = (nodes: (0, 1, ), style: (:))

  let config = (:)
  let options = config.at("options", default: (:))
  let options = (label-fill: palette.ink) + if options == none { (:) } else { options }
  if not ("mode", "amplitude-mode", "cross-section-mode").any(key => options.keys().contains(key)) {
    options += (mode: "amplitude",)
  }
  // Open only the sewn initial-state connections; final-state cut edges remain paired.
  if options.at("split-initial-state", default: true) {
    let initial = graph.edges(raw).filter(e => e.source != none and e.sink != none and e.data.at("is_cut", default: none) != none).map(e => e.edge)
    if initial.len() > 0 {
      raw = graph.cut(raw, left: subgraph.select(raw, sink: initial), right: subgraph.select(raw, source: initial))
      highlight-0 = subgraph.select(raw, hedges: graph.edges(raw).map(e => (e.source, e.sink).filter(h => h != none and subgraph.contains(highlight-0, h.origin)).map(h => h.hedge)).flatten())
      highlight-1 = subgraph.select(raw, hedges: graph.edges(raw).map(e => (e.source, e.sink).filter(h => h != none and subgraph.contains(highlight-1, h.origin)).map(h => h.hedge)).flatten())
      highlight-2 = subgraph.select(raw, hedges: graph.edges(raw).map(e => (e.source, e.sink).filter(h => h != none and subgraph.contains(highlight-2, h.origin)).map(h => h.hedge)).flatten())
    }
  }
  let layers = ((subgraph: highlight-0, edge-style: style-0), (subgraph: highlight-1, edge-style: style-1), (subgraph: highlight-2, edge-style: style-2), )
  let covered = subgraph.select(raw, hedges: layers.map(layer => subgraph.hedges(layer.subgraph)).flatten())
  let outside = subgraph.complement(raw, covered)
  layers = ((subgraph: outside, edge-style: (stroke: outside-stroke)),) + layers
  let style = config.at("style", default: (:))
  let node-label-default = if options.at("show-node-index", default: false) or options.at("debug", default: false) { (:) } else { (node-label: none) }
  let style = node-label-default + (node-style: node => { let matches = (node-0, node-1, node-2, ).filter(layer => node.vid in layer.nodes); physics.node-style + if matches.len() > 0 { matches.last().style } else { (stroke: outside-stroke) } }) + if style == none { (:) } else { style }
  let draw = config.at("draw", default: (:))
  let draw = if draw == none { (:) } else { draw }
  let drawing = (
    draw: (subgraph: layers, subgraph-edge-underlay: "casings"),
  )
  physics-layout.render-layout(
    config + (options: options, style: style, draw: drawing.at("draw", default: (:)) + draw),
    input: raw,
    graph: graph,
    renderer: renderer,
    physics: physics,
    edge-style: (map: particle-map, default-edge: physics.default-edge),
  )
}
