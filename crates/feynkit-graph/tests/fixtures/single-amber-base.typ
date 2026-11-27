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
  "phi": (source:source-stroke(c: palette.ink, thickness: massive, dash: dashed), sink:sink-stroke(c: palette.ink, thickness: massive, dash: dashed), label:$ "phi" $),
)

#context {
  let raw = build({
    node(<v0>, id: 0, inspection: (node: 0, edges: (0,1,2,)), feynkit-name: "left", interaction: "V_1", numerator: "1")
    node(<v1>, id: 1, inspection: (node: 1, edges: (1,2,3,)), feynkit-name: "right", interaction: "V_1", numerator: "1")
    edge(<e0>, sink(<v0>), id: 0, orientation: "undirected", particle: "phi", pdg: 25, directed: false, numerator: "1", inspection: (edge: 0, source: none, sink: 0), momentum: "p_0", momentum-signature: (loops: (0,), external: (1,0,)), external-state: "incoming", external-index: 0, external-name: "p1", external-connection: 0)
    edge(source(<v0>), <e1>, sink(<v1>), id: 1, orientation: "undirected", particle: "phi", pdg: 25, directed: false, numerator: "1", inspection: (edge: 1, source: 0, sink: 1), momentum: "k_0", momentum-signature: (loops: (1,), external: (0,0,)))
    edge(source(<v0>), <e2>, sink(<v1>), id: 2, orientation: "undirected", particle: "phi", pdg: 25, directed: false, numerator: "1", inspection: (edge: 2, source: 0, sink: 1), momentum: "-k_0+p_0", momentum-signature: (loops: (-1,), external: (1,0,)))
    edge(source(<v1>), <e3>, id: 3, orientation: "undirected", particle: "phi", pdg: 25, directed: false, numerator: "1", inspection: (edge: 3, source: 1, sink: none), momentum: "p_0", momentum-signature: (loops: (0,), external: (1,0,)), external-state: "outgoing", external-index: 1, external-name: "p2", external-connection: 1)
  }, name: "bubble", data: (symmetry-factor: 1, overall-factor: "1", numerator: "1", loop-count: 1))
  raw = graph.map(raw, edge: e => (inspection: e.data.inspection + (
    source-hedge: if e.source == none { none } else { e.source.hedge },
    sink-hedge: if e.sink == none { none } else { e.sink.hedge },
  )))
  let highlighted = subgraph.select(raw, source: (1, 3, ), sink: (0, 2, ))
  let selected-nodes = (0, 1, )
  let selected-stroke = rgb("#ffd166") + 1.2pt
  let selected-edge-style = edge => (stroke: physics.source-style(edge, map: particle-map).stroke + (paint: rgb("#ffd166")))
  let outside-stroke = (paint: rgb("#777777").transparentize(55%), thickness: 0.6pt, dash: "dotted")

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
      highlighted = subgraph.select(raw, hedges: graph.edges(raw).map(e => (e.source, e.sink).filter(h => h != none and subgraph.contains(highlighted, h.origin)).map(h => h.hedge)).flatten())
    }
  }
  let style = config.at("style", default: (:))
  let node-label-default = if options.at("show-node-index", default: false) or options.at("debug", default: false) { (:) } else { (node-label: none) }
  let style = node-label-default + (node-style: node => physics.node-style + if node.vid in selected-nodes { (stroke: selected-stroke) } else { (stroke: outside-stroke) }) + if style == none { (:) } else { style }
  let draw = config.at("draw", default: (:))
  let draw = if draw == none { (:) } else { draw }
  let drawing = (
    draw: (subgraph: ((subgraph: subgraph.complement(raw, highlighted), edge-style: (stroke: outside-stroke)), (subgraph: highlighted, edge-style: selected-edge-style)), subgraph-edge-underlay: false),
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
