// Decorate a snapshot of the owner's render configuration. Structural selections
// are graph-bound in Python; these half-edge/node IDs refer to that checked revision.
#let focus(config, hedges, nodes) = {
  let selected-stroke = rgb("#c58b13") + 1.2pt
  let outside-stroke = (paint: rgb("#777777").transparentize(55%), thickness: 0.6pt, dash: "dotted")
  let elements = config.elements
  elements.nodes = elements.nodes.enumerate().map(((i, record)) => {
    let original = record.at("node-style", default: (:))
    record + (node-style: node => {
      let style = if type(original) == function { original(node) } else { original }
      let style = if style == none { (:) } else { style }
      style + if i in nodes {
        (stroke: selected-stroke, fill: rgb("#ffd166").transparentize(65%))
      } else {
        (stroke: outside-stroke, fill: rgb("#777777").transparentize(85%))
      }
    })
  })
  let drawing = config.at("draw", default: (:))
  if drawing == none { drawing = (:) }
  config + (
    elements: elements,
    draw: drawing + (
      subgraph: (
        (subgraph: hedges.map(bit => not bit), edge-style: (stroke: outside-stroke)),
        (subgraph: hedges, edge-style: (stroke: selected-stroke)),
      ),
      subgraph-edge-underlay: false,
    ),
  )
}
