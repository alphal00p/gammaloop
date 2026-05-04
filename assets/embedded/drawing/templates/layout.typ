#import "linnest.typ": draw, graph, layout as apply-layout
#import "edge-style.typ" as edge-style

#let layout(
  input,
  split-edge: true,
  scope: (:),
  columns: 1fr,
  unit: 1,
  additional-data: (:),
) = {
  let graphs = graph.parse(input)
  let diags = ()
  for graph-bytes in graphs {
    let graph-bytes = apply-layout(graph-bytes, ..additional-data)
    diags.push(draw(
      graph-bytes,
      scope: scope,
      unit: unit,
      title: auto,
      source-style: edge-style.source-style,
      sink-style: edge-style.sink-style,
      edge-label: edge-style.edge-label,
    ))
  }
  for d in diags {
    d
  }
}
