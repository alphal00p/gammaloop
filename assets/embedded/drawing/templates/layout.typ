#import "linnest.typ": draw, graph, layout as apply-layout
#import "edge-style.typ" as edge-style

#let layout(
  input,
  split-edge: true,
  scope: (:),
  columns: 1fr,
  unit: 1,
  typst-fields: "plain",
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
      source-style: edge => edge-style.source-style(edge, typst-fields: typst-fields),
      sink-style: edge => edge-style.sink-style(edge, typst-fields: typst-fields),
      edge-label: edge => edge-style.edge-label(edge, typst-fields: typst-fields),
    ))
  }
  for d in diags {
    d
  }
}
