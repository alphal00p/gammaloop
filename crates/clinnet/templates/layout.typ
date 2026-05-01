#import "linnest.typ": draw, graph, layout as apply-layout

#let layout(
  input,
  split_edge: true,
  scope: (:),
  columns: 1fr,
  unit: 1,
  additional_data: (:),
) = {
  let graphs = graph.parse(input)
  let diags = ()
  for graph-bytes in graphs {
    let graph-bytes = apply-layout(graph-bytes, ..additional_data)
    diags.push(draw(graph-bytes, scope: scope, unit: unit, title: auto))
  }
  for d in diags{
    d
  }
}
