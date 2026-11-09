// The same bindings drive Tidy previews and standalone documentation cells.
#import "../src/lib.typ": draw, graph, layout, layouts, subgraph
#let options = layouts.options

#let api-link(prefix, name, display: none) = {
  raw(if display == none { name } else { display }, lang: none)
}

#let scope = (
  api-link: api-link,
  draw: draw,
  graph: graph,
  layout: layout,
  options: options,
  subgraph: subgraph,
)
