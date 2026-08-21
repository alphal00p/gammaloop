#import "crates/linnest/typst/src/lib.typ": (
  draw, graph, layout as apply-layout, physics, subgraph,
)
#import "edge-style.typ" as edge-style
#import "layout-core.typ": bind-layout, bind-render

// Save-dot adapter: bind the extracted Linnest package and the model-specific
// edge style to GammaLoop's shared layout algorithm.
#let layout = bind-layout(
  draw: draw,
  graph: graph,
  subgraph: subgraph,
  apply-layout: apply-layout,
  physics: physics,
  edge-style: edge-style,
)

#let render-layout = bind-render(
  draw: draw,
  graph: graph,
  subgraph: subgraph,
  apply-layout: apply-layout,
  physics: physics,
  edge-style: edge-style,
)
