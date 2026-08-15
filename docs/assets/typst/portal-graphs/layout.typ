#import "../../../../crates/linnest/typst/src/lib.typ": (
  draw, graph, layout as apply-layout, physics,
)
#import "../../../../assets/embedded/drawing/templates/layout-core.typ": (
  bind-layout,
)
#import "edge-style.typ" as edge-style
#import "../theme.typ": palette

// Website adapter: the graph algorithm is exactly GammaLoop's save-dot core;
// only the transparent, theme-aware presentation differs.
#let layout = bind-layout(
  draw: draw,
  graph: graph,
  apply-layout: apply-layout,
  physics: physics,
  edge-style: edge-style,
  diagram-options: (
    title: none,
    node-fill: none,
    node-stroke: palette.ink + 1.1pt,
  ),
)
