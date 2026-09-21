#import "crates/linnest/typst/src/lib.typ": graph
#import "crates/linnest/typst/src/render/layout.typ" as renderer
#import "physics-edge-style.typ" as physics
#import "edge-style.typ" as edge-style
#import "layout-core.typ": bind-layout, bind-render

// Save-dot adapter: bind the extracted Linnest package and the model-specific
// particle styles. Linnest owns the shared measurement and layout pipeline.
#let layout = bind-layout(
  graph: graph,
  renderer: renderer,
  physics: physics,
  edge-style: edge-style,
  diagram-options: (title: none),
)
#let render-layout = bind-render(
  graph: graph,
  renderer: renderer,
  physics: physics,
  edge-style: edge-style,
  diagram-options: (title: none),
)
