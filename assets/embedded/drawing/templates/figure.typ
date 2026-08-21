// Default V1 per-graph Typst template bundled with GammaLoop drawing state.
// Clinnet imports this module and calls render(config) with native Typst values.
#import "layout.typ": render-layout

#let render(config) = {
  set page(width: auto, height: auto, margin: (x: 2mm, y: 2mm))
  context render-layout(config)
}
