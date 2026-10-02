// Default V1 per-graph Typst renderer bundled with Linnest consumers.
// Renderers import this module and call render(config); configuration is native
// Typst data rather than command-line string inputs.
#import "layout.typ": layout-spec

#let render(config) = {
  set page(width: auto, height: auto, margin: (x: 2mm, y: 2mm))
  context layout-spec(config)
}
