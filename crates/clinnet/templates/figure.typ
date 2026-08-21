// Default V1 per-graph Typst template bundled with the Linnet CLI.
// Clinnet imports this module and calls render(config); configuration is native
// Typst data rather than command-line string inputs.
#import "layout.typ": layout

#let render(config) = {
  set page(width: auto, height: auto, margin: (x: 2mm, y: 2mm))
  context layout(config)
}
