// Clinnet's default figure delegates drawing to Linnest's generic renderer.
#import "layout.typ": layout

#let render(config) = {
  set page(width: auto, height: auto, margin: (x: 2mm, y: 2mm))
  context layout(config)
}
