// Clinnet's public template path adapts its DOT input to Linnest's generic renderer.
#import "crates/linnest/typst/src/render/layout.typ": *
#import "crates/linnest/typst/src/lib.typ": graph

#let layout(config) = {
  let path = config.at("data-path", default: none)
  if path == none {
    panic("render config requires data-path")
  }
  for parsed in graph.parse(read(path)) {
    layout-graph(config, parsed)
  }
}
