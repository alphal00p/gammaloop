#import "../shared.typ": *
#import "content/overview.typ": overview
#import "content/tutorial.typ": tutorial
#import "content/algorithms.typ": algorithms
#import "content/clinnet.typ": clinnet
#import "content/linnest.typ": linnest
#import "content/api.typ": api
#import "content/typst-reference.typ": typst-reference
#import "content/typst-graph.typ": typst-graph
#import "content/typst-layout.typ": typst-layout
#import "content/typst-drawing.typ": typst-drawing
#import "content/typst-physics.typ": typst-physics
#import "content/typst-subgraph.typ": typst-subgraph
#import "content/changelog.typ": changelog
#import "content/linnet-releases.typ": linnet-releases
#import "content/clinnet-releases.typ": clinnet-releases

#let manual = product-document(
  title: "Linnet",
  tagline: "Half-edge graphs and subgraph-first algorithms",
  version: "crates/linnet/Cargo.toml",
  owner: "Linnet project",
  body: [
    #overview
    #tutorial
    #algorithms
    #clinnet
    #linnest
    #api
    #typst-reference
    #typst-graph
    #typst-layout
    #typst-drawing
    #typst-physics
    #typst-subgraph
    #changelog
    #linnet-releases
    #clinnet-releases
  ],
)

#manual
