#import "../shared.typ": *
#import "content/overview.typ": overview
#import "content/tutorial.typ": tutorial
#import "content/algorithms.typ": algorithms
#import "content/api.typ": api
#import "content/changelog.typ": changelog

#let manual = product-document(
  title: "Linnet",
  tagline: "Half-edge graphs and subgraph-first algorithms",
  version: "crates/linnet/Cargo.toml",
  owner: "Linnet project",
  body: [
    #overview
    #tutorial
    #algorithms
    #api
    #changelog
  ],
)

#manual
