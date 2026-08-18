#import "../shared.typ": *
#import "content/overview.typ": overview
#import "content/tutorial.typ": tutorial
#import "content/algorithms.typ": algorithms
#import "content/clinnet.typ": clinnet
#import "content/linnest.typ": linnest
#import "content/api.typ": api
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
    #changelog
    #linnet-releases
    #clinnet-releases
  ],
)

#manual
