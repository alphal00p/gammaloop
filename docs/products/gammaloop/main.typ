#import "../shared.typ": *
#import "content/overview.typ": overview
#import "content/tutorial.typ": tutorial
#import "content/process-generation.typ": process-generation
#import "content/events.typ": events
#import "content/diagnostics.typ": diagnostics
#import "content/kurvst.typ": kurvst
#import "content/api.typ": api
#import "content/changelog.typ": changelog

#let manual = product-document(
  title: "GammaLoop",
  tagline: "Differential cross-sections with Local Unitarity",
  version: "pyproject.toml and crates/gammalooprs/Cargo.toml",
  owner: "GammaLoop project",
  body: [
    #overview
    #tutorial
    #process-generation
    #events
    #diagnostics
    #kurvst
    #api
    #changelog
  ],
)

#manual
