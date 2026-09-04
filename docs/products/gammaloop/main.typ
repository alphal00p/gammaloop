#import "../shared.typ": *
#import "content/overview.typ": overview
#import "content/quickstart.typ": quickstart
#import "content/quickstart-cli.typ": quickstart-cli
#import "content/quickstart-python.typ": quickstart-python
#import "content/quickstart-rust.typ": quickstart-rust
#import "content/physics-scope.typ": physics-scope
#import "content/conventions.typ": conventions
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
    #quickstart
    #quickstart-cli
    #quickstart-python
    #quickstart-rust
    #physics-scope
    #conventions
    #tutorial
    #process-generation
    #events
    #api
    #diagnostics
    #kurvst
    #changelog
  ],
)

#manual
