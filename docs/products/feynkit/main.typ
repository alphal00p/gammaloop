#import "../shared.typ": *
#import "content/overview.typ": overview
#import "content/quickstart.typ": quickstart
#import "content/quickstart-rust.typ": quickstart-rust
#import "content/quickstart-python.typ": quickstart-python
#import "content/tutorial.typ": tutorial
#import "content/tensor-reduction.typ": tensor-reduction
#import "content/notebooks.typ": notebooks
#import "content/showcases.typ": showcases
#import "content/showcase-first-diagram.typ": showcase-first-diagram
#import "content/showcase-models-and-diagrams.typ": showcase-models-and-diagrams
#import "content/showcase-cff.typ": showcase-cff
#import "content/showcase-kinematics.typ": showcase-kinematics
#import "content/showcase-ufo.typ": showcase-ufo
#import "content/showcase-tensor-reduction.typ": showcase-tensor-reduction
#import "content/community-host.typ": community-host
#import "content/api.typ": api
#import "content/changelog.typ": changelog

#let manual = product-document(
  title: "FeynKit",
  tagline: "Standalone model, diagram, and symbolic-physics tools",
  version: "crates/feynkit/Cargo.toml",
  owner: "FeynKit maintainers",
  body: [
    #overview
    #quickstart
    #quickstart-rust
    #quickstart-python
    #tutorial
    #tensor-reduction
    #notebooks
    #showcases
    #showcase-first-diagram
    #showcase-models-and-diagrams
    #showcase-cff
    #showcase-kinematics
    #showcase-ufo
    #showcase-tensor-reduction
    #community-host
    #api
    #changelog
  ],
)

#manual
