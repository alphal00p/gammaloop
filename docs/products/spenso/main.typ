#import "../shared.typ": *
#import "content/overview.typ": overview
#import "content/quickstart.typ": quickstart
#import "content/quickstart-python.typ": quickstart-python
#import "content/quickstart-rust.typ": quickstart-rust
#import "content/tutorial.typ": tutorial
#import "content/networks.typ": networks
#import "content/python.typ": python
#import "content/api.typ": api
#import "content/changelog.typ": changelog
#import "content/spenso-releases.typ": spenso-releases
#import "content/spenso-macros-releases.typ": spenso-macros-releases
#import "content/spenso-hep-lib-releases.typ": spenso-hep-lib-releases
#import "content/spynso3-releases.typ": spynso3-releases

#let manual = product-document(
  title: "Spenso",
  tagline: "Typed tensors, symbolic structures, and executable tensor networks",
  version: "crates/spenso/Cargo.toml",
  owner: "Spenso project",
  body: [
    #overview
    #quickstart
    #quickstart-python
    #quickstart-rust
    #tutorial
    #networks
    #python
    #api
    #changelog
    #spenso-releases
    #spenso-macros-releases
    #spenso-hep-lib-releases
    #spynso3-releases
  ],
)

#manual
