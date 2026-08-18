#import "../shared.typ": *
#import "content/overview.typ": overview
#import "content/tutorial.typ": tutorial
#import "content/networks.typ": networks
#import "content/api.typ": api
#import "content/changelog.typ": changelog
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
    #tutorial
    #networks
    #api
    #changelog
    #spenso-macros-releases
    #spenso-hep-lib-releases
    #spynso3-releases
  ],
)

#manual
