#import "../shared.typ": *
#import "content/overview.typ": overview
#import "content/tutorial.typ": tutorial
#import "content/networks.typ": networks
#import "content/api.typ": api
#import "content/changelog.typ": changelog

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
  ],
)

#manual
