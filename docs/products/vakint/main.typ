#import "../shared.typ": *
#import "content/overview.typ": overview
#import "content/tutorial.typ": tutorial
#import "content/evaluation.typ": evaluation
#import "content/api.typ": api
#import "content/changelog.typ": changelog

#let manual = product-document(
  title: "Vakint",
  tagline: "Single-scale vacuum-integral matching and evaluation",
  version: "crates/vakint/Cargo.toml",
  owner: "Vakint project",
  body: [
    #overview
    #tutorial
    #evaluation
    #api
    #changelog
  ],
)

#manual
