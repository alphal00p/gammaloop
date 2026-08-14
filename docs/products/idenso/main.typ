#import "../shared.typ": *
#import "content/overview.typ": overview
#import "content/algebra.typ": algebra
#import "content/api.typ": api
#import "content/changelog.typ": changelog

#let manual = product-document(
  title: "Idenso",
  tagline: "Symbolic tensor identities for Spenso and Symbolica",
  version: "crates/idenso/Cargo.toml",
  owner: "Idenso project",
  body: [
    #overview
    #algebra
    #api
    #changelog
  ],
)

#manual
