#import "../shared.typ": *
#import "content/overview.typ": overview
#import "content/tutorial.typ": tutorial
#import "content/algebra.typ": algebra
#import "content/form-symbolica-color-and-dirac.typ": form-color-dirac-content
#import "content/api.typ": api
#import "content/changelog.typ": changelog

#let manual = product-document(
  title: "Idenso",
  tagline: "Symbolic tensor identities for Spenso and Symbolica",
  version: "crates/idenso/Cargo.toml",
  owner: "Idenso project",
  body: [
    #overview
    #tutorial
    #algebra
    #form-color-dirac-content("manual")
    #api
    #changelog
  ],
)

#manual
