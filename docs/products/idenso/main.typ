#import "../shared.typ": *
#import "content/overview.typ": overview
#import "content/quickstart.typ": quickstart
#import "content/quickstart-python.typ": quickstart-python
#import "content/quickstart-rust.typ": quickstart-rust
#import "content/tutorial.typ": tutorial
#import "content/algebra.typ": algebra
#import "content/form-symbolica-color-and-dirac.typ": form-color-dirac-content
#import "content/api.typ": api
#import "content/changelog.typ": changelog
#import "content/idenso-releases.typ": idenso-releases

#let manual = product-document(
  title: "Idenso",
  tagline: "Symbolic tensor identities for Spenso and Symbolica",
  version: "crates/idenso/Cargo.toml",
  owner: "Idenso project",
  body: [
    #overview
    #quickstart
    #quickstart-python
    #quickstart-rust
    #tutorial
    #algebra
    #form-color-dirac-content("manual")
    #api
    #changelog
    #idenso-releases
  ],
)

#manual
