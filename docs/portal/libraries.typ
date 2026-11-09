#import "shared.typ": *
#import "components.typ": *

#let libraries-page = portal-document(
  "Libraries and APIs",
  "Guides, Rust crates, and Python APIs for the αLoop software libraries.",
  active: "libraries",
)[
  #page-hero[
    = Libraries for each part of the calculation

    FeynKit is a collection of tools for perturbative quantum field theory.
    Use its crates independently, alongside Linnet and Spenso, or together in GammaLoop.
    The #link("../#topics")[topic index] groups the documentation by what you want to calculate.
  ]
  #project-grid[
    #for product in registry.product {
      let root = product-root(product.id, root: "../")
      let components = product.rust_components + product.python_components
      let packages = components.map(component => component.package).sorted().dedup()
      let guide = product.pages.find(page => page.group == "Guides")
      project-card(product.id)[
        == #link(root)[#product.title]

        #project-summary[#product.tagline]

        #package-list(product.title, packages)

        #project-links(product.title)[
          #link(page-link((product: product.id, page: "quickstart"), root: "../"))[Get started]
          #link(root + guide.route)[Guides]
          #link(root + "reference/rust/")[Rust API]
          #link(root + "reference/python/")[Python API]
        ]
      ]
    }
  ]
]
