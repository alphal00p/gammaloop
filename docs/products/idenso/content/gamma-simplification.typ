#import "../../shared.typ": source-link, product-link

#let gamma-simplification = [
= Gamma simplification notebook

Explore the current Dirac simplifier through executable identities, configurable chain
ordering, and examples of expressions that deliberately remain explicit. Each identity
card checks its expected result and verifies that a second pass leaves it unchanged.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "gamma_simplification",
    "aria-label": "Idenso gamma simplification notebook",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This notebook requires JavaScript and the documentation build's notebook assets.
      Open the #source-link("examples/notebooks/gamma_simplification.py", label: "Marimo source")
      to run it locally.
    ]
  ]
} else {
  [Open the #source-link("examples/notebooks/gamma_simplification.py", label: "Marimo source")
  to run the interactive notebook.]
}

== Capabilities and conventions

The examples cover contracted open chains, four-dimensional Chisholm reduction, symbolic
Lorentz dimension with independent spinor trace normalization, slashed momenta, ordinary
traces, canonical ordering, gamma-five, chiral projectors, the optional three-gamma epsilon
expansion, and charge conjugation. A small adjustable trace illustrates the growth in the
number of metric pairings.

The default repeated-pair ordering avoids unnecessary open-chain expansion. Canonical
ordering can expose cancellations between different gamma orders. Gamma-five rules remain
strictly four-dimensional: the simplifier does not choose a dimensional-regularization
prescription. Mixed dimensions and unsupported transposed words can remain unevaluated.
The notebook checks these boundaries as well as successful reductions.

See #product-link("idenso", page: "reference/form-color-dirac/", label: "Shipped color and Dirac rules")
for the conventions and #product-link("idenso", page: "guides/showcase/", label: "the tensor display showcase")
for additional notation controls.

== Run locally

Use an interpreter containing the combined Symbolica community host with Spenso and Idenso,
Marimo 0.24.0, and Typst 0.15.0. Installing ordinary Symbolica alone does not supply the
community extension.

// docs-example: syntax
```sh
just notebook gamma_simplification /path/to/python
```
]
