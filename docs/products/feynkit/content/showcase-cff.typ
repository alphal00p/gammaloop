#import "../../shared.typ": source-link

#let showcase-cff = [
= Cross-Free Families

The `feynkit-cff` showcase starts from a one-loop scalar diagram. Edge IDs connect orientation constraints to the graph; each generated surface retains the energy and external-momentum data behind its Symbolica symbol.

The notebook runs in the browser when this site's notebook assets are available. Its cells
are editable; run a changed cell to update the dependent results. The first load downloads
Symbolica, the notebook runtime, and this example's data.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "02_cff_and_symbolica_marimo",
    "aria-label": "FeynKit: Cross-Free Families",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This notebook needs JavaScript and this documentation build's notebook assets.
      Open the #source-link("examples/notebooks/feynkit/02_cff_and_symbolica_marimo.py", label: "notebook source")
      to run the same example locally.
    ]
  ]
} else {
  [The online guide includes an interactive notebook. Open the
  #source-link("examples/notebooks/feynkit/02_cff_and_symbolica_marimo.py", label: "notebook source")
  to run it locally.]
}

== What to explore

Inspect the first orientation and its denominator products, then change the symbolic coupling multiplying the result. The surface arena remains available alongside the expression.

== Run locally

From a checkout, use a Python environment containing the shared Symbolica host and the
notebook dependencies. The #link("guides/showcases/")[showcase gallery] gives the setup.

// docs-example: syntax
```sh
just notebook feynkit/02_cff_and_symbolica_marimo /path/to/python
```

Continue with the #link("guides/showcases/")[other FeynKit showcases] or the
#link("reference/python/feynkit-community/")[Python API reference].
]
