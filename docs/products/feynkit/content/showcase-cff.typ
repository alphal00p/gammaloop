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

Displaying `diagram.cross_free_family()` returns the interactive display of
`symbolica.community.hepkit.CrossFreeFamily`. The compact view starts with the
factored expression and selected surface definitions. Expand "Explore graph and
families" to edit an orientation by clicking an arrowhead, choose a family
preview, or compare surface regions on the native Linnet graph. Shift-click
surface factors to select several surfaces. Shared
denominator factors stay outside the family sum; the active family's path is
highlighted in that expression. Selection only changes the display, never the
underlying CFF result.

The display's $C$ is shorthand for the denominator sum. Numerators, on-shell
energy prefactors and the spatial integration measure stay separate. Use
`to_expression()` to continue symbolically, or `to_expression(normalized=True)`
to include the generated energy product and loop measure. The surface arena
remains available alongside the expression.

== Run locally

From a checkout, use a Python environment containing the shared Symbolica host and the
notebook dependencies. The #link("guides/showcases/")[showcase gallery] gives the setup.

// docs-example: syntax
```sh
just notebook feynkit/02_cff_and_symbolica_marimo /path/to/python
```

Continue with #link("guides/showcases/numerical-integration/")[numerical integration with CFF and LTD], the #link("guides/showcases/")[other FeynKit showcases], or the
#link("reference/python/feynkit-community/")[Python API reference].
]
