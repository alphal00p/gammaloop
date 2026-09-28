#import "../../shared.typ": source-link

#let showcase-first-diagram = [
= A first diagram

This tour joins the model, generator, and graph APIs through `symbolica.community.feynkit`. Inspect the generation report, the diagram-wide factor, and the particle metadata on every edge.

The notebook runs in the browser when this site's notebook assets are available. Its cells
are editable; run a changed cell to update the dependent results. The first load downloads
Symbolica, the notebook runtime, and this example's data.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "00_quickstart_marimo",
    "aria-label": "FeynKit: A first diagram",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This notebook needs JavaScript and this documentation build's notebook assets.
      Open the #source-link("examples/notebooks/feynkit/00_quickstart_marimo.py", label: "notebook source")
      to run the same example locally.
    ]
  ]
} else {
  [The online guide includes an interactive notebook. Open the
  #source-link("examples/notebooks/feynkit/00_quickstart_marimo.py", label: "notebook source")
  to run it locally.]
}

== What to explore

Change the allowed interaction or loop order and compare the resulting diagrams. The model and graph stay typed throughout, while symbolic factors remain native Symbolica expressions.

== Run locally

From a checkout, use a Python environment containing the shared Symbolica host and the
notebook dependencies. The #link("guides/showcases/")[showcase gallery] gives the setup.

// docs-example: syntax
```sh
just notebook feynkit/00_quickstart_marimo /path/to/python
```

Continue with the #link("guides/showcases/")[other FeynKit showcases] or the
#link("reference/python/feynkit-community/")[Python API reference].
]
