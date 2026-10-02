#import "../../shared.typ": source-link

#let showcase-models-and-diagrams = [
= Models, generation, and graphs

This notebook covers `feynkit-model`, `feynkit-generator`, and `feynkit-graph` together. It shows immutable parameter updates, inclusive loop ranges, model-owned generation, JSON/DOT interchange, and signed momentum routing.

The notebook runs in the browser when this site's notebook assets are available. Its cells
are editable; run a changed cell to update the dependent results. The first load downloads
Symbolica, the notebook runtime, and this example's data.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "01_models_and_diagrams_marimo",
    "aria-label": "FeynKit: Models, generation, and graphs",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This notebook needs JavaScript and this documentation build's notebook assets.
      Open the #source-link("examples/notebooks/feynkit/01_models_and_diagrams_marimo.py", label: "notebook source")
      to run the same example locally.
    ]
  ]
} else {
  [The online guide includes an interactive notebook. Open the
  #source-link("examples/notebooks/feynkit/01_models_and_diagrams_marimo.py", label: "notebook source")
  to run it locally.]
}

== What to explore

Change the parameter card and inspect which dependent values require recomputation. Compare tree and one-loop graphs, then follow a momentum signature through a chosen basis.

== Run locally

From a checkout, use a Python environment containing the shared Symbolica host and the
notebook dependencies. The #link("guides/showcases/")[showcase gallery] gives the setup.

// docs-example: syntax
```sh
just notebook feynkit/01_models_and_diagrams_marimo /path/to/python
```

Continue with the #link("guides/showcases/")[other FeynKit showcases] or the
#link("reference/python/feynkit-community/")[Python API reference].
]
