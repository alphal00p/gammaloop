#import "../../shared.typ": source-link

#let showcase-four-loop-numerator = [
= Four-loop gluon propagator numerators

Generate all four-loop one-particle-irreducible gluon propagator graphs in QCD with one massless quark flavor, including gluon, ghost and quark loops. Select a graph and calculate its numerator.

The notebook runs in the browser when this site's notebook assets are available. Its cells
are editable; run a changed cell to update the dependent results. The first load downloads
Symbolica, the notebook runtime, and this example's data.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "05_four_loop_numerator_marimo",
    "aria-label": "FeynKit: Four-loop gluon propagator numerators",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This notebook needs JavaScript and this documentation build's notebook assets.
      Open the #source-link("examples/notebooks/feynkit/05_four_loop_numerator_marimo.py", label: "notebook source")
      to run the same example locally.
    ]
  ]
} else {
  [The online guide includes an interactive notebook. Open the
  #source-link("examples/notebooks/feynkit/05_four_loop_numerator_marimo.py", label: "notebook source")
  to run it locally.]
}

== Canonical example

Adapted from the #link("https://symbolica.io/gallery/notebooks/four_loop_numerator.html")[Symbolica Community notebook], pinned to #link("https://github.com/symbolica-dev/symbolica-community/blob/39b556c7d7dcc919f61d125f60ff83f8d49f7ac3/examples/hep/four_loop_numerator.py")[revision `39b556c7`]. Copyright 2024 Ruijl Research, MIT license; the source directory includes the license. The graph selector and scientific calculation cells are preserved.

== What to explore

Search the complete graph selector by name or internal particle content. Inspect tensor reduction, graph weights, the external-state projector, color and Dirac simplification, and the expanded term count. Generation and selected-graph algebra are separate computations; changing the selector reruns the dependent numerator cells. Quark traces and tensor reduction can take substantially longer than loading Python.

== Run locally

Use the shared Symbolica host and dependencies described in the #link("guides/showcases/")[showcase gallery].

// docs-example: syntax
```sh
just notebook feynkit/05_four_loop_numerator_marimo /path/to/python
```
]
