#import "../../shared.typ": source-link

#let playground = [
= Live Python playground

Edit DOT below and watch the layout update as the solver runs. The notebook loads automatically
when it comes into view.
Open the collapsible layout controls to change the solver settings. Editing DOT or a slider
restarts the layout; pause and resume keep the current positions.

The first visit downloads the browser Python runtime and notebook dependencies. Python and
Linnet then run locally in your browser through WebAssembly; no Python installation or remote
notebook server is needed. An internet connection is needed for the initial downloads.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "layout_stream",
    "aria-label": "Live DOT layout notebook",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This interactive notebook needs JavaScript and the notebook assets for this documentation
      build. The #source-link("crates/linnet-py/examples/layout_stream.py", label: "notebook source")
      is also available to run locally.
    ]
  ]
} else {
  [Open the #source-link("crates/linnet-py/examples/layout_stream.py", label: "streaming layout notebook source")
  to run this example locally. The online documentation embeds the interactive notebook.]
}

== Geometry and final drawings

The live SVG previews the solver's positions and edge geometry. Use
#link("guides/linnest/")[Linnest] or the Python rendering API for final Typst drawings with
measured labels, styles, and decorations.

== Try the Python API

The #link("quickstart/python/")[Python quickstart] has an editable notebook containing its graph
example. Edit the Python cell and use its play button to inspect the resulting graph.

The #link("guides/python-rendering/")[Python rendering guide] provides a complete editable
workflow with custom styles and arbitrary Python payloads.

The #link("reference/python/")[Python reference] documents the full API. Its examples can depend
on earlier definitions, local files, or rendering dependencies, so they are not all runnable
as independent notebook cells yet.
]
