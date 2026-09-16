#import "../../shared.typ": source-link

#let showcase-kinematics = [
= Kinematics and jets

The `feynkit-kinematics` showcase makes the mostly-minus metric and component order explicit. It covers on-shell momenta, boosts and their inverses, rotations, rapidity, angular distances, and generalized-kt jet clustering.

The notebook runs in the browser when this site's notebook assets are available. Its cells
are editable; run a changed cell to update the dependent results. The first load downloads
Symbolica, the notebook runtime, and this example's data.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "03_kinematics_and_jets_marimo",
    "aria-label": "FeynKit: Kinematics and jets",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This notebook needs JavaScript and this documentation build's notebook assets.
      Open the #source-link("examples/notebooks/feynkit/03_kinematics_and_jets_marimo.py", label: "notebook source")
      to run the same example locally.
    ]
  ]
} else {
  [The online guide includes an interactive notebook. Open the
  #source-link("examples/notebooks/feynkit/03_kinematics_and_jets_marimo.py", label: "notebook source")
  to run it locally.]
}

== What to explore

Change the momenta, boost velocity, or jet radius. Compare invariant masses before and after a transformation and follow each jet back to its input constituent indices.

== Run locally

From a checkout, use a Python environment containing the shared Symbolica host and the
notebook dependencies. The #link("guides/showcases/")[showcase gallery] gives the setup.

// docs-example: syntax
```sh
just notebook feynkit/03_kinematics_and_jets_marimo /path/to/python
```

Continue with the #link("guides/showcases/")[other FeynKit showcases] or the
#link("reference/python/feynkit-community/")[Python API reference].
]
