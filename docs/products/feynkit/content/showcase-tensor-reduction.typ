#import "../../shared.typ": source-link

#let showcase-tensor-reduction = [
= Vacuum tensor reduction

The `feynkit-tensor` showcase distinguishes integrated and external momenta, preserves the symbolic dimension, and uses contraction symmetry to keep high-rank results compact. It then reduces the numerator of a two-loop gluon vacuum graph.

The notebook runs in the browser when this site's notebook assets are available. Its cells
are editable; run a changed cell to update the dependent results. The first load downloads
Symbolica, the notebook runtime, and this example's data.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "07_tensor_reduction_marimo",
    "aria-label": "FeynKit: Vacuum tensor reduction",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This notebook needs JavaScript and this documentation build's notebook assets.
      Open the #source-link("examples/notebooks/feynkit/07_tensor_reduction_marimo.py", label: "notebook source")
      to run the same example locally.
    ]
  ]
} else {
  [The online guide includes an interactive notebook. Open the
  #source-link("examples/notebooks/feynkit/07_tensor_reduction_marimo.py", label: "notebook source")
  to run it locally.]
}

== What to explore

Compare symbolic projectors of different ranks, then inspect the scalar graph terms. The projection is valid under Lorentz-invariant vacuum integration; it is not an identity of the original unintegrated numerator.

== Run locally

From a checkout, use a Python environment containing the shared Symbolica host and the
notebook dependencies. The #link("guides/showcases/")[showcase gallery] gives the setup.

// docs-example: syntax
```sh
just notebook feynkit/07_tensor_reduction_marimo /path/to/python
```

Continue with the #link("guides/showcases/")[other FeynKit showcases] or the
#link("reference/python/feynkit-community/")[Python API reference].
]
