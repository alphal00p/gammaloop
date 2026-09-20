#import "../../shared.typ": source-link

#let showcase-ufo = [
= Loading UFO models

The `feynkit-ufo` showcase uses the repository's scalar UFO fixture and its default restriction card. After normalization, generation uses the same typed `Model` API as the other showcases.

The notebook runs in the browser when this site's notebook assets are available. Its cells
are editable; run a changed cell to update the dependent results. The first load downloads
Symbolica, the notebook runtime, and this example's data.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "04_ufo_loading_marimo",
    "aria-label": "FeynKit: Loading UFO models",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This notebook needs JavaScript and this documentation build's notebook assets.
      Open the #source-link("examples/notebooks/feynkit/04_ufo_loading_marimo.py", label: "notebook source")
      to run the same example locally.
    ]
  ]
} else {
  [The online guide includes an interactive notebook. Open the
  #source-link("examples/notebooks/feynkit/04_ufo_loading_marimo.py", label: "notebook source")
  to run it locally.]
}

The notebook pins the Symbolica 3-compatible UFO loader revision
#link("https://github.com/alphal00p/ufo_model_loader/tree/70ddee6b416f8c8b340e0d087646d77095c5d24b")[`70ddee6b`]
(version 0.1.8). Its square-root parsing and index-wrapping rules match the current Symbolica
API. The browser uses a wheel built from the same source revision.

== What to explore

Inspect the normalization counts, numerical parameter card, and tree-diagram generation report. Native execution needs Python 3.11 or newer and `ufo-model-loader`; the browser notebook installs the matching dependency.

== Run locally

From a checkout, use a Python environment containing the shared Symbolica host and the
notebook dependencies. The #link("guides/showcases/")[showcase gallery] gives the setup.

// docs-example: syntax
```sh
just notebook feynkit/04_ufo_loading_marimo /path/to/python
```

Continue with the #link("guides/showcases/")[other FeynKit showcases] or the
#link("reference/python/feynkit-community/")[Python API reference].
]
