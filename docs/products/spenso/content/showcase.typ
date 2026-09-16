#import "../../shared.typ": source-link, product-link

#let showcase = [
= Spenso + Idenso live showcase

Build a typed tensor expression, simplify its Dirac trace and contractions, and inspect
how its mathematical display changes. This Marimo notebook shares one Symbolica kernel
between Spenso and Idenso. Its concrete tensor network is drawn with Linnet and evaluated
to a checked scalar result.

== Explore the notebook

Choose Ports, Schoonschip, or Function call notation, then change dimensions, spacing,
parentheses, or symbol scripts. The cells are editable: changing a constructor or an
Idenso pass recomputes the affected results. The first visit downloads the browser Python
runtime and this documentation version's WebAssembly wheels.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "spenso_idenso_display",
    "aria-label": "Spenso and Idenso tensor display showcase",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This notebook requires JavaScript and the documentation build's notebook assets.
      Open the #source-link("examples/notebooks/spenso_idenso_display.py", label: "Marimo source")
      to run it locally.
    ]
  ]
} else {
  [Open the #source-link("examples/notebooks/spenso_idenso_display.py", label: "Marimo source")
  to run the interactive showcase.]
}

== What the examples establish

- Typed vector products and Dirac traces retain their representation and free-index structure.
- Idenso simplifies the underlying Symbolica expressions; Spenso presents the resulting tensor
  semantics using the selected display settings.
- A concrete two-vector network separates the tensor view, graph, DOT source, and scalar
  result. Contracting `(1, 2)` with `(3, 4)` gives `11`.

For the full workflows, see
#product-link("spenso", page: "guides/python/", label: "Python tensor workflows") and
#product-link("idenso", page: "guides/algebra/", label: "algebra and rewrites").
#product-link("feynkit", page: "guides/showcases/", label: "FeynKit's showcase gallery")
connects these symbolic tools to models, Feynman diagrams, CFF, and tensor reduction.

== Run and edit locally

From a checkout, use an interpreter with the combined Symbolica community host containing
Spenso and Idenso, plus Marimo 0.24.0, Typst 0.15.0, and this checkout's Linnet Python wheel:

// docs-example: syntax
```sh
just notebook spenso_idenso_display /path/to/python
```

The notebook source is shared by the browser export and the local editor. The browser
export installs the versioned wheels; a local editor uses the chosen interpreter's packages.
]
