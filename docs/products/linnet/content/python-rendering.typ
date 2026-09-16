#import "../../shared.typ": source-link, product-link

#let python-rendering = [
= Rendering graphs from Python

A Linnet graph can carry ordinary Python objects and render directly in a notebook. This guide
uses a document-processing workflow to explore the default drawing, layout selection, and
styles computed from application data. Begin with the #link("quickstart/python/")[Python
quickstart] for installation and basic graph queries.

The notebook below loads automatically when it comes into view. Python, Linnet, and the Typst
renderer run in your browser through WebAssembly. The first visit downloads their dependencies;
no local Python installation or notebook server is needed.

== Try the rendering API

Choose a layout and enable the custom line theme or half-edge IDs. Each change rebuilds the
example graph and redraws it. The cells are editable, so you can change the workflow or its
selectors and run the affected cell with its play button.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "rendering_api",
    "aria-label": "Linnet Python rendering notebook",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This interactive notebook requires JavaScript and this documentation build's notebook
      assets. The #source-link("crates/linnet-py/examples/rendering_api.py", label: "rendering notebook source")
      is also available to run locally.
    ]
  ]
} else {
  [The online guide includes an interactive notebook. Open the
  #source-link("crates/linnet-py/examples/rendering_api.py", label: "rendering notebook source")
  to run the same example locally.]
}

== What to explore

- Switch between *Force* and *Stable layered*. Both render the same graph; force layout tunes
  geometric positions, while the layered layout organizes the workflow in its chosen direction.
- Enable *Custom Python line theme*. Edge selectors inspect each channel's protocol to choose
  its color, label, and pattern. Endpoint selectors place direction marks using source/sink flow.
- Enable *Half-edge IDs* to see how the incoming request and outgoing publication differ from
  internal edges. The notebook reports the node, edge, and half-edge counts below the drawings.
- Switch the node store and inspect the payload checks. Node, edge, and endpoint dataclass
  instances remain attached by identity, including their shared runtime object.

== The rendering boundary

Returning a `Graph` from a cell invokes its SVG representation. `graph.to_svg()` explicitly
produces the same kind of drawing; the notebook shows both forms. The graph's `RenderConfig`
combines layout options, drawing defaults, and Python selectors. Only topology and the selectors'
typed drawing results pass to Typst; application payloads stay in Python.

The optional theme demonstrates that domain conventions belong in those selectors. It does not
change the workflow records or require them to be serializable. For the complete types and
method signatures, use the #link("reference/python/")[Python API reference].

To experiment with a DOT-backed physics example, continue to
#product-link("gammaloop", label: "GammaLoop DOT input", page: "guides/dot-input/").
For a continuously updating view of the force solver itself, use the
#link("playground/")[live DOT playground].
]
