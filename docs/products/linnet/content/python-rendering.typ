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

== Default layout and controls

`Graph.to_svg()`, `Graph.render()` and Feynman diagram rendering use the shared
EC-planarization and ImPrEd pipeline by default. The same implementation runs
in direct Typst and in CLI drawings. Layout controls remain explicit:

```python
config = linnet.RenderSettings(
    layouts=linnet.LayoutSettings(
        algorithm=linnet.LayoutAlgorithm.Impred,
        impred_parallel_balance=1.0,
        impred_pull=0.45,
        impred_external_max_points=2,
    )
)
svg = graph.to_svg(config=config)
```

The #link("reference/typst/layout/")[layout guide] describes every ImPrEd control.
Python uses underscores where the Typst names use hyphens. Explicit layouts
such as `LayoutAlgorithm.Force` continue to select their respective algorithms.

== The rendering boundary

Returning a `Graph` from a cell invokes its interactive HTML representation. `graph.to_svg()`
explicitly produces the same native drawing with hover labels, pinned click details, and local
selection using Shift-, Ctrl-, or Meta-click. These interactions do not update Python graph
selections automatically. In Marimo, direct rich display enables the embedded script; when
embedding the SVG explicitly, use `mo.iframe(graph.to_svg())`.
Drag the drawing to pan, or use Ctrl/Meta-scroll to zoom gently around the pointer.
At 100%, SVG points retain their physical size, so equally styled labels, vertices,
and strokes have the same size across drawings. The viewport grows vertically to
keep the whole drawing visible. Narrow columns or explicitly bounded previews
can shrink the drawing to fit; the zoom indicator reports that actual scale.
With the graph focused, `+` and `-` zoom in five-percent steps and `0` fits the
drawing without enlarging it above 100%. Hover previews details beside the graph when space permits,
and below it in narrow outputs. The preview disappears when the pointer leaves the element;
keyboard focus also previews details. Clicking pins the details until the panel is closed.
The graph's `RenderSettings` combines layout options, drawing defaults, and Python selectors.
Ordinary `Graph` and `Subgraph` drawings use Rust for layout and SVG geometry; Typst typesets
only their labels. Custom templates, layout sequences, partial-edge styles, and other options
outside the native renderer's supported subset use the complete Typst pipeline. Requested
settings are retained in that path; application payloads stay in Python.

`render()` snapshots the graph, settings, and selector results. `typst_source` retains the
configuration source for inspection. SVG, PDF, PNG, and Typst exports of a native result use
its rendered SVG snapshot, so exporting does not redraw the graph with a different renderer.
Feynman diagrams and tensor networks retain their own default styles; passing shared drawing
overrides is an explicit customization.

The optional theme demonstrates that domain conventions belong in those selectors. It does not
change the workflow records or require them to be serializable. For the complete types and
method signatures, use the #link("reference/python/")[Python API reference].

To experiment with a DOT-backed physics example, continue to
#product-link("gammaloop", label: "GammaLoop DOT input", page: "guides/dot-input/").
For a continuously updating view of the force solver itself, use the
#link("playground/")[live DOT playground].
]
