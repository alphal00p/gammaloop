#import "../../shared.typ": source-link, product-link

#let python-rendering = [
= Rendering graphs from Python

A Linnet graph can carry ordinary Python objects and render directly in a notebook. This guide
uses a document-processing workflow to explore the default drawing, layout selection, and
styles computed from application data. Begin with the #link("quickstart/python/")[Python
quickstart] for installation and basic graph queries.

The notebook below loads automatically when it comes into view. Python and Linnet
run through WebAssembly; the browser measures MathML labels for native graph drawing.
The first visit downloads their dependencies;
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

`Graph.render()` and Feynman diagram rendering use the shared
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
view = graph.render(config=config)
view
```

The #link("reference/typst/layout/")[layout guide] describes every ImPrEd control.
Python uses underscores where the Typst names use hyphens. Explicit layouts
such as `LayoutAlgorithm.Force` continue to select their respective algorithms.

== The rendering boundary

Returning a `Graph` from a cell invokes its notebook display representation.
`view = graph.render()` captures selectors and settings in a `DiagramView`;
displaying it shows a complete snapshot. When AnyWidget is installed, it also
measures MathML labels in the actual notebook browser and refreshes the drawing
through the existing native `Scene::render` path. Snapshot display does not wait
for font loading; exact measurements still require the bundled font.
Marimo isolates the drawing in a frame and reveals it after its theme, controls
and initial size are ready, including when AnyWidget is absent.
The completed drawing has hover labels, pinned click details and local selection
using Shift-, Ctrl-, or Meta-click. These interactions do not update Python graph
selections automatically. Display `view` directly in Marimo, IPython or Jupyter.
Drag the drawing to pan, or use Ctrl/Meta-scroll to zoom gently around the pointer.
At 100%, small standalone drawings use a bounded enlargement for readable labels.
Larger drawings and explicitly bounded or embedded previews retain their authored
scale. The viewport follows the drawing's height without reserving space for
hidden details or navigation. The zoom percentage and keyboard help float in a
compact corner control over the canvas. Narrow columns or explicitly bounded
previews can shrink the drawing to fit; the zoom indicator reports the scale
relative to its default size.
With the graph focused, `+` and `-` zoom in five-percent steps and `0` fits the
drawing up to that 100% baseline. Hover previews details below standalone drawings
without moving the graph; bounded or embedded previews show details in an overlay.
The preview disappears when the pointer leaves the element;
keyboard focus also previews details. Clicking pins the details until the panel is closed.
The graph's `RenderSettings` combines layout options, drawing defaults, and Python selectors.
`Graph` and `Subgraph` drawings use Rust for layout and SVG geometry, with direct
MathML labels and titles. Generic `TextLabel` styles and `MathSymbol` Greek
identifiers and scripts retain their typed meaning through `IndexDisplay`.
No Typst parser is involved. `LayoutSettings.then(...)` and `LayoutSettings.sequence(...)`
run native passes in order, continuing from the preceding geometry. Unsupported settings
raise an error; Python has no full-document Typst renderer or template fallback.

Keep the view when you want to export the notebook's drawing:

```python
from IPython.display import display

view = graph.render(config=config)
display(view)
```

After measurement completes, run a later cell:

```python
snapshot = view.measured()
svg = snapshot.to_svg()
html = snapshot.to_html()
snapshot.save("workflow.svg")
```

The immutable `DiagramRender` retains the measured boxes and embedded math font.
Its SVG uses XHTML `foreignObject` MathML and is a browser-oriented artifact.
For explicit Marimo embedding after measurement, use `mo.iframe(snapshot.to_html())`.

Scripts explicitly choose `view.measured(browser="chromium")`, requiring
Playwright and Chromium. Alternatively, a browser can measure the probes from
`view.measurement_request()` and return its receipt to
`view.measured(measurements=receipt)`. Receipts belong to that same view;
`graph.to_svg(browser="chromium")` creates a fresh view instead. Pending
`measured()` calls raise an error; strict measurement never uses estimates.

Default `view.snapshot()`, `view.to_svg()` and `view.to_html()` export complete
static output without launching a browser. They reuse completed measurements,
otherwise use optional typst-py outlining when available, otherwise estimated
STIX MathML label bounds. Estimates leave measurement pending. Installing
`anywidget` adds live notebook measurement; earlier snapshots remain immutable.


Use `outlined = view.outlined()` for vector-label SVG through optional Python
typst-py (distribution `typst`). PDF, PNG and Typst SVG import require this
route to retain labels. `outlined.save("workflow.pdf")` converts that SVG, and
`outlined.to_linnest()` embeds it in portable Typst source. The native extension
contains no Rust Typst compiler or graph plugin runtime.
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
