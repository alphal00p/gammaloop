#import "../../shared.typ": product-link

#let notebooks = [
= Notebook figures and symbolic output

Displaying an `Amplitude` shows compact rows pairing each diagram with its
weighted tensor contribution. Open a row to inspect the graph with Linnet's
pan, zoom, selection, and hover controls. `amplitude.expression()` remains the
full symbolic operator and uses Spenso's tensor printer.

Generated collections use a thumbnail strip and one shared graph viewer.
For `process.generate_cross_section(...)`, these are the actual cut
forward-scattering diagrams. Selecting a thumbnail changes the inspected graph;
it does not alter the diagrams or their momentum routing. Expanded collection
graphs keep a fixed viewport height, with hover previews overlaid and click-to-pin
details. Displays preview at most six diagrams; `.diagrams` exposes the full
collection. Single-diagram SVG and Typst exports remain available independently.


For complete executable examples, open the #link("guides/showcases/")[FeynKit showcase gallery].

Graph geometry, layout, and interaction are produced directly in Rust. The embedded
Typst compiler typesets labels and formulas with bundled fonts; it does not load
MiTeX, Linnest, or Kurvst Typst packages. NumPy is the host's only Python dependency
for these notebook displays.

The notebook compiler retains SVG output, HTML/MathML, text shaping, and math
layout. It omits PDF/PNG export, PDF image import, WebAssembly plugins,
syntax highlighting, bibliographies, JPEG/GIF/WebP decoding, and system font
discovery. Labels can use normal Typst math and text, including bold and italic;
raw text is displayed without syntax coloring. Fonts are bundled, so rendering
works offline without Typst packages or a system installation. General document
compilation remains available in the standalone tools.

== Model labels

Particles accept optional `typstname` and `antitypstname` fields; parameters accept
`typstname`. Supply native Typst math without dollar delimiters, for example
`"typstname": "alpha_s"` or `"antitypstname": "overline(u)"`. Built-in models
supply these labels. Missing labels display the ordinary model name as escaped text.
The separate `texname` and `antitexname` fields still control LaTeX output. Loading
parameter labels changes presentation only, preserving algebraic names and expressions.

== Native SVG configuration

`diagram.render(...)` returns a `DiagramRender` that displays the configured
figure directly in IPython, Jupyter, and Marimo. The result retains its rendered
SVG and labels, so later display and export reuse the same snapshot.
Use `to_svg()` or `to_html()` for text exports, or `to_linnest()` for a
self-contained Typst document embedding the SVG. The exported document needs
no graph packages.

// docs-example: compile
```python
from pathlib import Path

drawing = diagram.render()
Path("diagram.svg").write_text(drawing.to_svg(), encoding="utf-8")
Path("diagram.typ").write_text(drawing.to_linnest(), encoding="utf-8")
drawing
```

Pass `graph.RenderSettings` as `config`. Shared layout, stroke, drawing options, and
`DiagramRender` snapshots live in `symbolica.community.graph`. Named constructors support
editor completion and `help()`. Settings can be changed for future drawings; existing
render results retain their original inputs. Diagram names appear in surrounding captions;
set `title` explicitly to add a heading to the SVG. `DrawOptions.node_radius` uses graph
units; `Stroke.thickness` uses points and `Stroke.paint` accepts a typed `Color`.

`LayoutSettings` groups advanced layout controls such as `impred_steps`,
`impred_spacing`, `impred_repulsion`, and `impred_labels`. Use
`help(LayoutSettings)` for all options and their defaults. Invalid options and
numeric values are rejected when constructing settings.

Particle and momentum presentation uses boolean `hepkit.DiagramStyle` arguments:
`show_particle`, `show_momentum`, `show_edge_index`, `show_node_index`,
`momentum_arrows`, `split_initial_state`, and `debug` (node and edge indices).

// docs-example: compile
```python
from symbolica.community import graph, hepkit

settings = graph.RenderSettings(
    layouts=graph.LayoutSettings(impred_steps=100),
    drawing=graph.DrawOptions(
        edge_stroke=graph.Stroke(paint=graph.Color("#6f4d85"), thickness=1.2),
    ),
)
style = hepkit.DiagramStyle(show_particle=False)
drawing = diagram.render(momenta=True, config=settings, style=style)
Path("momenta.svg").write_text(
    drawing.to_svg(), encoding="utf-8"
)
basis = next(iter(diagram.loop_momentum_bases()))
Path("alternative-routing.svg").write_text(
    diagram.render(lmb=basis, config=settings).to_svg(), encoding="utf-8"
)
drawing
```

`momenta=True` displays the stored routing. Passing `lmb=basis` enables momentum
display with that basis without changing the diagram. Loop and external components
use zero-based `k_i` and `p_i`. Momentum arrows follow source to sink independently
of fermion arrows. Explicit physics options override these display defaults.
A basis from a different diagram is rejected.

Cross sections open initial-state connections into incoming and outgoing legs by
default, retaining final-state cut edges. Set
`RenderSettings(split_initial_state=False)` for the sewn view.
Both views preserve the original edge and half-edge IDs. Configuration is per-call;
it changes neither the physics graph nor its loop-momentum basis.

The example assumes `diagram` from the quickstart. Its rich representation draws the figure
automatically. Hover over a vertex or edge to identify it, or click to pin its details.
Shift-, Ctrl-, or Meta-click toggles that element in the displayed selection; an ordinary
click leaves the selection unchanged. Enter and Space provide the same inspection and selection
actions for a focused element. Escape dismisses the details without clearing the selection.
The original particle
line styles, labels, and transparent background are preserved. Model, generation, and CFF
objects also expose compact representations for inspection.

For explicit Marimo embedding, use `mo.iframe(diagram.to_html())` so the interaction script
runs. Displaying `diagram` directly already uses Marimo's interactive HTML path.
Without scripts, inline SVG still shows native hover labels; embedding it as a static image
shows the drawing. Figure selections are local browser state; they do not alter the Python
diagram or create a physics `Subgraph`. Copy the node and edge IDs from the details panel into
`diagram.subgraph(nodes=[...], edges=[...])` for subsequent symbolic calculations, or use
`diagram.filter(...)` to select by physics properties.

Display a FeynKit `Subgraph` directly to show the selected region in its original
diagram, with the remaining graph muted and dotted:

// docs-example: compile
```python
region = diagram.filter(edge=lambda edge: edge.particle_name == "b")
Path("highlighted-diagram.svg").write_text(
    region.render().to_svg(), encoding="utf-8"
)
region
```

`Subgraph` inherits the physics diagram's renderer and retains shared ownership of
the original diagram, including its particle styles and labels. Selections may identify
individual half-edges as well as complete edges. Use `diagram.subgraph(selection)`
to import a canonical Linnet selection; foreign or stale selections raise an error.
For a standalone drawing of the excised topology, display `region.excise()`.

`diagram.numerator_expression()` returns Spenso’s `TensorExpression`, retaining the tensor interface
and index display hooks. `diagram.build_cff().to_expression()` returns a native Symbolica
expression. Displaying those algebraic results is separate from rendering a graph. Use
#product-link("spenso", page: "guides/python/", label: "Spenso's display tools") for tensor-aware
algebra and #product-link("linnet", page: "guides/python-rendering/", label: "Linnet's rendering guide")
for the underlying graph renderer.

The #link("guides/community-host/")[host guide] describes embedding the renderer
in a distribution. Optional standalone Linnet graph interoperability is separate
from rendering and requires that package only when explicitly used.

== Configured amplitude collections

`amplitude.render()` returns a `graph.DiagramRender` snapshot that displays directly
in IPython, Jupyter, and Marimo. `config` accepts the same `RenderSettings` as
individual diagrams, while `term_settings` controls the weighted tensor notation.
The default preview shows six contributions; `max_diagrams=None` renders all.

```python
from symbolica.community.tensor import DisplaySettings
from symbolica.community.hepkit import RenderSettings

drawing = amplitude.render(
    config=RenderSettings(show_momentum=True),
    max_diagrams=3,
    term_settings=DisplaySettings(show_dimensions=True),
)
drawing
html = drawing.to_html()
svg = drawing.diagrams[0].to_svg()
```

The `diagrams` property contains only the rendered contributions, in amplitude
order. Rendering captures a presentation snapshot and leaves the amplitude intact.
]
