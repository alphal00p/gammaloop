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

Graph geometry, layout and interaction use the native Rust renderer. Labels and
formulas use direct MathML. The live notebook browser measures labels through
AnyWidget before the existing native `Scene::render` draws the graph; this path
has no Typst compiler. Graphs, subgraphs, process and generation collections,
amplitudes, CFF and LTD graph explorers share this `DiagramView` lifecycle.
Symbolic expression HTML is complete MathML; paged tensor output uses Symbolica's
shared viewer with Spenso's direct semantic leaves.

== Model labels

MathML particle labels display escaped literal `Particle.name` values. They do
not interpret Typst source. Particles also accept optional `typstname` and
`antitypstname` fields; parameters accept `typstname`. These apply to Typst and
vector exports. Supply native Typst math without dollar delimiters, for example
`"typstname": "alpha_s"` or `"antitypstname": "overline(u)"`. Built-in models
supply these labels. Missing labels display the ordinary model name as escaped text.
The separate `texname` and `antitexname` fields still control LaTeX output. Loading
parameter labels changes presentation only, preserving algebraic names and expressions.

== Native SVG configuration

`diagram.render(...)` returns a `DiagramView`. Displaying it in IPython, Jupyter
or Marimo uses static output by default; with `anywidget`, it also sends the
labels to the actual browser for measurement. Keep that view
when you want to export the drawing:

// docs-example: compile
```python
from IPython.display import display

view = diagram.render()
display(view)
```

After measurement completes, use a later cell:

```python
snapshot = view.measured()
snapshot.save("diagram.svg")
html = snapshot.to_html()
```

The immutable `DiagramRender` retains the measured boxes and embedded math font.
Its measured SVG contains XHTML `foreignObject` MathML and targets browsers.
`measured()` raises while measurement is pending. Scripts explicitly request
`view.measured(browser="chromium")` with optional Playwright and Chromium, or
return a browser receipt from `view.measurement_request()` to
`view.measured(measurements=receipt)`. Receipts must come from the same view;
strict `measured()` never uses estimates or launches a browser implicitly.

Default `view.snapshot()`, `view.to_svg()` and `view.to_html()` export complete
static output without launching a browser. They reuse completed measurements,
otherwise use optional typst-py outlining when available, otherwise estimated
STIX MathML label bounds. Estimates leave measurement pending. Installing
`anywidget` adds live notebook measurement; earlier snapshots remain immutable.


For vector labels, call `outlined = view.outlined()` with optional Python
typst-py (distribution `typst`). This route is required for PDF, PNG and Typst
SVG import to avoid missing MathML labels. Use `outlined.save("diagram.pdf")`
or `outlined.to_linnest()` for a self-contained Typst document embedding the
outlined SVG. The document needs no graph packages.

Pass `graph.RenderSettings` as `config`. Shared settings, `DiagramView` and
`DiagramRender` live in `symbolica.community.graph`. Named constructors support
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
This script example explicitly measures with Chromium before exporting:

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
view = diagram.render(momenta=True, config=settings, style=style)
snapshot = view.measured(browser="chromium")
snapshot.save("momenta.svg")
basis = next(iter(diagram.loop_momentum_bases()))
alternative_view = diagram.render(lmb=basis, config=settings)
alternative_view.measured(browser="chromium").save("alternative-routing.svg")
snapshot
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
line styles, labels, and transparent background are preserved. Model, generation, CFF and LTD
objects also expose compact representations for inspection.

For explicit Marimo embedding after measurement, use `mo.iframe(snapshot.to_html())`.
Displaying `diagram` or `view` directly uses the live widget MIME path.
`diagram.to_html(browser="chromium")` explicitly measures a fresh view for script
export. Its default `browser=None` does not start a browser and it accepts no
measurement receipt; retain `view` for notebook or external-browser receipts.
Inline measured SVG needs browser support for XHTML MathML labels; use an
outlined snapshot in consumers that require vector labels.
Figure selections are local browser state; they do not alter the Python
diagram or create a physics `Subgraph`. Copy the node and edge IDs from the details panel into
`diagram.subgraph(nodes=[...], edges=[...])` for subsequent symbolic calculations, or use
`diagram.filter(...)` to select by physics properties.

Display a FeynKit `Subgraph` directly to show the selected region in its original
diagram, with the remaining graph muted and dotted:

// docs-example: compile
```python
region = diagram.filter(edge=lambda edge: edge.particle_name == "b")
region.render().measured(browser="chromium").save("highlighted-diagram.svg")
region
```

`Subgraph` inherits the physics diagram's renderer and retains shared ownership of
the original diagram, including its particle styles and labels. Selections may identify
individual half-edges as well as complete edges. Use `diagram.subgraph(selection)`
to import a canonical Linnet selection; foreign or stale selections raise an error.
For a standalone drawing of the excised topology, display `region.excise()`.

`diagram.numerator_expression()` returns Spenso’s `TensorExpression`, retaining the tensor interface
and index display hooks. `diagram.integrate_energy(method="cff")` returns `CffRepresentation`;
`diagram.integrate_energy(method="ltd")` returns `LtdRepresentation`. Displaying either result
opens its expression and graph explorer. Calling `to_expression()` instead
returns a native Symbolica expression, with the normalization conventions
documented in the #link("reference/interfaces/")[interface guide]. Use
#product-link("spenso", page: "guides/python/", label: "Spenso's display tools") for tensor-aware
algebra and #product-link("linnet", page: "guides/python-rendering/", label: "Linnet's rendering guide")
for the underlying graph renderer.

The #link("guides/community-host/")[host guide] describes embedding the renderer
in a distribution. Optional standalone Linnet graph interoperability is separate
from rendering and requires that package only when explicitly used.

== Configured amplitude collections

`amplitude.render()` returns a `graph.DiagramView` that displays directly
in IPython, Jupyter and Marimo. `config` accepts the same `RenderSettings` as
individual diagrams, while `term_settings` controls the weighted tensor notation.
The default preview shows six contributions; `max_diagrams=None` renders all.

```python
from symbolica.community.tensor import DisplaySettings
from symbolica.community.hepkit import RenderSettings

view = amplitude.render(
    config=RenderSettings(show_momentum=True),
    max_diagrams=3,
    term_settings=DisplaySettings(show_dimensions=True),
)
snapshot = view.measured(browser="chromium")
html = snapshot.to_html()
svg = snapshot.diagrams[0].to_svg()
```

The view's `diagrams` property contains configured child views; the measured
snapshot contains child snapshots, in amplitude order. In notebooks, display
`view` and retrieve `view.measured()` after completion to use the notebook's
measurements. Rendering leaves the amplitude intact.
]
