#import "../../shared.typ": product-link

#let notebooks = [
= Notebook figures and symbolic output

For complete executable examples, open the #link("guides/showcases/")[FeynKit showcase gallery].

Use a diagram produced by the #link("quickstart/python/")[Python quickstart]. Its
`to_linnest()` method returns complete Typst source; it does not compile a figure. `render()`,
`to_html()`, `_repr_svg_()`, and `_repr_html_()` compile that source with Python's Typst package.
Rendering uses Linnet’s Python preparation and compilation pipeline. `render()` replaces
the earlier diagram `to_svg()` method and returns SVG text, including hover information.
SVG figures have a transparent background. Their palette follows the browser's
light/dark preference when opened separately; inline figures also follow explicit
Marimo and Jupyter notebook themes. Both modes use the website diagram palette,
including its neutral and charged-particle colours and lightened sink strokes.
The native module embeds the compatible Linnest/Kurvst source, Wasm files, and shared
physics styles. Particle metadata selects dashed scalar lines, fermion arrows,
photon waves, gluon coils, and mathematical labels from the model's TeX names.
FeynKit, the physics showcase, and `just draw` use the same physics layout template:
100 steps per epoch and 30 epochs by default. Both modes group incoming X coordinates on the
left and outgoing X coordinates on the right, and start Y coordinates in half-edge order.
Amplitudes keep Y independent; finalized cross-section diagrams share movable Y groups by
external connection ID (`is_cut`). Dangling-centroid repulsion spreads external endpoints,
with default strength `gamma-dangling-centroid=1.25`. A horizontal spring also biases
incoming endpoints left and outgoing endpoints right of the current node centroid,
with `external-centroid-bias=1.0` relative to the edge spring strength. Its target
distance is `external-centroid-distance=3.0` times the external edge's natural
spring length; Y stays free. The reaction
is shared over the nodes so this force does not translate the whole graph.
Use `just draw --input external-centroid-bias=0` to disable it, or increase that
value for a stronger bias. Generic Linnest layouts default to zero and expose the
same option through Python's `LayoutOptions(external_centroid_bias=...)` and
Typst's `constraints: (external-centroid-bias: ...)`.
Neither mode automatically pins X or Y;
explicit positions still take precedence.
The shared template owns label placement, force presets, and particle styles.
Internal and external label distances are independent: ordinary particle labels
use offsets of 0.60 and 0.45 times the graph spring length, respectively; with
momentum labels these become 0.75 and 0.60. These defaults also apply to `just draw`.
Linnest exposes `internal-label-length-scale` and `external-label-length-scale`
(or `labels: (internal-distance: ..., external-distance: ...)` in Typst).
For example, `just draw --input internal-label-length-scale=0.8` increases only
internal spacing; `--input external-label-length-scale=0.6` controls external labels.

// docs-example: syntax
```sh
python -m pip install "linnet==0.1.0" "typst>=0.15,<0.16"
```

Save the resulting SVG or leave the diagram as the last value in a notebook cell:

// docs-example: compile
```python
from pathlib import Path

Path("diagram.svg").write_text(diagram.render(), encoding="utf-8")
Path("diagram.typ").write_text(diagram.to_linnest(), encoding="utf-8")
diagram
```

Rendering controls use Linnet's existing typed groups. `layouts` contains solver and
spacing settings, `drawing` controls the canvas and geometry, and `style` supplies
node/edge styling. Physics controls belong in `template_options`, using the same
hyphenated names as `just draw --input`: `show-particle`, `show-edge-index`,
`show-node-index`, `show-half-edge-index`, `debug`, `momentum-arrows`,
`momentum-arrow-offset`, `momentum-arrow-length`, `momentum-arrow-side`,
`momentum-arrow-stroke`, and `momentum-label-gap`, among others. Catalogue pagination
(`rows` and `columns`) applies to `just draw`, rather than an individual SVG.

// docs-example: compile
```python
from linnet import DrawOptions, LayoutOptions, RenderConfig

settings = RenderConfig(
    layouts=LayoutOptions(
        seed=42,
        steps=100,
        epochs=30,
        external_centroid_bias=1.5,
        internal_label_length_scale=0.8,
        external_label_length_scale=0.7,
    ),
    drawing=DrawOptions(unit=1.5),
    template_options={"show-particle": False},
)
Path("momenta.svg").write_text(
    diagram.render(momenta=True, config=settings), encoding="utf-8"
)
basis = next(iter(diagram.loop_momentum_bases()))
Path("alternative-routing.svg").write_text(
    diagram.render(lmb=basis, config=settings), encoding="utf-8"
)
```

Internal particle labels use a uniform gap to the measured text box by default;
a deterministic annealing pass slides them along the rendered curves and can
switch sides to reduce overlaps without changing that gap. Explicit `label-side`
choices on path labels stay fixed. Adjust
`LayoutOptions(internal_label_length_scale=...)` to change that gap. External
labels keep their independent spacing. `LabelLayout.FixedGap` exposes the same
mode on generic Linnet layouts; `LabelLayout.DanglingTangent` retains the freely
relaxed internal labels.

Enable the notebook's *Collision boxes* toggle or set
`DrawOptions(debug_label_collisions=True)` to inspect the optimizer's padded
boxes. Dashed purple boxes repel other purple boxes; cyan label boxes repel
orange node/edge boxes. The orange rectangles sample the visible edge carriers,
including the width of waves and coils. Base padding is split equally across each
pair of boxes. `DrawOptions(label_collision_padding=0.15)` adds another 0.15
canvas units on every side of each label box by default. Increase it to encourage
more clearance from labels, nodes, and edges; use zero for the original base
sizes. The notebook's *Extra label collision padding* slider controls this value.
It changes collision avoidance while preserving the fixed normal label gap.
The debug overlay itself preserves layout, canvas size, and SVG interaction.
In Typst, pass `label-collision-padding: 0.15` and
`debug-label-collisions: true` to `draw`.

`momenta=True` draws the graph's stored routing. Passing `lmb=basis` also enables
momentum display, using that basis without changing the diagram. Loop and external
components use the same zero-based `k_i` and `p_i` conventions as
`MomentumSignature.format_momentum()`. Momentum arrows always follow source to sink,
independently of fermion-arrow orientation. A basis from a different diagram is
rejected; region highlighting can use its original diagram's complete basis.
Explicit physics settings override the display defaults enabled by `momenta` or
`lmb`: for example, `momentum-arrows: false` hides arrows while retaining labels,
and `show-momentum: false` hides momentum labels. Edge hover information still
contains the routing. `to_html()` and `to_linnest()` accept the same options.

Cross sections open their initial-state connections into matched incoming and outgoing
legs by default, retaining the physical final-state cut edges. Use
`RenderConfig(template_options={"split-initial-state": False})` for the sewn view.
Both views retain the original diagram's edge and half-edge IDs for hover, selection,
and highlighting; the physics graph and its loop-momentum basis are unchanged.

Configurations are per-call snapshots: neither the diagram nor the supplied
`RenderConfig` is changed. The embedded physics renderer accepts typed layout,
drawing, style and template options; custom templates and Python drawing selectors
belong on the generic Linnet graph's `prepare_render()` API. Imported Typst style
functions are supported through the usual `source_root` and module references.

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
region = diagram.filter(edge=lambda edge: edge.data.particle_name == "b")
Path("highlighted-diagram.svg").write_text(
    region.render(), encoding="utf-8"
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

If figure compilation reports a missing `linnet` or `typst` module, install it in the interpreter that runs
the Symbolica host. A model import or successful generation does not require that renderer.
The #link("guides/community-host/")[host guide] explains how a distribution can include this
optional dependency for notebook users.
]
