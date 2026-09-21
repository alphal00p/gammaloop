#import "../../shared.typ": product-link

#let notebooks = [
= Notebook figures and symbolic output

For complete executable examples, open the #link("guides/showcases/")[FeynKit showcase gallery].

Use a diagram produced by the #link("quickstart/python/")[Python quickstart]. Its
`to_linnest()` method returns complete Typst source; it does not compile a figure. `to_svg()`,
`to_html()`, `_repr_svg_()`, and `_repr_html_()` compile that source with Python's Typst package.
Rendering uses Linnet’s Python preparation and compilation pipeline.
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
with default strength `gamma-dangling-centroid=1.25`. Neither mode automatically pins X or Y;
explicit positions still take precedence.
The shared template owns label placement, force presets, and particle styles.
Internal and external label distances are independent: ordinary particle labels
use offsets of 0.60 and 0.30 times the graph spring length, respectively; with
momentum labels these become 0.75 and 0.45. These defaults also apply to `just draw`.
Linnest exposes `internal-label-length-scale` and `external-label-length-scale`
(or `labels: (internal-distance: ..., external-distance: ...)` in Typst).
For example, `just draw --input internal-label-length-scale=0.8` increases only
internal spacing; `--input external-label-length-scale=0.4` controls external labels.

// docs-example: syntax
```sh
python -m pip install "typst>=0.15,<0.16"
```

Save the resulting SVG or leave the diagram as the last value in a notebook cell:

// docs-example: compile
```python
from pathlib import Path

Path("diagram.svg").write_text(diagram.to_svg(), encoding="utf-8")
Path("diagram.typ").write_text(diagram.to_linnest(), encoding="utf-8")
diagram
```

The example assumes `diagram` from the quickstart. In Jupyter, its rich representation draws
the figure automatically. In another notebook frontend, display `diagram.to_html()` using that
frontend's HTML object when it does not consume the standard rich-display methods. Model,
generation, and CFF objects also expose compact representations for inspection.

Pass a Linnet selection as `highlight` to draw a region with Linnest's subgraph
highlighting while keeping the complete diagram visible:

// docs-example: compile
```python
graph = diagram.to_linnet()
selected = graph.filter(edge=lambda edge: edge.data.particle_name == "b")
Path("highlighted-diagram.svg").write_text(
    diagram.to_svg(highlight=selected), encoding="utf-8"
)
```

`to_linnest(highlight=selected)` and `to_html(highlight=selected)` accept the same
selection; omitting `highlight` keeps the usual drawing. Selections may identify
individual half-edges as well as complete edges. They must belong to this diagram's
exported graph at its current topology revision; foreign or stale selections raise
an error, as they do for the physics operations.

`diagram.numerator_expression()` returns Spenso’s `TensorExpression`, retaining the tensor interface
and index display hooks. `diagram.build_cff().to_expression()` returns a native Symbolica
expression. Displaying those algebraic results is separate from rendering a graph. Use
#product-link("spenso", page: "guides/python/", label: "Spenso's display tools") for tensor-aware
algebra and #product-link("linnet", page: "guides/python-rendering/", label: "Linnet's rendering guide")
for the underlying graph renderer.

If figure compilation reports a missing `typst` module, install it in the interpreter that runs
the Symbolica host. A model import or successful generation does not require that renderer.
The #link("guides/community-host/")[host guide] explains how a distribution can include this
optional dependency for notebook users.
]
