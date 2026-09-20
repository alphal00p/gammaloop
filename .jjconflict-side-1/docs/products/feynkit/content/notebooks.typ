#import "../../shared.typ": product-link

#let notebooks = [
= Notebook figures and symbolic output

For complete executable examples, open the #link("guides/showcases/")[FeynKit showcase gallery].

Use a diagram produced by the #link("quickstart/python/")[Python quickstart]. Its
`to_linnest()` method returns complete Typst source; it does not compile a figure. `to_svg()`,
`to_html()`, `_repr_svg_()`, and `_repr_html_()` compile that source with Python's Typst package.
Rendering uses Linnet’s Python preparation and compilation pipeline.
The native module embeds the compatible Linnest/Kurvst source, Wasm files, and shared
physics styles. Particle metadata selects dashed scalar lines, fermion arrows,
photon waves, gluon coils, and mathematical labels from the model's TeX names.
FeynKit, the physics showcase, and `just draw` use the same physics layout template:
100 steps per epoch and 30 epochs by default. Amplitudes place incoming legs on the
left and outgoing legs on the right using Linnest half-edge ordering. Finalized
cross-section diagrams instead pair both sides by their external connection IDs.
The shared template owns label placement, force presets, and particle styles.

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
