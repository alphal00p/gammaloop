#import "../../shared.typ": product-link

#let showcases = [
= FeynKit notebook showcases

Six executable Marimo notebooks cover every FeynKit component. Start with a small scalar
model, then choose the physics operation you want to explore. Diagrams, tables, and symbolic
expressions are computed by the same APIs documented in the reference.

== Choose a showcase

- #link("guides/showcases/first-diagram/")[A first diagram] follows the facade from a normalized
  model to generated one-loop graphs and notebook figures.
- #link("guides/showcases/models-and-diagrams/")[Models, generation, and graphs] explores
  parameter cards, process selectors, loop ranges, serialization, and momentum bases.
- #link("guides/showcases/cff/")[Cross-Free Families] inspects orientations, surfaces,
  denominator products, and native Symbolica expressions.
- #link("guides/showcases/kinematics/")[Kinematics and jets] covers momenta, Lorentz
  transformations, angular distances, and jet algorithms.
- #link("guides/showcases/ufo/")[Loading UFO models] normalizes a raw model and inspects
  its diagnostics before generating diagrams.
- #link("guides/showcases/tensor-reduction/")[Vacuum tensor reduction] projects low- and
  high-rank numerators and turns a vacuum graph into scalar contributions.

== Component coverage

#table(
  columns: (1fr, 2fr),
  table.header([Component], [Showcase]),
  [`feynkit` facade], [#link("guides/showcases/first-diagram/")[A first diagram]],
  [`feynkit-py`], [All six notebooks through `symbolica.community.feynkit`],
  [`feynkit-model`], [#link("guides/showcases/models-and-diagrams/")[Models and parameter cards]],
  [`feynkit-generator`], [#link("guides/showcases/models-and-diagrams/")[Processes and generation]],
  [`feynkit-graph`], [#link("guides/showcases/models-and-diagrams/")[Serialization and momentum bases]],
  [`feynkit-cff`], [#link("guides/showcases/cff/")[Cross-Free Families]],
  [`feynkit-kinematics`], [#link("guides/showcases/kinematics/")[Kinematics and jets]],
  [`feynkit-ufo`], [#link("guides/showcases/ufo/")[Loading UFO models]],
  [`feynkit-tensor`], [#link("guides/showcases/tensor-reduction/")[Vacuum tensor reduction]],
)

The Python notebooks exercise the focused Rust components through
`symbolica.community.feynkit`. The #link("reference/interfaces/")[interface guide] links the
Rust and Python references. For tensor-aware notation and algebraic identities, continue to
#product-link("spenso", page: "guides/showcase/", label: "the Spenso + Idenso showcase").

== Running a notebook locally

Use a shared Symbolica host containing FeynKit, as described in the
#link("guides/community-host/")[host guide]. Its Python environment needs `marimo==0.24.0`
and `typst==0.15.0`; the UFO example additionally needs Python 3.11 or newer and the
#link("guides/showcases/ufo/")[pinned Symbolica 3-compatible UFO loader]. Set
`SYMBOLICA_LICENSE` through your local environment when required by your Symbolica
distribution.

Run the following from the checkout, replacing the interpreter path with that environment's
Python. The command opens an editable Marimo notebook and reads the checkout's model fixtures.

// docs-example: syntax
```sh
just notebook feynkit/00_quickstart_marimo /path/to/python
```

Each showcase page links its source and gives the corresponding command. The
#link("guides/notebooks/")[notebook rendering guide] explains the SVG, HTML, Typst, and symbolic
display APIs used by the examples.
]
